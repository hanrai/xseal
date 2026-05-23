//! Isolated s-mer hash micro-benchmark (simd-minimizers / seq-hash): NtHash SIMD rolling.
//! Syncmer K=31, w=21 => s-mer length s=11. Protocol matches hash_isolate_xseal.

use std::env;
use std::fs::File;
use std::hint::black_box;
use std::io::Read;
use std::time::Instant;

use packed_seq::PackedSeq;

/// packed-seq requires 48 trailing bytes for SIMD lookahead (see PackedSeq::from_raw_parts).
const PACKED_SEQ_PAD: usize = 48;
use rayon::prelude::*;
use seq_hash::{KmerHasher, NtHasher};
use wide::u32x8;

const INTERNAL_LOOPS: usize = 13;
const WARMUP_LOOPS: usize = 3;
const RECORDED_LOOPS: usize = INTERNAL_LOOPS - WARMUP_LOOPS;

struct RunStats {
    throughput_gbp_s: f64,
    hash_count: usize,
    digest: u64,
}

fn touch_hot_pages(data: &[u8]) {
    let mut touch = 0u64;
    for chunk in data.chunks(4096) {
        touch = touch.wrapping_add(chunk[0] as u64);
    }
    black_box(touch);
}

fn load_or_generate(
    path: Option<&str>,
    num_bases: usize,
    seed: u64,
) -> (Vec<u8>, usize) {
    if let Some(p) = path {
        if p != "-" {
            let mut f = File::open(p).unwrap_or_else(|e| panic!("open {p}: {e}"));
            let mut raw = Vec::new();
            f.read_to_end(&mut raw).expect("read");
            let file_bases = raw.len() * 4;
            let use_bases = if num_bases == 0 {
                file_bases
            } else {
                num_bases.min(file_bases)
            };
            raw.truncate((use_bases + 3) / 4);
            raw.resize(raw.len() + PACKED_SEQ_PAD, 0);
            touch_hot_pages(&raw);
            return (raw, use_bases);
        }
    }

    assert!(num_bases > 0, "num_bases required for random input");
    let nbytes = (num_bases + 3) / 4;
    let mut raw = vec![0u8; nbytes + PACKED_SEQ_PAD];
    let mut rng = seed;
    for b in raw.iter_mut() {
        rng = rng.wrapping_mul(6364136223846793005).wrapping_add(1);
        *b = (rng >> 33) as u8;
    }
    touch_hot_pages(&raw);
    (raw, num_bases)
}

fn hash_scan_slice(data: &[u8], base_start: usize, num_bases: usize, s: usize) -> (usize, u64) {
    let hasher = NtHasher::<true>::new(s);
    let byte_off = base_start / 4;
    let base_off = base_start % 4;
    let seq = PackedSeq::from_raw_parts(&data[byte_off..], base_off, num_bases);
    let expected = num_bases.saturating_sub(s).saturating_add(1);
    let padded = hasher.hash_kmers_simd(seq, 0);

    let mut simd_digest = u32x8::default();
    let mut count = 0usize;
    let mut it = padded.it;

    let full_chunks = expected / 8;
    for _ in 0..full_chunks {
        if let Some(chunk) = it.next() {
            simd_digest ^= chunk;
            count += 8;
        } else {
            break;
        }
    }

    let arr = simd_digest.to_array();
    let mut digest = arr.iter().fold(0u64, |acc, &h| acc ^ (h as u64));

    let tail_len = expected - count;
    if tail_len > 0 {
        if let Some(chunk) = it.next() {
            let arr = chunk.to_array();
            for i in 0..tail_len {
                digest ^= arr[i] as u64;
                count += 1;
            }
        }
    }

    (count, digest)
}

fn run_once(data: &[u8], num_bases: usize, s: usize, num_threads: usize) -> RunStats {
    let num_hashes = num_bases.saturating_sub(s).saturating_add(1);
    let start = Instant::now();

    let partial: Vec<(usize, u64)> = (0..num_threads)
        .into_par_iter()
        .map(|t| {
            let per = num_hashes / num_threads;
            let start_h = t * per;
            let len_h = if t + 1 == num_threads {
                num_hashes - start_h
            } else {
                per
            };
            if len_h == 0 {
                return (0, 0);
            }
            let base_start = start_h;
            let base_len = len_h + s - 1;
            hash_scan_slice(data, base_start, base_len, s)
        })
        .collect();

    let sec = start.elapsed().as_secs_f64();
    let hash_count: usize = partial.iter().map(|(c, _)| c).sum();
    let digest = partial.iter().fold(0u64, |d, (_, h)| d ^ h);
    let gbp_s = num_bases as f64 / sec / 1e9;
    RunStats {
        throughput_gbp_s: gbp_s,
        hash_count,
        digest,
    }
}

fn main() {
    let args: Vec<String> = env::args().collect();
    if args.len() < 6 {
        eprintln!(
            "Usage: {} <2bit.bin|-> <num_bases|0> <threads> <sync_k> <s> [seed]",
            args[0]
        );
        eprintln!("  sync_k is ignored (hash length = s). NtHash SIMD rolling.");
        std::process::exit(1);
    }

    let path_arg = &args[1];
    let num_bases_arg: usize = args[2].parse().expect("num_bases");
    let num_threads: usize = args[3].parse().expect("threads");
    let _sync_k: usize = args[4].parse().expect("sync_k");
    let s: usize = args[5].parse().expect("s");
    let seed: u64 = args.get(6).map(|x| x.parse().expect("seed")).unwrap_or(42);

    let path_opt = if path_arg == "-" {
        None
    } else {
        Some(path_arg.as_str())
    };
    let (data, num_bases) = load_or_generate(path_opt, num_bases_arg, seed);

    println!(
        "NOTE|tool|SimdHashIsolate|s|{s}|num_bases|{num_bases}|threads|{num_threads}|seed|{seed}"
    );

    rayon::ThreadPoolBuilder::new()
        .num_threads(num_threads)
        .build_global()
        .expect("thread pool");

    for loop_idx in 0..INTERNAL_LOOPS {
        let st = run_once(&data, num_bases, s, num_threads);
        if loop_idx >= WARMUP_LOOPS {
            println!(
                "INNER|loop|{loop_idx}|gbp_s|{:.6}|hashes|{}|digest|0x{:x}",
                st.throughput_gbp_s, st.hash_count, st.digest
            );
        }
    }

    let mut recorded = Vec::with_capacity(RECORDED_LOOPS);
    let mut final_digest = 0u64;
    let mut final_hashes = 0usize;
    for _ in 0..RECORDED_LOOPS {
        let st = run_once(&data, num_bases, s, num_threads);
        recorded.push(st.throughput_gbp_s);
        final_digest = st.digest;
        final_hashes = st.hash_count;
    }
    let mean = recorded.iter().sum::<f64>() / recorded.len() as f64;
    let var = recorded
        .iter()
        .map(|v| {
            let d = v - mean;
            d * d
        })
        .sum::<f64>()
        / recorded.len() as f64;
    let stddev = var.sqrt();

    println!(
        "RESULT|SimdHashIsolate|{num_threads}|{mean:.6}|{stddev:.6}|{num_bases}|{final_hashes}|0x{final_digest:x}"
    );
}
