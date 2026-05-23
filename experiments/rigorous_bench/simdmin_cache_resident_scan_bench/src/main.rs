//! Hot-RAM 2-bit packed benchmark: ClosedSync, OpenSync, Minimizer (simd-minimizers).
//! Protocol: 3 warmup + 10 recorded trials per scheme (loops 4..13 of 13).
//! Syncmer K=31, s=11, w=21  =>  canonical_*_syncmers(11, 21).
//! Minimizer k=31, w=21       =>  canonical_minimizers(31, 21).

use std::env;
use std::hint::black_box;
use std::time::Instant;

use packed_seq::PackedSeq;
use rayon::prelude::*;
use simd_minimizers::{
    canonical_closed_syncmers, canonical_minimizers, canonical_open_syncmers,
};

const INTERNAL_LOOPS: usize = 13;
const WARMUP_LOOPS: usize = 3;
const RECORDED_LOOPS: usize = INTERNAL_LOOPS - WARMUP_LOOPS;

#[derive(Clone, Copy)]
enum Scheme {
    ClosedSync,
    OpenSync,
    Minimizer,
}

impl Scheme {
    fn name(self) -> &'static str {
        match self {
            Scheme::ClosedSync => "ClosedSync",
            Scheme::OpenSync => "OpenSync",
            Scheme::Minimizer => "Minimizer",
        }
    }
}

struct ThreadStats {
    hits: usize,
    digest: u64,
}

struct RunStats {
    throughput_gbp_s: f64,
    hits: usize,
    digest: u64,
}

fn load_hot_packed(path: &str) -> (Vec<u8>, usize) {
    use std::fs::File;
    use std::io::Read;

    let mut file = File::open(path).unwrap_or_else(|e| panic!("open {path}: {e}"));
    let mut raw = Vec::new();
    file.read_to_end(&mut raw)
        .unwrap_or_else(|e| panic!("read {path}: {e}"));

    let num_bases = raw.len() * 4;
    let mut touch = 0u64;
    for chunk in raw.chunks(4096) {
        touch = touch.wrapping_add(chunk[0] as u64);
    }
    black_box(touch);
    (raw, num_bases)
}

fn pos_capacity_for_thread(num_bases: usize, num_threads: usize) -> usize {
    std::cmp::max(
        1024 * 1024 * 64,
        num_bases / num_threads / 4 + 1024,
    )
}

fn digest_positions(pos: &[u32]) -> u64 {
    let mut d = 0u64;
    for &p in pos {
        d = d.wrapping_add(p as u64);
        d = d.rotate_left(17);
    }
    black_box(d)
}

fn run_scheme(
    scheme: Scheme,
    data: &[u8],
    num_bases: usize,
    lib_s: usize,
    w: usize,
    min_k: usize,
    min_w: usize,
    thread_bufs: &mut [Vec<u32>],
) -> RunStats {
    let num_threads = thread_bufs.len();
    let chunk_size_bases = (num_bases / num_threads / 128) * 128;

    let start = Instant::now();

    let partial: Vec<ThreadStats> = thread_bufs
        .par_iter_mut()
        .enumerate()
        .map(|(t, pos)| {
            pos.clear();

            let start_base = t * chunk_size_bases;
            let thread_num_bases = if t == num_threads - 1 {
                num_bases - start_base
            } else {
                chunk_size_bases
            };

            if thread_num_bases < 512 {
                return ThreadStats {
                    hits: 0,
                    digest: 0,
                };
            }

            let safe_len = thread_num_bases - 256;
            let end_base = start_base + thread_num_bases;
            let slice_end = (end_base / 4).min(data.len());
            let thread_data = &data[(start_base / 4)..slice_end];
            let packed_seq = PackedSeq::from_raw_parts(thread_data, 0, safe_len);

            match scheme {
                Scheme::ClosedSync => {
                    canonical_closed_syncmers(lib_s, w).run(packed_seq, pos);
                }
                Scheme::OpenSync => {
                    canonical_open_syncmers(lib_s, w).run(packed_seq, pos);
                }
                Scheme::Minimizer => {
                    canonical_minimizers(min_k, min_w).run(packed_seq, pos);
                }
            }

            let hits = pos.len();
            black_box(hits);
            black_box(pos.as_ptr());
            black_box(pos.len());

            ThreadStats { hits, digest: 0 }
        })
        .collect();

    let duration = start.elapsed().as_secs_f64();

    let per_thread: Vec<ThreadStats> = thread_bufs
        .par_iter()
        .zip(partial.par_iter())
        .map(|(pos, p)| {
            let digest = digest_positions(pos);
            black_box(digest);
            ThreadStats {
                hits: p.hits,
                digest,
            }
        })
        .collect();

    let hits: usize = per_thread.iter().map(|s| s.hits).sum();
    let digest: u64 = per_thread
        .iter()
        .map(|s| s.digest)
        .fold(0u64, |a, b| a.wrapping_add(b));

    black_box(hits);
    black_box(digest);

    RunStats {
        throughput_gbp_s: num_bases as f64 / duration / 1e9,
        hits,
        digest,
    }
}

fn mean_std(samples: &[f64]) -> (f64, f64) {
    let n = samples.len() as f64;
    let mean = samples.iter().sum::<f64>() / n;
    let var = samples.iter().map(|x| (x - mean).powi(2)).sum::<f64>() / n;
    (mean, var.sqrt())
}

fn bench_scheme(
    scheme: Scheme,
    data: &[u8],
    num_bases: usize,
    lib_s: usize,
    w: usize,
    min_k: usize,
    min_w: usize,
    thread_bufs: &mut [Vec<u32>],
) -> (f64, f64, usize, u64) {
    let mut recorded = Vec::with_capacity(RECORDED_LOOPS);
    let mut last_hits = 0usize;
    let mut last_digest = 0u64;

    for iter in 0..INTERNAL_LOOPS {
        let stats = run_scheme(
            scheme, data, num_bases, lib_s, w, min_k, min_w, thread_bufs,
        );
        if iter >= WARMUP_LOOPS {
            recorded.push(stats.throughput_gbp_s);
        }
        last_hits = stats.hits;
        last_digest = stats.digest;
    }

    let (mean, std) = mean_std(&recorded);
    (mean, std, last_hits, last_digest)
}

fn main() {
    let args: Vec<String> = env::args().collect();
    if args.len() < 3 {
        eprintln!(
            "Usage: {} <hg38_2bit.bin> <threads> [sync_k=31] [s=11] [min_w=21] [mode=0 (All)|1 (Closed)|2 (Open)|3 (Minimizer)]",
            args.first().map(String::as_str).unwrap_or("simd_min_hot_bench")
        );
        std::process::exit(1);
    }

    let path = &args[1];
    let num_threads: usize = args[2].parse().expect("threads");
    let sync_k: usize = args.get(3).map(|s| s.parse().expect("sync_k")).unwrap_or(31);
    let s: usize = args.get(4).map(|s| s.parse().expect("s")).unwrap_or(11);
    let min_w: usize = args.get(5).map(|s| s.parse().expect("min_w")).unwrap_or(21);
    let run_mode: usize = args.get(6).map(|s| s.parse().expect("run_mode")).unwrap_or(0);
    let w = sync_k - s + 1;
    let lib_s = s;
    let min_k = sync_k;

    rayon::ThreadPoolBuilder::new()
        .num_threads(num_threads)
        .build_global()
        .expect("rayon pool");

    let (data, num_bases) = load_hot_packed(path);
    let cap = pos_capacity_for_thread(num_bases, num_threads);
    let mut thread_bufs: Vec<Vec<u32>> = (0..num_threads)
        .map(|_| Vec::with_capacity(cap))
        .collect();

    println!(
        "NOTE|SimdMinHot|sync_k={}|s={}|w={}|min_k={}|min_w={}|bases={}|loops={}|warmup={}|recorded={}",
        sync_k, s, w, min_k, min_w, num_bases, INTERNAL_LOOPS, WARMUP_LOOPS, RECORDED_LOOPS
    );

    let schemes = match run_mode {
        1 => vec![Scheme::ClosedSync],
        2 => vec![Scheme::OpenSync],
        3 => vec![Scheme::Minimizer],
        _ => vec![Scheme::ClosedSync, Scheme::OpenSync, Scheme::Minimizer],
    };

    for scheme in schemes {
        let (mean, std, hits, digest) = bench_scheme(
            scheme,
            &data,
            num_bases,
            lib_s,
            w,
            min_k,
            min_w,
            &mut thread_bufs,
        );
        let coverage = hits as f64 / num_bases as f64;
        println!(
            "INNER|SimdMinHot|{}|{}|{:.6}|{:.6}|{}|{:.8}|{:016x}",
            scheme.name(),
            num_threads,
            mean,
            std,
            hits,
            coverage,
            digest
        );
        println!(
            "RESULT|SimdMinHot|{}|{}|{:.4}|{:.4}|{}|{:.6}|{:x}",
            scheme.name(),
            num_threads,
            mean,
            std,
            hits,
            coverage,
            digest
        );
    }
}
