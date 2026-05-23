//! 128 KiB 2-bit packed byte blocks -> rayon workers.
//! Per block: scan (overwrite) -> values_u64 harvest (hash) or hit count only.
//! Thread-local tallies cache-line padded; buffers allocated once and reused.

use std::env;
use std::hint::black_box;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::time::Instant;

use packed_seq::PackedSeq;
use rayon::prelude::*;
use simd_minimizers::{
    canonical_closed_syncmers, canonical_minimizers, canonical_open_syncmers,
};

const BLOCK_PACKED_BYTES: usize = 128 * 1024;
const SIMD_CHUNK_BASES: usize = 128;
/// packed-seq requires extra bytes after each slice (see packed_seq::PADDING).
const PACKED_SEQ_PAD_BYTES: usize = 48;
const CACHE_LINE: usize = 64;
const POS_CAP: usize = 512 * 1024;

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

struct BlockGeom {
    num_blocks: usize,
    tail_bytes: usize,
}

struct RunStats {
    throughput_gbp_s: f64,
    hits: usize,
    digest: u64,
}

/// Cache-line-isolated per-thread counters (sum after the timed region).
#[repr(align(64))]
struct ThreadTally {
    hits: usize,
    digest: u64,
    _pad: [u8; CACHE_LINE - std::mem::size_of::<usize>() - std::mem::size_of::<u64>()],
}

struct WorkerSlot {
    pos: Vec<u32>,
    hashes: Vec<u64>,
}

impl WorkerSlot {
    fn new() -> Self {
        Self {
            pos: Vec::with_capacity(POS_CAP),
            hashes: Vec::with_capacity(POS_CAP),
        }
    }
}

fn block_geometry(file_bytes: usize) -> BlockGeom {
    let tail_bytes = file_bytes % BLOCK_PACKED_BYTES;
    let num_blocks = file_bytes / BLOCK_PACKED_BYTES + usize::from(tail_bytes > 0);
    BlockGeom {
        num_blocks,
        tail_bytes,
    }
}

fn block_layout(geom: &BlockGeom, block_id: usize, file_bytes: usize) -> (usize, usize, usize) {
    let byte_off = block_id * BLOCK_PACKED_BYTES;
    let mut block_bytes = if block_id + 1 < geom.num_blocks {
        BLOCK_PACKED_BYTES
    } else if geom.tail_bytes > 0 {
        geom.tail_bytes
    } else {
        BLOCK_PACKED_BYTES
    };
    let avail = file_bytes.saturating_sub(byte_off);
    block_bytes = block_bytes.min(avail);
    let scan_bases = (block_bytes * 4 / SIMD_CHUNK_BASES) * SIMD_CHUNK_BASES;
    (byte_off, block_bytes, scan_bases)
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
    // packed-seq: need PADDING bytes after last accessed base (see from_raw_parts assert).
    raw.resize(raw.len() + PACKED_SEQ_PAD_BYTES + 64, 0);
    (raw, num_bases)
}

fn fold_digest_pos(d: &mut u64, pos: &[u32]) {
    for &p in pos {
        *d = d.wrapping_add(p as u64);
        *d = d.rotate_left(17);
    }
}

fn fold_digest_pos_hash(d: &mut u64, pos: &[u32], hash: &[u64]) {
    for (&p, &h) in pos.iter().zip(hash.iter()) {
        *d = d.wrapping_add(p as u64 ^ h);
        *d = d.rotate_left(17);
    }
}

fn process_block(
    scheme: Scheme,
    data: &[u8],
    byte_off: usize,
    scan_bases: usize,
    lib_s: usize,
    w: usize,
    min_k: usize,
    min_w: usize,
    harvest_hash: bool,
    slot: &mut WorkerSlot,
) -> usize {
    // packed-seq needs tail room inside the logical length (see simd_min_hot_bench).
    const TAIL_MARGIN_BASES: usize = 256;
    if scan_bases < TAIL_MARGIN_BASES + 256 {
        return 0;
    }
    let safe_len = scan_bases - TAIL_MARGIN_BASES;
    let offset_bp = byte_off * 4;
    debug_assert!(offset_bp + safe_len + PACKED_SEQ_PAD_BYTES * 4 <= data.len() * 4);
    let packed_seq = PackedSeq::from_raw_parts(data, offset_bp, safe_len);

    // Reuse capacity; next scan overwrites prefix [0..hits). No clear / memset.
    slot.pos.truncate(0);
    if harvest_hash {
        slot.hashes.truncate(0);
    }

    let hits = match scheme {
        Scheme::ClosedSync => {
            let out = canonical_closed_syncmers(lib_s, w).run(packed_seq, &mut slot.pos);
            if harvest_hash {
                slot.hashes.extend(out.values_u64());
            }
            slot.pos.len()
        }
        Scheme::OpenSync => {
            let out = canonical_open_syncmers(lib_s, w).run(packed_seq, &mut slot.pos);
            if harvest_hash {
                slot.hashes.extend(out.values_u64());
            }
            slot.pos.len()
        }
        Scheme::Minimizer => {
            let out = canonical_minimizers(min_k, min_w).run(packed_seq, &mut slot.pos);
            if harvest_hash {
                slot.hashes.extend(out.values_u64());
            }
            slot.pos.len()
        }
    };

    black_box(slot.pos.as_ptr());
    black_box(slot.hashes.as_ptr());
    hits
}

fn run_scheme_once(
    scheme: Scheme,
    data: &[u8],
    file_bytes: usize,
    num_bases: usize,
    geom: &BlockGeom,
    slots: &mut [WorkerSlot],
    tallies: &mut [ThreadTally],
    harvest_hash: bool,
    lib_s: usize,
    w: usize,
    min_k: usize,
    min_w: usize,
) -> RunStats {
    let next_block = AtomicUsize::new(0);

    for tl in tallies.iter_mut() {
        tl.hits = 0;
        tl.digest = 0;
    }

    let start = Instant::now();

    slots
        .par_iter_mut()
        .zip(tallies.par_iter_mut())
        .for_each(|(slot, tally)| {
            loop {
                let block_id = next_block.fetch_add(1, Ordering::Relaxed);
                if block_id >= geom.num_blocks {
                    break;
                }
                let (byte_off, _, scan_bases) = block_layout(geom, block_id, file_bytes);
                let n = process_block(
                    scheme,
                    data,
                    byte_off,
                    scan_bases,
                    lib_s,
                    w,
                    min_k,
                    min_w,
                    harvest_hash,
                    slot,
                );
                if n == 0 {
                    continue;
                }
                tally.hits += n;
                if harvest_hash {
                    fold_digest_pos_hash(
                        &mut tally.digest,
                        &slot.pos[..n],
                        &slot.hashes[..n],
                    );
                } else {
                    fold_digest_pos(&mut tally.digest, &slot.pos[..n]);
                }
            }
        });

    let duration = start.elapsed().as_secs_f64();

    let hits: usize = tallies.iter().map(|t| t.hits).sum();
    let digest: u64 = tallies
        .iter()
        .map(|t| t.digest)
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
    file_bytes: usize,
    num_bases: usize,
    geom: &BlockGeom,
    slots: &mut [WorkerSlot],
    tallies: &mut [ThreadTally],
    harvest_hash: bool,
    lib_s: usize,
    w: usize,
    min_k: usize,
    min_w: usize,
) -> (f64, f64, usize, u64) {
    let mut recorded = Vec::with_capacity(RECORDED_LOOPS);
    let mut last_hits = 0usize;
    let mut last_digest = 0u64;

    for iter in 0..INTERNAL_LOOPS {
        let stats = run_scheme_once(
            scheme,
            data,
            file_bytes,
            num_bases,
            geom,
            slots,
            tallies,
            harvest_hash,
            lib_s,
            w,
            min_k,
            min_w,
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
    if args.len() < 4 {
        eprintln!(
            "Usage: {} <hg38_2bit.bin> <threads> <hash_mode=0|1> [sync_k=31] [s=11] [min_w=21]",
            args.first()
                .map(String::as_str)
                .unwrap_or("simd_min_block128k_bench")
        );
        std::process::exit(1);
    }

    let path = &args[1];
    let num_threads: usize = args[2].parse().expect("threads");
    let hash_mode: usize = args[3].parse().expect("hash_mode");
    let harvest_hash = hash_mode != 0;
    let sync_k: usize = args.get(4).map(|s| s.parse().expect("sync_k")).unwrap_or(31);
    let s: usize = args.get(5).map(|s| s.parse().expect("s")).unwrap_or(11);
    let min_w: usize = args.get(6).map(|s| s.parse().expect("min_w")).unwrap_or(21);
    let w = sync_k - s + 1;
    let lib_s = s;
    let min_k = sync_k;

    let tool = if harvest_hash {
        "SimdMinBlock128k"
    } else {
        "SimdMinBlock128kCnt"
    };

    rayon::ThreadPoolBuilder::new()
        .num_threads(num_threads)
        .build_global()
        .expect("rayon pool");

    let (data, num_bases) = load_hot_packed(path);
    let file_bytes = data.len();
    let geom = block_geometry(file_bytes);

    let mut slots: Vec<WorkerSlot> = (0..num_threads).map(|_| WorkerSlot::new()).collect();
    let mut tallies: Vec<ThreadTally> = (0..num_threads)
        .map(|_| ThreadTally {
            hits: 0,
            digest: 0,
            _pad: [0u8; CACHE_LINE - std::mem::size_of::<usize>() - std::mem::size_of::<u64>()],
        })
        .collect();

    println!(
        "NOTE|{tool}|protocol=block128k_dispatch|block_packed_bytes={BLOCK_PACKED_BYTES}|scan_bases_per_full_block={}|num_blocks={}|hash_mode={hash_mode}|sync_k={sync_k}|s={s}|w={w}|min_k={min_k}|min_w={min_w}|bases={num_bases}|loops={INTERNAL_LOOPS}|warmup={WARMUP_LOOPS}|recorded={RECORDED_LOOPS}|tally_stride_bytes={}",
        BLOCK_PACKED_BYTES * 4,
        geom.num_blocks,
        std::mem::size_of::<ThreadTally>()
    );

    for scheme in [
        Scheme::ClosedSync,
        Scheme::OpenSync,
        Scheme::Minimizer,
    ] {
        let (mean, std, hits, digest) = bench_scheme(
            scheme,
            &data,
            file_bytes,
            num_bases,
            &geom,
            &mut slots,
            &mut tallies,
            harvest_hash,
            lib_s,
            w,
            min_k,
            min_w,
        );
        let coverage = hits as f64 / num_bases as f64;
        println!(
            "INNER|{tool}|{}|{num_threads}|{mean:.6}|{std:.6}|{hits}|{coverage:.8}|{digest:016x}",
            scheme.name()
        );
        println!(
            "RESULT|{tool}|{}|{num_threads}|{mean:.4}|{std:.4}|{hits}|{coverage:.6}|{digest:x}",
            scheme.name()
        );
    }
}
