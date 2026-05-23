use std::env;
use std::fs::File;
use std::time::Instant;
use memmap2::MmapOptions;
use simd_minimizers::{canonical_minimizers, canonical_closed_syncmers, canonical_open_syncmers};
use packed_seq::PackedSeq;
use rayon::prelude::*;

struct BenchResult {
    throughput: f64,
    hits: usize,
}

fn run_benchmark(data: &[u8], num_bases: usize, k: usize, w: usize, mode: &str, num_threads: usize) -> BenchResult {
    let pool = rayon::ThreadPoolBuilder::new().num_threads(num_threads).build().unwrap();
    let start = Instant::now();
    let result: usize = pool.install(|| {
        let chunk_size_bases = (num_bases / num_threads / 128) * 128;
        (0..num_threads).into_par_iter().map(|t| {
            let start_base = t * chunk_size_bases;
            let thread_num_bases = if t == num_threads - 1 { num_bases - start_base } else { chunk_size_bases };
            if thread_num_bases < 512 { return 0; }
            let safe_len = thread_num_bases - 256;
            let end_base = start_base + thread_num_bases;
            let slice_end = (end_base / 4).min(data.len());
            let thread_data = &data[(start_base / 4)..slice_end];
            let packed_seq = PackedSeq::from_raw_parts(thread_data, 0, safe_len);
            let mut pos = Vec::with_capacity(1024 * 1024 * 32);
            match mode {
                "ClosedSync" => { canonical_closed_syncmers(k, w).run(packed_seq, &mut pos); },
                "OpenSync" => { canonical_open_syncmers(k, w).run(packed_seq, &mut pos); },
                "Minimizer" => { canonical_minimizers(k, w).run(packed_seq, &mut pos); },
                _ => panic!("Unknown mode"),
            }
            pos.len()
        }).sum()
    });
    let duration = start.elapsed();
    BenchResult {
        throughput: num_bases as f64 / duration.as_secs_f64() / 1e9,
        hits: result,
    }
}

fn main() {
    let args: Vec<String> = env::args().collect();
    if args.len() < 6 { return; }

    let file = File::open(&args[1]).expect("failed to open file");
    let mmap = unsafe { MmapOptions::new().map(&file).expect("failed to mmap") };
    let data: &[u8] = &mmap;
    let num_bases = data.len() * 4;
    let num_threads: usize = args[2].parse().unwrap();
    let k: usize = args[3].parse().unwrap();
    let s: usize = args[4].parse().unwrap();
    let w: usize = args[5].parse().unwrap();

    // Warmup
    for _ in 0..3 { run_benchmark(data, num_bases, k, k - s + 1, "ClosedSync", num_threads); }

    let mut closed_t = Vec::new();
    let mut open_t = Vec::new();
    let mut min_t = Vec::new();
    let (mut c_hits, mut o_hits, mut m_hits) = (0, 0, 0);

    for _ in 0..10 {
        let r = run_benchmark(data, num_bases, k, k - s + 1, "ClosedSync", num_threads);
        closed_t.push(r.throughput); c_hits = r.hits;
        let r = run_benchmark(data, num_bases, k, k - s + 1, "OpenSync", num_threads);
        open_t.push(r.throughput); o_hits = r.hits;
        let r = run_benchmark(data, num_bases, k, w, "Minimizer", num_threads);
        min_t.push(r.throughput); m_hits = r.hits;
    }

    let stats = |v: &Vec<f64>| {
        let mean = v.iter().sum::<f64>() / v.len() as f64;
        let std = (v.iter().map(|&x| (x - mean).powi(2)).sum::<f64>() / v.len() as f64).sqrt();
        (mean, std)
    };

    let (c_m, c_s) = stats(&closed_t);
    let (o_m, o_s) = stats(&open_t);
    let (m_m, m_s) = stats(&min_t);

    // Match rigorous_bench/xseal_bench.cpp RESULT lines for downstream parsers.
    println!(
        "RESULT|SimdMin|Syncmer|{}|{:.4}|{:.4}|{}|0",
        num_threads, c_m, c_s, c_hits
    );
    println!(
        "RESULT|SimdMin|OpenSync|{}|{:.4}|{:.4}|{}|0",
        num_threads, o_m, o_s, o_hits
    );
    println!(
        "RESULT|SimdMin|Minimizer|{}|{:.4}|{:.4}|{}|0",
        num_threads, m_m, m_s, m_hits
    );

    if args.len() >= 7 {
        use std::io::Write;
        let mut f = File::create(&args[6]).unwrap();
        // Dump Minimizer hits (r3) for conservation comparison
        // We need to run it again or store it. Let's just run it one more time to dump.
        let pool = rayon::ThreadPoolBuilder::new().num_threads(num_threads).build().unwrap();
        pool.install(|| {
            let mut pos = Vec::new();
            let packed_seq = PackedSeq::from_raw_parts(data, 0, num_bases - 256);
            simd_minimizers::canonical_minimizers(k, w).run(packed_seq, &mut pos);
            for p in pos {
                f.write_all(&(p as u32).to_le_bytes()).unwrap();
            }
        });
        println!("INFO|Dumped positions to {}", args[6]);
    }
}
