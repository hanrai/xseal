//! End-to-end baseline: Needletail (parse + per-record `seq` materialization) + simd-minimizers.
//!
//! Use `--preload` to read the FASTA into RAM once, then time parsing from a `Cursor` so
//! throughput is not dominated by disk (fair comparison to memory-resident micro-benchmarks).
//! Use `--rigorous` for 3 warmup + 10 trials and mean ± std on Gbp/s (requires `--preload`).
//! Machine-readable lines: `PHASE|...`, `RESULT|...`, `SPLIT|...`.

use anyhow::{bail, Result};
use clap::Parser;
use needletail::{parse_fastx_file, parse_fastx_reader};
use simd_minimizers::packed_seq::AsciiSeq;
use std::fs;
use std::io::Cursor;
use std::time::Instant;
use rayon::prelude::*;

#[derive(Parser, Debug)]
struct Args {
    path: String,
    #[arg(short, long, default_value_t = 31)]
    k: usize,
    #[arg(short, long, default_value_t = 11)]
    s: usize,
    #[arg(short, long, default_value_t = 21)]
    w: usize,
    #[arg(short, long, default_value_t = 12)]
    threads: usize,
    /// Read entire FASTA into RAM once, then parse from memory (excludes disk read from timed parse).
    #[arg(long, default_value_t = false)]
    preload: bool,
    /// 3 warmup + 10 timed trials; mean ± std for Gbp/s (implies `--preload`).
    #[arg(long, default_value_t = false)]
    rigorous: bool,
}

fn parse_materialize_from_bytes(data: &[u8]) -> Result<(Vec<Vec<u8>>, usize)> {
    let mut reader = parse_fastx_reader(Cursor::new(data))?;
    let mut records = Vec::new();
    let mut total_bases = 0usize;
    while let Some(record) = reader.next() {
        let rec = record?;
        let seq = rec.seq().into_owned();
        total_bases += seq.len();
        records.push(seq);
    }
    Ok((records, total_bases))
}

fn parse_materialize_from_path(path: &str) -> Result<(Vec<Vec<u8>>, usize)> {
    let mut reader = parse_fastx_file(path)?;
    let mut records = Vec::new();
    let mut total_bases = 0usize;
    while let Some(record) = reader.next() {
        let rec = record?;
        let seq = rec.seq().into_owned();
        total_bases += seq.len();
        records.push(seq);
    }
    Ok((records, total_bases))
}

fn scan_only(records: &[Vec<u8>], s: usize, w: usize) -> usize {
    records
        .par_iter()
        .map(|seq_bytes| {
            if seq_bytes.len() < s + w - 1 {
                return 0;
            }
            let mut min_pos = Vec::new();
            simd_minimizers::canonical_closed_syncmers(s, w)
                .run(AsciiSeq(seq_bytes), &mut min_pos);
            min_pos.len()
        })
        .sum()
}

fn mean_std(v: &[f64]) -> (f64, f64) {
    if v.is_empty() {
        return (0.0, 0.0);
    }
    let n = v.len() as f64;
    let mean = v.iter().sum::<f64>() / n;
    let var = v.iter().map(|x| (x - mean).powi(2)).sum::<f64>() / n;
    (mean, var.sqrt())
}

fn main() -> Result<()> {
    let args = Args::parse();
    if args.rigorous && !args.preload {
        bail!("--rigorous requires --preload (each trial re-parses from RAM; disk would dominate noise)");
    }

    rayon::ThreadPoolBuilder::new()
        .num_threads(args.threads)
        .build_global()?;

    println!(
        "NOTE|rust_sota_bench|Needletail parse + seq materialize + simd-minimizers ({} threads, k/s/w={}/{}/{})",
        args.threads, args.k, args.s, args.w
    );

    if args.rigorous {
        let t_disk = Instant::now();
        let bytes = fs::read(&args.path)?;
        let disk_secs = t_disk.elapsed().as_secs_f64();
        let gbp = bytes.len() as f64 / (1024.0 * 1024.0 * 1024.0) / disk_secs;
        println!(
            "PHASE|read_file_to_ram|s|{:.4}|bytes|{}|GiB/s|{:.4}",
            disk_secs,
            bytes.len(),
            gbp
        );

        const WARMUP: usize = 3;
        const TRIALS: usize = 10;
        let mut e2e_gbps = Vec::with_capacity(TRIALS);
        let mut parse_gbps = Vec::with_capacity(TRIALS);
        let mut scan_gbps = Vec::with_capacity(TRIALS);
        let mut frac_parse_wall = Vec::with_capacity(TRIALS);
        let mut total_bases = 0usize;
        let mut total_hits = 0usize;

        for i in 0..WARMUP + TRIALS {
            let t0 = Instant::now();
            let (records, bases) = parse_materialize_from_bytes(&bytes)?;
            let t_parse = t0.elapsed();
            let t1 = Instant::now();
            let hits = scan_only(&records, args.s, args.w);
            let t_scan = t1.elapsed();

            total_bases = bases;
            total_hits = hits;

            if i >= WARMUP {
                let b = bases as f64 / 1e9;
                let ps = t_parse.as_secs_f64();
                let ss = t_scan.as_secs_f64();
                e2e_gbps.push(b / (ps + ss));
                parse_gbps.push(b / ps);
                scan_gbps.push(b / ss);
                frac_parse_wall.push(ps / (ps + ss));
            }
        }

        let (m_e2e, s_e2e) = mean_std(&e2e_gbps);
        let (m_p, s_p) = mean_std(&parse_gbps);
        let (m_s, s_s) = mean_std(&scan_gbps);
        let (m_fp, s_fp) = mean_std(&frac_parse_wall);

        println!("CALIB|logical_bases|{}", total_bases);
        println!("RESULT|TotalHits|{}", total_hits);
        println!(
            "RESULT|E2E_Gbp_s|mean|{:.4}|std|{:.4}|trials|{}",
            m_e2e, s_e2e, TRIALS
        );
        println!(
            "RESULT|ParseMaterialize_Gbp_s|mean|{:.4}|std|{:.4}",
            m_p, s_p
        );
        println!("RESULT|SyncmerScan_Gbp_s|mean|{:.4}|std|{:.4}", m_s, s_s);
        println!(
            "SPLIT|wall_time_fraction_in_parse_mean|{:.4}|std|{:.4}",
            m_fp, s_fp
        );
        return Ok(());
    }

    // Single-run modes (e2e Gbp/s = parse+scan only; disk read for --preload is reported separately)
    let (_, total_bases, t_parse, t_scan, total_hits) = if args.preload {
        let t_disk = Instant::now();
        let bytes = fs::read(&args.path)?;
        let disk_secs = t_disk.elapsed().as_secs_f64();
        println!(
            "PHASE|read_file_to_ram|s|{:.4}|bytes|{}",
            disk_secs,
            bytes.len()
        );

        let t0 = Instant::now();
        let (records, bases) = parse_materialize_from_bytes(&bytes)?;
        let t_parse = t0.elapsed();

        let t1 = Instant::now();
        let hits = scan_only(&records, args.s, args.w);
        let t_scan = t1.elapsed();
        ((), bases, t_parse, t_scan, hits)
    } else {
        let t0 = Instant::now();
        let (records, bases) = parse_materialize_from_path(&args.path)?;
        let t_parse = t0.elapsed();

        let t1 = Instant::now();
        let hits = scan_only(&records, args.s, args.w);
        let t_scan = t1.elapsed();
        ((), bases, t_parse, t_scan, hits)
    };

    let parse_secs = t_parse.as_secs_f64();
    let scan_secs = t_scan.as_secs_f64();
    let total_secs = parse_secs + scan_secs;
    let b = total_bases as f64 / 1e9;

    println!("----------------------------------------------------");
    println!(" Total Bases      : {}", total_bases);
    println!(" Total Hits       : {}", total_hits);
    println!(" Phase parse+mat. : {:.4} s  ({:.4} Gbp/s)", parse_secs, b / parse_secs);
    println!(" Phase scan       : {:.4} s  ({:.4} Gbp/s)", scan_secs, b / scan_secs);
    println!(" Wall parse+scan : {:.4} s  ({:.4} Gbp/s e2e)", total_secs, b / total_secs);
    println!(
        "RESULT|E2E_Gbp_s|single_run|{:.4}|preload|{}",
        b / total_secs,
        args.preload
    );
    println!("----------------------------------------------------");

    Ok(())
}
