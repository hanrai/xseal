//! Marker semantics audit for hot 2-bit benchmarks.
//!
//! Validates simd-minimizers (scalar==SIMD) and checks each reported position
//! against closed / open (strict) / open (relaxed, XSeal-style) predicates.

use std::fs::File;
use std::io::Read;
use std::path::PathBuf;

use clap::Parser;
use packed_seq::{PackedSeq, Seq, SeqVec};
use seq_hash::{KmerHasher, NtHasher};
use simd_minimizers::{
    canonical_closed_syncmers, canonical_minimizers, canonical_open_syncmers,
};

#[derive(Parser, Debug)]
#[command(name = "marker_audit")]
struct Args {
    /// 2-bit packed hg38 (or slice).
    bin_path: PathBuf,
    /// Only audit this many bases from the start (0 = entire file).
    #[arg(long, default_value_t = 50_000_000)]
    max_bases: usize,
    #[arg(long, default_value_t = 31)]
    sync_k: usize,
    #[arg(long, default_value_t = 11)]
    s: usize,
    #[arg(long, default_value_t = 21)]
    min_w: usize,
    /// Optional u32-le position dump from XSeal (same slice / coordinates).
    #[arg(long)]
    xseal_pos: Option<PathBuf>,
    #[arg(long, default_value = "ClosedSync")]
    xseal_scheme: String,
}

fn load_packed(path: &PathBuf, max_bases: usize) -> (Vec<u8>, usize) {
    let mut f = File::open(path).expect("open bin");
    let mut raw = Vec::new();
    f.read_to_end(&mut raw).expect("read bin");
    let file_bases = raw.len() * 4;
    let num_bases = if max_bases == 0 {
        file_bases
    } else {
        max_bases.min(file_bases)
    };
    let bytes = (num_bases + 3) / 4;
    raw.truncate(bytes);
    (raw, num_bases)
}

fn read_pos_dump(path: &PathBuf) -> Vec<u32> {
    let mut f = File::open(path).expect("open pos dump");
    let mut buf = Vec::new();
    f.read_to_end(&mut buf).expect("read pos dump");
    assert_eq!(buf.len() % 4, 0, "pos dump must be u32-le");
    buf.chunks_exact(4)
        .map(|c| u32::from_le_bytes([c[0], c[1], c[2], c[3]]))
        .collect()
}

fn all_smer_hashes(seq: PackedSeq<'_>, s: usize) -> Vec<u32> {
    let hasher = NtHasher::<true>::new(s);
    hasher.hash_kmers_scalar(seq).collect()
}

fn window_smer_hashes(all: &[u32], start: usize, w: usize) -> &[u32] {
    &all[start..start + w]
}

/// simd-minimizers open: canonical minimizer index must be exactly w/2.
fn open_strict_at(all: &[u32], start: usize, w: usize) -> bool {
    let hs = window_smer_hashes(all, start, w);
    let min_h = hs.iter().copied().min().unwrap();
    let mut min_idx = None;
    for (i, &h) in hs.iter().enumerate() {
        if h == min_h {
            if min_idx.is_some() {
                // tie: strict simd uses strand-aware tie break — we approximate by leftmost min hash index
                // only for audit triage; full check uses library output set.
            }
            if min_idx.is_none() {
                min_idx = Some(i);
            }
        }
    }
    min_idx == Some(w / 2)
}

/// XSeal-style open: middle s-mer hash equals window minimum (no tie-break index).
fn open_relaxed_at(all: &[u32], start: usize, w: usize) -> bool {
    let hs = window_smer_hashes(all, start, w);
    let min_h = hs.iter().copied().min().unwrap();
    hs[w / 2] == min_h
}

fn closed_at(all: &[u32], start: usize, w: usize) -> bool {
    let hs = window_smer_hashes(all, start, w);
    let min_h = hs.iter().copied().min().unwrap();
    hs[0] == min_h || hs[w - 1] == min_h
}

fn audit_positions(
    label: &str,
    seq_len: usize,
    all_hashes: &[u32],
    positions: &[u32],
    w: usize,
    k_sync: usize,
) {
    let max_start = seq_len.saturating_sub(k_sync - 1);
    let mut bad_closed = 0usize;
    let mut bad_open_strict = 0usize;
    let mut bad_open_relaxed = 0usize;
    let mut oob = 0usize;
    let mut gaps = Vec::new();

    let mut prev: Option<u32> = None;
    for &p in positions {
        if p as usize >= max_start {
            oob += 1;
            continue;
        }
        if let Some(pr) = prev {
            gaps.push(p.saturating_sub(pr));
        }
        prev = Some(p);

        if !closed_at(all_hashes, p as usize, w) {
            bad_closed += 1;
        }
        if !open_strict_at(all_hashes, p as usize, w) {
            bad_open_strict += 1;
        }
        if !open_relaxed_at(all_hashes, p as usize, w) {
            bad_open_relaxed += 1;
        }
    }

    let n = positions.len();
    let coverage = n as f64 / seq_len as f64;
    let adjacent_lt_k = if gaps.is_empty() {
        0.0
    } else {
        gaps.iter().filter(|&&g| (g as usize) < k_sync).count() as f64 / gaps.len() as f64
    };

    println!(
        "AUDIT|{}|hits|{}|coverage|{:.8}|oob|{}|fail_closed|{}|fail_open_strict|{}|fail_open_relaxed|{}|frac_gap_lt_K|{:.6}",
        label,
        n,
        coverage,
        oob,
        bad_closed,
        bad_open_strict,
        bad_open_relaxed,
        adjacent_lt_k
    );
}

fn compare_sets(a: &[u32], b: &[u32]) -> (usize, usize, usize, f64) {
    use std::collections::HashSet;
    let sa: HashSet<u32> = a.iter().copied().collect();
    let sb: HashSet<u32> = b.iter().copied().collect();
    let inter = sa.intersection(&sb).count();
    let uni = sa.union(&sb).count();
    let jaccard = if uni == 0 { 1.0 } else { inter as f64 / uni as f64 };
    (inter, sa.len(), sb.len(), jaccard)
}

fn main() {
    let args = Args::parse();
    let w = args.sync_k - args.s + 1;
    let lib_s = args.s;
    let k_sync = args.sync_k;

    let (raw, num_bases) = load_packed(&args.bin_path, args.max_bases);
    let safe_len = num_bases.saturating_sub(256);
    let seq = PackedSeq::from_raw_parts(&raw, 0, safe_len);
    let all_hashes = all_smer_hashes(seq, lib_s);

    println!(
        "NOTE|marker_audit|bases|{}|sync_k|{}|s|{}|w|{}|min_k|{}|min_w|{}",
        num_bases, k_sync, lib_s, w, args.sync_k, args.min_w
    );

    // --- simd-minimizers: scalar vs SIMD must match (how the library proves itself) ---
    for (name, scalar, simd) in [
        (
            "ClosedSync",
            canonical_closed_syncmers(lib_s, w).run_scalar_once(seq),
            canonical_closed_syncmers(lib_s, w).run_once(seq),
        ),
        (
            "OpenSync",
            canonical_open_syncmers(lib_s, w).run_scalar_once(seq),
            canonical_open_syncmers(lib_s, w).run_once(seq),
        ),
        (
            "Minimizer",
            canonical_minimizers(args.sync_k, args.min_w).run_scalar_once(seq),
            canonical_minimizers(args.sync_k, args.min_w).run_once(seq),
        ),
    ] {
        assert_eq!(
            scalar, simd,
            "{name}: scalar vs SIMD mismatch (len {} vs {})",
            scalar.len(),
            simd.len()
        );
        println!("OK|SimdMin|{name}|scalar_eq_simd|hits|{}", scalar.len());
        if name == "ClosedSync" {
            audit_positions(
                "SimdMin_ClosedSync",
                safe_len,
                &all_hashes,
                &scalar,
                w,
                k_sync,
            );
        }
        if name == "OpenSync" {
            audit_positions(
                "SimdMin_OpenSync",
                safe_len,
                &all_hashes,
                &scalar,
                w,
                k_sync,
            );
            let relaxed_only = scalar
                .iter()
                .filter(|&&p| {
                    open_relaxed_at(&all_hashes, p as usize, w)
                        && !open_strict_at(&all_hashes, p as usize, w)
                })
                .count();
            println!(
                "AUDIT|SimdMin_OpenSync|relaxed_not_strict|{}|frac|{:.6}",
                relaxed_only,
                relaxed_only as f64 / scalar.len().max(1) as f64
            );
        }
    }

    // RC symmetry (canonical minimizers): same multiset on revcomp.
    let rc = seq.to_revcomp();
    let mut fwd = Vec::new();
    let mut rc_pos = Vec::new();
    canonical_minimizers(args.sync_k, args.min_w).run(seq, &mut fwd);
    canonical_minimizers(args.sync_k, args.min_w).run(rc.as_slice(), &mut rc_pos);
    let len = seq.len();
    let mut rc_ok = 0usize;
    for (&p, &q) in fwd.iter().zip(rc_pos.iter().rev()) {
        if (p as usize) + (q as usize) == len - args.sync_k {
            rc_ok += 1;
        }
    }
    println!(
        "AUDIT|SimdMin|Minimizer|rc_symmetric_pairs|{}/{}",
        rc_ok,
        fwd.len().min(rc_pos.len())
    );

    if let Some(path) = args.xseal_pos {
        let xpos = read_pos_dump(&path);
        let simd_closed = canonical_closed_syncmers(lib_s, w).run_once(seq);
        let (inter, na, nb, j) = compare_sets(&xpos, &simd_closed);
        println!(
            "COMPARE|Xseal_vs_Simd|{}|inter|{}|xseal|{}|simd|{}|jaccard|{:.6}",
            args.xseal_scheme, inter, na, nb, j
        );
        audit_positions(
            &format!("Xseal_{}", args.xseal_scheme),
            safe_len,
            &all_hashes,
            &xpos,
            w,
            k_sync,
        );
        if args.xseal_scheme == "OpenSync" {
            let relaxed_only = xpos
                .iter()
                .filter(|&&p| {
                    open_relaxed_at(&all_hashes, p as usize, w)
                        && !open_strict_at(&all_hashes, p as usize, w)
                })
                .count();
            println!(
                "AUDIT|Xseal_OpenSync|relaxed_not_strict|{}|frac|{:.6}",
                relaxed_only,
                relaxed_only as f64 / xpos.len().max(1) as f64
            );
        }
    }
}
