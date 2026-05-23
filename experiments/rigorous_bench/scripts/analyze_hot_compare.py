#!/usr/bin/env python3

# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2026 hanrai. All Rights Reserved.

"""Aggregate hot-bench logs and compare XSeal vs simd-minimizers."""

from __future__ import annotations

import argparse
import glob
import math
import os
import re
from collections import defaultdict
from dataclasses import dataclass
from typing import Dict, List, Tuple

SCHEMES = ("ClosedSync", "OpenSync", "Minimizer")
INNER_RE = re.compile(
    r"^INNER\|(SimdMinHot|XsealHot)\|(ClosedSync|OpenSync|Minimizer)\|(\d+)\|"
    r"([\d.]+)\|([\d.]+)\|(\d+)\|([\d.eE+-]+)\|([0-9a-fA-F]+)"
)


@dataclass
class InnerRow:
    tool: str
    scheme: str
    threads: int
    gbp_mean: float
    gbp_std: float
    hits: int
    coverage: float
    digest: str


def parse_log(path: str) -> List[InnerRow]:
    rows: List[InnerRow] = []
    with open(path, encoding="utf-8", errors="replace") as f:
        for line in f:
            m = INNER_RE.match(line.strip())
            if not m:
                continue
            tool, scheme, threads, gmean, gstd, hits, cov, digest = m.groups()
            rows.append(
                InnerRow(
                    tool=tool,
                    scheme=scheme,
                    threads=int(threads),
                    gbp_mean=float(gmean),
                    gbp_std=float(gstd),
                    hits=int(hits),
                    coverage=float(cov),
                    digest=digest,
                )
            )
    return rows


def mean_std(vals: List[float]) -> Tuple[float, float]:
    n = len(vals)
    if n == 0:
        return 0.0, 0.0
    mean = sum(vals) / n
    var = sum((x - mean) ** 2 for x in vals) / n
    return mean, math.sqrt(var)


def load_runs(glob_pat: str) -> Dict[Tuple[str, str], List[InnerRow]]:
    by_key: Dict[Tuple[str, str], List[List[InnerRow]]] = defaultdict(list)
    for path in sorted(glob.glob(glob_pat)):
        rows = parse_log(path)
        if not rows:
            continue
        key = (rows[0].tool, rows[0].scheme)
        # one file may contain all schemes
        per_scheme: Dict[str, List[InnerRow]] = defaultdict(list)
        for r in rows:
            per_scheme[r.scheme].append(r)
        for scheme, rs in per_scheme.items():
            by_key[(rs[0].tool, scheme)].append(rs)

    flat: Dict[Tuple[str, str], List[InnerRow]] = {}
    for key, runs in by_key.items():
        merged: List[InnerRow] = []
        for run_rows in runs:
            merged.extend(run_rows)
        flat[key] = merged
    return flat


def aggregate(rows: List[InnerRow]) -> dict:
  # rows from multiple outer runs: group by (tool, scheme) — each outer run contributes 3 schemes
    outer: Dict[int, List[InnerRow]] = defaultdict(list)
    # infer outer run id by order: every 3 lines same tool is one process invocation
    # safer: one log file = one outer run with up to 3 schemes
    return {}


def aggregate_from_run_files(glob_pat: str) -> Dict[Tuple[str, str], dict]:
    """Each log file is one outer sandbox run (3 schemes inside)."""
    stats: Dict[Tuple[str, str], dict] = {}
    buckets: Dict[Tuple[str, str], List[InnerRow]] = defaultdict(list)

    for path in sorted(glob.glob(glob_pat)):
        for row in parse_log(path):
            buckets[(row.tool, row.scheme)].append(row)

    for key, rows in buckets.items():
        gbp_means = [r.gbp_mean for r in rows]
        hits_vals = [r.hits for r in rows]
        cov_vals = [r.coverage for r in rows]
        om, os = mean_std(gbp_means)
        hm = hits_vals[-1] if hits_vals else 0
        hits_stable = len(set(hits_vals)) == 1
        stats[key] = {
            "outer_runs": len(rows),
            "outer_gbp_mean": om,
            "outer_gbp_std": os,
            "inner_gbp_mean_avg": sum(r.gbp_mean for r in rows) / len(rows),
            "inner_gbp_std_avg": sum(r.gbp_std for r in rows) / len(rows),
            "hits": hm,
            "hits_stable": hits_stable,
            "coverage": cov_vals[-1] if cov_vals else 0.0,
            "coverage_pct": (cov_vals[-1] if cov_vals else 0.0) * 100.0,
            "digest_last": rows[-1].digest if rows else "",
            "threads": rows[0].threads if rows else 0,
        }
    return stats


def print_report(simd_stats, xseal_stats, out_path: str | None) -> None:
    lines: List[str] = []
    lines.append("# Hot 2-bit bench comparison (K=31, s=11, w=21)\n")
    lines.append(
        "| Scheme | Tool | Outer Gbp/s μ±σ | Inner σ̄ | Hits | Coverage % | vs peer hits | vs peer speed |\n"
    )
    lines.append("|---|---|---:|---:|---:|---:|---:|---:|\n")

    for scheme in SCHEMES:
        sk = ("SimdMinHot", scheme)
        xk = ("XsealHot", scheme)
        s = simd_stats.get(sk)
        x = xseal_stats.get(xk)
        for tool, st in (("SimdMin", s), ("XSeal", x)):
            if not st:
                continue
            hit_note = "stable" if st["hits_stable"] else "VARIES"
            lines.append(
                f"| {scheme} | {tool} | {st['outer_gbp_mean']:.4f}±{st['outer_gbp_std']:.4f} | "
                f"{st['inner_gbp_std_avg']:.4f} | {st['hits']:,} ({hit_note}) | "
                f"{st['coverage_pct']:.4f} | — | — |\n"
            )

    lines.append("\n## Pairwise ratios (XSeal / SimdMin)\n\n")
    lines.append("| Scheme | Speed ratio | Hits ratio | Δ coverage (pp) |\n")
    lines.append("|---|---:|---:|---:|\n")
    for scheme in SCHEMES:
        s = simd_stats.get(("SimdMinHot", scheme))
        x = xseal_stats.get(("XsealHot", scheme))
        if not s or not x:
            continue
        spd = x["outer_gbp_mean"] / s["outer_gbp_mean"] if s["outer_gbp_mean"] else 0
        hr = x["hits"] / s["hits"] if s["hits"] else 0
        dpp = (x["coverage_pct"] - s["coverage_pct"])
        lines.append(
            f"| {scheme} | {spd:.3f}× | {hr:.4f} | {dpp:+.4f} |\n"
        )

    lines.append(
        "\n## Biological / semantic validation (hash-agnostic)\n\n"
        "1. **Definition audit**: For each hit position `p`, verify the scheme predicate on the extracted "
        "substring (closed/open syncmer or window minimizer) using *that tool's* hash — not cross-compare hashes.\n"
        "2. **Reference oracle**: Run a slow scalar reference on a subset (e.g. 1 Mb) and require 100% "
        "predicate match for XSeal and simd-minimizers separately.\n"
        "3. **Set overlap**: On identical partitioning, report Jaccard = |A∩B|/|A∪B| for positions; "
        "low Jaccard with similar coverage implies different tie-break / canonical rules, not necessarily bugs.\n"
        "4. **RC symmetry**: Check canonical syncmers on `(seq, revcomp(seq))` share the same multiset of "
        "canonical k-mers (positions mirrored).\n"
        "5. **Downstream invariance**: Index recall / mapping sensitivity vs a gold aligner — the biology-relevant "
        "metric when hashes differ.\n"
    )

    text = "".join(lines)
    print(text)
    if out_path:
        with open(out_path, "w", encoding="utf-8") as f:
            f.write(text)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--simd-glob", required=True, help="e.g. simd_min_hot_out/run_*_t12.txt")
    ap.add_argument("--xseal-glob", required=True)
    ap.add_argument("-o", "--output", default=None)
    args = ap.parse_args()

    simd_stats = aggregate_from_run_files(args.simd_glob)
    xseal_stats = aggregate_from_run_files(args.xseal_glob)
    print_report(simd_stats, xseal_stats, args.output)


if __name__ == "__main__":
    main()
