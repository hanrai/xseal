#!/usr/bin/env python3

# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2026 hanrai. All Rights Reserved.

"""Mutation robustness (Figure 7): compare minimizer position conservation under random SNPs.

Usage:
  # Plot paper reference anchors only (no benchmark execution):
  python3 mutation_robustness.py --reference-only

  # Live run on hg38 2-bit (requires built scan throughput benches):
  python3 mutation_robustness.py --run --bin data/hg38_2bit.bin
"""
from __future__ import annotations

import argparse
import os
import random
import subprocess
import sys
from pathlib import Path

# Published Figure 7 anchors (manuscript); not recomputed unless --run succeeds.
PAPER_FIG7_ANCHORS = {
    "xseal": [1.0, 0.8602, 0.4960, 0.2927],
    "simd-min": [1.0, 0.8258, 0.4218, 0.2422],
}
RATES = [0.0, 0.01, 0.05, 0.10]


def repo_root() -> Path:
    return Path(__file__).resolve().parents[2]


def simulate_mutations(input_2bit: Path, output_2bit: Path, rate: float, seed: int) -> None:
    data = bytearray(input_2bit.read_bytes())
    num_bases = len(data) * 4
    num_mutations = int(num_bases * rate)
    rng = random.Random(seed)
    for idx in rng.sample(range(num_bases), num_mutations):
        byte_idx = idx // 4
        bit_offset = (idx % 4) * 2
        data[byte_idx] ^= 1 << bit_offset
    output_2bit.write_bytes(data)


def load_pos(path: Path) -> set[int]:
    if not path.is_file():
        return set()
    import numpy as np

    return set(np.fromfile(path, dtype=np.uint32).tolist())


def run_bench(repo: Path, bin_file: Path, tool: str, out_file: Path, k: int, s: int, w: int) -> set[int]:
    if tool == "xseal":
        exe = repo / "build" / "bin" / "xseal_scan_throughput_bench"
    else:
        exe = repo / "build" / "bin" / "simdmin_scan_throughput_bench"
    if not exe.is_file():
        raise FileNotFoundError(f"missing {exe}; build scan throughput benches first")
    cmd = [str(exe), str(bin_file), "1", str(k), str(s), str(w), str(out_file)]
    subprocess.run(cmd, check=True, capture_output=True)
    return load_pos(out_file)


def conservation(ref: set[int], mut: set[int]) -> float:
    if not ref:
        return 0.0
    return len(ref & mut) / len(ref)


def live_curve(repo: Path, bin_file: Path, k: int, s: int, w: int) -> dict[str, list[float]]:
    tmp = repo / "experiments" / "results" / "mutation_robustness_tmp"
    tmp.mkdir(parents=True, exist_ok=True)
    ref_x = tmp / "ref_xseal.pos"
    ref_s = tmp / "ref_simd.pos"
    ref_x_set = run_bench(repo, bin_file, "xseal", ref_x, k, s, w)
    ref_s_set = run_bench(repo, bin_file, "simd-min", ref_s, k, s, w)
    out: dict[str, list[float]] = {"xseal": [], "simd-min": []}
    for i, rate in enumerate(RATES):
        mut_bin = tmp / f"mut_{i}.2bit"
        simulate_mutations(bin_file, mut_bin, rate, seed=42 + i)
        mx = tmp / f"mut_x_{i}.pos"
        ms = tmp / f"mut_s_{i}.pos"
        mx_set = run_bench(repo, mut_bin, "xseal", mx, k, s, w)
        ms_set = run_bench(repo, mut_bin, "simd-min", ms, k, s, w)
        out["xseal"].append(conservation(ref_x_set, mx_set))
        out["simd-min"].append(conservation(ref_s_set, ms_set))
    return out


def plot_curves(results: dict[str, list[float]], out_pdf: Path, subtitle: str) -> None:
    import matplotlib.pyplot as plt

    plt.figure(figsize=(8, 5))
    plt.plot(RATES, results["xseal"], "o-", label="XSEAL (Independent Hash)", color="#1f77b4", linewidth=2)
    plt.plot(RATES, results["simd-min"], "s--", label="simd-min (Murmur Hash)", color="#ff7f0e", linewidth=2)
    plt.xlabel("Mutation Rate", fontsize=12, fontweight="bold")
    plt.ylabel("Minimizer Conservation", fontsize=12, fontweight="bold")
    plt.title(f"Figure 7: Mutation Robustness of Sampling Hashes\n{subtitle}", fontsize=12, fontweight="bold")
    plt.grid(True, linestyle="--", alpha=0.7)
    plt.legend(fontsize=11)
    out_pdf.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(out_pdf, bbox_inches="tight")
    print(f"Wrote {out_pdf}")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--run", action="store_true", help="Execute live benchmark (slow; needs hg38_2bit.bin)")
    ap.add_argument("--reference-only", action="store_true", help="Plot manuscript anchor curves only")
    ap.add_argument("--bin", type=Path, default=None, help="2-bit packed genome for --run")
    ap.add_argument("-o", "--output", type=Path, default=None)
    args = ap.parse_args()
    repo = repo_root()
    out_pdf = args.output or (repo / "experiments" / "results" / "fig7_robustness.pdf")

    if args.run:
        bin_file = args.bin or (repo / "data" / "hg38_2bit.bin")
        if not bin_file.is_file():
            print(f"error: missing {bin_file}", file=sys.stderr)
            return 1
        print("Running live mutation robustness (K=21 S=11 W=11) ...", file=sys.stderr)
        results = live_curve(repo, bin_file, k=21, s=11, w=11)
        plot_curves(results, out_pdf, subtitle="(measured)")
        return 0

    if args.reference_only:
        plot_curves(PAPER_FIG7_ANCHORS, out_pdf, subtitle="(manuscript reference anchors)")
        return 0

    print("Specify --run (live) or --reference-only (published anchors).", file=sys.stderr)
    print("Hard-coded curves without flags are disabled in the cleanroom release.", file=sys.stderr)
    return 2


if __name__ == "__main__":
    sys.exit(main())
