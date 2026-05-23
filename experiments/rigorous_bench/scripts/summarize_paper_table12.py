#!/usr/bin/env python3

# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2026 hanrai. All Rights Reserved.

"""Summarize repro_paper_table12_out -> Table 1 / Table 2 with μ±σ (10 outer runs)."""

from __future__ import annotations

import os
import re
import statistics as st
from pathlib import Path

SCHEMES = ("Minimizer", "ClosedSync", "OpenSync")
SCHEME_LABEL = {
    "Minimizer": "Minimizer",
    "ClosedSync": "Sync. (Closed)",
    "OpenSync": "Sync. (Open)",
}

PAPER_TAB1 = {
    ("SimdMinHot", "Minimizer"): 0.59,
    ("SimdMinHot", "ClosedSync"): 0.56,
    ("SimdMinHot", "OpenSync"): 0.71,
    ("XsealHot", "Minimizer"): 1.17,
    ("XsealHot", "ClosedSync"): 1.08,
    ("XsealHot", "OpenSync"): 1.36,
}

PAPER_TAB2 = {
    1: {
        "SimdMinHot": {"Minimizer": 0.59, "ClosedSync": 0.56, "OpenSync": 0.71},
        "XsealHot": {"Minimizer": 1.17, "ClosedSync": 1.08, "OpenSync": 1.36},
    },
    6: {
        "SimdMinHot": {"Minimizer": 2.65, "ClosedSync": 2.45, "OpenSync": 3.12},
        "XsealHot": {"Minimizer": 7.30, "ClosedSync": 6.69, "OpenSync": 7.85},
    },
    12: {
        "SimdMinHot": {"Minimizer": 2.58, "ClosedSync": 2.79, "OpenSync": 3.25},
        "XsealHot": {"Minimizer": 7.30, "ClosedSync": 6.69, "OpenSync": 8.41},
    },
}


def parse_log(path: Path, tool: str) -> dict[tuple[str, int], list[float]]:
    pat = re.compile(
        rf"^INNER\|{re.escape(tool)}\|(ClosedSync|OpenSync|Minimizer)\|(\d+)\|"
        r"([\d.]+)\|"
    )
    out: dict[tuple[str, int], list[float]] = {}
    if not path.is_file():
        return out
    for line in path.read_text(errors="replace").splitlines():
        m = pat.match(line.strip())
        if not m:
            continue
        scheme, threads, gbp = m.group(1), int(m.group(2)), float(m.group(3))
        out.setdefault((scheme, threads), []).append(gbp)
    return out


def mu_sigma(vals: list[float]) -> tuple[float, float]:
    if not vals:
        return float("nan"), float("nan")
    if len(vals) == 1:
        return vals[0], 0.0
    return st.mean(vals), st.stdev(vals)


def fmt_mu_sigma(mu: float, sig: float) -> str:
    if mu != mu:
        return "N/A"
    return f"{mu:.2f}±{sig:.2f}"


def pct_delta(measured: float, paper: float) -> str:
    if paper == 0 or measured != measured:
        return "N/A"
    return f"{100.0 * (measured - paper) / paper:+.1f}%"


def main() -> None:
    out_dir = Path(os.environ.get("OUT_DIR", "repro_paper_table12_out"))
    outer_runs = int(os.environ.get("OUTER_RUNS", "10"))

    tab1_s = parse_log(out_dir / "table1_t1.log", "SimdMinHot")
    tab1_x = parse_log(out_dir / "table1_t1.log", "XsealHot")

    lines: list[str] = []
    lines.append("=" * 72)
    lines.append("Paper Table 1 — Single-thread (hg38 hot RAM, pos-only)")
    lines.append(f"Protocol: 3 warmup + 10 recorded loops × {outer_runs} outer runs")
    lines.append("=" * 72)
    lines.append(
        f"{'Implementation':<22} {'Mode':<16} {'Measured':>12} {'Paper':>10} {'Δ%':>8}"
    )
    lines.append("-" * 72)

    rows_tab1 = [
        ("SimdMinHot", "simd-minimizers"),
        ("XsealHot", "XSeal (genome-part.)"),
    ]
    for tool, label in rows_tab1:
        tab = tab1_s if tool == "SimdMinHot" else tab1_x
        for scheme in SCHEMES:
            vals = tab.get((scheme, 1), [])
            mu, sig = mu_sigma(vals)
            paper = PAPER_TAB1.get((tool, scheme), float("nan"))
            lines.append(
                f"{label:<22} {SCHEME_LABEL[scheme]:<16} {fmt_mu_sigma(mu, sig):>12} "
                f"{paper:>10.2f} {pct_delta(mu, paper):>8}"
            )
    lines.append("")

    # block128k @ 1T (supplementary)
    blk1 = parse_log(out_dir / "t1/block128k.log", "XsealBlock128kCnt")
    lines.append("Supplementary 1T — XSeal block128k dispatch (pos-only):")
    for scheme in SCHEMES:
        vals = blk1.get((scheme, 1), [])
        mu, sig = mu_sigma(vals)
        paper = PAPER_TAB1.get(("XsealHot", scheme), float("nan"))
        lines.append(
            f"  {SCHEME_LABEL[scheme]:<16} {fmt_mu_sigma(mu, sig):>12}  "
            f"(vs paper XSeal {paper:.2f}: {pct_delta(mu, paper)})"
        )
    lines.append("")

    lines.append("=" * 72)
    lines.append("Paper Table 2 — Multi-thread aggregate Gbp/s (pos-only)")
    lines.append("=" * 72)

    for t in (1, 6, 12):
        gp_log = out_dir / f"t{t}/genome_partition.log"
        blk_log = out_dir / f"t{t}/block128k.log"
        gp_s = parse_log(gp_log, "SimdMinHot")
        gp_x = parse_log(gp_log, "XsealHot")
        blk_x = parse_log(blk_log, "XsealBlock128kCnt")
        blk_s = parse_log(blk_log, "SimdMinBlock128kCnt")

        lines.append(f"\n--- {t} thread(s) ---")
        lines.append(
            f"{'Backend':<28} {'Minimizer':>14} {'Closed':>14} {'Open':>14}"
        )
        lines.append("-" * 72)

        def row(label: str, data: dict, tool: str) -> None:
            cells = []
            for scheme in SCHEMES:
                mu, sig = mu_sigma(data.get((scheme, t), []))
                cells.append(fmt_mu_sigma(mu, sig))
            lines.append(f"{label:<28} {cells[0]:>14} {cells[1]:>14} {cells[2]:>14}")

        row("simd-min (genome-part.)", gp_s, "SimdMinHot")
        row("XSeal (genome-part.)", gp_x, "XsealHot")
        row("simd-min (block128k)", blk_s, "SimdMinBlock128kCnt")
        row("XSeal (block128k)", blk_x, "XsealBlock128kCnt")

        if t in PAPER_TAB2:
            p = PAPER_TAB2[t]
            lines.append(
                f"{'paper simd-min':<28} {p['SimdMinHot']['Minimizer']:>14.2f} "
                f"{p['SimdMinHot']['ClosedSync']:>14.2f} "
                f"{p['SimdMinHot']['OpenSync']:>14.2f}"
            )
            lines.append(
                f"{'paper XSeal':<28} {p['XsealHot']['Minimizer']:>14.2f} "
                f"{p['XsealHot']['ClosedSync']:>14.2f} "
                f"{p['XsealHot']['OpenSync']:>14.2f}"
            )

    lines.append("\n" + "=" * 72)
    lines.append("Table 2 @ 12T — XSeal genome-part. vs block128k vs paper")
    lines.append("=" * 72)
    t = 12
    gp_x = parse_log(out_dir / f"t{t}/genome_partition.log", "XsealHot")
    blk_x = parse_log(out_dir / f"t{t}/block128k.log", "XsealBlock128kCnt")
    for scheme in SCHEMES:
        mu_gp, sig_gp = mu_sigma(gp_x.get((scheme, t), []))
        mu_blk, sig_blk = mu_sigma(blk_x.get((scheme, t), []))
        paper = PAPER_TAB2[12]["XsealHot"][scheme]
        lines.append(f"\n{SCHEME_LABEL[scheme]}:")
        lines.append(
            f"  genome-part.  {fmt_mu_sigma(mu_gp, sig_gp)}  (paper {paper:.2f}, "
            f"{pct_delta(mu_gp, paper)})"
        )
        lines.append(
            f"  block128k     {fmt_mu_sigma(mu_blk, sig_blk)}  (vs paper "
            f"{pct_delta(mu_blk, paper)})"
        )
        if mu_gp == mu_gp and mu_blk == mu_blk and mu_gp > 0:
            lines.append(
                f"  128k vs coarse +{100.0 * (mu_blk - mu_gp) / mu_gp:.1f}%"
            )

    text = "\n".join(lines) + "\n"
    (out_dir / "PAPER_TABLE12_SUMMARY.txt").write_text(text, encoding="utf-8")
    print(text)


if __name__ == "__main__":
    main()
