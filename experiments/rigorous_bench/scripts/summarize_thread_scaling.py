#!/usr/bin/env python3

# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2026 hanrai. All Rights Reserved.

"""Summarize thread-scaling logs (XsealHot / SimdMinHot) + optional 1T baseline."""

from __future__ import annotations

import argparse
import re
import statistics as st
from pathlib import Path

SCHEMES = ("ClosedSync", "OpenSync", "Minimizer")
SCHEME_SHORT = {"ClosedSync": "Closed", "OpenSync": "Open", "Minimizer": "Min"}


def parse_inner_lines(text: str, tool_prefix: str) -> dict[tuple[str, int], list[float]]:
    """(scheme, threads) -> list of inner means from each outer run."""
    pat = re.compile(
        rf"^INNER\|{re.escape(tool_prefix)}\|(ClosedSync|OpenSync|Minimizer)\|(\d+)\|"
        r"([\d.]+)\|"
    )
    out: dict[tuple[str, int], list[float]] = {}
    for line in text.splitlines():
        m = pat.match(line.strip())
        if not m:
            continue
        scheme, threads, gbp = m.group(1), int(m.group(2)), float(m.group(3))
        out.setdefault((scheme, threads), []).append(gbp)
    return out


def load_1t_baseline(
    baseline_dir: Path, tool_x: str, tool_s: str, hash_mode: str
) -> dict[tuple[str, str], float]:
    """From table1_vs_paper.txt or hash summary."""
    vals: dict[tuple[str, str], float] = {}

    paper_path = baseline_dir / "table1_vs_paper.txt"
    if paper_path.is_file():
        for line in paper_path.read_text(errors="replace").splitlines():
            m = re.match(
                r"^(SimdMinHot|XsealHot|SimdMinHotHash|XsealHotHash)\s+"
                r"(ClosedSync|OpenSync|Minimizer)\s+([\d.]+)",
                line.strip(),
            )
            if m:
                vals[(m.group(1), m.group(2))] = float(m.group(3))
        if vals:
            return vals

    if hash_mode == "1":
        summ = baseline_dir / "table1_hash_summary.txt"
        if summ.is_file():
            pat = re.compile(
                rf"^(SimdMinHotHash|XsealHotHash)\s+(ClosedSync|OpenSync|Minimizer)\s+"
                r"([\d.]+)\s+([\d.]+)"
            )
            for line in summ.read_text(errors="replace").splitlines():
                m = pat.match(line.strip())
                if m:
                    vals[(m.group(1), m.group(2))] = float(m.group(3))
            if vals:
                return vals

    return vals


def mean_std(vals: list[float]) -> tuple[float, float]:
    if not vals:
        return 0.0, 0.0
    if len(vals) == 1:
        return vals[0], 0.0
    return st.mean(vals), st.stdev(vals)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out-dir", required=True, help="repro_scaling_out/pos_only")
    ap.add_argument("--baseline-1t-dir", default=None, help="repro_tab1_tab3_out for 1T")
    ap.add_argument("--hash-mode", default="0")
    ap.add_argument("--tool-x", default=None, help="INNER| tool prefix for XSeal")
    ap.add_argument("--tool-s", default=None, help="INNER| tool prefix for simd")
    ap.add_argument("--title", default=None, help="Report heading")
    ap.add_argument(
        "--no-paper-row",
        action="store_true",
        help="Omit fixed paper 12T comparison row",
    )
    ap.add_argument("-o", "--markdown", default=None)
    args = ap.parse_args()

    out_dir = Path(args.out_dir)
    if args.tool_x:
        tool_x = args.tool_x
    else:
        tool_x = "XsealHotHash" if args.hash_mode == "1" else "XsealHot"
    if args.tool_s:
        tool_s = args.tool_s
    else:
        tool_s = "SimdMinHotHash" if args.hash_mode == "1" else "SimdMinHot"

    all_x: dict[tuple[str, int], list[float]] = {}
    all_s: dict[tuple[str, int], list[float]] = {}
    thread_set: set[int] = set()

    for log in sorted(out_dir.glob("t*/scaling_*t.log")):
        text = log.read_text(errors="replace")
        for key, vals in parse_inner_lines(text, tool_x).items():
            all_x.setdefault(key, []).extend(vals)
            thread_set.add(key[1])
        for key, vals in parse_inner_lines(text, tool_s).items():
            all_s.setdefault(key, []).extend(vals)
            thread_set.add(key[1])

    threads_sorted = sorted(thread_set)
    baseline_1t: dict[tuple[str, str], float] = {}
    if args.baseline_1t_dir:
        baseline_1t = load_1t_baseline(
            Path(args.baseline_1t_dir), tool_x, tool_s, args.hash_mode
        )

    lines: list[str] = []
    mode = "pos+hash" if args.hash_mode == "1" else "pos-only"
    title = args.title or f"Thread scaling summary ({mode})"
    lines.append(f"# {title}\n\n")
    if args.baseline_1t_dir:
        lines.append(
            f"Threads tested: {', '.join(map(str, threads_sorted))} "
            f"(1T from baseline dir)\n\n"
        )
    else:
        lines.append(f"Threads tested: {', '.join(map(str, threads_sorted))}\n\n")

    # Table per scheme
    for scheme in SCHEMES:
        short = SCHEME_SHORT[scheme]
        lines.append(f"## {scheme}\n\n")
        lines.append(
            "| Threads | simd Gbp/s (μ±σ) | XSeal Gbp/s (μ±σ) | XSeal/simd | "
            "XSeal eff vs 1T | simd eff vs 1T |\n"
        )
        lines.append("|---:|---:|---:|---:|---:|---:|\n")

        bx = baseline_1t.get((tool_x, scheme), 0.0)
        bs = baseline_1t.get((tool_s, scheme), 0.0)

        # Fallback to measured 1T values if not loaded from baseline
        if bx == 0.0 and 1 in threads_sorted:
            vx1 = all_x.get((scheme, 1), [])
            if vx1:
                bx = st.mean(vx1)
        if bs == 0.0 and 1 in threads_sorted:
            vs1 = all_s.get((scheme, 1), [])
            if vs1:
                bs = st.mean(vs1)

        if bx > 0:
            lines.append(
                f"| **1** (baseline) | {bs:.4f}* | {bx:.4f}* | "
                f"{bx/bs:.3f}× | 1.00 | 1.00 |\n"
            )

        for t in threads_sorted:
            if t == 1 and bx > 0:
                # Skip duplicate 1T row if printed as baseline above
                continue
            vx = all_x.get((scheme, t), [])
            vs = all_s.get((scheme, t), [])
            if not vx and not vs:
                continue
            mx, sx = mean_std(vx)
            ms, ss = mean_std(vs)
            ratio = mx / ms if ms else 0.0
            eff_x = mx / (bx * t) if bx and t else 0.0
            eff_s = ms / (bs * t) if bs and t else 0.0
            lines.append(
                f"| {t} | {ms:.4f}±{ss:.4f} | {mx:.4f}±{sx:.4f} | {ratio:.3f}× | "
                f"{eff_x:.3f} | {eff_s:.3f} |\n"
            )
        lines.append("\n*1T from prior `repro_tab1_tab3_out` run or direct 1T scaling measurement.\n\n")

    # Paper tab:scaling style (Min / Closed / Open rows as columns)
    lines.append("## Paper-style aggregate (Gbp/s μ)\n\n")
    lines.append("| Threads | simd Min | simd Closed | simd Open | XSeal Min | XSeal Closed | XSeal Open |\n")
    lines.append("|---:|---:|---:|---:|---:|---:|---:|\n")

    def cell(tool_data: dict, tool: str, scheme: str, t: int) -> str:
        if t == 1:
            v = baseline_1t.get((tool, scheme))
            if v:
                return f"{v:.4f}"
            vals = tool_data.get((scheme, t), [])
            return f"{st.mean(vals):.4f}" if vals else "—"
        vals = tool_data.get((scheme, t), [])
        return f"{st.mean(vals):.4f}" if vals else "—"

    if baseline_1t or 1 in threads_sorted:
        lines.append(
            f"| 1 | {cell(all_s, tool_s, 'Minimizer', 1)} | "
            f"{cell(all_s, tool_s, 'ClosedSync', 1)} | {cell(all_s, tool_s, 'OpenSync', 1)} | "
            f"{cell(all_x, tool_x, 'Minimizer', 1)} | "
            f"{cell(all_x, tool_x, 'ClosedSync', 1)} | "
            f"{cell(all_x, tool_x, 'OpenSync', 1)} |\n"
        )

    for t in threads_sorted:
        if t == 1:
            continue
        lines.append(
            f"| {t} | {cell(all_s, tool_s, 'Minimizer', t)} | "
            f"{cell(all_s, tool_s, 'ClosedSync', t)} | {cell(all_s, tool_s, 'OpenSync', t)} | "
            f"{cell(all_x, tool_x, 'Minimizer', t)} | "
            f"{cell(all_x, tool_x, 'ClosedSync', t)} | "
            f"{cell(all_x, tool_x, 'OpenSync', t)} |\n"
        )

    if not args.no_paper_row:
        lines.append("\n| 12 (paper) | 2.58 | 2.79 | 3.25 | 7.30 | 6.69 | 8.41 |\n")

    lines.append("\n## Notes\n\n")
    lines.append("- **Efficiency** = throughput / (threads × 1T throughput); 1.0 = perfect linear scaling.\n")
    lines.append("- Ryzen 5600X: 6C/12T; threads > 12 oversubscribe.\n")
    lines.append("- Compare paper 12T row to local t=12 for validation.\n")

    text = "".join(lines)
    md_path = Path(args.markdown) if args.markdown else out_dir / "SCALING_REPORT.md"
    md_path.write_text(text, encoding="utf-8")
    print(text)


if __name__ == "__main__":
    main()
