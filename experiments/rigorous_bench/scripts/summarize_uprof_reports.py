#!/usr/bin/env python3

# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2026 hanrai. All Rights Reserved.

"""Parse AMD uProf report.csv (assess / custom PMU) and print xseal vs simd-min comparison."""
from __future__ import annotations

import csv
import re
import sys
from pathlib import Path


def parse_report(path: Path) -> dict[str, dict[str, float]]:
    text = path.read_text(errors="replace")
    rows: dict[str, dict[str, float]] = {}

    # Process-level row: first line after '"10 HOTTEST PROCESSES"'
    lines = text.splitlines()
    for i, line in enumerate(lines):
        if line.startswith("PROCESS,") and i > 0 and "HOTTEST PROCESSES" in lines[i - 1]:
            headers = [h.strip().strip('"') for h in line[len("PROCESS,") :].split(",")]
            if i + 1 >= len(lines):
                break
            values = next(csv.reader([lines[i + 1]]))
            if not values:
                continue
            proc = values[0]
            rows[proc] = {}
            for h, v in zip(headers, values[1:], strict=False):
                try:
                    rows[proc][h] = float(v)
                except ValueError:
                    pass
            break
    return rows


def main() -> None:
    root = Path(__file__).resolve().parents[1] / "uprof_assess_out"
    pairs = [
        ("assess", root / "xseal_st/report.csv", root / "simd_st/report.csv"),
        ("avx_pmu", root / "xseal_avx/report.csv", root / "simd_avx/report.csv"),
    ]
    metrics = [
        "RETIRED_INST",
        "CYCLES_NOT_IN_HALT",
        "IPC",
        "CPI",
        "RETIRED_BR_INST_MISP (PTI)",
        "%RETIRED_BR_INST_MISP",
        "L1_DC_ACCESSES (PTI)",
        "L1_DC_MISSES (PTI)",
        "%L1_DC_MISSES",
        "RETIRED_SSE_AVX_FLOPS",
        "RETIRED_MACRO_OPS",
    ]

    for label, xs, sim in pairs:
        if not xs.is_file() or not sim.is_file():
            print(f"[{label}] skip (missing csv)")
            continue
        dx = next(iter(parse_report(xs).values()), {})
        ds = next(iter(parse_report(sim).values()), {})
        if not dx or not ds:
            print(f"[{label}] no process row in csv")
            continue
        print(f"\n=== {label} (PMU sample weights, same hg38 2bit 1T hot bench) ===")
        print(f"{'metric':<32} {'XSeal':>12} {'simd-min':>12} {'simd/xseal':>10}")
        for k in metrics:
            if k not in dx or k not in ds:
                continue
            xv, sv = dx[k], ds[k]
            ratio = sv / xv if xv else float("nan")
            print(f"{k:<32} {xv:12.4f} {sv:12.4f} {ratio:10.3f}")
        if "RETIRED_INST" in dx and "RETIRED_SSE_AVX_FLOPS" in dx:
            fx = dx["RETIRED_SSE_AVX_FLOPS"] / dx["RETIRED_INST"]
            fs = ds["RETIRED_SSE_AVX_FLOPS"] / ds["RETIRED_INST"]
            print(f"{'AVX_FLOPs per 1k INST':<32} {fx*1000:12.2f} {fs*1000:12.2f} {fs/fx if fx else 0:10.3f}")


if __name__ == "__main__":
    main()
