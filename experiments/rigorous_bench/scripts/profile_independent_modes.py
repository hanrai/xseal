#!/usr/bin/env python3

# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2026 hanrai. All Rights Reserved.

import subprocess
import re
import os
import sys

BASES = 3209286112
LOOPS = 13
TOTAL_BP = BASES * LOOPS

commands = {
    ("XSeal", "Minimizer"): ["experiments/rigorous_bench/xseal_hot_bench", "data/hg38_2bit.bin", "1", "31", "11", "21", "3"],
    ("XSeal", "ClosedSync"): ["experiments/rigorous_bench/xseal_hot_bench", "data/hg38_2bit.bin", "1", "31", "11", "21", "1"],
    ("XSeal", "OpenSync"): ["experiments/rigorous_bench/xseal_hot_bench", "data/hg38_2bit.bin", "1", "31", "11", "21", "2"],
    ("simd-minimizers", "Minimizer"): ["experiments/rigorous_bench/simd_min_hot_bench/target/release/simd_min_hot_bench", "data/hg38_2bit.bin", "1", "31", "11", "21", "3"],
    ("simd-minimizers", "ClosedSync"): ["experiments/rigorous_bench/simd_min_hot_bench/target/release/simd_min_hot_bench", "data/hg38_2bit.bin", "1", "31", "11", "21", "1"],
    ("simd-minimizers", "OpenSync"): ["experiments/rigorous_bench/simd_min_hot_bench/target/release/simd_min_hot_bench", "data/hg38_2bit.bin", "1", "31", "11", "21", "2"]
}

results = {}

for (tool, mode), cmd in commands.items():
    print(f"Running profile for {tool} - {mode}...", flush=True)
    perf_cmd = ["perf", "stat", "-e", "instructions,cycles,branches,branch-misses"] + cmd
    
    res = subprocess.run(perf_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    
    # Parse perf output from stderr
    stderr = res.stderr
    
    inst = 0
    cyc = 0
    branches = 0
    misses = 0
    time_s = 0.0
    
    for line in stderr.splitlines():
        line = line.strip()
        if "instructions" in line:
            inst = int(re.sub(r"[^\d]", "", line.split("instructions")[0]))
        elif "cycles" in line:
            cyc = int(re.sub(r"[^\d]", "", line.split("cycles")[0]))
        elif "branch-misses" in line:
            misses = int(re.sub(r"[^\d]", "", line.split("branch-misses")[0]))
        elif "branches" in line:
            branches = int(re.sub(r"[^\d]", "", line.split("branches")[0]))
        elif "seconds time elapsed" in line:
            time_s = float(re.search(r"([\d.]+)\s+seconds time elapsed", line).group(1))
            
    ipc = inst / cyc if cyc > 0 else 0.0
    insn_bp = inst / TOTAL_BP
    b_miss_pct = (misses / branches * 100.0) if branches > 0 else 0.0
    b_miss_mbp = (misses / TOTAL_BP * 1000000.0)
    
    results[(tool, mode)] = {
        "instructions": inst,
        "cycles": cyc,
        "branches": branches,
        "branch_misses": misses,
        "time_s": time_s,
        "ipc": ipc,
        "insn_bp": insn_bp,
        "b_miss_pct": b_miss_pct,
        "b_miss_mbp": b_miss_mbp
    }
    
    print(f"  Instructions: {inst:,}")
    print(f"  Cycles: {cyc:,}")
    print(f"  IPC: {ipc:.2f}")
    print(f"  Insn/bp: {insn_bp:.2f}")
    print(f"  Branch Misses: {misses:,} ({b_miss_pct:.2f}%)")
    print(f"  B-Miss / Mbp: {b_miss_mbp:.2f}")
    print(f"  Time: {time_s:.2f}s")
    print()

# Print markdown table
print("\n" + "="*80)
print("INDEPENDENT PMU PROFILING RESULTS FOR TABLE 1")
print("="*80)
print(f"{'Tool':<18} | {'Mode':<12} | {'insn/bp':<8} | {'IPC':<6} | {'B-Miss %':<8} | {'B-Miss/Mbp':<10} | {'Time (s)':<8}")
print("-" * 80)
for (tool, mode), r in results.items():
    print(f"{tool:<18} | {mode:<12} | {r['insn_bp']:8.2f} | {r['ipc']:6.2f} | {r['b_miss_pct']:7.2f}% | {r['b_miss_mbp']:10.2f} | {r['time_s']:8.2f}")
print("="*80)
