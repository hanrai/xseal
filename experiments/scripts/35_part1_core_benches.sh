#!/usr/bin/env bash

# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2026 hanrai. All Rights Reserved.

set -euo pipefail
# shellcheck source=experiments/scripts/common.sh
source "$(cd "$(dirname "$0")" && pwd)/common.sh"
load_config

mkdir -p "$REPO_ROOT/build/bin" "$EXPERIMENTS_ROOT/results"
OUT_JSON="$EXPERIMENTS_ROOT/results/metrics_part1.json"
LOG="$EXPERIMENTS_ROOT/results/part1_core.log"
SKIP="$EXPERIMENTS_ROOT/results/part1_core.skip"
: >"$LOG"

TWO_BIT="${TWO_BIT_HG38:-$REPO_ROOT/build/hg38_2bit.bin}"
K=31
S=11
W=21

if [[ ! -f "${HG38_FASTA:-}" ]] && [[ -f "$REPO_ROOT/data/hg38.fa" ]]; then
  HG38_FASTA="$REPO_ROOT/data/hg38.fa"
fi

if [[ ! -f "$TWO_BIT" ]]; then
  if [[ -f "${HG38_FASTA:-}" ]]; then
    echo "[part1] packing 2-bit from $HG38_FASTA -> $TWO_BIT" | tee -a "$LOG"
    python3 "$REPO_ROOT/experiments/rigorous_bench/scripts/pack_hg38.py" "$HG38_FASTA" "$TWO_BIT" | tee -a "$LOG"
  else
    echo "No hg38 FASTA for 2-bit pack; set HG38_FASTA or TWO_BIT_HG38" | tee "$SKIP"
    echo '{}' >"$OUT_JSON"
    exit 0
  fi
fi

if [[ ! -f "$REPO_ROOT/external/simd-minimizers/Cargo.toml" ]]; then
  echo "simd-minimizers missing; run 05_fetch_vendor.sh" | tee "$SKIP"
  echo '{}' >"$OUT_JSON"
  exit 0
fi

XR="$REPO_ROOT/build/bin/xseal_scan_throughput_bench"
clang++ -O3 -march=native -fopenmp=libgomp -I"$REPO_ROOT/include" \
  "$REPO_ROOT/experiments/rigorous_bench/xseal_scan_throughput_bench.cpp" -o "$XR" 2>>"$LOG" || {
  echo "xseal_scan_throughput_bench compile failed" | tee "$SKIP"
  echo '{}' >"$OUT_JSON"
  exit 0
}

if ! (cd "$REPO_ROOT/experiments/rigorous_bench/simdmin_scan_throughput_bench" && env -u CARGO_TARGET_DIR RUSTFLAGS="-C target-cpu=native" cargo build --release -q --target-dir "$REPO_ROOT/build/rust_target"); then
  echo "simdmin_scan_throughput_bench cargo build failed" | tee "$SKIP"
  echo '{}' >"$OUT_JSON"
  exit 0
fi
cp "$REPO_ROOT/build/rust_target/release/simdmin_scan_throughput_bench" "$REPO_ROOT/build/bin/"
SM="$REPO_ROOT/build/bin/simdmin_scan_throughput_bench"

IFS=',' read -r -a PT <<<"${PART1_THREADS_LIST:-1,12}"
for t in "${PT[@]}"; do
  t="$(echo "$t" | tr -d ' ')"
  [[ -z "$t" ]] && continue
  echo "=== threads=$t XSeal ===" | tee -a "$LOG"
  if ! "$XR" "$TWO_BIT" "$t" "$K" "$S" "$W" 2>&1 | tee -a "$LOG"; then
    echo "part1: xseal_scan failed threads=$t" | tee -a "$LOG"
    exit 1
  fi
  echo "=== threads=$t SimdMin ===" | tee -a "$LOG"
  if ! "$SM" "$TWO_BIT" "$t" "$K" "$S" "$W" 2>&1 | tee -a "$LOG"; then
    echo "part1: simdmin_scan failed threads=$t" | tee -a "$LOG"
    exit 1
  fi
done

PERF_OUT="$EXPERIMENTS_ROOT/results/part1_perf_onethread.txt"
if [[ "${PART1_PERF:-0}" == "1" ]] && command -v perf >/dev/null 2>&1; then
  echo "[part1] perf stat (1 thread)" | tee -a "$LOG"
  perf stat -e instructions,cycles,branches,branch-misses -- "$XR" "$TWO_BIT" 1 "$K" "$S" "$W" \
    >"$PERF_OUT" 2>&1 || echo "perf failed" >>"$PERF_OUT"
else
  echo "[part1] skip perf (set PART1_PERF=1 and install perf)" | tee -a "$LOG"
fi

export LOG OUT_JSON
python3 <<'PY'
import json, re, statistics as st, os

log_path = os.environ["LOG"]
out_path = os.environ["OUT_JSON"]
text = open(log_path, errors="ignore").read().splitlines()


def collect(tool: str, mode: str, threads: int):
    vals = []
    prefix = f"RESULT|{tool}|{mode}|{threads}|"
    for line in text:
        if not line.startswith(prefix):
            continue
        parts = line.split("|")
        if len(parts) < 6:
            continue
        try:
            vals.append(float(parts[4]))
        except ValueError:
            continue
    if not vals:
        return None
    return {
        "mean": st.mean(vals),
        "stdev": st.stdev(vals) if len(vals) > 1 else 0.0,
        "n": len(vals),
        "unit": "Gbp/s",
    }


obj = {}
modes = [("OpenSync", "open"), ("Syncmer", "closed_sync"), ("Minimizer", "minimizer")]
for t in (1, 12):
    for mode, slug in modes:
        for tool, keyp in (("Xseal", "xseal"), ("SimdMin", "simdmin")):
            mid = f"{keyp}_{slug}_{t}t_gbp_s"
            c = collect(tool, mode, t)
            if c:
                obj[mid] = c

open(out_path, "w").write(json.dumps(obj, indent=2))
print(json.dumps(obj, indent=2)[:4000])
PY

echo "Wrote $OUT_JSON"
