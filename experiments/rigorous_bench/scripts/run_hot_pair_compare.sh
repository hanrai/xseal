#!/usr/bin/env bash

# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2026 hanrai. All Rights Reserved.

# Build + run simd-minimizers and XSeal hot benches (3 schemes each), 10 outer runs.
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/../../.." && pwd)"
BIN_FILE="${BIN_FILE:-$ROOT/data/hg38_2bit.bin}"
THREADS="${THREADS:-12}"
SYNC_K="${SYNC_K:-31}"
S="${S:-11}"
MIN_W="${MIN_W:-21}"
OUTER_RUNS="${OUTER_RUNS:-10}"

SIMD_OUT="${SIMD_OUT:-$ROOT/experiments/rigorous_bench/simd_min_hot_out}"
XSEAL_OUT="${XSEAL_OUT:-$ROOT/experiments/rigorous_bench/xseal_hot_out}"
COMPARE_OUT="${COMPARE_OUT:-$ROOT/experiments/rigorous_bench/hot_compare_out}"

if [[ ! -f "$BIN_FILE" ]]; then
  echo "Error: missing $BIN_FILE" >&2
  exit 1
fi

mkdir -p "$SIMD_OUT" "$XSEAL_OUT" "$COMPARE_OUT"

# shellcheck source=build_hot_benches.sh
source "$(dirname "$0")/build_hot_benches.sh"

run_tool() {
  local name=$1 bin=$2 out_dir=$3
  local summary="$out_dir/summary_t${THREADS}.txt"
  : >"$summary"
  echo "NOTE|outer|tool|$name|runs|$OUTER_RUNS|threads|$THREADS" | tee -a "$summary"
  for run in $(seq 1 "$OUTER_RUNS"); do
    local log="$out_dir/run_$(printf '%02d' "$run")_t${THREADS}.txt"
    "$bin" "$BIN_FILE" "$THREADS" "$SYNC_K" "$S" "$MIN_W" | tee "$log"
    echo "OUTER_RUN|$run|tool|$name|log|$log" >>"$summary"
  done
}

echo "=== simd-minimizers ($OUTER_RUNS outer) ===" >&2
run_tool SimdMin "$SIMD_BIN" "$SIMD_OUT"

echo "=== XSeal ($OUTER_RUNS outer) ===" >&2
run_tool Xseal "$XSEAL_BIN" "$XSEAL_OUT"

REPORT="$COMPARE_OUT/report_t${THREADS}.md"
python3 "$ROOT/experiments/rigorous_bench/scripts/analyze_hot_compare.py" \
  --simd-glob "$SIMD_OUT/run_*_t${THREADS}.txt" \
  --xseal-glob "$XSEAL_OUT/run_*_t${THREADS}.txt" \
  -o "$REPORT"

echo "Wrote $REPORT"
