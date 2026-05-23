#!/usr/bin/env bash

# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2026 hanrai. All Rights Reserved.

# Outer sandbox: 10 process runs. Each run: 3 warmup + 10 timed loops (13 total).
# Aggregates mean/std across the 10 inner means (Gbp/s).
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/../../.." && pwd)"
BIN_FILE="${BIN_FILE:-$ROOT/data/hg38_2bit.bin}"
THREADS="${THREADS:-12}"
SYNC_K="${SYNC_K:-31}"
S="${S:-11}"
MIN_W="${MIN_W:-21}"
OUTER_RUNS="${OUTER_RUNS:-10}"
OUT_DIR="${OUT_DIR:-$ROOT/experiments/rigorous_bench/simd_min_hot_out}"

if [[ ! -f "$BIN_FILE" ]]; then
  echo "Error: 2-bit packed input not found: $BIN_FILE" >&2
  echo "Pack with: python3 $ROOT/experiments/rigorous_bench/scripts/pack_hg38.py <hg38.fa> $BIN_FILE" >&2
  exit 1
fi

mkdir -p "$OUT_DIR"
# shellcheck source=build_hot_benches.sh
source "$(dirname "$0")/build_hot_benches.sh"
BENCH_BIN="$SIMD_BIN"

if [[ ! -x "$BENCH_BIN" ]]; then
  echo "Error: binary missing after build: $BENCH_BIN" >&2
  exit 1
fi

SUMMARY="$OUT_DIR/summary_t${THREADS}.txt"
: >"$SUMMARY"

echo "NOTE|outer|runs|$OUTER_RUNS|bin|$BIN_FILE|threads|$THREADS|sync_k|$SYNC_K|s|$S|min_w|$MIN_W" | tee -a "$SUMMARY"

for run in $(seq 1 "$OUTER_RUNS"); do
  log="$OUT_DIR/run_$(printf '%02d' "$run")_t${THREADS}.txt"
  "$BENCH_BIN" "$BIN_FILE" "$THREADS" "$SYNC_K" "$S" "$MIN_W" | tee "$log"
  echo "OUTER_RUN|$run|log|$log" | tee -a "$SUMMARY"
done

python3 "$ROOT/experiments/rigorous_bench/scripts/analyze_hot_compare.py" \
  --simd-glob "$OUT_DIR/run_*_t${THREADS}.txt" \
  --xseal-glob "/dev/null" \
  -o "$OUT_DIR/report_t${THREADS}.md" 2>/dev/null || true

echo "Wrote $SUMMARY and per-run logs"
