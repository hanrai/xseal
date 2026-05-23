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
OUT_DIR="${OUT_DIR:-$ROOT/experiments/rigorous_bench/xseal_hot_out}"

if [[ ! -f "$BIN_FILE" ]]; then
  echo "Error: 2-bit packed input not found: $BIN_FILE" >&2
  exit 1
fi

mkdir -p "$OUT_DIR"
# shellcheck source=build_hot_benches.sh
source "$(dirname "$0")/build_hot_benches.sh"
BENCH_BIN="$XSEAL_BIN"

SUMMARY="$OUT_DIR/summary_t${THREADS}.txt"
: >"$SUMMARY"

echo "NOTE|outer|runs|$OUTER_RUNS|bin|$BIN_FILE|threads|$THREADS|sync_k|$SYNC_K|s|$S" | tee -a "$SUMMARY"

for run in $(seq 1 "$OUTER_RUNS"); do
  log="$OUT_DIR/run_$(printf '%02d' "$run")_t${THREADS}.txt"
  "$BENCH_BIN" "$BIN_FILE" "$THREADS" "$SYNC_K" "$S" "$MIN_W" | tee "$log"
  echo "OUTER_RUN|$run|log|$log" | tee -a "$SUMMARY"
done

echo "Wrote $SUMMARY and per-run logs"
