#!/usr/bin/env bash

# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2026 hanrai. All Rights Reserved.

# Paper Table 1 (tab:throughput) + Table 2 (tab:scaling) reproduction.
# Protocol: each process run = 3 warmup + 10 recorded loops (13 total);
#           10 outer process runs -> report μ±σ across outer inner-means.
# pos-only (no hash materialization); K=31 S=11 W=21; hg38_2bit.bin.
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/../../.." && pwd)"
BIN_FILE="${BIN_FILE:-$ROOT/data/hg38_2bit.bin}"
SYNC_K="${SYNC_K:-31}"
S="${S:-11}"
MIN_W="${MIN_W:-21}"
OUTER_RUNS="${OUTER_RUNS:-10}"
THREAD_LIST_TAB2="${THREAD_LIST_TAB2:-1,6,12}"
OUT_DIR="${OUT_DIR:-$ROOT/experiments/rigorous_bench/repro_paper_table12_out}"

XSEAL_HOT="$ROOT/experiments/rigorous_bench/xseal_hot_bench"
XSEAL_BLK="$ROOT/experiments/rigorous_bench/xseal_block128k_bench"
SIMD_HOT="$ROOT/experiments/rigorous_bench/simd_min_hot_bench/target/release/simd_min_hot_bench"
SIMD_BLK="$ROOT/experiments/rigorous_bench/simd_min_block128k_bench/target/release/simd_min_block128k_bench"

CXX_FLAGS=(
  -O3 -DNDEBUG -mavx2 -mbmi2 -mfma -ffast-math -march=native -std=c++20 -pthread
  -I"$ROOT/include"
)

if [[ ! -f "$BIN_FILE" ]]; then
  echo "error: missing $BIN_FILE" >&2
  exit 1
fi

mkdir -p "$OUT_DIR"

check_bins() {
  local missing=0
  for b in "$XSEAL_HOT" "$XSEAL_BLK" "$SIMD_HOT" "$SIMD_BLK"; do
    if [[ ! -x "$b" ]]; then
      echo "error: missing executable: $b" >&2
      missing=1
    fi
  done
  return "$missing"
}

if [[ "${SKIP_BUILD:-0}" != "1" ]]; then
  # shellcheck source=build_hot_benches.sh
  source "$(dirname "$0")/build_hot_benches.sh"
  echo "[build] xseal_block128k_bench" >&2
  clang++ "${CXX_FLAGS[@]}" \
    "$ROOT/experiments/rigorous_bench/xseal_fixed_block_scan_bench.cpp" -o "$XSEAL_BLK"
  echo "[build] simd_min_block128k_bench" >&2
  (
    cd "$ROOT/experiments/rigorous_bench/simdmin_fixed_block_scan_bench"
    env -u CARGO_TARGET_DIR RUSTFLAGS="-C target-cpu=native" cargo build --release -q
  )
  mkdir -p "$(dirname "$SIMD_BLK")"
  cp -f "$ROOT/experiments/rigorous_bench/simdmin_fixed_block_scan_bench/target/release/simdmin_fixed_block_scan_bench" \
    "$SIMD_BLK"
else
  echo "[skip build] SKIP_BUILD=1 — using existing binaries" >&2
  check_bins || exit 1
fi

run_hot_pair() {
  local t=$1 log=$2
  echo "=== genome_partition threads=$t ($OUTER_RUNS outer) ===" | tee -a "$log"
  echo "NOTE|benchmark=genome_partition|threads=$t|outer=$OUTER_RUNS|warmup=3|recorded=10" | tee -a "$log"
  for run in $(seq 1 "$OUTER_RUNS"); do
    echo "--- outer run $run ---" | tee -a "$log"
    "$XSEAL_HOT" "$BIN_FILE" "$t" "$SYNC_K" "$S" "$MIN_W" | tee -a "$log"
    "$SIMD_HOT" "$BIN_FILE" "$t" "$SYNC_K" "$S" "$MIN_W" | tee -a "$log"
  done
}

run_block128k_pair() {
  local t=$1 log=$2
  echo "=== block128k threads=$t ($OUTER_RUNS outer) ===" | tee -a "$log"
  echo "NOTE|benchmark=block128k|threads=$t|outer=$OUTER_RUNS|hash_mode=0|warmup=3|recorded=10" | tee -a "$log"
  for run in $(seq 1 "$OUTER_RUNS"); do
    echo "--- outer run $run ---" | tee -a "$log"
    "$XSEAL_BLK" "$BIN_FILE" "$t" 0 "$SYNC_K" "$S" "$MIN_W" | tee -a "$log"
    "$SIMD_BLK" "$BIN_FILE" "$t" 0 "$SYNC_K" "$S" "$MIN_W" | tee -a "$log"
  done
}

TAB1_LOG="$OUT_DIR/table1_t1.log"
: >"$TAB1_LOG"
run_hot_pair 1 "$TAB1_LOG"

IFS=',' read -r -a TAB2_THREADS <<<"$THREAD_LIST_TAB2"
for t in "${TAB2_THREADS[@]}"; do
  t="$(echo "$t" | tr -d ' ')"
  [[ -z "$t" ]] && continue
  TDIR="$OUT_DIR/t${t}"
  mkdir -p "$TDIR"
  run_hot_pair "$t" "$TDIR/genome_partition.log"
  run_block128k_pair "$t" "$TDIR/block128k.log"
done

{
  echo "bin_file=$BIN_FILE"
  echo "K=$SYNC_K S=$S MIN_W=$MIN_W"
  echo "outer_runs=$OUTER_RUNS"
  echo "thread_list_tab2=$THREAD_LIST_TAB2"
  echo "finished=$(date -Iseconds)"
} >"$OUT_DIR/run_meta.txt"

export OUT_DIR OUTER_RUNS
python3 "$ROOT/experiments/rigorous_bench/scripts/summarize_paper_table12.py" | tee "$OUT_DIR/PAPER_TABLE12_SUMMARY.txt"

echo "Done -> $OUT_DIR/PAPER_TABLE12_SUMMARY.txt"
