#!/usr/bin/env bash

# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2026 hanrai. All Rights Reserved.

# Genome-partition baseline (legacy hot bench): each worker owns a large contiguous
# slice of hg38 (~num_bases/threads). Preserves prior scaling experiments that showed
# chunk size vs throughput effects. Does NOT use 128 KiB block dispatch.
#
# Outputs: repro_scaling_out/genome_partition/{pos_only,hash} (or legacy pos_only/hash paths)
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/../../.." && pwd)"
BIN_FILE="${BIN_FILE:-$ROOT/data/hg38_2bit.bin}"
SYNC_K="${SYNC_K:-31}"
S="${S:-11}"
MIN_W="${MIN_W:-21}"
OUTER_RUNS="${OUTER_RUNS:-10}"
THREAD_LIST="${THREAD_LIST:-2,4,5,6,8,11,12,16}"
HASH_MODE="${HASH_MODE:-0}"

# Default new layout; set OUT_DIR explicitly to keep writing legacy repro_scaling_out/pos_only
OUT_DIR="${OUT_DIR:-$ROOT/experiments/rigorous_bench/repro_scaling_out/genome_partition/pos_only}"

if [[ "$HASH_MODE" == "1" ]]; then
  OUT_DIR="${OUT_DIR%/pos_only}/hash"
  [[ "$OUT_DIR" == *genome_partition* ]] || OUT_DIR="$ROOT/experiments/rigorous_bench/repro_scaling_out/genome_partition/hash"
  XSEAL_SRC="$ROOT/experiments/rigorous_bench/xseal_cache_resident_hash_bench.cpp"
  XSEAL_BIN="$ROOT/experiments/rigorous_bench/xseal_hot_hash_bench"
  SIMD_CRATE_DIR="$ROOT/experiments/rigorous_bench/simdmin_cache_resident_hash_bench"
  SIMD_BIN="$ROOT/experiments/rigorous_bench/simd_min_hot_hash_bench/target/release/simd_min_hot_hash_bench"
  TOOL_TAG_X="XsealHotHash"
  TOOL_TAG_S="SimdMinHotHash"
else
  XSEAL_SRC="$ROOT/experiments/rigorous_bench/xseal_cache_resident_scan_bench.cpp"
  XSEAL_BIN="$ROOT/experiments/rigorous_bench/xseal_hot_bench"
  SIMD_CRATE_DIR="$ROOT/experiments/rigorous_bench/simdmin_cache_resident_scan_bench"
  SIMD_BIN="$ROOT/experiments/rigorous_bench/simd_min_hot_bench/target/release/simd_min_hot_bench"
  TOOL_TAG_X="XsealHot"
  TOOL_TAG_S="SimdMinHot"
fi

if [[ ! -f "$BIN_FILE" ]]; then
  echo "error: missing $BIN_FILE" >&2
  exit 1
fi

mkdir -p "$OUT_DIR"

echo "[build] xseal genome-partition ($HASH_MODE hash mode)" >&2
clang++ -O3 -DNDEBUG -mavx2 -mbmi2 -mfma -ffast-math -march=native -std=c++20 -pthread \
  -I"$ROOT/include" "$XSEAL_SRC" -o "$XSEAL_BIN"

echo "[build] simd hot hash/scan bench" >&2
(
  cd "$SIMD_CRATE_DIR"
  env -u CARGO_TARGET_DIR RUSTFLAGS="-C target-cpu=native" cargo build --release -q
)
mkdir -p "$(dirname "$SIMD_BIN")"
cp -f "$SIMD_CRATE_DIR/target/release/$(basename "$SIMD_CRATE_DIR")" "$SIMD_BIN"

META="$OUT_DIR/run_meta.txt"
{
  echo "benchmark=genome_partition"
  echo "protocol=$([[ "$HASH_MODE" == "1" ]] && echo pos_and_hash || echo pos_only)"
  echo "thread_list=$THREAD_LIST"
  echo "outer_runs=$OUTER_RUNS"
  echo "K=$SYNC_K S=$S MIN_W=$MIN_W"
  echo "bin_file=$BIN_FILE"
  echo "started=$(date -Iseconds)"
} >"$META"

IFS=',' read -r -a THREADS_ARR <<<"$THREAD_LIST"

for t in "${THREADS_ARR[@]}"; do
  t="$(echo "$t" | tr -d ' ')"
  [[ -z "$t" ]] && continue
  TDIR="$OUT_DIR/t${t}"
  mkdir -p "$TDIR"
  LOG="$TDIR/scaling_${t}t.log"
  : >"$LOG"

  echo "=== threads=$t ($TOOL_TAG_X / $TOOL_TAG_S) genome_partition ===" | tee -a "$LOG"
  echo "NOTE|scaling|benchmark=genome_partition|threads=$t|outer=$OUTER_RUNS|hash_mode=$HASH_MODE" | tee -a "$LOG"

  for run in $(seq 1 "$OUTER_RUNS"); do
    echo "--- outer run $run ---" | tee -a "$LOG"
    "$XSEAL_BIN" "$BIN_FILE" "$t" "$SYNC_K" "$S" "$MIN_W" | tee -a "$LOG"
    "$SIMD_BIN" "$BIN_FILE" "$t" "$SYNC_K" "$S" "$MIN_W" | tee -a "$LOG"
  done
done

echo "finished=$(date -Iseconds)" >>"$META"

if [[ "$HASH_MODE" == "1" ]]; then
  BASELINE_DIR="$ROOT/experiments/rigorous_bench/repro_tab1_hash_out"
else
  BASELINE_DIR="$ROOT/experiments/rigorous_bench/repro_tab1_tab3_out"
fi

export OUT_DIR HASH_MODE TOOL_TAG_X TOOL_TAG_S
python3 "$ROOT/experiments/rigorous_bench/scripts/summarize_thread_scaling.py" \
  --out-dir "$OUT_DIR" \
  --baseline-1t-dir "$BASELINE_DIR" \
  --hash-mode "$HASH_MODE" \
  --title "Genome-partition scaling" \
  | tee "$OUT_DIR/scaling_summary.txt"

echo "Done. Results: $OUT_DIR/scaling_summary.txt"
