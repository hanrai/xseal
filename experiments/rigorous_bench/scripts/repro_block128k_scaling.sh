#!/usr/bin/env bash

# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2026 hanrai. All Rights Reserved.

# 128 KiB 2-bit block dispatch: work queue hands exact 128k packed blocks to all workers.
# Per block: scan -> harvest (hash_mode=1) or hit count only (hash_mode=0).
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/../../.." && pwd)"
BIN_FILE="${BIN_FILE:-$ROOT/data/hg38_2bit.bin}"
SYNC_K="${SYNC_K:-31}"
S="${S:-11}"
MIN_W="${MIN_W:-21}"
OUTER_RUNS="${OUTER_RUNS:-10}"
THREAD_LIST="${THREAD_LIST:-1,2,4,5,6,8,11,12,16}"
HASH_MODE="${HASH_MODE:-0}"

OUT_DIR="${OUT_DIR:-$ROOT/experiments/rigorous_bench/repro_scaling_out/block128k/pos_only}"

if [[ "$HASH_MODE" == "1" ]]; then
  OUT_DIR="${OUT_DIR%/pos_only}/hash"
  TOOL_TAG_X="XsealBlock128k"
  TOOL_TAG_S="SimdMinBlock128k"
else
  TOOL_TAG_X="XsealBlock128kCnt"
  TOOL_TAG_S="SimdMinBlock128kCnt"
fi

XSEAL_BIN="$ROOT/experiments/rigorous_bench/xseal_block128k_bench"
SIMD_CRATE="simd_min_block128k_bench"
SIMD_BIN="$ROOT/experiments/rigorous_bench/${SIMD_CRATE}/target/release/${SIMD_CRATE}"

if [[ ! -f "$BIN_FILE" ]]; then
  echo "error: missing $BIN_FILE" >&2
  exit 1
fi

mkdir -p "$OUT_DIR"

echo "[build] xseal_block128k_bench (hash_mode=$HASH_MODE)" >&2
clang++ -O3 -DNDEBUG -mavx2 -mbmi2 -mfma -ffast-math -march=native -std=c++20 -pthread -I"$ROOT/include" \
  "$ROOT/experiments/rigorous_bench/xseal_block128k_bench.cpp" -o "$XSEAL_BIN"

echo "[build] $SIMD_CRATE" >&2
(
  cd "$ROOT/experiments/rigorous_bench/${SIMD_CRATE}"
  env -u CARGO_TARGET_DIR RUSTFLAGS="-C target-cpu=native" cargo build --release -q
)

META="$OUT_DIR/run_meta.txt"
{
  echo "benchmark=block128k_dispatch"
  echo "block_bytes=131072"
  echo "protocol=$([[ "$HASH_MODE" == "1" ]] && echo block128k_pos_and_hash || echo block128k_pos_count)"
  echo "thread_list=$THREAD_LIST"
  echo "outer_runs=$OUTER_RUNS"
  echo "hash_mode=$HASH_MODE"
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

  echo "=== threads=$t ($TOOL_TAG_X / $TOOL_TAG_S) block128k ===" | tee -a "$LOG"
  echo "NOTE|scaling|benchmark=block128k|threads=$t|outer=$OUTER_RUNS|hash_mode=$HASH_MODE" | tee -a "$LOG"

  for run in $(seq 1 "$OUTER_RUNS"); do
    echo "--- outer run $run ---" | tee -a "$LOG"
    "$XSEAL_BIN" "$BIN_FILE" "$t" "$HASH_MODE" "$SYNC_K" "$S" "$MIN_W" | tee -a "$LOG"
    "$SIMD_BIN" "$BIN_FILE" "$t" "$HASH_MODE" "$SYNC_K" "$S" "$MIN_W" | tee -a "$LOG"
  done
done

echo "finished=$(date -Iseconds)" >>"$META"

export OUT_DIR HASH_MODE TOOL_TAG_X TOOL_TAG_S
python3 "$ROOT/experiments/rigorous_bench/scripts/summarize_thread_scaling.py" \
  --out-dir "$OUT_DIR" \
  --hash-mode "$HASH_MODE" \
  --tool-x "$TOOL_TAG_X" \
  --tool-s "$TOOL_TAG_S" \
  --title "128 KiB block dispatch scaling" \
  --no-paper-row \
  | tee "$OUT_DIR/scaling_summary.txt"

echo "Done. Results: $OUT_DIR/scaling_summary.txt"
