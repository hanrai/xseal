#!/usr/bin/env bash

# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2026 hanrai. All Rights Reserved.

# Build + run isolated s-mer hash benches (XSeal Murmur vs seq-hash NtHash SIMD).
# Syncmer K=31, w=21 => s=11. Optional perf stat for insn/bp.
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/../../.." && pwd)"
BIN_FILE="${BIN_FILE:--}"
NUM_BASES="${NUM_BASES:-134217728}"
THREADS="${THREADS:-1}"
SYNC_K="${SYNC_K:-31}"
S="${S:-11}"
SEED="${SEED:-42}"
OUTER_RUNS="${OUTER_RUNS:-3}"
RUN_PERF="${RUN_PERF:-1}"

XSEAL_BIN="${XSEAL_BIN:-$ROOT/experiments/rigorous_bench/hash_isolate_xseal}"
SIMD_BIN="${SIMD_BIN:-$ROOT/experiments/rigorous_bench/hash_isolate_simd/target/release/hash_isolate_simd}"
OUT_DIR="${OUT_DIR:-$ROOT/experiments/rigorous_bench/hash_isolate_out}"

mkdir -p "$OUT_DIR"

echo "[build] hash_isolate_xseal" >&2
clang++ -O3 -DNDEBUG -mavx2 -mbmi2 -mfma -ffast-math -march=native -std=c++20 -pthread \
  -I"$ROOT/include" \
  "$ROOT/experiments/rigorous_bench/hash_isolate_xseal.cpp" -o "$XSEAL_BIN"

echo "[build] hash_isolate_simd" >&2
(
  cd "$ROOT/experiments/rigorous_bench/hash_isolate_simd"
  env -u CARGO_TARGET_DIR RUSTFLAGS="-C target-cpu=native" cargo build --release -q
)

run_bench() {
  local tag=$1 bin=$2
  local log="$OUT_DIR/${tag}_t${THREADS}.txt"
  : >"$log"
  echo "NOTE|config|num_bases|$NUM_BASES|threads|$THREADS|sync_k|$SYNC_K|s|$S|seed|$SEED" | tee -a "$log"
  for run in $(seq 1 "$OUTER_RUNS"); do
    echo "=== $tag run $run/$OUTER_RUNS ===" | tee -a "$log"
    "$bin" "$BIN_FILE" "$NUM_BASES" "$THREADS" "$SYNC_K" "$S" "$SEED" | tee -a "$log"
  done
  echo "log|$log" >&2
}

run_bench xseal "$XSEAL_BIN"
run_bench simd "$SIMD_BIN"

if [[ "$RUN_PERF" == "1" ]]; then
  PERF_LOG="$OUT_DIR/perf_insn_t${THREADS}.txt"
  : >"$PERF_LOG"
  echo "=== perf stat (insn/bp denominator = NUM_BASES=$NUM_BASES) ===" | tee -a "$PERF_LOG"
  for tag in xseal simd; do
  bin="$XSEAL_BIN"
  [[ "$tag" == "simd" ]] && bin="$SIMD_BIN"
  echo "--- $tag ---" | tee -a "$PERF_LOG"
  perf stat -e instructions,cycles,branches,branch-misses -r 3 -- \
    "$bin" "$BIN_FILE" "$NUM_BASES" "$THREADS" "$SYNC_K" "$S" "$SEED" \
    2>&1 | tee -a "$PERF_LOG"
  python3 - "$PERF_LOG" "$tag" "$NUM_BASES" <<'PY' | tee -a "$PERF_LOG"
import re, sys
log_path, tag, num_bases = sys.argv[1], sys.argv[2], int(sys.argv[3])
text = open(log_path).read()
# last perf block for this tag
blocks = text.split(f"--- {tag} ---")
if len(blocks) < 2:
    sys.exit(0)
blk = blocks[-1].split("---")[0]
m = re.search(r"([\d,]+)\s+instructions", blk)
if not m:
    sys.exit(0)
ins = int(m.group(1).replace(",", ""))
insn_bp = ins / num_bases
print(f"SUMMARY|{tag}|instructions|{ins}|insn_per_bp|{insn_bp:.4f}")
PY
  done
fi

echo "" >&2
echo "Done. Logs: $OUT_DIR" >&2
echo "Parse RESULT lines for Gbp/s; SUMMARY lines for insn/bp." >&2
