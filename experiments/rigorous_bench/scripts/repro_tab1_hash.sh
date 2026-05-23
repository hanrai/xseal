#!/usr/bin/env bash

# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2026 hanrai. All Rights Reserved.

# Table 1 style hot-RAM benchmark with explicit hash materialization.
# XSeal: scan + xseal_harvest_kmers (harvest layer).
# simd-minimizers: scan + values_u64().
# Does NOT touch repro_tab1_tab3_pos_only.sh or repro_tab1_tab3_out/.
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/../../.." && pwd)"
BIN_FILE="${BIN_FILE:-$ROOT/data/hg38_2bit.bin}"
THREADS="${THREADS:-1}"
SYNC_K="${SYNC_K:-31}"
S="${S:-11}"
MIN_W="${MIN_W:-21}"
OUTER_RUNS="${OUTER_RUNS:-10}"
OUT_DIR="${OUT_DIR:-$ROOT/experiments/rigorous_bench/repro_tab1_hash_out}"

XSEAL_BIN="$ROOT/experiments/rigorous_bench/xseal_hot_hash_bench"
SIMD_BIN="$ROOT/experiments/rigorous_bench/simd_min_hot_hash_bench/target/release/simd_min_hot_hash_bench"

if [[ ! -f "$BIN_FILE" ]]; then
  echo "error: missing $BIN_FILE" >&2
  exit 1
fi

mkdir -p "$OUT_DIR"

echo "[build] xseal_hot_hash_bench (scan + harvest layer)" >&2
clang++ -O3 -DNDEBUG -mavx2 -mbmi2 -mfma -ffast-math -march=native -std=c++20 -pthread \
  -I"$ROOT/include" \
  "$ROOT/experiments/rigorous_bench/xseal_hot_hash_bench.cpp" -o "$XSEAL_BIN"

echo "[build] simd_min_hot_hash_bench (scan + values_u64)" >&2
(
  cd "$ROOT/experiments/rigorous_bench/simd_min_hot_hash_bench"
  env -u CARGO_TARGET_DIR RUSTFLAGS="-C target-cpu=native" cargo build --release -q
)

LOG="$OUT_DIR/table1_hash_${THREADS}t.log"
: >"$LOG"

echo "=== Table 1 hash experiment (pos+hash, threads=$THREADS) ===" | tee -a "$LOG"
echo "NOTE|protocol|pos_and_hash|xseal_harvest|simd_values_u64|warmup=3|recorded=10|loops=13|K=$SYNC_K|S=$S|W=$MIN_W" | tee -a "$LOG"

for run in $(seq 1 "$OUTER_RUNS"); do
  echo "--- outer run $run ---" | tee -a "$LOG"
  "$XSEAL_BIN" "$BIN_FILE" "$THREADS" "$SYNC_K" "$S" "$MIN_W" | tee -a "$LOG"
  "$SIMD_BIN" "$BIN_FILE" "$THREADS" "$SYNC_K" "$S" "$MIN_W" | tee -a "$LOG"
done

export LOG THREADS OUTER_RUNS OUT_DIR
python3 <<'PY' | tee "$OUT_DIR/table1_hash_summary.txt"
import os, re, statistics as st
from pathlib import Path

log = Path(os.environ["LOG"]).read_text(errors="replace").splitlines()
threads = os.environ["THREADS"]
outer_runs = os.environ["OUTER_RUNS"]

inner_re = re.compile(
    r"^INNER\|(XsealHotHash|SimdMinHotHash)\|(ClosedSync|OpenSync|Minimizer)\|"
    + re.escape(threads)
    + r"\|([\d.]+)\|([\d.]+)\|(\d+)\|([\d.eE+-]+)\|([0-9a-fA-F]+)"
)

vals = {}
digests = {}
for line in log:
    m = inner_re.match(line)
    if not m:
        continue
    key = (m.group(1), m.group(2))
    vals.setdefault(key, []).append(float(m.group(3)))
    digests[key] = m.group(6)

print(f"Hash experiment summary (threads={threads}, {outer_runs} outer runs)")
print(f"{'Tool':<16} {'Mode':<12} {'Gbp/s':>10} {'±σ':>8} {'Hits':>12} {'Digest':>18}")
print("-" * 72)
for key in sorted(vals.keys()):
    mean = st.mean(vals[key])
    std = st.stdev(vals[key]) if len(vals[key]) > 1 else 0.0
  # last hits from final matching line
    hits = 0
    for line in log:
        m = inner_re.match(line)
        if m and (m.group(1), m.group(2)) == key:
            hits = int(m.group(5))
    print(
        f"{key[0]:<16} {key[1]:<12} {mean:10.4f} {std:8.4f} {hits:12d} {digests.get(key, '?'):>18}"
    )

# Speedup vs simd per mode
print()
print("XSeal/simd speedup (Minimizer):")
xs = vals.get(("XsealHotHash", "Minimizer"), [])
sm = vals.get(("SimdMinHotHash", "Minimizer"), [])
if xs and sm:
    print(f"  {st.mean(xs)/st.mean(sm):.3f}x")
PY

echo "Done. Log: $LOG"
echo "Summary: $OUT_DIR/table1_hash_summary.txt"
