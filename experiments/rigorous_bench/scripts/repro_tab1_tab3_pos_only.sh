#!/usr/bin/env bash

# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2026 hanrai. All Rights Reserved.

# Reproduce paper Table 1 (tab:throughput) and verify Table 3 (tab:microarch).
# Both scanners: position-only output (no hash buffers; digest outside timed region).
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/../../.." && pwd)"
BIN_FILE="${BIN_FILE:-$ROOT/data/hg38_2bit.bin}"
THREADS="${THREADS:-1}"
SYNC_K="${SYNC_K:-31}"
S="${S:-11}"
MIN_W="${MIN_W:-21}"
OUTER_RUNS="${OUTER_RUNS:-10}"
OUT_DIR="${OUT_DIR:-$ROOT/experiments/rigorous_bench/repro_tab1_tab3_out}"

# shellcheck source=build_hot_benches.sh
source "$(dirname "$0")/build_hot_benches.sh"

if [[ ! -f "$BIN_FILE" ]]; then
  echo "error: missing $BIN_FILE (run: python3 experiments/rigorous_bench/scripts/pack_hg38.py data/hg38.fa data/hg38_2bit.bin)" >&2
  exit 1
fi

mkdir -p "$OUT_DIR"

TAB1_LOG="$OUT_DIR/table1_hot_${THREADS}t.log"
TAB1_SUM="$OUT_DIR/table1_summary.md"
: >"$TAB1_LOG"

echo "=== Table 1 reproduction (hot RAM, pos-only, threads=$THREADS) ===" | tee -a "$TAB1_LOG"
echo "NOTE|protocol|pos_only|warmup=3|recorded=10|loops=13|K=$SYNC_K|S=$S|W=$MIN_W" | tee -a "$TAB1_LOG"

for run in $(seq 1 "$OUTER_RUNS"); do
  echo "--- outer run $run ---" | tee -a "$TAB1_LOG"
  "$XSEAL_BIN" "$BIN_FILE" "$THREADS" "$SYNC_K" "$S" "$MIN_W" | tee -a "$TAB1_LOG"
  "$SIMD_BIN" "$BIN_FILE" "$THREADS" "$SYNC_K" "$S" "$MIN_W" | tee -a "$TAB1_LOG"
done

# Paper Table 1 anchors (experiments.tex)
export TAB1_LOG THREADS OUTER_RUNS OUT_DIR
python3 <<'PY' | tee "$OUT_DIR/table1_vs_paper.txt"
import os, re, statistics as st
from pathlib import Path

log_path = os.environ["TAB1_LOG"]
threads = os.environ["THREADS"]
outer_runs = os.environ["OUTER_RUNS"]
log = Path(log_path).read_text(errors="replace").splitlines()
paper = {
    ("SimdMinHot", "Minimizer"): 0.82,
    ("SimdMinHot", "ClosedSync"): 0.81,
    ("SimdMinHot", "OpenSync"): 0.84,
    ("XsealHot", "Minimizer"): 1.45,
    ("XsealHot", "ClosedSync"): 1.70,
    ("XsealHot", "OpenSync"): 1.74,
}
inner_re = re.compile(
    rf"^INNER\|(SimdMinHot|XsealHot)\|(ClosedSync|OpenSync|Minimizer)\|{threads}\|"
    r"([\d.]+)\|([\d.]+)\|"
)
vals = {}
for line in log:
    m = inner_re.match(line)
    if not m:
        continue
    key = (m.group(1), m.group(2))
    vals.setdefault(key, []).append(float(m.group(3)))

print(f"Table 1 vs paper (threads={threads}, mean Gbp/s over {outer_runs} outer runs)")
print(f"{'Tool':<12} {'Mode':<12} {'Measured':>10} {'Paper':>8} {'Δ%':>8}")
print("-" * 54)
for key in sorted(paper.keys()):
    if key not in vals:
        print(f"{key[0]:<12} {key[1]:<12} {'N/A':>10}")
        continue
    mean = st.mean(vals[key])
    exp = paper[key]
    delta = 100.0 * (mean - exp) / exp
    print(f"{key[0]:<12} {key[1]:<12} {mean:10.4f} {exp:8.2f} {delta:+7.1f}%")
PY

TAB3_PERF="$OUT_DIR/table3_perf.txt"
: >"$TAB3_PERF"

if command -v perf >/dev/null 2>&1; then
  echo "=== Table 3 PMU (perf stat, full hot bench, pos-only) ===" | tee "$TAB3_PERF"
  echo "--- XSeal ---" | tee -a "$TAB3_PERF"
  perf stat -r 3 -e instructions,cycles,branches,branch-misses \
    "$XSEAL_BIN" "$BIN_FILE" "$THREADS" "$SYNC_K" "$S" "$MIN_W" \
    >>"$TAB3_PERF" 2>&1 || echo "perf xseal failed" | tee -a "$TAB3_PERF"
  echo "--- simd-minimizers ---" | tee -a "$TAB3_PERF"
  perf stat -r 3 -e instructions,cycles,branches,branch-misses \
    "$SIMD_BIN" "$BIN_FILE" "$THREADS" "$SYNC_K" "$S" "$MIN_W" \
    >>"$TAB3_PERF" 2>&1 || echo "perf simd failed" | tee -a "$TAB3_PERF"

  export TAB3_PERF
  python3 <<'PY' | tee "$OUT_DIR/table3_vs_paper.txt"
import os, re
from pathlib import Path

text = Path(os.environ["TAB3_PERF"]).read_text(errors="replace")
paper = {
    "xseal": {"inst_B": 217.17, "ipc": 2.80, "time_s": 18.49},
    "simd": {"inst_B": 460.97, "ipc": 3.37, "time_s": 32.64},
}

def parse_block(label):
    if label not in text:
        return None
    blk = text.split(label)[-1]
    m = re.search(
        r"([\d,]+)\s+instructions\s+([\d.eE+-]+)\s+seconds.*?([\d,]+)\s+cycles",
        blk,
        re.S,
    )
    if not m:
        return None
    inst = float(m.group(1).replace(",", ""))
    secs = float(m.group(2))
    cyc = float(m.group(3).replace(",", ""))
    ipc = inst / cyc if cyc else 0.0
    return {"inst": inst, "secs": secs, "ipc": ipc}


def fmt_row(name, got, ref):
    if not got:
        print(f"{name}: no perf data")
        return
    inst_b = got["inst"] / 1e9
    print(f"{name}:")
    print(f"  instructions: {inst_b:.2f} B  (paper {ref['inst_B']:.2f} B, "
          f"ratio {inst_b/ref['inst_B']:.3f})")
    print(f"  IPC:          {got['ipc']:.2f}     (paper {ref['ipc']:.2f})")
    print(f"  wall time:    {got['secs']:.2f} s  (paper {ref['time_s']:.2f} s)")


fmt_row("XSeal", parse_block("--- XSeal ---"), paper["xseal"])
print()
fmt_row("simd-min", parse_block("--- simd-minimizers ---"), paper["simd"])
xs = parse_block("--- XSeal ---")
sm = parse_block("--- simd-minimizers ---")
if xs and sm:
    print(f"\ninst ratio simd/xseal: {sm['inst']/xs['inst']:.3f}  (paper {paper['simd']['inst_B']/paper['xseal']['inst_B']:.3f})")
PY
else
  echo "perf not found; skip Table 3 PMU" | tee "$OUT_DIR/table3_vs_paper.txt"
fi

echo "Done. Outputs:"
echo "  $TAB1_LOG"
echo "  $OUT_DIR/table1_vs_paper.txt"
echo "  $OUT_DIR/table3_vs_paper.txt"
