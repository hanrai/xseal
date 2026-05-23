#!/usr/bin/env bash

# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2026 hanrai. All Rights Reserved.

# Vanilla Minimap2 mm_sketch only — perf instructions for insn/bp (no idx sort/post).
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/../../.." && pwd)"
MM="$ROOT/external/minimap2"
BIN="$MM/mm_sketch_bench"
FASTA="${FASTA:-$ROOT/data/hg38.fa}"
K="${K:-31}"
W="${W:-21}"
THREADS="${THREADS:-1}"
INNER_LOOPS="${INNER_LOOPS:-3}"
OUT="${OUT:-$ROOT/experiments/rigorous_bench/reports/mm_sketch_perf}"

mkdir -p "$OUT"

if [[ ! -x "$BIN" ]]; then
  echo "[mm_sketch_perf] building $BIN ..." >&2
  make -C "$MM" mm_sketch_bench
fi

if [[ ! -f "$FASTA" ]]; then
  echo "FASTA not found: $FASTA" >&2
  exit 1
fi

if [[ "$K" -gt 28 ]]; then
  echo "[mm_sketch_perf] WARNING: vanilla mm_sketch requires k<=28; got K=$K, using K=21 W=11" >&2
  K=21
  W=11
fi

CACHE="${CACHE:-$OUT/hg38.mmsk}"
LOG="$OUT/run_meta.txt"
{
  echo "# mm_sketch_perf"
  echo "fasta=$FASTA"
  echo "k=$K w=$W threads=$THREADS inner_loops=$INNER_LOOPS"
  echo "date=$(date -Is)"
} | tee "$LOG"

if ! command -v perf >/dev/null 2>&1; then
  echo "perf not found; running bench without PMU" | tee -a "$LOG"
  "$BIN" "$FASTA" "$K" "$W" "$THREADS" "$INNER_LOOPS" | tee "$OUT/bench.txt"
  exit 0
fi

echo "[mm_sketch_perf] dumping cache $CACHE (outside perf) ..." | tee -a "$LOG"
"$BIN" --dump "$CACHE" "$FASTA"

echo "[mm_sketch_perf] perf on --from-dump (sketch only, no FASTA parse) ..." | tee -a "$LOG"
perf stat -e instructions,cycles,branches,branch-misses -- \
  "$BIN" --from-dump "$CACHE" "$K" "$W" "$THREADS" "$INNER_LOOPS" \
  2>"$OUT/perf_stat.txt" | tee "$OUT/bench.txt"

python3 - "$OUT" <<'PY'
import re, sys
from pathlib import Path

out = Path(sys.argv[1])
perf_txt = (out / "perf_stat.txt").read_text()
bench_txt = (out / "bench.txt").read_text()

m_ins = re.search(r"([\d,]+)\s+instructions", perf_txt)
m_cyc = re.search(r"([\d,]+)\s+cycles", perf_txt)
m_res = re.search(
    r"RESULT\|mm_sketch_bench\|bases\|(\d+)\|hits\|(\d+)\|loops\|(\d+)\|", bench_txt
)
if not m_ins or not m_res:
    print("parse failed; see", out, file=sys.stderr)
    sys.exit(1)

ins = int(m_ins.group(1).replace(",", ""))
bases = int(m_res.group(1))
loops = int(m_res.group(3))
denom = bases * loops
insn_bp = ins / denom if denom else 0.0
ipc = ins / int(m_cyc.group(1).replace(",", "")) if m_cyc else float("nan")

summary = (
    f"instructions={ins:,}\n"
    f"bases={bases:,} loops={loops} denom_bases={denom:,}\n"
    f"insn_per_bp={insn_bp:.4f}\n"
    f"ipc={ipc:.4f}\n"
)
(out / "summary.txt").write_text(summary)
print(summary, end="")
PY

echo "Wrote $OUT/summary.txt"
