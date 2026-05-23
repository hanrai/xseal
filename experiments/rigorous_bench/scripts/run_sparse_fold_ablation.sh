#!/usr/bin/env bash

# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2026 hanrai. All Rights Reserved.

# Reproducible fold vs sparse-table window reduction + perf counters (Zen 3).
# Run from repo root or any cwd; uses absolute paths below.
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
OUT="${ROOT}/experiments/rigorous_bench/reports/sparse_table_fold_ablation_20260515"
BIN="${ROOT}/experiments/rigorous_bench/window_reduce_bench"
ITERS="${1:-5000000}"

mkdir -p "$OUT"
{
  echo "# sparse_table_fold_ablation"
  echo "# date: $(date -Iseconds)"
  echo "# host: $(uname -n)"
  echo "# cpu: $(grep -m1 'model name' /proc/cpuinfo || true)"
  echo "# outer_iters: $ITERS"
  echo "# binary: $BIN"
  echo "# compile: clang++ -O3 -march=native -std=c++20 experiments/rigorous_bench/window_reduce_bench.cpp"
  echo ""
} | tee "$OUT/run_meta.txt"

echo "== build ==" | tee -a "$OUT/run_meta.txt"
clang++ -O3 -march=native -std=c++20 "${ROOT}/experiments/rigorous_bench/window_reduce_bench.cpp" -o "$BIN" 2>&1 | tee -a "$OUT/compile.log"

echo "" | tee -a "$OUT/run_meta.txt"
echo "== correctness (checksums must match) ==" | tee -a "$OUT/run_meta.txt"
"$BIN" fold 1 | tee "$OUT/checksum_fold.txt"
"$BIN" sparse 1 | tee "$OUT/checksum_sparse.txt"

echo "" | tee -a "$OUT/run_meta.txt"
echo "== perf fold (5 runs) ==" | tee -a "$OUT/run_meta.txt"
perf stat -r 5 -e instructions,cycles,ls_bad_status2.stli_other,ls_dc_accesses -- \
  "$BIN" fold "$ITERS" 2>&1 | tee "$OUT/perf_fold.txt"

echo "" | tee -a "$OUT/run_meta.txt"
echo "== perf sparse (5 runs) ==" | tee -a "$OUT/run_meta.txt"
perf stat -r 5 -e instructions,cycles,ls_bad_status2.stli_other,ls_dc_accesses -- \
  "$BIN" sparse "$ITERS" 2>&1 | tee "$OUT/perf_sparse.txt"

OUT_ESC="$OUT" python3 -c "
import re, pathlib, os
out = pathlib.Path(os.environ['OUT_ESC'])
def parse(path):
    t = pathlib.Path(path).read_text()
    ins = int(re.search(r'([\d,]+)\s+instructions', t).group(1).replace(',',''))
    cyc = int(re.search(r'([\d,]+)\s+cycles', t).group(1).replace(',',''))
    stl = int(re.search(r'([\d,]+)\s+ls_bad_status2\.stli_other', t).group(1).replace(',',''))
    return ins, cyc, stl
insf, cyf, stf = parse(out/'perf_fold.txt')
inss, cys, sts = parse(out/'perf_sparse.txt')
print(f'fold:   insn={insf:,}  cycles={cyf:,}  stli_other={stf:,}  stli/insn={stf/insf:.6g}')
print(f'sparse: insn={inss:,}  cycles={cys:,}  stli_other={sts:,}  stli/insn={sts/inss:.6g}')
print(f'ratio sparse/fold: cycles={cys/cyf:.3f}  stli={sts/stf:.3f}  insn={inss/insf:.3f}')
" | tee "$OUT/ratios.txt"

echo "Done. Artifacts under: $OUT"
