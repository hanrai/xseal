#!/usr/bin/env bash

# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2026 hanrai. All Rights Reserved.

set -euo pipefail
# shellcheck source=experiments/scripts/common.sh
source "$(cd "$(dirname "$0")" && pwd)/common.sh"
load_config

mkdir -p "$REPO_ROOT/build/bin" "$EXPERIMENTS_ROOT/results"
OUT="$EXPERIMENTS_ROOT/results/metrics_part2.json"
LOG="$EXPERIMENTS_ROOT/results/part2_pmu.log"
SKIP="$EXPERIMENTS_ROOT/results/part2_pmu.skip"
: >"$LOG"

STLF="$REPO_ROOT/build/bin/stlf_penalty_microbench"
WR="$REPO_ROOT/build/bin/window_reduction_microbench"
clang++ -O3 -march=native -std=c++20 "$REPO_ROOT/experiments/rigorous_bench/stlf_penalty_microbench.cpp" -o "$STLF" 2>>"$LOG" || true
clang++ -O3 -march=native -std=c++20 "$REPO_ROOT/experiments/rigorous_bench/window_reduction_microbench.cpp" -o "$WR" 2>>"$LOG" || true

STLF_ITERS="${STLF_ITERS:-500000}"
FOLD_ITERS="${FOLD_SPARSE_ITERS:-2000000}"

if [[ -x "$STLF" ]]; then
  echo "== stlf_penalty_microbench ==" | tee -a "$LOG"
  "$STLF" good "$STLF_ITERS" | tee -a "$LOG" || true
  "$STLF" bad "$STLF_ITERS" | tee -a "$LOG" || true
else
  echo "stlf compile failed" | tee "$SKIP"
fi

if [[ "${PART2_PMU:-0}" == "1" ]] && command -v perf >/dev/null 2>&1 && [[ -x "$WR" ]]; then
  echo "== perf fold ==" | tee -a "$LOG"
  if ! perf stat -r 2 -e instructions,cycles,ls_bad_status2.stli_other -- "$WR" fold "$FOLD_ITERS" >>"$LOG" 2>&1; then
    echo "perf fold failed" | tee -a "$SKIP"
  fi
  echo "== perf sparse ==" | tee -a "$LOG"
  perf stat -r 2 -e instructions,cycles,ls_bad_status2.stli_other -- "$WR" sparse "$FOLD_ITERS" >>"$LOG" 2>&1 || true
else
  echo "[part2] skip perf fold/sparse (PART2_PMU=1 needs perf + window_reduction_microbench)" | tee -a "$LOG"
fi

export LOG OUT
python3 <<'PY'
import json, os, re
from pathlib import Path

log = Path(os.environ["LOG"]).read_text(errors="ignore")
out: dict = {}


def perf_block(label: str, block: str) -> None:
    m = re.search(r"([\d,]+)\s+instructions", block)
    ins = float(m.group(1).replace(",", "")) if m else None
    m = re.search(r"([\d,]+)\s+cycles", block)
    cyc = float(m.group(1).replace(",", "")) if m else None
    m = re.search(r"([\d,]+)\s+ls_bad_status2\.stli_other", block)
    st = float(m.group(1).replace(",", "")) if m else None
    if ins is not None:
        out[label] = {"instructions": ins, "cycles": cyc, "stli_other": st, "unit": "counts"}


if "== perf fold ==" in log:
    blk = log.split("== perf fold ==")[1].split("== perf sparse ==")[0]
    perf_block("pmu_fold_window_reduce", blk)
if "== perf sparse ==" in log:
    blk = log.split("== perf sparse ==")[-1]
    perf_block("pmu_sparse_window_reduce", blk)

Path(os.environ["OUT"]).write_text(json.dumps(out, indent=2))
print(json.dumps(out, indent=2)[:2000])
PY

echo "Wrote $OUT"
