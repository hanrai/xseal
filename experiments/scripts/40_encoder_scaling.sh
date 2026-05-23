#!/usr/bin/env bash

# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2026 hanrai. All Rights Reserved.

set -euo pipefail
# shellcheck source=experiments/scripts/common.sh
source "$(cd "$(dirname "$0")" && pwd)/common.sh"
load_config

cd "$REPO_ROOT"
OUT="$EXPERIMENTS_ROOT/results/encoder_scaling.log"
ENC_BIN="$REPO_ROOT/build/bin/encoder_throughput_scaling_bench"
mkdir -p "$REPO_ROOT/build/bin"

if [[ ! -f "$REPO_ROOT/experiments/rigorous_bench/encoder_throughput_scaling_bench.cpp" ]]; then
  echo "encoder_throughput_scaling_bench.cpp missing; skip" | tee "$EXPERIMENTS_ROOT/results/encoder_scaling.skip"
  exit 0
fi

clang++ -O3 -march=native -std=c++20 -I"$REPO_ROOT/include" -I"$REPO_ROOT/include/xseal" \
  "$REPO_ROOT/experiments/rigorous_bench/encoder_throughput_scaling_bench.cpp" \
  -o "$ENC_BIN" 2>"$EXPERIMENTS_ROOT/results/encoder_scaling.compile.log" || {
  echo "encoder compile failed; see encoder_scaling.compile.log" | tee "$EXPERIMENTS_ROOT/results/encoder_scaling.skip"
  exit 0
}

: >"$OUT"
THREADS=(1 6 12)
MODES=("minimap2" "xseal")
FREQS=(0 80 10)
SIZE=1024
for f in "${FREQS[@]}"; do
  for t in "${THREADS[@]}"; do
    for m in "${MODES[@]}"; do
      echo "=== size=$SIZE threads=$t mode=$m freq=$f ===" | tee -a "$OUT"
      "$ENC_BIN" "$SIZE" "$t" "$m" "$f" 2>&1 | tee -a "$OUT" || true
    done
  done
done

export ENC_LOG="$OUT" ENC_JSON="$EXPERIMENTS_ROOT/results/metrics_encoder.json"
python3 <<'PY'
import json, re
from pathlib import Path
import os

log_path = Path(os.environ["ENC_LOG"])
out_path = Path(os.environ["ENC_JSON"])
text = log_path.read_text(errors="ignore")
pat = re.compile(
    r"RESULT\|(\w+)\|(\d+)\|(\d+)\|([\d.]+)\|([\d.]+)\|GiB/s"
)
obj = {}
for mode, threads, freq, mean, std in pat.findall(text):
    t, f = int(threads), int(freq)
    if (t, f) in ((6, 0), (6, 80), (1, 0), (1, 80)):
        key = f"encoder_{mode}_{t}t_freq{f}_gib_s"
        obj[key] = {
            "mean": float(mean),
            "stdev": float(std),
            "n": 1,
            "unit": "GiB/s",
        }

out_path.write_text(json.dumps(obj, indent=2))
print(json.dumps(obj, indent=2)[:2000])
PY

echo "Wrote $OUT and metrics_encoder.json"
