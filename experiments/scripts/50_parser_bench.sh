#!/usr/bin/env bash

# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2026 hanrai. All Rights Reserved.

set -euo pipefail
# shellcheck source=experiments/scripts/common.sh
source "$(cd "$(dirname "$0")" && pwd)/common.sh"
load_config

cd "$REPO_ROOT"
SKIP="$EXPERIMENTS_ROOT/results/parser_bench.skip"
mkdir -p "$REPO_ROOT/build/bin"
PB="$REPO_ROOT/build/bin/parser_throughput_bench"

if [[ ! -f "$REPO_ROOT/experiments/rigorous_bench/parser_throughput_bench.cpp" ]]; then
  echo "parser_throughput_bench.cpp missing; skip" | tee "$SKIP"
  exit 0
fi

if ! clang++ -O3 -march=native -std=c++20 -I"$REPO_ROOT/include" -I"$REPO_ROOT/external/minimap2" \
  "$REPO_ROOT/experiments/rigorous_bench/parser_throughput_bench.cpp" \
  -o "$PB" -lz 2>"$EXPERIMENTS_ROOT/results/parser_bench.build.log"; then
  echo "parser_throughput_bench target failed to build (see parser_bench.build.log)" | tee "$SKIP"
  exit 0
fi

INPUT="${HG38_FASTA:-$REPO_ROOT/data/hg38.fa}"
if [[ ! -f "$INPUT" ]]; then
  echo "no hg38 for parser bench" | tee "$SKIP"
  exit 0
fi

"$PB" "$INPUT" 2>&1 | tee "$EXPERIMENTS_ROOT/results/parser_bench.log" || {
  echo "parser_bench run failed" | tee "$SKIP"
  exit 0
}

python3 <<PY
import json, re
text = open("$EXPERIMENTS_ROOT/results/parser_bench.log", errors="ignore").read()
# parser_bench prints decimal GB/s (file_size/1e9 / time); convert to GiB/s for paper anchor.
m = re.search(r"Throughput\s*:\s*([\d.]+)\s*GB/s", text)
gb = float(m.group(1)) if m else None
gib = gb * (1e9 / (1024**3)) if gb is not None else None
open("$EXPERIMENTS_ROOT/results/metrics_parser.json", "w").write(
    json.dumps(
        {
            "parser_micro_xseal_gib_s": {
                "mean": gib,
                "stdev": 0.0,
                "n": 1,
                "unit": "GiB/s",
                "note": "converted from printed decimal GB/s",
            }
            if gib is not None
            else {"error": "unparsed"}
        },
        indent=2,
    )
)
PY

echo "Wrote parser_bench.log and metrics_parser.json (if parsed)"
