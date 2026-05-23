#!/usr/bin/env bash

# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2026 hanrai. All Rights Reserved.

set -euo pipefail
# shellcheck source=experiments/scripts/common.sh
source "$(cd "$(dirname "$0")" && pwd)/common.sh"
load_config

OUT="$EXPERIMENTS_ROOT/results/metrics_minimap2_index.json"
SKIP="$EXPERIMENTS_ROOT/results/minimap2_index.skip"
LOG="$EXPERIMENTS_ROOT/results/minimap2_index.log"

if [[ "${RUN_MINIMAP2_E2E:-0}" != "1" ]]; then
  echo "RUN_MINIMAP2_E2E!=1; skip minimap2 index timing" | tee "$SKIP"
  echo '{}' >"$OUT"
  exit 0
fi

INPUT="${HG38_FASTA:-$REPO_ROOT/data/hg38.fa}"
if [[ ! -f "$INPUT" ]]; then
  echo "no hg38" | tee "$SKIP"
  echo '{}' >"$OUT"
  exit 0
fi

VENV="$EXPERIMENTS_ROOT/vendor_versions.env"
# shellcheck source=/dev/null
set -a && source "$VENV" && set +a

MM_DIR="$REPO_ROOT/external/minimap2"
if [[ ! -f "$MM_DIR/Makefile" ]]; then
  mkdir -p "$(dirname "$MM_DIR")"
  rm -rf "$MM_DIR"
  if ! git clone --depth 1 --branch "${MINIMAP2_GIT_REF:-v2.28}" \
    "${MINIMAP2_REPO_URL:-https://github.com/lh3/minimap2.git}" "$MM_DIR" 2>>"$EXPERIMENTS_ROOT/results/minimap2_clone.log"; then
    echo "minimap2 clone failed" | tee "$SKIP"
    echo '{}' >"$OUT"
    exit 0
  fi
fi

(
  cd "$MM_DIR" && make -j"$(nproc)" 2>>"$EXPERIMENTS_ROOT/results/minimap2_build.log"
) || {
  echo "minimap2 make failed" | tee "$SKIP"
  echo '{}' >"$OUT"
  exit 0
}

M2="$MM_DIR/minimap2"
IDX="${MINIMAP2_INDEX_OUT:-/tmp/hg38_paper_repro_index.mmi}"
rm -f "$IDX"
THREADS="${MINIMAP2_THREADS:-12}"

export M2="$M2" IDX="$IDX" INPUT="$INPUT" THREADS="$THREADS" LOG="$LOG" OUT="$OUT" SKIP="$SKIP"
python3 <<'PY'
import json
import os
import subprocess
import sys
import time
from pathlib import Path

m2 = os.environ["M2"]
idx = os.environ["IDX"]
inp = os.environ["INPUT"]
threads = int(os.environ["THREADS"])
log = Path(os.environ["LOG"])
out = Path(os.environ["OUT"])
skip = Path(os.environ["SKIP"])

t0 = time.perf_counter()
with log.open("w") as fp:
    r = subprocess.run(
        [m2, "-t", str(threads), "-d", idx, inp],
        stdout=fp,
        stderr=subprocess.STDOUT,
    )
elapsed = time.perf_counter() - t0
if r.returncode != 0:
    skip.write_text(f"minimap2 exit {r.returncode}\n", encoding="utf-8")
    out.write_text("{}", encoding="utf-8")
    sys.exit(0)
obj = {
    "minimap2_index_total_s": {
        "mean": elapsed,
        "stdev": 0.0,
        "n": 1,
        "unit": "s",
        "non_portable": True,
    }
}
out.write_text(json.dumps(obj, indent=2), encoding="utf-8")
print(json.dumps(obj, indent=2))
PY

echo "Wrote $OUT"
