#!/usr/bin/env bash

# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2026 hanrai. All Rights Reserved.

set -euo pipefail
# shellcheck source=experiments/scripts/common.sh
source "$(cd "$(dirname "$0")" && pwd)/common.sh"
load_config

OUT="$EXPERIMENTS_ROOT/results/metrics_minimap2_xseal_e2e.json"
SKIP="$EXPERIMENTS_ROOT/results/minimap2_xseal_e2e.skip"
LOG="$EXPERIMENTS_ROOT/results/minimap2_xseal_e2e.log"

if [[ "${RUN_MINIMAP2_XSEAL_E2E:-0}" != "1" ]]; then
  echo "RUN_MINIMAP2_XSEAL_E2E!=1; skip minimap2+XSeal segmented repro" | tee "$SKIP"
  echo '{}' >"$OUT"
  exit 0
fi

INPUT="${HG38_FASTA:-$REPO_ROOT/data/hg38.fa}"
if [[ ! -f "$INPUT" ]]; then
  echo "no hg38 at $INPUT" | tee "$SKIP"
  echo '{}' >"$OUT"
  exit 0
fi

MM_DIR="$REPO_ROOT/integration/minimap2_xseal_hardened"
BIN="${MINIMAP2_XSEAL_BIN:-$MM_DIR/minimap2_xseal}"
THREADS="${MINIMAP2_THREADS:-12}"
IDX="${MINIMAP2_XSEAL_INDEX_OUT:-/tmp/hg38_xseal_paper_repro.mmi}"

if [[ ! -x "$BIN" ]]; then
  echo "[66] building minimap2_xseal in $MM_DIR ..." >&2
  if ! make -C "$MM_DIR" -j"$(nproc)" 2>>"$EXPERIMENTS_ROOT/results/minimap2_xseal_build.log"; then
    echo "make minimap2_xseal failed (see minimap2_xseal_build.log)" | tee "$SKIP"
    echo '{}' >"$OUT"
    exit 0
  fi
fi

if [[ ! -x "$BIN" ]]; then
  echo "no executable at $BIN" | tee "$SKIP"
  echo '{}' >"$OUT"
  exit 0
fi

rm -f "$IDX"
export BIN INPUT IDX THREADS LOG OUT SKIP
python3 <<'PY'
import json
import os
import re
import subprocess
import sys
from pathlib import Path

bin_p = os.environ["BIN"]
inp = os.environ["INPUT"]
idx = os.environ["IDX"]
threads = int(os.environ["THREADS"])
log = Path(os.environ["LOG"])
out = Path(os.environ["OUT"])
skip = Path(os.environ["SKIP"])

log.parent.mkdir(parents=True, exist_ok=True)
with log.open("w") as fp:
    r = subprocess.run(
        [bin_p, "-t", str(threads), "-d", idx, inp],
        stdout=fp,
        stderr=subprocess.STDOUT,
        text=True,
    )
text = log.read_text(encoding="utf-8", errors="replace")
if r.returncode != 0:
    skip.write_text(f"minimap2_xseal exit {r.returncode}\n", encoding="utf-8")
    out.write_text("{}", encoding="utf-8")
    sys.exit(0)

pat = re.compile(
    r"^RESULT\|([^|]+)\|mean\|([\d.]+)\|std\|([\d.]+)\|trials\|(\d+)\s*$",
    re.MULTILINE,
)
obj = {}
for m in pat.finditer(text):
    key, mean_s, std_s, trials = m.group(1), m.group(2), m.group(3), m.group(4)
    obj[key] = {
        "mean": float(mean_s),
        "stdev": float(std_s),
        "n": int(trials),
        "unit": "GiB" if key.endswith("_gib") else "s",
    }
if not obj:
    skip.write_text("no RESULT lines in log (see minimap2_xseal_e2e.log)\n", encoding="utf-8")
out.write_text(json.dumps(obj, indent=2), encoding="utf-8")
print(json.dumps(obj, indent=2))
PY

echo "Wrote $OUT"
