#!/usr/bin/env bash

# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2026 hanrai. All Rights Reserved.

set -euo pipefail
# shellcheck source=experiments/scripts/common.sh
source "$(cd "$(dirname "$0")" && pwd)/common.sh"
load_config

OUT="$EXPERIMENTS_ROOT/results/metrics_rust_sota.json"
SKIP="$EXPERIMENTS_ROOT/results/rust_sota.skip"
LOG="$EXPERIMENTS_ROOT/results/rust_sota.log"
FETCH_LOG="$EXPERIMENTS_ROOT/results/rust_sota_fetch.log"

if [[ "${RUN_RUST_SOTA:-0}" != "1" ]]; then
  echo "RUN_RUST_SOTA!=1; skip rust modular baseline" | tee "$SKIP"
  echo '{}' >"$OUT"
  exit 0
fi

INPUT="${HG38_FASTA:-$REPO_ROOT/data/hg38.fa}"
if [[ ! -f "$INPUT" ]]; then
  echo "no hg38 for rust_sota" | tee "$SKIP"
  echo '{}' >"$OUT"
  exit 0
fi

if [[ ! -f "$REPO_ROOT/external/simd-minimizers/Cargo.toml" ]]; then
  echo "simd-minimizers missing; run 05_fetch_vendor.sh (FETCH_VENDOR=1)" | tee "$SKIP"
  echo '{}' >"$OUT"
  exit 0
fi

if ! command -v cargo >/dev/null 2>&1; then
  echo "cargo not found in PATH; install Rust toolchain to run rust_sota_bench" | tee "$SKIP"
  echo '{}' >"$OUT"
  exit 0
fi

RS_DIR="$REPO_ROOT/external/rust_sota_bench"
if ! (cd "$RS_DIR" && cargo fetch 2>&1 | tee "$FETCH_LOG"); then
  {
    echo "cargo fetch failed (see rust_sota_fetch.log); attempting cargo build anyway."
  } | tee -a "$SKIP"
fi

if ! (cd "$RS_DIR" && env -u CARGO_TARGET_DIR RUSTFLAGS="-C target-cpu=native" cargo build --release -q); then
  echo "rust_sota_bench cargo build failed" | tee "$SKIP"
  echo '{}' >"$OUT"
  exit 0
fi

BIN="$RS_DIR/target/release/rust_sota_bench"
TH="${RUST_SOTA_THREADS:-12}"
"$BIN" "$INPUT" -t "$TH" --preload --rigorous 2>&1 | tee "$LOG"

export RUST_SOTA_LOG="$LOG" RUST_SOTA_OUT="$OUT"
python3 <<'PY'
import json
import os
import re

log_path = os.environ["RUST_SOTA_LOG"]
out_path = os.environ["RUST_SOTA_OUT"]
text = open(log_path, errors="ignore").read()

obj = {}

m = re.search(r"RESULT\|E2E_Gbp_s\|mean\|([\d.]+)\|std\|([\d.]+)\|trials\|(\d+)", text)
trials = 10
if m:
    trials = int(m.group(3))
    obj["integration_rust_modular_e2e_gbp_s"] = {
        "mean": float(m.group(1)),
        "stdev": float(m.group(2)),
        "n": trials,
        "unit": "Gbp/s",
    }

m2 = re.search(r"RESULT\|ParseMaterialize_Gbp_s\|mean\|([\d.]+)\|std\|([\d.]+)", text)
if m2:
    obj["integration_rust_modular_parse_gbp_s"] = {
        "mean": float(m2.group(1)),
        "stdev": float(m2.group(2)),
        "n": trials,
        "unit": "Gbp/s",
    }

m3 = re.search(r"RESULT\|SyncmerScan_Gbp_s\|mean\|([\d.]+)\|std\|([\d.]+)", text)
if m3:
    obj["integration_rust_modular_scan_gbp_s"] = {
        "mean": float(m3.group(1)),
        "stdev": float(m3.group(2)),
        "n": trials,
        "unit": "Gbp/s",
    }

open(out_path, "w").write(json.dumps(obj, indent=2))
print(json.dumps(obj, indent=2))
PY

echo "Wrote $OUT"
