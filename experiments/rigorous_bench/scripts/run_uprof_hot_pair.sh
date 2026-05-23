#!/usr/bin/env bash

# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2026 hanrai. All Rights Reserved.

# AMD uProf: assess + overview (overview/BPF needs root on this host).
# Outputs under experiments/rigorous_bench/uprof_assess_out/
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/../../.." && pwd)"
UPROF="${UPROF:-/opt/AMDuProf_5.3-518/bin/AMDuProfCLI}"
BIN_FILE="${BIN_FILE:-$ROOT/data/hg38_2bit.bin}"
THREADS="${THREADS:-1}"
SYNC_K="${SYNC_K:-31}"
S="${S:-11}"
MIN_W="${MIN_W:-21}"
OUT_ROOT="${OUT_ROOT:-$ROOT/experiments/rigorous_bench/uprof_assess_out}"
# overview = TBP + AVX FLOPs + function trace (NOT Top-Down pipeline_util).
# Pipeline Utilization (Retiring/Backend_Bound) needs GUI "Profile" with full PMU set;
# CLI --config assess|overview does not collect those TDCycles events.
# overview / BPF tracing: 1 = sudo collect (required if setcap-only failed)
USE_SUDO_OVERVIEW="${USE_SUDO_OVERVIEW:-1}"

XSEAL_BIN="${XSEAL_BIN:-$ROOT/experiments/rigorous_bench/xseal_hot_bench_g}"
SIMD_BIN="${SIMD_BIN:-$ROOT/experiments/rigorous_bench/simd_min_hot_bench/target/release/simd_min_hot_bench}"

if [[ ! -x "$UPROF" ]]; then
  echo "error: AMDuProfCLI not found at $UPROF" >&2
  exit 1
fi
if [[ ! -f "$BIN_FILE" ]]; then
  echo "error: missing 2bit input: $BIN_FILE" >&2
  exit 1
fi

BUILD_HOT_WITH_G=1 XSEAL_BIN="$XSEAL_BIN"
# shellcheck source=build_hot_benches.sh
source "$(dirname "$0")/build_hot_benches.sh"
(
  cd "$ROOT/experiments/rigorous_bench/simdmin_cache_resident_scan_bench"
  env -u CARGO_TARGET_DIR RUSTFLAGS="-C target-cpu=native -C debuginfo=1" \
    cargo build --release -q
)
cp -f "$ROOT/experiments/rigorous_bench/simdmin_cache_resident_scan_bench/target/release/simdmin_cache_resident_scan_bench" \
  "$SIMD_BIN"

run_collect() {
  local tag=$1 bin=$2
  local assess_dir="$OUT_ROOT/${tag}_assess"
  local overview_dir="$OUT_ROOT/${tag}_overview"
  mkdir -p "$assess_dir" "$overview_dir"

  echo "== $tag: assess (no sudo) =="
  "$UPROF" collect --config assess -o "$assess_dir" \
    "$bin" "$BIN_FILE" "$THREADS" "$SYNC_K" "$S" "$MIN_W"
  "$UPROF" report -i "$assess_dir" --category cpu

  echo "== $tag: overview (BPF; sudo if USE_SUDO_OVERVIEW=1) =="
  if [[ "$USE_SUDO_OVERVIEW" == "1" ]]; then
    sudo "$UPROF" collect --config overview -o "$overview_dir" \
      "$bin" "$BIN_FILE" "$THREADS" "$SYNC_K" "$S" "$MIN_W"
  else
    "$UPROF" collect --config overview -o "$overview_dir" \
      "$bin" "$BIN_FILE" "$THREADS" "$SYNC_K" "$S" "$MIN_W"
  fi
  # Report: overview sessions use a timestamped subdir; GUI is most reliable.
  local session
  session="$(find "$overview_dir" -maxdepth 1 -type d -name 'AMDuProf-*' | sort | tail -1)"
  if [[ -n "$session" ]]; then
    echo "overview session: $session"
    echo "  GUI (use --session on the directory, not session.uprof):"
    echo "    /opt/AMDuProf_5.3-518/bin/AMDuProf --session \"$session\" \\"
    echo "      --bin-path \"$ROOT/rigorous_bench;$ROOT/experiments/rigorous_bench/simd_min_hot_bench/target/release\""
    echo "  Black-window workaround: export AMDUPROF_GUI_WORKAROUND=1"
    echo "  CLI report (needs sudo on this host): sudo $UPROF report -i \"$session\" --category cpu"
  fi
}

run_collect xseal "$XSEAL_BIN"
run_collect simd "$SIMD_BIN"

echo "Done. Assess CSV: $OUT_ROOT/{xseal,simd}_assess/report.csv"
