#!/usr/bin/env bash

# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2026 hanrai. All Rights Reserved.

# shellcheck shell=bash
# Loaded by other scripts in this directory.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")" && pwd)"
# EXPERIMENTS_ROOT = experiments/ directory (parent of scripts/)
export EXPERIMENTS_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"

load_config() {
  if [[ -f "$EXPERIMENTS_ROOT/config.env" ]]; then
    # shellcheck source=/dev/null
    source "$EXPERIMENTS_ROOT/config.env"
  elif [[ -f "$EXPERIMENTS_ROOT/config.env.example" ]]; then
    echo "[experiments] Using config.env.example (copy to config.env to customize)" >&2
    # shellcheck source=/dev/null
    source "$EXPERIMENTS_ROOT/config.env.example"
  else
    echo "Missing experiments/config.env and config.env.example" >&2
    exit 1
  fi

  # REPO_ROOT = repository root (parent of experiments/)
  : "${REPO_ROOT:=$(cd "$EXPERIMENTS_ROOT/.." && pwd)}"
  : "${BUILD_DIR:=build}"
  : "${WARMUP:=3}"
  : "${TRIALS:=10}"
  : "${USE_SHM:=0}"
  : "${E2E_PROFILE_ONCE:=1}"
  : "${GOVERNOR_CHECK:=0}"
  : "${PURE_SCAN_THREADS_LIST:=1,6,12}"
  : "${FETCH_VENDOR:=1}"
  : "${USE_HYPERFINE:=0}"
  : "${RUN_RUST_SOTA:=0}"
  : "${RUST_SOTA_THREADS:=12}"
  : "${RUN_MINIMAP2_E2E:=0}"
  : "${RUN_MINIMAP2_XSEAL_E2E:=0}"
  : "${PART2_PMU:=0}"
  : "${PART1_PERF:=0}"
  : "${RUN_PART1:=1}"
  : "${PART1_THREADS_LIST:=1,12}"

  if [[ -z "${XSEAL_E2E_THREADS:-}" ]]; then
    local n
    n=$(nproc)
    XSEAL_E2E_THREADS=$((n > 1 ? n - 1 : 1))
  fi

  : "${HG38_FASTA:=$REPO_ROOT/data/hg38.fa}"
  if [[ ! -f "$HG38_FASTA" && -f "$REPO_ROOT/data/hg38.fa" ]]; then
    HG38_FASTA="$REPO_ROOT/data/hg38.fa"
  fi

  export REPO_ROOT BUILD_DIR WARMUP TRIALS USE_SHM E2E_PROFILE_ONCE GOVERNOR_CHECK
  export XSEAL_E2E_THREADS PURE_SCAN_THREADS_LIST HG38_FASTA TWO_BIT_HG38
  export FETCH_VENDOR USE_HYPERFINE RUN_RUST_SOTA RUST_SOTA_THREADS RUN_MINIMAP2_E2E RUN_MINIMAP2_XSEAL_E2E PART2_PMU PART1_PERF
  export RUN_PART1 PART1_THREADS_LIST
}

check_governor() {
  [[ "${GOVERNOR_CHECK:-0}" == "1" ]] || return 0
  if [[ -f /sys/devices/system/cpu/cpu0/cpufreq/scaling_governor ]]; then
    g=$(cat /sys/devices/system/cpu/cpu0/cpufreq/scaling_governor 2>/dev/null || echo unknown)
    if [[ "$g" != "performance" ]]; then
      echo "[experiments] WARNING: cpufreq governor is '$g' (not performance). Set GOVERNOR_CHECK=0 to silence." >&2
    fi
  fi
}

ensure_repo_build_bin() {
  local exe="$REPO_ROOT/$BUILD_DIR/xseal_bench"
  if [[ ! -x "$exe" ]]; then
    echo "Missing $exe — run scripts/10_build.sh first." >&2
    exit 1
  fi
}
