#!/usr/bin/env bash

# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2026 hanrai. All Rights Reserved.

# 2bit mmap + xseal_harvest_kmers<true> (xseal/xseal_kmer_harvester.hpp); warmup/runs via hyperfine for wall time.
#
# Defaults: 512 KiB synthetic 2bit (urandom) at data/harvest_synth_512k.bin, N=524288 (512k positions).
# Override: BIN_FILE=... N=...  e.g. full hg38: BIN_FILE=data/hg38_2bit.bin N=50000000
# Example:
#   K=31 WARMUP=3 RUNS=10 ./experiments/rigorous_bench/scripts/run_harvest_hyperfine.sh
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
cd "${ROOT}"

SYNTH_512K="${ROOT}/data/harvest_synth_512k.bin"
BIN_FILE="${BIN_FILE:-$SYNTH_512K}"
K="${K:-31}"
N="${N:-524288}"
WARMUP="${WARMUP:-3}"
RUNS="${RUNS:-10}"
HARVEST_BIN="${HARVEST_BIN:-${ROOT}/experiments/rigorous_bench/harvest_bench}"

if [[ ! -f "${BIN_FILE}" ]]; then
  if [[ "${BIN_FILE}" == "${SYNTH_512K}" ]]; then
    mkdir -p "$(dirname "${SYNTH_512K}")"
    dd if=/dev/urandom of="${SYNTH_512K}" bs=1024 count=512 status=none
  else
    echo "error: 2bit file not found: ${BIN_FILE}" >&2
    exit 1
  fi
fi
if [[ ! -x "${HARVEST_BIN}" ]]; then
  echo "error: harvest_bench not executable: ${HARVEST_BIN}" >&2
  echo "build: clang++ -O3 -DNDEBUG -mavx2 -mbmi2 -mfma -ffast-math -march=native -std=c++20 -pthread -Iinclude experiments/rigorous_bench/harvest_bench.cpp -o experiments/rigorous_bench/harvest_bench" >&2
  exit 1
fi
if ! command -v hyperfine >/dev/null 2>&1; then
  echo "error: hyperfine not in PATH" >&2
  exit 1
fi

THREADS_LIST=(1 6 12 16)
for t in "${THREADS_LIST[@]}"; do
  echo "=== hyperfine threads=${t} K=${K} N=${N} warmup=${WARMUP} runs=${RUNS} ===" >&2
  hyperfine --warmup "${WARMUP}" --runs "${RUNS}" --show-output \
    "${HARVEST_BIN} ${BIN_FILE} ${K} ${N} ${t}"
done
