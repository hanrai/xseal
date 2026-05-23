#!/usr/bin/env bash

# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2026 hanrai. All Rights Reserved.

# run_chunk_sweep — parameter sweep: xseal_bench chunk/batch MiB vs Raw File Throughput (GB/s).
#
# Uses the CMake-built xseal_bench (benchmark/xseal_bench.cpp), NOT rigorous_bench/xseal_bench.cpp.
#
# Usage:
#   benchmark/run_chunk_sweep.sh [PATH_TO_xseal_bench]
# Env:
#   XSEAL_BENCH_EXE     — override bench binary (else first argv)
#   XSEAL_BENCH_INPUT   — FASTA path; if unset, use <repo>/data/hg38.fa when present, else ~32MiB temp synthetic FASTA
#   XSEAL_BENCH_THREADS — worker threads (default: nproc-1, min 1; leaves one logical CPU for reader/OS)
#   XSEAL_BENCH_READER  — mmap | uring | pipe (default: mmap)
#   XSEAL_BENCH_EXTRA_ARGS — extra CLI tokens for each run (e.g. --profile); word-split on $IFS
#   CHUNK_SWEEP_MAX_I    — max grid exponent i (default 10 → 8 MiB); use 9 to stop at 4 MiB
#
# Grid includes 1/512 and 1/256 MiB below the former 1/128 minimum; xseal_bench clamps -B to
# at least one 4KiB page, so sub-4KiB requests may map to the same effective I/O size as 1/256 MiB.

set -euo pipefail

RUNS_PER_POINT="${RUNS_PER_POINT:-10}"
SYNTH_BYTES="${SYNTH_BYTES:-$((32 * 1024 * 1024))}"
CHUNK_SWEEP_MAX_I="${CHUNK_SWEEP_MAX_I:-10}"

die() { echo "run_chunk_sweep: $*" >&2; exit 1; }

ncpus() {
  if command -v nproc >/dev/null 2>&1; then
    nproc
  elif [[ -r /proc/cpuinfo ]]; then
    grep -c '^processor' /proc/cpuinfo
  else
    echo 4
  fi
}

default_worker_threads() {
  local n
  n=$(ncpus)
  if [[ "$n" -gt 1 ]]; then
    echo $((n - 1))
  else
    echo 1
  fi
}

BENCH="${XSEAL_BENCH_EXE:-${1:-}}"
[[ -n "$BENCH" ]] || die "set XSEAL_BENCH_EXE or pass path to xseal_bench as argv1"
[[ -x "$BENCH" ]] || die "not executable: $BENCH"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
DEFAULT_HG38="$REPO_ROOT/data/hg38.fa"

INPUT="${XSEAL_BENCH_INPUT:-}"
TMPFASTA=""
if [[ -z "$INPUT" ]]; then
  if [[ -f "$DEFAULT_HG38" ]]; then
    INPUT=$DEFAULT_HG38
  else
    TMPFASTA=$(mktemp "${TMPDIR:-/tmp}/xseal_run_chunk_sweep_XXXXXX.fa")
    trap '[[ -n "$TMPFASTA" ]] && rm -f "$TMPFASTA"' EXIT
    set +o pipefail
    {
      echo '>synth'
      yes ACGT | tr -d '\n' | head -c "$SYNTH_BYTES"
      echo
    } >"$TMPFASTA"
    set -o pipefail
    INPUT=$TMPFASTA
  fi
elif [[ ! -f "$INPUT" ]]; then
  die "XSEAL_BENCH_INPUT is not a file: $INPUT"
fi

THREADS="${XSEAL_BENCH_THREADS:-$(default_worker_threads)}"
[[ "$THREADS" =~ ^[0-9]+$ ]] && [[ "$THREADS" -ge 1 ]] || THREADS=1

READER="${XSEAL_BENCH_READER:-mmap}"
case "$READER" in
mmap | uring | pipe) ;;
*)
  die "XSEAL_BENCH_READER must be mmap|uring|pipe, got: $READER"
  ;;
esac

# shellcheck disable=SC2086 # intentional word-split for user-supplied flags
EXTRA_ARGS=${XSEAL_BENCH_EXTRA_ARGS-}

extract_raw_gbs() {
  sed -n '/Raw File Throughput/s/.*: *\([0-9.eE+-]*\) *GB\/s.*/\1/p' | head -1
}

one_run_gbs() {
  local mib="$1"
  local out v
  out=$("$BENCH" "$INPUT" -r "$READER" -t "$THREADS" -B "$mib" $EXTRA_ARGS 2>&1) || return 1
  v=$(echo "$out" | extract_raw_gbs)
  [[ -n "$v" ]] || { echo "$out" >&2; return 1; }
  echo "$v"
}

mean_stddev() {
  awk '
    { x[NR] = $1; s += $1; ss += $1 * $1 }
    END {
      n = NR
      if (n < 1) exit 1
      m = s / n
      if (n < 2) { printf "%.17g 0\n", m; exit }
      printf "%.17g %.17g\n", m, sqrt((ss - n * m * m) / (n - 1))
    }'
}

echo "run_chunk_sweep: bench=$BENCH input=$INPUT threads=$THREADS runs_per_mib=$RUNS_PER_POINT reader=$READER extra_args=${EXTRA_ARGS:-<none>}"
echo "chunk_batch_mib mean_raw_gbs stddev_raw_gbs"

# MiB grid: (1/128) * 2^i for i = -2 .. CHUNK_SWEEP_MAX_I (default 10 => 8 MiB)
for i in $(seq -2 "$CHUNK_SWEEP_MAX_I"); do
  mib=$(
    awk -v i="$i" 'BEGIN {
      m = 1.0 / 128.0
      if (i >= 0) { for (j = 0; j < i; j++) m *= 2.0 }
      else { for (j = 0; j < -i; j++) m /= 2.0 }
      printf "%.17g\n", m
    }'
  )
  samples=()
  for ((r = 0; r < RUNS_PER_POINT; r++)); do
    g=$(one_run_gbs "$mib") || die "run failed mib=$mib trial=$r"
    samples+=("$g")
  done
  read -r mean stddev < <(printf '%s\n' "${samples[@]}" | mean_stddev)
  echo "$mib $mean $stddev"
done
