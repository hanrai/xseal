#!/usr/bin/env bash

# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2026 hanrai. All Rights Reserved.

# Build paper hot-RAM benches from cache-resident sources (stable names for repro scripts).
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/../../.." && pwd)"
XSEAL_SRC="$ROOT/experiments/rigorous_bench/xseal_cache_resident_scan_bench.cpp"
XSEAL_BIN="${XSEAL_BIN:-$ROOT/experiments/rigorous_bench/xseal_hot_bench}"
SIMD_CRATE="$ROOT/experiments/rigorous_bench/simdmin_cache_resident_scan_bench"
SIMD_BIN="${SIMD_BIN:-$ROOT/experiments/rigorous_bench/simd_min_hot_bench/target/release/simd_min_hot_bench}"

XSEAL_EXTRA_FLAGS=()
if [[ "${BUILD_HOT_WITH_G:-0}" == "1" ]]; then
  XSEAL_EXTRA_FLAGS+=(-g)
fi

echo "[build] xseal_hot_bench <- xseal_cache_resident_scan_bench.cpp" >&2
clang++ -O3 -DNDEBUG -mavx2 -mbmi2 -mfma -ffast-math -march=native -std=c++20 -pthread \
  "${XSEAL_EXTRA_FLAGS[@]}" \
  -I"$ROOT/include" \
  "$XSEAL_SRC" -o "$XSEAL_BIN"

echo "[build] simd_min_hot_bench <- simdmin_cache_resident_scan_bench" >&2
(
  cd "$SIMD_CRATE"
  env -u CARGO_TARGET_DIR RUSTFLAGS="-C target-cpu=native" cargo build --release -q
)
mkdir -p "$(dirname "$SIMD_BIN")"
cp -f "$SIMD_CRATE/target/release/simdmin_cache_resident_scan_bench" "$SIMD_BIN"

echo "Built: $XSEAL_BIN" >&2
echo "Built: $SIMD_BIN" >&2

build_hot_hash_benches() {
  local xseal_src="$ROOT/experiments/rigorous_bench/xseal_cache_resident_hash_bench.cpp"
  local xseal_bin="${1:-$ROOT/experiments/rigorous_bench/xseal_hot_hash_bench}"
  local simd_crate="$ROOT/experiments/rigorous_bench/simdmin_cache_resident_hash_bench"
  local simd_bin="${2:-$ROOT/experiments/rigorous_bench/simd_min_hot_hash_bench/target/release/simd_min_hot_hash_bench}"
  clang++ -O3 -DNDEBUG -mavx2 -mbmi2 -mfma -ffast-math -march=native -std=c++20 -pthread \
    "${XSEAL_EXTRA_FLAGS[@]}" -I"$ROOT/include" "$xseal_src" -o "$xseal_bin"
  (
    cd "$simd_crate"
    env -u CARGO_TARGET_DIR RUSTFLAGS="-C target-cpu=native" cargo build --release -q
  )
  mkdir -p "$(dirname "$simd_bin")"
  cp -f "$simd_crate/target/release/simdmin_cache_resident_hash_bench" "$simd_bin"
  echo "Built: $xseal_bin" >&2
  echo "Built: $simd_bin" >&2
}
