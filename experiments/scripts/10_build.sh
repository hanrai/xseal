#!/usr/bin/env bash

# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2026 hanrai. All Rights Reserved.

set -euo pipefail
# shellcheck source=experiments/scripts/common.sh
source "$(cd "$(dirname "$0")" && pwd)/common.sh"
load_config
check_governor

mkdir -p "$REPO_ROOT/build/bin"
cd "$REPO_ROOT"

if [[ ! -d "$BUILD_DIR" ]]; then
  echo "Creating $REPO_ROOT/$BUILD_DIR (Release)..."
  cmake -S "$REPO_ROOT" -B "$BUILD_DIR" -DCMAKE_BUILD_TYPE=Release -G Ninja 2>/dev/null \
    || cmake -S "$REPO_ROOT" -B "$BUILD_DIR" -DCMAKE_BUILD_TYPE=Release
fi

cmake --build "$BUILD_DIR" -j"$(nproc)" --target xseal_bench

clang++ -O3 -march=native -std=c++20 \
  -I"$REPO_ROOT/include" \
  "$EXPERIMENTS_ROOT/src/pure_syncmer_micro.cpp" \
  -o "$REPO_ROOT/build/bin/pure_syncmer_micro" \
  -pthread

echo "Built: $REPO_ROOT/$BUILD_DIR/xseal_bench"
echo "Built: $REPO_ROOT/build/bin/pure_syncmer_micro"
