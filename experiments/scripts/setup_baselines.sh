#!/usr/bin/env bash

# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2026 hanrai. All Rights Reserved.

set -euo pipefail

# setup_baselines.sh — unified baseline dependency setup script
# Clones and compiles third-party baseline dependencies into root-level git-ignored 'external/' directory

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
REPO_ROOT=$(cd "$SCRIPT_DIR/../.." && pwd)
EXTERNAL_DIR="$REPO_ROOT/external"

SIMD_MINIMIZERS_URL="https://github.com/rust-seq/simd-minimizers.git"
SIMD_MINIMIZERS_REF="v2.3.1"
MINIMAP2_URL="https://github.com/lh3/minimap2.git"
MINIMAP2_REF="v2.28"

mkdir -p "$EXTERNAL_DIR"

echo "========================================================================="
echo " Setting up Cleanroom Baseline Dependencies under: external/"
echo "========================================================================="

# 1. Setup simd-minimizers
SIMD_MIN_DIR="$EXTERNAL_DIR/simd-minimizers"
if [ -d "$SIMD_MIN_DIR" ]; then
    echo "[vendor] simd-minimizers is already present in $SIMD_MIN_DIR"
else
    echo "[vendor] Cloning simd-minimizers at pinned ref: $SIMD_MINIMIZERS_REF..."
    git clone --depth 1 --branch "$SIMD_MINIMIZERS_REF" "$SIMD_MINIMIZERS_URL" "$SIMD_MIN_DIR"
fi

# 2. Setup minimap2
MINIMAP2_DIR="$EXTERNAL_DIR/minimap2"
if [ -d "$MINIMAP2_DIR" ]; then
    echo "[vendor] minimap2 is already present in $MINIMAP2_DIR"
else
    echo "[vendor] Cloning minimap2 at pinned ref: $MINIMAP2_REF..."
    git clone --depth 1 --branch "$MINIMAP2_REF" "$MINIMAP2_URL" "$MINIMAP2_DIR"
fi

# Build minimap2
echo "[vendor] Compiling minimap2 in release mode..."
cd "$MINIMAP2_DIR"
make -j"$(nproc)"

# Build the sketch bench tool if it exists
if [ -f mm_sketch_bench.c ]; then
    echo "[vendor] Building mm_sketch_bench tool..."
    cc -O3 -I. mm_sketch_bench.c -o mm_sketch_bench -L. -lminimap2 -lz -lpthread -lm
fi

# Symlink internal Rust baselines into external/ for compatibility
echo "[vendor] Symlinking internal Rust baselines to external/..."
ln -sfns "$REPO_ROOT/experiments/baselines/rust_sota_bench" "$EXTERNAL_DIR/rust_sota_bench"

echo "========================================================================="
echo " Cleanroom Baseline Dependencies Setup Complete!"
echo "========================================================================="
