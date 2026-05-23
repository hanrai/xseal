#!/usr/bin/env bash

# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2026 hanrai. All Rights Reserved.

set -euo pipefail

ROOT="$(cd "$(dirname "$0")/../../.." && pwd)"
BIN_FILE="${BIN_FILE:-$ROOT/data/hg38_2bit.bin}"
MAX_BASES="${MAX_BASES:-50000000}"
SYNC_K="${SYNC_K:-31}"
S="${S:-11}"
MIN_W="${MIN_W:-21}"

AUDIT_BIN="$ROOT/experiments/rigorous_bench/marker_audit/target/release/marker_audit"

echo "[build] marker_audit" >&2
(
  cd "$ROOT/experiments/rigorous_bench/marker_audit"
  env -u CARGO_TARGET_DIR RUSTFLAGS="-C target-cpu=native" cargo build --release -q
)

echo "=== simd-minimizers semantics audit (first ${MAX_BASES} bases) ==="
"$AUDIT_BIN" "$BIN_FILE" --max-bases "$MAX_BASES" --sync-k "$SYNC_K" --s "$S" --min-w "$MIN_W"
