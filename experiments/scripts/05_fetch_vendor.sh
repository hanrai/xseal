#!/usr/bin/env bash

# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2026 hanrai. All Rights Reserved.

set -euo pipefail
# shellcheck source=experiments/scripts/common.sh
source "$(cd "$(dirname "$0")" && pwd)/common.sh"
load_config

if [[ "${FETCH_VENDOR:-1}" != "1" ]]; then
  echo "[vendor] FETCH_VENDOR!=1; skip fetch" >&2
  exit 0
fi

mkdir -p "$EXPERIMENTS_ROOT/results"
SKIP="$EXPERIMENTS_ROOT/results/vendor.skip"
rm -f "$SKIP" "$EXPERIMENTS_ROOT/results/vendor.ok"

VENV="$EXPERIMENTS_ROOT/vendor_versions.env"
if [[ ! -f "$VENV" ]]; then
  echo "Missing $VENV" | tee "$SKIP"
  exit 0
fi
# shellcheck source=/dev/null
set -a
# shellcheck source=/dev/null
source "$VENV"
set +a

ensure_simd_minimizers() {
  local dest="$REPO_ROOT/external/simd-minimizers"
  local url="$SIMD_MINIMIZERS_REPO_URL"
  local ref="$SIMD_MINIMIZERS_GIT_REF"

  if [[ -f "$dest/Cargo.toml" ]]; then
    echo "[vendor] simd-minimizers already present: $dest" >&2
    (cd "$dest" && git rev-parse --short HEAD 2>/dev/null || true)
    return 0
  fi

  echo "[vendor] fetching simd-minimizers -> $dest (ref=$ref)" >&2
  mkdir -p "$(dirname "$dest")"
  rm -rf "$dest"
  if git clone --depth 1 --branch "$ref" "$url" "$dest" 2>/dev/null; then
    :
  else
    rm -rf "$dest"
    git clone "$url" "$dest"
    (cd "$dest" && git checkout -q "$ref") || {
      echo "checkout $ref failed in $dest" | tee "$SKIP"
      rm -rf "$dest"
      return 1
    }
  fi
  [[ -f "$dest/Cargo.toml" ]] || {
    echo "invalid simd-minimizers tree" | tee "$SKIP"
    return 1
  }
  return 0
}

ensure_simd_minimizers || exit 0

echo "[vendor] simd-minimizers OK: $(git -C "$REPO_ROOT/external/simd-minimizers" rev-parse --short HEAD 2>/dev/null || echo unknown)" >&2
echo ok >"$EXPERIMENTS_ROOT/results/vendor.ok"
