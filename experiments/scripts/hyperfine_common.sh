#!/usr/bin/env bash

# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2026 hanrai. All Rights Reserved.

# Optional hyperfine helpers (sourced by e2e / README examples).
# No-op if hyperfine is not installed unless USE_HYPERFINE=1 (then warn once).

paper_repro_has_hyperfine() {
  command -v hyperfine >/dev/null 2>&1
}

# Run hyperfine with JSON export; caller must parse JSON (see 20_e2e_xseal_hg38.sh).
paper_repro_hyperfine_json() {
  local out_json=$1
  shift
  if ! paper_repro_has_hyperfine; then
    [[ "${USE_HYPERFINE:-0}" == "1" ]] && echo "[experiments] USE_HYPERFINE=1 but hyperfine not in PATH" >&2
    return 1
  fi
  local w="${WARMUP:-3}"
  local r="${TRIALS:-10}"
  hyperfine --warmup "$w" --runs "$r" --show-output --export-json "$out_json" "$@"
}
