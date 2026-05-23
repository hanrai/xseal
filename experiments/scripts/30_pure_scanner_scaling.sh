#!/usr/bin/env bash

# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2026 hanrai. All Rights Reserved.

set -euo pipefail
# shellcheck source=experiments/scripts/common.sh
source "$(cd "$(dirname "$0")" && pwd)/common.sh"
load_config

BIN="$REPO_ROOT/build/bin/pure_syncmer_micro"
if [[ ! -x "$BIN" ]]; then
  echo "Run 10_build.sh first (missing $BIN)" >&2
  exit 1
fi

OUT="$EXPERIMENTS_ROOT/results/pure_scan_scaling.csv"
echo "threads,aggregate_gbp_s,wall_s" >"$OUT"

IFS=',' read -r -a TLIST <<<"${PURE_SCAN_THREADS_LIST}"
for t in "${TLIST[@]}"; do
  t="$(echo "$t" | tr -d ' ')"
  [[ -z "$t" ]] && continue
  echo "[pure_scan] threads=$t" >&2
  run_out=$("$BIN" "$t" 2>&1 || true)
  g=$(echo "$run_out" | grep '^RESULT aggregate_gbp_s' | awk '{print $3}')
  w=$(echo "$run_out" | grep '^RESULT wall_seconds' | awk '{print $3}')
  echo "$t,${g:-nan},${w:-nan}" >>"$OUT"
done

python3 <<PY
import csv, json, statistics as st
rows = []
with open("$OUT") as f:
    r = csv.DictReader(f)
    for row in r:
        try:
            rows.append(float(row["aggregate_gbp_s"]))
        except ValueError:
            pass
# map 12-thread row to paper metric id
by_t = {}
with open("$OUT") as f:
    r = csv.DictReader(f)
    for row in r:
        try:
            by_t[int(row["threads"])] = float(row["aggregate_gbp_s"])
        except (ValueError, KeyError):
            pass
obj = {"pure_scan_by_threads": by_t}
if 12 in by_t:
    v = {"mean": by_t[12], "stdev": 0.0, "n": 1, "unit": "Gbp/s"}
    obj["paper_repro_pure_xseal_open_12t_gbp_s"] = v
    obj["e2e_xseal_mt_scaling_12t_open_gbp_s"] = v  # legacy alias for older reports
open("$EXPERIMENTS_ROOT/results/metrics_pure_scan.json", "w").write(json.dumps(obj, indent=2))
print(json.dumps(obj, indent=2))
PY

echo "Wrote $OUT and results/metrics_pure_scan.json"
