#!/usr/bin/env bash

# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2026 hanrai. All Rights Reserved.

set -euo pipefail
# shellcheck source=experiments/scripts/common.sh
source "$(cd "$(dirname "$0")" && pwd)/common.sh"
load_config
check_governor
ensure_repo_build_bin

mkdir -p "$EXPERIMENTS_ROOT/results"
BIN="$REPO_ROOT/$BUILD_DIR/xseal_bench"
INPUT="${HG38_FASTA:-$REPO_ROOT/data/hg38.fa}"

if [[ ! -f "$INPUT" ]]; then
  echo "Input FASTA not found: $INPUT" >&2
  exit 1
fi

if [[ "${USE_SHM}" == "1" ]]; then
  SHM_PATH=/dev/shm/hg38.fa
  echo "Copying to $SHM_PATH ..."
  cp -f "$INPUT" "$SHM_PATH"
  INPUT="$SHM_PATH"
fi

RAW="$EXPERIMENTS_ROOT/results/e2e_xseal_raw.txt"
HF_JSON="$EXPERIMENTS_ROOT/results/hyperfine_e2e.json"
: >"$RAW"

for ((i = 1; i <= WARMUP; i++)); do
  echo "[e2e] warmup $i/$WARMUP" >&2
  "$BIN" "$INPUT" -r mmap -m syncmer -t "$XSEAL_E2E_THREADS" >/dev/null 2>&1 || true
done

for ((i = 1; i <= TRIALS; i++)); do
  echo "[e2e] trial $i/$TRIALS" >&2
  "$BIN" "$INPUT" -r mmap -m syncmer -t "$XSEAL_E2E_THREADS" | tee -a "$RAW"
done

if [[ "${USE_HYPERFINE:-0}" == "1" ]] && command -v hyperfine >/dev/null 2>&1; then
  echo "[e2e] hyperfine (wall clock)" >&2
  hf_cmd=( "$BIN" "$INPUT" -r mmap -m syncmer -t "$XSEAL_E2E_THREADS" )
  hyperfine --warmup "$WARMUP" --runs "$TRIALS" --export-json "$HF_JSON" "${hf_cmd[@]}" || true
else
  [[ "${USE_HYPERFINE:-0}" == "1" ]] && echo "[e2e] USE_HYPERFINE=1 but hyperfine not installed" >&2
  rm -f "$HF_JSON"
fi

if [[ "${E2E_PROFILE_ONCE}" == "1" ]]; then
  echo "[e2e] one --profile run" >&2
  "$BIN" "$INPUT" -r mmap -m syncmer -t "$XSEAL_E2E_THREADS" --profile \
    >"$EXPERIMENTS_ROOT/results/e2e_profile_once.txt" 2>&1 || true
fi

python3 <<PY
import json, re, statistics as st
from pathlib import Path
raw_path = "$RAW"
out_path = "$EXPERIMENTS_ROOT/results/metrics_e2e.json"
hf_path = Path("$HF_JSON")
text = open(raw_path, "r", errors="ignore").read().splitlines()
gbp, tot = [], []
for line in text:
    if "Base Throughput" in line:
        m = re.search(r":\s*([\d.]+)\s*Gbp/s", line)
        if m:
            gbp.append(float(m.group(1)))
    if "Total Time" in line:
        m = re.search(r":\s*([\d.]+)\s*s\s*$", line)
        if m:
            tot.append(float(m.group(1)))
n = min(len(gbp), len(tot))
if n == 0:
    open(out_path, "w").write(json.dumps({"error": "no_parse", "raw": raw_path}, indent=2))
    raise SystemExit(0)
gbp, tot = gbp[-n:], tot[-n:]
obj = {
    "e2e_xseal_integration_gbp_s": {
        "mean": st.mean(gbp),
        "stdev": st.stdev(gbp) if len(gbp) > 1 else 0.0,
        "n": len(gbp),
        "unit": "Gbp/s",
    },
    "e2e_total_time_s": {"mean": st.mean(tot), "stdev": st.stdev(tot) if len(tot) > 1 else 0.0, "n": len(tot)},
}
if hf_path.is_file():
    try:
        hf = json.loads(hf_path.read_text(encoding="utf-8"))
        rs = hf.get("results") or []
        if rs and "mean" in rs[0]:
            obj["e2e_hyperfine_wall_mean_s"] = {
                "mean": float(rs[0]["mean"]),
                "stdev": float(rs[0].get("stddev", 0.0) or 0.0),
                "n": int(rs[0].get("runs", len(rs[0].get("times") or [])) or 1),
                "unit": "s",
            }
    except (json.JSONDecodeError, KeyError, TypeError, ValueError):
        pass
open(out_path, "w").write(json.dumps(obj, indent=2))
print(json.dumps(obj["e2e_xseal_integration_gbp_s"], indent=2))
PY

echo "Wrote $EXPERIMENTS_ROOT/results/metrics_e2e.json"
