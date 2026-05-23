#!/usr/bin/env bash
# Benchmark minimap2_xseal fused indexer on hg38: disk mmap vs memory-preload.
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
MM_DIR="$(cd "$(dirname "$0")" && pwd)"
BIN="${MINIMAP2_XSEAL_BIN:-$MM_DIR/minimap2_xseal}"
INPUT="${HG38_FASTA:-$ROOT/data/hg38.fa}"
THREADS="${MINIMAP2_THREADS:-12}"
WARMUP="${WARMUP:-1}"
TRIALS="${TRIALS:-3}"
OUT_DIR="${OUT_DIR:-$MM_DIR/bench_hg38_out}"

if [[ ! -f "$INPUT" ]]; then
  echo "hg38 not found: $INPUT" >&2
  exit 1
fi

if [[ ! -x "$BIN" ]]; then
  echo "Building $BIN ..." >&2
  make -C "$MM_DIR" -j"$(nproc)"
fi

mkdir -p "$OUT_DIR"
SUMMARY="$OUT_DIR/summary.txt"
: >"$SUMMARY"

run_mode() {
  local mode="$1"
  local io_env="$2"
  local drop_env="${3:-0}"
  local run_input="${4:-$INPUT}"
  local log="$OUT_DIR/${mode}.log"
  : >"$log"

  echo "=== mode=$mode input=$run_input io=$io_env drop_cache=$drop_env ===" | tee -a "$SUMMARY"
  for ((t = 1; t <= WARMUP + TRIALS; t++)); do
    local idx="$OUT_DIR/.idx_${mode}_${t}.mmi"
    rm -f "$idx"
  done

  for ((t = 1; t <= WARMUP + TRIALS; t++)); do
    local idx="$OUT_DIR/.idx_${mode}_${t}.mmi"
    rm -f "$idx"
    local label="trial"
    [[ $t -le $WARMUP ]] && label="warmup"
    echo "[$mode] $label $t/$((WARMUP + TRIALS))" >&2
    if ! XSEAL_IO_MODE="$io_env" XSEAL_DROP_PAGE_CACHE="$drop_env" \
      "$BIN" -t "$THREADS" -d "$idx" "$run_input" >>"$log" 2>&1; then
      echo "ERROR: $mode trial $t failed (see $log)" >&2
      return 1
    fi
    rm -f "$idx"
  done

  python3 - "$log" "$mode" "$TRIALS" >>"$SUMMARY" <<'PY'
import re, sys, statistics as st
log_path, mode, trials_s = sys.argv[1], sys.argv[2], int(sys.argv[3])
text = open(log_path, encoding="utf-8", errors="replace").read()
pat = re.compile(r"^RESULT\|([^|]+)\|mean\|([\d.]+)\|", re.M)
runs = []
cur = {}
for line in text.splitlines():
    m = pat.match(line)
    if not m:
        continue
    cur[m.group(1)] = float(m.group(2))
    if "minimap2_xseal_total_s" in cur and len(cur) >= 4:
        runs.append(dict(cur))
        cur = {}
if len(runs) <= trials_s:
    print(f"{mode}: insufficient RESULT blocks ({len(runs)})")
    sys.exit(0)
runs = runs[-trials_s:]
keys = ["minimap2_xseal_frontend_s", "minimap2_xseal_frontend_gbp_s",
        "minimap2_xseal_backend_s", "minimap2_xseal_total_s", "minimap2_xseal_total_gbp_s"]
for k in keys:
    vals = [r[k] for r in runs if k in r]
    if not vals:
        continue
    unit = "Gbp/s" if k.endswith("_gbp_s") else "s"
    mu = st.mean(vals)
    sd = st.stdev(vals) if len(vals) > 1 else 0.0
    print(f"  {k}: {mu:.4f} ± {sd:.4f} {unit} (n={len(vals)})")
PY
  echo "" | tee -a "$SUMMARY"
}

# Direct disk I/O: mmap from SSD path; evict page cache before frontend pass
run_mode "disk_io" "mmap" "1" "$INPUT"

# Memory-resident: FASTA on tmpfs (/dev/shm), mmap without forced re-read from block device
SHM_INPUT="/dev/shm/hg38_xseal_bench.fa"
if [[ ! -f "$SHM_INPUT" ]] || [[ "$INPUT" -nt "$SHM_INPUT" ]]; then
  echo "Copying hg38 to $SHM_INPUT for memory-resident runs ..." >&2
  cp -f "$INPUT" "$SHM_INPUT"
fi
run_mode "memory_resident" "mmap" "0" "$SHM_INPUT"

echo "Full logs: $OUT_DIR/{disk_io,memory_resident}.log"
echo "Summary: $SUMMARY"
cat "$SUMMARY"
