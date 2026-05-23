#!/usr/bin/env bash
# hg38: native vs XSeal (k=21,w=11,-t12,11 workers). Sub-metrics for paper frontend claims.
# Safe for repeated runs: one case per invocation segment, cooldown, optional skip shm.
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
MM_DIR="$(cd "$(dirname "$0")" && pwd)"
INPUT="${HG38_FASTA:-$ROOT/data/hg38.fa}"
THREADS="${MINIMAP2_THREADS:-12}"
WORKERS="${XSEAL_WORKERS:-11}"
K="${MM_K:-21}"
W="${MM_W:-11}"
WARMUP="${WARMUP:-0}"
TRIALS="${TRIALS:-1}"
COOLDOWN_SEC="${COOLDOWN_SEC:-45}"
OUT_DIR="${OUT_DIR:-$MM_DIR/bench_hg38_rerun}"
# Space-separated: native_disk xseal_min_disk xseal_syn_disk (default: disk only, 3 cases)
CASES="${CASES:-native_disk xseal_min_disk xseal_syn_disk}"
USE_SHM="${USE_SHM:-0}"

if [[ ! -f "$INPUT" ]]; then
  echo "hg38 not found: $INPUT" >&2
  exit 1
fi

if [[ ! -x "$MM_DIR/minimap2_xseal" ]]; then
  make -C "$MM_DIR" -j"$(nproc)"
fi

BIN_MM="$MM_DIR/minimap2_xseal"
BIN_SYN="$MM_DIR/minimap2_xseal_syncmer"
mkdir -p "$OUT_DIR"
SUMMARY="$OUT_DIR/summary.txt"
: >"$SUMMARY"

mem_ok() {
  local avail_kb
  avail_kb=$(awk '/^MemAvailable:/ {print $2}' /proc/meminfo)
  local swap_used_kb
  swap_used_kb=$(awk '/^SwapTotal:/ {t=$2} /^SwapFree:/ {print t-$2}' /proc/meminfo)
  if [[ "${swap_used_kb:-0}" -gt 6291456 ]]; then
    echo "WARN: swap used >6GiB (${swap_used_kb} KiB); stopping benchmark loop." >&2
    return 1
  fi
  if [[ "${avail_kb:-0}" -lt 4194304 ]]; then
    echo "WARN: MemAvailable <4GiB (${avail_kb} KiB); stopping." >&2
    return 1
  fi
  return 0
}

cooldown() {
  echo "Cooldown ${COOLDOWN_SEC}s (sync)..." >&2
  sync
  sleep "$COOLDOWN_SEC"
}

run_case() {
  local name="$1"
  local bin="$2"
  local extra_env="$3"
  local run_input="${4:-$INPUT}"
  local log="$OUT_DIR/${name}.log"
  : >"$log"

  echo "=== $name (input=$run_input) ===" | tee -a "$SUMMARY"
  for ((t = 1; t <= WARMUP + TRIALS; t++)); do
    mem_ok || return 1
    local idx="$OUT_DIR/.${name}_${t}.mmi"
    rm -f "$idx"
    local label=trial
    [[ $t -le $WARMUP ]] && label=warmup
    echo "[$name] $label $t/$((WARMUP + TRIALS))" >&2
    env $extra_env XSEAL_WORKERS="$WORKERS" \
      "$bin" -t "$THREADS" -k "$K" -w "$W" -d "$idx" "$run_input" >>"$log" 2>&1 || {
        echo "  RUN FAILED (see $log)" | tee -a "$SUMMARY"
        return 1
      }
    rm -f "$idx"
    [[ $t -lt $((WARMUP + TRIALS)) ]] && cooldown
  done

  python3 - "$log" "$name" "$TRIALS" >>"$SUMMARY" <<'PY'
import re, sys, statistics as st
log_path, name, trials_s = sys.argv[1], sys.argv[2], int(sys.argv[3])
text = open(log_path, encoding="utf-8", errors="replace").read()
pat = re.compile(r"^RESULT\|([^|]+)\|mean\|([\d.]+)\|", re.M)
runs, cur = [], {}
for line in text.splitlines():
    m = pat.match(line)
    if m:
        cur[m.group(1)] = float(m.group(2))
    if "[M::main] Real time:" in line and cur:
        runs.append(dict(cur))
        cur = {}
if cur:
    runs.append(dict(cur))
if len(runs) < trials_s:
    print(f"  insufficient trials ({len(runs)})")
    sys.exit(0)
runs = runs[-trials_s:]
keys = sorted({k for r in runs for k in r})
prio_suffix = (
    "_io_s", "_file_map_s", "_parse_s", "_compute_wall_s", "_worker_compute_s",
    "_index_insert_wall_s", "_index_insert_accum_s", "_index_insert_s",
    "_frontend_s", "_compute_gbp_s", "_frontend_gbp_s",
    "_backend_s", "_e2e_s", "_e2e_gbp_s",
    "_distinct_minimizers", "_singleton_pct", "_seed_hits", "_peak_rss_gib",
)
prio = [k for k in keys if any(k.endswith(s) for s in prio_suffix)]
for k in prio:
    vals = [r[k] for r in runs if k in r]
    if not vals:
        continue
    unit = ""
    if k.endswith("_gbp_s"):
        unit = "Gbp/s"
    elif k.endswith("_pct"):
        unit = "%"
    elif k.endswith("_gib"):
        unit = "GiB"
    elif k.endswith("_s"):
        unit = "s"
    mu = st.mean(vals)
    sd = st.stdev(vals) if len(vals) > 1 else 0.0
    print(f"  {k}: {mu:.4f} ± {sd:.4f} {unit} (n={len(vals)})")
print()
PY
  cooldown
}

SHM_INPUT="/dev/shm/hg38_xseal_bench.fa"
if [[ "$USE_SHM" == "1" ]]; then
  if [[ ! -f "$SHM_INPUT" ]] || [[ "$INPUT" -nt "$SHM_INPUT" ]]; then
    echo "Staging hg38 to $SHM_INPUT ..." >&2
    cp -f "$INPUT" "$SHM_INPUT"
  fi
fi

for case in $CASES; do
  mem_ok || break
  case "$case" in
    native_disk)
      run_case "native_disk" "$BIN_MM" "XSEAL_USE_NATIVE=1 XSEAL_IO_MODE=mmap XSEAL_DROP_PAGE_CACHE=1" "$INPUT"
      ;;
    native_ram)
      run_case "native_ram" "$BIN_MM" "XSEAL_USE_NATIVE=1 XSEAL_IO_MODE=mmap XSEAL_DROP_PAGE_CACHE=0" "$SHM_INPUT"
      ;;
    xseal_min_disk)
      run_case "xseal_min_disk" "$BIN_MM" "XSEAL_IO_MODE=mmap XSEAL_DROP_PAGE_CACHE=1" "$INPUT"
      ;;
    xseal_min_ram)
      run_case "xseal_min_ram" "$BIN_MM" "XSEAL_IO_MODE=mmap XSEAL_DROP_PAGE_CACHE=0" "$SHM_INPUT"
      ;;
    xseal_syn_disk)
      run_case "xseal_syn_disk" "$BIN_SYN" "XSEAL_IO_MODE=mmap XSEAL_DROP_PAGE_CACHE=1" "$INPUT"
      ;;
    xseal_syn_ram)
      run_case "xseal_syn_ram" "$BIN_SYN" "XSEAL_IO_MODE=mmap XSEAL_DROP_PAGE_CACHE=0" "$SHM_INPUT"
      ;;
    *)
      echo "unknown case: $case" >&2
      exit 1
      ;;
  esac
done

echo "Summary: $SUMMARY"
cat "$SUMMARY"
