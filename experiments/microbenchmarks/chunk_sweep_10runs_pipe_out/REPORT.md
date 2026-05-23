# Chunk sweep report — reader `pipe`

## Summary

- **Reader**: `pipe` (`xseal_bench -r pipe`)
- **Input**: `/home/hanrai/Projects/xseal/data/hg38.fa`
- **Threads**: `12`
- **Trials per MiB point**: `10` (inner); **outer sweeps**: `10`
- **Metric**: Raw File Throughput (GB/s, decimal `1e9` bytes/s)
- **Grid span**: `0.00195312` … `8` MiB (13 points)
- **Peak mean throughput**: **5.1773 GB/s** at `chunk_batch_mib=0.5`
- **Reference lines (plot)**: L1d/L2 full capacity (solid); **4× cache capacity** (dashed) at 0.125 / 2 MiB (from sysfs L1d=32 KiB, L2=512 KiB)

Batch log: `sweep_batch.log` (795 bytes).

### Method note

- Full grid **1/512 … 8 MiB** (`CHUNK_SWEEP_MAX_I=10`, 13 points).
- Requires `xseal_bench` built after fixing `get_file_size()` on pipe (uses `total_offset`, not `0`).
- ~**194 s** wall time per outer sweep on this machine.

## Results (mean ± std across outer sweeps)

| chunk_batch_mib | mean_raw_gbs | std_across_sweeps |
|-----------------|-------------:|------------------:|
| 0.00195312 | 1.6966 | 0.0438 |
| 0.00390625 | 1.6909 | 0.0364 |
| 0.0078125 | 2.7969 | 0.0843 |
| 0.015625 | 4.2000 | 0.0691 |
| 0.03125 | 4.9969 | 0.1381 |
| 0.0625 | 4.9531 | 0.1218 |
| 0.125 | 5.1610 | 0.1169 |
| 0.25 | 5.1184 | 0.1571 |
| 0.5 | 5.1773 | 0.0951 |
| 1 | 4.3796 | 0.1321 |
| 2 | 3.7657 | 0.0956 |
| 4 | 3.5782 | 0.0426 |
| 8 | 3.4692 | 0.0554 |

## Artifacts

- `run_01.txt` … `run_10.txt` — per-sweep tables
- `aggregate_mean_10runs.csv` — table above
- `chunk_sweep_avg_10runs.png` / `.pdf` — plot

## Reproduce

```bash
XSEAL_BENCH_READER=pipe benchmark/run_chunk_sweep.sh build/xseal_bench
python3 benchmark/plot_chunk_sweep_avg_chart.py --out-dir chunk_sweep_10runs_pipe_out --reader-label pipe --write-csv --write-report
```
