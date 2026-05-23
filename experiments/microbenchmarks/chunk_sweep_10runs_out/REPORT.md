# Chunk sweep report — reader `mmap`

## Summary

- **Reader**: `mmap` (`xseal_bench -r mmap`)
- **Input**: `/home/hanrai/Projects/xseal/data/hg38.fa`
- **Threads**: `12`
- **Trials per MiB point**: `10` (inner); **outer sweeps**: `10`
- **Metric**: Raw File Throughput (GB/s, decimal `1e9` bytes/s)
- **Grid span**: `0.00195312` … `8` MiB (13 points)
- **Peak mean throughput**: **6.2967 GB/s** at `chunk_batch_mib=0.25`
- **Reference lines (plot)**: L1d/L2 full capacity (solid); **4× cache capacity** (dashed) at 0.125 / 2 MiB (from sysfs L1d=32 KiB, L2=512 KiB)

Batch log: `sweep_batch.log` (500 bytes).

## Results (mean ± std across outer sweeps)

| chunk_batch_mib | mean_raw_gbs | std_across_sweeps |
|-----------------|-------------:|------------------:|
| 0.00195312 | 4.5861 | 0.1988 |
| 0.00390625 | 4.7472 | 0.1653 |
| 0.0078125 | 5.6712 | 0.1991 |
| 0.015625 | 6.0200 | 0.1066 |
| 0.03125 | 6.1654 | 0.1169 |
| 0.0625 | 6.1565 | 0.1104 |
| 0.125 | 6.2202 | 0.1171 |
| 0.25 | 6.2967 | 0.0726 |
| 0.5 | 6.1970 | 0.1605 |
| 1 | 6.1310 | 0.0951 |
| 2 | 5.5575 | 0.1888 |
| 4 | 4.8459 | 0.1038 |
| 8 | 4.7480 | 0.0947 |

## Artifacts

- `run_01.txt` … `run_10.txt` — per-sweep tables
- `aggregate_mean_10runs.csv` — table above
- `chunk_sweep_avg_10runs.png` / `.pdf` — plot

## Reproduce

```bash
XSEAL_BENCH_READER=mmap benchmark/run_chunk_sweep.sh build/xseal_bench
python3 benchmark/plot_chunk_sweep_avg_chart.py --out-dir chunk_sweep_10runs_out --reader-label mmap --write-csv --write-report
```
