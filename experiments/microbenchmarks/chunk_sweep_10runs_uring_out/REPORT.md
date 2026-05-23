# Chunk sweep report — reader `uring`

## Summary

- **Reader**: `uring` (`xseal_bench -r uring`)
- **Input**: `/home/hanrai/Projects/xseal/data/hg38.fa`
- **Threads**: `12`
- **Trials per MiB point**: `10` (inner); **outer sweeps**: `10`
- **Metric**: Raw File Throughput (GB/s, decimal `1e9` bytes/s)
- **Grid span**: `0.00195312` … `4` MiB (12 points)
- **Peak mean throughput**: **5.4710 GB/s** at `chunk_batch_mib=0.5`
- **Reference lines (plot)**: L1d/L2 full capacity (solid); **4× cache capacity** (dashed) at 0.125 / 2 MiB (from sysfs L1d=32 KiB, L2=512 KiB)

Batch log: `sweep_batch.log` (885 bytes).

### Grid note (vs mmap experiment)

- The **mmap** run in [`../chunk_sweep_10runs_out/`](../chunk_sweep_10runs_out/) used the full grid through **8 MiB** (`CHUNK_SWEEP_MAX_I=10`, 13 points).
- This **uring** batch used `CHUNK_SWEEP_MAX_I=9` (stops at **4 MiB**, 12 points) because a single probe at **8 MiB** aborts:

  ```text
  ./build/xseal_bench data/hg38.fa -r uring -t 12 -B 8
  → std::runtime_error: io_uring fatal error: -14 (EFAULT)
  ```

- Outer sweep wall time was ~**252 s** per pass (vs ~**77 s** for mmap on the same machine/input).

## Results (mean ± std across outer sweeps)

| chunk_batch_mib | mean_raw_gbs | std_across_sweeps |
|-----------------|-------------:|------------------:|
| 0.00195312 | 0.5817 | 0.0050 |
| 0.00390625 | 0.5864 | 0.0099 |
| 0.0078125 | 1.1842 | 0.0195 |
| 0.015625 | 2.1786 | 0.0529 |
| 0.03125 | 3.3585 | 0.0950 |
| 0.0625 | 4.3582 | 0.0762 |
| 0.125 | 5.1503 | 0.0903 |
| 0.25 | 5.2648 | 0.0939 |
| 0.5 | 5.4710 | 0.1703 |
| 1 | 5.3220 | 0.1183 |
| 2 | 5.0223 | 0.1208 |
| 4 | 5.2424 | 0.0813 |

## Artifacts

- `run_01.txt` … `run_10.txt` — per-sweep tables
- `aggregate_mean_10runs.csv` — table above
- `chunk_sweep_avg_10runs.png` / `.pdf` — plot

## Reproduce

```bash
XSEAL_BENCH_READER=uring CHUNK_SWEEP_MAX_I=9 benchmark/run_chunk_sweep.sh build/xseal_bench
python3 benchmark/plot_chunk_sweep_avg_chart.py --out-dir benchmark/chunk_sweep_10runs_uring_out \
  --reader-label uring --write-csv --write-report
```

## Comparison with mmap (same input, 12 threads)

| | mmap (`chunk_sweep_10runs_out`) | uring (this run) |
|--|--------------------------------|------------------|
| Peak mean Raw GB/s | ~5.94 @ 0.25 MiB | ~5.47 @ 0.5 MiB |
| Grid max | 8 MiB | 4 MiB (8 MiB crashes) |
| ~Wall / outer sweep | ~77 s | ~252 s |
