# Chunk sweep — io_uring reader

**io_uring** (`-r uring`) experiment: 10 outer sweeps on `data/hg38.fa`, 2026-05-16.

- **Report**: [`REPORT.md`](REPORT.md)
- **Plot (uring, distinct name)**: [`chunk_sweep_avg_10runs_uring.png`](chunk_sweep_avg_10runs_uring.png)
- Legacy name in this folder: `chunk_sweep_avg_10runs.png` (same data; mmap chart lives under `../chunk_sweep_10runs_out/`)
- **Batch log**: [`sweep_batch.log`](sweep_batch.log)

The archived **mmap** run is in [`../chunk_sweep_10runs_out/`](../chunk_sweep_10runs_out/).

**Note:** `-B 8` (8 MiB) is excluded (`CHUNK_SWEEP_MAX_I=9`) due to `io_uring` EFAULT; see REPORT.
