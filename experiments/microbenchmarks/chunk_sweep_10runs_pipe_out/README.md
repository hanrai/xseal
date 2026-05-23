# Chunk sweep — pipe reader

**pipe** (`-r pipe`, file opened via `open()` on `data/hg38.fa`) experiment: 10 outer sweeps, 2026-05-16.

- **Report**: [`REPORT.md`](REPORT.md)
- **Plot**: [`chunk_sweep_avg_10runs_pipe.png`](chunk_sweep_avg_10runs_pipe.png)
- **Batch log**: [`sweep_batch.log`](sweep_batch.log)

Sibling runs: [`../chunk_sweep_10runs_out/`](../chunk_sweep_10runs_out/) (mmap), [`../chunk_sweep_10runs_uring_out/`](../chunk_sweep_10runs_uring_out/) (uring).

**Note:** `XsealPipeReader::get_file_size()` was fixed to return bytes delivered (`total_offset`) so Raw File Throughput is meaningful (rebuild `xseal_bench` required).
