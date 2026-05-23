# XSeal Benchmark Reproduction Guide

This directory contains the complete set of C++ and Rust benchmarks, microbenchmarks, and reproducibility harnesses used in the XSeal manuscript. Every file and crate has been renamed to align perfectly with the microarchitectural experiments and performance findings in the paper.

---

## 0. Reference genome (hg38)

The repository does not ship hg38. Download and place at `data/hg38.fa`:

```bash
./experiments/scripts/download_hg38.sh
cp experiments/config.env.example experiments/config.env
# set HG38_FASTA=$REPO_ROOT/data/hg38.fa if needed
```

Pack 2-bit for Part I hot benches:

```bash
python3 experiments/rigorous_bench/scripts/pack_hg38.py data/hg38.fa data/hg38_2bit.bin
```

Paper Table 1 hot-RAM reproduction uses `experiments/rigorous_bench/scripts/build_hot_benches.sh` (builds `xseal_hot_bench` / `simd_min_hot_bench` from cache-resident sources).

## 1. Setup Baseline Dependencies

Before running any experiments, execute the unified baseline setup script. This script automatically downloads, compiles, and pins the SOTA baselines into the root-level git-ignored `external/` directory:

```bash
./experiments/scripts/setup_baselines.sh
```

This will set up:
- **`simd-minimizers`** (v2.3.1) at `external/simd-minimizers/`
- **`minimap2`** (v2.28) at `external/minimap2/`
- Symlinks for `rust_sota_bench` into `external/` for compatibility

---

## 2. Benchmark Suite Directory (Rigorous Benchmark Alignment)

All rigorous microarchitectural and throughput benchmarks reside in [`experiments/rigorous_bench/`](file:///home/hanrai/Projects/xseal/experiments/rigorous_bench/). Below is the exact mapping of source files to the paper's figures and tables.

### 2.1. Genomic Scanning Throughput & Thread Scaling (Table 1 / Table 2 / Figure 3)
Measures the syncmer and minimizer scanning throughput on packed 2-bit hg38 genomes across 1 to 12 threads.
- **XSeal C++ Scanner**:
  - Source: [`experiments/rigorous_bench/xseal_scan_throughput_bench.cpp`](file:///home/hanrai/Projects/xseal/experiments/rigorous_bench/xseal_scan_throughput_bench.cpp)
  - Compilation:
    ```bash
    clang++ -O3 -march=native -fopenmp=libgomp -Iinclude experiments/rigorous_bench/xseal_scan_throughput_bench.cpp -o build/bin/xseal_scan_throughput_bench
    ```
- **SIMD-Minimizers Rust Scanner**:
  - Source: [`experiments/rigorous_bench/simdmin_scan_throughput_bench`](file:///home/hanrai/Projects/xseal/experiments/rigorous_bench/simdmin_scan_throughput_bench/)
  - Compilation:
    ```bash
    cargo build --manifest-path experiments/rigorous_bench/simdmin_scan_throughput_bench/Cargo.toml --release --target-dir build/rust_target
    cp build/rust_target/release/simdmin_scan_throughput_bench build/bin/
    ```

### 2.2. Microarchitectural Hazards: Store-to-Load Forwarding (STLF) & Window Reduction (Section 5.3 / Figure 4)
Examines the CPU performance penalties of STLF memory hazards and evaluates AVX2 window reduction efficiency.
- **STLF Penalty Microbenchmark**:
  - Source: [`experiments/rigorous_bench/stlf_penalty_microbench.cpp`](file:///home/hanrai/Projects/xseal/experiments/rigorous_bench/stlf_penalty_microbench.cpp)
  - Compilation:
    ```bash
    clang++ -O3 -march=native -std=c++20 experiments/rigorous_bench/stlf_penalty_microbench.cpp -o build/bin/stlf_penalty_microbench
    ```
- **Window Reduction Microbenchmark**:
  - Source: [`experiments/rigorous_bench/window_reduction_microbench.cpp`](file:///home/hanrai/Projects/xseal/experiments/rigorous_bench/window_reduction_microbench.cpp)
  - Compilation:
    ```bash
    clang++ -O3 -march=native -std=c++20 experiments/rigorous_bench/window_reduction_microbench.cpp -o build/bin/window_reduction_microbench
    ```

### 2.3. Base Encoder Performance (Table 4)
Measures the scaling and raw performance of two-bit base encoding algorithms on synthetic fasta buffers.
- **Encoder Scaling Throughput**:
  - Source: [`experiments/rigorous_bench/encoder_throughput_scaling_bench.cpp`](file:///home/hanrai/Projects/xseal/experiments/rigorous_bench/encoder_throughput_scaling_bench.cpp)
  - Compilation:
    ```bash
    clang++ -O3 -march=native -std=c++20 -Iinclude -Iinclude/xseal experiments/rigorous_bench/encoder_throughput_scaling_bench.cpp -o build/bin/encoder_throughput_scaling_bench
    ```
- **Encoder Microbenchmark**:
  - Source: [`experiments/rigorous_bench/encoder_microbench.cpp`](file:///home/hanrai/Projects/xseal/experiments/rigorous_bench/encoder_microbench.cpp)
  - Compilation:
    ```bash
    clang++ -O3 -march=native -std=c++20 -Iinclude experiments/rigorous_bench/encoder_microbench.cpp -o build/bin/encoder_microbench
    ```
- **Encoder Single-Threaded PMU Profiler**:
  - Source: [`experiments/rigorous_bench/encoder_pmu_profiler.cpp`](file:///home/hanrai/Projects/xseal/experiments/rigorous_bench/encoder_pmu_profiler.cpp)
  - Compilation:
    ```bash
    clang++ -O3 -march=native -std=c++20 -Iinclude experiments/rigorous_bench/encoder_pmu_profiler.cpp -o build/bin/encoder_pmu_profiler
    ```
- **Encoder Multi-Threaded PMU Profiler**:
  - Source: [`experiments/rigorous_bench/encoder_multithread_pmu_profiler.cpp`](file:///home/hanrai/Projects/xseal/experiments/rigorous_bench/encoder_multithread_pmu_profiler.cpp)
  - Compilation:
    ```bash
    clang++ -O3 -march=native -std=c++20 -Iinclude experiments/rigorous_bench/encoder_multithread_pmu_profiler.cpp -o build/bin/encoder_multithread_pmu_profiler
    ```

### 2.4. FASTA Parser Throughput & PMU Counters (Table 5)
Measures DNA sequence FASTA parsing performance compared to LH3's reference `kseq.h` library.
- **Parser Throughput Benchmark**:
  - Source: [`experiments/rigorous_bench/parser_throughput_bench.cpp`](file:///home/hanrai/Projects/xseal/experiments/rigorous_bench/parser_throughput_bench.cpp)
  - Compilation:
    ```bash
    clang++ -O3 -march=native -std=c++20 -Iinclude -Iexternal/minimap2 experiments/rigorous_bench/parser_throughput_bench.cpp -o build/bin/parser_throughput_bench -lz
    ```
- **Parser PMU Profiler**:
  - Source: [`experiments/rigorous_bench/parser_pmu_profiler.cpp`](file:///home/hanrai/Projects/xseal/experiments/rigorous_bench/parser_pmu_profiler.cpp)
  - Compilation:
    ```bash
    clang++ -O3 -march=native -std=c++20 -Iinclude -Iexternal/minimap2 experiments/rigorous_bench/parser_pmu_profiler.cpp -o build/bin/parser_pmu_profiler -lz
    ```

### 2.5. Hashing & Cache-Resident Syncmer Loops (Figure 2)
Isolates cache effects by running multiple syncmer iterations on a hot L3/L2 cache-resident sequence.
- **XSeal Cache-Resident Scan**:
  - Source: [`experiments/rigorous_bench/xseal_cache_resident_scan_bench.cpp`](file:///home/hanrai/Projects/xseal/experiments/rigorous_bench/xseal_cache_resident_scan_bench.cpp)
  - Compilation:
    ```bash
    clang++ -O3 -march=native -std=c++20 -Iinclude experiments/rigorous_bench/xseal_cache_resident_scan_bench.cpp -o build/bin/xseal_cache_resident_scan_bench
    ```
- **XSeal Cache-Resident Hashing**:
  - Source: [`experiments/rigorous_bench/xseal_cache_resident_hash_bench.cpp`](file:///home/hanrai/Projects/xseal/experiments/rigorous_bench/xseal_cache_resident_hash_bench.cpp)
  - Compilation:
    ```bash
    clang++ -O3 -march=native -std=c++20 -Iinclude experiments/rigorous_bench/xseal_cache_resident_hash_bench.cpp -o build/bin/xseal_cache_resident_hash_bench
    ```
- **SimdMin Cache-Resident Scan (Rust)**:
  - Source: [`experiments/rigorous_bench/simdmin_cache_resident_scan_bench`](file:///home/hanrai/Projects/xseal/experiments/rigorous_bench/simdmin_cache_resident_scan_bench/)
  - Compilation:
    ```bash
    cargo build --manifest-path experiments/rigorous_bench/simdmin_cache_resident_scan_bench/Cargo.toml --release --target-dir build/rust_target
    cp build/rust_target/release/simdmin_cache_resident_scan_bench build/bin/
    ```
- **SimdMin Cache-Resident Hashing (Rust)**:
  - Source: [`experiments/rigorous_bench/simdmin_cache_resident_hash_bench`](file:///home/hanrai/Projects/xseal/experiments/rigorous_bench/simdmin_cache_resident_hash_bench/)
  - Compilation:
    ```bash
    cargo build --manifest-path experiments/rigorous_bench/simdmin_cache_resident_hash_bench/Cargo.toml --release --target-dir build/rust_target
    cp build/rust_target/release/simdmin_cache_resident_hash_bench build/bin/
    ```

### 2.6. Isolated Hashing Throughput
Isolates the raw processing throughput of the syncmer hashing function without scanning memory.
- **XSeal SIMD Hashing Isolation**:
  - Source: [`experiments/rigorous_bench/xseal_simd_hash_isolation_bench.cpp`](file:///home/hanrai/Projects/xseal/experiments/rigorous_bench/xseal_simd_hash_isolation_bench.cpp)
  - Compilation:
    ```bash
    clang++ -O3 -march=native -std=c++20 -Iinclude experiments/rigorous_bench/xseal_simd_hash_isolation_bench.cpp -o build/bin/xseal_simd_hash_isolation_bench
    ```
- **SimdMin Hashing Isolation (Rust)**:
  - Source: [`experiments/rigorous_bench/simdmin_hash_isolation_bench`](file:///home/hanrai/Projects/xseal/experiments/rigorous_bench/simdmin_hash_isolation_bench/)
  - Compilation:
    ```bash
    cargo build --manifest-path experiments/rigorous_bench/simdmin_hash_isolation_bench/Cargo.toml --release --target-dir build/rust_target
    cp build/rust_target/release/simdmin_hash_isolation_bench build/bin/
    ```

### 2.7. Decoupled K-mer Harvesting Throughput
Measures sequence traversal speed when extracting sequences at syncmer coordinate positions.
- **Decoupled K-mer Harvesting**:
  - Source: [`experiments/rigorous_bench/kmer_harvest_throughput_bench.cpp`](file:///home/hanrai/Projects/xseal/experiments/rigorous_bench/kmer_harvest_throughput_bench.cpp)
  - Compilation:
    ```bash
    clang++ -O3 -march=native -std=c++20 -Iinclude experiments/rigorous_bench/kmer_harvest_throughput_bench.cpp -o build/bin/kmer_harvest_throughput_bench
    ```

### 2.8. Fixed 128 KiB Block Buffer Scaling
Benchmarks scanning throughput limited to a rigid 128 KiB sequence buffer block size.
- **XSeal Fixed Block Scan**:
  - Source: [`experiments/rigorous_bench/xseal_fixed_block_scan_bench.cpp`](file:///home/hanrai/Projects/xseal/experiments/rigorous_bench/xseal_fixed_block_scan_bench.cpp)
  - Compilation:
    ```bash
    clang++ -O3 -march=native -std=c++20 -Iinclude experiments/rigorous_bench/xseal_fixed_block_scan_bench.cpp -o build/bin/xseal_fixed_block_scan_bench
    ```
- **SimdMin Fixed Block Scan (Rust)**:
  - Source: [`experiments/rigorous_bench/simdmin_fixed_block_scan_bench`](file:///home/hanrai/Projects/xseal/experiments/rigorous_bench/simdmin_fixed_block_scan_bench/)
  - Compilation:
    ```bash
    cargo build --manifest-path experiments/rigorous_bench/simdmin_fixed_block_scan_bench/Cargo.toml --release --target-dir build/rust_target
    cp build/rust_target/release/simdmin_fixed_block_scan_bench build/bin/
    ```

### 2.9. Coordinate Auditor & Fused Integration
- **Coordinate Audit (Rust)**:
  - Source: [`experiments/rigorous_bench/genomic_marker_coordinate_audit`](file:///home/hanrai/Projects/xseal/experiments/rigorous_bench/genomic_marker_coordinate_audit/)
  - Compilation:
    ```bash
    cargo build --manifest-path experiments/rigorous_bench/genomic_marker_coordinate_audit/Cargo.toml --release --target-dir build/rust_target
    cp build/rust_target/release/genomic_marker_coordinate_audit build/bin/
    ```
- **Fused Integration Wall Bench (C++)**:
  - Source: [`experiments/rigorous_bench/fused_integration_wall_bench.cpp`](file:///home/hanrai/Projects/xseal/experiments/rigorous_bench/fused_integration_wall_bench.cpp)
  - Compilation:
    ```bash
    clang++ -O3 -march=native -std=c++20 -Iinclude experiments/rigorous_bench/fused_integration_wall_bench.cpp -o build/bin/fused_integration_wall_bench
    ```

---

## 3. Per-table reproduction (`experiments/scripts/`)

Run only the scripts for the tables or sections you need. Shared setup (§0) and anchors in [`expected_metrics.yaml`](expected_metrics.yaml). Hot-RAM Table 1 scripts `source` [`build_hot_benches.sh`](rigorous_bench/scripts/build_hot_benches.sh).

| Target | Script | Output |
|--------|--------|--------|
| Table 1 + 3 (hot RAM, pos-only) | `rigorous_bench/scripts/repro_tab1_tab3_pos_only.sh` | `repro_tab1_tab3_out/` |
| Tables 1 + 2 (hot RAM + 128 KiB scaling) | `rigorous_bench/scripts/repro_paper_table12.sh` | `repro_paper_table12_out/` |
| Part I throughput (scan bench) | `scripts/35_part1_core_benches.sh` | `results/metrics_part1.json` |
| Table 4 encoder | `scripts/40_encoder_scaling.sh` | `results/metrics_encoder.json` |
| Table 5 parser | `scripts/50_parser_bench.sh` | `results/metrics_parser.json` |
| PMU / STLF (optional) | `PART2_PMU=1 scripts/55_part2_pmu.sh` | `results/metrics_part2.json` |
| minimap2 integration | `integration/minimap2_xseal_hardened/run_hg38_bench.sh` | `bench_hg38_out/` |
| Compare to yaml anchors | `python3 scripts/90_compare_report.py` (after metrics JSON exist) | `results/report.md` |
