# XSeal: eXtreme Sequence Exact Alignment Locator (X-SEAL)

<p align="center">
  <img src="docs/xseal_logo.png" alt="XSeal Logo" width="300px">
</p>

[![License](https://img.shields.io/badge/License-Apache_2.0-blue.svg)](https://opensource.org/licenses/Apache-2.0)
[![Platform](https://img.shields.io/badge/Platform-Linux-lightgrey.svg)](https://linux.org)
[![Language](https://img.shields.io/badge/Language-C%2B%2B20-blue.svg)](https://isocpp.org/)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.20338868.svg)](https://doi.org/10.5281/zenodo.20338868)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.20351426.svg)](https://doi.org/10.5281/zenodo.20351426)

**XSeal** is a state-of-the-art, header-only C++20 library designed for ultra-high-speed sequence analysis. By leveraging modern micro-architectural optimizations and SIMD-BMI2 hybrid pipelines, XSeal pushes SIMD-resident k-mer sampling throughput toward the practical limits observed on modern AVX2 CPUs.

---

## 🚀 Flagship Performance

XSeal is engineered for performance-critical bioinformatics applications. Profiled under realistic execution conditions (warm page-cache, memory-resident sequence data) on an AMD Ryzen 5 5600X (Zen 3, 6 physical cores, 12 SMT threads), it delivers:

- **Ultra-Fast FASTA/FASTQ Ingestion**: DFA sequential parser ingestion rate exceeding **20 GiB/s** per core, bypassing intermediate string allocations.
- **Extreme Single-Threaded Scanning (AVX2)**: Sustained core scanning speedup of **1.77$\times$--2.10$\times$** over recent SIMD rolling-hash baselines:
  - **Minimizer Scanning**: **~1.45 Gbp/s** (vs. ~0.82 Gbp/s for `simd-minimizers`).
  - **Closed Syncmer Scanning**: **~1.70 Gbp/s** (vs. ~0.81 Gbp/s for `simd-minimizers`).
  - **Open Syncmer Scanning**: **~1.74 Gbp/s** (vs. ~0.84 Gbp/s for `simd-minimizers`).
- **Efficient Multi-Thread Scaling**: Plateaus at **10.58--10.60 Gbp/s** aggregate throughput at 12–16 threads, suggesting a transition from memory pressure toward vector execution saturation.
- **Lower Instruction Footprint**: Reduces retired instructions per base (instruction density) down to **7.8--9.3 insn/bp** (compared to 21.5--23.3 insn/bp in existing rolling-hash SIMD implementations).
- **Substantial End-to-End Acceleration**: Accelerates `minimap2`'s sequence ingestion and sketching stage by **44.73$\times$** (from 29.97s to 0.67s), leading to a **2.27$\times$ overall end-to-end index-building speedup** (from 51.59s to 22.71s).

---

## 🛠️ Key Architectural Innovations

XSeal's massive performance gains are achieved through several key micro-architectural optimizations detailed in our preprint:

- **VSB (Vector-Stream-Buffer) Pipeline Pattern**: A cache-aware design that structures genomic data flow onto the CPU's memory hierarchy:
  - *Stream (Cache-Level Staging)*: Tasks are chunked to $\sim$128 KB (2-bit encoded), maintaining strict L2 cache residency to mitigate OS scheduling jitter.
  - *Buffer (L1-Level Staging)*: Decouples phase 1 s-mer extraction and phase 2 reduction using a 640-byte internal scratch buffer (`hash_buf[160]`), avoiding register spill.
  - *Vector (Register-Level Coupling)*: Binds hot loop iterations directly to 256-bit YMM registers.
- **Stateless SIMD Hashing (vs. Stateful NtHash)**: Rather than maintaining stateful sliding rolling hash recurrence relations across vector lanes (which causes branch and register bookkeeping bloat), XSeal computes canonical s-mer hashes independently from scratch inside YMM registers, slashing instruction density by more than $2.6\times$.
- **4-Way Accumulator ILP Forcing via Fold Expressions**: Unrolls window reduction loops using C++17 fold expressions into four independent dependency chains, allowing modern Out-of-Order engines to schedule execution port cycles concurrently.
- **PHI (Packed Hash-Index) Reduction**: Combines the 24-bit MSB hash with its 8-bit local position offset into a single `uint32_t`. A single `vpminud` vector lane reduction simultaneously extracts the minimum hash and its coordinate, eliminating conditional branching.
- **Zero-Copy Boundary-Safe SIMD Streaming**: Streamlines direct memory mapping (`mmap`) processing using SIMD lookahead and an automated scalar fallback tail, avoiding segmentation faults at memory page boundaries.

## 💡 Initial Context: An Accidental Breakthrough

XSeal did not start as a targeted bioinformatics research project, but rather as a weekend escape. Disillusioned by the capital-dominated "dark forest" of Solana MEV (Maximal Extractable Value) trading, I asked an AI to recommend a completely unrelated field where I could just write C++ and squeeze out raw CPU cycles for fun. It suggested genomic k-mer sampling. 

What began as a therapeutic coding exercise to decompress unexpectedly resulted in shattering the performance ceiling of a foundational algorithm. For the philosophical and personal background of this journey, see [`docs/JOURNEY.md`](docs/JOURNEY.md).

## 📂 Project Structure

- **`include/`**: Header files (Header-only core).
- **`tests/verify/`**: Correctness verification suite.
- **`experiments/`**: Scripts, microbenchmarks, and reproducibility pipelines for our paper.
- **`docs/`**: API documentation and the developer's journey (`JOURNEY.md`).

---

## 🔬 Paper Reproducibility

For per-table reproduction scripts, see **[`experiments/README.md`](experiments/README.md)** §3.

## 📦 Quick Start

### Installation (C++ Header-only)
XSeal is a header-only library. Simply include the `include/` directory in your project.
```bash
git clone https://github.com/hanrai/xseal.git
# Add -I/path/to/xseal/include to your compiler flags
```
> [!WARNING]
> **Compiler Recommendation**: While XSeal supports GCC, extensive micro-architectural profiling indicates that **GCC-compiled binaries suffer a >30% performance degradation** compared to Clang. This is due to differences in AVX2 loop unrolling, instruction scheduling, and register allocation. To achieve the paper-grade optimal performance, we strongly recommend compiling with **Clang**.

### CLI Tool (High-Performance Benchmark)
Build the project and run the robust showcase benchmark tool:
```bash
mkdir build && cd build
cmake ..
make -j

# Run Syncmer scanning with 12 threads and mmap reader
./xseal_bench -m syncmer -r mmap -t 12 <input.fasta>
```

### Correctness Verification
XSeal has been validated to produce bit-perfect parity against a scalar PHI reference on **chromosome 1** of hg38 (not the full multi-contig scan). Download hg38 with `./experiments/scripts/download_hg38.sh`, then:
```bash
# Minimizer / syncmer streaming checks on chr1 only
./verify_streaming data/hg38.fa m 31 21

# Run open syncmer streaming verification on chr1
./verify_streaming data/hg38.fa so 31 11 0

# Run closed syncmer streaming verification on chr1
./verify_streaming data/hg38.fa sc 31 11 5
```

### Basic Library Usage (Minimizer)
```cpp
#include <xseal/xseal_minimizer.hpp>

// Initialize minimizer scanner: K-mer size (31), Window size (11)
xseal::XSealMinimizer scanner(31, 11);

// Allocate output buffers for positions and hashes
std::vector<uint32_t> out_pos(1000000);
std::vector<uint32_t> out_hash(1000000);

// Scan pre-encoded 2-bit chunks using AVX2 (Window size = 11, CANONICAL = true)
size_t num_minimizers = scanner.scan<11, true>(
    seq_chunks,      // const __m256i* encoded chunks
    num_chunks,      // number of 32-byte chunks
    out_pos.data(),
    out_hash.data()  // optional hash output buffer
);
```

---

## 🔬 Academic Citation

If you use XSeal in your research, please cite our work as:

### APA
```text
Fan, L. (2026). XSeal: Microarchitectural Optimization of SIMD-Resident Genomic Sampling Pipelines. Zenodo. https://doi.org/10.5281/zenodo.20338868
```

### BibTeX
```bibtex
@misc{fan_2026_20338868,
  author       = {Fan, Lei},
  title        = {XSeal: Microarchitectural Optimization of SIMD-
                   Resident Genomic Sampling Pipelines
                  },
  month        = may,
  year         = 2026,
  publisher    = {Zenodo},
  doi          = {10.5281/zenodo.20338868},
  url          = {https://doi.org/10.5281/zenodo.20338868},
}
```

---

## ⚖️ License
 
 XSeal is released under the **Apache License 2.0**.
 
 The Apache License 2.0 ensures that XSeal remains free for academic, research, commercial, and open-source projects. You are free to use, modify, and distribute the code, provided that you include the original copyright and license notice in any copy of the software/source.
 
 ### Third-Party Code
 This repository includes benchmark integrations and modified code derived from other open-source projects:
 
 - **Minimap2**: Unmodified and lightly-hardened components from Minimap2 (developed by Heng Li) are located under the `integration/minimap2_xseal_hardened/` directory.
 - **simd-minimizers**: Automatically downloaded and compiled for baseline comparative benchmarking scripts in `experiments/rigorous_bench/`.
 
 Both of these components remain under their original **MIT License**. Please see the `THIRD_PARTY_LICENSES/` directory for full attribution and license details.

---
© 2026 hanrai, XSeal Project. Developed for extreme performance.
