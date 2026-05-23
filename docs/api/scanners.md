# XSeal Scanners: Design & Usage

XSeal provides two state-of-the-art, high-performance sequence scanners designed to locate structural sampling points (Minimizers and Syncmers) inside 2-bit encoded genomic sequences.

---

## 1. `XSealSyncmer`: AVX2-Accelerated Syncmer Scanner

`XSealSyncmer` identifies open or closed syncmers over a sequence using highly optimized AVX2 SIMD operations.

### Design & Optimizations
- **Parallel Hash Extraction**: Processes 8 s-mers in parallel. It uses vectorised 2-bit extraction and applies Murmur3-style mixed hashing directly in vector registers.
- **Ultra-Parallel Window Min**: Utilizes 4-way interleaved accumulators to find the minimum hash in a sliding window of size $W$, avoiding standard comparative branches.
- **Branchless Position Collection**: Uses vectorized look-up tables (`offset_lut`) to store found positions sequentially in memory without branching.

### API & Usage
```cpp
#include <xseal/xseal_syncmer.hpp>

// Initialize syncmer scanner: K-mer size, S-mer size
xseal::XSealSyncmer scanner(31, 11);

// seq_chunks: pointer to __m256i aligned 2-bit encoded sequence chunks
// num_chunks: number of 32-byte (128 bases) sequence chunks
// out_pos: target array to store the resulting absolute base positions of syncmers
// W: sliding window size
// T: target s-mer position within the k-mer (T=0 for open/closed syncmers)
// CLOSED: whether to find closed syncmers (T or W-1-T is the minimum)
// CANONICAL: whether to scan canonical (minimum of forward and reverse-complement) s-mers
size_t num_syncmers = scanner.scan<11, 0, true, true>(seq_chunks, num_chunks, out_pos);
```

---

## 2. `XSealMinimizer`: AVX2-Accelerated Minimizer Scanner

`XSealMinimizer` computes sequence minimizers using a vectorized sliding-window approach.

### Design & Optimizations
- **Vectorized PHI-Reduction**: Simultaneously tracks both the minimum hash values and their relative offsets inside vector registers using bit-packed 32-bit fields, avoiding separate tracking overhead.
- **Redundancy Elimination**: Directly compares each window's minimum position with the previously selected minimizer, emitting only new unique minimizers branchlessly.

### API & Usage
```cpp
#include <xseal/xseal_minimizer.hpp>

// Initialize minimizer scanner: K-mer size, W (window size, passed as template in scan)
xseal::XSealMinimizer scanner(31, 11);

// seq_chunks: pointer to 2-bit encoded chunks
// num_chunks: number of sequence chunks
// out_pos: target array to store absolute base positions of minimizers
// out_hash (optional): target array to store the minimizers' hash values
// WindowSize: sliding window size (e.g. 11)
// CANONICAL: whether to compile canonical minimizer scanning
size_t num_minimizers = scanner.scan<11, true>(seq_chunks, num_chunks, out_pos, out_hash);
```
