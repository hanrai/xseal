# XSeal Harvesters: Design & Usage

XSeal utilizes a decoupled architecture where sampling location discovery (scanning) is separated from actual k-mer sequence extraction (harvesting). The harvester subsystem is responsible for collecting and hashing k-mer sequences at selected sampling positions.

---

## 1. Overview & Decoupled Architecture

In traditional sequence sampling, k-mer sequence extraction and hash calculations are tightly coupled inside the main scanning loop. This introduces significant instruction overhead and cache pollution.

XSeal's decoupled pipeline model solves this:
1. **Scanner**: Scans the 2-bit DNA stream using register-resident sliding windows, emitting *only* raw absolute positions (`out_pos`) of target syncmers/minimizers.
2. **Harvester**: In a separate parallel stage, reads the 2-bit DNA stream at the selected positions, extracts the full $K$-mer sequences, and computes their 64-bit canonical hashes.

---

## 2. API: `xseal_harvest_kmers`

The harvester is implemented as a highly optimized, template-driven free function:

```cpp
#include <xseal/xseal_kmer_harvester.hpp>

template <bool CANONICAL, bool HASH = true>
void xseal_harvest_kmers(const uint8_t *__restrict enc,
                         const uint32_t *__restrict pos, size_t n, int k,
                         uint64_t *__restrict out_hashes,
                         size_t prefetch_ahead = 0);
```

### Parameters
- **`enc`**: Pointer to the 2-bit encoded DNA sequence buffer.
- **`pos`**: Array containing absolute base coordinates (positions) of selected k-mers.
- **`n`**: Total number of positions to harvest.
- **`k`**: The k-mer size (up to 32 bases).
- **`out_hashes`**: Output array where the harvested 64-bit hashes will be written.
- **`prefetch_ahead`** (optional): Instruction pre-fetch threshold to hide memory access latency.

---

## 3. Micro-architectural Optimizations

- **4-Way SIMD Extraction**: Extracts four $K$-mers from the 2-bit bitstream simultaneously, resolving shift boundaries in parallel.
- **AVX2 Reverse Complement**: Employs an ultra-fast vectorised reverse-complement logic using raw byte shuffles (`_mm256_shuffle_epi8`) with pre-calculated look-up tables (`V_RC_LUT_L` and `V_RC_LUT_H`) to obtain canonical k-mers directly in vector registers.
- **AVX2 Murmur3 64-bit Mixer**: Computes 64-bit Murmur3 mixing hashes across 4-ways in parallel using vectorised multiplications.

---

## 4. Usage Example

```cpp
#include <xseal/xseal_kmer_harvester.hpp>

// enc_seq: 2-bit encoded genomic stream
// syncmer_positions: absolute positions returned by XSealSyncmer::scan
size_t num_positions = 50000;
std::vector<uint64_t> harvested_hashes(num_positions);

// Harvest 31-mers at the selected syncmer positions:
// CANONICAL = true (canonical k-mers)
// HASH = true (compute Murmur3 64-bit mixing hashes)
xseal::xseal_harvest_kmers<true, true>(
    enc_seq.data(),
    syncmer_positions.data(),
    num_positions,
    31,
    harvested_hashes.data()
);
```
