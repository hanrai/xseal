// SPDX-License-Identifier: Apache-2.0
// Copyright (c) 2026 hanrai. All Rights Reserved.

#pragma once

#include <cstdint>
#include <vector>

// Scalar reference matching XSealMinimizer SIMD PHI protocol (see minimizer_test_scalar.cpp).
class MinimizerScalarRef {
public:
  int K;
  explicit MinimizerScalarRef(int k) : K(k) {}

  uint32_t reverse_complement(uint32_t kmer) const {
    int k_eff = (K > 16) ? 16 : K;
    uint32_t res = 0;
    for (int i = 0; i < k_eff; ++i) {
      uint32_t base = (kmer >> (2 * i)) & 3;
      uint32_t comp = base ^ 3;
      res |= (comp << (2 * (k_eff - 1 - i)));
    }
    return res;
  }

  static uint32_t murmur_hash(uint32_t h) {
    h ^= (h >> 13);
    h *= 0x85ebca6b;
    h ^= (h >> 16);
    return h;
  }

  struct Hit {
    uint32_t pos;
    uint32_t hash;
  };

  template <int WindowSize, bool CANONICAL = true>
  std::vector<Hit> scan_hits(const uint8_t *encoded_seq, size_t num_chunks,
                             uint32_t seg_offset) const {
    std::vector<Hit> results;
    if (num_chunks <= 1)
      return results;

    uint32_t k_mask = (K >= 16) ? 0xFFFFFFFFu : (1U << (2 * K)) - 1;
    uint32_t last_abs_pos = 0xFFFFFFFFu;

    for (size_t c = 0; c < num_chunks - 1; ++c) {
      size_t chunk_start = c * 128;
      uint32_t phi_buf[256];
      int total_pos = ((128 + WindowSize - 1) + 7) & ~7;

      for (int i = 0; i < total_pos; ++i) {
        size_t kmer_start = chunk_start + i;
        uint32_t kmer = 0;
        for (int j = 0; j < 16; ++j) {
          size_t pos = kmer_start + j;
          uint8_t b = (encoded_seq[pos >> 2] >> ((pos & 3) << 1)) & 3;
          kmer |= (uint32_t)b << (j * 2);
        }
        kmer &= k_mask;

        uint32_t h = kmer;
        if constexpr (CANONICAL) {
          uint32_t r = reverse_complement(kmer);
          if (r < h)
            h = r;
        }
        uint32_t hash_val = murmur_hash(h);
        phi_buf[i] = (hash_val & 0xFFFFFF00u) | (static_cast<uint32_t>(i) & 0xFFu);
      }

      for (int batch = 0; batch < 128; ++batch) {
        uint32_t min_phi = 0xFFFFFFFFu;
        for (int j = 0; j < WindowSize; ++j) {
          if (phi_buf[batch + j] < min_phi)
            min_phi = phi_buf[batch + j];
        }
        uint32_t abs_pos =
            static_cast<uint32_t>(seg_offset + chunk_start + (min_phi & 0xFFu));
        if (abs_pos != last_abs_pos) {
          results.push_back({abs_pos, min_phi & 0xFFFFFF00u});
          last_abs_pos = abs_pos;
        }
      }
    }
    return results;
  }
};
