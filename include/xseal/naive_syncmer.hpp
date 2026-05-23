// SPDX-License-Identifier: Apache-2.0
// Copyright (c) 2026 hanrai. All Rights Reserved.

#pragma once

#include <cstdint>
#include <vector>
#include <algorithm>

namespace xseal {

/**
 * NaiveSyncmer: A scalar, non-vectorized implementation of the syncmer scanner.
 * Used as a baseline to isolate the benefits of SIMD and the VSB design pattern.
 * Uses independent hashing (NtHash-like but scalar).
 */
class NaiveSyncmer {
public:
    NaiveSyncmer(int k, int s) : K(k), S(s) {
        k_mask = (k >= 16) ? 0xFFFFFFFF : (1U << (k * 2)) - 1;
        s_mask = (s >= 16) ? 0xFFFFFFFF : (1U << (s * 2)) - 1;
    }

    template <int W, int T = 0>
    size_t scan(const uint8_t* encoded_seq, size_t num_bases, uint32_t* out_pos) {
        if (num_bases < (size_t)K) return 0;
        size_t total_hits = 0;
        size_t num_windows = num_bases - K + 1;

        std::vector<uint32_t> window_hashes(W);

        for (size_t i = 0; i < num_windows; ++i) {
            // Compute hashes for the current window if needed, or maintain a sliding window of hashes
            // To be "fair" but "naive", we'll compute independent hashes for the window
            // Or we can be slightly better and just compute one new hash and slide.
            // Let's do the "Naive Independent" approach as suggested: compute S-mer hash independently.
            
            uint32_t min_h = 0xFFFFFFFF;
            uint32_t first_h = 0;
            uint32_t last_h = 0;

            for (int j = 0; j < W; ++j) {
                uint32_t h = get_smer_hash(encoded_seq, i + j);
                if (h < min_h) min_h = h;
                if (j == T) first_h = h;
                if (j == W - 1) last_h = h;
            }

            if (first_h == min_h || last_h == min_h) {
                out_pos[total_hits++] = (uint32_t)i;
            }
        }
        return total_hits;
    }

private:
    int K, S;
    uint32_t k_mask, s_mask;

    inline uint32_t get_smer_hash(const uint8_t* seq, size_t pos) const {
        uint32_t smer = 0;
        for (int j = 0; j < S; ++j) {
            size_t p = pos + j;
            uint8_t b = (seq[p >> 2] >> ((p & 3) << 1)) & 3;
            smer |= (uint32_t)b << (j * 2);
        }
        smer &= s_mask;

        // Simple Murmur-like mix
        uint32_t h = smer;
        h ^= (h >> 13);
        h *= 0x85ebca6b;
        h ^= (h >> 16);
        return h;
    }
};

} // namespace xseal
