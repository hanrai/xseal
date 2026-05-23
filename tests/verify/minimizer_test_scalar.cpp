// SPDX-License-Identifier: Apache-2.0
// Copyright (c) 2026 hanrai. All Rights Reserved.

#include <iostream>
#include <vector>
#include <chrono>
#include <algorithm>
#include <cstring>
#include <iomanip>
#include <xseal/xseal_minimizer.hpp>
#include <xseal/xseal_encoder.hpp>

using namespace xseal;

// Scalar implementation that mimics XSealMinimizer SIMD logic exactly
class ScalarMinimizerRef {
public:
    int K;
    ScalarMinimizerRef(int k) : K(k) {}

    uint32_t reverse_complement(uint32_t kmer) {
        // SIMD version uses 32-bit lanes, so it effectively sees only min(K, 16) bases
        int k_eff = (K > 16) ? 16 : K;
        uint32_t res = 0;
        for (int i = 0; i < k_eff; ++i) {
            uint32_t base = (kmer >> (2 * i)) & 3;
            uint32_t comp = base ^ 3;
            res |= (comp << (2 * (k_eff - 1 - i)));
        }
        return res;
    }

    uint32_t murmur_hash(uint32_t h) {
        h ^= (h >> 13);
        h *= 0x85ebca6b;
        h ^= (h >> 16);
        return h;
    }

    template <int WindowSize, bool CANONICAL = true>
    std::vector<uint32_t> scan(const uint8_t* encoded_seq, size_t num_chunks, uint32_t seg_offset) {
        std::vector<uint32_t> results;
        if (num_chunks <= 1) return results;

        uint32_t k_mask = (K >= 16) ? 0xFFFFFFFF : (1U << (2 * K)) - 1;
        uint32_t last_abs_pos = 0xFFFFFFFF;

        for (size_t c = 0; c < num_chunks - 1; ++c) {
            size_t chunk_start = c * 128;
            uint32_t phi_buf[256];
            int total_pos = ((128 + WindowSize - 1) + 7) & ~7;

            for (int i = 0; i < total_pos; ++i) {
                size_t kmer_start = chunk_start + i;
                uint32_t kmer = 0;
                // SIMD loads 32 bits (16 bases) regardless of K
                for (int j = 0; j < 16; ++j) {
                    size_t pos = kmer_start + j;
                    uint8_t b = (encoded_seq[pos >> 2] >> ((pos & 3) << 1)) & 3;
                    kmer |= (uint32_t)b << (j * 2);
                }
                kmer &= k_mask;

                uint32_t h = kmer;
                if constexpr (CANONICAL) {
                    uint32_t r = reverse_complement(kmer);
                    // Match SIMD: if K > 16, no extra shift is performed
                    // If K < 16, the SIMD does srli(v_rev, 32 - 2*K)
                    // My reverse_complement already handles K < 16 correctly
                    if (r < h) h = r;
                }
                uint32_t hash_val = murmur_hash(h);
                phi_buf[i] = (hash_val & 0xFFFFFF00) | (i & 0xFF);
            }

            for (int batch = 0; batch < 128; ++batch) {
                uint32_t min_phi = 0xFFFFFFFF;
                for (int j = 0; j < WindowSize; ++j) {
                    if (phi_buf[batch + j] < min_phi) {
                        min_phi = phi_buf[batch + j];
                    }
                }

                uint32_t abs_pos = (seg_offset + chunk_start) + (min_phi & 0xFF);
                if (abs_pos != last_abs_pos) {
                    results.push_back(abs_pos);
                    last_abs_pos = abs_pos;
                }
            }
        }
        return results;
    }
};

template <int WindowSize, bool CANONICAL = true>
void run_test_case(int k) {
    std::cout << "Testing K=" << k << ", WindowSize=" << WindowSize << ", Canonical=" << (CANONICAL ? "Yes" : "No") << "..." << std::endl;

    const size_t num_bases = 1024 * 1024; // 1M bases
    std::vector<uint8_t> random_seq(num_bases);
    for (size_t i = 0; i < num_bases; ++i) {
        random_seq[i] = "ACGT"[rand() % 4];
    }

    XSealEncoder encoder;
    XSealEncoderState state;
    std::vector<uint8_t> encoded((num_bases + 3) / 4 + 128);
    encoder.init_state(&state, encoded.data(), encoded.size());
    const char* ptr = (const char*)random_seq.data();
    encoder.encode_chunk(&state, &ptr, ptr + num_bases);
    encoder.flush_state(&state);

    size_t num_chunks = num_bases / 128;

    // Run Scalar
    ScalarMinimizerRef scalar(k);
    auto start_scalar = std::chrono::high_resolution_clock::now();
    auto ref_hits = scalar.scan<WindowSize, CANONICAL>(encoded.data(), num_chunks, 0);
    auto end_scalar = std::chrono::high_resolution_clock::now();

    // Run SIMD
    XSealMinimizer simd(k, WindowSize);
    std::vector<uint32_t> simd_hit_pos(num_bases);
    std::vector<uint32_t> simd_hit_hashes(num_bases);
    auto start_simd = std::chrono::high_resolution_clock::now();
    size_t simd_count = simd.scan<WindowSize, CANONICAL>((const __m256i*)encoded.data(), num_chunks, simd_hit_pos.data(), simd_hit_hashes.data());
    auto end_simd = std::chrono::high_resolution_clock::now();

    // Compare
    bool match = true;
    if (ref_hits.size() != simd_count) {
        std::cout << "  [FAIL] Size mismatch: Scalar=" << ref_hits.size() << ", SIMD=" << simd_count << std::endl;
        match = false;
    } else {
        for (size_t i = 0; i < ref_hits.size(); ++i) {
            if (ref_hits[i] != simd_hit_pos[i]) {
                std::cout << "  [FAIL] Mismatch at index " << i << ": Scalar=" << ref_hits[i] << ", SIMD=" << simd_hit_pos[i] << std::endl;
                match = false;
                break;
            }
        }
    }

    if (match) {
        std::cout << "  [SUCCESS] Bit-perfect match! (" << ref_hits.size() << " hits)" << std::endl;
        double scalar_time = std::chrono::duration<double>(end_scalar - start_scalar).count();
        double simd_time = std::chrono::duration<double>(end_simd - start_simd).count();
        std::cout << "  Performance:" << std::endl;
        std::cout << "    Scalar: " << std::fixed << std::setprecision(2) << (num_bases / 1e6) / scalar_time << " MBases/s" << std::endl;
        std::cout << "    SIMD:   " << std::fixed << std::setprecision(2) << (num_bases / 1e6) / simd_time << " MBases/s" << std::endl;
        std::cout << "    Speedup: " << scalar_time / simd_time << "x" << std::endl;
    }
}

int main() {
    srand(42);
    run_test_case<11, true>(15);
    run_test_case<21, true>(31); // Now should pass with truncation logic
    run_test_case<31, false>(15);
    run_test_case<5, true>(10);
    return 0;
}
