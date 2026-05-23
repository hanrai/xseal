// SPDX-License-Identifier: Apache-2.0
// Copyright (c) 2026 hanrai. All Rights Reserved.

#include <iostream>
#include <vector>
#include <chrono>
#include <algorithm>
#include <cstring>
#include <iomanip>
#include <xseal/xseal_syncmer.hpp>
#include <xseal/xseal_encoder.hpp>

using namespace xseal;

class ScalarSyncmerRef {
public:
    int K;
    int S;
    ScalarSyncmerRef(int k, int s) : K(k), S(s) {}

    uint32_t reverse_complement(uint32_t smer) {
        int s_eff = (S > 16) ? 16 : S;
        uint32_t res = 0;
        for (int i = 0; i < s_eff; ++i) {
            uint32_t base = (smer >> (2 * i)) & 3;
            uint32_t comp = base ^ 3;
            res |= (comp << (2 * (s_eff - 1 - i)));
        }
        return res;
    }

    uint32_t murmur_hash(uint32_t h) {
        h ^= (h >> 13);
        h *= 0x85ebca6b;
        h ^= (h >> 16);
        return h;
    }

    template <int WindowSize, int T = 0, bool CLOSED = true, bool CANONICAL = true>
    std::vector<uint32_t> scan(const uint8_t* encoded_seq, size_t num_chunks, uint32_t seg_offset, size_t actual_bases) {
        std::vector<uint32_t> results;
        if (num_chunks == 0) return results;

        size_t total_bases_to_use = (actual_bases > 0) ? actual_bases : (num_chunks * 128);
        size_t num_windows = (total_bases_to_use >= (size_t)K) ? (total_bases_to_use - K + 1) : 0;
        size_t num_simd_chunks = num_windows / 128;

        uint32_t s_mask = (S >= 16) ? 0xFFFFFFFF : (1U << (2 * S)) - 1;

        for (size_t c = 0; c < num_simd_chunks; ++c) {
            size_t chunk_start = c * 128;
            uint32_t phi_buf[256];
            constexpr int total_pos = ((128 + WindowSize - 1) + 7) & ~7;

            // Generate hashes for all s-mers starting in this extended chunk
            for (int i = 0; i < total_pos; ++i) {
                size_t smer_start = chunk_start + i;
                uint32_t smer = 0;
                for (int j = 0; j < 16; ++j) {
                    size_t pos = smer_start + j;
                    uint8_t b = (encoded_seq[pos >> 2] >> ((pos & 3) << 1)) & 3;
                    smer |= (uint32_t)b << (j * 2);
                }
                smer &= s_mask;

                uint32_t h = smer;
                if constexpr (CANONICAL) {
                    uint32_t r = reverse_complement(smer);
                    if (r < h) h = r;
                }
                phi_buf[i] = murmur_hash(h);
            }

            // Perform Window sliding over 128 positions
            for (int batch = 0; batch < 128; ++batch) {
                uint32_t min_phi = 0xFFFFFFFF;
                for (int j = 0; j < WindowSize; ++j) {
                    if (phi_buf[batch + j] < min_phi) {
                        min_phi = phi_buf[batch + j];
                    }
                }

                // Check syncmer condition
                bool is_syncmer = false;
                if constexpr (CLOSED) {
                    if (phi_buf[batch + T] == min_phi || phi_buf[batch + WindowSize - 1] == min_phi) {
                        is_syncmer = true;
                    }
                } else {
                    if (phi_buf[batch + T] == min_phi) {
                        is_syncmer = true;
                    }
                }

                if (is_syncmer) {
                    uint32_t abs_pos = (seg_offset + chunk_start) + batch;
                    results.push_back(abs_pos);
                }
            }
        }
        return results;
    }
};

template <int K, int S, int T = 0, bool CLOSED = true, bool CANONICAL = true>
void run_syncmer_test_case() {
    constexpr int WindowSize = K - S + 1;
    std::cout << "Testing K=" << K << ", S=" << S << ", W=" << WindowSize 
              << ", T=" << T << ", Canonical=" << (CANONICAL ? "Yes" : "No") << "..." << std::endl;

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

    size_t num_windows = (num_bases - K + 1) / 128 * 128;
    size_t test_bases = num_windows + K - 1;

    // Run Scalar
    ScalarSyncmerRef scalar(K, S);
    auto start_scalar = std::chrono::high_resolution_clock::now();
    auto ref_hits = scalar.scan<WindowSize, T, CLOSED, CANONICAL>(encoded.data(), num_chunks, 0, test_bases);
    auto end_scalar = std::chrono::high_resolution_clock::now();

    // Run SIMD
    XSealSyncmer simd(K, S);
    std::vector<uint32_t> simd_hit_pos(num_bases);
    auto start_simd = std::chrono::high_resolution_clock::now();
    size_t simd_count = simd.scan<WindowSize, T, CLOSED, CANONICAL>((const __m256i*)encoded.data(), num_chunks, simd_hit_pos.data(), test_bases);
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
    // Standard configurations
    run_syncmer_test_case<31, 11, 0, true, true>();
    run_syncmer_test_case<21, 11, 0, true, true>();
    run_syncmer_test_case<15, 5, 0, true, true>();
    run_syncmer_test_case<31, 11, 0, false, true>(); // Open Syncmer
    run_syncmer_test_case<21, 11, 0, true, false>(); // Non-canonical
    return 0;
}
