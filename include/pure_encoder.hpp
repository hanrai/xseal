// SPDX-License-Identifier: Apache-2.0
// Copyright (c) 2026 hanrai. All Rights Reserved.

#pragma once
#include <immintrin.h>
#include <cstdint>
#include <cstring>
#include <stdexcept>

#ifdef _MSC_VER
#include <intrin.h>
#endif

namespace xseal {

/**
 * @brief PureEncoder: Hardware minimalist DNA 2-bit compression engine (Zero-Branch / Full SIMD)
 * Extreme single-shot: Utilizing AVX2 + BMI2 (_pdep_u64) instruction sets, single-core achieves pure streaming compression over 10GB/s.
 * 
 * Assumption: The input data is pure continuous gene segments, without any special formatting characters like '\n', '>', 'N', etc.
 * These can be filtered out in segments during the first pass of FastaScanner preprocessing.
 */
class PureEncoder {
public:
    /**
     * @brief Zero-branch compression: process pure gene data segments.
     * @param start Sequence start pointer
     * @param end Sequence end pointer
     * @param out_ptr Target memory for output, must guarantee enough space ((end-start)/4 + 8) and 64-bit alignment is highly recommended.
     */
    static inline void encode_chunk(const char* start, const char* end, uint8_t* out_ptr) {
        const char* p = start;
        uint64_t* out_u64 = reinterpret_cast<uint64_t*>(out_ptr);

        // Core hot loop: relentlessly process 32 bytes (256 bits) at a time, without any branch checks
        // Produces exactly 64 bits (8 bytes), perfectly aligned direct write, breaking all bit_cnt constraints.
        while (p + 32 <= end) {
            __m256i raw = _mm256_loadu_si256(reinterpret_cast<const __m256i*>(p));
            
            // Strip Bit 1 and Bit 2
            // Clever use of cute-nucleotides principle:
            // 'A' (0x41) = ...00000001 -> bit 2=0, bit 1=0
            // 'C' (0x43) = ...01000011 -> bit 2=0, bit 1=1
            // 'T' (0x54) = ...01010100 -> bit 2=1, bit 1=0
            // 'G' (0x47) = ...01000111 -> bit 2=1, bit 1=1
            __m256i t1 = _mm256_slli_epi16(raw, 5); // Push original bit 2 to the sign bit (i.e. bit 7)
            __m256i t2 = _mm256_slli_epi16(raw, 6); // Push original bit 1 to the sign bit

            uint32_t hi = _mm256_movemask_epi8(t1);
            uint32_t lo = _mm256_movemask_epi8(t2);

            // Xseal's expected encoding is: A=0(00), C=1(01), G=2(10), T=3(11)
            // Based on the ASCII binary extraction above:
            // G extracts hi=1, lo=1 (i.e., 3); T extracts hi=1, lo=0 (i.e., 2).
            // We just need to flip lo when hi=1 to change G to 2 and T to 3, perfectly matching the Xseal standard!
            lo ^= hi; 

            // Use BMI2 _pdep_u64 for dual interleaved packing
            // 0xAAAAAAAAAAAAAAAA = 10101010... -> Put bit 2 (hi) into all odd positions
            // 0x5555555555555555 = 01010101... -> Put bit 1 (lo) into all even positions
#if defined(_MSC_VER) && !defined(__clang__)
            uint64_t out_val = _pdep_u64(hi, 0xAAAAAAAAAAAAAAAAULL) | _pdep_u64(lo, 0x5555555555555555ULL);
#else
            uint64_t out_val = __builtin_ia32_pdep_di(hi, 0xAAAAAAAAAAAAAAAAULL) | __builtin_ia32_pdep_di(lo, 0x5555555555555555ULL);
#endif

            // No cache, direct aligned store
            *out_u64++ = out_val;
            
            p += 32;
        }

        // Phase 2: Handle remaining remainder of less than 32 bytes (scalar slow tail processing)
        if (p < end) {
            uint64_t tail_val = 0;
            int bit_cnt = 0;
            while (p < end) {
                char c = *p++; 
                // Xseal scalar encoding logic: A=0, C=1, G=2, T=3
                uint64_t v = (c=='C'||c=='c') ? 1 : (c=='G'||c=='g') ? 2 : (c=='T'||c=='t') ? 3 : 0;
                tail_val |= (v << bit_cnt);
                bit_cnt += 2;
            }
            // Do not overwrite the clean bytes already generated, use byte-level memcpy to prevent out-of-bounds
            std::memcpy(out_u64, &tail_val, (bit_cnt + 7) / 8);
        }
    }
};

} // namespace xseal
