#pragma once

#include <cassert>
#include <cstdint>

#if defined(__GNUC__) || defined(__clang__)
#include <immintrin.h>
#else
#include <intrin.h>
#endif

namespace xseal {

// --- Static Constants for Harvester ---
namespace detail {
    static inline const __m256i V_MURMUR_C = _mm256_set1_epi64x(0xff51afd7ed558ccdULL);
    static inline const __m256i V_M32 = _mm256_set1_epi64x(0xFFFFFFFFULL);
    static inline const __m256i V_MURMUR_C_LO = _mm256_and_si256(V_MURMUR_C, V_M32);
    static inline const __m256i V_MURMUR_C_HI = _mm256_srli_epi64(V_MURMUR_C, 32);

    static inline const __m256i V_RC_LUT_L = _mm256_setr_epi8(
        0xF0, 0xB0, 0x70, 0x30, 0xE0, 0xA0, 0x60, 0x20, 0xD0, 0x90, 0x50, 0x10, 0xC0, 0x80, 0x40,
        0x00, 0xF0, 0xB0, 0x70, 0x30, 0xE0, 0xA0, 0x60, 0x20, 0xD0, 0x90, 0x50, 0x10, 0xC0, 0x80,
        0x40, 0x00);
    static inline const __m256i V_RC_LUT_H = _mm256_setr_epi8(
        0x0F, 0x0B, 0x07, 0x03, 0x0E, 0x0A, 0x06, 0x02, 0x0D, 0x09, 0x05, 0x01, 0x0C, 0x08, 0x04,
        0x00, 0x0F, 0x0B, 0x07, 0x03, 0x0E, 0x0A, 0x06, 0x02, 0x0D, 0x09, 0x05, 0x01, 0x0C, 0x08,
        0x04, 0x00);
    static inline const __m256i V_M0F = _mm256_set1_epi8(0x0F);
    static inline const __m256i V_REV_BYTES = _mm256_setr_epi8(
        7, 6, 5, 4, 3, 2, 1, 0, 15, 14, 13, 12, 11, 10, 9, 8,
        7, 6, 5, 4, 3, 2, 1, 0, 15, 14, 13, 12, 11, 10, 9, 8);
    static inline const __m256i V_H_BIT = _mm256_set1_epi64x(0x8000000000000000ULL);
}

/**
 * @brief Murmur3-style 64-bit Mixer (AVX2 optimized)
 * Matches the style of the syncmer scanner but for 64-bit k-mers.
 */
static inline __m256i xseal_murmur3_mix64_avx2(__m256i x) {
    // x ^= x >> 33
    x = _mm256_xor_si256(x, _mm256_srli_epi64(x, 33));
    
    __m256i x_hi = _mm256_srli_epi64(x, 32);
    __m256i x_lo = _mm256_and_si256(x, detail::V_M32);
    
    // Exact 64-bit mullo: (lo*lo) + (lo*hi + hi*lo) << 32
    __m256i r_ll = _mm256_mul_epu32(x_lo, detail::V_MURMUR_C_LO);
    __m256i r_lh = _mm256_mul_epu32(x_lo, detail::V_MURMUR_C_HI);
    __m256i r_hl = _mm256_mul_epu32(x_hi, detail::V_MURMUR_C_LO);
    __m256i r_mid = _mm256_add_epi64(r_lh, r_hl);
    x = _mm256_add_epi64(r_ll, _mm256_slli_epi64(r_mid, 32));
    
    // x ^= x >> 33
    x = _mm256_xor_si256(x, _mm256_srli_epi64(x, 33));
    return x;
}

/**
 * @brief Ultra-Fast DNA Reverse Complement (AVX2 PSHUFB)
 */
static inline __m256i xseal_revcomp_avx2_fast(__m256i vr, uint32_t r_shift) {
    __m256i low = _mm256_and_si256(vr, detail::V_M0F);
    __m256i high = _mm256_and_si256(_mm256_srli_epi16(vr, 4), detail::V_M0F);
    vr = _mm256_or_si256(_mm256_shuffle_epi8(detail::V_RC_LUT_L, low),
                         _mm256_shuffle_epi8(detail::V_RC_LUT_H, high));
    vr = _mm256_shuffle_epi8(vr, detail::V_REV_BYTES);
    return _mm256_srli_epi64(vr, r_shift);
}

/**
 * @brief XSeal Harvester - Ultimate Modular Edition
 * @tparam CANONICAL Whether to find the canonical k-mer min.
 * @tparam HASH Whether to compute the 64-bit Wyhash mixer.
 */
template <bool CANONICAL, bool HASH = true>
__attribute__((noinline)) void xseal_harvest_kmers(const uint8_t *__restrict enc,
                                                   const uint32_t *__restrict pos, size_t n, int k,
                                                   uint64_t *__restrict out_hashes,
                                                   size_t /*prefetch_ahead*/ = 0) {
    const uint32_t r_shift = 64u - (static_cast<unsigned>(k) << 1u);
    const uint64_t k_mask = (k == 32) ? ~0ULL : ((1ULL << (2 * k)) - 1);
    const __m256i v_k_mask = _mm256_set1_epi64x(k_mask);

    size_t i = 0;

    for (; i + 4 <= n; i += 4) {
        uint32_t p0 = pos[i];
        uint32_t p1 = pos[i+1];
        uint32_t p2 = pos[i+2];
        uint32_t p3 = pos[i+3];

        auto extract = [&](uint32_t p) -> uint64_t {
            const uint64_t *ptr = reinterpret_cast<const uint64_t *>(enc + (p >> 2));
            uint32_t shift = (p & 3) << 1;
            return ((ptr[0] >> shift) | (ptr[1] << (64 - shift))) & k_mask;
        };

        __m256i vf = _mm256_setr_epi64x(extract(p0), extract(p1), extract(p2), extract(p3));

        if constexpr (CANONICAL) {
            __m256i vr = xseal_revcomp_avx2_fast(vf, r_shift);
            __m256i v_mask = _mm256_cmpgt_epi64(_mm256_xor_si256(vf, detail::V_H_BIT), 
                                                _mm256_xor_si256(vr, detail::V_H_BIT));
            vf = _mm256_blendv_epi8(vf, vr, v_mask);
        }

        if constexpr (HASH) {
            _mm256_storeu_si256((__m256i *)&out_hashes[i], xseal_murmur3_mix64_avx2(vf));
        } else {
            _mm256_storeu_si256((__m256i *)&out_hashes[i], vf);
        }
    }

    // Tail
    for (; i < n; ++i) {
        uint32_t p = pos[i];
        const uint64_t *ptr = reinterpret_cast<const uint64_t *>(enc + (p >> 2));
        uint32_t shift = (p & 3) << 1;
        uint64_t c = ((ptr[0] >> shift) | (ptr[1] << (64 - shift))) & k_mask;
        if constexpr (CANONICAL) {
            uint64_t rev = ~c & k_mask;
            rev = ((rev >> 2) & 0x3333333333333333ULL) | ((rev & 0x3333333333333333ULL) << 2);
            rev = ((rev >> 4) & 0x0F0F0F0F0F0F0F0FULL) | ((rev & 0x0F0F0F0F0F0F0F0FULL) << 4);
            rev = ((rev >> 8) & 0x00FF00FF00FF00FFULL) | ((rev & 0x00FF00FF00FF00FFULL) << 8);
            rev = ((rev >> 16) & 0x0000FFFF0000FFFFULL) | ((rev & 0x0000FFFF0000FFFFULL) << 16);
            rev = (rev >> 32) | (rev << 32);
            rev >>= r_shift;
            if (rev < c) c = rev;
        }
        if constexpr (HASH) {
            c ^= c >> 33;
            c *= 0xff51afd7ed558ccdULL;
            c ^= c >> 33;
            out_hashes[i] = c;
        } else {
            out_hashes[i] = c;
        }
    }
}

} // namespace xseal
