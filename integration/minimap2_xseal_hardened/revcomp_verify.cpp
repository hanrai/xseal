// Verify harvester revcomp: compare SIMD vs scalar reference.
#include <cstdio>
#include <cstdint>
#include <xseal/xseal_kmer_harvester.hpp>

static uint64_t scalar_revcomp_canon(uint64_t c, int k)
{
    const uint64_t k_mask = (1ULL << (2 * k)) - 1;
    const uint32_t r_shift = 64u - (k << 1u);
    uint64_t rev = ~c & k_mask;
    rev = ((rev >> 2) & 0x3333333333333333ULL) | ((rev & 0x3333333333333333ULL) << 2);
    rev = ((rev >> 4) & 0x0F0F0F0F0F0F0F0FULL) | ((rev & 0x0F0F0F0F0F0F0F0FULL) << 4);
    rev = ((rev >> 8) & 0x00FF00FF00FF00FFULL) | ((rev & 0x00FF00FF00FF00FFULL) << 8);
    rev = ((rev >> 16) & 0x0000FFFF0000FFFFULL) | ((rev & 0x0000FFFF0000FFFFULL) << 16);
    rev = (rev >> 32) | (rev << 32);
    rev >>= r_shift;
    return rev < c ? rev : c;
}

// Broken SIMD path (old harvester logic)
static inline __m256i broken_revcomp(__m256i vr, uint32_t r_shift)
{
    const __m256i V_M0F = _mm256_set1_epi8(0x0F);
    const __m256i V_RC_LUT_L = _mm256_setr_epi8(
        15, 11, 7, 3, 14, 10, 6, 2, 13, 9, 5, 1, 12, 8, 4, 0, 15, 11, 7, 3, 14, 10, 6, 2, 13, 9,
        5, 1, 12, 8, 4, 0);
    const __m256i V_RC_LUT_H = V_RC_LUT_L;
    const __m256i V_REV_BYTES = _mm256_setr_epi8(7, 6, 5, 4, 3, 2, 1, 0, 15, 14, 13, 12, 11, 10, 9, 8,
                                                 7, 6, 5, 4, 3, 2, 1, 0, 15, 14, 13, 12, 11, 10, 9, 8);
    __m256i low = _mm256_and_si256(vr, V_M0F);
    __m256i high = _mm256_and_si256(_mm256_srli_epi16(vr, 4), V_M0F);
    vr = _mm256_or_si256(_mm256_shuffle_epi8(V_RC_LUT_L, low), _mm256_shuffle_epi8(V_RC_LUT_H, high));
    vr = _mm256_shuffle_epi8(vr, V_REV_BYTES);
    return _mm256_srli_epi64(vr, r_shift);
}

int main()
{
    const int k = 21;
    const uint64_t k_mask = (1ULL << (2 * k)) - 1;
    const uint32_t r_shift = 64u - (k << 1u);
    const __m256i v_k_mask = _mm256_set1_epi64x(k_mask);

    uint64_t mism_old = 0, mism_new = 0, unique_old = 0, unique_new = 0;
    uint64_t seen_old[1 << 16] = {}, seen_new[1 << 16] = {};

    for (uint64_t t = 0; t < (1ULL << 16); ++t) {
        uint64_t c = (t * 0x9e3779b97f4a7c15ULL) & k_mask;
        uint64_t ref = scalar_revcomp_canon(c, k);

        __m256i vf = _mm256_set1_epi64x(c);
        __m256i vr_old = broken_revcomp(_mm256_andnot_si256(vf, v_k_mask), r_shift);
        __m256i vr_new = xseal::xseal_revcomp_avx2_fast(vf, r_shift);
        uint64_t old_c = _mm256_extract_epi64(vr_old, 0);
        uint64_t new_c = _mm256_extract_epi64(vr_new, 0);
        if ((old_c < c ? old_c : c) != ref) mism_old++;
        if ((new_c < c ? new_c : c) != ref) mism_new++;

        uint64_t ho = (old_c < c ? old_c : c) >> (64 - 22);
        uint64_t hn = (new_c < c ? new_c : c) >> (64 - 22);
        if (!seen_old[ho]) {
            seen_old[ho] = 1;
            unique_old++;
        }
        if (!seen_new[hn]) {
            seen_new[hn] = 1;
            unique_new++;
        }
    }

    fprintf(stderr, "sample 65536 random k-mers (k=21)\n");
    fprintf(stderr, "canonical mismatches: broken_simd=%llu fixed_simd=%llu\n",
            (unsigned long long)mism_old, (unsigned long long)mism_new);
    fprintf(stderr, "unique high-22-bit canonical (proxy): broken=%llu fixed=%llu\n",
            (unsigned long long)unique_old, (unsigned long long)unique_new);
    return mism_new > 0 ? 1 : 0;
}
