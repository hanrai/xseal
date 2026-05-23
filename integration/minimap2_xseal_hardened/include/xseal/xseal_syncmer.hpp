#pragma once

#include <cstdint>
#include <cstring>
#include <immintrin.h>
#include <utility>

#include <xseal/xseal_common.hpp>

namespace xseal {

struct LutInit {
  uint8_t table[256][8];
  LutInit() {
    for (int i = 0; i < 256; ++i) {
      int k = 0;
      for (int j = 0; j < 8; ++j) {
        if ((i >> j) & 1)
          table[i][k++] = j;
      }
      while (k < 8)
        table[i][k++] = 0;
    }
  }
};
static inline const LutInit offset_lut;

struct SyncmerShuffleLut {
  alignas(64) uint32_t shuffle_table[256][8];
  SyncmerShuffleLut() {
    for (int i = 0; i < 256; ++i) {
      int k = 0;
      for (int j = 0; j < 8; ++j) {
        if ((i >> j) & 1)
          shuffle_table[i][k++] = j;
      }
      while (k < 8)
        shuffle_table[i][k++] = 0;
    }
  }
};
static inline const SyncmerShuffleLut syncmer_shuffle_lut;

class XSealSyncmer {
private:
  int K, S;

  static inline const __m256i V_REV_BYTES_32 =
      _mm256_setr_epi8(3, 2, 1, 0, 7, 6, 5, 4, 11, 10, 9, 8, 15, 14, 13, 12, 3,
                       2, 1, 0, 7, 6, 5, 4, 11, 10, 9, 8, 15, 14, 13, 12);
  static inline const __m256i V_MASK_0F = _mm256_set1_epi8(0x0F);
  static inline const __m256i V_PACK =
      _mm256_setr_epi32(0, 2, 4, 6, 0, 2, 4, 6);
  static inline const __m256i V_RC_LUT_L = _mm256_setr_epi8(
      0xF0, 0xB0, 0x70, 0x30, 0xE0, 0xA0, 0x60, 0x20, 0xD0, 0x90, 0x50, 0x10,
      0xC0, 0x80, 0x40, 0x00, 0xF0, 0xB0, 0x70, 0x30, 0xE0, 0xA0, 0x60, 0x20,
      0xD0, 0x90, 0x50, 0x10, 0xC0, 0x80, 0x40, 0x00);
  static inline const __m256i V_RC_LUT_H = _mm256_setr_epi8(
      0x0F, 0x0B, 0x07, 0x03, 0x0E, 0x0A, 0x06, 0x02, 0x0D, 0x09, 0x05, 0x01,
      0x0C, 0x08, 0x04, 0x00, 0x0F, 0x0B, 0x07, 0x03, 0x0E, 0x0A, 0x06, 0x02,
      0x0D, 0x09, 0x05, 0x01, 0x0C, 0x08, 0x04, 0x00);
  static inline const __m256i V_SH_FWD = _mm256_setr_epi64x(0, 2, 4, 6);

  __m256i v_smer_mask;

  static inline __m256i fast_hash(__m256i x, __m256i v_murmur_c) {
    x = _mm256_xor_si256(x, _mm256_srli_epi32(x, 13));
    x = _mm256_mullo_epi32(x, v_murmur_c);
    x = _mm256_xor_si256(x, _mm256_srli_epi32(x, 16));
    return x;
  }

  template <bool CANONICAL>
  inline void extract_8(const uint8_t *fwd_ptr, int base_idx, int out_idx,
                        __m256i vsh_f, __m256i vmask, int s, __m256i v_murmur_c,
                        __m256i v_rev_bytes, __m256i v_m0f, __m256i v_lut_l,
                        __m256i v_lut_h, uint32_t *hash_buf) const {
    int fwd_byte = base_idx >> 2;
    __m256i fv1 = _mm256_castpd_si256(
        _mm256_broadcast_sd((const double *)(fwd_ptr + fwd_byte)));
    __m256i fv2 = _mm256_castpd_si256(
        _mm256_broadcast_sd((const double *)(fwd_ptr + fwd_byte + 1)));
    __m256i fs1 = _mm256_srlv_epi64(fv1, vsh_f);
    __m256i fs2 = _mm256_srlv_epi64(fv2, vsh_f);
    __m256i v_fwd_full =
        _mm256_blend_epi32(_mm256_permutevar8x32_epi32(fs1, V_PACK),
                           _mm256_permutevar8x32_epi32(fs2, V_PACK), 0xF0);
    __m256i v_fwd = _mm256_and_si256(v_fwd_full, vmask);

    __m256i v_h;
    if constexpr (CANONICAL) {
      __m256i vr = _mm256_shuffle_epi8(v_fwd_full, v_rev_bytes);
      __m256i lo = _mm256_and_si256(vr, v_m0f);
      __m256i hi = _mm256_and_si256(_mm256_srli_epi16(vr, 4), v_m0f);
      __m256i v_rev = _mm256_or_si256(_mm256_shuffle_epi8(v_lut_l, lo),
                                      _mm256_shuffle_epi8(v_lut_h, hi));
      v_rev = _mm256_and_si256(_mm256_srli_epi32(v_rev, 32 - 2 * s), vmask);
      v_h = _mm256_min_epu32(v_fwd, v_rev);
    } else
      v_h = v_fwd;
    _mm256_storeu_si256((__m256i *)&hash_buf[out_idx],
                        fast_hash(v_h, v_murmur_c));
  }

  template <bool CANONICAL>
  inline __m256i extract_8_vector(const uint8_t *fwd_ptr, int base_idx,
                                  __m256i vsh_f, __m256i vmask, int s, __m256i v_murmur_c,
                                  __m256i v_rev_bytes, __m256i v_m0f, __m256i v_lut_l,
                                  __m256i v_lut_h) const {
    int fwd_byte = base_idx >> 2;
    __m256i fv1 = _mm256_castpd_si256(
        _mm256_broadcast_sd((const double *)(fwd_ptr + fwd_byte)));
    __m256i fv2 = _mm256_castpd_si256(
        _mm256_broadcast_sd((const double *)(fwd_ptr + fwd_byte + 1)));
    __m256i fs1 = _mm256_srlv_epi64(fv1, vsh_f);
    __m256i fs2 = _mm256_srlv_epi64(fv2, vsh_f);
    __m256i v_fwd_full =
        _mm256_blend_epi32(_mm256_permutevar8x32_epi32(fs1, V_PACK),
                           _mm256_permutevar8x32_epi32(fs2, V_PACK), 0xF0);
    __m256i v_fwd = _mm256_and_si256(v_fwd_full, vmask);

    __m256i v_h;
    if constexpr (CANONICAL) {
      __m256i vr = _mm256_shuffle_epi8(v_fwd_full, v_rev_bytes);
      __m256i lo = _mm256_and_si256(vr, v_m0f);
      __m256i hi = _mm256_and_si256(_mm256_srli_epi16(vr, 4), v_m0f);
      __m256i v_rev = _mm256_or_si256(_mm256_shuffle_epi8(v_lut_l, lo),
                                      _mm256_shuffle_epi8(v_lut_h, hi));
      v_rev = _mm256_and_si256(_mm256_srli_epi32(v_rev, 32 - 2 * s), vmask);
      v_h = _mm256_min_epu32(v_fwd, v_rev);
    } else
      v_h = v_fwd;
    return fast_hash(v_h, v_murmur_c);
  }

  template <int W, size_t... J>
  static inline __m256i reduce_window_impl(const uint32_t *hash_buf, int batch,
                                           std::index_sequence<J...>) {
    __m256i va = _mm256_loadu_si256((const __m256i *)&hash_buf[batch]);
    __m256i vb = _mm256_loadu_si256((const __m256i *)&hash_buf[batch + 1]);
    __m256i vc = _mm256_loadu_si256((const __m256i *)&hash_buf[batch + 2]);
    __m256i vd = _mm256_loadu_si256((const __m256i *)&hash_buf[batch + 3]);
    (
        [&]() __attribute__((always_inline)) {
          if constexpr (J % 4 == 0)
            va = _mm256_min_epu32(
                va,
                _mm256_loadu_si256((const __m256i *)&hash_buf[batch + J + 4]));
          else if constexpr (J % 4 == 1)
            vb = _mm256_min_epu32(
                vb,
                _mm256_loadu_si256((const __m256i *)&hash_buf[batch + J + 4]));
          else if constexpr (J % 4 == 2)
            vc = _mm256_min_epu32(
                vc,
                _mm256_loadu_si256((const __m256i *)&hash_buf[batch + J + 4]));
          else
            vd = _mm256_min_epu32(
                vd,
                _mm256_loadu_si256((const __m256i *)&hash_buf[batch + J + 4]));
        }(),
        ...);
    return _mm256_min_epu32(_mm256_min_epu32(va, vb), _mm256_min_epu32(vc, vd));
  }

  template <int W>
  static inline __m256i reduce_window(const uint32_t *hash_buf, int batch) {
    if constexpr (W <= 1)
      return _mm256_loadu_si256((const __m256i *)&hash_buf[batch]);
    else if constexpr (W == 2)
      return _mm256_min_epu32(
          _mm256_loadu_si256((const __m256i *)&hash_buf[batch]),
          _mm256_loadu_si256((const __m256i *)&hash_buf[batch + 1]));
    else if constexpr (W == 3)
      return _mm256_min_epu32(
          _mm256_min_epu32(
              _mm256_loadu_si256((const __m256i *)&hash_buf[batch]),
              _mm256_loadu_si256((const __m256i *)&hash_buf[batch + 1])),
          _mm256_loadu_si256((const __m256i *)&hash_buf[batch + 2]));
    else
      return reduce_window_impl<W>(hash_buf, batch,
                                   std::make_index_sequence<W - 4>{});
  }

public:
  XSealSyncmer(int k, int s) : K(k), S(s) {
    uint32_t sm = (s >= 16) ? 0xFFFFFFFF : (1U << (s * 2)) - 1;
    v_smer_mask = _mm256_set1_epi32(sm);
  }

  /// Hash-only micro-benchmark: 2-bit extract + Murmur at every s-mer start (no window).
  template <bool CANONICAL = true>
  size_t hash_smers_only(const uint8_t *enc, size_t num_bases,
                         uint64_t *digest_out) const {
    const size_t num_hashes =
        (num_bases >= static_cast<size_t>(S)) ? num_bases - S + 1 : 0;
    if (num_hashes == 0) {
      if (digest_out)
        *digest_out = 0;
      return 0;
    }
    const __m256i vmask = v_smer_mask;
    const __m256i v_murmur_c = _mm256_set1_epi32(0x85ebca6b);
    const int s_val = S;
    
    __m256i v_digest = _mm256_setzero_si256();
    size_t pos = 0;
#pragma unroll 8
    while (pos + 8 <= num_hashes) {
      __m256i v_h = extract_8_vector<CANONICAL>(enc, static_cast<int>(pos), V_SH_FWD, vmask, s_val,
                                                v_murmur_c, V_REV_BYTES_32, V_MASK_0F, V_RC_LUT_L,
                                                V_RC_LUT_H);
      v_digest = _mm256_xor_si256(v_digest, v_h);
      pos += 8;
    }

    alignas(32) uint32_t hash_buf[8];
    _mm256_storeu_si256((__m256i *)hash_buf, v_digest);
    uint64_t digest = hash_buf[0] ^ hash_buf[1] ^ hash_buf[2] ^ hash_buf[3] ^
                      hash_buf[4] ^ hash_buf[5] ^ hash_buf[6] ^ hash_buf[7];

    for (; pos < num_hashes; ++pos) {
      extract_8<CANONICAL>(enc, static_cast<int>(pos), 0, V_SH_FWD, vmask, s_val,
                           v_murmur_c, V_REV_BYTES_32, V_MASK_0F, V_RC_LUT_L,
                           V_RC_LUT_H, hash_buf);
      digest ^= hash_buf[0];
    }
    if (digest_out)
      *digest_out = digest;
    asm volatile("" : "+r"(digest) : : "memory");
    return num_hashes;
  }

  template <int W, int T = 0, bool CLOSED = true, bool CANONICAL = true>
  size_t scan(const __m256i *seq_chunks, size_t num_chunks, uint32_t *out_pos,
              size_t actual_bases = 0) {
    if (num_chunks == 0)
      return 0;
    alignas(32) uint32_t hash_buf[160];
    size_t total_out_idx = 0;
    uint32_t chunk_pos = 0;
    constexpr int total_pos = ((128 + W - 1) + 7) & ~7;
    const __m256i vmask = v_smer_mask;
    const __m256i v_murmur_c = _mm256_set1_epi32(0x85ebca6b);
    const int s_val = S;
    static const __m256i v_lane_idx = _mm256_setr_epi32(0, 1, 2, 3, 4, 5, 6, 7);

    size_t total_bases_to_use =
        (actual_bases > 0) ? actual_bases : (num_chunks * 128);
    size_t num_windows =
        (total_bases_to_use >= (size_t)K) ? (total_bases_to_use - K + 1) : 0;

    // V10.0 Elite-Scan: Manual Inlining & Port Optimization
    size_t num_simd_chunks = num_windows / 128;
    const __m256i v_sh_f = V_SH_FWD;
    const __m256i v_pack = V_PACK;

    for (size_t c = 0; c < num_simd_chunks; ++c) {
      const uint8_t *fwd_ptr = (const uint8_t *)&seq_chunks[c];

// Inline Hash Extraction
#pragma unroll 19
      for (int i = 0; i < total_pos; i += 8) {
        int fwd_byte = i >> 2;
        __m256i fv1 = _mm256_castpd_si256(
            _mm256_broadcast_sd((const double *)(fwd_ptr + fwd_byte)));
        __m256i fv2 = _mm256_castpd_si256(
            _mm256_broadcast_sd((const double *)(fwd_ptr + fwd_byte + 1)));
        __m256i fs1 = _mm256_srlv_epi64(fv1, v_sh_f);
        __m256i fs2 = _mm256_srlv_epi64(fv2, v_sh_f);
        __m256i v_fwd_full =
            _mm256_blend_epi32(_mm256_permutevar8x32_epi32(fs1, v_pack),
                               _mm256_permutevar8x32_epi32(fs2, v_pack), 0xF0);
        __m256i v_fwd = _mm256_and_si256(v_fwd_full, vmask);
        __m256i v_h;
        if constexpr (CANONICAL) {
          __m256i vr = _mm256_shuffle_epi8(v_fwd_full, V_REV_BYTES_32);
          __m256i lo = _mm256_and_si256(vr, V_MASK_0F);
          __m256i hi = _mm256_and_si256(_mm256_srli_epi16(vr, 4), V_MASK_0F);
          __m256i v_rev = _mm256_or_si256(_mm256_shuffle_epi8(V_RC_LUT_L, lo),
                                          _mm256_shuffle_epi8(V_RC_LUT_H, hi));
          v_rev =
              _mm256_and_si256(_mm256_srli_epi32(v_rev, 32 - 2 * s_val), vmask);
          v_h = _mm256_min_epu32(v_fwd, v_rev);
        } else
          v_h = v_fwd;
        _mm256_storeu_si256((__m256i *)&hash_buf[i],
                            fast_hash(v_h, v_murmur_c));
      }

      // Phase 2: Ultra-Parallel Window Min (4-way interleaved accumulators)
#pragma unroll 16
      for (int batch = 0; batch < 128; batch += 8) {
        __m256i v_m0 = _mm256_loadu_si256((const __m256i *)&hash_buf[batch]);
        __m256i v_m1 =
            _mm256_loadu_si256((const __m256i *)&hash_buf[batch + 1]);
        __m256i v_m2 =
            _mm256_loadu_si256((const __m256i *)&hash_buf[batch + 2]);
        __m256i v_m3 =
            _mm256_loadu_si256((const __m256i *)&hash_buf[batch + 3]);

        auto unroll_min = [&]<size_t... Is>(std::index_sequence<Is...>) {
          (
              [&]() {
                __m256i vn = _mm256_loadu_si256(
                    (const __m256i *)&hash_buf[batch + Is + 4]);
                if constexpr (Is % 4 == 0)
                  v_m0 = _mm256_min_epu32(v_m0, vn);
                else if constexpr (Is % 4 == 1)
                  v_m1 = _mm256_min_epu32(v_m1, vn);
                else if constexpr (Is % 4 == 2)
                  v_m2 = _mm256_min_epu32(v_m2, vn);
                else
                  v_m3 = _mm256_min_epu32(v_m3, vn);
              }(),
              ...);
        };
        unroll_min(std::make_index_sequence<W - 4>{});
        __m256i v_min = _mm256_min_epu32(_mm256_min_epu32(v_m0, v_m1),
                                         _mm256_min_epu32(v_m2, v_m3));

        __m256i v_ht =
            _mm256_loadu_si256((const __m256i *)&hash_buf[batch + T]);
        __m256i v_eq = _mm256_cmpeq_epi32(v_ht, v_min);
        if constexpr (CLOSED) {
          __m256i v_hl =
              _mm256_loadu_si256((const __m256i *)&hash_buf[batch + W - 1]);
          v_eq = _mm256_or_si256(v_eq, _mm256_cmpeq_epi32(v_hl, v_min));
        }

        int mask = _mm256_movemask_ps((__m256)v_eq) & 0xFF;

        __m256i v_lut32 = _mm256_cvtepu8_epi32(
            _mm_loadl_epi64((const __m128i *)offset_lut.table[mask]));
        __m256i v_res_pos =
            _mm256_add_epi32(_mm256_set1_epi32(chunk_pos + batch), v_lut32);

        // THE BRANCHLESS STORE: Always write, advance index only by popcount
        _mm256_storeu_si256((__m256i *)&out_pos[total_out_idx], v_res_pos);
        total_out_idx += _mm_popcnt_u32(mask);
      }
      chunk_pos += 128;
    }

    // 2. MASKED TAIL (Single SIMD block for the remaining windows)
    size_t remaining = num_windows % 128;
    if (remaining > 0) {
      const uint8_t *fwd_ptr = (const uint8_t *)&seq_chunks[num_simd_chunks];
      for (int i = 0; i < total_pos; i += 8) {
        extract_8<CANONICAL>(fwd_ptr, i, i, V_SH_FWD, vmask, s_val, v_murmur_c,
                             V_REV_BYTES_32, V_MASK_0F, V_RC_LUT_L, V_RC_LUT_H,
                             hash_buf);
      }
      for (int batch = 0; batch < 128; batch += 8) {
        __m256i v_min = reduce_window<W>(hash_buf, batch);
        __m256i v_ht =
            _mm256_loadu_si256((const __m256i *)&hash_buf[batch + T]);
        __m256i v_eq = _mm256_cmpeq_epi32(v_ht, v_min);
        if constexpr (CLOSED) {
          __m256i v_hl =
              _mm256_loadu_si256((const __m256i *)&hash_buf[batch + W - 1]);
          v_eq = _mm256_or_si256(v_eq, _mm256_cmpeq_epi32(v_hl, v_min));
        }
        int mask = _mm256_movemask_ps((__m256)v_eq) & 0xFF;
        int batch_limit = remaining - batch;
        if (batch_limit <= 0)
          mask = 0;
        else if (batch_limit < 8)
          mask &= (1 << batch_limit) - 1;

        if (mask == 0)
          continue;

        __m256i v_lut32 = _mm256_cvtepu8_epi32(
            _mm_loadl_epi64((const __m128i *)offset_lut.table[mask]));
        __m256i v_res_pos =
            _mm256_add_epi32(_mm256_set1_epi32(chunk_pos + batch), v_lut32);
        _mm256_storeu_si256((__m256i *)&out_pos[total_out_idx], v_res_pos);
        total_out_idx += _mm_popcnt_u32(mask);
      }
    }

    return total_out_idx;
  }
};
} // namespace xseal