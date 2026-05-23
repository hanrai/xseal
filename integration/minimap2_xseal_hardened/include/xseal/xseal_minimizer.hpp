#pragma once

#include <cstdint>
#include <cstring>
#include <immintrin.h>
#include <utility>
#include <vector>

#include <xseal/xseal_common.hpp>

namespace xseal {

struct MinimizerShuffleLut {
  alignas(64) uint32_t shuffle_table[256][8];
  MinimizerShuffleLut() {
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

static inline const MinimizerShuffleLut minimizer_shuffle_lut;

class XSealMinimizer {
private:
  int K;

  template <int WindowSize, size_t... J>
  static inline __m256i reduce_phi_impl(const uint32_t *phi_buf, int batch,
                                        std::index_sequence<J...>) {
    __m256i va = _mm256_loadu_si256((const __m256i *)&phi_buf[batch]);
    __m256i vb = _mm256_loadu_si256((const __m256i *)&phi_buf[batch + 1]);
    __m256i vc = _mm256_loadu_si256((const __m256i *)&phi_buf[batch + 2]);
    __m256i vd = _mm256_loadu_si256((const __m256i *)&phi_buf[batch + 3]);
    (
        [&]() __attribute__((always_inline)) {
          if constexpr (J % 4 == 0)
            va = _mm256_min_epu32(
                va,
                _mm256_loadu_si256((const __m256i *)&phi_buf[batch + J + 4]));
          else if constexpr (J % 4 == 1)
            vb = _mm256_min_epu32(
                vb,
                _mm256_loadu_si256((const __m256i *)&phi_buf[batch + J + 4]));
          else if constexpr (J % 4 == 2)
            vc = _mm256_min_epu32(
                vc,
                _mm256_loadu_si256((const __m256i *)&phi_buf[batch + J + 4]));
          else
            vd = _mm256_min_epu32(
                vd,
                _mm256_loadu_si256((const __m256i *)&phi_buf[batch + J + 4]));
        }(),
        ...);
    return _mm256_min_epu32(_mm256_min_epu32(va, vb), _mm256_min_epu32(vc, vd));
  }

  template <int WindowSize>
  static inline __m256i reduce_phi(const uint32_t *phi_buf, int batch) {
    if constexpr (WindowSize <= 4) {
      __m256i v = _mm256_loadu_si256((const __m256i *)&phi_buf[batch]);
      if constexpr (WindowSize >= 2)
        v = _mm256_min_epu32(
            v, _mm256_loadu_si256((const __m256i *)&phi_buf[batch + 1]));
      if constexpr (WindowSize >= 3)
        v = _mm256_min_epu32(
            v, _mm256_loadu_si256((const __m256i *)&phi_buf[batch + 2]));
      if constexpr (WindowSize >= 4)
        v = _mm256_min_epu32(
            v, _mm256_loadu_si256((const __m256i *)&phi_buf[batch + 3]));
      return v;
    } else
      return reduce_phi_impl<WindowSize>(
          phi_buf, batch, std::make_index_sequence<WindowSize - 4>{});
  }

public:
  XSealMinimizer(int k, [[maybe_unused]] int w) : K(k) {}

  static inline uint32_t calculate_kmer_hash(const uint8_t* encoded_seq, size_t pos, int k, bool canonical = true) {
    uint32_t k_mask = (k >= 16) ? 0xFFFFFFFF : (1U << (k * 2)) - 1;
    uint32_t kmer = 0;
    for (int j = 0; j < 16; ++j) { // SIMD limit
        size_t p = pos + j;
        uint8_t b = (encoded_seq[p >> 2] >> ((p & 3) << 1)) & 3;
        kmer |= (uint32_t)b << (j * 2);
    }
    kmer &= k_mask;

    uint32_t h = kmer;
    if (canonical) {
        uint32_t r = 0;
        int k_eff = (k > 16) ? 16 : k;
        for (int i = 0; i < k_eff; ++i) {
            uint32_t base = (kmer >> (2 * i)) & 3;
            uint32_t comp = base ^ 3;
            r |= (comp << (2 * (k_eff - 1 - i)));
        }
        if (r < h) h = r;
    }

    h ^= (h >> 13);
    h *= 0x85ebca6b;
    h ^= (h >> 16);
    return h;
  }

  template <int WindowSize, bool CANONICAL = true>
  size_t scan(const __m256i *seq_chunks, size_t num_chunks, uint32_t *out_pos, uint32_t *out_hash = nullptr) {
    if (num_chunks <= 1)
      return 0;
    alignas(32) uint32_t phi_buf[256];
    size_t total_out_idx = 0;
    uint32_t chunk_pos = 0;
    constexpr int total_pos = ((128 + WindowSize - 1) + 7) & ~7;
    const __m256i v_murmur_c = _mm256_set1_epi32(0x85ebca6b);
    const __m256i v_rev_bytes =
        _mm256_setr_epi8(3, 2, 1, 0, 7, 6, 5, 4, 11, 10, 9, 8, 15, 14, 13, 12,
                         3, 2, 1, 0, 7, 6, 5, 4, 11, 10, 9, 8, 15, 14, 13, 12);
    const __m256i v_m0f = _mm256_set1_epi8(0x0F);
    const __m256i v_lut_l = _mm256_setr_epi8(
        0xF0, 0xB0, 0x70, 0x30, 0xE0, 0xA0, 0x60, 0x20, 0xD0, 0x90, 0x50, 0x10,
        0xC0, 0x80, 0x40, 0x00, 0xF0, 0xB0, 0x70, 0x30, 0xE0, 0xA0, 0x60, 0x20,
        0xD0, 0x90, 0x50, 0x10, 0xC0, 0x80, 0x40, 0x00);
    const __m256i v_lut_h = _mm256_setr_epi8(
        0x0F, 0x0B, 0x07, 0x03, 0x0E, 0x0A, 0x06, 0x02, 0x0D, 0x09, 0x05, 0x01,
        0x0C, 0x08, 0x04, 0x00, 0x0F, 0x0B, 0x07, 0x03, 0x0E, 0x0A, 0x06, 0x02,
        0x0D, 0x09, 0x05, 0x01, 0x0C, 0x08, 0x04, 0x00);
    const __m256i vsh_fwd = _mm256_setr_epi64x(0, 2, 4, 6);
    const __m256i v_pack = _mm256_setr_epi32(0, 2, 4, 6, 0, 2, 4, 6);
    const __m256i v_phi_offsets = _mm256_setr_epi32(0, 1, 2, 3, 4, 5, 6, 7);
    const __m256i v_phi_mask = _mm256_set1_epi32(0xFFFFFF00);
    const int k_val = K;
    uint32_t sm = (k_val >= 16) ? 0xFFFFFFFF : (1U << (k_val * 2)) - 1;
    const __m256i vmask = _mm256_set1_epi32(sm);
    uint32_t last_min_abs_pos = 0xFFFFFFFF;
    const __m256i v_perm_prev = _mm256_setr_epi32(7, 0, 1, 2, 3, 4, 5, 6);

    for (size_t c = 0; c < num_chunks; ++c) {
      const uint8_t *fwd_ptr = (const uint8_t *)&seq_chunks[c];
#pragma unroll 19
      for (int i = 0; i < total_pos; i += 8) {
        int fwd_byte = i >> 2;
        __m256i fv1 = _mm256_castpd_si256(
            _mm256_broadcast_sd((const double *)(fwd_ptr + fwd_byte)));
        __m256i fv2 = _mm256_castpd_si256(
            _mm256_broadcast_sd((const double *)(fwd_ptr + fwd_byte + 1)));
        __m256i fs1 = _mm256_srlv_epi64(fv1, vsh_fwd);
        __m256i fs2 = _mm256_srlv_epi64(fv2, vsh_fwd);
        __m256i v_fwd_full =
            _mm256_blend_epi32(_mm256_permutevar8x32_epi32(fs1, v_pack),
                               _mm256_permutevar8x32_epi32(fs2, v_pack), 0xF0);
        __m256i v_fwd = _mm256_and_si256(v_fwd_full, vmask);

        __m256i v_h;
        if constexpr (CANONICAL) {
          __m256i vr = _mm256_shuffle_epi8(v_fwd_full, v_rev_bytes);
          __m256i lo = _mm256_and_si256(vr, v_m0f);
          __m256i hi = _mm256_and_si256(_mm256_srli_epi16(vr, 4), v_m0f);
          __m256i v_rev = _mm256_or_si256(_mm256_shuffle_epi8(v_lut_l, lo),
                                          _mm256_shuffle_epi8(v_lut_h, hi));
          int shift_amt = 32 - 2 * k_val;
          if (shift_amt > 0)
            v_rev = _mm256_srli_epi32(v_rev, shift_amt);
          v_rev = _mm256_and_si256(v_rev, vmask);
          v_h = _mm256_min_epu32(v_fwd, v_rev);
        } else
          v_h = v_fwd;
        v_h = _mm256_xor_si256(v_h, _mm256_srli_epi32(v_h, 13));
        v_h = _mm256_mullo_epi32(v_h, v_murmur_c);
        v_h = _mm256_xor_si256(v_h, _mm256_srli_epi32(v_h, 16));
        __m256i v_phi = _mm256_or_si256(
            _mm256_and_si256(v_h, v_phi_mask),
            _mm256_add_epi32(_mm256_set1_epi32(i), v_phi_offsets));
        _mm256_storeu_si256((__m256i *)&phi_buf[i], v_phi);
      }

#pragma unroll 16
      for (int batch = 0; batch < 128; batch += 8) {
        __m256i v_min_phi = reduce_phi<WindowSize>(phi_buf, batch);
        __m256i v_min_h = _mm256_and_si256(v_min_phi, v_phi_mask);
        __m256i v_local_idx =
            _mm256_and_si256(v_min_phi, _mm256_set1_epi32(0xFF));
        __m256i v_abs_p =
            _mm256_add_epi32(_mm256_set1_epi32(chunk_pos), v_local_idx);

        __m256i v_prev_abs_p = _mm256_blend_epi32(
            _mm256_permutevar8x32_epi32(v_abs_p, v_perm_prev),
            _mm256_set1_epi32(last_min_abs_pos), 0x01);
        __m256i v_is_new =
            _mm256_xor_si256(_mm256_cmpeq_epi32(v_abs_p, v_prev_abs_p),
                             _mm256_set1_epi32(0xFFFFFFFF));
        int mask = _mm256_movemask_ps((__m256)v_is_new) & 0xFF;

        __m256i v_shuf_mask = _mm256_loadu_si256(
            (const __m256i *)minimizer_shuffle_lut.shuffle_table[mask]);
        __m256i v_p_p = _mm256_permutevar8x32_epi32(v_abs_p, v_shuf_mask);

        _mm256_storeu_si256((__m256i *)&out_pos[total_out_idx], v_p_p);
        if (out_hash) {
            __m256i v_h_p = _mm256_permutevar8x32_epi32(v_min_h, v_shuf_mask);
            _mm256_storeu_si256((__m256i *)&out_hash[total_out_idx], v_h_p);
        }
        total_out_idx += _mm_popcnt_u32(mask);
        last_min_abs_pos =
            _mm_extract_epi32(_mm256_extractf128_si256(v_abs_p, 1), 3);
      }
      chunk_pos += 128;
    }
    return total_out_idx;
  }
};
} // namespace xseal
