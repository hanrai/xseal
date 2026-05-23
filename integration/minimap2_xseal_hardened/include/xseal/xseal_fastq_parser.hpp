#pragma once

#include <cstdint>
#include <cstring>
#include <immintrin.h>
#include <xseal/xseal_common.hpp>
#include <xseal/xseal_slot_pool.hpp>

namespace xseal {

template <typename BasesHandler, typename HeaderHandler, typename SkipHandler,
          uint32_t NumSlots = 64>
class XsealFastqParserT {
public:
  static constexpr uint64_t EMIT_THRESHOLD = 16 * 1024 * 1024; // 16MB

  XsealFastqParserT(BasesHandler b, HeaderHandler h, SkipHandler s)
      : on_bases(std::move(b)), on_header(std::move(h)), on_skip(std::move(s)) {
    init_luts();
  }

  void parse(const uint8_t *data, size_t size, uint64_t absolute_offset,
             const XsealSlot<NumSlots> &slot = XsealSlot<NumSlots>()) {
    current_slot = slot;
    if (size == 0)
      return;
    size_t ptr = 0;

    while (ptr < size) {
      switch (state) {
      case State::HEADER:
      case State::PLUS:
      case State::QUAL: {
        const uint8_t *base = data + ptr;
        size_t rem = size - ptr;
        const uint8_t *nl = (const uint8_t *)memchr(base, '\n', rem);
        if (nl) {
          size_t span = nl - base + 1;
          active_segment_logical_len += span;
          active_segment_physical_bytes += span;
          ptr += span;
          push_segment(false);
          advance_state();
          active_segment_start_offset = absolute_offset + ptr;
        } else {
          active_segment_logical_len += rem;
          active_segment_physical_bytes += rem;
          ptr += rem;
        }
        break;
      }

      case State::BASES: {
        while (ptr < size) {
          // 4-way SIMD unrolling (128 bytes)
          while (ptr + 128 <= size) {
            __m256i c0 = _mm256_loadu_si256((const __m256i *)(data + ptr));
            __m256i c1 = _mm256_loadu_si256((const __m256i *)(data + ptr + 32));
            __m256i c2 = _mm256_loadu_si256((const __m256i *)(data + ptr + 64));
            __m256i c3 = _mm256_loadu_si256((const __m256i *)(data + ptr + 96));

            __m256i cl0 = classify(c0);
            __m256i cl1 = classify(c1);
            __m256i cl2 = classify(c2);
            __m256i cl3 = classify(c3);

            __m256i stop0 = _mm256_cmpeq_epi8(
                _mm256_and_si256(cl0, v_class_base_nl), _mm256_setzero_si256());
            __m256i stop1 = _mm256_cmpeq_epi8(
                _mm256_and_si256(cl1, v_class_base_nl), _mm256_setzero_si256());
            __m256i stop2 = _mm256_cmpeq_epi8(
                _mm256_and_si256(cl2, v_class_base_nl), _mm256_setzero_si256());
            __m256i stop3 = _mm256_cmpeq_epi8(
                _mm256_and_si256(cl3, v_class_base_nl), _mm256_setzero_si256());

            __m256i combined_stop = _mm256_or_si256(
                _mm256_or_si256(stop0, stop1), _mm256_or_si256(stop2, stop3));

            if (__builtin_expect(
                    _mm256_testz_si256(combined_stop, combined_stop), 1)) {
              active_segment_logical_len += 128;
              active_segment_physical_bytes += 128;
              ptr += 128;
            } else {
              break;
            }
          }

          // 1-way SIMD fallback
          while (ptr + 32 <= size) {
            __m256i chunk = _mm256_loadu_si256((const __m256i *)(data + ptr));
            __m256i char_class = classify(chunk);

            __m256i is_stop =
                _mm256_cmpeq_epi8(_mm256_and_si256(char_class, v_class_base_nl),
                                  _mm256_setzero_si256());
            uint32_t stop_mask = (uint32_t)_mm256_movemask_epi8(is_stop);

            if (stop_mask == 0) {
              active_segment_logical_len += 32;
              active_segment_physical_bytes += 32;
              ptr += 32;
            } else {
              uint32_t off = __builtin_ctz(stop_mask);
              active_segment_logical_len += off;
              active_segment_physical_bytes += off;
              ptr += off;
              goto bases_transition;
            }
          }
          if (ptr >= size)
            break;

        bases_transition:
          uint8_t c_class = scalar_class[data[ptr]];
          if (c_class == 0x01) {
            active_segment_logical_len++;
            active_segment_physical_bytes++;
            ptr++;
          } else if (c_class == 0x02) {
            active_segment_physical_bytes++;
            ptr++;
            if (data[ptr - 1] == '\n') {
              push_segment(false);
              advance_state();
              active_segment_start_offset = absolute_offset + ptr;
              break; // out of BASES state
            }
          } else { // N or Ambiguous
            push_segment(false);
            state = State::N_GAP;
            current_type = SegType::N_GAP;
            active_segment_start_offset = absolute_offset + ptr;
            active_segment_logical_len = 0;
            active_segment_physical_bytes = 0;
            break;
          }
        }
        break;
      }

      case State::N_GAP: {
        while (ptr < size) {
          // 1-way SIMD fallback for N_GAP
          while (ptr + 32 <= size) {
            __m256i chunk = _mm256_loadu_si256((const __m256i *)(data + ptr));
            __m256i char_class = classify(chunk);

            __m256i is_stop = _mm256_and_si256(
                char_class, _mm256_set1_epi8(0x03)); // Base or NL
            uint32_t stop_mask = (uint32_t)_mm256_movemask_epi8(
                _mm256_cmpgt_epi8(is_stop, _mm256_setzero_si256()));

            if (stop_mask == 0) {
              active_segment_logical_len += 32;
              active_segment_physical_bytes += 32;
              ptr += 32;
            } else {
              uint32_t off = __builtin_ctz(stop_mask);
              active_segment_logical_len += off;
              active_segment_physical_bytes += off;
              ptr += off;
              goto gap_transition;
            }
          }
          if (ptr >= size)
            break;
        gap_transition:
          uint8_t c_class = scalar_class[data[ptr]];
          if (c_class == 0x08) {
            active_segment_logical_len++;
            active_segment_physical_bytes++;
            ptr++;
          } else if (c_class == 0x01) {
            push_segment(false);
            state = State::BASES;
            current_type = SegType::BASES;
            active_segment_start_offset = absolute_offset + ptr;
            active_segment_logical_len = 0;
            active_segment_physical_bytes = 0;
            break;
          } else if (c_class == 0x02) { // \n or \r
            active_segment_physical_bytes++;
            ptr++;
            if (data[ptr - 1] == '\n') {
              push_segment(false);
              advance_state(); // N_GAP to PLUS
              active_segment_start_offset = absolute_offset + ptr;
              break;
            }
          }
        }
        break;
      }
      }
    }

    if (active_segment_physical_bytes > 0) {
      push_segment(true);
      active_segment_start_offset = absolute_offset + size;
    }
  }

  void finish() { push_segment(false); }

private:
  void init_luts() {
    v_0F = _mm256_set1_epi8(0x0F);
    v_class_base_nl = _mm256_set1_epi8(0x01);
    v_class_base = _mm256_set1_epi8(0x01);
    v_class_nl = _mm256_set1_epi8(0x02);
    v_class_gap = _mm256_set1_epi8(0x08);

    for (int i = 0; i < 256; ++i)
      scalar_class[i] = 0x08;
    scalar_class['\n'] = 0x02;
    scalar_class['\r'] = 0x02;
    const char *bases = "ACGTacgt";
    for (const char *p = bases; *p; ++p)
      scalar_class[(uint8_t)*p] = 0x01;

    uint8_t low_tbl[16] = {0};
    uint8_t high_tbl[16] = {0};
    for (int i = 0; i < 256; ++i) {
      uint8_t cls = scalar_class[i];
      low_tbl[i & 0x0F] |= cls;
      high_tbl[i >> 4] |= cls;
    }

    uint8_t low_32[32], high_32[32];
    for (int i = 0; i < 16; ++i) {
      low_32[i] = low_32[i + 16] = low_tbl[i];
      high_32[i] = high_32[i + 16] = high_tbl[i];
    }
    lut_low = _mm256_loadu_si256((const __m256i *)low_32);
    lut_high = _mm256_loadu_si256((const __m256i *)high_32);
  }

  inline __m256i classify(__m256i v) {
    __m256i low = _mm256_and_si256(v, v_0F);
    __m256i high = _mm256_and_si256(_mm256_srli_epi16(v, 4), v_0F);
    return _mm256_and_si256(_mm256_shuffle_epi8(lut_low, low),
                            _mm256_shuffle_epi8(lut_high, high));
  }

  void advance_state() {
    // N_GAP implicitly advances to PLUS when hitting \n
    if (state == State::HEADER) {
      state = State::BASES;
      current_type = SegType::BASES;
    } else if (state == State::BASES || state == State::N_GAP) {
      state = State::PLUS;
      current_type = SegType::SKIP;
    } else if (state == State::PLUS) {
      state = State::QUAL;
      current_type = SegType::SKIP;
    } else if (state == State::QUAL) {
      state = State::HEADER;
      current_type = SegType::HEADER;
    }
  }

  void push_segment(bool is_physical) {
    if (active_segment_logical_len > 0 || current_type == SegType::HEADER ||
        current_type == SegType::SKIP) {
      uint32_t log_len =
          (current_type == SegType::BASES || current_type == SegType::N_GAP)
              ? (uint32_t)active_segment_logical_len
              : 1;
      SeqSegment seg = {active_segment_start_offset, log_len,
                        (uint32_t)active_segment_physical_bytes, current_type,
                        is_physical};

      if (current_type == SegType::BASES || current_type == SegType::N_GAP) {
        if (current_type == SegType::BASES)
          on_bases(seg, current_slot);
        else
          on_skip(seg, current_slot);
      } else if (current_type == SegType::HEADER)
        on_header(seg, current_slot);
      else
        on_skip(seg, current_slot);
    }
    active_segment_logical_len = 0;
    active_segment_physical_bytes = 0;
  }

  BasesHandler on_bases;
  HeaderHandler on_header;
  SkipHandler on_skip;
  XsealSlot<NumSlots> current_slot;

  enum class State { HEADER, BASES, N_GAP, PLUS, QUAL } state = State::HEADER;
  SegType current_type = SegType::HEADER;
  uint64_t active_segment_start_offset = 0;
  uint64_t active_segment_logical_len = 0;
  uint64_t active_segment_physical_bytes = 0;

  __m256i v_0F;
  __m256i v_class_base_nl;
  __m256i v_class_base;
  __m256i v_class_nl;
  __m256i v_class_gap;
  __m256i lut_low;
  __m256i lut_high;
  uint8_t scalar_class[256];
};

} // namespace xseal
