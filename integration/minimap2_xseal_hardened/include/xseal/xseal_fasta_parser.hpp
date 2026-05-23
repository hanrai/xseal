#pragma once

#include <cstdint>
#include <cstring>
#include <immintrin.h>
#include <vector>

#include <utility>
#include <xseal/xseal_common.hpp>
#include <xseal/xseal_slot_pool.hpp>
namespace xseal {

using XsealRunTable = std::vector<SeqSegment>;

/**
 * @brief 核心 FASTA 扫描引擎
 * 负责识别序列段并生成结构表（Run-Length Table）。
 * 无模板，极致性能，易于扩展 QC 和 Profiling 功能。
 */
class XsealFastaParser {
public:
    XsealFastaParser() { init_luts(); }

    void parse(const uint8_t *data, size_t size, uint64_t absolute_offset,
              XsealRunTable &run_table) {
        if (size == 0) return;
        size_t ptr = 0;

        while (ptr < size) {
            switch (state) {
                case State::HEADER:
                    handle_header(data, size, ptr, absolute_offset, run_table);
                    break;
                case State::BASES:
                    handle_bases(data, size, ptr, absolute_offset, run_table);
                    break;
                case State::N_GAP:
                    handle_gap(data, size, ptr, absolute_offset, run_table);
                    break;
            }
        }

        // 处理块末尾可能被截断的物理段
        if (active_segment_physical_bytes > 0) {
            push_segment(run_table, true);
            active_segment_start_offset = absolute_offset + size;
        }
    }

    void finish(XsealRunTable &run_table) {
        push_segment(run_table, false);
    }

    uint32_t get_max_id_len() const { return max_id_len; }

private:
    enum class State { HEADER, BASES, N_GAP } state = State::HEADER;
    SegType current_type = SegType::HEADER;
    uint64_t active_segment_start_offset = 0;
    uint64_t active_segment_logical_len = 0;
    uint64_t active_segment_physical_bytes = 0;

    // Profiler metrics
    uint32_t max_id_len = 0;
    bool in_header_id = false;
    uint32_t current_id_scan_len = 0;

    inline void handle_header(const uint8_t *data, size_t size, size_t &ptr,
                              uint64_t absolute_offset, XsealRunTable &run_table) {
        // 如果物理字节为 0 或 1(仅有'>')，说明这是一个新 Header 的开始
        if (active_segment_physical_bytes <= 1) {
            in_header_id = true;
            current_id_scan_len = 0;
        }

        while (ptr < size) {
            const uint8_t *base = data + ptr;
            size_t rem = size - ptr;
            const uint8_t *nl = (const uint8_t *)memchr(base, '\n', rem);
            
            size_t span = nl ? (size_t)(nl - base) : rem;

            if (in_header_id) {
                size_t i = 0;
                // 如果是本段的最开始，且第一个字符是 '>'，则跳过它不计入 ID 长度
                if (active_segment_physical_bytes == 0 && span > 0 && base[0] == '>') {
                    i = 1;
                }
                
                for (; i < span; ++i) {
                    if (base[i] <= ' ') { 
                        in_header_id = false;
                        if (current_id_scan_len > max_id_len) max_id_len = current_id_scan_len;
                        break;
                    }
                    current_id_scan_len++;
                }
                // 如果到了行尾还没遇到空格，此时 ID 结束（或者本段 ID 扫描暂告一段落）
                if (nl && in_header_id) {
                    if (current_id_scan_len > max_id_len) max_id_len = current_id_scan_len;
                    in_header_id = false;
                }
            }

            if (nl) {
                active_segment_logical_len += span;
                active_segment_physical_bytes += span;
                ptr += span;

                if (span > 0 && base[span - 1] == '\r') {
                    active_segment_logical_len--;
                }

                push_segment(run_table, false);
                state = State::BASES;
                current_type = SegType::BASES;
                ptr++; // Skip '\n'
                active_segment_start_offset = absolute_offset + ptr;
                active_segment_logical_len = 0;
                active_segment_physical_bytes = 0;
                return;
            } else {
                active_segment_logical_len += rem;
                active_segment_physical_bytes += rem;
                ptr += rem;
                return;
            }
        }
    }

    inline void handle_bases(const uint8_t *data, size_t size, size_t &ptr,
                             uint64_t absolute_offset, std::vector<SeqSegment> &run_table) {
        while (ptr < size) {
            while (ptr + 128 <= size) {
                __m256i c0 = _mm256_loadu_si256((const __m256i *)(data + ptr));
                __m256i c1 = _mm256_loadu_si256((const __m256i *)(data + ptr + 32));
                __m256i c2 = _mm256_loadu_si256((const __m256i *)(data + ptr + 64));
                __m256i c3 = _mm256_loadu_si256((const __m256i *)(data + ptr + 96));

                __m256i cl0 = classify(c0);
                __m256i cl1 = classify(c1);
                __m256i cl2 = classify(c2);
                __m256i cl3 = classify(c3);

                __m256i stop0 = _mm256_cmpeq_epi8(_mm256_and_si256(cl0, v_class_base_nl), _mm256_setzero_si256());
                __m256i stop1 = _mm256_cmpeq_epi8(_mm256_and_si256(cl1, v_class_base_nl), _mm256_setzero_si256());
                __m256i stop2 = _mm256_cmpeq_epi8(_mm256_and_si256(cl2, v_class_base_nl), _mm256_setzero_si256());
                __m256i stop3 = _mm256_cmpeq_epi8(_mm256_and_si256(cl3, v_class_base_nl), _mm256_setzero_si256());

                __m256i combined_stop = _mm256_or_si256(_mm256_or_si256(stop0, stop1),
                                                        _mm256_or_si256(stop2, stop3));

                if (__builtin_expect(_mm256_testz_si256(combined_stop, combined_stop), 1)) {
                    uint32_t nl0 = _mm256_movemask_epi8(_mm256_cmpeq_epi8(_mm256_and_si256(cl0, v_class_nl), v_class_nl));
                    uint32_t nl1 = _mm256_movemask_epi8(_mm256_cmpeq_epi8(_mm256_and_si256(cl1, v_class_nl), v_class_nl));
                    uint32_t nl2 = _mm256_movemask_epi8(_mm256_cmpeq_epi8(_mm256_and_si256(cl2, v_class_nl), v_class_nl));
                    uint32_t nl3 = _mm256_movemask_epi8(_mm256_cmpeq_epi8(_mm256_and_si256(cl3, v_class_nl), v_class_nl));

                    active_segment_logical_len += 128 - (__builtin_popcount(nl0) + __builtin_popcount(nl1) +
                                                         __builtin_popcount(nl2) + __builtin_popcount(nl3));
                    active_segment_physical_bytes += 128;
                    ptr += 128;
                } else {
                    break;
                }
            }

            while (ptr + 32 <= size) {
                __m256i chunk = _mm256_loadu_si256((const __m256i *)(data + ptr));
                __m256i char_class = classify(chunk);
                __m256i is_stop = _mm256_cmpeq_epi8(_mm256_and_si256(char_class, v_class_base_nl), _mm256_setzero_si256());
                uint32_t stop_mask = (uint32_t)_mm256_movemask_epi8(is_stop);

                if (stop_mask == 0) {
                    uint32_t nl_mask = _mm256_movemask_epi8(_mm256_cmpeq_epi8(_mm256_and_si256(char_class, v_class_nl), v_class_nl));
                    active_segment_logical_len += 32 - __builtin_popcount(nl_mask);
                    active_segment_physical_bytes += 32;
                    ptr += 32;
                } else {
                    uint32_t off = __builtin_ctz(stop_mask);
                    uint32_t nl_mask = _mm256_movemask_epi8(_mm256_cmpeq_epi8(_mm256_and_si256(char_class, v_class_nl), v_class_nl)) & ((1u << off) - 1);
                    active_segment_logical_len += off - __builtin_popcount(nl_mask);
                    active_segment_physical_bytes += off;
                    ptr += off;
                    break; 
                }
            }
            if (ptr >= size) return;

            // Transition to scalar logic
            uint8_t c_class = scalar_class[data[ptr]];
            if (c_class == 0x01) {
                active_segment_logical_len++;
                active_segment_physical_bytes++;
                ptr++;
            } else if (c_class == 0x02) {
                active_segment_physical_bytes++;
                ptr++;
            } else if (c_class == 0x04) { // '>'
                push_segment(run_table, false);
                state = State::HEADER;
                current_type = SegType::HEADER;
                active_segment_start_offset = absolute_offset + ptr;
                active_segment_logical_len = 1;
                active_segment_physical_bytes = 1;
                ptr++;
                return;
            } else { // N or Ambiguous
                push_segment(run_table, false);
                state = State::N_GAP;
                current_type = SegType::N_GAP;
                active_segment_start_offset = absolute_offset + ptr;
                active_segment_logical_len = 0;
                active_segment_physical_bytes = 0;
                return;
            }
        }
    }

    inline void handle_gap(const uint8_t *data, size_t size, size_t &ptr,
                           uint64_t absolute_offset, std::vector<SeqSegment> &run_table) {
        while (ptr < size) {
            while (ptr + 128 <= size) {
                __m256i c0 = _mm256_loadu_si256((const __m256i *)(data + ptr));
                __m256i c1 = _mm256_loadu_si256((const __m256i *)(data + ptr + 32));
                __m256i c2 = _mm256_loadu_si256((const __m256i *)(data + ptr + 64));
                __m256i c3 = _mm256_loadu_si256((const __m256i *)(data + ptr + 96));

                __m256i cl0 = classify(c0);
                __m256i cl1 = classify(c1);
                __m256i cl2 = classify(c2);
                __m256i cl3 = classify(c3);

                __m256i stop0 = _mm256_and_si256(cl0, _mm256_set1_epi8(0x05));
                __m256i stop1 = _mm256_and_si256(cl1, _mm256_set1_epi8(0x05));
                __m256i stop2 = _mm256_and_si256(cl2, _mm256_set1_epi8(0x05));
                __m256i stop3 = _mm256_and_si256(cl3, _mm256_set1_epi8(0x05));

                __m256i combined_stop = _mm256_or_si256(_mm256_or_si256(stop0, stop1),
                                                        _mm256_or_si256(stop2, stop3));

                if (__builtin_expect(_mm256_testz_si256(combined_stop, combined_stop), 1)) {
                    uint32_t nl0 = _mm256_movemask_epi8(_mm256_cmpeq_epi8(_mm256_and_si256(cl0, v_class_nl), v_class_nl));
                    uint32_t nl1 = _mm256_movemask_epi8(_mm256_cmpeq_epi8(_mm256_and_si256(cl1, v_class_nl), v_class_nl));
                    uint32_t nl2 = _mm256_movemask_epi8(_mm256_cmpeq_epi8(_mm256_and_si256(cl2, v_class_nl), v_class_nl));
                    uint32_t nl3 = _mm256_movemask_epi8(_mm256_cmpeq_epi8(_mm256_and_si256(cl3, v_class_nl), v_class_nl));

                    active_segment_logical_len += 128 - (__builtin_popcount(nl0) + __builtin_popcount(nl1) +
                                                         __builtin_popcount(nl2) + __builtin_popcount(nl3));
                    active_segment_physical_bytes += 128;
                    ptr += 128;
                } else {
                    break;
                }
            }

            while (ptr + 32 <= size) {
                __m256i chunk = _mm256_loadu_si256((const __m256i *)(data + ptr));
                __m256i char_class = classify(chunk);
                __m256i is_stop = _mm256_and_si256(char_class, _mm256_set1_epi8(0x05));
                uint32_t stop_mask = (uint32_t)_mm256_movemask_epi8(_mm256_cmpgt_epi8(is_stop, _mm256_setzero_si256()));

                if (stop_mask == 0) {
                    uint32_t nl_mask = _mm256_movemask_epi8(_mm256_cmpeq_epi8(_mm256_and_si256(char_class, v_class_nl), v_class_nl));
                    active_segment_logical_len += 32 - __builtin_popcount(nl_mask);
                    active_segment_physical_bytes += 32;
                    ptr += 32;
                } else {
                    uint32_t off = __builtin_ctz(stop_mask);
                    uint32_t nl_mask = _mm256_movemask_epi8(_mm256_cmpeq_epi8(_mm256_and_si256(char_class, v_class_nl), v_class_nl)) & ((1u << off) - 1);
                    active_segment_logical_len += off - __builtin_popcount(nl_mask);
                    active_segment_physical_bytes += off;
                    ptr += off;
                    break;
                }
            }
            if (ptr >= size) return;

            // Transition to scalar logic
            uint8_t c_class = scalar_class[data[ptr]];
            if (c_class == 0x08) {
                active_segment_logical_len++;
                active_segment_physical_bytes++;
                ptr++;
            } else if (c_class == 0x02) {
                active_segment_physical_bytes++;
                ptr++;
            } else if (c_class == 0x01) {
                push_segment(run_table, false);
                state = State::BASES;
                current_type = SegType::BASES;
                active_segment_start_offset = absolute_offset + ptr;
                active_segment_logical_len = 0;
                active_segment_physical_bytes = 0;
                return;
            } else { // '>'
                push_segment(run_table, false);
                state = State::HEADER;
                current_type = SegType::HEADER;
                active_segment_start_offset = absolute_offset + ptr;
                active_segment_logical_len = 1;
                active_segment_physical_bytes = 1;
                ptr++;
                return;
            }
        }
    }

    void init_luts() {
        v_0F = _mm256_set1_epi8(0x0F);
        v_class_base_nl = _mm256_set1_epi8(0x03);
        v_class_base = _mm256_set1_epi8(0x01);
        v_class_nl = _mm256_set1_epi8(0x02);
        v_class_header = _mm256_set1_epi8(0x04);
        v_class_gap = _mm256_set1_epi8(0x08);

        for (int i = 0; i < 256; ++i) scalar_class[i] = 0x08;
        scalar_class['>'] = 0x04;
        scalar_class['\n'] = 0x02;
        scalar_class['\r'] = 0x02;
        const char *bases = "ACGTacgt";
        for (const char *p = bases; *p; ++p) scalar_class[(uint8_t)*p] = 0x01;

        uint8_t low_tbl[16] = {0}, high_tbl[16] = {0};
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

    void push_segment(XsealRunTable &run_table, bool is_physical) {
        if (active_segment_logical_len > 0 || current_type == SegType::HEADER) {
            run_table.push_back({active_segment_start_offset, (uint32_t)active_segment_logical_len,
                                 (uint32_t)active_segment_physical_bytes, current_type, is_physical});
        }
        active_segment_logical_len = 0;
        active_segment_physical_bytes = 0;
    }

    __m256i v_0F, v_class_base_nl, v_class_base, v_class_nl, v_class_header, v_class_gap;
    __m256i lut_low, lut_high;
    uint8_t scalar_class[256];
};

/**
 * @brief 兼容层：将基于表的扫描转换为原有的回调模式
 */
template <typename BasesHandler, typename HeaderHandler, typename GapHandler, uint32_t NumSlots = 64>
class XsealFastaParserT {
public:
    XsealFastaParserT(BasesHandler b, HeaderHandler h, GapHandler g)
        : on_bases(std::move(b)), on_header(std::move(h)), on_gap(std::move(g)) {
        run_table_cache.reserve(1024); // 预分配一些空间
    }

    void parse(const uint8_t *data, size_t size, uint64_t absolute_offset,
              const XsealSlot<NumSlots> &slot = XsealSlot<NumSlots>()) {
        run_table_cache.clear();
        scanner.parse(data, size, absolute_offset, run_table_cache);

        // 兼容层：模拟回调
        for (const auto &seg : run_table_cache) {
            if (seg.type == SegType::BASES)
                on_bases(seg, slot);
            else if (seg.type == SegType::HEADER)
                on_header(seg, slot);
            else
                on_gap(seg, slot);
        }
    }

    void finish() {
        run_table_cache.clear();
        scanner.finish(run_table_cache);
        for (const auto &seg : run_table_cache) {
            if (seg.type == SegType::BASES)
                on_bases(seg, XsealSlot<NumSlots>());
            else if (seg.type == SegType::HEADER)
                on_header(seg, XsealSlot<NumSlots>());
            else
                on_gap(seg, XsealSlot<NumSlots>());
        }
    }

    uint32_t get_max_id_len() const { return scanner.get_max_id_len(); }

    private:
    XsealFastaParser scanner;
    XsealRunTable run_table_cache;
    BasesHandler on_bases;
    HeaderHandler on_header;
    GapHandler on_gap;
};

} // namespace xseal
