// SPDX-License-Identifier: Apache-2.0
// Copyright (c) 2026 hanrai. All Rights Reserved.

#pragma once

#include <cstdint>
#include <cstring>
#include <immintrin.h>

namespace xseal {

enum XSealEncoderStatus {
    XSEAL_ENC_OK          = 0,
    XSEAL_ENC_BUFFER_FULL = 1,  
    XSEAL_ENC_EOF         = 2
};

struct XSealEncoderState {
    uint8_t* out_buffer;
    size_t   out_byte_idx;
    uint64_t bit_buffer;
    int      bit_count;
    size_t   out_capacity;
    size_t   total_bases_encoded; 
};

class XSealEncoder {
private:
    uint64_t pext_mask_lut[256];
    int bit_len_lut[256];
    uint8_t seq_nt4_table[256];

public:
    XSealEncoder() {
        std::memset(seq_nt4_table, 0, 256);
        seq_nt4_table[(uint8_t)'A'] = 0; seq_nt4_table[(uint8_t)'C'] = 1;
        seq_nt4_table[(uint8_t)'G'] = 2; seq_nt4_table[(uint8_t)'T'] = 3;
        seq_nt4_table[(uint8_t)'a'] = 0; seq_nt4_table[(uint8_t)'c'] = 1;
        seq_nt4_table[(uint8_t)'g'] = 2; seq_nt4_table[(uint8_t)'t'] = 3;

        for (int i = 0; i < 256; i++) {
            uint64_t mask = 0; int bits = 0;
            for (int j = 0; j < 8; j++) {
                if ((i >> j) & 1) { mask |= (3ULL << (j * 8)); bits += 2; }
            }
            pext_mask_lut[i] = mask; bit_len_lut[i] = bits;
        }
    }

    void init_state(XSealEncoderState* state, uint8_t* target_buffer, size_t capacity) {
        state->out_buffer = target_buffer; state->out_byte_idx = 0;
        state->bit_buffer = 0; state->bit_count = 0;
        state->out_capacity = capacity; state->total_bases_encoded = 0;
    }

    void flush_state(XSealEncoderState* state) {
        if (state->bit_count > 0) {
            int bytes = (state->bit_count + 7) / 8;
            std::memcpy(state->out_buffer + state->out_byte_idx, &state->bit_buffer, bytes);
            state->out_byte_idx += bytes;
            state->bit_count = 0; state->bit_buffer = 0;
        }
    }

    XSealEncoderStatus encode_chunk(XSealEncoderState* state, const char** data_in_out, const char* end) {
        const char* ptr = *data_in_out;
        uint8_t* local_out_buf = state->out_buffer; size_t local_out_idx = state->out_byte_idx;
        uint64_t local_bit_buf = state->bit_buffer; int local_bit_count = state->bit_count;
        size_t local_cap = state->out_capacity; size_t local_bases = state->total_bases_encoded;

        // V5.0 Constant Engine (Ultra Slim)
        const __m256i v_lookup = _mm256_setr_epi8(
            0x80,0x00,0x80,0x01,0x03,0x80,0x80,0x02,0x80,0x80,0x80,0x80,0x80,0x80,0x80,0x80,
            0x80,0x00,0x80,0x01,0x03,0x80,0x80,0x02,0x80,0x80,0x80,0x80,0x80,0x80,0x80,0x80);
        const __m256i v_mask_0F = _mm256_set1_epi8(0x0F);
        const __m256i v_nt4_lut = _mm256_setr_epi8(0, 0, 0, 1, 3, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 3, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0);
        const __m256i v_mul1 = _mm256_set1_epi32(0x40100401); const __m256i v_mul2 = _mm256_set1_epi32(0x00010001);
        const __m256i v_shuf_mask = _mm256_setr_epi8(0,4,8,12,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,0,4,8,12,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1);
        const __m256i v_nl = _mm256_set1_epi8('\n');
        const __m256i v_cr = _mm256_set1_epi8('\r');

        XSealEncoderStatus final_status = XSEAL_ENC_EOF;

        // V6.0 Turbo Engine (4-way Unrolled Madden)
        while (ptr + 128 <= end) {
            __m256i v0 = _mm256_loadu_si256((const __m256i*)ptr);
            __m256i v1 = _mm256_loadu_si256((const __m256i*)(ptr + 32));
            __m256i v2 = _mm256_loadu_si256((const __m256i*)(ptr + 64));
            __m256i v3 = _mm256_loadu_si256((const __m256i*)(ptr + 96));

            __m256i nl0 = _mm256_or_si256(_mm256_cmpeq_epi8(v0, v_nl), _mm256_cmpeq_epi8(v0, v_cr));
            __m256i nl1 = _mm256_or_si256(_mm256_cmpeq_epi8(v1, v_nl), _mm256_cmpeq_epi8(v1, v_cr));
            __m256i nl2 = _mm256_or_si256(_mm256_cmpeq_epi8(v2, v_nl), _mm256_cmpeq_epi8(v2, v_cr));
            __m256i nl3 = _mm256_or_si256(_mm256_cmpeq_epi8(v3, v_nl), _mm256_cmpeq_epi8(v3, v_cr));
            
            __m256i combined_nl = _mm256_or_si256(_mm256_or_si256(nl0, nl1), _mm256_or_si256(nl2, nl3));
            
            if (__builtin_expect(_mm256_testz_si256(combined_nl, combined_nl), 1)) {
                if (__builtin_expect(local_out_idx + 32 > local_cap, 0)) break; 

                // 4-way Madden Pipeline
                __m256i t0 = _mm256_shuffle_epi8(v_nt4_lut, v0);
                __m256i t1 = _mm256_shuffle_epi8(v_nt4_lut, v1);
                __m256i t2 = _mm256_shuffle_epi8(v_nt4_lut, v2);
                __m256i t3 = _mm256_shuffle_epi8(v_nt4_lut, v3);

                t0 = _mm256_maddubs_epi16(t0, v_mul1); t1 = _mm256_maddubs_epi16(t1, v_mul1);
                t2 = _mm256_maddubs_epi16(t2, v_mul1); t3 = _mm256_maddubs_epi16(t3, v_mul1);

                t0 = _mm256_madd_epi16(t0, v_mul2); t1 = _mm256_madd_epi16(t1, v_mul2);
                t2 = _mm256_madd_epi16(t2, v_mul2); t3 = _mm256_madd_epi16(t3, v_mul2);

                t0 = _mm256_shuffle_epi8(t0, v_shuf_mask); t1 = _mm256_shuffle_epi8(t1, v_shuf_mask);
                t2 = _mm256_shuffle_epi8(t2, v_shuf_mask); t3 = _mm256_shuffle_epi8(t3, v_shuf_mask);

                auto pack_and_store = [&](__m256i v) {
                    uint64_t bits = (uint32_t)_mm_cvtsi128_si32(_mm256_castsi256_si128(v)) | ((uint64_t)(uint32_t)_mm_cvtsi128_si32(_mm256_extracti128_si256(v, 1)) << 32);
                    *((uint64_t*)(local_out_buf + local_out_idx)) = local_bit_buf | (bits << local_bit_count);
                    local_out_idx += 8;
                    local_bit_buf = local_bit_count ? (bits >> (64 - local_bit_count)) : 0;
                };

                pack_and_store(t0); pack_and_store(t1); pack_and_store(t2); pack_and_store(t3);
                local_bases += 128; ptr += 128;
            } else break; 
        }

        while (ptr + 32 <= end) {
            __m256i v_seq = _mm256_loadu_si256((const __m256i*)ptr);
            __m256i is_nl = _mm256_cmpeq_epi8(v_seq, v_nl);
            __m256i is_cr = _mm256_cmpeq_epi8(v_seq, v_cr);
            __m256i any_newline_v = _mm256_or_si256(is_nl, is_cr);
            uint32_t nl_mask = _mm256_movemask_epi8(any_newline_v);

            if (__builtin_expect(nl_mask == 0, 1)) {
                // Madden Path (God Speed): 100% ACGT
                if (__builtin_expect(local_out_idx + 8 > local_cap, 0)) { final_status = XSEAL_ENC_BUFFER_FULL; goto exit_encode; }
                __m256i trans = _mm256_shuffle_epi8(v_nt4_lut, v_seq);
                __m256i t1 = _mm256_maddubs_epi16(trans, v_mul1);
                __m256i t2 = _mm256_madd_epi16(t1, v_mul2);
                __m256i t3 = _mm256_shuffle_epi8(t2, v_shuf_mask);
                uint32_t lo32 = _mm_cvtsi128_si32(_mm256_castsi256_si128(t3));
                uint32_t hi32 = _mm_cvtsi128_si32(_mm256_extracti128_si256(t3, 1));
                uint64_t chunk_bits = lo32 | ((uint64_t)hi32 << 32);

                *((uint64_t*)(local_out_buf + local_out_idx)) = local_bit_buf | (chunk_bits << local_bit_count);
                local_out_idx += 8;
                local_bit_buf = local_bit_count ? (chunk_bits >> (64 - local_bit_count)) : 0;
                local_bases += 32; ptr += 32;
            } else {
                // PEXT Path (Elite Speed): Contains newlines, but guaranteed only bases otherwise
                uint32_t valid_mask = ~nl_mask;
                __m256i lo = _mm256_and_si256(v_seq, v_mask_0F);
                __m256i v_mapped = _mm256_shuffle_epi8(v_lookup, lo);
                
                uint64_t ext0 = _mm256_extract_epi64(v_mapped, 0); uint64_t ext1 = _mm256_extract_epi64(v_mapped, 1);
                uint64_t ext2 = _mm256_extract_epi64(v_mapped, 2); uint64_t ext3 = _mm256_extract_epi64(v_mapped, 3);
                uint8_t m0 = valid_mask & 0xFF; uint8_t m1 = (valid_mask >> 8) & 0xFF;
                uint8_t m2 = (valid_mask >> 16) & 0xFF; uint8_t m3 = (valid_mask >> 24) & 0xFF;
                uint64_t p0 = _pext_u64(ext0, pext_mask_lut[m0]); uint64_t p1 = _pext_u64(ext1, pext_mask_lut[m1]);
                uint64_t p2 = _pext_u64(ext2, pext_mask_lut[m2]); uint64_t p3 = _pext_u64(ext3, pext_mask_lut[m3]);
                int c0 = bit_len_lut[m0]; int c1 = bit_len_lut[m1]; int c2 = bit_len_lut[m2]; int c3 = bit_len_lut[m3];
                uint64_t combined = p0 | (p1 << c0) | (p2 << (c0 + c1)) | (p3 << (c0 + c1 + c2));
                int total_bits = c0 + c1 + c2 + c3;

                local_bit_buf |= (combined << local_bit_count);
                int new_count = local_bit_count + total_bits;
                if (new_count >= 64) {
                    if (__builtin_expect(local_out_idx + 8 > local_cap, 0)) { final_status = XSEAL_ENC_BUFFER_FULL; goto exit_encode; }
                    *((uint64_t*)(local_out_buf + local_out_idx)) = local_bit_buf;
                    local_out_idx += 8; local_bit_count = new_count - 64;
                    local_bit_buf = (local_bit_count > 0) ? (combined >> (total_bits - local_bit_count)) : 0;
                } else { local_bit_count = new_count; }
                local_bases += (total_bits / 2); ptr += 32;
            }
        }

        while (ptr < end) {
            if (__builtin_expect(local_out_idx + 8 > local_cap, 0)) { final_status = XSEAL_ENC_BUFFER_FULL; goto exit_encode; }
            char c = *ptr;
            if (c == '\n' || c == '\r') { ptr++; continue; }
            
            uint8_t val = seq_nt4_table[(uint8_t)c];
            local_bit_buf |= ((uint64_t)val << local_bit_count);
            local_bit_count += 2;
            if (local_bit_count >= 64) {
                *((uint64_t*)(local_out_buf + local_out_idx)) = local_bit_buf;
                local_out_idx += 8; local_bit_buf = 0; local_bit_count = 0;
            }
            local_bases++; ptr++;
        }

    exit_encode:
        state->out_byte_idx = local_out_idx; state->bit_buffer = local_bit_buf;
        state->bit_count = local_bit_count; state->total_bases_encoded = local_bases;
        *data_in_out = ptr; return final_status;
    }
};

} // namespace xseal
