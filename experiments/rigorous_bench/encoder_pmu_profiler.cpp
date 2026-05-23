// SPDX-License-Identifier: Apache-2.0
// Copyright (c) 2026 hanrai. All Rights Reserved.

#include <iostream>
#include <vector>
#include <chrono>
#include <iomanip>
#include <numeric>
#include <algorithm>
#include <cstring>
#include <cmath>
#include <x86intrin.h>
#include "xseal/xseal_encoder.hpp"

using namespace xseal;

unsigned char minimap2_nt4_table[256] = {
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  3, 3, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  3, 3, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
};

size_t minimap2_encode(const char* start, const char* end, uint8_t* out) {
    size_t out_idx = 0;
    uint64_t bit_buf = 0;
    int bit_cnt = 0;
    for (const char* p = start; p < end; ++p) {
        uint8_t c = minimap2_nt4_table[(uint8_t)*p];
        if (c < 4) {
            bit_buf |= ((uint64_t)c << bit_cnt);
            bit_cnt += 2;
            if (bit_cnt >= 64) {
                std::memcpy(out + out_idx, &bit_buf, 8);
                out_idx += 8; bit_buf = 0; bit_cnt = 0;
            }
        }
    }
    if (bit_cnt > 0) std::memcpy(out + out_idx, &bit_buf, (bit_cnt + 7) / 8);
    return out_idx;
}

int main(int argc, char** argv) {
    if (argc < 3) {
        std::cerr << "Usage: " << argv[0] << " <size_mb> <mode: minimap2|xseal>\n";
        return 1;
    }
    size_t data_size = std::stoull(argv[1]) * 1024 * 1024;
    std::string mode = argv[2];

    char* input = (char*)_mm_malloc(data_size, 64);
    uint8_t* output = (uint8_t*)_mm_malloc(data_size / 4 + 1024, 64);

    const char bases[] = "ACGT";
    for (size_t i = 0; i < data_size; ++i) input[i] = bases[rand() % 4];

    if (mode == "minimap2") {
        for (int i = 0; i < 10; ++i) {
            minimap2_encode(input, input + data_size, output);
        }
    } else if (mode == "xseal") {
        XSealEncoder encoder;
        for (int i = 0; i < 10; ++i) {
            XSealEncoderState state;
            encoder.init_state(&state, output, data_size / 4 + 1024);
            const char* ptr = input;
            encoder.encode_chunk(&state, &ptr, input + data_size);
            encoder.flush_state(&state);
        }
    }

    _mm_free(input); _mm_free(output);
    return 0;
}
