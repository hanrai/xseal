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
#include "pure_encoder.hpp"

using namespace xseal;

// Minimap2's standard nt4 table
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

std::pair<double, double> stats(const std::vector<double>& v) {
    double sum = std::accumulate(v.begin(), v.end(), 0.0);
    double mean = sum / v.size();
    double sq_sum = std::inner_product(v.begin(), v.end(), v.begin(), 0.0);
    double stdev = std::sqrt(std::abs(sq_sum / v.size() - mean * mean));
    return {mean, stdev};
}

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
    size_t data_size = 1024 * 1024 * 1024; // 1GB
    if (argc > 1) data_size = std::stoull(argv[1]) * 1024 * 1024;

    char* input = (char*)_mm_malloc(data_size, 64);
    uint8_t* output = (uint8_t*)_mm_malloc(data_size / 4 + 1024, 64);

    const char bases[] = "ACGT";
    for (size_t i = 0; i < data_size; ++i) input[i] = bases[rand() % 4];

    auto run_bench = [&](const std::string& name, auto func) {
        std::vector<double> throughputs;
        std::vector<double> cycles_per_base;
        
        for (int i = 0; i < 13; ++i) { // 3 warmup + 10 trials
            std::memset(output, 0, data_size / 4 + 1024);
            
            uint64_t t0 = __rdtsc();
            auto start = std::chrono::high_resolution_clock::now();
            func();
            auto end = std::chrono::high_resolution_clock::now();
            uint64_t t1 = __rdtsc();
            
            if (i >= 3) {
                std::chrono::duration<double> diff = end - start;
                throughputs.push_back((double)data_size / (1024.0 * 1024.0 * 1024.0) / diff.count());
                cycles_per_base.push_back((double)(t1 - t0) / data_size);
            }
        }
        
        auto s_tp = stats(throughputs);
        auto s_cy = stats(cycles_per_base);
        std::cout << std::fixed << std::setprecision(4);
        std::cout << "RESULT|" << name << "|" << s_tp.first << "|" << s_tp.second << "|" << s_cy.first << " GiB/s, " << s_cy.first << " cycles/base\n";
    };

    std::cout << "Benchmarking Front-end Preprocessing (Data size: " << data_size / (1024*1024) << " MB)...\n";

    run_bench("Minimap2 (nt4)", [&]() { minimap2_encode(input, input + data_size, output); });

    XSealEncoder prod_enc;
    run_bench("XSeal (Production)", [&]() {
        XSealEncoderState state;
        prod_enc.init_state(&state, output, data_size / 4 + 1024);
        const char* ptr = input;
        prod_enc.encode_chunk(&state, &ptr, input + data_size);
        prod_enc.flush_state(&state);
    });

    run_bench("XSeal (Pure SIMD)", [&]() { PureEncoder::encode_chunk(input, input + data_size, output); });

    _mm_free(input); _mm_free(output);
    return 0;
}
