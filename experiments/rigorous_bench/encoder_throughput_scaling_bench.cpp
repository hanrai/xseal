// SPDX-License-Identifier: Apache-2.0
// Copyright (c) 2026 hanrai. All Rights Reserved.

#include <iostream>
#include <vector>
#include <chrono>
#include <iomanip>
#include <numeric>
#include <algorithm>
#include <cstring>
#include <thread>
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

void minimap2_encode_worker(const char* start, const char* end, uint8_t* out) {
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
}

void xseal_encode_worker(const char* start, const char* end, uint8_t* out) {
    XSealEncoder encoder;
    XSealEncoderState state;
    encoder.init_state(&state, out, (end - start) / 4 + 1024);
    const char* ptr = start;
    encoder.encode_chunk(&state, &ptr, end);
    encoder.flush_state(&state);
}

std::pair<double, double> get_stats(const std::vector<double>& v) {
    double sum = std::accumulate(v.begin(), v.end(), 0.0);
    double mean = sum / v.size();
    double sq_sum = std::inner_product(v.begin(), v.end(), v.begin(), 0.0);
    double stdev = std::sqrt(std::abs(sq_sum / v.size() - mean * mean));
    return {mean, stdev};
}

int main(int argc, char** argv) {
    if (argc < 5) {
        std::cerr << "Usage: " << argv[0] << " <size_mb> <threads> <mode: minimap2|xseal> <newline_freq>\n";
        return 1;
    }
    size_t data_size = std::stoull(argv[1]) * 1024 * 1024;
    int num_threads = std::stoi(argv[2]);
    std::string mode = argv[3];
    int newline_freq = std::stoi(argv[4]);

    char* input = (char*)_mm_malloc(data_size, 64);
    uint8_t* output = (uint8_t*)_mm_malloc(data_size / 4 + 1024 * num_threads, 64);

    srand(42);
    for (size_t i = 0; i < data_size; ++i) {
        if (newline_freq > 0 && i > 0 && i % newline_freq == 0) {
            input[i] = '\n';
        } else {
            input[i] = "ACGT"[rand() % 4];
        }
    }

    std::vector<double> tps;
    for (int trial = 0; trial < 13; ++trial) {
        std::vector<std::thread> threads;
        size_t chunk_size = data_size / num_threads;
        
        auto start_t = std::chrono::high_resolution_clock::now();
        for (int i = 0; i < num_threads; ++i) {
            const char* s = input + i * chunk_size;
            const char* e = (i == num_threads - 1) ? (input + data_size) : (s + chunk_size);
            uint8_t* o = output + i * (chunk_size / 4 + 1024);
            if (mode == "minimap2") {
                threads.emplace_back(minimap2_encode_worker, s, e, o);
            } else {
                threads.emplace_back(xseal_encode_worker, s, e, o);
            }
        }
        for (auto& t : threads) t.join();
        auto end_t = std::chrono::high_resolution_clock::now();
        
        if (trial >= 3) {
            std::chrono::duration<double> diff = end_t - start_t;
            tps.push_back((data_size / (1024.0 * 1024.0 * 1024.0)) / diff.count());
        }
    }

    auto s = get_stats(tps);
    std::cout << "RESULT|" << mode << "|" << num_threads << "|" << newline_freq << "|" << s.first << "|" << s.second << "|GiB/s\n";

    _mm_free(input); _mm_free(output);
    return 0;
}
