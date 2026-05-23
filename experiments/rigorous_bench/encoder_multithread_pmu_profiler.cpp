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

void minimap2_encode(const char* start, const char* end, uint8_t* out) {
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

int main(int argc, char** argv) {
    if (argc < 4) {
        std::cerr << "Usage: " << argv[0] << " <size_mb> <threads> <mode: minimap2|xseal>\n";
        return 1;
    }
    size_t data_size = std::stoull(argv[1]) * 1024 * 1024;
    int num_threads = std::stoi(argv[2]);
    std::string mode = argv[3];

    char* input = (char*)_mm_malloc(data_size, 64);
    uint8_t* output = (uint8_t*)_mm_malloc(data_size / 4 + 1024 * num_threads, 64);

    const char bases[] = "ACGT";
    for (size_t i = 0; i < data_size; ++i) input[i] = bases[rand() % 4];

    size_t chunk_size = data_size / num_threads;
    std::vector<std::thread> threads;

    auto start_t = std::chrono::high_resolution_clock::now();
    
    for (int i = 0; i < num_threads; ++i) {
        const char* s = input + i * chunk_size;
        const char* e = (i == num_threads - 1) ? (input + data_size) : (s + chunk_size);
        uint8_t* o = output + i * (chunk_size / 4 + 1024);
        if (mode == "minimap2") {
            threads.emplace_back(minimap2_encode, s, e, o);
        } else {
            threads.emplace_back(xseal_encode_worker, s, e, o);
        }
    }

    for (auto& t : threads) t.join();

    auto end_t = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> diff = end_t - start_t;
    
    std::cout << "Mode: " << mode << " | Threads: " << num_threads 
              << " | Throughput: " << (data_size / 1e9) / diff.count() << " GB/s\n";

    _mm_free(input); _mm_free(output);
    return 0;
}
