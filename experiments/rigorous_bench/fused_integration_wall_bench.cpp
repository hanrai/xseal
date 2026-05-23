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
#include "xseal/xseal_syncmer.hpp"
#include "xseal/xseal_fasta_parser.hpp"

using namespace xseal;

// Mock work for Scanner to avoid empty loops
static uint64_t dummy_sum = 0;

void mock_scan(const uint8_t* twobit, size_t num_bases, XSealSyncmer& scanner) {
    size_t num_chunks = (num_bases / 128);
    if (num_chunks == 0) return;
    std::vector<uint32_t> hit_pos(num_bases / 5 + 1024); // More space
    size_t hits = scanner.scan<21, 10, false, true>(
        reinterpret_cast<const __m256i*>(twobit),
        num_chunks, hit_pos.data(), num_bases);
    dummy_sum += hits;
}

int main(int argc, char** argv) {
    size_t data_size = 512 * 1024 * 1024; // 512MB
    if (argc > 1) data_size = std::stoull(argv[1]) * 1024 * 1024;

    char* input = (char*)_mm_malloc(data_size, 64);
    uint8_t* twobit_buffer = (uint8_t*)_mm_malloc(data_size / 4 + 1024, 64);

    const char bases[] = "ACGT";
    srand(42);
    for (size_t i = 0; i < data_size; ++i) {
        if (i % 81 == 80) input[i] = '\n';
        else input[i] = bases[rand() % 4];
    }

    XsealFastaParser parser;
    XSealEncoder encoder;
    XSealSyncmer scanner(31, 11);

    auto run_bench = [&](const std::string& name, auto func) {
        std::cout << "Starting " << name << "..." << std::endl;
        std::vector<double> times;
        for (int i = 0; i < 5; ++i) {
            dummy_sum = 0;
            auto start = std::chrono::high_resolution_clock::now();
            func();
            auto end = std::chrono::high_resolution_clock::now();
            if (i >= 2) times.push_back(std::chrono::duration<double>(end - start).count());
            std::cout << "  Trial " << i << " done." << std::endl;
        }
        if (times.empty()) return;
        double avg_time = std::accumulate(times.begin(), times.end(), 0.0) / times.size();
        double throughput = (double)data_size / (1024.0*1024.0*1024.0) / avg_time;
        std::cout << "RESULT|" << name << "|" << throughput << " GiB/s | " << avg_time << " s\n";
    };

    // 1. Traditional Modular (Full Materialization)
    run_bench("Traditional Modular", [&]() {
        XsealRunTable run_table;
        parser.parse((const uint8_t*)input, data_size, 0, run_table);
        parser.finish(run_table);

        XSealEncoderState enc_state;
        encoder.init_state(&enc_state, twobit_buffer, data_size / 4 + 1024);
        for (auto& seg : run_table) {
            if (seg.type == SegType::BASES) {
                const char* p = input + seg.file_offset;
                encoder.encode_chunk(&enc_state, &p, p + seg.byte_len);
            }
        }
        encoder.flush_state(&enc_state);

        mock_scan(twobit_buffer, enc_state.total_bases_encoded, scanner);
    });

    _mm_free(input); _mm_free(twobit_buffer);
    return 0;
}
