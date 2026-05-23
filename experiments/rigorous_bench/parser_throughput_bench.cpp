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
#include <zlib.h>
#include <cstdio>
#include <sys/stat.h>
#include "xseal/xseal_fasta_parser.hpp"
#include "xseal/xseal_encoder.hpp"

// Modified parser to count transitions if needed, but we can just check run_table.size()
typedef struct {
    const char *buf;
    size_t size;
    size_t pos;
} mem_buf_t;

static int mem_read(mem_buf_t *mb, uint8_t *buf, int len) {
    if (mb->pos >= mb->size) return 0;
    size_t rem = mb->size - mb->pos;
    size_t n = (size_t)len < rem ? (size_t)len : rem;
    std::memcpy(buf, mb->buf + mb->pos, n);
    mb->pos += n;
    return (int)n;
}

#include "kseq.h"
KSEQ_INIT(mem_buf_t*, mem_read)

using namespace xseal;

std::pair<double, double> stats(const std::vector<double>& v) {
    if (v.empty()) return {0.0, 0.0};
    double sum = std::accumulate(v.begin(), v.end(), 0.0);
    double mean = sum / v.size();
    double sq_sum = std::inner_product(v.begin(), v.end(), v.begin(), 0.0);
    double stdev = std::sqrt(std::abs(sq_sum / v.size() - mean * mean));
    return {mean, stdev};
}

int main(int argc, char** argv) {
    const char* filename = nullptr;
    bool is_random = false;
    if (argc > 1) filename = argv[1];

    char* input = nullptr;
    size_t actual_size = 0;

    if (filename && (std::string(filename).find(".fa") != std::string::npos || 
                     std::string(filename).find(".fasta") != std::string::npos)) {
        struct stat st;
        if (stat(filename, &st) != 0) {
            std::cerr << "File not found: " << filename << "\n";
            return 1;
        }
        actual_size = st.st_size;
        input = (char*)_mm_malloc(actual_size, 64);
        FILE* f = fopen(filename, "rb");
        fread(input, 1, actual_size, f);
        fclose(f);
        std::cout << "DATA_TYPE|Real|" << filename << "|" << actual_size << "\n";
    } else {
        is_random = true;
        size_t data_size = 512 * 1024 * 1024; // 512MB
        if (argc > 1) data_size = std::stoull(argv[1]) * 1024 * 1024;
        input = (char*)_mm_malloc(data_size, 64);
        srand(42);
        size_t p = 0;
        while (p < data_size - 1000) {
            input[p++] = '>';
            for (int i = 0; i < 50; ++i) input[p++] = 'A' + (rand() % 26);
            input[p++] = '\n';
            for (int i = 0; i < 800; ++i) {
                input[p++] = "ACGT"[rand() % 4];
                if (i % 80 == 79) input[p++] = '\n';
            }
            if (p > 0 && input[p-1] != '\n') input[p++] = '\n';
        }
        actual_size = p;
        std::cout << "DATA_TYPE|Random|Mock|" << actual_size << "\n";
    }

    uint8_t* output = (uint8_t*)_mm_malloc(actual_size / 4 + 1024, 64);
    XsealRunTable run_table;
    run_table.reserve(20000000);

    auto run_bench = [&](const std::string& name, auto func) {
        std::vector<double> tps;
        std::vector<double> cpb;
        size_t seg_count = 0;
        const int WARMUP = 3;
        const int TRIALS = 10;
        for (int i = 0; i < WARMUP + TRIALS; ++i) {
            run_table.clear();
            uint64_t t0 = __rdtsc();
            auto start = std::chrono::high_resolution_clock::now();
            func();
            auto end = std::chrono::high_resolution_clock::now();
            uint64_t t1 = __rdtsc();
            seg_count = run_table.size();
            if (i >= WARMUP) {
                std::chrono::duration<double> diff = end - start;
                tps.push_back((double)actual_size / (1024.0*1024.0*1024.0) / diff.count());
                cpb.push_back((double)(t1 - t0) / actual_size);
            }
        }
        auto s_tp = stats(tps);
        auto s_cp = stats(cpb);
        std::cout << "METRIC|" << name << "|" << s_tp.first << "|" << s_tp.second << "|" << s_cp.first << "|" << seg_count << "\n";
    };

    // 1. kseq.h
    run_bench("kseq.h", [&]() {
        mem_buf_t mb = {input, actual_size, 0};
        kseq_t *ks = kseq_init(&mb);
        while (kseq_read(ks) >= 0) { }
        kseq_destroy(ks);
    });

    // 2. Xseal Parser Only
    XsealFastaParser parser;
    run_bench("XSeal_Parser", [&]() {
        parser.parse((const uint8_t*)input, actual_size, 0, run_table);
        parser.finish(run_table);
    });

    // 3. Xseal Fused
    XSealEncoder encoder;
    run_bench("XSeal_Fused", [&]() {
        parser.parse((const uint8_t*)input, actual_size, 0, run_table);
        parser.finish(run_table);
        XSealEncoderState enc_state;
        encoder.init_state(&enc_state, output, actual_size / 4 + 1024);
        for (const auto& seg : run_table) {
            if (seg.type == SegType::BASES) {
                const char* ptr = input + seg.file_offset;
                encoder.encode_chunk(&enc_state, &ptr, ptr + seg.byte_len);
            }
        }
        encoder.flush_state(&enc_state);
    });

    _mm_free(input); _mm_free(output);
    return 0;
}
