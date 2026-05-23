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

int main(int argc, char** argv) {
    if (argc < 3) {
        std::cerr << "Usage: " << argv[0] << " <size_mb> <mode: kseq|xseal>\n";
        return 1;
    }
    char* input = nullptr;
    size_t actual_size = 0;
    std::string path_or_size = argv[1];
    std::string mode = argv[2];

    if (path_or_size.find(".fa") != std::string::npos || path_or_size.find(".fasta") != std::string::npos) {
        struct stat st;
        if (stat(path_or_size.c_str(), &st) != 0) return 1;
        actual_size = st.st_size;
        input = (char*)_mm_malloc(actual_size, 64);
        FILE* f = fopen(path_or_size.c_str(), "rb");
        fread(input, 1, actual_size, f);
        fclose(f);
    } else {
        size_t data_size = std::stoull(path_or_size) * 1024 * 1024;
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
    }

    XsealRunTable run_table;
    run_table.reserve(20000000);

    auto start_t = std::chrono::high_resolution_clock::now();
    if (mode == "kseq") {
        for (int i = 0; i < 10; ++i) {
            mem_buf_t mb = {input, actual_size, 0};
            kseq_t *ks = kseq_init(&mb);
            while (kseq_read(ks) >= 0) { }
            kseq_destroy(ks);
        }
    } else if (mode == "xseal") {
        XsealFastaParser parser;
        for (int i = 0; i < 10; ++i) {
            run_table.clear();
            parser.parse((const uint8_t*)input, actual_size, 0, run_table);
            parser.finish(run_table);
        }
    }
    auto end_t = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> diff = end_t - start_t;
    std::cout << "Mode: " << mode << " | Throughput: " << (actual_size * 10 / 1e9) / diff.count() << " GB/s\n";

    _mm_free(input);
    return 0;
}
