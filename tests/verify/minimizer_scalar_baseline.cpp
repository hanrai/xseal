// SPDX-License-Identifier: Apache-2.0
// Copyright (c) 2026 hanrai. All Rights Reserved.

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <algorithm>
#include <cstdint>

inline uint32_t fast_hash(uint32_t x) {
    const uint32_t murmur_c = 0x85ebca6b;
    x ^= x >> 13;
    x *= murmur_c;
    x ^= x >> 16;
    return x;
}

inline uint32_t reverse_complement(uint32_t x, int k) {
    uint32_t res = 0;
    for (int i = 0; i < k; ++i) {
        uint32_t b = (x >> (2 * i)) & 3;
        res |= ((b ^ 3) << (2 * (k - 1 - i)));
    }
    return res;
}

void process_sequence(const std::string& name, const std::string& seq, int k, int w) {
    std::cerr << "Processing " << name << " (" << seq.length() << " bp)..." << std::endl;
    uint32_t mask = (k >= 16) ? 0xFFFFFFFF : (1U << (2 * k)) - 1;

    auto get_h = [&](size_t pos) -> uint32_t {
        uint32_t kmer = 0;
        for (int i = 0; i < k; ++i) {
            char c = std::toupper(seq[pos + i]);
            int b = (c == 'A') ? 0 : (c == 'C') ? 1 : (c == 'G') ? 2 : (c == 'T') ? 3 : -1;
            if (b == -1) return 0xFFFFFFFF;
            kmer |= ((uint32_t)b << (2 * i));
        }
        kmer &= mask;
        uint32_t canon = std::min(kmer, reverse_complement(kmer, k));
        return fast_hash(canon);
    };

    size_t global_abs_offset = 0;
    const char* ptr = seq.c_str();
    const char* end = ptr + seq.length();

    while (ptr < end) {
        while (ptr < end && (*ptr == 'N' || *ptr == 'n')) { ptr++; global_abs_offset++; }
        if (ptr == end) break;
        const char* seg_start = ptr;
        while (ptr < end && !(*ptr == 'N' || *ptr == 'n')) ptr++;
        size_t seg_len = ptr - seg_start;
        
        if (seg_len >= (size_t)(k + w - 1)) {
            uint32_t last_abs_pos = 0xFFFFFFFF;
            // Block-wise processing (128-bp chunks)
            for (size_t chunk_start = 0; chunk_start + 128 + w - 1 <= seg_len; chunk_start += 128) {
                std::vector<uint32_t> hashes(128 + w - 1);
                for (int i = 0; i < 128 + w - 1; ++i) hashes[i] = get_h(seg_start - seq.c_str() + chunk_start + i);

                for (int batch = 0; batch < 128; batch += 8) {
                    for (int i = 0; i < 8; ++i) {
                        int win_idx = batch + i;
                        uint32_t min_h = 0xFFFFFFFF;
                        int min_rel = 0;
                        for (int j = 0; j < w; ++j) {
                            uint32_t h = hashes[win_idx + j];
                            if (h < min_h) { min_h = h; min_rel = win_idx + j; }
                        }
                        uint32_t abs_pos = global_abs_offset + chunk_start + min_rel;
                        if (min_h != 0xFFFFFFFF && abs_pos != last_abs_pos) {
                            std::cout << abs_pos << "\t" << (min_h & 0xFFFFFF00) << "\n";
                            last_abs_pos = abs_pos;
                        }
                    }
                }
            }
        }
        global_abs_offset += seg_len;
    }
}

int main(int argc, char** argv) {
    if (argc < 4) return 1;
    std::string fa_path = argv[1];
    int k = std::stoi(argv[2]); int w = std::stoi(argv[3]);
    std::ifstream fa(fa_path); std::string line, seq;
    while (std::getline(fa, line)) {
        if (line.empty() || line[0] == '>') continue;
        seq += line;
    }
    process_sequence("chr1", seq, k, w);
    return 0;
}
