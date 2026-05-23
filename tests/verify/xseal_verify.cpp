// SPDX-License-Identifier: Apache-2.0
// Copyright (c) 2026 hanrai. All Rights Reserved.

#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <fstream>
#include <cstdint>

/**
 * XSeal Official Verification Tool (Naive Reference Implementation)
 * 
 * This tool implements the syncmer sampling logic in a plain, scalar, and
 * obviously-correct manner. It is used as the "Ground Truth" to verify the
 * bit-perfect accuracy of the high-performance SIMD-optimized XSeal core.
 * 
 * Usage: ./xseal_verify <input.fa> [K] [S] [T]
 */

uint32_t fast_hash(uint32_t x) {
    x ^= x >> 13;
    x *= 0x85ebca6b;
    x ^= x >> 16;
    return x;
}

uint32_t get_smer_bits(const std::string& smer) {
    uint32_t bits = 0;
    for (int i = 0; i < (int)smer.length(); ++i) {
        uint32_t val = 0;
        if (smer[i] == 'C' || smer[i] == 'c') val = 1;
        else if (smer[i] == 'G' || smer[i] == 'g') val = 2;
        else if (smer[i] == 'T' || smer[i] == 't') val = 3;
        bits |= (val << (2 * i));
    }
    return bits;
}

uint32_t reverse_complement_bits(uint32_t bits, int s) {
    uint32_t rc = 0;
    for (int i = 0; i < s; ++i) {
        uint32_t val = (bits >> (2 * i)) & 3;
        uint32_t rc_val = val ^ 3;
        rc |= (rc_val << (2 * (s - 1 - i)));
    }
    return rc;
}

void process_segment(const std::string& segment, int k, int s, int t, uint64_t start_pos) {
    if (segment.length() < (size_t)k) return;
    int w = k - s + 1;
    std::vector<uint32_t> hashes;
    hashes.reserve(segment.length() - s + 1);
    
    for (int i = 0; i <= (int)segment.length() - s; ++i) {
        uint32_t fwd = get_smer_bits(segment.substr(i, s));
        uint32_t rev = reverse_complement_bits(fwd, s);
        hashes.push_back(fast_hash(std::min(fwd, rev)));
    }

    for (int i = 0; i <= (int)hashes.size() - w; ++i) {
        uint32_t min_h = hashes[i];
        for (int j = 1; j < w; ++j) {
            if (hashes[i+j] < min_h) min_h = hashes[i+j];
        }
        if (hashes[i + t] == min_h || hashes[i + w - 1] == min_h) {
            std::cout << (start_pos + i) << "\n";
        }
    }
}

int main(int argc, char** argv) {
    if (argc < 2) {
        std::cerr << "XSeal Naive Verifier\nUsage: " << argv[0] << " <input.fa> [K] [S] [T]\n";
        return 1;
    }
    
    std::string path = argv[1];
    int k = (argc > 2) ? std::stoi(argv[2]) : 31;
    int s = (argc > 3) ? std::stoi(argv[3]) : 11;
    int t = (argc > 4) ? std::stoi(argv[4]) : 0;

    std::ifstream f(path);
    if (!f.is_open()) {
        std::cerr << "Error: Cannot open file " << path << "\n";
        return 1;
    }

    std::string line, seq;
    auto flush_seq = [&]() {
        if (seq.empty()) return;
        uint64_t pos = 0;
        while (pos < seq.length()) {
            size_t start = seq.find_first_of("ACGTacgt", pos);
            if (start == std::string::npos) break;
            size_t end = seq.find_first_not_of("ACGTacgt", start);
            if (end == std::string::npos) end = seq.length();
            process_segment(seq.substr(start, end - start), k, s, t, start);
            pos = end;
        }
        seq.clear();
    };

    while (std::getline(f, line)) {
        if (line.empty()) continue;
        if (line[0] == '>') {
            flush_seq();
        } else {
            seq += line;
        }
    }
    flush_seq();
    
    return 0;
}
