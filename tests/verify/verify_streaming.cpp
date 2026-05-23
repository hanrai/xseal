// SPDX-License-Identifier: Apache-2.0
// Copyright (c) 2026 hanrai. All Rights Reserved.

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <algorithm>
#include <cstdint>
#include <immintrin.h>
#include <xseal/xseal_syncmer.hpp>
#include <xseal/xseal_minimizer.hpp>
#include <xseal/xseal_encoder.hpp>
#include "minimizer_scalar_ref.hpp"

inline uint32_t fast_hash(uint32_t x) {
    const uint32_t murmur_c = 0x85ebca6b;
    x ^= x >> 13; x *= murmur_c; x ^= x >> 16;
    return x;
}

inline uint32_t reverse_complement(uint32_t x, int s) {
    uint32_t res = 0;
    for (int i = 0; i < s; ++i) {
        uint32_t b = (x >> (2 * i)) & 3;
        res |= ((b ^ 3) << (2 * (s - 1 - i)));
    }
    return res;
}

struct Hit { uint32_t pos; uint32_t hash; };

void get_syncmer_baseline(const std::string& seq, int k, int s, int t, bool open, size_t abs_start, std::vector<Hit>& hits) {
    int w = k - s + 1;
    uint32_t mask = (s >= 16) ? 0xFFFFFFFF : (1U << (2 * s)) - 1;
    for (size_t chunk_start = 0; chunk_start + 128 + w - 1 <= seq.length(); chunk_start += 128) {
        std::vector<uint32_t> h_buf(128 + w - 1);
        for (int i = 0; i < 128 + w - 1; ++i) {
            uint32_t smer = 0;
            for (int j = 0; j < s; ++j) {
                int b = (seq[chunk_start + i + j] == 'A') ? 0 : (seq[chunk_start + i + j] == 'C') ? 1 : (seq[chunk_start + i + j] == 'G') ? 2 : 3;
                smer |= ((uint32_t)b << (2 * j));
            }
            smer &= mask;
            h_buf[i] = fast_hash(std::min(smer, reverse_complement(smer, s)));
        }
        for (int i = 0; i < 128; ++i) {
            uint32_t min_h = 0xFFFFFFFF;
            for (int j = 0; j < w; ++j) if (h_buf[i + j] < min_h) min_h = h_buf[i + j];
            bool match = (h_buf[i + t] == min_h);
            if (open) match |= (h_buf[i + w - 1] == min_h);
            if (match) hits.push_back({(uint32_t)(abs_start + chunk_start + i), h_buf[i + t]});
        }
    }
}

int main(int argc, char** argv) {
    if (argc < 5) { std::cerr << "Usage: " << argv[0] << " <fa> <mode:m/so/sc> <k> <w/s> [t]" << std::endl; return 1; }
    std::string mode = argv[2]; int k = std::stoi(argv[3]); int ws = std::stoi(argv[4]); int t = (argc > 5) ? std::stoi(argv[5]) : 0;
    std::ifstream fa(argv[1]); std::string line, segment;
    xseal::XSealEncoder encoder; xseal::XSealEncoderState enc_state;
    xseal::XSealMinimizer m_scanner(k, ws); xseal::XSealSyncmer s_scanner(k, ws);
    std::vector<uint8_t> enc_buf; std::vector<xseal::SyncmerHitFull> sim_hits; std::vector<Hit> base_hits;
    size_t abs_offset = 0; bool chr1_found = false;

    auto verify = [&](const std::string& s, size_t start) {
        int w = (mode == "m") ? ws : (k - ws + 1);
        if (s.length() < (size_t)(k + ws + 128)) return;
        base_hits.clear();
        if (mode != "m")
            get_syncmer_baseline(s, k, ws, t, (mode == "so"), start, base_hits);

        enc_buf.assign(s.length() / 4 + 256, 0);
        encoder.init_state(&enc_state, enc_buf.data(), enc_buf.size());
        const char* p = s.c_str(); const char* e = p + s.length();
        encoder.encode_chunk(&enc_state, &p, e);
        encoder.flush_state(&enc_state);
        size_t n_chunks = (enc_state.out_byte_idx + 31) / 32;
        if (n_chunks < 2) return;

        if (mode == "m") {
            MinimizerScalarRef ref(k);
            auto ref_hits = ref.scan_hits<21, true>(enc_buf.data(), n_chunks,
                                                    static_cast<uint32_t>(start));
            base_hits.reserve(ref_hits.size());
            for (const auto &h : ref_hits)
                base_hits.push_back({h.pos, h.hash});
        }

        sim_hits.assign(s.length() + 1024, {0,0});
        size_t count = 0;

        auto get_smer_hash_at = [&](const std::string& seq, size_t pos) -> uint32_t {
            uint32_t mask = (ws >= 16) ? 0xFFFFFFFF : (1U << (2 * ws)) - 1;
            uint32_t smer = 0;
            for (int j = 0; j < ws; ++j) {
                int b = (seq[pos + j] == 'A') ? 0 : (seq[pos + j] == 'C') ? 1 : (seq[pos + j] == 'G') ? 2 : 3;
                smer |= ((uint32_t)b << (2 * j));
            }
            smer &= mask;
            return fast_hash(std::min(smer, reverse_complement(smer, ws)));
        };

        if (mode == "m") {
            std::vector<uint32_t> sim_hit_pos(s.length() + 1024);
            std::vector<uint32_t> sim_hit_hashes(s.length() + 1024);
            count = m_scanner.scan<21, true>((const __m256i*)enc_buf.data(), n_chunks, sim_hit_pos.data(), sim_hit_hashes.data());
            for (size_t i = 0; i < count; ++i) {
                sim_hits[i].pos = sim_hit_pos[i] + start;
                sim_hits[i].hash = sim_hit_hashes[i] & 0xFFFFFF00;
            }
        } else if (mode == "so") {
            std::vector<uint32_t> sim_hit_pos(s.length() + 1024);
            count = s_scanner.scan<21, 0, true, true>((const __m256i*)enc_buf.data(), n_chunks, sim_hit_pos.data(), s.length());
            for (size_t i = 0; i < count; ++i) {
                sim_hits[i].pos = sim_hit_pos[i] + start;
                sim_hits[i].hash = get_smer_hash_at(s, sim_hit_pos[i] + 0);
            }
        } else {
            std::vector<uint32_t> sim_hit_pos(s.length() + 1024);
            count = s_scanner.scan<21, 5, false, true>((const __m256i*)enc_buf.data(), n_chunks, sim_hit_pos.data(), s.length());
            for (size_t i = 0; i < count; ++i) {
                sim_hits[i].pos = sim_hit_pos[i] + start;
                sim_hits[i].hash = get_smer_hash_at(s, sim_hit_pos[i] + 5);
            }
        }
        
        size_t s_ptr = 0, b_ptr = 0, w_end = start + s.length() - (k + w);
        while (s_ptr < count && b_ptr < base_hits.size()) {
            if (sim_hits[s_ptr].pos >= (uint32_t)w_end) { s_ptr++; continue; }
            if (base_hits[b_ptr].pos >= (uint32_t)w_end) { b_ptr++; continue; }
            if (sim_hits[s_ptr].pos != base_hits[b_ptr].pos || sim_hits[s_ptr].hash != base_hits[b_ptr].hash) {
                std::cerr << "Mismatch at " << base_hits[b_ptr].pos << ": SIMD(p=" << sim_hits[s_ptr].pos << ",h=" << sim_hits[s_ptr].hash << ") vs Base(p=" << base_hits[b_ptr].pos << ",h=" << base_hits[b_ptr].hash << ")" << std::endl;
                exit(1);
            }
            s_ptr++; b_ptr++;
        }
    };

    while (std::getline(fa, line)) {
        if (line.empty()) continue;
        if (line[0] == '>') {
            if (chr1_found) break;
            if (line.find("chr1") != std::string::npos) chr1_found = true;
            abs_offset = 0; segment.clear(); continue;
        }
        if (!chr1_found) continue;
        for (char c : line) {
            char uc = std::toupper(c);
            if (uc == 'A' || uc == 'C' || uc == 'G' || uc == 'T') segment += uc;
            else { verify(segment, abs_offset - segment.length()); segment.clear(); }
            abs_offset++;
        }
    }
    verify(segment, abs_offset - segment.length());
    if (chr1_found) std::cout << mode << " Verification PASSED!" << std::endl;
    return 0;
}
