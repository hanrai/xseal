// SPDX-License-Identifier: Apache-2.0
// Copyright (c) 2026 hanrai. All Rights Reserved.

#include <iostream>
#include <vector>
#include <chrono>
#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>
#include <immintrin.h>
#include <iomanip>
#include <thread>
#include <numeric>
#include <cmath>
#include <algorithm>
#include <cstring>

#include <xseal/xseal_syncmer.hpp>
#include <xseal/xseal_minimizer.hpp>

using namespace xseal;

struct BenchResult {
    double throughput;
    size_t hits;
    uint64_t digest;
};

// Position-only anti-DCE digest (matches xseal_hot_bench; no hash materialization).
inline uint64_t digest_positions(const uint32_t* pos, size_t hits) {
    uint64_t d = 0;
    for (size_t i = 0; i < hits; ++i) {
        d += static_cast<uint64_t>(pos[i]);
        d = (d << 17) | (d >> (64 - 17));
    }
    asm volatile("" : "+r"(d) : : "memory");
    return d;
}

uint32_t** global_out_pos = nullptr;
size_t per_thread_buf_size = 0;

template<int W>
BenchResult run_closed_sync_benchmark_impl(const uint8_t* data, size_t num_bases, int k, int s, int num_threads) {
    size_t bases_per_thread = (num_bases / num_threads / 128) * 128;
    std::vector<std::thread> threads;
    std::vector<size_t> thread_hits(num_threads, 0);
    XSealSyncmer scanner(k, s);
    
    // Clear buffer outside timing
    for (int i = 0; i < num_threads; ++i) std::memset(global_out_pos[i], 0, per_thread_buf_size * sizeof(uint32_t));

    auto t0 = std::chrono::high_resolution_clock::now();
    for (int t = 0; t < num_threads; ++t) {
        threads.emplace_back([&, t]() {
            size_t start_base = t * bases_per_thread;
            size_t thread_num_bases = (t == num_threads - 1) ? (num_bases - start_base) : bases_per_thread;
            if (thread_num_bases < 256) return;
            size_t num_chunks = thread_num_bases / 128;
            const __m256i* chunk_ptr = (const __m256i*)(data + (start_base / 4));
            uint32_t* out_pos = global_out_pos[t];
            thread_hits[t] = scanner.scan<W, 0, true, true>(chunk_ptr, num_chunks, out_pos, thread_num_bases);
            asm volatile("" : : "g"(thread_hits[t]), "g"(out_pos) : "memory");
        });
    }
    for (auto& t : threads) t.join();
    auto t1 = std::chrono::high_resolution_clock::now();
    double secs = std::chrono::duration<double>(t1 - t0).count();
    uint64_t digest = 0;
    for (int t = 0; t < num_threads; ++t) {
        digest += digest_positions(global_out_pos[t], thread_hits[t]);
    }
    return { (double)num_bases / secs / 1e9, std::accumulate(thread_hits.begin(), thread_hits.end(), 0ULL), digest };
}

template<int W, int T>
BenchResult run_open_sync_benchmark_impl(const uint8_t* data, size_t num_bases, int k, int s, int num_threads) {
    size_t bases_per_thread = (num_bases / num_threads / 128) * 128;
    std::vector<std::thread> threads;
    std::vector<size_t> thread_hits(num_threads, 0);
    XSealSyncmer scanner(k, s);
    for (int i = 0; i < num_threads; ++i) std::memset(global_out_pos[i], 0, per_thread_buf_size * sizeof(uint32_t));
    auto t0 = std::chrono::high_resolution_clock::now();
    for (int t = 0; t < num_threads; ++t) {
        threads.emplace_back([&, t]() {
            size_t start_base = t * bases_per_thread;
            size_t thread_num_bases = (t == num_threads - 1) ? (num_bases - start_base) : bases_per_thread;
            if (thread_num_bases < 256) return;
            size_t num_chunks = thread_num_bases / 128;
            const __m256i* chunk_ptr = (const __m256i*)(data + (start_base / 4));
            thread_hits[t] = scanner.scan<W, T, false, true>(chunk_ptr, num_chunks, global_out_pos[t], thread_num_bases);
        });
    }
    for (auto& t : threads) t.join();
    auto t1 = std::chrono::high_resolution_clock::now();
    double secs = std::chrono::duration<double>(t1 - t0).count();
    uint64_t digest = 0;
    for (int t = 0; t < num_threads; ++t) {
        digest += digest_positions(global_out_pos[t], thread_hits[t]);
    }
    return { (double)num_bases / secs / 1e9, std::accumulate(thread_hits.begin(), thread_hits.end(), 0ULL), digest };
}

template<int W>
BenchResult run_minimizer_benchmark_impl(const uint8_t* data, size_t num_bases, int k, int w, int num_threads) {
    size_t bases_per_thread = (num_bases / num_threads / 128) * 128;
    std::vector<std::thread> threads;
    std::vector<size_t> thread_hits(num_threads, 0);
    XSealMinimizer scanner(k, w);
    for (int i = 0; i < num_threads; ++i) std::memset(global_out_pos[i], 0, per_thread_buf_size * sizeof(uint32_t));
    auto t0 = std::chrono::high_resolution_clock::now();
    for (int t = 0; t < num_threads; ++t) {
        threads.emplace_back([&, t]() {
            size_t start_base = t * bases_per_thread;
            size_t thread_num_bases = (t == num_threads - 1) ? (num_bases - start_base) : bases_per_thread;
            if (thread_num_bases < 256) return;
            size_t num_chunks = thread_num_bases / 128;
            const __m256i* chunk_ptr = (const __m256i*)(data + (start_base / 4));
            // out_hash omitted: position-only output (no hash materialization).
            thread_hits[t] = scanner.scan<W, true>(chunk_ptr, num_chunks, global_out_pos[t]);
            asm volatile("" : : "g"(thread_hits[t]), "g"(global_out_pos[t]) : "memory");
        });
    }
    for (auto& t : threads) t.join();
    auto t1 = std::chrono::high_resolution_clock::now();
    double secs = std::chrono::duration<double>(t1 - t0).count();
    uint64_t digest = 0;
    for (int t = 0; t < num_threads; ++t) {
        digest += digest_positions(global_out_pos[t], thread_hits[t]);
    }
    return { (double)num_bases / secs / 1e9, std::accumulate(thread_hits.begin(), thread_hits.end(), 0ULL), digest };
}

BenchResult dispatch_sync(int w, const uint8_t* data, size_t num_bases, int k, int s, int num_threads) {
    switch(w) {
        case 5: return run_closed_sync_benchmark_impl<5>(data, num_bases, k, s, num_threads);
        case 9: return run_closed_sync_benchmark_impl<9>(data, num_bases, k, s, num_threads);
        case 11: return run_closed_sync_benchmark_impl<11>(data, num_bases, k, s, num_threads);
        case 13: return run_closed_sync_benchmark_impl<13>(data, num_bases, k, s, num_threads);
        case 17: return run_closed_sync_benchmark_impl<17>(data, num_bases, k, s, num_threads);
        case 21: return run_closed_sync_benchmark_impl<21>(data, num_bases, k, s, num_threads);
        default: return run_closed_sync_benchmark_impl<21>(data, num_bases, k, s, num_threads);
    }
}

BenchResult dispatch_open(int w, const uint8_t* data, size_t num_bases, int k, int s, int num_threads) {
    if (w == 5) return run_open_sync_benchmark_impl<5, 2>(data, num_bases, k, s, num_threads);
    if (w == 9) return run_open_sync_benchmark_impl<9, 4>(data, num_bases, k, s, num_threads);
    if (w == 11) return run_open_sync_benchmark_impl<11, 5>(data, num_bases, k, s, num_threads);
    if (w == 13) return run_open_sync_benchmark_impl<13, 6>(data, num_bases, k, s, num_threads);
    if (w == 17) return run_open_sync_benchmark_impl<17, 8>(data, num_bases, k, s, num_threads);
    return run_open_sync_benchmark_impl<21, 10>(data, num_bases, k, s, num_threads);
}

BenchResult dispatch_min(int w, const uint8_t* data, size_t num_bases, int k, int num_threads) {
    switch(w) {
        case 11: return run_minimizer_benchmark_impl<11>(data, num_bases, k, w, num_threads);
        case 21: return run_minimizer_benchmark_impl<21>(data, num_bases, k, w, num_threads);
        case 31: return run_minimizer_benchmark_impl<31>(data, num_bases, k, w, num_threads);
        case 51: return run_minimizer_benchmark_impl<51>(data, num_bases, k, w, num_threads);
        default: return run_minimizer_benchmark_impl<21>(data, num_bases, k, w, num_threads);
    }
}

int main(int argc, char** argv) {
    if (argc < 6) return 1;
    int fd = open(argv[1], O_RDONLY);
    struct stat st; fstat(fd, &st);
    uint8_t* mmapped_ptr = (uint8_t*)mmap(NULL, st.st_size, PROT_READ, MAP_PRIVATE, fd, 0);
    size_t padded_size = (st.st_size + 1024 * 1024 + 63) & ~63;
    void* aligned_buf = nullptr;
    posix_memalign(&aligned_buf, 64, padded_size);
    std::memset(aligned_buf, 0, padded_size);
    std::memcpy(aligned_buf, mmapped_ptr, st.st_size);
    munmap(mmapped_ptr, st.st_size); close(fd);
    
    uint8_t* data = (uint8_t*)aligned_buf;
    size_t num_bases = (st.st_size * 4);
    int num_threads = std::stoi(argv[2]);
    int K = std::stoi(argv[3]), S = std::stoi(argv[4]), W = std::stoi(argv[5]);
    
    // Safer buffer size: 3.2GB per thread
    per_thread_buf_size = std::max((size_t)(1024 * 1024 * 64), (size_t)(num_bases / num_threads / 4 + 1024));
    global_out_pos = new uint32_t*[num_threads];
    for (int i = 0; i < num_threads; ++i) {
        global_out_pos[i] = (uint32_t*)aligned_alloc(64, (per_thread_buf_size + 256) * sizeof(uint32_t));
    }
    
    int W_sync = K - S + 1;
    for (int i = 0; i < 3; ++i) dispatch_sync(W_sync, data, num_bases, K, S, num_threads);

    auto stats = [](const std::vector<double>& v) {
        double sum = std::accumulate(v.begin(), v.end(), 0.0);
        double mean = sum / v.size();
        double sq_sum = std::inner_product(v.begin(), v.end(), v.begin(), 0.0);
        double stdev = std::sqrt(std::abs(sq_sum / v.size() - mean * mean));
        return std::make_pair(mean, stdev);
    };

    std::vector<double> sync_t, open_t, min_t;
    size_t s_hits, o_hits, m_hits; uint64_t s_dg, o_dg, m_dg;
    for (int i = 0; i < 10; ++i) {
        auto r1 = dispatch_sync(W_sync, data, num_bases, K, S, num_threads);
        sync_t.push_back(r1.throughput); s_hits = r1.hits; s_dg = r1.digest;
        auto r2 = dispatch_open(W_sync, data, num_bases, K, S, num_threads);
        open_t.push_back(r2.throughput); o_hits = r2.hits; o_dg = r2.digest;
        auto r3 = dispatch_min(W, data, num_bases, K, num_threads);
        min_t.push_back(r3.throughput); m_hits = r3.hits; m_dg = r3.digest;
    }

    auto s1 = stats(sync_t); auto s2 = stats(open_t); auto s3 = stats(min_t);
    std::cout << std::fixed << std::setprecision(4);
    std::cout << "RESULT|Xseal|Syncmer|" << num_threads << "|" << s1.first << "|" << s1.second << "|" << s_hits << "|" << std::hex << s_dg << std::dec << "\n";
    std::cout << "RESULT|Xseal|OpenSync|" << num_threads << "|" << s2.first << "|" << s2.second << "|" << o_hits << "|" << std::hex << o_dg << std::dec << "\n";
    std::cout << "RESULT|Xseal|Minimizer|" << num_threads << "|" << s3.first << "|" << s3.second << "|" << m_hits << "|" << std::hex << m_dg << std::dec << "\n";
    if (argc >= 7) {
        FILE* fp = fopen(argv[6], "wb");
        if (fp) {
            // We use m_hits (Minimizer) for mutation robustness comparison
            fwrite(global_out_pos[0], sizeof(uint32_t), m_hits, fp);
            fclose(fp);
            std::cout << "INFO|Dumped " << m_hits << " positions to " << argv[6] << "\n";
        }
    }

    for (int i = 0; i < num_threads; ++i) free(global_out_pos[i]);
    delete[] global_out_pos; free(aligned_buf);
    return 0;
}
