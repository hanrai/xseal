// SPDX-License-Identifier: Apache-2.0
// Copyright (c) 2026 hanrai. All Rights Reserved.

#include <xseal/xseal_syncmer.hpp>
#include <vector>
#include <chrono>
#include <iostream>
#include <random>
#include <thread>
#include <iomanip>

using namespace xseal;

void bench_thread(int thread_id, size_t num_chunks, int iterations) {
    XSealSyncmer scanner(31, 11);
    // Pre-allocate and fill random 2-bit data (simulate encoded sequence)
    std::vector<__m256i> data(num_chunks + 16);
    std::mt19937 gen(42 + thread_id);
    uint32_t* ptr = (uint32_t*)data.data();
    for (size_t i = 0; i < num_chunks * 8; ++i) ptr[i] = gen();

    // Pre-allocate result space
    std::vector<uint32_t> results(num_chunks * 32);
    
    auto start = std::chrono::high_resolution_clock::now();
    size_t total_hits = 0;
    for (int it = 0; it < iterations; ++it) {
        // W=21, T=0, CLOSED=true, CANONICAL=true, FULL_OUTPUT=false
        size_t hits = scanner.scan<21, 0, true, true, false>(data.data(), num_chunks, 0, results.data());
        total_hits += hits;
    }
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> diff = end - start;
    
    double total_gb = (double)num_chunks * 128 * iterations / 1e9;
    std::cout << "  🧵 Thread " << std::setw(2) << thread_id << ": " 
              << std::fixed << std::setprecision(2) << total_gb / diff.count() << " GB/s"
              << " (Density: " << (double)total_hits / ((double)num_chunks * 128 * iterations) * 100.0 << "%)\n";
}

int main(int argc, char** argv) {
    int num_threads = (argc > 1) ? std::stoi(argv[1]) : 1;
    size_t chunks_per_thread = 100000; // 12.8 MB (Fits in L3)
    int iterations = 100;

    std::cout << "====================================================\n";
    std::cout << " 🚀 XSeal PURE SYNC-SCAN MICRO-BENCHMARK\n";
    std::cout << "----------------------------------------------------\n";
    std::cout << " 🛠️  Config: W=21, K=31, S=11, Canonical=True\n";
    std::cout << " 🧵 Threads: " << num_threads << "\n";
    std::cout << "====================================================\n";

    std::vector<std::thread> threads;
    for (int i = 0; i < num_threads; ++i) {
        threads.emplace_back(bench_thread, i, chunks_per_thread, iterations);
    }
    for (auto& t : threads) t.join();
    
    std::cout << "====================================================\n";

    return 0;
}
