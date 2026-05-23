// SPDX-License-Identifier: Apache-2.0
// Copyright (c) 2026 hanrai. All Rights Reserved.

// In-memory multi-thread syncmer micro-benchmark (paper Part I style, not FASTA e2e).
#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <random>
#include <thread>
#include <vector>
#include <xseal/xseal_syncmer.hpp>

using namespace xseal;

static void bench_thread(int thread_id, size_t num_chunks, int iterations) {
  XSealSyncmer scanner(31, 11);
  std::vector<__m256i> data(num_chunks + 16);
  std::mt19937 gen(42 + thread_id);
  uint32_t *ptr = reinterpret_cast<uint32_t *>(data.data());
  for (size_t i = 0; i < num_chunks * 8; ++i)
    ptr[i] = gen();

  std::vector<uint32_t> results(num_chunks * 32 + 1024);

  size_t total_hits = 0;
  for (int it = 0; it < iterations; ++it) {
    const size_t bases = num_chunks * 128;
    size_t hits = scanner.scan<21, 0, true, true>(
        data.data(), num_chunks, results.data(), bases);
    total_hits += hits;
  }
  (void)total_hits;
  (void)thread_id;
}

int main(int argc, char **argv) {
  int num_threads = (argc > 1) ? std::stoi(argv[1]) : 1;
  size_t chunks_per_thread = 100000;
  int iterations = 100;

  std::vector<std::thread> threads;
  auto t0 = std::chrono::high_resolution_clock::now();
  for (int i = 0; i < num_threads; ++i) {
    threads.emplace_back([&, i]() {
      bench_thread(i, chunks_per_thread, iterations);
    });
  }
  for (auto &t : threads)
    t.join();
  auto t1 = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> wall = t1 - t0;

  double total_bases =
      (double)num_threads * (double)chunks_per_thread * 128.0 * (double)iterations;
  double aggregate = (total_bases / 1e9) / wall.count();
  std::cout << std::fixed << std::setprecision(4);
  std::cout << "RESULT aggregate_gbp_s " << aggregate << "\n";
  std::cout << "RESULT wall_seconds " << wall.count() << "\n";
  std::cout << "RESULT threads " << num_threads << "\n";
  return 0;
}
