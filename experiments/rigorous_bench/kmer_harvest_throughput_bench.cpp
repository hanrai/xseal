// SPDX-License-Identifier: Apache-2.0
// Copyright (c) 2026 hanrai. All Rights Reserved.

// Micro-benchmark: mmap 2-bit stream + xseal_harvest_kmers<true> (see xseal/xseal_kmer_harvester.hpp).
// Optional outer timing: hyperfine; this binary also does a small fixed warmup/trial loop for console stats.

#include <algorithm>
#include <chrono>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <random>
#include <string>
#include <thread>
#include <vector>

#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>

#include <xseal/xseal_kmer_harvester.hpp>

static void usage(const char *argv0) {
  std::cerr << "Usage: " << argv0 << " <2bit.bin> <k> <N_positions> <threads>\n"
            << "  k in [1,31], N = total positions (split across threads).\n"
            << "  Positions are uniform random in [0, num_bases-k) with seed=1.\n";
}

int main(int argc, char **argv) {
  if (argc != 5) {
    usage(argv[0]);
    return 1;
  }

  const char *path = argv[1];
  const int k = std::stoi(argv[2]);
  const size_t N = static_cast<size_t>(std::stoull(argv[3]));
  const int num_threads = std::stoi(argv[4]);

  if (k < 1 || k > 31 || num_threads < 1) {
    usage(argv[0]);
    return 1;
  }
  if (N == 0) {
    std::cerr << "N_positions must be > 0\n";
    return 1;
  }

  int fd = open(path, O_RDONLY);
  if (fd < 0) {
    perror("open");
    return 1;
  }
  struct stat st {};
  if (fstat(fd, &st) != 0) {
    perror("fstat");
    close(fd);
    return 1;
  }
  if (st.st_size <= 0) {
    std::cerr << "empty file\n";
    close(fd);
    return 1;
  }

  void *mapped =
      mmap(nullptr, static_cast<size_t>(st.st_size), PROT_READ, MAP_PRIVATE, fd, 0);
  close(fd);
  if (mapped == MAP_FAILED) {
    perror("mmap");
    return 1;
  }

  const uint8_t *enc = static_cast<const uint8_t *>(mapped);
  const size_t num_bases = static_cast<size_t>(st.st_size) * 4u;
  if (num_bases <= static_cast<size_t>(k)) {
    std::cerr << "genome too short for k\n";
    munmap(mapped, static_cast<size_t>(st.st_size));
    return 1;
  }
  if (num_bases > static_cast<size_t>(UINT32_MAX)) {
    std::cerr << "num_bases > UINT32_MAX: position API is uint32_t; use a smaller 2bit slice\n";
    munmap(mapped, static_cast<size_t>(st.st_size));
    return 1;
  }

  const uint32_t span = static_cast<uint32_t>(num_bases - static_cast<size_t>(k));
  if (span == 0) {
    munmap(mapped, static_cast<size_t>(st.st_size));
    return 1;
  }

  std::vector<uint32_t> pos(N);
  std::vector<uint64_t> out(N);
  {
    std::mt19937_64 gen(1);
    std::uniform_int_distribution<uint32_t> dist(0, span > 0 ? span - 1 : 0);
    for (size_t i = 0; i < N; ++i)
      pos[i] = dist(gen);
    std::sort(pos.begin(), pos.end());
  }

  const int warmups = 3;
  for (int w = 0; w < warmups; ++w) {
    std::vector<std::thread> workers;
    workers.reserve(static_cast<size_t>(num_threads));
    for (int t = 0; t < num_threads; ++t) {
      const size_t chunk = N / static_cast<size_t>(num_threads);
      const size_t start = static_cast<size_t>(t) * chunk;
      const size_t len = (t == num_threads - 1) ? (N - start) : chunk;
      workers.emplace_back([=, &pos, &out]() {
        xseal::xseal_harvest_kmers<true>(enc, pos.data() + start, len, k, out.data() + start);
      });
    }
    for (auto &worker : workers) worker.join();
  }

  const int trials = 10;
  std::vector<double> results_sec;
  results_sec.reserve(trials);

  for (int trial = 0; trial < trials; ++trial) {
    std::vector<std::thread> workers;
    workers.reserve(static_cast<size_t>(num_threads));

    const auto t0 = std::chrono::high_resolution_clock::now();

    for (int t = 0; t < num_threads; ++t) {
      const size_t chunk = N / static_cast<size_t>(num_threads);
      const size_t start = static_cast<size_t>(t) * chunk;
      const size_t len = (t == num_threads - 1) ? (N - start) : chunk;
      workers.emplace_back([=, &pos, &out]() {
        xseal::xseal_harvest_kmers<true>(enc, pos.data() + start, len, k, out.data() + start);
      });
    }
    for (auto &w : workers)
      w.join();

    const auto t1 = std::chrono::high_resolution_clock::now();
    results_sec.push_back(std::chrono::duration<double>(t1 - t0).count());
    
    // Optional: verification of some results to prevent optimization
    volatile uint64_t sink = 0;
    for (size_t i = 0; i < N; i += N / 100 + 1) sink ^= out[i];
  }

  double sum = 0;
  for (double s : results_sec) sum += s;
  double avg_sec = sum / trials;

  std::cout << std::fixed << std::setprecision(6);
  const double g_hashes = static_cast<double>(N) / avg_sec / 1e9;
  uint64_t xor_ck = 0;
  for (size_t i = 0; i < N; ++i)
    xor_ck ^= out[i];
  std::cout << "\n[HARVESTER BENCHMARK RESULTS]\n";
  std::cout << "Threads:    " << num_threads << "\n";
  std::cout << "K:          " << k << "\n";
  std::cout << "Positions:  " << N << "\n";
  std::cout << "Avg Time:   " << avg_sec << " s\n";
  std::cout << "Throughput: " << g_hashes << " G-hashes/s\n";
  std::cout << "----------------------------\n";
  // Machine-parseable; xor_ck should match across thread counts for same N/k/seed/sorted pos.
  std::cout << "RESULT|HarvestSplitmix|" << num_threads << "|" << avg_sec << "|" << N << "|" << g_hashes << "|"
            << std::hex << xor_ck << std::dec << "\n";

  munmap(mapped, static_cast<size_t>(st.st_size));
  return 0;
}
