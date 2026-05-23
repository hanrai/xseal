// SPDX-License-Identifier: Apache-2.0
// Copyright (c) 2026 hanrai. All Rights Reserved.

// Isolated s-mer hash micro-benchmark (XSeal): 2-bit extract + independent Murmur.
// Syncmer K=31, w=21 => s-mer length s=11. Compare with hash_isolate_simd (NtHash SIMD).

#include <algorithm>
#include <chrono>
#include <cstdint>
#include <cstdio>
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

#include <xseal/xseal_syncmer.hpp>

using namespace xseal;

namespace {

constexpr int kInternalLoops = 13;
constexpr int kWarmupLoops = 3;
constexpr int kRecordedLoops = kInternalLoops - kWarmupLoops;

struct RunStats {
  double throughput_gbp_s;
  size_t hash_count;
  uint64_t digest;
};

void touch_hot_pages(const uint8_t *data, size_t nbytes) {
  uint64_t touch = 0;
  for (size_t i = 0; i < nbytes; i += 4096)
    touch += data[i];
  asm volatile("" : "+r"(touch) : : "memory");
}

uint8_t *load_or_generate(const char *path, size_t num_bases, unsigned seed,
                           size_t *num_bases_out) {
  if (path && path[0] != '\0' && std::strcmp(path, "-") != 0) {
    int fd = open(path, O_RDONLY);
    if (fd < 0) {
      std::perror("open");
      std::exit(1);
    }
    struct stat st {};
    if (fstat(fd, &st) != 0) {
      std::perror("fstat");
      std::exit(1);
    }
    void *mmapped =
        mmap(nullptr, static_cast<size_t>(st.st_size), PROT_READ, MAP_PRIVATE, fd, 0);
    close(fd);
    if (mmapped == MAP_FAILED) {
      std::perror("mmap");
      std::exit(1);
    }
    size_t file_bases = static_cast<size_t>(st.st_size) * 4;
    size_t use_bases = num_bases > 0 ? std::min(num_bases, file_bases) : file_bases;
    size_t nbytes = (use_bases + 3) / 4;
    size_t padded = (nbytes + 1024 * 1024 + 63) & ~size_t(63);
    void *aligned_buf = nullptr;
    if (posix_memalign(&aligned_buf, 64, padded) != 0) {
      std::cerr << "posix_memalign failed\n";
      std::exit(1);
    }
    std::memset(aligned_buf, 0, padded);
    std::memcpy(aligned_buf, mmapped, nbytes);
    munmap(mmapped, static_cast<size_t>(st.st_size));
    auto *data = static_cast<uint8_t *>(aligned_buf);
    touch_hot_pages(data, nbytes);
    *num_bases_out = use_bases;
    return data;
  }

  size_t nbytes = (num_bases + 3) / 4;
  size_t padded = (nbytes + 63) & ~size_t(63);
  void *aligned_buf = nullptr;
  if (posix_memalign(&aligned_buf, 64, padded) != 0) {
    std::cerr << "posix_memalign failed\n";
    std::exit(1);
  }
  std::memset(aligned_buf, 0, padded);
  auto *data = static_cast<uint8_t *>(aligned_buf);
  std::mt19937_64 rng(seed);
  std::uniform_int_distribution<unsigned> dist(0, 255);
  for (size_t i = 0; i < nbytes; ++i)
    data[i] = static_cast<uint8_t>(dist(rng));
  touch_hot_pages(data, nbytes);
  *num_bases_out = num_bases;
  return data;
}

RunStats run_once(const uint8_t *data, size_t num_bases, int sync_k, int s,
                  int num_threads) {
  const size_t num_hashes =
      num_bases >= static_cast<size_t>(s) ? num_bases - s + 1 : 0;
  std::vector<uint64_t> digests(static_cast<size_t>(num_threads), 0);
  std::vector<size_t> counts(static_cast<size_t>(num_threads), 0);

  const auto t0 = std::chrono::high_resolution_clock::now();

  std::vector<std::thread> threads;
  threads.reserve(static_cast<size_t>(num_threads));
  for (int t = 0; t < num_threads; ++t) {
    threads.emplace_back([&, t]() {
      const size_t per =
          num_hashes / static_cast<size_t>(num_threads);
      const size_t start = static_cast<size_t>(t) * per;
      const size_t len =
          (t == num_threads - 1) ? (num_hashes - start) : per;
      if (len == 0)
        return;
      const size_t base_start = start;
      const size_t base_len = len + static_cast<size_t>(s) - 1;
      const size_t byte_off = base_start / 4;
      const size_t base_off = base_start % 4;
      XSealSyncmer scanner(sync_k, s);
      counts[static_cast<size_t>(t)] = scanner.hash_smers_only<true>(
          data + byte_off, base_len, &digests[static_cast<size_t>(t)]);
    });
  }
  for (auto &th : threads)
    th.join();

  const auto t1 = std::chrono::high_resolution_clock::now();
  double sec = std::chrono::duration<double>(t1 - t0).count();
  size_t total_hashes = 0;
  uint64_t digest = 0;
  for (int t = 0; t < num_threads; ++t) {
    total_hashes += counts[static_cast<size_t>(t)];
    digest ^= digests[static_cast<size_t>(t)];
  }
  double gbp_s = static_cast<double>(num_bases) / sec / 1e9;
  return {gbp_s, total_hashes, digest};
}

void usage(const char *argv0) {
  std::cerr
      << "Usage: " << argv0
      << " <2bit.bin|-> <num_bases|0> <threads> <sync_k> <s> [seed]\n"
      << "  '-' or omit file: generate random 2-bit (num_bases required).\n"
      << "  num_bases=0: use entire file. Default K=31 syncmer => s=11.\n"
      << "  Measures s-mer hash only (2-bit extract + Murmur), not window logic.\n";
}

} // namespace

int main(int argc, char **argv) {
  if (argc < 6) {
    usage(argv[0]);
    return 1;
  }
  const char *path = argv[1];
  const size_t num_bases_arg = static_cast<size_t>(std::stoull(argv[2]));
  const int num_threads = std::stoi(argv[3]);
  const int sync_k = std::stoi(argv[4]);
  const int s = std::stoi(argv[5]);
  const unsigned seed =
      argc >= 7 ? static_cast<unsigned>(std::stoul(argv[6])) : 42u;

  if (num_threads < 1 || s < 1 || sync_k < s) {
    usage(argv[0]);
    return 1;
  }

  size_t num_bases = num_bases_arg;
  uint8_t *data = load_or_generate(
      (path[0] == '-' && path[1] == '\0') ? nullptr : path, num_bases, seed,
      &num_bases);

  std::cout << std::fixed << std::setprecision(6);
  std::cout << "NOTE|tool|XsealHashIsolate|sync_k|" << sync_k << "|s|" << s
            << "|num_bases|" << num_bases << "|threads|" << num_threads
            << "|seed|" << seed << "\n";

  for (int loop = 0; loop < kInternalLoops; ++loop) {
    RunStats st = run_once(data, num_bases, sync_k, s, num_threads);
    if (loop >= kWarmupLoops) {
      std::cout << "INNER|loop|" << loop << "|gbp_s|" << st.throughput_gbp_s
                << "|hashes|" << st.hash_count << "|digest|0x" << std::hex
                << st.digest << std::dec << "\n";
    }
  }

  std::vector<double> recorded;
  recorded.reserve(kRecordedLoops);
  uint64_t final_digest = 0;
  size_t final_hashes = 0;
  for (int loop = 0; loop < kRecordedLoops; ++loop) {
    RunStats st = run_once(data, num_bases, sync_k, s, num_threads);
    recorded.push_back(st.throughput_gbp_s);
    final_digest = st.digest;
    final_hashes = st.hash_count;
  }
  double mean = 0;
  for (double v : recorded)
    mean += v;
  mean /= static_cast<double>(recorded.size());
  double var = 0;
  for (double v : recorded) {
    double d = v - mean;
    var += d * d;
  }
  var /= static_cast<double>(recorded.size());
  double stddev = std::sqrt(var);

  std::cout << "RESULT|XsealHashIsolate|" << num_threads << "|" << mean << "|"
            << stddev << "|" << num_bases << "|" << final_hashes << "|0x"
            << std::hex << final_digest << std::dec << "\n";

  free(data);
  return 0;
}
