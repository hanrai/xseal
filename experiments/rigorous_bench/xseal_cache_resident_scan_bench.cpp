// SPDX-License-Identifier: Apache-2.0
// Copyright (c) 2026 hanrai. All Rights Reserved.

// Hot-RAM 2-bit packed benchmark: ClosedSync, OpenSync, Minimizer (XSeal).
// Protocol matches simd_min_hot_bench.

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <string>
#include <barrier>
#include <thread>
#include <vector>

#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>

#include <xseal/xseal_minimizer.hpp>
#include <xseal/xseal_syncmer.hpp>

using namespace xseal;

namespace {

constexpr int kInternalLoops = 13;
constexpr int kWarmupLoops = 3;
constexpr int kRecordedLoops = kInternalLoops - kWarmupLoops;

enum class Scheme { ClosedSync, OpenSync, Minimizer };

const char *scheme_name(Scheme s) {
  switch (s) {
    case Scheme::ClosedSync:
      return "ClosedSync";
    case Scheme::OpenSync:
      return "OpenSync";
    case Scheme::Minimizer:
      return "Minimizer";
  }
  return "Unknown";
}

struct RunStats {
  double throughput_gbp_s;
  size_t hits;
  uint64_t digest;
};

uint64_t digest_positions(const uint32_t *pos, size_t hits) {
  uint64_t d = 0;
  for (size_t i = 0; i < hits; ++i) {
    d += static_cast<uint64_t>(pos[i]);
    d = (d << 17) | (d >> (64 - 17));
  }
  asm volatile("" : "+r"(d) : : "memory");
  return d;
}

void touch_hot_pages(const uint8_t *data, size_t nbytes) {
  uint64_t touch = 0;
  for (size_t i = 0; i < nbytes; i += 4096) {
    touch += data[i];
  }
  asm volatile("" : "+r"(touch) : : "memory");
}

uint8_t *load_hot_packed(const char *path, size_t *num_bases_out) {
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

  uint8_t *mmapped_ptr =
      static_cast<uint8_t *>(mmap(nullptr, st.st_size, PROT_READ, MAP_PRIVATE, fd, 0));
  close(fd);
  if (mmapped_ptr == MAP_FAILED) {
    std::perror("mmap");
    std::exit(1);
  }

  size_t padded_size = (static_cast<size_t>(st.st_size) + 1024 * 1024 + 63) & ~size_t(63);
  void *aligned_buf = nullptr;
  if (posix_memalign(&aligned_buf, 64, padded_size) != 0) {
    std::cerr << "posix_memalign failed\n";
    std::exit(1);
  }
  std::memset(aligned_buf, 0, padded_size);
  std::memcpy(aligned_buf, mmapped_ptr, st.st_size);
  munmap(mmapped_ptr, st.st_size);

  uint8_t *data = static_cast<uint8_t *>(aligned_buf);
  touch_hot_pages(data, st.st_size);
  *num_bases_out = static_cast<size_t>(st.st_size) * 4;
  return data;
}

size_t pos_capacity_for_thread(size_t num_bases, int num_threads) {
  return std::max(static_cast<size_t>(1024 * 1024 * 64),
                  num_bases / static_cast<size_t>(num_threads) / 4 + 1024);
}

struct RunResult {
  std::vector<double> recorded_throughputs;
  size_t last_hits;
  uint64_t last_digest;
};

template <int W>
RunResult run_closed_sync(const uint8_t *data, size_t num_bases, int sync_k, int s,
                           int num_threads, uint32_t **out_pos, int internal_loops, int warmup_loops) {
  const size_t chunk_size_bases = (num_bases / num_threads / 128) * 128;
  std::barrier<> sync_barrier(num_threads + 1);
  std::vector<size_t> thread_hits(num_threads, 0);

  std::vector<std::thread> threads;
  threads.reserve(num_threads);
  for (int t = 0; t < num_threads; ++t) {
    threads.emplace_back([&, t]() {
      XSealSyncmer scanner(sync_k, s);
      const size_t start_base = static_cast<size_t>(t) * chunk_size_bases;
      const size_t thread_num_bases =
          (t == num_threads - 1) ? (num_bases - start_base) : chunk_size_bases;

      for (int iter = 0; iter < internal_loops; ++iter) {
        sync_barrier.arrive_and_wait();

        if (thread_num_bases >= 256) {
          const size_t num_chunks = thread_num_bases / 128;
          const __m256i *chunk_ptr =
              reinterpret_cast<const __m256i *>(data + (start_base / 4));
          thread_hits[t] = scanner.scan<W, 0, true, true>(
              chunk_ptr, num_chunks, out_pos[t], thread_num_bases);
        } else {
          thread_hits[t] = 0;
        }
        asm volatile("" : : "g"(thread_hits[t]), "g"(out_pos[t]) : "memory");

        sync_barrier.arrive_and_wait();
      }
    });
  }

  std::vector<double> recorded;
  recorded.reserve(internal_loops - warmup_loops);
  size_t last_hits = 0;
  uint64_t last_digest = 0;

  for (int iter = 0; iter < internal_loops; ++iter) {
    const size_t cap = pos_capacity_for_thread(num_bases, num_threads);
    for (int t = 0; t < num_threads; ++t) {
      std::memset(out_pos[t], 0, cap * sizeof(uint32_t));
    }

    auto t0 = std::chrono::high_resolution_clock::now();
    sync_barrier.arrive_and_wait();

    sync_barrier.arrive_and_wait();
    auto t1 = std::chrono::high_resolution_clock::now();

    const double secs = std::chrono::duration<double>(t1 - t0).count();
    const double throughput = static_cast<double>(num_bases) / secs / 1e9;

    if (iter >= warmup_loops) {
      recorded.push_back(throughput);
    }

    size_t total_hits = 0;
    uint64_t digest = 0;
    for (int t = 0; t < num_threads; ++t) {
      total_hits += thread_hits[t];
      digest += digest_positions(out_pos[t], thread_hits[t]);
    }
    last_hits = total_hits;
    last_digest = digest;
  }

  for (auto &th : threads) {
    th.join();
  }

  return {recorded, last_hits, last_digest};
}

template <int W, int T>
RunResult run_open_sync(const uint8_t *data, size_t num_bases, int sync_k, int s,
                         int num_threads, uint32_t **out_pos, int internal_loops, int warmup_loops) {
  const size_t chunk_size_bases = (num_bases / num_threads / 128) * 128;
  std::barrier<> sync_barrier(num_threads + 1);
  std::vector<size_t> thread_hits(num_threads, 0);

  std::vector<std::thread> threads;
  threads.reserve(num_threads);
  for (int t = 0; t < num_threads; ++t) {
    threads.emplace_back([&, t]() {
      XSealSyncmer scanner(sync_k, s);
      const size_t start_base = static_cast<size_t>(t) * chunk_size_bases;
      const size_t thread_num_bases =
          (t == num_threads - 1) ? (num_bases - start_base) : chunk_size_bases;

      for (int iter = 0; iter < internal_loops; ++iter) {
        sync_barrier.arrive_and_wait();

        if (thread_num_bases >= 256) {
          const size_t num_chunks = thread_num_bases / 128;
          const __m256i *chunk_ptr =
              reinterpret_cast<const __m256i *>(data + (start_base / 4));
          thread_hits[t] = scanner.scan<W, T, false, true>(
              chunk_ptr, num_chunks, out_pos[t], thread_num_bases);
        } else {
          thread_hits[t] = 0;
        }
        asm volatile("" : : "g"(thread_hits[t]), "g"(out_pos[t]) : "memory");

        sync_barrier.arrive_and_wait();
      }
    });
  }

  std::vector<double> recorded;
  recorded.reserve(internal_loops - warmup_loops);
  size_t last_hits = 0;
  uint64_t last_digest = 0;

  for (int iter = 0; iter < internal_loops; ++iter) {
    const size_t cap = pos_capacity_for_thread(num_bases, num_threads);
    for (int t = 0; t < num_threads; ++t) {
      std::memset(out_pos[t], 0, cap * sizeof(uint32_t));
    }

    auto t0 = std::chrono::high_resolution_clock::now();
    sync_barrier.arrive_and_wait();

    sync_barrier.arrive_and_wait();
    auto t1 = std::chrono::high_resolution_clock::now();

    const double secs = std::chrono::duration<double>(t1 - t0).count();
    const double throughput = static_cast<double>(num_bases) / secs / 1e9;

    if (iter >= warmup_loops) {
      recorded.push_back(throughput);
    }

    size_t total_hits = 0;
    uint64_t digest = 0;
    for (int t = 0; t < num_threads; ++t) {
      total_hits += thread_hits[t];
      digest += digest_positions(out_pos[t], thread_hits[t]);
    }
    last_hits = total_hits;
    last_digest = digest;
  }

  for (auto &th : threads) {
    th.join();
  }

  return {recorded, last_hits, last_digest};
}

template <int W>
RunResult run_minimizer(const uint8_t *data, size_t num_bases, int min_k, int min_w,
                         int num_threads, uint32_t **out_pos, int internal_loops, int warmup_loops) {
  const size_t chunk_size_bases = (num_bases / num_threads / 128) * 128;
  std::barrier<> sync_barrier(num_threads + 1);
  std::vector<size_t> thread_hits(num_threads, 0);

  std::vector<std::thread> threads;
  threads.reserve(num_threads);
  for (int t = 0; t < num_threads; ++t) {
    threads.emplace_back([&, t]() {
      XSealMinimizer scanner(min_k, min_w);
      const size_t start_base = static_cast<size_t>(t) * chunk_size_bases;
      const size_t thread_num_bases =
          (t == num_threads - 1) ? (num_bases - start_base) : chunk_size_bases;

      for (int iter = 0; iter < internal_loops; ++iter) {
        sync_barrier.arrive_and_wait();

        if (thread_num_bases >= 256) {
          const size_t num_chunks = thread_num_bases / 128;
          const __m256i *chunk_ptr =
              reinterpret_cast<const __m256i *>(data + (start_base / 4));
          thread_hits[t] = scanner.scan<W, true>(chunk_ptr, num_chunks, out_pos[t]);
        } else {
          thread_hits[t] = 0;
        }
        asm volatile("" : : "g"(thread_hits[t]), "g"(out_pos[t]) : "memory");

        sync_barrier.arrive_and_wait();
      }
    });
  }

  std::vector<double> recorded;
  recorded.reserve(internal_loops - warmup_loops);
  size_t last_hits = 0;
  uint64_t last_digest = 0;

  for (int iter = 0; iter < internal_loops; ++iter) {
    const size_t cap = pos_capacity_for_thread(num_bases, num_threads);
    for (int t = 0; t < num_threads; ++t) {
      std::memset(out_pos[t], 0, cap * sizeof(uint32_t));
    }

    auto t0 = std::chrono::high_resolution_clock::now();
    sync_barrier.arrive_and_wait();

    sync_barrier.arrive_and_wait();
    auto t1 = std::chrono::high_resolution_clock::now();

    const double secs = std::chrono::duration<double>(t1 - t0).count();
    const double throughput = static_cast<double>(num_bases) / secs / 1e9;

    if (iter >= warmup_loops) {
      recorded.push_back(throughput);
    }

    size_t total_hits = 0;
    uint64_t digest = 0;
    for (int t = 0; t < num_threads; ++t) {
      total_hits += thread_hits[t];
      digest += digest_positions(out_pos[t], thread_hits[t]);
    }
    last_hits = total_hits;
    last_digest = digest;
  }

  for (auto &th : threads) {
    th.join();
  }

  return {recorded, last_hits, last_digest};
}

RunResult dispatch_scheme(Scheme scheme, int w_sync, const uint8_t *data, size_t num_bases,
                         int sync_k, int s, int min_k, int min_w, int num_threads,
                         uint32_t **out_pos, int internal_loops, int warmup_loops) {
  switch (scheme) {
    case Scheme::ClosedSync:
      return run_closed_sync<21>(data, num_bases, sync_k, s, num_threads, out_pos, internal_loops, warmup_loops);
    case Scheme::OpenSync:
      return run_open_sync<21, 10>(data, num_bases, sync_k, s, num_threads, out_pos, internal_loops, warmup_loops);
    case Scheme::Minimizer:
      return run_minimizer<21>(data, num_bases, min_k, min_w, num_threads, out_pos, internal_loops, warmup_loops);
  }
  return {};
}

std::pair<double, double> mean_std(const std::vector<double> &v) {
  const double n = static_cast<double>(v.size());
  const double mean = std::accumulate(v.begin(), v.end(), 0.0) / n;
  double var = 0.0;
  for (double x : v) {
    var += (x - mean) * (x - mean);
  }
  var /= n;
  return {mean, std::sqrt(var)};
}

void print_scheme(Scheme scheme, int num_threads, size_t num_bases, double mean, double inner_std,
                  const RunStats &stats) {
  const double coverage = static_cast<double>(stats.hits) / static_cast<double>(num_bases);
  std::cout << std::fixed << std::setprecision(6);
  std::cout << "INNER|XsealHot|" << scheme_name(scheme) << "|" << num_threads << "|" << mean
            << "|" << inner_std << "|" << stats.hits << "|"
            << coverage << "|" << std::hex << std::setw(16) << std::setfill('0')
            << stats.digest << std::dec << "\n";
  std::cout << std::setprecision(4);
  std::cout << "RESULT|XsealHot|" << scheme_name(scheme) << "|" << num_threads << "|" << mean
            << "|" << inner_std << "|" << stats.hits << "|"
            << std::setprecision(6) << coverage << "|" << std::hex << stats.digest << std::dec
            << "\n";
}

}  // namespace

int main(int argc, char **argv) {
  if (argc < 3) {
    std::cerr << "Usage: " << (argc > 0 ? argv[0] : "xseal_hot_bench")
              << " <hg38_2bit.bin> <threads> [sync_k=31] [s=11] [min_w=21] [mode=0 (All)|1 (Closed)|2 (Open)|3 (Minimizer)]\n";
    return 1;
  }

  const char *path = argv[1];
  const int num_threads = std::stoi(argv[2]);
  const int sync_k = (argc >= 4) ? std::stoi(argv[3]) : 31;
  const int s = (argc >= 5) ? std::stoi(argv[4]) : 11;
  const int min_w = (argc >= 6) ? std::stoi(argv[5]) : 21;
  const int run_mode = (argc >= 7) ? std::stoi(argv[6]) : 0;
  const int w_sync = sync_k - s + 1;
  const int min_k = sync_k;

  size_t num_bases = 0;
  uint8_t *data = load_hot_packed(path, &num_bases);

  uint32_t **out_pos = new uint32_t *[num_threads];
  const size_t per_thread_buf_size = pos_capacity_for_thread(num_bases, num_threads);
  for (int t = 0; t < num_threads; ++t) {
    out_pos[t] = static_cast<uint32_t *>(
        aligned_alloc(64, (per_thread_buf_size + 256) * sizeof(uint32_t)));
    if (!out_pos[t]) {
      std::cerr << "aligned_alloc failed\n";
      return 1;
    }
  }

  std::cout << std::fixed;
  std::cout << "NOTE|XsealHot|sync_k=" << sync_k << "|s=" << s << "|w_sync=" << w_sync
            << "|min_k=" << min_k << "|min_w=" << min_w << "|bases=" << num_bases
            << "|loops=" << kInternalLoops << "|warmup=" << kWarmupLoops
            << "|recorded=" << kRecordedLoops << "\n";

  std::vector<Scheme> schemes;
  if (run_mode == 0) {
    schemes = {Scheme::ClosedSync, Scheme::OpenSync, Scheme::Minimizer};
  } else if (run_mode == 1) {
    schemes = {Scheme::ClosedSync};
  } else if (run_mode == 2) {
    schemes = {Scheme::OpenSync};
  } else if (run_mode == 3) {
    schemes = {Scheme::Minimizer};
  } else {
    std::cerr << "Invalid mode\n";
    return 1;
  }

  for (Scheme scheme : schemes) {
    RunResult result = dispatch_scheme(scheme, w_sync, data, num_bases, sync_k, s, min_k,
                                       min_w, num_threads, out_pos, kInternalLoops, kWarmupLoops);

    const auto [mean, std] = mean_std(result.recorded_throughputs);
    RunStats out{mean, result.last_hits, result.last_digest};
    print_scheme(scheme, num_threads, num_bases, mean, std, out);
  }

  for (int t = 0; t < num_threads; ++t) {
    free(out_pos[t]);
  }
  delete[] out_pos;
  free(data);
  return 0;
}
