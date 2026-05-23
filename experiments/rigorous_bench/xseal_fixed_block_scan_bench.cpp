// SPDX-License-Identifier: Apache-2.0
// Copyright (c) 2026 hanrai. All Rights Reserved.

// Hot-RAM benchmark: 128 KiB 2-bit packed byte blocks -> worker pool.
// Per block: scan (overwrite thread-local pos) -> harvest immediately (hash mode) or hit count.
// No cross-block accumulation of hits; thread-local counters padded to avoid false sharing.
//
// Usage:
//   xseal_block128k_bench <hg38_2bit.bin> <threads> <hash_mode=0|1> [sync_k=31] [s=11] [min_w=21]

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <new>
#include <numeric>
#include <string>
#include <barrier>
#include <thread>
#include <vector>

#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>

#include <xseal/xseal_kmer_harvester.hpp>
#include <xseal/xseal_minimizer.hpp>
#include <xseal/xseal_syncmer.hpp>

using namespace xseal;

namespace {

constexpr size_t kBlockPackedBytes = 128u * 1024u;
constexpr size_t kSimdChunkBases = 128u;
constexpr size_t kCacheLine = 64u;
// Upper bound on hits per 128 KiB block (~0.14 density @ minimizer).
constexpr size_t kPosCapPerBlock = 512u * 1024u;

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

struct BlockGeom {
  size_t num_blocks;
  size_t tail_bytes;
};

struct SharedJob {
  const uint8_t *data;
  size_t file_bytes;
  size_t num_bases;
  BlockGeom geom;
  Scheme scheme;
  int sync_k;
  int s;
  int min_k;
  int min_w;
  int harvest_k;
  bool harvest_hash;
  std::atomic<size_t> next_block{0};
};

// Per-thread stats in separate cache lines (avoid false sharing on hit/digest updates).
struct alignas(kCacheLine) ThreadTally {
  size_t hits;
  uint64_t digest;
};

struct WorkerLocal {
  uint32_t *pos = nullptr;
  uint64_t *hash = nullptr;
  XSealSyncmer syncmer{31, 11};
  XSealMinimizer minimizer{31, 21};
};

void fold_digest_pos(uint64_t &d, const uint32_t *pos, size_t n) {
  for (size_t i = 0; i < n; ++i) {
    d += static_cast<uint64_t>(pos[i]);
    d = (d << 17) | (d >> (64 - 17));
  }
}

void fold_digest_pos_hash(uint64_t &d, const uint32_t *pos, const uint64_t *hash, size_t n) {
  for (size_t i = 0; i < n; ++i) {
    d += static_cast<uint64_t>(pos[i]) ^ hash[i];
    d = (d << 17) | (d >> (64 - 17));
  }
}

void touch_hot_pages(const uint8_t *data, size_t nbytes) {
  uint64_t touch = 0;
  for (size_t i = 0; i < nbytes; i += 4096) {
    touch += data[i];
  }
  asm volatile("" : "+r"(touch) : : "memory");
}

uint8_t *load_hot_packed(const char *path, size_t *file_bytes_out, size_t *num_bases_out) {
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
  *file_bytes_out = static_cast<size_t>(st.st_size);
  *num_bases_out = static_cast<size_t>(st.st_size) * 4;
  return data;
}

BlockGeom block_geometry(size_t file_bytes) {
  BlockGeom g{};
  g.tail_bytes = file_bytes % kBlockPackedBytes;
  g.num_blocks = file_bytes / kBlockPackedBytes + (g.tail_bytes ? 1u : 0u);
  return g;
}

void block_layout(const BlockGeom &geom, size_t block_id, size_t file_bytes, size_t *byte_off,
                  size_t *block_bytes, size_t *scan_bases) {
  *byte_off = block_id * kBlockPackedBytes;
  if (block_id + 1 < geom.num_blocks) {
    *block_bytes = kBlockPackedBytes;
  } else {
    *block_bytes = geom.tail_bytes ? geom.tail_bytes : kBlockPackedBytes;
  }
  const size_t avail_bytes = (*byte_off < file_bytes) ? (file_bytes - *byte_off) : 0;
  *block_bytes = std::min(*block_bytes, avail_bytes);
  const size_t raw_bases = (*block_bytes) * 4;
  *scan_bases = (raw_bases / kSimdChunkBases) * kSimdChunkBases;
}

template <int W>
size_t scan_block(Scheme scheme, XSealSyncmer &syncmer, XSealMinimizer &minimizer,
                  const __m256i *chunk_ptr, size_t num_chunks, size_t scan_bases,
                  uint32_t *out_pos) {
  switch (scheme) {
    case Scheme::ClosedSync:
      return syncmer.scan<W, 0, true, true>(chunk_ptr, num_chunks, out_pos, scan_bases);
    case Scheme::OpenSync:
      return syncmer.scan<W, 10, false, true>(chunk_ptr, num_chunks, out_pos, scan_bases);
    case Scheme::Minimizer:
      return minimizer.scan<W, true>(chunk_ptr, num_chunks, out_pos);
  }
  return 0;
}

// Scan then harvest on the same block buffer; pos/hash storage is overwritten each block.
size_t process_one_block(const SharedJob &job, uint32_t *pos, uint64_t *hash,
                         XSealSyncmer &syncmer, XSealMinimizer &minimizer, size_t block_id) {
  size_t byte_off = 0;
  size_t block_bytes = 0;
  size_t scan_bases = 0;
  block_layout(job.geom, block_id, job.file_bytes, &byte_off, &block_bytes, &scan_bases);
  if (scan_bases < 256) {
    return 0;
  }

  const size_t num_chunks = scan_bases / kSimdChunkBases;
  const __m256i *chunk_ptr = reinterpret_cast<const __m256i *>(job.data + byte_off);
  const uint8_t *enc = job.data + byte_off;

  size_t n = 0;
  if (job.scheme != Scheme::Minimizer) {
    n = scan_block<21>(job.scheme, syncmer, minimizer, chunk_ptr, num_chunks, scan_bases,
                       pos);
  } else {
    n = minimizer.scan<21, true>(chunk_ptr, num_chunks, pos);
  }
  if (n == 0) {
    return 0;
  }
  if (job.harvest_hash) {
    xseal_harvest_kmers<true>(enc, pos, n, job.harvest_k, hash);
  }
  asm volatile("" : : "g"(n), "g"(pos), "g"(hash) : "memory");
  return n;
}

struct RunResult {
  std::vector<double> recorded_throughputs;
  size_t last_hits;
  uint64_t last_digest;
};

RunResult run_scheme(SharedJob &job, std::vector<WorkerLocal> &workers,
                     std::vector<ThreadTally> &tallies, int internal_loops, int warmup_loops) {
  const int num_threads = static_cast<int>(workers.size());
  std::barrier<> sync_barrier(num_threads + 1);

  std::vector<std::thread> threads;
  threads.reserve(num_threads);
  for (int t = 0; t < num_threads; ++t) {
    threads.emplace_back([&job, &workers, &tallies, t, &sync_barrier, internal_loops]() {
      // Create stack-local variables to avoid heap/pointer-escape aliasing fences
      uint32_t *const local_pos = workers[t].pos;
      uint64_t *const local_hash = workers[t].hash;
      XSealSyncmer local_syncmer(job.sync_k, job.s);
      XSealMinimizer local_minimizer(job.min_k, job.min_w);

      for (int iter = 0; iter < internal_loops; ++iter) {
        // Wait for start signal
        sync_barrier.arrive_and_wait();

        // Thread resets its tally
        tallies[t].hits = 0;
        tallies[t].digest = 0;

        // Do work
        while (true) {
          const size_t block_id = job.next_block.fetch_add(1, std::memory_order_relaxed);
          if (block_id >= job.geom.num_blocks) {
            break;
          }
          const size_t n = process_one_block(job, local_pos, local_hash, local_syncmer, local_minimizer, block_id);
          if (n == 0) {
            continue;
          }
          tallies[t].hits += n;
          if (job.harvest_hash) {
            fold_digest_pos_hash(tallies[t].digest, local_pos, local_hash, n);
          } else {
            fold_digest_pos(tallies[t].digest, local_pos, n);
          }
        }

        // Wait for end signal
        sync_barrier.arrive_and_wait();
      }
    });
  }

  std::vector<double> recorded;
  recorded.reserve(internal_loops - warmup_loops);
  size_t last_hits = 0;
  uint64_t last_digest = 0;

  for (int iter = 0; iter < internal_loops; ++iter) {
    job.next_block.store(0, std::memory_order_relaxed);

    // Wait for all workers to be ready to start
    auto t0 = std::chrono::high_resolution_clock::now();
    sync_barrier.arrive_and_wait();

    // Wait for all workers to finish processing
    sync_barrier.arrive_and_wait();
    auto t1 = std::chrono::high_resolution_clock::now();

    const double secs = std::chrono::duration<double>(t1 - t0).count();
    const double throughput = static_cast<double>(job.num_bases) / secs / 1e9;

    if (iter >= warmup_loops) {
      recorded.push_back(throughput);
    }

    size_t total_hits = 0;
    uint64_t digest = 0;
    for (const ThreadTally &tl : tallies) {
      total_hits += tl.hits;
      digest += tl.digest;
    }
    last_hits = total_hits;
    last_digest = digest;
  }

  for (auto &th : threads) {
    th.join();
  }

  return {recorded, last_hits, last_digest};
}

bool init_workers(int num_threads, bool harvest_hash, std::vector<WorkerLocal> &workers,
                  std::vector<ThreadTally> &tallies) {
  workers.assign(static_cast<size_t>(num_threads), WorkerLocal{});
  tallies.assign(static_cast<size_t>(num_threads), ThreadTally{});

  for (int t = 0; t < num_threads; ++t) {
    workers[static_cast<size_t>(t)].pos = static_cast<uint32_t *>(
        aligned_alloc(64, kPosCapPerBlock * sizeof(uint32_t)));
    if (!workers[static_cast<size_t>(t)].pos) {
      return false;
    }
    if (harvest_hash) {
      workers[static_cast<size_t>(t)].hash = static_cast<uint64_t *>(
          aligned_alloc(64, kPosCapPerBlock * sizeof(uint64_t)));
      if (!workers[static_cast<size_t>(t)].hash) {
        return false;
      }
    }
  }
  return true;
}

void free_workers(std::vector<WorkerLocal> &workers) {
  for (auto &w : workers) {
    free(w.pos);
    free(w.hash);
    w.pos = nullptr;
    w.hash = nullptr;
  }
  workers.clear();
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

void print_scheme(const char *tool, Scheme scheme, int num_threads, size_t num_bases,
                  double mean, double inner_std, const RunStats &stats) {
  const double coverage = static_cast<double>(stats.hits) / static_cast<double>(num_bases);
  std::cout << std::fixed << std::setprecision(6);
  std::cout << "INNER|" << tool << "|" << scheme_name(scheme) << "|" << num_threads << "|" << mean
            << "|" << inner_std << "|" << stats.hits << "|" << coverage << "|"
            << std::hex << std::setw(16) << std::setfill('0') << stats.digest << std::dec << "\n";
  std::cout << std::setprecision(4);
  std::cout << "RESULT|" << tool << "|" << scheme_name(scheme) << "|" << num_threads << "|"
            << mean << "|" << inner_std << "|" << stats.hits << "|" << std::setprecision(6)
            << coverage << "|" << std::hex << stats.digest << std::dec << "\n";
}

}  // namespace

int main(int argc, char **argv) {
  if (argc < 4) {
    std::cerr << "Usage: " << (argc > 0 ? argv[0] : "xseal_block128k_bench")
              << " <hg38_2bit.bin> <threads> <hash_mode=0|1> [sync_k=31] [s=11] [min_w=21]\n";
    return 1;
  }

  const char *path = argv[1];
  const int num_threads = std::stoi(argv[2]);
  const int hash_mode = std::stoi(argv[3]);
  const int sync_k = (argc >= 5) ? std::stoi(argv[4]) : 31;
  const int s = (argc >= 6) ? std::stoi(argv[5]) : 11;
  const int min_w = (argc >= 7) ? std::stoi(argv[6]) : 21;
  const int w_sync = sync_k - s + 1;
  const int min_k = sync_k;
  const int harvest_k = sync_k;
  const bool harvest_hash = hash_mode != 0;

  size_t file_bytes = 0;
  size_t num_bases = 0;
  uint8_t *data = load_hot_packed(path, &file_bytes, &num_bases);
  const BlockGeom geom = block_geometry(file_bytes);

  const char *tool = harvest_hash ? "XsealBlock128k" : "XsealBlock128kCnt";

  std::vector<WorkerLocal> workers;
  std::vector<ThreadTally> tallies;
  if (!init_workers(num_threads, harvest_hash, workers, tallies)) {
    std::cerr << "worker buffer alloc failed\n";
    free(data);
    return 1;
  }

  std::cout << std::fixed;
  std::cout << "NOTE|" << tool
            << "|protocol=block128k_dispatch|block_packed_bytes=" << kBlockPackedBytes
            << "|scan_bases_per_full_block=" << (kBlockPackedBytes * 4)
            << "|num_blocks=" << geom.num_blocks << "|hash_mode=" << hash_mode
            << "|sync_k=" << sync_k << "|s=" << s << "|w_sync=" << w_sync
            << "|min_k=" << min_k << "|min_w=" << min_w << "|harvest_k=" << harvest_k
            << "|bases=" << num_bases << "|loops=" << kInternalLoops
            << "|warmup=" << kWarmupLoops << "|recorded=" << kRecordedLoops
            << "|tally_stride_bytes=" << sizeof(ThreadTally) << "\n";

  const Scheme schemes[] = {Scheme::ClosedSync, Scheme::OpenSync, Scheme::Minimizer};
  for (Scheme scheme : schemes) {
    SharedJob job{};
    job.data = data;
    job.file_bytes = file_bytes;
    job.num_bases = num_bases;
    job.geom = geom;
    job.scheme = scheme;
    job.sync_k = sync_k;
    job.s = s;
    job.min_k = min_k;
    job.min_w = min_w;
    job.harvest_k = harvest_k;
    job.harvest_hash = harvest_hash;

    RunResult result = run_scheme(job, workers, tallies, kInternalLoops, kWarmupLoops);
    const auto [mean, std] = mean_std(result.recorded_throughputs);
    RunStats out{mean, result.last_hits, result.last_digest};
    print_scheme(tool, scheme, num_threads, num_bases, mean, std, out);
  }

  free_workers(workers);
  free(data);
  return 0;
}
