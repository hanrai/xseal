// SPDX-License-Identifier: Apache-2.0
// Copyright (c) 2026 hanrai. All Rights Reserved.

#pragma once

#include <atomic>
#include <chrono>
#include <cstring>
#include <immintrin.h>
#include <thread>
#include <vector>
#include <xseal/xseal_encoder.hpp>
#include <xseal/xseal_fasta_parser.hpp>
#include <xseal/xseal_minimizer.hpp>
#include <xseal/xseal_syncmer.hpp>
#include <xseal/xseal_slot_pool.hpp>
#include <xseal/xseal_kmer_harvester.hpp>

namespace xseal {

/// Default raw-byte batch before a distributor task is pushed (historically 2MiB).
inline constexpr size_t kXsealDefaultDistributorBatchBytes = 2 * 1024 * 1024;

enum class XsealScannerType { SYNCMER, MINIMIZER };

struct XsealSegmentRef {
  const uint8_t *ptr;
  size_t len;
  uint64_t absolute_offset;
};

template <uint32_t NumSlots = 768> struct XsealTask {
  std::vector<XsealSegmentRef> segments;
  std::vector<uint8_t> carry_over;
  size_t total_raw_len = 0;
  uint64_t task_id = 0;
  std::vector<XsealSlot<NumSlots>> slots; // ✅ Upgraded to slot array
};

template <uint32_t NumSlots = 768> class XsealDistributor {
  /// Packed encoder output → syncmer/minimizer working buffer.
  static constexpr size_t kEncoderToScannerBufferCap = 16 * 1024 * 1024;

public:
  XsealDistributor(int num_threads, int k, int s,
                   XsealScannerType type = XsealScannerType::SYNCMER,
                   size_t batch_threshold = kXsealDefaultDistributorBatchBytes)
      : K(k), S(s), scanner_type(type), BATCH_THRESHOLD(batch_threshold) {
    finished.store(false, std::memory_order_relaxed);
    for (size_t i = 0; i < RING_SIZE; ++i) {
      task_ring[i].seq.store(i, std::memory_order_relaxed);
      task_done_flags[i].store(0, std::memory_order_relaxed);
    }
    for (int i = 0; i < num_threads; ++i) {
      workers.emplace_back(&XsealDistributor::worker_loop, this);
    }
  }

  ~XsealDistributor() { stop(); }

  void stop() {
    flush_current_batch();
    finished.store(true, std::memory_order_release);
    for (auto &t : workers) {
      if (t.joinable())
        t.join();
    }
  }

  void wait_for_task(uint64_t target_id) {
    if (target_id == 0)
      return;
    while (task_done_flags[target_id % RING_SIZE].load(
               std::memory_order_acquire) != target_id &&
           !finished.load(std::memory_order_acquire)) {
      _mm_pause();
    }
  }

  uint64_t get_tasks_pushed_count() const { return tasks_pushed_count.load(); }

  void process_bases(const SeqSegment &seg, const uint8_t *chunk_base,
                     uint64_t chunk_off, const XsealSlot<NumSlots> &slot = XsealSlot<NumSlots>()) {
    // ✅ Core fix: Hold the lifecycle of all cross-chunk memory
    if (slot.valid()) {
      bool already_have = false;
      for (const auto &s : current_task.slots) {
        if (s.id() == slot.id()) {
          already_have = true;
          break;
        }
      }
      if (!already_have)
        current_task.slots.push_back(slot);
    }
    XsealSegmentRef ref;
    ref.ptr = chunk_base + (seg.file_offset - chunk_off);
    ref.len = seg.byte_len;
    ref.absolute_offset = seg.file_offset;

    if (current_task.segments.empty() && !saved_tail.empty()) {
      current_task.carry_over = std::move(saved_tail);
      saved_tail.clear();
    }

    // 2. Handle physical truncation with "Reverse Extraction" to avoid the
    // Dirty Tail Problem
    if (seg.is_physically_truncated) {
      size_t valid_bases_collected = 0;
      size_t tail_bytes_to_extract = 0;
      const uint8_t *p = ref.ptr + ref.len - 1;

      while (p >= ref.ptr && valid_bases_collected < (size_t)(K - 1)) {
        if (*p != '\n' && *p != '\r') {
          valid_bases_collected++;
        }
        tail_bytes_to_extract++;
        p--;
      }

      if (valid_bases_collected == (size_t)(K - 1)) {
        saved_tail.assign(ref.ptr + ref.len - tail_bytes_to_extract,
                          ref.ptr + ref.len);
      } else {
        saved_tail.insert(saved_tail.end(), ref.ptr, ref.ptr + ref.len);
      }
    }

    if (ref.len > 0) {
      current_task.total_raw_len += ref.len;
      current_task.segments.push_back(ref);
    }

    if (current_task.total_raw_len >= BATCH_THRESHOLD) {
        flush_current_batch();
    }
  }

  void reset_sequence() {
    flush();
    saved_tail.clear();
  }

  void flush() { flush_current_batch(); }

  size_t get_total_hits() const { return total_hits.load(); }
  size_t get_total_bases() const { return total_bases.load(); }
  uint64_t get_total_scan_ns() const { return total_scan_ns.load(); }
  uint64_t get_total_encode_ns() const { return total_encode_ns.load(); }
  uint64_t get_total_harvest_ns() const { return total_harvest_ns.load(); }

protected:
  void flush_current_batch() {
    if (current_task.segments.empty() && current_task.carry_over.empty())
      return;

    uint64_t idx = ring_push_idx.load(std::memory_order_relaxed);
    while (task_ring[idx % RING_SIZE].seq.load(std::memory_order_acquire) !=
           idx) {
      _mm_pause();
    }

    current_task.task_id = ++tasks_pushed_count;
    // ✅ Copy all currently held slots, pass to the next Task to continue lifecycle

    task_ring[idx % RING_SIZE].task = std::move(current_task);
    task_ring[idx % RING_SIZE].seq.store(idx + 1, std::memory_order_release);
    ring_push_idx.store(idx + 1, std::memory_order_relaxed);

    current_task = XsealTask<NumSlots>();
  }

  void worker_loop() {
    XSealEncoder encoder;
    XSealEncoderState enc_state;
    XSealSyncmer syncmer_scanner(K, S);
    XSealMinimizer minimizer_scanner(K, K - S + 1); // W = K - S + 1

    std::vector<uint8_t> out_buf(kEncoderToScannerBufferCap + 64);
    std::vector<uint32_t> hit_pos_pool(1 * 1024 * 1024); // 4MB enough
    std::vector<uint64_t> hit_hash_pool(1 * 1024 * 1024); // 8MB enough
    
    uint64_t local_scan_ns = 0;
    uint64_t local_encode_ns = 0;
    uint64_t local_harvest_ns = 0;
    size_t local_hits = 0;
    size_t local_bases = 0;

    while (true) {
      uint64_t idx = ring_consume_idx.fetch_add(1, std::memory_order_relaxed);
      while (task_ring[idx % RING_SIZE].seq.load(std::memory_order_acquire) !=
             idx + 1) {
        if (finished.load(std::memory_order_relaxed) &&
            idx >= tasks_pushed_count.load(std::memory_order_relaxed)) {
          // Summary update before exit
          total_scan_ns += local_scan_ns;
          total_encode_ns += local_encode_ns;
          total_harvest_ns += local_harvest_ns;
          total_hits += local_hits;
          total_bases += local_bases;
          return;
        }
        _mm_pause();
      }
      XsealTask<NumSlots> task = std::move(task_ring[idx % RING_SIZE].task);
      task_ring[idx % RING_SIZE].seq.store(idx + RING_SIZE,
                                           std::memory_order_release);

      encoder.init_state(&enc_state, out_buf.data(), kEncoderToScannerBufferCap);
      auto t_encode_0 = std::chrono::high_resolution_clock::now();

      // Calculate actual bases in carry_over for coordinate adjustment
      uint64_t carry_offset_adj = 0;
      for (uint8_t c : task.carry_over) {
        if (c != '\n' && c != '\r')
          carry_offset_adj++;
      }

      // 1. Encode carry-over
      if (!task.carry_over.empty()) {
        const char *p = (const char *)task.carry_over.data();
        const char *end = p + task.carry_over.size();
        encoder.encode_chunk(&enc_state, &p, end);
      }

      // 2. Encode all segments in batch
      for (const auto &seg : task.segments) {
        const char *p = (const char *)seg.ptr;
        const char *end = p + seg.len;
        encoder.encode_chunk(&enc_state, &p, end);
      }

      encoder.flush_state(&enc_state);
      local_encode_ns += std::chrono::duration_cast<std::chrono::nanoseconds>(
                             std::chrono::high_resolution_clock::now() - t_encode_0)
                             .count();

      // 3. Scan
      if (enc_state.total_bases_encoded >= (size_t)K) {
        _mm256_storeu_si256(
            (__m256i *)(enc_state.out_buffer + enc_state.out_byte_idx),
            _mm256_setzero_si256());

        size_t num_chunks = (enc_state.out_byte_idx + 31) / 32;
        size_t hits = 0;

        auto t_scan_0 = std::chrono::high_resolution_clock::now();
        if (scanner_type == XsealScannerType::SYNCMER) {
          hits = syncmer_scanner.scan<21, 0, true, true>(
              reinterpret_cast<const __m256i *>(enc_state.out_buffer),
              num_chunks, hit_pos_pool.data(), enc_state.total_bases_encoded);
        } else {
          hits = minimizer_scanner.scan<21, true>(
              reinterpret_cast<const __m256i *>(enc_state.out_buffer),
              num_chunks, hit_pos_pool.data(), nullptr);
        }
        auto t_scan_1 = std::chrono::high_resolution_clock::now();
        local_scan_ns += std::chrono::duration_cast<std::chrono::nanoseconds>(t_scan_1 - t_scan_0).count();

        // 4. Harvest 64-bit hashes
        if (hits > 0) {
          auto t_harvest_0 = std::chrono::high_resolution_clock::now();
          xseal_harvest_kmers<true>(enc_state.out_buffer, hit_pos_pool.data(), hits, K, hit_hash_pool.data());
          auto t_harvest_1 = std::chrono::high_resolution_clock::now();
          local_harvest_ns += std::chrono::duration_cast<std::chrono::nanoseconds>(t_harvest_1 - t_harvest_0).count();
          
          // 🛡️ Prevent optimization: Sink the results
          static volatile uint64_t sink;
          sink = hit_hash_pool[hits - 1]; 
        }

        local_hits += hits;
        local_bases += enc_state.total_bases_encoded - carry_offset_adj;
      } else {
        local_bases += enc_state.total_bases_encoded - carry_offset_adj;
      }

      // Mark task as completed in order
      task_done_flags[task.task_id % RING_SIZE].store(
          task.task_id, std::memory_order_release);
    }
  }

  int K, S;
  XsealScannerType scanner_type;
  size_t BATCH_THRESHOLD;
  std::vector<std::thread> workers;
  static constexpr size_t RING_SIZE = 4096;
  struct TaskSlot {
    std::atomic<uint64_t> seq;
    XsealTask<NumSlots> task;
    char pad[64]; // Padding to avoid false sharing
  };

  alignas(64) TaskSlot task_ring[RING_SIZE];
  alignas(64) std::atomic<uint64_t> ring_push_idx{0};
  alignas(64) std::atomic<uint64_t> ring_consume_idx{0};
  alignas(64) std::atomic<uint64_t> task_done_flags[RING_SIZE];

  std::atomic<bool> finished;

  std::vector<uint8_t> saved_tail;
  XsealTask<NumSlots> current_task;

  std::atomic<uint64_t> tasks_pushed_count{0};
  std::atomic<uint64_t> tasks_completed_upto{0};

  std::atomic<size_t> total_hits{0};
  std::atomic<size_t> total_bases{0};
  std::atomic<uint64_t> total_scan_ns{0};
  std::atomic<uint64_t> total_encode_ns{0};
  std::atomic<uint64_t> total_harvest_ns{0};
};

} // namespace xseal
