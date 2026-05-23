// SPDX-License-Identifier: Apache-2.0
// Copyright (c) 2026 hanrai. All Rights Reserved.

#pragma once

#include <algorithm>
#include <cstdint>
#include <cstring>
#include <fcntl.h>
#include <immintrin.h> // for _mm_pause()
#include <liburing.h>
#include <stdexcept>
#include <string>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>
#include <xseal/xseal_slot_pool.hpp>

namespace xseal {

/**
 * XsealUringReader: High-performance io_uring reader for files/NVMe.
 * V4.0 Architecture: Zero-blocking, SeqID ROB, Direct I/O, Registered Buffers.
 */
template <uint32_t NumSlots = 768> class XsealUringReaderT {
public:
  static constexpr uint32_t NUM_SLOTS = NumSlots;

  explicit XsealUringReaderT(const std::string &filepath,
                             size_t chunk_bytes = 2 * 1024 * 1024)
      : fd(-1), file_size(0), use_direct_io(true), chunk_size_(0) {

    if (chunk_bytes == 0)
      throw std::runtime_error("❌ XsealUringReader: chunk_bytes must be > 0");

    // Attempt to open file with O_DIRECT
    fd = open(filepath.c_str(), O_RDONLY | O_DIRECT);
    if (fd < 0) {
      fd = open(filepath.c_str(), O_RDONLY);
      use_direct_io = false;
    }

    struct stat st;
    if (fstat(fd, &st) == 0) {
      file_size = st.st_size;
    } else {
      throw std::runtime_error(
          "❌ XsealUringReader: Failed to stat file size.");
    }

    size_t cs = chunk_bytes;
    if (use_direct_io) {
      cs = std::max(cs, size_t(4096));
      cs = (cs + 4095) & ~size_t(4095);
    }
    chunk_size_ = cs;

    // Initialize ROB, queue depth must accommodate the concurrency limit
    if (io_uring_queue_init(ROB_SIZE, &ring, 0) < 0) {
      if (fd >= 0)
        close(fd);
      throw std::runtime_error(
          "❌ XsealUringReader: Failed to initialize io_uring.");
    }

    // 4K aligned allocation required by O_DIRECT
    if (posix_memalign((void **)&pool_base, 4096,
                       NUM_SLOTS * chunk_size_) != 0) {
      io_uring_queue_exit(&ring);
      if (fd >= 0)
        close(fd);
      throw std::bad_alloc();
    }
    std::memset(pool_base, 0, NUM_SLOTS * chunk_size_);

    // Core acceleration: register buffers and FDs, the key to squeezing NVMe
    struct iovec iovs[NUM_SLOTS];
    for (uint32_t i = 0; i < NUM_SLOTS; ++i) {
      iovs[i].iov_base = pool_base + i * chunk_size_;
      iovs[i].iov_len = chunk_size_;
    }
    io_uring_register_buffers(&ring, iovs, NUM_SLOTS);
    io_uring_register_files(&ring, &fd, 1);

    _fill_pipeline(); // Warm up pipeline
  }

  ~XsealUringReaderT() {
    if (fd >= 0)
      close(fd);
    io_uring_queue_exit(&ring);
    if (pool_base)
      free(pool_base);
  }

  uint64_t get_file_size() const { return file_size; }
  const uint8_t *get_pool_base() const { return pool_base; }
  size_t chunk_size() const noexcept { return chunk_size_; }

  /**
   * getNext: V4.0 Pull-based API. Zero OS blocking, purely hardware spinning.
   */
  XsealSlot<NumSlots> getNext(size_t &out_len, uint64_t &out_offset) {
    if (processed_offset >= file_size)
      return XsealSlot<NumSlots>(); // Completely read

    while (true) {
      // 1. Attempt to furiously fill the underlying NVMe queue
      _fill_pipeline();

      // 2. Check ROB: is the expected SeqID block completed?
      uint32_t expected_idx = expected_seq & ROB_MASK;
      if (rob[expected_idx].ready) {
        ssize_t res = rob[expected_idx].res;

        // 🛡️ Defense 1: Handle kernel-level resource exhaustion and errors
        if (res < 0) {
          if (res == -EAGAIN || res == -EINTR) {
            // Kernel saturated by 7GB/s burst concurrency, retrying requested
            rob[expected_idx].ready = false;
            _retry_submit(expected_idx, expected_seq);
            continue; // Kick back to reread, continue polling
          } else {
            throw std::runtime_error("❌ io_uring fatal error: " +
                                     std::to_string(res));
          }
        }

        // 🛡️ Defense 2: Strictly guard against silent corruption from short reads
        if (res >= 0 && (size_t)res < rob[expected_idx].actual_len) {
          throw std::runtime_error(
              "❌ Unhandled Short Read! Expected: " +
              std::to_string(rob[expected_idx].actual_len) +
              ", Got: " + std::to_string(res));
        }

        // Grab the goods!
        out_len = rob[expected_idx].actual_len;
        out_offset = rob[expected_idx].offset;
        XsealSlot<NumSlots> result = std::move(rob[expected_idx].slot);

        // Advance state
        rob[expected_idx].ready = false;
        expected_seq++;
        active_kernel_reads--;
        processed_offset += out_len;

        return result;
      }

      // 3. Out of stock? Pure user-space peek at CQE
      struct io_uring_cqe *cqe;
      unsigned head;
      unsigned peeked_count = 0;

      io_uring_for_each_cqe(&ring, head, cqe) {
        // Extract the SeqID we hid in SQE
        uint64_t seq = (uint64_t)io_uring_cqe_get_data(cqe);
        uint32_t idx = seq & ROB_MASK;

        // 🎯 Core fix: Record the true kernel receipt
        rob[idx].res = cqe->res;
        rob[idx].ready = true;
        peeked_count++;
      }

      // 4. If data was peeked, bulk clear CQE and retry immediately!
      if (peeked_count > 0) {
        io_uring_cq_advance(&ring, peeked_count);
        continue;
      }

      // 5. Really out of stock, nanosecond-level spin-yield, never context yield!
      if (processed_offset >= file_size) {
        return XsealSlot<NumSlots>();
      }
      for (int i = 0; i < 16; ++i)
        _mm_pause();
    }
  }

private:
  // ROB size must be a power of 2 and larger than NUM_SLOTS to form an anti-aliasing barrier
  static constexpr uint32_t ROB_SIZE = 1024;
  static constexpr uint32_t ROB_MASK = ROB_SIZE - 1;

  struct RobEntry {
    XsealSlot<NumSlots> slot;
    size_t actual_len = 0;
    uint64_t offset = 0;
    ssize_t res = 0; // Added: save kernel return result
    bool ready = false;
  };

  int fd;
  struct io_uring ring;
  uint8_t *pool_base = nullptr;
  uint64_t file_size;
  bool use_direct_io;
  size_t chunk_size_;
  XsealSlotPool<NumSlots> slot_pool;

  RobEntry rob[ROB_SIZE];

  uint64_t expected_seq = 0;
  uint64_t next_submit_seq = 0;
  uint64_t submitted_offset = 0;
  uint64_t processed_offset = 0;
  uint32_t active_kernel_reads = 0;

  // Core filling engine
  inline void _fill_pipeline() {
    if (submitted_offset >= file_size)
      return;

    bool submitted_any = false;
    while (active_kernel_reads < NUM_SLOTS && submitted_offset < file_size) {
      XsealSlot<NumSlots> slot = slot_pool.acquire();
      if (!slot.valid())
        break; // Worker not done, pool empty, trigger natural backpressure

      uint64_t seq = next_submit_seq++;
      uint32_t idx = seq & ROB_MASK;

      // Calculate the true length of this read
      uint32_t actual_len = (submitted_offset + chunk_size_ > file_size)
                                ? (file_size - submitted_offset)
                                : static_cast<uint32_t>(chunk_size_);

      rob[idx].slot = std::move(slot);
      rob[idx].actual_len = actual_len;
      rob[idx].offset = submitted_offset;
      rob[idx].ready = false;

      // Calculate aligned request length for O_DIRECT
      uint32_t req_len = actual_len;
      if (use_direct_io) {
        req_len = (actual_len + 4095) & ~4095;
      }

      struct io_uring_sqe *sqe = io_uring_get_sqe(&ring);
      if (!sqe) {
        // I/O queue full, rollback state
        rob[idx].slot.reset();
        next_submit_seq--;
        break;
      }

      // Top performance esoteric: use FIXED_FILE and registered memory buffer index (slot.id)
      uint32_t buf_index = rob[idx].slot.id();
      io_uring_prep_read_fixed(sqe, 0, pool_base + buf_index * chunk_size_,
                               req_len, submitted_offset, buf_index);
      sqe->flags |= IOSQE_FIXED_FILE;

      // Stuff SeqID as an amulet into the kernel
      io_uring_sqe_set_data(sqe, (void *)(uintptr_t)seq);

      submitted_offset += actual_len;
      active_kernel_reads++;
      submitted_any = true;
    }

    if (submitted_any) {
      io_uring_submit(&ring);
    }
  }

  // 🛠️ Retry engine specifically designed to combat -EAGAIN
  inline void _retry_submit(uint32_t idx, uint64_t seq) {
    struct io_uring_sqe *sqe = io_uring_get_sqe(&ring);
    if (!sqe) {
      io_uring_submit(&ring);
      sqe = io_uring_get_sqe(&ring);
    }

    uint32_t buf_index = rob[idx].slot.id();
    uint32_t req_len = rob[idx].actual_len;
    if (use_direct_io) {
      req_len = (req_len + 4095) & ~4095;
    }

    io_uring_prep_read_fixed(sqe, 0, pool_base + buf_index * chunk_size_,
                             req_len, rob[idx].offset, buf_index);
    sqe->flags |= IOSQE_FIXED_FILE;
    io_uring_sqe_set_data(sqe, (void *)(uintptr_t)seq);

    io_uring_submit(&ring); // Immediately kick back into kernel
  }
};

using XsealUringReader = XsealUringReaderT<>;

} // namespace xseal
