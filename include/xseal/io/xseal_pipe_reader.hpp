// SPDX-License-Identifier: Apache-2.0
// Copyright (c) 2026 hanrai. All Rights Reserved.

#pragma once

#include <cstring>
#include <fcntl.h>
#include <immintrin.h>
#include <liburing.h>
#include <stdexcept>
#include <unistd.h>
#include <xseal/xseal_slot_pool.hpp>

namespace xseal {

/**
 * XsealPipeReader: High-performance stream reader for Pipes/FIFOs.
 * Fixes "Short Read" fragmentation and "Interleaving" by using QD=1 and
 * Read-Until-Full.
 */
template <uint32_t NumSlots = 768> class XsealPipeReaderT {
public:
  static constexpr uint32_t NUM_SLOTS = NumSlots;

  explicit XsealPipeReaderT(size_t chunk_bytes = 2 * 1024 * 1024)
      : fd(STDIN_FILENO), chunk_size_(chunk_bytes) {
    if (chunk_size_ == 0)
      throw std::runtime_error("❌ XsealPipeReader: chunk_bytes must be > 0");
    _init();
  }

  explicit XsealPipeReaderT(int input_fd,
                            size_t chunk_bytes = 2 * 1024 * 1024)
      : fd(input_fd), chunk_size_(chunk_bytes) {
    if (chunk_size_ == 0)
      throw std::runtime_error("❌ XsealPipeReader: chunk_bytes must be > 0");
    _init();
  }

  ~XsealPipeReaderT() {
    io_uring_queue_exit(&ring);
    if (pool_base)
      free(pool_base);
  }

  size_t chunk_size() const noexcept { return chunk_size_; }

  XsealSlot<NumSlots> getNext(size_t &out_len, uint64_t &out_offset) {
    if (eof_reached && !is_read_in_flight && head == tail) {
      return XsealSlot<NumSlots>();
    }

    while (true) {
      // 1. Is there a fully loaded 24MB box on the shelf? Dispatch it!
      if (head != tail) {
        uint32_t idx = head % NUM_SLOTS;
        out_len = ready_queue[idx].filled_bytes;
        XsealSlot<NumSlots> result = std::move(ready_queue[idx].slot);

        head++;
        out_offset = total_offset;
        total_offset += out_len;
        return result;
      }

      // 2. Shelf is empty, urge packing at the bottom
      _fill_pipeline();

      // 3. Peek at the kernel: is this scoop of water fully received?
      struct io_uring_cqe *cqe;
      if (io_uring_peek_cqe(&ring, &cqe) == 0) {
        ssize_t res = cqe->res;
        io_uring_cqe_seen(&ring, cqe);
        is_read_in_flight = false; // This request has landed

        if (res <= 0) {
          // EOF: pipe is closed
          eof_reached = true;
          if (active_filled_bytes > 0) {
            // Push the last partially filled box onto the shelf
            uint32_t idx = tail % NUM_SLOTS;
            ready_queue[idx].slot = std::move(active_filling_slot);
            ready_queue[idx].filled_bytes = active_filled_bytes;
            tail++;
          }
        } else {
          // Received water! Accumulate into the current box
          active_filled_bytes += res;

          if (active_filled_bytes == chunk_size_) {
            // Finally filled 24MB! Seal the box and push to shelf!
            uint32_t idx = tail % NUM_SLOTS;
            ready_queue[idx].slot = std::move(active_filling_slot);
            ready_queue[idx].filled_bytes = active_filled_bytes;
            tail++;
          }
        }
        continue; // State changed, restart immediately, don't stop!
      }

      // 4. No water received, no completion events, nap for a few nanoseconds
      if (eof_reached && !is_read_in_flight && head == tail) {
        return XsealSlot<NumSlots>();
      }
      for (int i = 0; i < 16; ++i)
        _mm_pause();
    }
  }

  /// Bytes delivered to the parser (for Raw File Throughput after EOF).
  uint64_t get_file_size() const { return total_offset; }
  const uint8_t *get_pool_base() const { return pool_base; }

private:
  struct ReadySlot {
    XsealSlot<NumSlots> slot;
    size_t filled_bytes = 0;
  };

  int fd;
  const size_t chunk_size_;
  struct io_uring ring;
  uint8_t *pool_base = nullptr;
  XsealSlotPool<NumSlots> slot_pool;

  ReadySlot ready_queue[NUM_SLOTS];
  uint32_t head = 0;
  uint32_t tail = 0;

  // Reservoir state machine (specialized for short reads)
  XsealSlot<NumSlots> active_filling_slot;
  size_t active_filled_bytes = 0;
  bool is_read_in_flight = false;

  uint64_t total_offset = 0;
  bool eof_reached = false;

  void _init() {
    if (fd < 0)
      throw std::runtime_error("❌ XsealPipeReader: Invalid FD");
    if (io_uring_queue_init(16, &ring, 0) < 0) { // Set depth small since we are queuing
      throw std::runtime_error("❌ XsealPipeReader: io_uring init failed");
    }
    if (posix_memalign((void **)&pool_base, 4096, NUM_SLOTS * chunk_size_) !=
        0) {
      throw std::bad_alloc();
    }
    std::memset(pool_base, 0, NUM_SLOTS * chunk_size_);
  }

  // Single-thread precise water lifting
  inline void _fill_pipeline() {
    if (eof_reached || is_read_in_flight)
      return; // Absolutely no concurrency!

    if (!active_filling_slot.valid()) {
      active_filling_slot = slot_pool.acquire();
      if (!active_filling_slot.valid())
        return;                // Out of boxes, creating natural backpressure
      active_filled_bytes = 0; // Got a new box, reset water level
    }

    struct io_uring_sqe *sqe = io_uring_get_sqe(&ring);
    if (!sqe)
      return;

    // Essence: determine write pointer offset and bytes to read based on current water level
    uint8_t *target_buf =
        pool_base + active_filling_slot.id() * chunk_size_ + active_filled_bytes;
    size_t bytes_to_request = chunk_size_ - active_filled_bytes;

    io_uring_prep_read(sqe, fd, target_buf, bytes_to_request, -1);
    io_uring_sqe_set_data(sqe, nullptr); // No SeqID needed, only one in flight

    is_read_in_flight = true;
    io_uring_submit(&ring);
  }
};

using XsealPipeReader = XsealPipeReaderT<>;

} // namespace xseal
