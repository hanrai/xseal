// SPDX-License-Identifier: Apache-2.0
// Copyright (c) 2026 hanrai. All Rights Reserved.

#pragma once

#include <liburing.h>
#include <fcntl.h>
#include <unistd.h>
#include <sys/stat.h>
#include <sys/mman.h>
#include <stdexcept>
#include <string>
#include <vector>
#include <cstring>
#include <atomic>
#include <immintrin.h>

namespace xseal {

/**
 * XsealUringWriter: High-performance asynchronous writer using io_uring.
 * Supports both sequential file writing and pipe streaming.
 */
template <size_t ChunkSize = 24 * 1024 * 1024>
class XsealUringWriterT {
public:
    static constexpr size_t CHUNK_SIZE = ChunkSize;
    static constexpr uint32_t NUM_SLOTS = 16;

    XsealUringWriterT(int fd, bool is_pipe = true) : fd(fd), is_pipe(is_pipe) {
        // 1. Pipe Optimization
        if (is_pipe) {
#ifdef __linux__
            int current_pipe_size = fcntl(fd, F_GETPIPE_SZ);
            if (current_pipe_size > 0) {
                fcntl(fd, F_SETPIPE_SZ, 1048576); 
            }
#endif
        }

        // 2. Initialize io_uring
        if (io_uring_queue_init(NUM_SLOTS * 2, &ring, 0) < 0) {
            throw std::runtime_error("❌ XsealUringWriter: Failed to initialize io_uring.");
        }

        // 3. Allocate and Register Fixed Buffers
        if (posix_memalign((void**)&pool_base, 4096, NUM_SLOTS * CHUNK_SIZE) != 0) {
            io_uring_queue_exit(&ring);
            throw std::bad_alloc();
        }

        for (uint32_t i = 0; i < NUM_SLOTS; ++i) {
            iovs[i].iov_base = pool_base + i * CHUNK_SIZE;
            iovs[i].iov_len = CHUNK_SIZE;
            slot_active[i] = false;
            slot_completed[i] = true; // Initially all slots are ready to be filled
        }
        io_uring_register_buffers(&ring, iovs, NUM_SLOTS);
    }

    ~XsealUringWriterT() {
        flush();
        io_uring_queue_exit(&ring);
        if (pool_base) free(pool_base);
    }

    /**
     * @brief Get a pointer to the next available slot for writing.
     * Blocks if all slots are currently being written to disk.
     */
    uint8_t* get_next_slot() {
        uint32_t slot_id = next_user_slot;
        if (slot_active[slot_id]) {
            _wait_for_slot(slot_id);
        }
        return (uint8_t*)iovs[slot_id].iov_base;
    }

    /**
     * @brief Copy data to the next available slot using non-temporal streaming stores.
     * This avoids cache pollution and is ideal for data that won't be read back by the CPU.
     * Note: The next_user_slot is automatically committed after copying.
     */
    void write_data(const void* src, size_t len) {
        uint8_t* slot = get_next_slot();
        
        const __m256i* s = (const __m256i*)src;
        __m256i* d = (__m256i*)slot;
        size_t n = len / 32;
        
        for (size_t i = 0; i < n; ++i) {
            _mm256_stream_si256(d + i, _mm256_loadu_si256(s + i));
        }
        
        // Handle tail with standard memcpy
        size_t rem = len % 32;
        if (rem > 0) {
            std::memcpy(slot + n * 32, (const uint8_t*)src + n * 32, rem);
        }
        
        _mm_sfence(); // Memory barrier to ensure streaming stores are visible
        commit_slot(len);
    }

    /**
     * @brief Submit the current slot for asynchronous writing.
     */
    void commit_slot(size_t len) {
        uint32_t slot_id = next_user_slot;
        
        struct io_uring_sqe* sqe = io_uring_get_sqe(&ring);
        if (!sqe) {
            io_uring_submit(&ring);
            sqe = io_uring_get_sqe(&ring);
        }

        int64_t offset = is_pipe ? -1 : current_file_offset;
        io_uring_prep_write_fixed(sqe, fd, iovs[slot_id].iov_base, len, offset, slot_id);
        io_uring_sqe_set_data(sqe, (void*)(uintptr_t)slot_id);
        
        slot_active[slot_id] = true;
        slot_completed[slot_id] = false;
        
        if (!is_pipe) current_file_offset += len;
        
        io_uring_submit(&ring);
        
        next_user_slot = (next_user_slot + 1) % NUM_SLOTS;
        inflight_count++;
    }

    /**
     * @brief Wait for all pending writes to complete.
     */
    void flush() {
        while (inflight_count > 0) {
            struct io_uring_cqe* cqe;
            int ret = io_uring_wait_cqe(&ring, &cqe);
            if (ret < 0) break;
            
            uint32_t slot_id = (uint32_t)(uintptr_t)io_uring_cqe_get_data(cqe);
            slot_active[slot_id] = false;
            slot_completed[slot_id] = true;
            inflight_count--;
            io_uring_cqe_seen(&ring, cqe);
        }
    }

private:
    int fd;
    bool is_pipe;
    struct io_uring ring;
    uint8_t* pool_base = nullptr;
    struct iovec iovs[NUM_SLOTS];
    
    bool slot_active[NUM_SLOTS];
    bool slot_completed[NUM_SLOTS];
    
    uint32_t next_user_slot = 0;
    uint64_t current_file_offset = 0;
    std::atomic<uint32_t> inflight_count{0};

    void _wait_for_slot(uint32_t slot_id) {
        while (!slot_completed[slot_id]) {
            struct io_uring_cqe* cqe;
            int ret = io_uring_wait_cqe(&ring, &cqe);
            if (ret < 0) {
                if (ret == -EINTR) continue;
                throw std::runtime_error("❌ XsealUringWriter: io_uring_wait_cqe failed");
            }
            
            uint32_t comp_slot = (uint32_t)(uintptr_t)io_uring_cqe_get_data(cqe);
            slot_active[comp_slot] = false;
            slot_completed[comp_slot] = true;
            inflight_count--;
            io_uring_cqe_seen(&ring, cqe);
        }
    }
};

using XsealUringWriter = XsealUringWriterT<>;

} // namespace xseal
