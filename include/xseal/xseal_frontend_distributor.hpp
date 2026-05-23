// SPDX-License-Identifier: Apache-2.0
// Copyright (c) 2026 hanrai. All Rights Reserved.

#pragma once

#include <atomic>
#include <cstring>
#include <immintrin.h>
#include <thread>
#include <vector>
#include <functional>
#include <memory>
#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>
#include <chrono>
#include <xseal/xseal_encoder.hpp>
#include <xseal/xseal_fasta_parser.hpp>
#include <xseal/xseal_minimizer.hpp>
#include <xseal/xseal_syncmer.hpp>

namespace xseal {

struct XsealHit {
    uint32_t hash;
    uint32_t pos;
    uint32_t rid;
    uint32_t span;
    uint32_t strand;
};

struct XsealSegmentRef {
    size_t buffer_offset;
    size_t len;
    uint32_t rid;
};

template <uint32_t NumSlots = 768> struct XsealFrontendTask {
    std::vector<uint8_t> ascii_buffer; 
    std::vector<XsealSegmentRef> segments;
    std::vector<std::pair<size_t, uint32_t>> rid_boundaries; // end_pos_in_out_buf -> rid
    uint64_t task_id = 0;
};

template <uint32_t NumSlots = 768> class XsealFrontendDistributor {
public:
    XsealFrontendDistributor(int num_threads, int k, int w) : k_val(k), w_val(w), num_workers(num_threads) {
        finished.store(false, std::memory_order_relaxed);
        ring_push_idx.store(0);
        ring_consume_idx.store(0);
        tasks_pushed_count.store(0);
        tasks_completed_count.store(0);
        total_bases_processed.store(0);
        total_hits.store(0);
        current_rid.store(0xFFFFFFFF);
        start_time = std::chrono::high_resolution_clock::now();

        for (size_t i = 0; i < RING_SIZE; ++i) {
            task_ring[i].seq.store(i, std::memory_order_relaxed);
        }
        for (int i = 0; i < num_threads; ++i) {
            workers.emplace_back(&XsealFrontendDistributor::worker_loop, this);
        }
        current_task.ascii_buffer.reserve(16 * 1024 * 1024 + 256);
    }

    ~XsealFrontendDistributor() { stop(); }

    void stop() {
        if (finished.load()) return;
        flush_current_batch();
        finished.store(true, std::memory_order_release);
        for (auto &t : workers) {
            if (t.joinable()) t.join();
        }
        
        auto end_time = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> diff = end_time - start_time;
        
        fprintf(stderr, "\n[XSeal FINAL BENCHMARK]\n");
        fprintf(stderr, "File: %s\n", last_fn.c_str());
        fprintf(stderr, "Total Hits: %zu\n", total_hits.load());
        fprintf(stderr, "Total Bases: %zu\n", total_bases_processed.load());
        fprintf(stderr, "End-to-End Time: %.3fs\n", diff.count());
        if (diff.count() > 0) {
            fprintf(stderr, "Throughput: %.3f MiB/s\n", (total_bases_processed.load() / (1024.0 * 1024.0)) / diff.count());
        }
        fprintf(stderr, "------------------------\n");
    }

    void process_file(const char* fn, std::function<void(const std::vector<XsealHit>&)> callback) {
        hit_callback = callback;
        last_fn = fn;
        int fd = open(fn, O_RDONLY);
        if (fd < 0) return;
        struct stat st;
        fstat(fd, &st);
        size_t size = st.st_size;
        uint8_t* data = (uint8_t*)mmap(NULL, size, PROT_READ, MAP_SHARED, fd, 0);
        if (data == MAP_FAILED) { close(fd); return; }

        XsealFastaParser parser;
        XsealRunTable run_table;
        run_table.reserve(4096);

        size_t chunk_size = 64 * 1024 * 1024;
        for (size_t offset = 0; offset < size; offset += chunk_size) {
            size_t current_chunk = std::min(chunk_size, size - offset);
            run_table.clear();
            parser.parse(data + offset, current_chunk, offset, run_table);
            
            for (const auto &seg : run_table) {
                if (seg.type == SegType::HEADER) {
                    flush_current_batch();
                    current_rid.fetch_add(1);
                } else if (seg.type == SegType::BASES) {
                    size_t buf_off = current_task.ascii_buffer.size();
                    current_task.ascii_buffer.insert(current_task.ascii_buffer.end(), data + seg.file_offset, data + seg.file_offset + seg.byte_len);
                    
                    XsealSegmentRef ref;
                    ref.buffer_offset = buf_off;
                    ref.len = seg.byte_len;
                    ref.rid = current_rid.load();
                    current_task.segments.push_back(ref);
                    
                    if (current_task.ascii_buffer.size() >= 12 * 1024 * 1024) flush_current_batch();
                }
            }
        }
        
        munmap(data, size);
        close(fd);
        stop();
    }

protected:
    void flush_current_batch() {
        if (current_task.segments.empty()) return;
        current_task.ascii_buffer.insert(current_task.ascii_buffer.end(), 128, 'N');

        uint64_t idx = ring_push_idx.load(std::memory_order_relaxed);
        while (task_ring[idx % RING_SIZE].seq.load(std::memory_order_acquire) != idx) {
            _mm_pause();
        }
        task_ring[idx % RING_SIZE].task = std::move(current_task);
        task_ring[idx % RING_SIZE].seq.store(idx + 1, std::memory_order_release);
        ring_push_idx.store(idx + 1, std::memory_order_relaxed);
        tasks_pushed_count.fetch_add(1, std::memory_order_release);
        
        current_task = XsealFrontendTask<NumSlots>();
        current_task.ascii_buffer.reserve(16 * 1024 * 1024 + 256);
    }

    void worker_loop() {
        XSealEncoder encoder;
        XSealEncoderState enc_state;
        XSealMinimizer scanner(k_val, w_val);
        
        size_t buffer_cap = 16 * 1024 * 1024 + 1024;
        void* enc_ptr = nullptr;
        posix_memalign(&enc_ptr, 32, buffer_cap);
        std::unique_ptr<uint8_t, void(*)(void*)> out_buf((uint8_t*)enc_ptr, free);
        
        size_t hit_cap = 16 * 1024 * 1024;
        void* hit_ptr = nullptr;
        posix_memalign(&hit_ptr, 32, hit_cap * sizeof(uint32_t));
        std::unique_ptr<uint32_t, void(*)(void*)> hit_pos((uint32_t*)hit_ptr, free);
        
        std::vector<XsealHit> local_hits;
        local_hits.reserve(1024 * 1024);

        while (true) {
            uint64_t idx = ring_consume_idx.fetch_add(1, std::memory_order_relaxed);
            if (idx >= tasks_pushed_count.load(std::memory_order_acquire) && finished.load()) return;
            
            while (task_ring[idx % RING_SIZE].seq.load(std::memory_order_acquire) != idx + 1) {
                if (finished.load() && idx >= tasks_pushed_count.load()) return;
                _mm_pause();
            }
            XsealFrontendTask<NumSlots> task = std::move(task_ring[idx % RING_SIZE].task);
            task_ring[idx % RING_SIZE].seq.store(idx + RING_SIZE, std::memory_order_release);

            local_hits.clear();
            memset(&enc_state, 0, sizeof(enc_state));
            encoder.init_state(&enc_state, out_buf.get(), buffer_cap);
            
            // RID mapping setup
            struct RidMap { uint32_t end_pos; uint32_t rid; };
            std::vector<RidMap> rid_map;
            rid_map.reserve(task.segments.size());

            const char *p_start = (const char *)task.ascii_buffer.data();
            const char *p = p_start;
            const char *end = p + task.ascii_buffer.size();
            
            // Encode segments and build RID map
            for (const auto& seg : task.segments) {
                const char* s_p = (const char*)task.ascii_buffer.data() + seg.buffer_offset;
                const char* s_end = s_p + seg.len;
                encoder.encode_chunk(&enc_state, &s_p, s_end);
                rid_map.push_back({(uint32_t)enc_state.total_bases_encoded, seg.rid});
            }
            encoder.flush_state(&enc_state);

            size_t num_chunks = enc_state.total_bases_encoded / 128;
            if (num_chunks > 1) {
                size_t n = 0;
                if (w_val == 10) n = scanner.scan<10, true>((const __m256i*)out_buf.get(), num_chunks, hit_pos.get());
                else if (w_val == 20) n = scanner.scan<20, true>((const __m256i*)out_buf.get(), num_chunks, hit_pos.get());
                else if (w_val == 50) n = scanner.scan<50, true>((const __m256i*)out_buf.get(), num_chunks, hit_pos.get());
                else if (w_val == 100) n = scanner.scan<100, true>((const __m256i*)out_buf.get(), num_chunks, hit_pos.get());
                else n = scanner.scan<10, true>((const __m256i*)out_buf.get(), num_chunks, hit_pos.get());
                
                size_t rid_idx = 0;
                for (size_t i = 0; i < n; ++i) {
                    uint32_t pos = hit_pos.get()[i];
                    while (rid_idx < rid_map.size() && pos >= rid_map[rid_idx].end_pos) rid_idx++;
                    if (rid_idx >= rid_map.size()) break;

                    XsealHit h;
                    h.hash = hit_pos.get()[i]; // Note: XSealMinimizer scan doesn't output hashes, this is a placeholder
                    h.pos = pos;
                    h.rid = rid_map[rid_idx].rid;
                    h.span = (uint32_t)k_val;
                    h.strand = 0; 
                    local_hits.push_back(h);
                }
            }
            
            if (hit_callback && !local_hits.empty()) {
                total_hits.fetch_add(local_hits.size());
                hit_callback(local_hits);
            }
            total_bases_processed.fetch_add(enc_state.total_bases_encoded);
            tasks_completed_count.fetch_add(1, std::memory_order_release);
        }
    }

    static constexpr size_t RING_SIZE = 1024;
    struct TaskSlot {
        std::atomic<uint64_t> seq;
        XsealFrontendTask<NumSlots> task;
    };

    TaskSlot task_ring[RING_SIZE];
    std::atomic<uint64_t> ring_push_idx;
    std::atomic<uint64_t> ring_consume_idx;
    std::atomic<bool> finished;
    XsealFrontendTask<NumSlots> current_task;
    std::atomic<uint64_t> tasks_pushed_count;
    std::atomic<uint64_t> tasks_completed_count;
    std::atomic<uint64_t> total_bases_processed;
    std::atomic<uint32_t> current_rid;
    std::atomic<size_t> total_hits;
    std::function<void(const std::vector<XsealHit>&)> hit_callback;
    std::vector<std::thread> workers;
    std::chrono::time_point<std::chrono::high_resolution_clock> start_time;
    std::string last_fn;
    int num_workers;
    int k_val, w_val;
};

} // namespace xseal
