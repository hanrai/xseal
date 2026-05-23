#pragma once

#include <atomic>
#include <cstring>
#include <cstdio>
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
#include <xseal/xseal_kmer_harvester.hpp>
#include <xseal/xseal_index_insert_tls.hpp>
#include <xseal/xseal_scanner_type.hpp>

namespace xseal {

/** mmap: map FASTA from filesystem; preload: read entire file into RAM first. */
enum class XsealIoMode { MMAP, PRELOAD };

struct XsealHit {
    uint64_t hash;
    uint32_t pos;
    uint32_t rid;
    uint32_t span;
    uint32_t strand;
};

struct XsealSegmentRef {
    const uint8_t *ptr;
    size_t len;
    uint32_t rid;
};

template <uint32_t NumSlots = 768> struct XsealFrontendTask {
    std::vector<XsealSegmentRef> segments;
    size_t total_len = 0;
    uint64_t task_id = 0;
};

template <uint32_t NumSlots = 768> class XsealFrontendDistributor {
public:
    XsealFrontendDistributor(int num_threads, int k, int w, XsealScannerType scanner_type,
                             int s_syncmer)
        : num_workers(num_threads), k_val(k), w_val(w), scanner_type_(scanner_type),
          s_syncmer_(s_syncmer) {
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
    }

    XsealFrontendDistributor(int num_threads, int k, int w)
        : XsealFrontendDistributor(num_threads, k, w, XsealScannerType::MINIMIZER, 0) {}

    void set_io_mode(XsealIoMode mode) { io_mode_ = mode; }
    void set_drop_page_cache(bool drop) { drop_page_cache_ = drop; }
    size_t bases_processed() const { return total_bases_processed.load(); }
    size_t hits_found() const { return total_hits.load(); }

    double file_map_seconds() const
    {
        return static_cast<double>(file_map_ns.load(std::memory_order_relaxed)) / 1e9;
    }
    double parse_seconds() const
    {
        return static_cast<double>(parse_ns.load(std::memory_order_relaxed)) / 1e9;
    }
    double compute_wall_seconds() const
    {
        return static_cast<double>(compute_wall_ns.load(std::memory_order_relaxed)) / 1e9;
    }
    double worker_compute_seconds() const
    {
        return static_cast<double>(worker_compute_ns.load(std::memory_order_relaxed)) / 1e9;
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

        fprintf(stderr, "\n[XSeal HARDENED BENCHMARK]\n");
        fprintf(stderr, "Scanner: %s\n",
                scanner_type_ == XsealScannerType::SYNCMER ? "syncmer" : "minimizer");
        fprintf(stderr, "Hash: canonical k-mer -> minimap2 hash64 (index keys)\n");
        fprintf(stderr, "IO mode: %s\n",
                io_mode_ == XsealIoMode::PRELOAD ? "preload (memory-resident)" : "mmap");
        fprintf(stderr, "File: %s\n", last_fn.c_str());
        fprintf(stderr, "Total Hits: %zu\n", total_hits.load());
        fprintf(stderr, "Total Bases: %zu\n", total_bases_processed.load());
        fprintf(stderr, "End-to-End Time: %.3fs\n", diff.count());
        if (diff.count() > 0) {
            fprintf(stderr, "Throughput: %.3f MiB/s\n",
                    (total_bases_processed.load() / (1024.0 * 1024.0)) / diff.count());
        }
        fprintf(stderr, "------------------------\n");
    }

    void process_file(const char *fn, std::function<void(const std::vector<XsealHit> &)> callback) {
        hit_callback = callback;
        last_fn = fn;
        file_map_ns.store(0, std::memory_order_relaxed);
        parse_ns.store(0, std::memory_order_relaxed);
        compute_wall_ns.store(0, std::memory_order_relaxed);
        worker_compute_ns.store(0, std::memory_order_relaxed);

        using clock = std::chrono::steady_clock;
        const auto t_file0 = clock::now();

        std::vector<uint8_t> owned;
        const uint8_t *data = nullptr;
        size_t size = 0;
        int fd = -1;
        uint8_t *mapped = nullptr;

        fd = open(fn, O_RDONLY);
        if (fd < 0)
            return;
        struct stat st;
        if (fstat(fd, &st) != 0 || st.st_size == 0) {
            close(fd);
            return;
        }
        size = static_cast<size_t>(st.st_size);

        if (io_mode_ == XsealIoMode::PRELOAD) {
            owned.resize(size);
            ssize_t off = 0;
            while (off < static_cast<ssize_t>(size)) {
                ssize_t n = read(fd, owned.data() + off, size - static_cast<size_t>(off));
                if (n <= 0)
                    break;
                off += n;
            }
            close(fd);
            fd = -1;
            if (static_cast<size_t>(off) != size)
                return;
            data = owned.data();
        } else {
            if (drop_page_cache_) {
#if defined(__linux__)
                posix_fadvise(fd, 0, static_cast<off_t>(size), POSIX_FADV_DONTNEED);
#endif
            }
            mapped = (uint8_t *)mmap(NULL, size, PROT_READ, MAP_SHARED, fd, 0);
            close(fd);
            fd = -1;
            if (mapped == MAP_FAILED)
                return;
            data = mapped;
        }

        file_map_ns.store(
            std::chrono::duration_cast<std::chrono::nanoseconds>(clock::now() - t_file0).count(),
            std::memory_order_relaxed);

        const auto t_parse0 = clock::now();
        ingest_buffer(data, size);
        parse_ns.store(
            std::chrono::duration_cast<std::chrono::nanoseconds>(clock::now() - t_parse0).count(),
            std::memory_order_relaxed);

        const auto t_compute0 = clock::now();
        stop();
        compute_wall_ns.store(
            std::chrono::duration_cast<std::chrono::nanoseconds>(clock::now() - t_compute0).count(),
            std::memory_order_relaxed);

        if (mapped)
            munmap(mapped, size);
    }

    void ingest_buffer(const uint8_t *data, size_t size) {
        XsealFastaParser parser;
        XsealRunTable run_table;
        run_table.reserve(4096);

        const size_t chunk_size = 64 * 1024 * 1024;
        for (size_t offset = 0; offset < size; offset += chunk_size) {
            const size_t current_chunk = std::min(chunk_size, size - offset);
            run_table.clear();
            parser.parse(data + offset, current_chunk, offset, run_table);

            for (const auto &seg : run_table) {
                if (seg.type == SegType::HEADER) {
                    flush_current_batch();
                    current_rid.fetch_add(1);
                } else if (seg.type == SegType::BASES) {
                    if (seg.file_offset + seg.byte_len > size)
                        continue;

                    XsealSegmentRef ref;
                    ref.ptr = data + seg.file_offset;
                    ref.len = seg.byte_len;
                    ref.rid = current_rid.load();
                    current_task.segments.push_back(ref);
                    current_task.total_len += seg.byte_len;

                    if (current_task.segments.size() >= NumSlots ||
                        current_task.total_len >= 12 * 1024 * 1024) {
                        flush_current_batch();
                    }
                }
            }
        }
    }

protected:
    void flush_current_batch() {
        if (current_task.segments.empty()) return;

        uint64_t idx = ring_push_idx.load(std::memory_order_relaxed);
        while (task_ring[idx % RING_SIZE].seq.load(std::memory_order_acquire) != idx) {
            _mm_pause();
        }
        task_ring[idx % RING_SIZE].task = std::move(current_task);
        task_ring[idx % RING_SIZE].seq.store(idx + 1, std::memory_order_release);
        ring_push_idx.store(idx + 1, std::memory_order_relaxed);
        tasks_pushed_count.fetch_add(1, std::memory_order_release);

        current_task = XsealFrontendTask<NumSlots>();
    }

    void worker_loop() {
        XSealEncoder encoder;
        XSealEncoderState enc_state;
        XSealMinimizer minimizer_scanner(k_val, w_val);
        XSealSyncmer syncmer_scanner(k_val, s_syncmer_);

        size_t buffer_cap = 32 * 1024 * 1024;
        void *enc_ptr = nullptr;
        posix_memalign(&enc_ptr, 64, buffer_cap);
        std::unique_ptr<uint8_t, void (*)(void *)> out_buf((uint8_t *)enc_ptr, free);

        size_t hit_cap = 16 * 1024 * 1024;
        void *hit_ptr = nullptr;
        posix_memalign(&hit_ptr, 64, hit_cap * sizeof(uint32_t));
        std::unique_ptr<uint32_t, void (*)(void *)> hit_pos((uint32_t *)hit_ptr, free);

        void *hash_ptr = nullptr;
        posix_memalign(&hash_ptr, 64, hit_cap * sizeof(uint64_t));
        std::unique_ptr<uint64_t, void (*)(void *)> hit_hash64((uint64_t *)hash_ptr, free);

        std::vector<XsealHit> local_hits;
        local_hits.reserve(1024 * 1024);

        const int eff_w = w_val;

        while (true) {
            uint64_t idx = ring_consume_idx.fetch_add(1, std::memory_order_relaxed);
            if (idx >= tasks_pushed_count.load(std::memory_order_acquire) && finished.load())
                return;

            while (task_ring[idx % RING_SIZE].seq.load(std::memory_order_acquire) != idx + 1) {
                if (finished.load() && idx >= tasks_pushed_count.load()) return;
                _mm_pause();
            }
            XsealFrontendTask<NumSlots> task = std::move(task_ring[idx % RING_SIZE].task);
            task_ring[idx % RING_SIZE].seq.store(idx + RING_SIZE, std::memory_order_release);

            const auto t_task0 = std::chrono::steady_clock::now();

            local_hits.clear();
            memset(&enc_state, 0, sizeof(enc_state));
            encoder.init_state(&enc_state, out_buf.get(), buffer_cap - 2048);

            struct RidMap {
                uint32_t end_pos;
                uint32_t rid;
            };
            std::vector<RidMap> rid_map;
            rid_map.reserve(task.segments.size());

            for (const auto &seg : task.segments) {
                if (seg.ptr == nullptr || seg.len == 0) continue;

                const char *p = (const char *)seg.ptr;
                const char *end = p + seg.len;
                encoder.encode_chunk(&enc_state, &p, end);
                rid_map.push_back({(uint32_t)enc_state.total_bases_encoded, seg.rid});
            }
            encoder.flush_state(&enc_state);

            size_t num_chunks = enc_state.total_bases_encoded / 128;
            if (num_chunks > 0) {
                size_t n = 0;
                const __m256i *buf = (const __m256i *)out_buf.get();
                const uint8_t *enc_bytes = out_buf.get();

                if (scanner_type_ == XsealScannerType::MINIMIZER) {
                    if (eff_w == 11)
                        n = minimizer_scanner.scan<11, true>(buf, num_chunks, hit_pos.get(), nullptr);
                    else if (eff_w == 10)
                        n = minimizer_scanner.scan<10, true>(buf, num_chunks, hit_pos.get(), nullptr);
                    else if (eff_w == 20)
                        n = minimizer_scanner.scan<20, true>(buf, num_chunks, hit_pos.get(), nullptr);
                    else if (eff_w == 21)
                        n = minimizer_scanner.scan<21, true>(buf, num_chunks, hit_pos.get(), nullptr);
                    else if (eff_w == 50)
                        n = minimizer_scanner.scan<50, true>(buf, num_chunks, hit_pos.get(), nullptr);
                    else if (eff_w == 100)
                        n = minimizer_scanner.scan<100, true>(buf, num_chunks, hit_pos.get(), nullptr);
                    else {
                        fprintf(stderr, "[XSeal] unsupported minimizer w=%d\n", eff_w);
                        n = 0;
                    }
                } else {
                    if (enc_state.total_bases_encoded < (size_t)k_val) {
                        /* skip scan */
                    } else if (eff_w == 11)
                        n = syncmer_scanner.scan<11, 0, true, true>(
                            buf, num_chunks, hit_pos.get(), enc_state.total_bases_encoded);
                    else if (eff_w == 10)
                        n = syncmer_scanner.scan<10, 0, true, true>(
                            buf, num_chunks, hit_pos.get(), enc_state.total_bases_encoded);
                    else if (eff_w == 20)
                        n = syncmer_scanner.scan<20, 0, true, true>(
                            buf, num_chunks, hit_pos.get(), enc_state.total_bases_encoded);
                    else if (eff_w == 21)
                        n = syncmer_scanner.scan<21, 0, true, true>(
                            buf, num_chunks, hit_pos.get(), enc_state.total_bases_encoded);
                    else if (eff_w == 50)
                        n = syncmer_scanner.scan<50, 0, true, true>(
                            buf, num_chunks, hit_pos.get(), enc_state.total_bases_encoded);
                    else if (eff_w == 100)
                        n = syncmer_scanner.scan<100, 0, true, true>(
                            buf, num_chunks, hit_pos.get(), enc_state.total_bases_encoded);
                    else {
                        fprintf(stderr, "[XSeal] unsupported syncmer w=%d\n", eff_w);
                        n = 0;
                    }
                }

                if (n > 0) {
                    xseal_harvest_kmers<true, false>(enc_bytes, hit_pos.get(), n, k_val,
                                                     hit_hash64.get(), 16);
                    size_t rid_idx = 0;
                    for (size_t i = 0; i < n; ++i) {
                        uint32_t pos = hit_pos.get()[i];
                        while (rid_idx < rid_map.size() && pos >= rid_map[rid_idx].end_pos) rid_idx++;
                        if (rid_idx >= rid_map.size()) break;

                        XsealHit h;
                        h.hash = hit_hash64.get()[i];
                        h.pos = pos;
                        h.rid = rid_map[rid_idx].rid;
                        h.span = (uint32_t)k_val;
                        h.strand = 0;
                        local_hits.push_back(h);
                    }
                }
            }

            if (hit_callback && !local_hits.empty()) {
                total_hits.fetch_add(local_hits.size());
                hit_callback(local_hits);
            }
            total_bases_processed.fetch_add(enc_state.total_bases_encoded);
            tasks_completed_count.fetch_add(1, std::memory_order_release);

            worker_compute_ns.fetch_add(
                static_cast<uint64_t>(std::chrono::duration_cast<std::chrono::nanoseconds>(
                                          std::chrono::steady_clock::now() - t_task0)
                                          .count()),
                std::memory_order_relaxed);
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
    std::atomic<size_t> total_bases_processed;
    std::atomic<uint32_t> current_rid;
    std::atomic<size_t> total_hits;
    std::atomic<uint64_t> file_map_ns{0};
    std::atomic<uint64_t> parse_ns{0};
    std::atomic<uint64_t> compute_wall_ns{0};
    std::atomic<uint64_t> worker_compute_ns{0};
    std::function<void(const std::vector<XsealHit> &)> hit_callback;
    std::vector<std::thread> workers;
    std::chrono::time_point<std::chrono::high_resolution_clock> start_time;
    std::string last_fn;
    int num_workers;
    int k_val, w_val;
    XsealScannerType scanner_type_;
    int s_syncmer_;
    XsealIoMode io_mode_{XsealIoMode::MMAP};
    bool drop_page_cache_{false};
};

} // namespace xseal
