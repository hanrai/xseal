#pragma once

#include <cstdint>
#include <vector>

namespace xseal {

struct XsealHit;

/** Per-worker TLS bucket buffers; flush once per fused task in worker_loop. */
namespace index_insert {

void push_hits(const std::vector<XsealHit> &hits, uint64_t kmer_mask);
void flush();
void note_task_timing(uint64_t t0_ns, uint64_t t1_ns);
void flush_all_registered(struct mm_idx_s *mi);
void cleanup_all_registered();

} // namespace index_insert
} // namespace xseal
