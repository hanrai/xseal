#pragma once

#include <cstdint>

namespace xseal {

enum class SegType : uint8_t { HEADER = 0, BASES = 1, N_GAP = 2, SKIP = 3 };

struct SeqSegment {
    uint64_t file_offset;      
    uint32_t logical_len; 
    uint32_t byte_len; 
    SegType  type; 
    bool     is_physically_truncated; 
};

/**
 * @brief Standard hit structure for position-only output.
 */
using SyncmerHit = uint64_t;

/**
 * @brief Full hit structure including position and 32-bit hash.
 * Total size: 8 bytes (32-bit pos + 32-bit hash).
 * This fits perfectly into a 64-bit register for vectorized collection.
 */
struct SyncmerHitFull {
  uint64_t pos;
  uint64_t hash;

  bool operator==(const SyncmerHitFull &other) const {
    return pos == other.pos && hash == other.hash;
  }
};

} // namespace xseal
