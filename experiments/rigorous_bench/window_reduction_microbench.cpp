// SPDX-License-Identifier: Apache-2.0
// Copyright (c) 2026 hanrai. All Rights Reserved.

/**
 * window_reduce_bench — reproducible fold vs sparse-table window minimum (SIMD)
 *
 * Mirrors the Phase-2 pattern in include/xseal/xseal_syncmer.hpp (W=21, batches
 * 0..127 step 8): for each batch, reduce W=21 lane-wise minima over __m256i
 * hashes in a scratch hash_buf[].
 *
 * Usage:
 *   ./window_reduce_bench fold|sparse <outer_iters>
 *
 * "fold"  : register-heavy O(W) fold (copied from production reduce_window<21>).
 * "sparse": O(W log W) sparse table with explicit __m256i st[levels][width]
 *           written then read (intended to stress store-load forwarding).
 */
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <immintrin.h>
#include <utility>

namespace {

// ---- Copied from xseal_syncmer.hpp (reduce_window<21>) — production fold path ----

template <int W, size_t... J>
static inline __m256i reduce_window_impl(const uint32_t* hash_buf, int batch,
                                         std::index_sequence<J...>) {
  __m256i va = _mm256_loadu_si256((const __m256i*)&hash_buf[batch]);
  __m256i vb = _mm256_loadu_si256((const __m256i*)&hash_buf[batch + 1]);
  __m256i vc = _mm256_loadu_si256((const __m256i*)&hash_buf[batch + 2]);
  __m256i vd = _mm256_loadu_si256((const __m256i*)&hash_buf[batch + 3]);
  ([&]() __attribute__((always_inline)) {
    if constexpr (J % 4 == 0)
      va = _mm256_min_epu32(va, _mm256_loadu_si256((const __m256i*)&hash_buf[batch + J + 4]));
    else if constexpr (J % 4 == 1)
      vb = _mm256_min_epu32(vb, _mm256_loadu_si256((const __m256i*)&hash_buf[batch + J + 4]));
    else if constexpr (J % 4 == 2)
      vc = _mm256_min_epu32(vc, _mm256_loadu_si256((const __m256i*)&hash_buf[batch + J + 4]));
    else
      vd = _mm256_min_epu32(vd, _mm256_loadu_si256((const __m256i*)&hash_buf[batch + J + 4]));
  }(),
   ...);
  return _mm256_min_epu32(_mm256_min_epu32(va, vb), _mm256_min_epu32(vc, vd));
}

template <int W>
static inline __m256i reduce_window_fold(const uint32_t* hash_buf, int batch) {
  if constexpr (W <= 1)
    return _mm256_loadu_si256((const __m256i*)&hash_buf[batch]);
  else if constexpr (W == 2)
    return _mm256_min_epu32(_mm256_loadu_si256((const __m256i*)&hash_buf[batch]),
                            _mm256_loadu_si256((const __m256i*)&hash_buf[batch + 1]));
  else if constexpr (W == 3)
    return _mm256_min_epu32(
        _mm256_min_epu32(_mm256_loadu_si256((const __m256i*)&hash_buf[batch]),
                         _mm256_loadu_si256((const __m256i*)&hash_buf[batch + 1])),
        _mm256_loadu_si256((const __m256i*)&hash_buf[batch + 2]));
  else
    return reduce_window_impl<W>(hash_buf, batch, std::make_index_sequence<W - 4>{});
}

// ---- Sparse table RMQ on W SIMD minima (lane-wise), explicit memory levels ----

template <int W>
static inline int sparse_log2_floor(int n) {
  int k = 0;
  while ((1 << (k + 1)) <= n) ++k;
  return k;
}

template <int W>
static inline __m256i reduce_window_sparse(const uint32_t* hash_buf, int batch) {
  constexpr int kMax = 6; // enough for W <= 64
  alignas(64) __m256i st[kMax][W];

  for (int i = 0; i < W; ++i)
    st[0][i] = _mm256_loadu_si256((const __m256i*)&hash_buf[batch + i]);

  int max_k = sparse_log2_floor<W>(W);
  for (int k = 1; k <= max_k; ++k) {
    int half = 1 << (k - 1);
    int span = 1 << k;
    for (int i = 0; i + span <= W; ++i)
      st[k][i] = _mm256_min_epu32(st[k - 1][i], st[k - 1][i + half]);
  }

  int len = W;
  int lg = sparse_log2_floor<W>(len);
  int p2 = 1 << lg;
  int r = len - 1;
  if (p2 == len)
    return st[lg][0];
  return _mm256_min_epu32(st[lg][0], st[lg][r - p2 + 1]);
}

// Deterministic fill (no I/O): mimics spread of 32-bit hashes in hash_buf.
static void fill_hash_buf(uint32_t* buf, int words, uint64_t seed) {
  uint64_t x = seed;
  for (int i = 0; i < words; ++i) {
    x ^= x << 7;
    x ^= x >> 9;
    x *= 0xD6E8FEB86659FD93ULL;
    buf[i] = (uint32_t)(x >> 32) ^ (uint32_t)x;
  }
}

enum class Mode { Fold, Sparse };

static __m256i run_kernel(Mode m, const uint32_t* hash_buf, int64_t outer) {
  __m256i acc = _mm256_setzero_si256();
  constexpr int W = 21;
  for (int64_t t = 0; t < outer; ++t) {
    for (int batch = 0; batch < 128; batch += 8) {
      __m256i v =
          (m == Mode::Fold) ? reduce_window_fold<W>(hash_buf, batch)
                            : reduce_window_sparse<W>(hash_buf, batch);
      acc = _mm256_add_epi32(acc, v);
    }
  }
  return acc;
}

} // namespace

int main(int argc, char** argv) {
  if (argc < 3) {
    std::fprintf(stderr, "usage: %s fold|sparse <outer_iters>\n", argv[0]);
    return 2;
  }
  Mode mode = std::strcmp(argv[1], "fold") == 0   ? Mode::Fold
              : std::strcmp(argv[1], "sparse") == 0 ? Mode::Sparse
                                                    : Mode::Fold;
  if (std::strcmp(argv[1], "sparse") != 0 && std::strcmp(argv[1], "fold") != 0) {
    std::fprintf(stderr, "mode must be fold or sparse\n");
    return 2;
  }
  int64_t outer = std::atoll(argv[2]);
  if (outer < 1) return 2;

  constexpr int W = 21;
  constexpr int words = (128 + W + 8) * 8; // slop past tail batches
  alignas(64) static uint32_t hash_buf[4096];
  std::memset(hash_buf, 0, sizeof(hash_buf));
  fill_hash_buf(hash_buf, words, 0xC0FFEEULL);

  __m256i r = run_kernel(mode, hash_buf, outer);
  alignas(64) uint32_t sink[8];
  _mm256_storeu_si256((__m256i*)sink, r);
  volatile uint32_t vsink = sink[0] ^ sink[7];
  (void)vsink;
  std::fprintf(stdout, "CHECKSUM");
  for (int i = 0; i < 8; ++i) std::fprintf(stdout, "|%08x", sink[i]);
  std::fprintf(stdout, "\n");
  return 0;
}
