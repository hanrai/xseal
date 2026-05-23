// SPDX-License-Identifier: Apache-2.0
// Copyright (c) 2026 hanrai. All Rights Reserved.

// Minimal micro-benchmark: classic size-mismatch store-then-load vs register-only.
// Used with `perf stat -e ls_bad_status2.stli_other,...` on AMD Zen 3 (Family 19h).
#include <cinttypes>
#include <cstdint>
#include <cstdio>
#include <cstring>

alignas(64) static volatile unsigned char buf[128];
static volatile std::uint64_t sink;

static void run_bad(std::uint64_t iters) {
    // Agner-style pattern: 32-bit store then 64-bit load overlapping same line.
    for (std::uint64_t i = 0; i < iters; ++i) {
        *reinterpret_cast<volatile std::uint32_t*>(buf + 0) =
            static_cast<std::uint32_t>(i);
        sink = *reinterpret_cast<volatile std::uint64_t*>(buf + 0);
    }
}

static void run_good(std::uint64_t iters) {
    std::uint64_t x = 1;
    for (std::uint64_t i = 0; i < iters; ++i) {
        x ^= (i * 0x9E3779B185EBCA87ULL) + (x >> 7);
    }
    sink = x;
}

int main(int argc, char** argv) {
    if (argc < 3) {
        std::fprintf(stderr, "usage: %s {bad|good} <iters>\n", argv[0]);
        return 2;
    }
    const char* mode = argv[1];
    std::uint64_t iters = 0;
    if (std::sscanf(argv[2], "%" SCNu64, &iters) != 1 || iters == 0)
        return 2;
    std::memset(const_cast<unsigned char*>(buf), 0, sizeof(buf));
    if (std::strcmp(mode, "bad") == 0)
        run_bad(iters);
    else if (std::strcmp(mode, "good") == 0)
        run_good(iters);
    else
        return 2;
    return 0;
}
