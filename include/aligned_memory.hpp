// SPDX-License-Identifier: Apache-2.0
// Copyright (c) 2026 hanrai. All Rights Reserved.

#pragma once

#include <cstddef>
#include <cstdlib>
#include <memory>
#include <stdexcept>

#ifdef _WIN32
#include <malloc.h>
#endif

namespace xseal {

    /**
     * 🛡️ AlignedDeleter: Ultra-stable cross-platform memory deallocator
     * Uses _aligned_free (Windows) or std::free (Unix) to ensure full compatibility with the CRT heap,
     * avoiding NT kernel-level segmentation faults potentially caused by VirtualFree.
     */
    struct AlignedDeleter {
        void operator()(void* ptr) const {
            if (!ptr) return;
#ifdef _WIN32
            _aligned_free(ptr);
#else
            std::free(ptr);
#endif
        }
    };

    using AlignedMemPtr = std::unique_ptr<uint8_t[], AlignedDeleter>;

    /**
     * 🚀 allocate_aligned_buyout: Physical memory "buyout" allocation
     * Uses 64-byte alignment to perfectly fit AVX-512 / AVX2 Cache Line boundaries,
     * eliminating CPU pipeline stalls caused by unaligned accesses.
     */
    inline AlignedMemPtr allocate_aligned_buyout(size_t size, size_t alignment = 64) {
        size_t padded_size = (size + alignment - 1) & ~(alignment - 1);
        void* ptr = nullptr;

#ifdef _WIN32
        ptr = _aligned_malloc(padded_size, alignment);
#else
        ptr = std::aligned_alloc(alignment, padded_size);
#endif
        if (!ptr) throw std::bad_alloc();
        return AlignedMemPtr(static_cast<uint8_t*>(ptr));
    }

} // namespace xseal