// SPDX-License-Identifier: Apache-2.0
// Copyright (c) 2026 hanrai. All Rights Reserved.

#pragma once

#include <algorithm>
#include <cstdint>
#include <type_traits>
#include <fcntl.h>
#include <immintrin.h>
#include <stdexcept>
#include <string>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>
#include <xseal/xseal_slot_pool.hpp>

namespace xseal {

/**
 * XsealMmapReader: Fast mmap-based sequential reader.
 * V4.0 Architecture: Unified Pull-based API with Slot backpressure.
 */
template <uint32_t NumSlots = 768> class XsealMmapReaderT {
public:
  static constexpr uint32_t NUM_SLOTS = NumSlots;

  explicit XsealMmapReaderT(const std::string &path,
                            size_t chunk_bytes = 2 * 1024 * 1024)
      : chunk_size_(chunk_bytes) {
    if (chunk_size_ == 0)
      throw std::runtime_error("XsealMmapReader: chunk_bytes must be > 0");

    fd = open(path.c_str(), O_RDONLY);
    if (fd == -1)
      throw std::runtime_error("Cannot open file: " + path);

    struct stat st;
    if (fstat(fd, &st) != 0) {
      close(fd);
      throw std::runtime_error("fstat failed");
    }
    file_size = st.st_size;

    posix_fadvise(fd, 0, file_size, POSIX_FADV_SEQUENTIAL);
    posix_fadvise(fd, 0, file_size, POSIX_FADV_WILLNEED);

    if (file_size > 0) {
      mapped_data =
          (const char *)mmap(NULL, file_size, PROT_READ, MAP_PRIVATE, fd, 0);
      if (mapped_data == MAP_FAILED) {
        close(fd);
        throw std::runtime_error("mmap failed");
      }
      madvise((void *)mapped_data, file_size, MADV_SEQUENTIAL);
      madvise((void *)mapped_data, file_size, MADV_HUGEPAGE);
      madvise((void *)mapped_data, file_size, MADV_WILLNEED);
    }
  }

  ~XsealMmapReaderT() {
    if (mapped_data != MAP_FAILED && mapped_data != nullptr)
      munmap((void *)mapped_data, file_size);
    if (fd != -1)
      close(fd);
  }

  uint64_t get_file_size() const { return file_size; }
  const char *get_mapped_data() const { return mapped_data; }
  size_t chunk_size() const noexcept { return chunk_size_; }

  /**
   * V4.0 Pull API: Provides chunks naturally limited by SlotPool backpressure.
   */
  XsealSlot<NumSlots> getNext(size_t &out_len, uint64_t &out_offset) {
    if (current_offset >= file_size)
      return XsealSlot<NumSlots>(); // EOF

    while (true) {
      XsealSlot<NumSlots> slot = slot_pool.acquire();
      if (slot.valid()) {
        out_offset = current_offset;
        out_len =
            std::min(chunk_size_, (size_t)(file_size - current_offset));
        current_offset += out_len;
        return slot;
      }
      // Natural backpressure: if all 64 boxes are with the Worker, main thread spin-yields for nanoseconds, waiting for return
      for (int i = 0; i < 16; ++i)
        _mm_pause();
    }
  }

private:
  const size_t chunk_size_;
  int fd = -1;
  size_t file_size = 0;
  uint64_t current_offset = 0;
  const char *mapped_data = (const char *)MAP_FAILED;

  XsealSlotPool<NumSlots> slot_pool;
};

using XsealMmapReader = XsealMmapReaderT<>;

template <class T> struct is_xseal_mmap_reader_t : std::false_type {};

template <uint32_t NumSlots>
struct is_xseal_mmap_reader_t<XsealMmapReaderT<NumSlots>> : std::true_type {};

template <class T>
inline constexpr bool is_xseal_mmap_reader_t_v = is_xseal_mmap_reader_t<T>::value;

} // namespace xseal
