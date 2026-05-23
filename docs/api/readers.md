# XSeal Reader Subsystem: Design & Usage

This document provides an overview of the XSeal V4.0 high-performance reader subsystem, detailing its design principles, architecture, usage, and performance characteristics.

## 1. Overview & Design Principles

XSeal implements a suite of high-performance sequential readers designed to feed genomic data into SIMD-optimized processing pipelines.

### Core Principles:
- **Unified Pull-Based API**: Instead of push-based callbacks, all readers provide a unified pull-based method using slot pool backpressure:
  ```cpp
  XsealSlot<NumSlots> getNext(size_t &out_len, uint64_t &out_offset);
  ```
  Returns a valid `XsealSlot` containing a pre-allocated aligned buffer slot when data is available. Returns an invalid/empty slot on EOF.
- **Natural Backpressure**: Working threads and reader threads are synchronized via a compile-time bounded `XsealSlotPool<NumSlots>`. If all slots in the pool are currently occupied by downstream workers, the reader naturally blocks or spins with nanosecond-level sleep (`_mm_pause()`) until slots are returned.
- **Zero-Copy / Registered Buffers**: Uses pre-allocated memory aligned to 4096-byte page boundaries, optimizing for `O_DIRECT` and registered buffers (`io_uring_register_buffers`) to eliminate OS copy overhead.

---

## 2. Reader Architectures

### `XsealMmapReaderT<uint32_t NumSlots>`
Uses the Linux `mmap` system call to map the entire file into virtual memory.
- **Mechanism**: Employs `madvise` hints (`MADV_SEQUENTIAL`, `MADV_HUGEPAGE`, `MADV_WILLNEED`) to trigger aggressive kernel-level page-table pre-fetching.
- **Pros**: Highest throughput for cached files (~13.5 GB/s).
- **Cons**: Performance is dependent on the host page cache.

### `XsealUringReaderT<uint32_t NumSlots>`
Leverages the modern Linux `io_uring` interface for fully asynchronous and zero-blocking NVMe/SSD storage operations.
- **Mechanism**: Registers the pre-allocated slot buffers (`io_uring_register_buffers`) and direct file descriptors to completely bypass the kernel page cache via `O_DIRECT`. A sliding ring of submit sequences (ROB) handles out-of-order kernel completions.
- **Pros**: Extremely low CPU overhead and high concurrent throughput (~6.3 GB/s) on NVMe storage.

### `XsealPipeReaderT<uint32_t NumSlots>`
Designed for stream-processing from standard input (`stdin`) or Unix pipelines.
- **Mechanism**: Uses `io_uring` with a queue depth of 1 (`QD=1`) and a custom "Read-Until-Full" sliding state machine to systematically resolve short-read fragmentation and interleaving over pipe boundaries.
- **Pros**: Perfectly handles non-seekable streaming pipelines.

---

## 3. Usage Guide

### Basic Pull Loop Usage
```cpp
#include <xseal/io/xseal_mmap_reader.hpp>

// Create a reader (defaults to 768 slots and 2MB chunk size)
xseal::XsealMmapReader reader("data/hg38.fa", 2 * 1024 * 1024);

size_t len = 0;
uint64_t offset = 0;

while (true) {
    // getNext pulls the next available data slot
    auto slot = reader.getNext(len, offset);
    if (!slot.valid()) {
        break; // EOF reached
    }

    // Access raw chunk data using slot.id() into the reader's buffer pool base
    const uint8_t* chunk_data = reader.get_mapped_data() + offset; 
    
    // Process chunk_data (len bytes)...
    
    // Slot goes out of scope and automatically returns to the pool
}
```

### Configuring custom slots
```cpp
// Allocate a smaller slot pool (e.g., 256 slots) for restricted memory environments
xseal::XsealUringReaderT<256> reader("data/hg38.fa", 2 * 1024 * 1024);
```

---

## 4. Performance Benchmarks

Measured on `data/hg38.fa` (3.27 GB) using **24MB** chunks:

| Reader | Throughput | Bottleneck |
| :--- | :--- | :--- |
| **MMAP** | **13.53 GB/s** | RAM Bandwidth |
| **IO_URING** | **6.33 GB/s** | Disk I/O / Latency |
| **PIPE (RAM)**| **2.13 GB/s** | Pipe Copy / Context Switch |

---

## 5. Testing & Verification

Unit and stress tests verify:
- **Natural Backpressure**: Slot pools block correctly when consumers are artificially stalled.
- **Direct I/O Boundary Integrity**: Page and size alignments function correctly across arbitrary file lengths.
- **Short-Read Robustness**: Pipe readers reconstruct streamed data bit-perfectly without fragmenting records.
