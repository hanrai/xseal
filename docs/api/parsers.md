# XSeal Parser Subsystem: Design & Usage

This document details the XSeal high-performance parsing subsystem, focusing on the SIMD-accelerated FASTA/FASTQ state machines.

## 1. Overview & Design Principles

XSeal parsers bridge the gap between high-speed readers (up to 13 GB/s) and SIMD-optimized processing pipelines (e.g., k-mer encoders).

### Core Principles:
- **Zero-Copy Streaming**: Data is never copied out of the reader's buffer. The parser emits lightweight `SeqSegment` references (pointers/offsets) linked to active `XsealSlot` buffers.
- **SIMD Scanning**: Uses AVX2 to scan for record boundaries, newlines, and special characters 32 or 128 bytes at a time.
- **Extreme Throughput**: Reaches **21.87 GB/s** single-core throughput on genomic data by minimizing micro-instructions and memory overhead.
- **Unified Slot Integration**: Designed to be integrated directly with the unified slot pool architecture, passing memory slot references through to the next stages of the pipeline.

---

## 2. Parser Implementations

### A. Core FASTA Engine: `XsealFastaParser`
A non-templated, ultra-low-overhead FASTA scanner that processes a memory block and appends metadata records to a run-length table (`XsealRunTable` or `std::vector<SeqSegment>`).
```cpp
class XsealFastaParser {
public:
    XsealFastaParser();
    void parse(const uint8_t *data, size_t size, uint64_t absolute_offset, XsealRunTable &run_table);
    void finish(XsealRunTable &run_table);
};
```

### B. Callback Wrapper: `XsealFastaParserT`
A template-based compatibility layer wrapping `XsealFastaParser`. It decodes sequence segments and fires optimized callback handlers, forwarding their corresponding active slot.
```cpp
template <typename BasesHandler, typename HeaderHandler, typename GapHandler, uint32_t NumSlots = 64>
class XsealFastaParserT {
public:
    XsealFastaParserT(BasesHandler b, HeaderHandler h, GapHandler g);
    void parse(const uint8_t *data, size_t size, uint64_t absolute_offset, const XsealSlot<NumSlots> &slot = XsealSlot<NumSlots>());
    void finish();
};
```

### C. FASTQ Engine: `XsealFastqParserT`
A template-based, AVX2-accelerated parser for FASTQ records.
```cpp
template <typename BasesHandler, typename HeaderHandler, typename SkipHandler, uint32_t NumSlots = 64>
class XsealFastqParserT {
public:
    XsealFastqParserT(BasesHandler b, HeaderHandler h, SkipHandler s);
    void parse(const uint8_t *data, size_t size, uint64_t absolute_offset, const XsealSlot<NumSlots> &slot = XsealSlot<NumSlots>());
    void finish();
};
```

---

## 3. Segment Representation

The parser emits sequence segments as `SeqSegment` structures:
```cpp
struct SeqSegment {
    uint64_t file_offset;      // Start of segment in the file
    uint32_t logical_len;      // Length of actual sequence data (excluding \n)
    uint32_t byte_len;         // Physical length in buffer (including \n)
    SegType  type;             // SegType::HEADER, SegType::BASES, SegType::N_GAP, or SegType::SKIP
    bool     is_physically_truncated; // True if segment was cut by a chunk boundary
};
```

---

## 4. Usage Guide

### Basic FASTA Parsing Loop
```cpp
#include <xseal/xseal_fasta_parser.hpp>
#include <xseal/io/xseal_mmap_reader.hpp>

// Define callback handlers accepting segments and slots
auto on_header = [](const xseal::SeqSegment& seg, const xseal::XsealSlot<768>& slot) {
    // Process FASTA Header...
};
auto on_bases = [](const xseal::SeqSegment& seg, const xseal::XsealSlot<768>& slot) {
    // Process DNA Bases...
};
auto on_gap = [](const xseal::SeqSegment& seg, const xseal::XsealSlot<768>& slot) {
    // Process Gaps (N/Ambiguous)...
};

// Initialize FASTA parser with handlers and slot pool size matching the reader
xseal::XsealFastaParserT<decltype(on_bases), decltype(on_header), decltype(on_gap), 768> parser(on_bases, on_header, on_gap);

xseal::XsealMmapReader reader("data/hg38.fa");
size_t len = 0;
uint64_t offset = 0;

while (true) {
    auto slot = reader.getNext(len, offset);
    if (!slot.valid()) {
        break; // EOF
    }

    const uint8_t* chunk_data = reader.get_pool_base() ? (reader.get_pool_base() + slot.id() * reader.chunk_size()) : (const uint8_t*)(reader.get_mapped_data() + offset);

    // Parse the chunk
    parser.parse(chunk_data, len, offset, slot);
}
parser.finish(); // Flush final segment
```

---

## 5. Performance Benchmarks

Measured on a single core (Intel/AMD, hg38.fa, 3.1 GB):

| Parser Implementation | Throughput (GB/s) |
| :--- | :--- |
| **XsealFastaParser (AVX2)** | **21.87 GB/s** |
| Scalar Baseline | 2.63 GB/s |

*Note: Achieving this speed requires pre-loading data into RAM to avoid OS page fault overhead.*
