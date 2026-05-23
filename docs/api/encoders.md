# XSeal Encoders: Design & Usage

XSeal provides two state-of-the-art DNA sequence encoders designed to convert ASCII DNA bases (`A`, `C`, `G`, `T`, and their lowercase forms) into a compact 2-bit representation (A=0, C=1, G=2, T=3).

---

## 1. `XSealEncoder`: Production-Grade State Machine Encoder

`XSealEncoder` is designed to process physical FASTA/FASTQ sequence chunks that contain newlines or formatting breaks, ignoring those breaks while streaming 2-bit output into the target buffer.

### Processing Paths
- **Turbo Madden Path (God Speed)**: When a block of 32 or 128 bytes contains only pure ACGT bases (no newlines), the encoder uses `_mm256_maddubs_epi16` to perform high-speed vertical bit-packing, reaching up to 40 GB/s on a single core.
- **PEXT Path (Elite Speed)**: If a block contains newlines, it falls back to the PEXT path. It uses `_pext_u64` BMI2 instructions to squeeze out 2-bit values for valid bases while ignoring newlines.
- **Scalar Fallback**: Handles remaining bytes at the end of a chunk or small segments.

### Usage Example
```cpp
#include <xseal/xseal_encoder.hpp>

xseal::XSealEncoder encoder;
xseal::XSealEncoderState state;
std::vector<uint8_t> out_buffer(1024 * 1024);

// Initialize state
encoder.init_state(&state, out_buffer.data(), out_buffer.size());

const char* data = "ACGTACGT\nACGTACGT";
const char* p = data;
const char* end = data + strlen(data);

// Encode chunk (updates pointer p)
xseal::XSealEncoderStatus status = encoder.encode_chunk(&state, &p, end);

// Flush remaining bits (padding to nearest byte boundary)
encoder.flush_state(&state);
```

### Status Codes
- `XSEAL_ENC_OK`: Normal operation (intermediate).
- `XSEAL_ENC_EOF`: Reached the end of the input buffer.
- `XSEAL_ENC_BUFFER_FULL`: Output buffer capacity reached.

---

## 2. `PureEncoder`: Hardware-Optimized Stream Encoder

`PureEncoder` is a zero-branch, full-SIMD streaming encoder designed for pre-cleaned genomic datasets (where special characters like `\n` or `N` have already been filtered out).

### Design & Optimizations
- **Zero-Branch Hot Loop**: Processes 32-byte (256-bit) blocks aggressively without branch checks, bit counters, or state preservation.
- **Cute-Nucleotides Extraction**: Leverages the binary representation of ASCII characters (`A=0x41`, `C=0x43`, `G=0x47`, `T=0x54`) to extract the second and first bits using shift and movemask (`_mm256_movemask_epi8`) primitives.
- **BMI2 Interleaving**: Uses double-interleaved `_pdep_u64` to pack the extracted bits into a 64-bit unsigned integer in a single operation.

### Usage Example
```cpp
#include <pure_encoder.hpp>

const char* raw_seq = "ACGTACGTACGTACGTACGTACGTACGTACGT"; // 32 bases
const char* end = raw_seq + 32;

// Output space: requires at least (end - start) / 4 + 8 bytes
std::vector<uint8_t> encoded_seq(64);

xseal::PureEncoder::encode_chunk(raw_seq, end, encoded_seq.data());
```
