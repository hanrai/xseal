// SPDX-License-Identifier: Apache-2.0
#include <memory>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <vector>
#include <xseal/io/xseal_mmap_reader.hpp>
#include <xseal/xseal_encoder.hpp>
#include <xseal/xseal_fasta_parser.hpp>
#include <xseal/xseal_syncmer.hpp>

namespace py = pybind11;

struct FullResults {
  py::array_t<uint32_t> positions;
  py::array_t<uint64_t> hashes;
};

py::array_t<uint32_t> extract_positions(const std::string &path, int k, int s) {
  xseal::XsealMmapReader reader(path);
  xseal::XSealEncoder encoder;
  xseal::XSealSyncmer scanner(k, s);
  std::vector<uint32_t> all_pos;

  size_t enc_buf_size = 16 * 1024 * 1024 / 4 + 4096;
  void *enc_buf_raw = nullptr;
  if (posix_memalign(&enc_buf_raw, 32, enc_buf_size) != 0)
    throw std::runtime_error("memalign failed");
  std::unique_ptr<void, decltype(&free)> enc_ptr(enc_buf_raw, free);
  uint8_t *enc_buf = static_cast<uint8_t *>(enc_buf_raw);

  std::vector<xseal::SyncmerHit> hit_buf(1024 * 1024);

  reader.read_sequential([&](const uint8_t *data, size_t len, uint64_t off) {
    auto on_bases = [&](const xseal::SeqSegment &seg) {
      xseal::XSealEncoderState state;
      encoder.init_state(&state, enc_buf, enc_buf_size);
      const char *ptr = (const char *)data + (seg.file_offset - off);
      const char *end = ptr + seg.byte_len;

      while (ptr < end) {
        auto status = encoder.encode_chunk(&state, &ptr, end);
        encoder.flush_state(&state);
        size_t num_chunks = (state.out_byte_idx + 31) / 32;
        size_t num_hits =
            scanner.scan<21, true>(reinterpret_cast<const __m256i *>(enc_buf),
                                   num_chunks, 0, (uint32_t *)hit_buf.data());

        for (size_t i = 0; i < num_hits; ++i)
          all_pos.push_back(hit_buf[i].pos);
        if (status == xseal::XSEAL_ENC_BUFFER_FULL) {
          encoder.init_state(&state, enc_buf, enc_buf_size);
        } else {
          break;
        }
      }
    };
    auto on_header = [](const xseal::SeqSegment &) {};
    auto on_gap = [](const xseal::SeqSegment &) {};

    xseal::XsealFastaParserT parser(on_bases, on_header, on_gap);
    parser.process_chunk(data, len, off);
    parser.finish();
  });

  auto pos_arr = py::array_t<uint32_t>(all_pos.size());
  std::memcpy(pos_arr.mutable_data(), all_pos.data(),
              all_pos.size() * sizeof(uint32_t));
  return pos_arr;
}

FullResults extract_full(const std::string &path, int k, int s) {
  xseal::XsealMmapReader reader(path);
  xseal::XSealEncoder encoder;
  xseal::XSealSyncmer scanner(k, s);
  std::vector<uint32_t> all_pos;
  std::vector<uint64_t> all_hashes;

  size_t enc_buf_size = 16 * 1024 * 1024 / 4 + 4096;
  void *enc_buf_raw = nullptr;
  if (posix_memalign(&enc_buf_raw, 32, enc_buf_size) != 0)
    throw std::runtime_error("memalign failed");
  std::unique_ptr<void, decltype(&free)> enc_ptr(enc_buf_raw, free);
  uint8_t *enc_buf = static_cast<uint8_t *>(enc_buf_raw);

  std::vector<xseal::SyncmerHitFull> hit_buf(1024 * 1024);

  reader.read_sequential([&](const uint8_t *data, size_t len, uint64_t off) {
    auto on_bases = [&](const xseal::SeqSegment &seg) {
      xseal::XSealEncoderState state;
      encoder.init_state(&state, enc_buf, enc_buf_size);
      const char *ptr = (const char *)data + (seg.file_offset - off);
      const char *end = ptr + seg.byte_len;

      while (ptr < end) {
        auto status = encoder.encode_chunk(&state, &ptr, end);
        encoder.flush_state(&state);
        size_t num_chunks = (state.out_byte_idx + 31) / 32;
        size_t num_hits = scanner.scan_with_hash<21, true>(
            reinterpret_cast<const __m256i *>(enc_buf), num_chunks, 0,
            hit_buf.data());

        for (size_t i = 0; i < num_hits; ++i) {
          all_pos.push_back(hit_buf[i].pos);
          all_hashes.push_back(hit_buf[i].hash);
        }
        if (status == xseal::XSEAL_ENC_BUFFER_FULL) {
          encoder.init_state(&state, enc_buf, enc_buf_size);
        } else {
          break;
        }
      }
    };
    auto on_header = [](const xseal::SeqSegment &) {};
    auto on_gap = [](const xseal::SeqSegment &) {};

    xseal::XsealFastaParserT parser(on_bases, on_header, on_gap);
    parser.process_chunk(data, len, off);
    parser.finish();
  });

  auto pos_arr = py::array_t<uint32_t>(all_pos.size());
  auto hash_arr = py::array_t<uint64_t>(all_hashes.size());
  std::memcpy(pos_arr.mutable_data(), all_pos.data(),
              all_pos.size() * sizeof(uint32_t));
  std::memcpy(hash_arr.mutable_data(), all_hashes.data(),
              all_hashes.size() * sizeof(uint64_t));

  return {pos_arr, hash_arr};
}

PYBIND11_MODULE(xsealsyncmer, m) {
  m.doc() = "XSealSyncmer: Ultra-high throughput Syncmer extraction";

  py::class_<FullResults>(m, "FullResults")
      .def_readonly("positions", &FullResults::positions)
      .def_readonly("hashes", &FullResults::hashes);

  m.def("extract_positions", &extract_positions,
        "Extract only positions of Syncmers from a FASTA file", py::arg("path"),
        py::arg("k") = 31, py::arg("s") = 11);

  m.def("extract_full", &extract_full,
        "Extract positions and 64-bit hashes of Syncmers from a FASTA file",
        py::arg("path"), py::arg("k") = 31, py::arg("s") = 11);
}
