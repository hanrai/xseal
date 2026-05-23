// SPDX-License-Identifier: Apache-2.0
// Copyright (c) 2026 hanrai. All Rights Reserved.

#pragma once

#include <chrono>
#include <cstdint>
#include <iostream>
#include <string>
#include <type_traits>
#include <xseal/xseal_distributor.hpp>
#include <xseal/xseal_fasta_parser.hpp>
#include <xseal/xseal_fastq_parser.hpp>
#include <xseal/io/xseal_mmap_reader.hpp>

namespace xseal {

/// Optional wall-time breakdown (main thread segments are sequential; worker sums overlap).
struct PipelineProfileBreakdown {
  uint64_t distributor_construct_ns = 0;
  uint64_t get_next_ns = 0;
  uint64_t parse_ns = 0;
  uint64_t flush_ns = 0;
  uint64_t parser_finish_ns = 0;
  uint64_t distributor_stop_ns = 0;
};

// =========================================================================
// 1. Core Pipeline Executor
// The ReaderType and ParserType here are completely determined at compile time, without any polymorphism overhead.
// =========================================================================
template <typename ReaderType, template <typename, typename, typename, uint32_t> class ParserTemplate>
class PipelineExecutor {
public:
  static constexpr uint32_t NumSlots = ReaderType::NUM_SLOTS;

  PipelineExecutor(ReaderType &reader, XsealScannerType scanner_type,
                   int num_threads, const std::string &reader_name,
                   bool force_fastq, bool is_fastq, bool profile_breakdown = false,
                   size_t batch_threshold = kXsealDefaultDistributorBatchBytes)
      : reader_(reader), scanner_type_(scanner_type), num_threads_(num_threads),
        reader_name_(reader_name), force_fastq_(force_fastq),
        is_fastq_(is_fastq), profile_(profile_breakdown),
        batch_threshold_(batch_threshold) {}

  void execute() {
    int k = 31, s = 11;
    PipelineProfileBreakdown pb;
    uint64_t distributor_construct_ns = 0;
    std::chrono::high_resolution_clock::time_point ctor_start;
    if (profile_)
      ctor_start = std::chrono::high_resolution_clock::now();
    XsealDistributor<NumSlots> distributor(num_threads_, k, s, scanner_type_,
                                           batch_threshold_);
    if (profile_) {
      distributor_construct_ns =
          std::chrono::duration_cast<std::chrono::nanoseconds>(
              std::chrono::high_resolution_clock::now() - ctor_start)
              .count();
      pb.distributor_construct_ns = distributor_construct_ns;
    }

    const uint8_t *current_chunk_data = nullptr;
    uint64_t current_chunk_off = 0;

    // [Performance Key Point]: Use lambda to ensure these callbacks can be inlined into the Parser
    auto on_bases = [&](const SeqSegment &seg,
                        const XsealSlot<NumSlots> &slot = XsealSlot<NumSlots>()) {
      distributor.process_bases(seg, current_chunk_data, current_chunk_off,
                                 slot);
    };
    auto on_header = [&](const SeqSegment &, const XsealSlot<NumSlots> & = XsealSlot<NumSlots>()) {
      distributor.reset_sequence();
    };
    auto on_gap = [&](const SeqSegment &, const XsealSlot<NumSlots> & = XsealSlot<NumSlots>()) {
      distributor.reset_sequence();
    };

    ParserTemplate<decltype(on_bases), decltype(on_header), decltype(on_gap), NumSlots>
        parser(on_bases, on_header, on_gap);

    auto start = std::chrono::high_resolution_clock::now();

    // Core hot loop
    while (true) {
      size_t len;
      uint64_t off;
      XsealSlot<NumSlots> slot;
      if (profile_) {
        auto tn0 = std::chrono::high_resolution_clock::now();
        slot = reader_.getNext(len, off);
        pb.get_next_ns += std::chrono::duration_cast<std::chrono::nanoseconds>(
                              std::chrono::high_resolution_clock::now() - tn0)
                              .count();
      } else {
        slot = reader_.getNext(len, off);
      }
      if (!slot.valid())
        break;

      const uint8_t *data_ptr = get_data_ptr(slot, off);
      current_chunk_data = data_ptr;
      current_chunk_off = off;

      if (profile_) {
        auto tp0 = std::chrono::high_resolution_clock::now();
        parser.parse(data_ptr, len, off, slot);
        pb.parse_ns += std::chrono::duration_cast<std::chrono::nanoseconds>(
                           std::chrono::high_resolution_clock::now() - tp0)
                           .count();
        auto tf0 = std::chrono::high_resolution_clock::now();
        distributor.flush();
        pb.flush_ns += std::chrono::duration_cast<std::chrono::nanoseconds>(
                           std::chrono::high_resolution_clock::now() - tf0)
                           .count();
      } else {
        parser.parse(data_ptr, len, off, slot);
        distributor.flush();
      }
    }

    if (profile_) {
      auto tf0 = std::chrono::high_resolution_clock::now();
      parser.finish();
      pb.parser_finish_ns =
          std::chrono::duration_cast<std::chrono::nanoseconds>(
              std::chrono::high_resolution_clock::now() - tf0)
              .count();
      auto ts0 = std::chrono::high_resolution_clock::now();
      distributor.stop();
      pb.distributor_stop_ns =
          std::chrono::duration_cast<std::chrono::nanoseconds>(
              std::chrono::high_resolution_clock::now() - ts0)
              .count();
    } else {
      parser.finish();
      distributor.stop();
    }

    if (profile_) {
      print_stats(start, distributor.get_total_bases(),
                  distributor.get_total_hits(), distributor.get_total_scan_ns(),
                  distributor.get_total_encode_ns(), distributor.get_total_harvest_ns(), &pb);
    } else {
      print_stats(start, distributor.get_total_bases(),
                  distributor.get_total_hits(), distributor.get_total_scan_ns(),
                  distributor.get_total_harvest_ns());
    }
  }

private:
  ReaderType &reader_;
  XsealScannerType scanner_type_;
  int num_threads_;
  std::string reader_name_;
  bool force_fastq_;
  bool is_fastq_;
  bool profile_;
  size_t batch_threshold_;

  // Compile-time branch: Determine how to obtain the data pointer based on the Reader type (zero runtime overhead)
  const uint8_t *get_data_ptr(const XsealSlot<NumSlots> &slot, uint64_t off) {
    if constexpr (is_xseal_mmap_reader_t_v<ReaderType>) {
      return reinterpret_cast<const uint8_t *>(reader_.get_mapped_data()) + off;
    } else {
      return reader_.get_pool_base() + slot.id() * reader_.chunk_size();
    }
  }

  void print_stats(
      std::chrono::time_point<std::chrono::high_resolution_clock> start,
      uint64_t total_bases, uint64_t total_hits, uint64_t total_scan_ns,
      uint64_t total_harvest_ns) {
    print_stats(start, total_bases, total_hits, total_scan_ns, 0, total_harvest_ns, nullptr);
  }

  void print_stats(
      std::chrono::time_point<std::chrono::high_resolution_clock> start,
      uint64_t total_bases, uint64_t total_hits, uint64_t total_scan_ns,
      uint64_t total_encode_ns, uint64_t total_harvest_ns, const PipelineProfileBreakdown *pb) {
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> diff = end - start;
    double file_bytes = reader_.get_file_size();

    std::cout << "\n====================================================\n";
    std::cout << " 🚀 XSeal V4.0 Robust-Elite (New Paradigm Architecture)\n";
    std::cout << "----------------------------------------------------\n";
    std::cout << " 🛠️  Scanner Mode         : "
              << (scanner_type_ == XsealScannerType::SYNCMER ? "Syncmer"
                                                             : "Minimizer")
              << "\n";
    std::cout << " 📂 Input Format          : " << (is_fastq_ ? "FASTQ" : "FASTA")
              << (force_fastq_ ? " (Forced)" : "") << "\n";
    std::cout << " 📥 Reader Type           : " << reader_name_ << "\n";
    const double chunk_mib =
        static_cast<double>(reader_.chunk_size()) / (1024.0 * 1024.0);
    const double batch_mib =
        static_cast<double>(batch_threshold_) / (1024.0 * 1024.0);
    std::cout << " 📦 Chunk Size           : " << chunk_mib << " MiB\n";
    std::cout << " 📐 Distributor batch     : " << batch_mib << " MiB (raw)\n";
    std::cout << " 🧵 Threads               : " << num_threads_ << "\n";
    std::cout << " 🧬 Total Bases processed : " << total_bases << " bp\n";
    std::cout << " 🎯 Syncmer Hits found    : " << total_hits << "\n";
    std::cout << " 🏎️  Base Throughput      : "
              << (total_bases / 1e9) / diff.count() << " Gbp/s\n";
    std::cout << " 💾 Raw File Throughput  : "
              << (file_bytes / 1e9) / diff.count() << " GB/s\n";
    std::cout << " ⏱️  Total Time           : " << diff.count() << " s\n";
    std::cout << " ⚡ Core Scan Time (Wall) : " << (total_scan_ns / 1e9) / num_threads_ << " s\n";
    std::cout << " 📊 Core Scan Throughput : " << (total_bases / 1e9) / ((total_scan_ns / 1e9) / num_threads_) << " Gbp/s\n";
    std::cout << " 🌾 Harvest Throughput  : " << (total_hits / 1e6) / ((total_harvest_ns / 1e9) / num_threads_) << " M-hits/s\n";

    if (pb != nullptr) {
      std::cout << "----------------------------------------------------\n";
      std::cout << " Profile (--profile): wall-time sums (see note)\n";
      std::cout << "  Note: main-thread segments sum to serial pipeline time;\n";
      std::cout << "  worker encode/scan sums overlap with parse/flush and each other.\n";
      std::cout << "  distributor_construct_ns : " << pb->distributor_construct_ns << "\n";
      std::cout << "  get_next_ns (sum)          : " << pb->get_next_ns << "\n";
      std::cout << "  parse_ns (sum)             : " << pb->parse_ns << "\n";
      std::cout << "  flush_ns (sum)             : " << pb->flush_ns << "\n";
      std::cout << "  parser_finish_ns           : " << pb->parser_finish_ns << "\n";
      std::cout << "  distributor_stop_ns      : " << pb->distributor_stop_ns << "\n";
      std::cout << "  worker total_encode_ns (sum threads): " << total_encode_ns << "\n";
      std::cout << "  worker total_scan_ns (sum threads)  : " << total_scan_ns << "\n";
      std::cout << "  worker total_harvest_ns (sum threads): " << total_harvest_ns << "\n";
    }
    std::cout << "====================================================\n\n";
  }
};

// =========================================================================
// 2. Parser Dispatch Factory
// Convert runtime bool is_fastq into compile-time types
// =========================================================================
template <typename ReaderType>
void dispatch_parser(ReaderType &reader, bool is_fastq,
                     XsealScannerType scanner_type, int num_threads,
                     const std::string &reader_name, bool force_fastq,
                     bool profile_breakdown = false,
                     size_t batch_threshold = kXsealDefaultDistributorBatchBytes) {
  if (is_fastq) {
    PipelineExecutor<ReaderType, XsealFastqParserT> executor(
        reader, scanner_type, num_threads, reader_name, force_fastq, is_fastq,
        profile_breakdown, batch_threshold);
    executor.execute();
  } else {
    PipelineExecutor<ReaderType, XsealFastaParserT> executor(
        reader, scanner_type, num_threads, reader_name, force_fastq, is_fastq,
        profile_breakdown, batch_threshold);
    executor.execute();
  }
}

} // namespace xseal
