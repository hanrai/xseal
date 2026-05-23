// SPDX-License-Identifier: Apache-2.0
// Copyright (c) 2026 hanrai. All Rights Reserved.

#include <CLI/CLI.hpp>
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <fcntl.h>
#include <iostream>
#include <string>
#include <sys/stat.h>
#include <thread>

#include <unistd.h>
#include <xseal/xseal_pipeline_executor.hpp>
#include <xseal/io/xseal_mmap_reader.hpp>
#include <xseal/io/xseal_pipe_reader.hpp>
#include <xseal/io/xseal_uring_reader.hpp>

using namespace xseal;

namespace {

constexpr size_t kAlign4096(size_t x) noexcept {
  return (x + 4095) & ~size_t(4095);
}

/// MiB (1024²) → bytes: round from double, then at least one 4KiB page (aligned).
size_t chunk_batch_bytes_from_mib(double mib) {
  const size_t raw =
      static_cast<size_t>(std::llround(mib * 1024.0 * 1024.0));
  return kAlign4096(std::max(raw, size_t(4096)));
}

} // namespace

int main(int argc, char **argv) {
  CLI::App app{"🚀 XSeal Robust-Elite Benchmark Tool"};

  std::string path;
  int num_threads = std::thread::hardware_concurrency();
  std::string mode_str = "syncmer";
  std::string reader_type = "mmap";
  bool force_fastq = false;
  bool profile_breakdown = false;
  double chunk_batch_mib = 0.25;

  app.add_option("path", path,
                 "Input FASTA/FASTQ file path (use '-' for stdin)");
  app.add_option("-t,--threads", num_threads, "Number of worker threads")
      ->default_val(num_threads);
  app.add_option("-m,--mode", mode_str, "Scanner mode: syncmer | minimizer")
      ->transform(CLI::Transformer(
          std::map<std::string, std::string>{{"syncmer", "syncmer"},
                                             {"minimizer", "minimizer"}},
          CLI::ignore_case))
      ->default_val("syncmer");
  app.add_option("-r,--reader", reader_type,
                 "Reader type: uring | mmap | pipe")
      ->transform(CLI::Transformer(
          std::map<std::string, std::string>{{"uring", "uring"},
                                             {"mmap", "mmap"},
                                             {"pipe", "pipe"}},
          CLI::ignore_case))
      ->default_val("mmap");
  app.add_option(
         "-B,--chunk-batch-mib", chunk_batch_mib,
         "Reader chunk size and distributor raw batch (MiB, base 1024). "
         "Decimals allowed; converted to bytes with rounding, then aligned up "
         "to 4KiB (matches uring O_DIRECT constraints).")
      ->default_val(0.25);
  app.add_flag("--fastq", force_fastq, "Force FASTQ format parsing");
  app.add_flag("--profile", profile_breakdown,
               "Print wall-time breakdown (main thread vs worker sums)");

  CLI11_PARSE(app, argc, argv);

  const size_t io_bytes = chunk_batch_bytes_from_mib(chunk_batch_mib);

  // Auto-detect stdin if path is empty and not a terminal
  if (path.empty()) {
    if (!isatty(STDIN_FILENO)) {
      path = "-";
    } else {
      std::cout << app.help() << std::endl;
      return 1;
    }
  }

  bool is_pipe = (path == "-" || path == "/dev/stdin");
  if (!is_pipe) {
    struct stat st;
    if (stat(path.c_str(), &st) == 0 && S_ISFIFO(st.st_mode))
      is_pipe = true;
  }

  // Auto-switch to pipe reader if input is a pipe
  if (is_pipe && reader_type == "uring") {
    reader_type = "pipe";
  }

  XsealScannerType scanner_type = (mode_str == "minimizer")
                                      ? XsealScannerType::MINIMIZER
                                      : XsealScannerType::SYNCMER;

  // Detect FASTQ by extension or force flag
  bool is_fastq = force_fastq;
  if (!is_fastq && !is_pipe) {
    std::string lower_path = path;
    std::transform(lower_path.begin(), lower_path.end(), lower_path.begin(),
                   ::tolower);
    is_fastq = lower_path.ends_with(".fastq") || lower_path.ends_with(".fq") ||
               lower_path.ends_with(".fastq.gz") ||
               lower_path.ends_with(".fq.gz");
  }

  if (is_pipe || reader_type == "pipe") {
    int fd = (path == "-" || path == "/dev/stdin")
                 ? STDIN_FILENO
                 : open(path.c_str(), O_RDONLY);
    XsealPipeReader reader(fd, io_bytes);
    dispatch_parser(reader, is_fastq, scanner_type, num_threads, reader_type,
                    force_fastq, profile_breakdown, io_bytes);
  } else if (reader_type == "mmap") {
    XsealMmapReader reader(path, io_bytes);
    dispatch_parser(reader, is_fastq, scanner_type, num_threads, reader_type,
                    force_fastq, profile_breakdown, io_bytes);
  } else {
    XsealUringReader reader(path, io_bytes);
    dispatch_parser(reader, is_fastq, scanner_type, num_threads, reader_type,
                    force_fastq, profile_breakdown, io_bytes);
  }

  return 0;
}
