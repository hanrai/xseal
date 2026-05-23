#!/bin/bash

# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2026 hanrai. All Rights Reserved.

set -e

BIN_FILE="data/hg38_2bit.bin"
RESULTS_FILE="experiments/rigorous_bench/generality_results.txt"

echo "Compiling..."
clang++ -O3 -DNDEBUG -mavx2 -mbmi2 -mfma -ffast-math -march=native -std=c++20 -pthread -Iinclude experiments/rigorous_bench/xseal_bench.cpp -o experiments/rigorous_bench/xseal_bench
(cd experiments/rigorous_bench/simd_min_bench && RUSTFLAGS="-C target-cpu=native" cargo build --release)

echo "ParamSet|Tool|Mode|Throughput" > $RESULTS_FILE

# Set A: Small
echo "Running Set A (k=15, s=7, w=11)..."
./experiments/rigorous_bench/xseal_bench $BIN_FILE 1 15 7 11 | grep "RESULT" | awk -F'|' '{print "Small|" $2 "|" $3 "|" $5}' >> $RESULTS_FILE
./experiments/rigorous_bench/simd_min_bench/target/release/simd_min_bench $BIN_FILE 1 15 7 11 | grep "RESULT" | awk -F'|' '{print "Small|" $2 "|" $3 "|" $5}' >> $RESULTS_FILE

# Set B: Medium
echo "Running Set B (k=21, s=11, w=21)..."
./experiments/rigorous_bench/xseal_bench $BIN_FILE 1 21 11 21 | grep "RESULT" | awk -F'|' '{print "Medium|" $2 "|" $3 "|" $5}' >> $RESULTS_FILE
./experiments/rigorous_bench/simd_min_bench/target/release/simd_min_bench $BIN_FILE 1 21 11 21 | grep "RESULT" | awk -F'|' '{print "Medium|" $2 "|" $3 "|" $5}' >> $RESULTS_FILE

# Set C: Large
echo "Running Set C (k=31, s=15, w=51)..."
./experiments/rigorous_bench/xseal_bench $BIN_FILE 1 31 15 51 | grep "RESULT" | awk -F'|' '{print "Large|" $2 "|" $3 "|" $5}' >> $RESULTS_FILE
./experiments/rigorous_bench/simd_min_bench/target/release/simd_min_bench $BIN_FILE 1 31 15 51 | grep "RESULT" | awk -F'|' '{print "Large|" $2 "|" $3 "|" $5}' >> $RESULTS_FILE

echo "Generality test complete."
column -t -s '|' $RESULTS_FILE
