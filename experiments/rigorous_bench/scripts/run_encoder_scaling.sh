#!/bin/bash

# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2026 hanrai. All Rights Reserved.

THREADS=(1 6 12)
MODES=("minimap2" "xseal")
FREQS=(0 80 10)
SIZE=1024

echo "Running Encoder Scaling Benchmarks..."
for f in "${FREQS[@]}"; do
    for t in "${THREADS[@]}"; do
        for m in "${MODES[@]}"; do
            ./experiments/rigorous_bench/encoder_scaling_bench $SIZE $t $m $f
        done
    done
done
