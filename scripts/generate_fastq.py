# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2026 hanrai. All Rights Reserved.

import random
import sys

bases = ['A', 'C', 'G', 'T']
num_reads = 10000000  # 10 Million

with open('data/illumina_storm.fastq', 'w') as f:
    for i in range(num_reads):
        f.write(f"@SIMULATED_SEQ:{i} length=150\n")
        # generate simple repeating patterns to make it faster
        seq = "ATGC" * 37 + "AT"
        f.write(seq + '\n')
        f.write("+\n")
        f.write('I' * 150 + '\n')
