# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2026 hanrai. All Rights Reserved.

import random
import sys

num_reads = 1000000 
with open('data/illumina_chaos.fastq', 'wb') as f:
    for i in range(num_reads):
        # randomly decide newline char
        nl = b'\r\n' if random.random() < 0.5 else b'\n'
        
        # inject BOM on 1% of reads before the @
        if random.random() < 0.01:
            f.write(b'\xEF\xBB\xBF')
            
        header = f"@CHAOS_SEQ:{i} length=150".encode()
        seq = (b"ATGC" * 37 + b"AT")
        plus = b"+"
        qual = b"I" * 150
        
        f.write(header + nl)
        f.write(seq + nl)
        f.write(plus + nl)
        f.write(qual + nl)
