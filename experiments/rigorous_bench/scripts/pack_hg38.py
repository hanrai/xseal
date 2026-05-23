# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2026 hanrai. All Rights Reserved.

import sys
import os

def pack_fasta(input_path, output_path):
    # Mapping: A=0, C=1, G=2, T=3. Map N and others to A.
    mapping = {
        'A': 0, 'C': 1, 'G': 2, 'T': 3,
        'a': 0, 'c': 1, 'g': 2, 't': 3
    }
    
    with open(input_path, 'r') as f_in, open(output_path, 'wb') as f_out:
        bit_buffer = 0
        bit_count = 0
        total_bases = 0
        
        for line in f_in:
            if line.startswith('>'):
                continue
            for char in line.strip():
                val = mapping.get(char, 0)
                bit_buffer |= (val << bit_count)
                bit_count += 2
                total_bases += 1
                
                if bit_count == 64:
                    f_out.write(bit_buffer.to_bytes(8, 'little'))
                    bit_buffer = 0
                    bit_count = 0
        
        if bit_count > 0:
            f_out.write(bit_buffer.to_bytes(8, 'little'))
            
    print(f"Packed {total_bases} bases to {output_path}")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python pack_hg38.py <input.fa> <output.bin>")
    else:
        pack_fasta(sys.argv[1], sys.argv[2])
