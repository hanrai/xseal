# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2026 hanrai. All Rights Reserved.

import xsealsyncmer
import os

def create_dummy_fasta(path):
    with open(path, "w") as f:
        f.write(">seq1\n")
        f.write("ACGT" * 1000 + "\n")
        f.write(">seq2\n")
        f.write("TGC" * 500 + "\n")

def test_bindings():
    fasta_path = "dummy.fasta"
    create_dummy_fasta(fasta_path)
    
    print(f"Extracting syncmers from {fasta_path}...")
    try:
        results = xsealsyncmer.extract_from_fasta(fasta_path, k=31, s=11)
        print(f"Extracted {len(results.positions)} syncmers.")
        
        if len(results.positions) > 0:
            print(f"First pos: {results.positions[0]}")
            print(f"First hash: {hex(results.hashes[0])}")
            print("SUCCESS: Bindings are working!")
        else:
            print("WARNING: No syncmers found (check k/s or sequence length).")
            
    finally:
        if os.path.exists(fasta_path):
            os.remove(fasta_path)

if __name__ == "__main__":
    test_bindings()
