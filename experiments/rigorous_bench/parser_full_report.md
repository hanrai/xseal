# XSeal Parser Micro-architectural Performance Report

This report documents the performance of the XSeal DFA-based FASTA parser under two distinct data scenarios: **Real Genomic Data (hg38)** and **Random Mock Data**.

## Executive Summary
The XSeal Parser demonstrates near-peak micro-architectural efficiency on real biological sequences, achieving **20.57 GiB/s (0.19 cycles/base)** on a single thread. The performance gap between real data and random data (11.23 GiB/s) is primarily attributed to **State Transition Density** and **Branch Predictor Pressure**.

---

## 1. Quantitative Comparison

| Metric | Real hg38 (data/hg38.fa) | Random Mock (512MB) | Gap / Note |
| :--- | :--- | :--- | :--- |
| **Throughput (GiB/s)** | **20.57** | 11.23 | **1.83x Advantage** |
| **Cycles/base** | **0.190** | 0.341 | Lower is better |
| **State Transitions** | **3,346** | **670,896** | 200x Difference in density |
| **Instructions Per Cycle (IPC)** | **3.40** | 1.21 | Near-peak for AVX2 vs logic-bound |
| **Branch Miss Rate** | **0.32%** | **2.76%** | High branch pressure in random |

---

## 2. Micro-architectural Analysis

### 2.1 The "State Transition" Penalty
The parser uses a DFA (Deterministic Finite Automaton) to identify segments (Bases, Gaps, Headers). 
- In **hg38**, segments are extremely long (chromosomes). The parser stays in the **Vectorized Fast Path** for millions of cycles without interruption.
- In **Random Data**, we injected segments every 800 bytes. This forces the CPU to frequently exit the vectorized loop, enter scalar transition logic, and update metadata (`push_segment`). This "jitter" prevents the CPU from achieving its maximum IPC.

### 2.2 Branch Prediction and Pipelining
- **hg38**: The almost non-existent branch miss rate (0.32%) allows the out-of-order execution engine to look far ahead and hide memory latency.
- **Random Data**: The higher branch miss rate (2.76%) causes frequent pipeline flushes. Every time `testz` detects a stop/transition, the speculative execution must be re-evaluated, incurring a heavy penalty.

---

## 3. Fused Frontend Impact (hg38)
When encoding is integrated (Fused mode), XSeal still maintains high efficiency.

| Configuration | Throughput (GiB/s) | Cycles/base |
| :--- | :--- | :--- |
| **XSeal Parser Only** | 20.57 | 0.190 |
| **XSeal Fused (Parse+Encode)** | 5.75 | 0.680 |
| **kseq.h (Parse Only Baseline)** | 1.71 | 2.286 |

**Insight**: XSeal Fused is **3.4x faster** than the baseline parser *alone*, confirming that our streaming co-design effectively eliminates the "Integration Wall" overhead.

---
**Test Environment**: Ubuntu 24.04, Clang 18.1.3, CPU with AVX2 support.
**Reproducibility**: Run `bash rigorous_bench/profile_parser.sh` to regenerate this data.
