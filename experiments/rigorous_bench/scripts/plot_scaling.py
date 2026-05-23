# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2026 hanrai. All Rights Reserved.

import matplotlib.pyplot as plt
import numpy as np

# Balanced sequence for 6-core/12-thread CPU: 1, 2, 4, 6, 12, 16
threads = [1, 2, 4, 6, 12, 16]

# XSeal (Clang)
x_open = [1.36, 2.71, 5.35, 7.85, 8.41, 8.04]
x_min = [1.17, 2.33, 4.52, 7.30, 7.30, 7.32]
x_closed = [1.08, 2.15, 4.19, 6.69, 6.69, 6.52]

# Rust (simd-minimizers)
r_open = [0.71, 1.36, 2.29, 3.12, 3.25, 4.04]
r_min = [0.59, 1.15, 1.89, 2.65, 2.58, 3.09]
r_closed = [0.56, 1.06, 1.75, 2.45, 2.79, 2.91]

fig, ax = plt.subplots(figsize=(10, 6))

ax.plot(threads, x_open, marker='s', label='XSeal (Open Sync)', color='#1f77b4', linewidth=2, markersize=8)
ax.plot(threads, x_min, marker='o', label='XSeal (Minimizer)', color='#ff7f0e', linewidth=2, markersize=8)
ax.plot(threads, x_closed, marker='^', label='XSeal (Closed Sync)', color='#2ca02c', linewidth=2, markersize=8)

ax.plot(threads, r_open, marker='s', label='simd-minimizers (Open Sync)', color='#1f77b4', linestyle='--', alpha=0.5)
ax.plot(threads, r_min, marker='o', label='simd-minimizers (Minimizer)', color='#ff7f0e', linestyle='--', alpha=0.5)
ax.plot(threads, r_closed, marker='^', label='simd-minimizers (Closed Sync)', color='#2ca02c', linestyle='--', alpha=0.5)

ax.set_xlabel('Thread Count', fontsize=12, fontweight='bold')
ax.set_ylabel('Aggregate Throughput (Gbp/s)', fontsize=12, fontweight='bold')
ax.set_title('Figure 3: Throughput Scaling: XSeal vs. simd-minimizers', fontsize=14, fontweight='bold')
ax.set_xticks([1, 2, 4, 6, 12, 16])
ax.grid(True, linestyle='--', alpha=0.5)
ax.legend(fontsize=10, ncol=2, loc='lower right')

plt.tight_layout()
plt.savefig('paper/figures/fig3_scaling.pdf', bbox_inches='tight')
print("Figure 3 updated: removed incorrect memory wall label.")
