# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2026 hanrai. All Rights Reserved.

import matplotlib.pyplot as plt
import numpy as np

# Data from latest balanced Clang run
labels = ['Minimizer', 'Syncmer (Closed)', 'Syncmer (Open)']

# XSeal (Clang)
xseal_means = [1.1711, 1.0822, 1.3576]
xseal_stds = [0.0036, 0.0026, 0.0024]

# simd-minimizers (Rust)
rust_means = [0.5915, 0.5598, 0.7093]
rust_stds = [0.0177, 0.0071, 0.0173]

x = np.arange(len(labels))
width = 0.35

fig, ax = plt.subplots(figsize=(10, 6))

rects1 = ax.bar(x - width/2, xseal_means, width, yerr=xseal_stds, 
                label='XSeal (Clang)', color='#1f77b4', capsize=5, edgecolor='black', alpha=0.9)
rects2 = ax.bar(x + width/2, rust_means, width, yerr=rust_stds, 
                label='simd-minimizers', color='#d62728', capsize=5, edgecolor='black', alpha=0.9)

ax.set_ylabel('Throughput (Gbp/s)', fontsize=12, fontweight='bold')
ax.set_title('Single-threaded Throughput: XSeal vs. simd-minimizers', fontsize=14, fontweight='bold')
ax.set_xticks(x)
ax.set_xticklabels(labels, fontsize=11, fontweight='bold')
ax.legend(fontsize=11)

ax.grid(axis='y', linestyle='--', alpha=0.7)
ax.set_ylim(0, 1.6)

def autolabel(rects):
    for rect in rects:
        height = rect.get_height()
        ax.annotate(f'{height:.2f}',
                    xy=(rect.get_x() + rect.get_width() / 2, height),
                    xytext=(0, 3),
                    textcoords="offset points",
                    ha='center', va='bottom', fontsize=10, fontweight='bold')

autolabel(rects1)
autolabel(rects2)

fig.tight_layout()
plt.savefig('paper/figures/fig1_throughput_comparison.pdf', bbox_inches='tight')
print("Figure 1 updated with balanced 3x3 comparison.")
