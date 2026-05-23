# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2026 hanrai. All Rights Reserved.

import matplotlib.pyplot as plt
import numpy as np

labels = ['Minimap2 (nt4)', 'XSeal (Prod)', 'XSeal (Pure)']
throughput = [1.80, 15.15, 16.60]
cycles = [2.17, 0.258, 0.236]

x = np.arange(len(labels))
width = 0.35

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

# Plot 1: Throughput (GiB/s)
bars1 = ax1.bar(labels, throughput, width, color=['#7f7f7f', '#1f77b4', '#aec7e8'], alpha=0.8, edgecolor='black')
ax1.set_ylabel('Throughput (GiB/s)', fontsize=12, fontweight='bold')
ax1.set_title('Front-end Throughput', fontsize=14, fontweight='bold')
ax1.set_ylim(0, 20)
ax1.grid(axis='y', linestyle='--', alpha=0.7)

# Add values on top of bars
for bar in bars1:
    height = bar.get_height()
    ax1.text(bar.get_x() + bar.get_width()/2., height + 0.5, f'{height:.2f}', ha='center', va='bottom', fontsize=11, fontweight='bold')

# Plot 2: Cycles/base
bars2 = ax2.bar(labels, cycles, width, color=['#7f7f7f', '#d62728', '#ff9896'], alpha=0.8, edgecolor='black')
ax2.set_ylabel('Cycles per Base', fontsize=12, fontweight='bold')
ax2.set_title('Micro-architectural Efficiency', fontsize=14, fontweight='bold')
ax2.set_ylim(0, 2.5)
ax2.grid(axis='y', linestyle='--', alpha=0.7)

# Add values on top of bars
for bar in bars2:
    height = bar.get_height()
    ax2.text(bar.get_x() + bar.get_width()/2., height + 0.05, f'{height:.3f}', ha='center', va='bottom', fontsize=11, fontweight='bold')

plt.tight_layout()
plt.savefig('paper/figures/fig8_encoder.pdf', bbox_inches='tight')
print("Figure 8 (Encoder Efficiency) generated.")
