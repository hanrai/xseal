# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2026 hanrai. All Rights Reserved.

import matplotlib.pyplot as plt
import numpy as np

labels = ['XSeal', 'simd-minimizers']
simd_pct = [71.55, 57.56]
scalar_pct = [28.45, 42.44]
ipc = [2.80, 3.37]

x = np.arange(len(labels))
width = 0.35

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

# Plot 1: Instruction Mix (Stacked)
ax1.bar(labels, simd_pct, width, label='SIMD (SSE/AVX)', color='#1f77b4', alpha=0.8)
ax1.bar(labels, scalar_pct, width, bottom=simd_pct, label='Scalar/Other', color='#ff7f0e', alpha=0.8)
ax1.set_ylabel('Instruction Mix (%)', fontsize=12, fontweight='bold')
ax1.set_title('Vectorization Density', fontsize=14, fontweight='bold')
ax1.legend(loc='upper right')
ax1.set_ylim(0, 115)

# Plot 2: IPC Comparison
ax2.bar(labels, ipc, width, color=['#1f77b4', '#d62728'], alpha=0.8, edgecolor='black')
ax2.set_ylabel('IPC (Instructions Per Cycle)', fontsize=12, fontweight='bold')
ax2.set_title('Instruction Per Cycle (IPC)', fontsize=14, fontweight='bold')
ax2.set_ylim(0, 4.5)

# Add IPC values on bars
for i, v in enumerate(ipc):
    ax2.text(i, v + 0.1, f'{v:.2f}', ha='center', fontsize=12, fontweight='bold')

plt.tight_layout()
plt.savefig('paper/figures/fig5_comparison.pdf', bbox_inches='tight')
print("Figure 5 (Comparative Mix) generated.")
