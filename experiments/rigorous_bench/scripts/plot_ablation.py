# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2026 hanrai. All Rights Reserved.

import matplotlib.pyplot as plt
import numpy as np

# Data from experiments.tex Table 2 (Ablation)
versions = ['V6.1', 'V7.0', 'V7.1', 'V7.4', 'V7.5', 'V10.0']
throughput = [0.67, 1.58, 0.75, 1.52, 1.09, 1.17]
labels = [
    'Baseline', 
    'VSB (Fwd-only)', 
    '+RC (Memory)', 
    'Reg-RC+ILP', 
    'Sparse-Table', 
    'Final (Clang)'
]

fig, ax = plt.subplots(figsize=(10, 5))

# Plot line
ax.plot(versions, throughput, marker='o', linestyle='-', color='#2ca02c', linewidth=2.5, markersize=8)

# Annotate points with labels
for i, txt in enumerate(labels):
    ax.annotate(txt, (versions[i], throughput[i]), xytext=(0, 10), 
                textcoords='offset points', ha='center', fontsize=9, fontweight='bold')

ax.set_ylabel('Throughput (Gbp/s)', fontsize=12, fontweight='bold')
ax.set_title('Figure 2: XSeal Architectural Evolution (Ablation Study)', fontsize=14, fontweight='bold')
ax.grid(True, linestyle='--', alpha=0.6)
ax.set_ylim(0, 1.8)

# Highlight the "Correctness Gap"
ax.annotate('The "Correctness Gap"', xy=('V7.1', 0.75), xytext=('V6.1', 0.2),
            arrowprops=dict(facecolor='black', shrink=0.05, width=1, headwidth=8),
            fontsize=10, color='red', fontweight='bold')

plt.tight_layout()
plt.savefig('paper/figures/fig2_ablation.pdf', bbox_inches='tight')
print("Figure 2 (Ablation) generated.")
