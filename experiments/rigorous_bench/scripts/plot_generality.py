# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2026 hanrai. All Rights Reserved.

import matplotlib.pyplot as plt
import numpy as np

sets = ['Small', 'Medium', 'Large']
modes = ['Syncmer', 'OpenSync', 'Minimizer']

# Data from results
xseal = {
    'Small': [0.8654, 1.5645, 0.9542],
    'Medium': [0.7136, 1.6776, 0.9881],
    'Large': [0.8454, 1.4565, 0.7715]
}
rust = {
    'Small': [0.4557, 0.6233, 0.5307],
    'Medium': [0.5013, 0.6512, 0.6103],
    'Large': [0.5589, 0.7005, 0.6866]
}

fig, ax = plt.subplots(figsize=(10, 6))

x = np.arange(len(sets))
width = 0.25

# Combine modes into grouped bars
for i, mode in enumerate(modes):
    xseal_vals = [xseal[s][i] for s in sets]
    rust_vals = [rust[s][i] for s in sets]
    
    pos = x + (i - 1) * width
    ax.bar(pos, xseal_vals, width/2, label=f'XSeal ({mode})', color=['#1f77b4', '#2ca02c', '#d62728'][i], alpha=0.8)
    ax.bar(pos + width/2, rust_vals, width/2, label=f'simd-min ({mode})', hatch='//', color=['#1f77b4', '#2ca02c', '#d62728'][i], alpha=0.5)

ax.set_ylabel('Throughput (Gbp/s) - Single Thread', fontsize=12, fontweight='bold')
ax.set_title('Figure 6: Performance Generalization Across Parameter Sets', fontsize=14, fontweight='bold')
ax.set_xticks(x)
ax.set_xticklabels(['Small\n(15,7,11)', 'Medium\n(21,11,21)', 'Large\n(31,15,51)'])
ax.legend(ncol=3, loc='upper center', bbox_to_anchor=(0.5, -0.15))

plt.tight_layout()
plt.savefig('paper/figures/fig6_generality.pdf', bbox_inches='tight')
print("Figure 6 (Generality) generated.")
