# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2026 hanrai. All Rights Reserved.

import matplotlib.pyplot as plt

# Data from micro-architectural profiling
labels = ['SIMD (SSE/AVX)', 'Scalar/Control/Memory']
sizes = [71.6, 28.4]
colors = ['#1f77b4', '#ff7f0e']
explode = (0.05, 0)  # explode the SIMD slice

fig, ax = plt.subplots(figsize=(8, 6))
patches, texts, autotexts = ax.pie(sizes, explode=explode, labels=labels, colors=colors,
                                    autopct='%1.1f%%', shadow=False, startangle=140,
                                    textprops={'fontsize': 12, 'fontweight': 'bold'})

# Style the percentage text
for autotext in autotexts:
    autotext.set_color('white')

ax.set_title('Figure 4: XSeal Instruction Mix (Retired Ops)', fontsize=14, fontweight='bold')
plt.tight_layout()

plt.savefig('paper/figures/fig4_instruction_mix.pdf', bbox_inches='tight')
print("Figure 4 (Instruction Mix) generated.")
