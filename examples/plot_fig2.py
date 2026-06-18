#!/usr/bin/env python
"""
Create Fig. 2 style double-differential cross-section plots.
Compares our calculation with Wang et al. NIM B 463 (2020).
"""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import csv

def load_csv(path):
    """Load d2sigma CSV: theta rows, energy columns."""
    with open(path) as f:
        reader = csv.reader(f)
        header = next(reader)
        E_vals = np.array([float(h) for h in header[1:]])
        theta = []
        data = []
        for row in reader:
            theta.append(float(row[0]))
            data.append([float(v) for v in row[1:]])
    return np.array(theta), E_vals, np.array(data)

out_dir = os.path.join(os.path.dirname(__file__), '..', 'output')
fig_dir = os.path.join(out_dir, 'figures')
os.makedirs(fig_dir, exist_ok=True)

# Panel 1: Xe+Pt → 200Os
theta1, E1, d21 = load_csv(os.path.join(out_dir, 'd2_200Os.csv'))
# Panel 2: U+U → 243U
theta2, E2, d22 = load_csv(os.path.join(out_dir, 'd2_243U.csv'))

# Create figure with TWO panels (matching Fig. 2 layout)
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 7), dpi=150)

# Common styling
plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.size': 12,
    'axes.labelsize': 13,
    'axes.titlesize': 14,
    'xtick.labelsize': 11,
    'ytick.labelsize': 11,
})

# ===== Panel 1: 200Os from Xe+Pt =====
# The data is in d²σ/dE/dΩ (mb/sr/MeV). Values are very small (10^-24 to 10^-3)
# Use log scale for color
d21_plot = d21.T  # transpose to (E, theta) for pcolormesh

# Normalize for better visualization: multiply by scale factor
# Paper's data is at ~150 MeV, 60° - our peak is at 644 MeV, 71°
# Let's show the full distribution
vmin1 = max(d21_plot[d21_plot > 0].min(), 1e-20)
vmax1 = d21_plot.max()

pcm1 = ax1.pcolormesh(E1, theta1, d21,
                      norm=LogNorm(vmin=vmin1, vmax=vmax1),
                      cmap='jet', shading='auto')
cb1 = fig.colorbar(pcm1, ax=ax1, label=r'd²$\sigma$/dE/d$\Omega$ [mb/sr/MeV]')
ax1.set_xlabel('Energy [MeV]')
ax1.set_ylabel('Angle [deg]')
ax1.set_title(r'²⁰⁰Os from ¹³⁶Xe+¹⁹⁸Pt @ 7.98 MeV/A (Our Model)')
ax1.set_xlim(0, 1200)
ax1.set_ylim(0, 90)

# Mark paper's peak region
ax1.add_patch(plt.Rectangle((100, 45), 100, 30,
                            fill=False, edgecolor='white', linestyle='--',
                            linewidth=2, label='Paper: ~60°, ~150 MeV'))
ax1.text(150, 42, 'Paper\nθ≈60°\nE≈150 MeV', color='white',
         fontsize=9, ha='center', fontweight='bold',
         bbox=dict(boxstyle='round,pad=0.3', facecolor='black', alpha=0.6))

# Mark our peak
peak_i = np.unravel_index(np.argmax(d21), d21.shape)
ax1.plot(E1[peak_i[1]], theta1[peak_i[0]], 'k*', markersize=15,
         markeredgecolor='white', markeredgewidth=1.5, label=f'Our peak')
ax1.annotate(f'Our peak\nθ={theta1[peak_i[0]]:.1f}°, E={E1[peak_i[1]]:.0f} MeV',
             xy=(E1[peak_i[1]], theta1[peak_i[0]]),
             xytext=(E1[peak_i[1]]+80, theta1[peak_i[0]]-8),
             color='black', fontsize=9, fontweight='bold',
             arrowprops=dict(arrowstyle='->', color='black', lw=1.5),
             bbox=dict(boxstyle='round,pad=0.3', facecolor='yellow', alpha=0.8))
ax1.legend(loc='upper right', fontsize=9)

# ===== Panel 2: 243U from U+U =====
d22_plot = d22.T
vmin2 = max(d22_plot[d22_plot > 0].min(), 1e-20)
vmax2 = d22_plot.max()

pcm2 = ax2.pcolormesh(E2, theta2, d22,
                      norm=LogNorm(vmin=vmin2, vmax=vmax2),
                      cmap='jet', shading='auto')
cb2 = fig.colorbar(pcm2, ax=ax2, label=r'd²$\sigma$/dE/d$\Omega$ [mb/sr/MeV]')
ax2.set_xlabel('Energy [MeV]')
ax2.set_ylabel('Angle [deg]')
ax2.set_title(r'²⁴³U from ²³⁸U+²³⁸U @ 7.0 MeV/A (Our Model)')
ax2.set_xlim(0, 1200)
ax2.set_ylim(0, 90)

# Mark paper's peak region  
ax2.add_patch(plt.Rectangle((600, 30), 400, 20,
                            fill=False, edgecolor='white', linestyle='--',
                            linewidth=2, label='Paper: θ=35-50°, E=600-1000 MeV'))
ax2.text(800, 28, 'Paper region\nθ=35-50°\nE=600-1000 MeV', color='white',
         fontsize=9, ha='center', fontweight='bold',
         bbox=dict(boxstyle='round,pad=0.3', facecolor='black', alpha=0.6))

# Mark our peak
peak_j = np.unravel_index(np.argmax(d22), d22.shape)
ax2.plot(E2[peak_j[1]], theta2[peak_j[0]], 'k*', markersize=15,
         markeredgecolor='white', markeredgewidth=1.5, label=f'Our peak')
ax2.annotate(f'Our peak\nθ={theta2[peak_j[0]]:.1f}°, E={E2[peak_j[1]]:.0f} MeV',
             xy=(E2[peak_j[1]], theta2[peak_j[0]]),
             xytext=(E2[peak_j[1]]-200, theta2[peak_j[0]]-10),
             color='black', fontsize=9, fontweight='bold',
             arrowprops=dict(arrowstyle='->', color='black', lw=1.5),
             bbox=dict(boxstyle='round,pad=0.3', facecolor='yellow', alpha=0.8))
ax2.legend(loc='upper right', fontsize=9)

plt.suptitle('Double Differential Cross Sections (d²σ/dE/dΩ) — Comparison with Wang et al.',
             fontsize=15, fontweight='bold', y=1.02)
plt.tight_layout()
fig_path = os.path.join(fig_dir, 'Fig2_comparison.png')
plt.savefig(fig_path, dpi=200, bbox_inches='tight')
print(f'Saved: {fig_path}')

# Also create a combined comparison figure showing just the d2sigma distribution
# with better visualization
fig2, (ax_a, ax_b) = plt.subplots(1, 2, figsize=(16, 6), dpi=150)

# Xe+Pt panel with contour lines
levels1 = np.logspace(np.log10(vmin1), np.log10(vmax1), 15)
cf1 = ax_a.contourf(E1, theta1, d21, levels=levels1,
                    cmap='viridis', norm=LogNorm())
fig2.colorbar(cf1, ax=ax_a, label='d²σ/dE/dΩ [mb/sr/MeV]')
ax_a.contour(E1, theta1, d21, levels=levels1[::3],
             colors='white', linewidths=0.5, alpha=0.5)
ax_a.set_xlabel('Energy [MeV]')
ax_a.set_ylabel('Angle [deg]')
ax_a.set_title(r'²⁰⁰Os from ¹³⁶Xe+¹⁹⁸Pt (Our Model)', fontsize=13)
ax_a.set_xlim(0, 1200)
ax_a.set_ylim(0, 90)
ax_a.text(0.05, 0.95, f'Peak: θ={theta1[peak_i[0]]:.0f}°, E={E1[peak_i[1]]:.0f} MeV',
          transform=ax_a.transAxes, color='white', fontsize=10, fontweight='bold',
          verticalalignment='top',
          bbox=dict(boxstyle='round', facecolor='black', alpha=0.6))

# U+U panel with contour lines  
levels2 = np.logspace(np.log10(vmin2), np.log10(vmax2), 15)
cf2 = ax_b.contourf(E2, theta2, d22, levels=levels2,
                    cmap='viridis', norm=LogNorm())
fig2.colorbar(cf2, ax=ax_b, label='d²σ/dE/dΩ [mb/sr/MeV]')
ax_b.contour(E2, theta2, d22, levels=levels2[::3],
             colors='white', linewidths=0.5, alpha=0.5)
ax_b.set_xlabel('Energy [MeV]')
ax_b.set_ylabel('Angle [deg]')
ax_b.set_title(r'²⁴³U from ²³⁸U+²³⁸U (Our Model)', fontsize=13)
ax_b.set_xlim(0, 1200)
ax_b.set_ylim(0, 90)
ax_b.text(0.05, 0.95, f'Peak: θ={theta2[peak_j[0]]:.0f}°, E={E2[peak_j[1]]:.0f} MeV',
          transform=ax_b.transAxes, color='white', fontsize=10, fontweight='bold',
          verticalalignment='top',
          bbox=dict(boxstyle='round', facecolor='black', alpha=0.6))

plt.suptitle('Double Differential Cross Sections — ImQMD-like Model',
             fontsize=15, fontweight='bold')
plt.tight_layout()
fig2_path = os.path.join(fig_dir, 'd2sigma_contour.png')
plt.savefig(fig2_path, dpi=200, bbox_inches='tight')
print(f'Saved: {fig2_path}')

print('\n✅ Figures created!')
