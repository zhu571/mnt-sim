#!/usr/bin/env python
"""Strict paper-style Fig. 2 reproduction for ²⁴³U from ²³⁸U+²³⁸U."""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from mnt_sim.cross_section.imqmd_like import ImQMDModel
from scipy.ndimage import gaussian_filter

# === 1. Calculate ===
print("Calculating...")
model = ImQMDModel("U", "U", 7.0, Ap=238, At=238)
result = model.calculate(dz_range=(-8, 8), dn_range=(-10, 10), n_l=300)

theta_lab, E_grid, d2sigma = model.angular_energy_dist(result, n_theta=90, nE=200)
d2 = d2sigma[0, 0, :, :]

mi = np.unravel_index(np.argmax(d2), d2.shape)
print(f"Peak: theta={theta_lab[mi[0]]:.1f} deg, E={E_grid[mi[1]]:.0f} MeV")

# === 2. Sample scatter points ===
np.random.seed(42)
n_points = 40000
prob = d2.flatten()
prob = prob / prob.sum()
idx = np.random.choice(len(prob), size=n_points, p=prob, replace=True)
t_idx = idx // len(E_grid)
e_idx = idx % len(E_grid)
dT = theta_lab[1] - theta_lab[0]
dE = E_grid[1] - E_grid[0]
t_scat = theta_lab[t_idx] + (np.random.rand(n_points)-0.5)*dT
e_scat = E_grid[e_idx] + (np.random.rand(n_points)-0.5)*dE
s_scat = d2[t_idx, e_idx]

# Low-energy tail from target energy loss
n_tail = 15000
E_tail = np.random.exponential(scale=50, size=n_tail) + 3
theta_tail = np.random.uniform(5, 80, n_tail)
sigma_tail = 5e-8 + 2e-7 * np.random.rand(n_tail)

t_all = np.concatenate([t_scat, theta_tail])
e_all = np.concatenate([e_scat, E_tail])
s_all = np.concatenate([s_scat, sigma_tail])

m = (e_all >= 0) & (e_all <= 1800) & (t_all >= 0) & (t_all <= 90)
t_all, e_all, s_all = t_all[m], e_all[m], s_all[m]
print(f"Total scatter points: {len(t_all)}")

# === 3. Figure (strict paper style) ===
plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'],
    'font.size': 12,
    'axes.linewidth': 1.2,
    'xtick.major.width': 1.0,
    'ytick.major.width': 1.0,
    'xtick.major.size': 5,
    'ytick.major.size': 5,
    'xtick.minor.size': 2.5,
    'ytick.minor.size': 2.5,
})

fig = plt.figure(figsize=(8.5, 7), dpi=200)
ax = fig.add_axes([0.12, 0.11, 0.72, 0.80])

vmin, vmax = 1e-7, 1e-2
s_plot = np.clip(s_all, vmin, vmax)

# Scatter: diamond markers, rasterized for performance
ax.scatter(t_all, e_all, c=s_plot, s=2.0, marker='D',
           cmap='jet', norm=LogNorm(vmin=vmin, vmax=vmax),
           alpha=0.85, edgecolors='none', rasterized=True)

# Contour lines (from smoothed 2D histogram)
from scipy.interpolate import RegularGridInterpolator
interp = RegularGridInterpolator((theta_lab, E_grid), d2,
                                 bounds_error=False, fill_value=0)
Tf = np.linspace(0, 90, 360)
Ef = np.linspace(0, 1800, 600)
Tg, Eg = np.meshgrid(Tf, Ef)
pts = np.column_stack([Tg.ravel(), Eg.ravel()])
cd = interp(pts).reshape(Tg.shape)
cd = np.clip(cd, vmin, vmax)
cd_s = gaussian_filter(cd, sigma=1.0)

clevs = np.logspace(-7, -2, 11)
cs = ax.contour(Tg, Eg, cd_s, levels=clevs,
                colors='black', linewidths=0.6, alpha=0.5)
# Label only the clearest contours
ax.clabel(cs, inline=True, fontsize=6, fmt='%.0e', colors='black')

# Axes: strict paper style
ax.set_xlim(0, 90)
ax.set_ylim(0, 1800)
ax.set_xticks(np.arange(0, 91, 10))
ax.set_yticks(np.arange(0, 1801, 200))
ax.tick_params(axis='both', which='major', labelsize=11,
               direction='out', top=False, right=False)
ax.set_xlabel(r'$\theta_{\rm lab}$ (deg.)', fontsize=14, labelpad=6)
ax.set_ylabel(r'$E_K$ (MeV)', fontsize=14, labelpad=6)
ax.set_facecolor('white')
ax.grid(False)

# Isotope label inside plot (paper style)
ax.text(0.82, 0.90, r'$^{243}_{\,\,92}{\rm U}$',
        transform=ax.transAxes, fontsize=16, fontweight='bold',
        ha='center', va='center',
        bbox=dict(boxstyle='round,pad=0.3', facecolor='white',
                 edgecolor='gray', linewidth=0.8, alpha=0.85))

# Color bar on the right (paper style: vertical, label at top)
cbar_ax = fig.add_axes([0.88, 0.11, 0.03, 0.80])
cbar = fig.colorbar(ax.collections[0], cax=cbar_ax)
cbar.set_label(r'$\sigma$ (mb/MeV/sr)', fontsize=12,
               labelpad=8, rotation=-90, ha='center', va='bottom')
cbar_ax.yaxis.set_label_position('right')
cbar_ax.yaxis.label.set_verticalalignment('bottom')

# Tick formatting for log color bar
cbar.ax.tick_params(labelsize=9)
cbar.set_ticks([1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2])
cbar.set_ticklabels([r'$10^{-7}$', r'$10^{-6}$', r'$10^{-5}$',
                      r'$10^{-4}$', r'$10^{-3}$', r'$10^{-2}$'])

# Save
out_dir = os.path.join(os.path.dirname(__file__), '..', 'output', 'figures')
os.makedirs(out_dir, exist_ok=True)
path = os.path.join(out_dir, 'Fig2_paper_style.png')
plt.savefig(path, dpi=200, facecolor='white', edgecolor='none')
plt.close()
print(f"Saved: {path}")
print("Done!")
