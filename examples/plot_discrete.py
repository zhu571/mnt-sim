#!/usr/bin/env python
"""Paper-style Fig. 2 with discrete color bands (no continuous gradient)."""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm, LogNorm
from matplotlib.cm import ScalarMappable
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from mnt_sim.cross_section.imqmd_like import ImQMDModel

print("Calculating...")
model = ImQMDModel("U", "U", 7.0, Ap=238, At=238)
result = model.calculate(dz_range=(-8, 8), dn_range=(-10, 10), n_l=300)
theta, E, d2 = model.angular_energy_dist(result, n_theta=90, nE=200)
d2d = d2[0, 0, :, :]

# Sample scatter points
np.random.seed(42)
n_p = 50000
prob = d2d.flatten()
prob = prob / prob.sum()
idx = np.random.choice(len(prob), size=n_p, p=prob, replace=True)
t_i = idx // len(E)
e_i = idx % len(E)
dT, dE = theta[1]-theta[0], E[1]-E[0]
EK = E[e_i] + (np.random.rand(n_p)-0.5)*dE
theta_scat = theta[t_i] + (np.random.rand(n_p)-0.5)*dT
sigma = d2d[t_i, e_i]

mask = (EK>=0)&(EK<=1800)&(theta_scat>=0)&(theta_scat<=90)
theta_scat, EK, sigma = theta_scat[mask], EK[mask], sigma[mask]
sigma = np.clip(sigma, 1e-7, 1e-2)

print(f"Points: {len(theta_scat)}")

# === Discrete color boundaries (paper style) ===
boundaries = [1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2]
ncolors = len(boundaries) - 1

# Colors matching paper description:
# 10^-2→orange/red, 10^-3→yellow, 10^-4→yellow-green,
# 10^-5→green, 10^-6→cyan, 10^-7→blue
colors = ['#00008B', '#00CED1', '#228B22', '#9ACD32', '#FFD700', '#FF4500']

cmap = matplotlib.colors.ListedColormap(colors)
cmap.set_under('white')  # values < 1e-7 → white (no color)
norm = BoundaryNorm(boundaries, cmap.N, clip=False, extend='min')

# === Plot ===
fig, ax = plt.subplots(figsize=(8, 6.5))

# 1. hist2d with discrete colors
counts, xedges, yedges, im = ax.hist2d(
    theta_scat, EK, bins=[80, 160], weights=sigma,
    cmap=cmap, norm=norm
)

# 2. Contour lines from hist2d counts
levels = boundaries
X, Y = np.meshgrid((xedges[:-1]+xedges[1:])/2, (yedges[:-1]+yedges[1:])/2)
cs = ax.contour(X, Y, counts.T, levels=levels,
                colors='black', linewidths=0.6)
ax.clabel(cs, inline=True, fontsize=6, fmt='%.0e', colors='black')

# 3. Scatter overlay (same discrete colors)
ax.scatter(theta_scat, EK, c=sigma, norm=norm, cmap=cmap,
           s=6, marker='D', edgecolors='none', alpha=0.6, rasterized=True)

# Axes
ax.set_xlim(0, 90)
ax.set_ylim(0, 1800)
ax.set_xticks(np.arange(0, 91, 10))
ax.set_yticks(np.arange(0, 1801, 200))
ax.tick_params(labelsize=10)
ax.set_xlabel(r'$\theta_{\rm lab}$ (deg.)', fontsize=13)
ax.set_ylabel(r'$E_K$ (MeV)', fontsize=13)
ax.set_facecolor('white')
ax.grid(False)

# 4. Colorbar inset
cax = inset_axes(ax, width="15%", height="40%", loc='upper right',
                  bbox_to_anchor=(0.05, 0.05, 1, 1),
                  bbox_transform=ax.transAxes)
# Create a colorbar with discrete colors
mappable = ScalarMappable(norm=norm, cmap=cmap)
mappable.set_array([])
cb = fig.colorbar(mappable, cax=cax, ticks=boundaries)
cb.set_label(r'$\sigma$ (mb/MeV/sr)', fontsize=9, labelpad=3)
cax.tick_params(labelsize=7)
cb.set_ticks(boundaries)
labels = [r'$10^{-7}$', r'$10^{-6}$', r'$10^{-5}$',
           r'$10^{-4}$', r'$10^{-3}$', r'$10^{-2}$']
cb.set_ticklabels(labels)

# 5. Annotations
ax.text(55, 1550, r'$^{243}_{\;\;92}\mathrm{U}$',
        fontsize=14, fontweight='bold',
        bbox=dict(boxstyle='round,pad=0.2', facecolor='white',
                 edgecolor='gray', linewidth=0.5, alpha=0.8))
ax.text(72, 1650, r'$\sigma$ (mb/MeV/sr)',
        fontsize=9, fontweight='normal', color='black')

od = os.path.join(os.path.dirname(__file__), '..', 'output', 'figures')
os.makedirs(od, exist_ok=True)
p = os.path.join(od, 'Fig2_discrete.png')
plt.savefig(p, dpi=300, facecolor='white', edgecolor='none')
plt.close()
print(f"Saved: {p}")
print("Done!")
