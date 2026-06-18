#!/usr/bin/env python
"""Run MC model and plot paper-style figure."""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm, LogNorm
from matplotlib.cm import ScalarMappable
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from mnt_sim.cross_section.imqmd_mc import ImQMD_MC_Model

print("Generating MC events...")
model = ImQMD_MC_Model()
events = model.generate_events(n_events=80000, target_mgcm2=95.0)

theta_all = events[:, 0]
E_all = events[:, 1]
sigma_all = events[:, 2]

print(f"Events: {len(theta_all)}")
print(f"Theta range: {theta_all.min():.1f}-{theta_all.max():.1f}")
print(f"E range: {E_all.min():.0f}-{E_all.max():.0f}")

# Create figure
boundaries = [1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2]
colors = ['#00008B', '#00CED1', '#228B22', '#9ACD32', '#FFD700', '#FF4500']
cmap = matplotlib.colors.ListedColormap(colors)
cmap.set_under('white')
norm = BoundaryNorm(boundaries, cmap.N, clip=False, extend='min')

fig, ax = plt.subplots(figsize=(8, 6.5))

# hist2d
counts, xedges, yedges, im = ax.hist2d(
    theta_all, E_all, bins=[80, 160], weights=sigma_all,
    cmap=cmap, norm=norm
)

# Contour
levels = boundaries
X = (xedges[:-1]+xedges[1:])/2
Y = (yedges[:-1]+yedges[1:])/2
Xg, Yg = np.meshgrid(X, Y)
cs = ax.contour(Xg, Yg, counts.T, levels=levels,
                colors='black', linewidths=0.6)
ax.clabel(cs, inline=True, fontsize=6, fmt='%.0e', colors='black')

# Scatter overlay
ax.scatter(theta_all, E_all, c=np.clip(sigma_all, 1e-7, 1e-2),
           norm=norm, cmap=cmap, s=4, marker='D',
           edgecolors='none', alpha=0.5, rasterized=True)

ax.set_xlim(0, 90)
ax.set_ylim(0, 1800)
ax.set_xticks(np.arange(0, 91, 10))
ax.set_yticks(np.arange(0, 1801, 200))
ax.tick_params(labelsize=10)
ax.set_xlabel(r'$\theta_{\rm lab}$ (deg.)', fontsize=13)
ax.set_ylabel(r'$E_K$ (MeV)', fontsize=13)
ax.set_facecolor('white')
ax.grid(False)

# Colorbar
cax = inset_axes(ax, width="15%", height="40%", loc='upper right',
                  bbox_to_anchor=(0.05, 0.05, 1, 1),
                  bbox_transform=ax.transAxes)
mappable = ScalarMappable(norm=norm, cmap=cmap)
mappable.set_array([])
cb = fig.colorbar(mappable, cax=cax, ticks=boundaries)
cb.set_label(r'$\sigma$ (mb/MeV/sr)', fontsize=9, labelpad=3)
cax.tick_params(labelsize=7)
cb.set_ticks(boundaries)
cb.set_ticklabels([r'$10^{-7}$', r'$10^{-6}$', r'$10^{-5}$',
                    r'$10^{-4}$', r'$10^{-3}$', r'$10^{-2}$'])

# Annotations
ax.text(55, 1550, r'$^{243}_{\;\;92}\mathrm{U}$',
        fontsize=14, fontweight='bold',
        bbox=dict(boxstyle='round,pad=0.2', facecolor='white',
                 edgecolor='gray', linewidth=0.5, alpha=0.8))
ax.text(72, 1650, r'$\sigma$ (mb/MeV/sr)',
        fontsize=9, fontweight='normal', color='black')

od = os.path.join(os.path.dirname(__file__), '..', 'output', 'figures')
os.makedirs(od, exist_ok=True)
p = os.path.join(od, 'Fig2_mc.png')
plt.savefig(p, dpi=300, facecolor='white', edgecolor='none')
plt.close()
print(f"Saved: {p}")

# Print stats
print(f"\n10^-2 region (sigma >= 0.01):")
for ti in range(28, 58, 3):
    t = X[ti]
    vals = [(Y[ei], counts.T[ti, ei]) for ei in range(len(Y)) if counts.T[ti, ei] >= 0.01]
    if vals:
        print(f"  th={t:5.1f}: E=({vals[0][0]:.0f}-{vals[-1][0]:.0f})")
