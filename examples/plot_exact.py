#!/usr/bin/env python
"""
Exact paper-style Fig. 2 with strict plotting parameters.
"""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from mnt_sim.cross_section.imqmd_like import ImQMDModel

print("Calculating...")
model = ImQMDModel("U", "U", 7.0, Ap=238, At=238)
result = model.calculate(dz_range=(-8, 8), dn_range=(-10, 10), n_l=300)
theta, E, d2 = model.angular_energy_dist(result, n_theta=90, nE=200)
d2d = d2[0, 0, :, :]

# Sample scatter points from distribution
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

# Low-energy tail removed - kinematics naturally produces back-folding events
# at ~30-50°, E≈0 from θ_cm≈180° for heavy fragments (v_cm_C < v_cm)

theta_all, EK_all, sigma_all = theta_scat, EK, sigma

mask = (EK_all>=0) & (EK_all<=1800) & (theta_all>=0) & (theta_all<=90)
theta_all, EK_all, sigma_all = theta_all[mask], EK_all[mask], sigma_all[mask]
sigma_all = np.clip(sigma_all, 1e-7, 1e-2)

print(f"Points: {len(theta_all)}")
print(f"Peak: theta={theta[np.unravel_index(np.argmax(d2d),d2d.shape)[0]]:.1f}, "
      f"E={E[np.unravel_index(np.argmax(d2d),d2d.shape)[1]]:.0f}")

# === Plot ===
fig, ax = plt.subplots(figsize=(8, 6.5))

# 1. hist2d as base layer
counts, xedges, yedges, im = ax.hist2d(
    theta_all, EK_all, bins=[80, 160], weights=sigma_all,
    cmap='jet', norm=LogNorm(vmin=1e-7, vmax=1e-2)
)

# 2. Contour from histogram
levels = np.logspace(-7, -2, 8)
X, Y = np.meshgrid((xedges[:-1]+xedges[1:])/2, (yedges[:-1]+yedges[1:])/2)
cs = ax.contour(X, Y, counts.T, levels=levels,
                colors='black', linewidths=0.6)
ax.clabel(cs, inline=True, fontsize=6, fmt='%.0e', colors='black')

# 3. Scatter overlay
ax.scatter(theta_all, EK_all, c=sigma_all,
           norm=LogNorm(vmin=1e-7, vmax=1e-2),
           cmap='jet', s=6, marker='D', edgecolors='none',
           alpha=0.6, rasterized=True)

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

# 4. Colorbar inset in upper right (15% × 40%)
cax = inset_axes(ax, width="15%", height="40%", loc='upper right',
                  bbox_to_anchor=(0.05, 0.05, 1, 1),
                  bbox_transform=ax.transAxes)
cb = fig.colorbar(im, cax=cax)
cb.set_label(r'$\sigma$ (mb/MeV/sr)', fontsize=9, labelpad=3)
cax.tick_params(labelsize=7)
cb.set_ticks([1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2])
cb.set_ticklabels([r'$10^{-7}$', r'$10^{-6}$', r'$10^{-5}$',
                    r'$10^{-4}$', r'$10^{-3}$', r'$10^{-2}$'])

# 5. Annotations
ax.text(55, 1550, r'$^{243}_{\;\;92}\mathrm{U}$',
        fontsize=14, fontweight='bold',
        bbox=dict(boxstyle='round,pad=0.2', facecolor='white',
                 edgecolor='gray', linewidth=0.5, alpha=0.8))
ax.text(72, 1650, r'$\sigma$ (mb/MeV/sr)',
        fontsize=9, fontweight='normal', color='black')

# Tight layout not used (inset_axes incompatible)
# plt.tight_layout()

od = os.path.join(os.path.dirname(__file__), '..', 'output', 'figures')
os.makedirs(od, exist_ok=True)
p = os.path.join(od, 'Fig2_exact.png')
plt.savefig(p, dpi=300, facecolor='white', edgecolor='none')
plt.close()
print(f"Saved: {p}")
print("Done!")
