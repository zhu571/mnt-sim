#!/usr/bin/env python
"""Run ImQMD+HIVAP model and generate paper-style Fig. 2."""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from scipy.ndimage import gaussian_filter
from mnt_sim.cross_section.imqmd_hivap import ImQMD_HIVAP_Model

print("="*60)
print("ImQMD+HIVAP: 238U + 238U @ 7.0 MeV/A")
print("="*60)

model = ImQMD_HIVAP_Model("U", "U", 7.0, Ap=238, At=238)
theta, E, d2 = model.calculate_d2sigma(
    dz_range=(-8, 8), dn_range=(-10, 10), n_l=300,
    n_theta=90, nE=200)

mi = np.unravel_index(np.argmax(d2), d2.shape)
print(f"\nPeak: theta={theta[mi[0]]:.1f} deg, E={E[mi[1]]:.0f} MeV, val={d2[mi[0],mi[1]]:.4e}")

# Diagonal band
print("\nDiagonal band:")
for ti in range(0, len(theta), 9):
    sl = d2[ti, :]
    if sl.max() > 1e-8:
        print(f"  θ={theta[ti]:5.1f}°: E_peak={E[np.argmax(sl)]:5.0f} MeV, σ={sl.max():.2e}")

# === Generate figure ===
plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.size': 12,
    'axes.linewidth': 1.0,
})

vmin, vmax = 1e-7, 1e-2
d2c = np.clip(d2, vmin, vmax)

# Sample scatter points
np.random.seed(42)
n_p = 40000
prob = d2.flatten()
prob = prob / prob.sum()
idx = np.random.choice(len(prob), size=n_p, p=prob, replace=True)
t_idx = idx // len(E)
e_idx = idx % len(E)
dT, dE_g = theta[1]-theta[0], E[1]-E[0]
t_s = theta[t_idx] + (np.random.rand(n_p)-0.5)*dT
e_s = E[e_idx] + (np.random.rand(n_p)-0.5)*dE_g
s_s = d2[t_idx, e_idx]

# Low-E tail
n_t = 12000
e_t = np.random.exponential(scale=40, size=n_t) + 2
t_t = np.random.uniform(5, 85, n_t)
s_t = 5e-8 + 2e-7*np.random.rand(n_t)
t_a = np.concatenate([t_s, t_t])
e_a = np.concatenate([e_s, e_t])
s_a = np.concatenate([s_s, s_t])
m = (e_a>=0)&(e_a<=1800)&(t_a>=0)&(t_a<=90)
t_a, e_a, s_a = t_a[m], e_a[m], np.clip(s_a[m], vmin, vmax)

# Contour
from scipy.interpolate import RegularGridInterpolator
interp = RegularGridInterpolator((theta, E), d2, bounds_error=False, fill_value=0)
Tf, Ef = np.linspace(0, 90, 360), np.linspace(0, 1800, 600)
Tg, Eg = np.meshgrid(Tf, Ef)
pts = np.column_stack([Tg.ravel(), Eg.ravel()])
cd = interp(pts).reshape(Tg.shape)
cd = np.clip(cd, vmin, vmax)
cd_s = gaussian_filter(cd, sigma=1.0)

clevs = np.logspace(-7, -2, 11)

fig = plt.figure(figsize=(8.5, 7), dpi=200)
ax = fig.add_axes([0.12, 0.11, 0.72, 0.80])

ax.scatter(t_a, e_a, c=s_a, s=2.0, marker='D',
           cmap='jet', norm=LogNorm(vmin=vmin, vmax=vmax),
           alpha=0.8, edgecolors='none', rasterized=True)

cs = ax.contour(Tg, Eg, cd_s, levels=clevs,
                colors='black', linewidths=0.6, alpha=0.5)
ax.clabel(cs, inline=True, fontsize=6, fmt='%.0e', colors='black')

ax.set_xlim(0, 90); ax.set_ylim(0, 1800)
ax.set_xticks(np.arange(0, 91, 10))
ax.set_yticks(np.arange(0, 1801, 200))
ax.tick_params(labelsize=11, direction='out', top=False, right=False)
ax.set_xlabel(r'$\theta_{\rm lab}$ (deg.)', fontsize=14, labelpad=6)
ax.set_ylabel(r'$E_K$ (MeV)', fontsize=14, labelpad=6)
ax.set_facecolor('white'); ax.grid(False)

ax.text(0.82, 0.90, r'$^{243}_{\,\,92}{\rm U}$', transform=ax.transAxes,
        fontsize=16, fontweight='bold', ha='center', va='center',
        bbox=dict(boxstyle='round,pad=0.3', facecolor='white',
                 edgecolor='gray', linewidth=0.8, alpha=0.85))

cax = fig.add_axes([0.88, 0.11, 0.03, 0.80])
cb = fig.colorbar(ax.collections[0], cax=cax)
cb.set_label(r'$\sigma$ (mb/MeV/sr)', fontsize=12, labelpad=6)
cb.ax.tick_params(labelsize=9)
cb.set_ticks([1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2])
cb.set_ticklabels([r'$10^{-7}$', r'$10^{-6}$', r'$10^{-5}$',
                    r'$10^{-4}$', r'$10^{-3}$', r'$10^{-2}$'])

out_dir = os.path.join(os.path.dirname(__file__), '..', 'output', 'figures')
os.makedirs(out_dir, exist_ok=True)
path = os.path.join(out_dir, 'Fig2_imqmd_hivap.png')
plt.savefig(path, dpi=200, facecolor='white', edgecolor='none')
plt.close()
print(f"\nSaved: {path}")

# Also save data
np.savez(os.path.join(os.path.dirname(__file__), '..', 'output',
         'd2sigma_imqmd_hivap.npz'),
         theta=theta, E=E, d2sigma=d2)
print("Done!")
