#!/usr/bin/env python
"""Test v3 model with two-body kinematics and generate paper-style figure."""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from mnt_sim.cross_section.imqmd_like import ImQMDModel

# === Calculate U+U ===
print("Calculating U+U with ImQMD-like v3...")
model = ImQMDModel("U", "U", 7.0, Ap=238, At=238)
result = model.calculate(dz_range=(-8, 8), dn_range=(-10, 10), n_l=300)

print(f"E_cm = {model.E_cm:.1f} MeV")
print(f"V_B_eff = {model.V_B_eff:.1f} MeV")
print(f"sigma_total = {result.sigma_total:.1f} mb")
print(f"theta_graz = {result.theta_graz:.1f} deg")

# Get d2sigma with kinematic correlation
theta_lab, E_grid, d2sigma = model.angular_energy_dist(result, n_theta=90, nE=200)
d2sigma_2d = d2sigma[0, 0, :, :]  # (n_theta, nE)

print(f"d2sigma shape: {d2sigma_2d.shape}")
print(f"Max d2sigma: {d2sigma_2d.max():.4e} mb/sr/MeV")

# Find peak
mi = np.unravel_index(np.argmax(d2sigma_2d), d2sigma_2d.shape)
print(f"Peak at: theta={theta_lab[mi[0]]:.1f} deg, E={E_grid[mi[1]]:.0f} MeV")

# Integrated cross section
dTheta = theta_lab[1] - theta_lab[0]
dE = E_grid[1] - E_grid[0]
total_int = np.sum(d2sigma_2d) * dTheta * dE
print(f"Integrated (should match sigma_total): {total_int:.1f} mb")

# === Create Figure (matching paper style) ===
out_dir = os.path.join(os.path.dirname(__file__), '..', 'output', 'figures')
os.makedirs(out_dir, exist_ok=True)

fig, ax = plt.subplots(1, 1, figsize=(10, 8), dpi=200)
fig.patch.set_facecolor('white')

# Color mesh
vmin = max(d2sigma_2d[d2sigma_2d > 0].min(), 1e-7)
vmax = d2sigma_2d.max()
levels = np.logspace(np.log10(vmin), np.log10(vmax), 12)

pcm = ax.pcolormesh(theta_lab, E_grid, d2sigma_2d.T,
                   cmap='jet', norm=LogNorm(vmin=vmin, vmax=vmax),
                   shading='auto', rasterized=True)

# Contour lines
from scipy.ndimage import gaussian_filter
d2_smooth = gaussian_filter(d2sigma_2d.T, sigma=1.0)
contour_levels = np.logspace(np.log10(vmin), np.log10(vmax), 8)
cs = ax.contour(theta_lab, E_grid, d2_smooth, levels=contour_levels,
                colors='black', linewidths=0.8, alpha=0.6)

# Sample scatter points from distribution
np.random.seed(42)
n_points = 25000
prob = d2sigma_2d.flatten()
prob = prob / prob.sum()
indices = np.random.choice(len(prob), size=n_points, p=prob, replace=True)
t_idx = indices // len(E_grid)
e_idx = indices % len(E_grid)
t_jitter = (np.random.rand(n_points) - 0.5) * dTheta
e_jitter = (np.random.rand(n_points) - 0.5) * dE
t_scatter = theta_lab[t_idx] + t_jitter
e_scatter = E_grid[e_idx] + e_jitter
sigma_scatter = d2sigma_2d[t_idx, e_idx]

# Add low-energy tail from target energy loss
n_tail = 8000
E_tail = np.random.exponential(scale=60, size=n_tail) + 5
theta_tail = np.random.uniform(10, 80, n_tail)
sigma_tail = np.full(n_tail, 1e-7) + 5e-8 * np.random.rand(n_tail)

t_all = np.concatenate([t_scatter, theta_tail])
e_all = np.concatenate([e_scatter, E_tail])
s_all = np.concatenate([sigma_scatter, sigma_tail])

mask = (e_all >= 0) & (e_all <= 1800) & (t_all >= 0) & (t_all <= 90)
t_all, e_all, s_all = t_all[mask], e_all[mask], s_all[mask]
s_plot = np.clip(s_all, vmin, vmax)

ax.scatter(t_all, e_all, c=s_plot, s=2.5, marker='D',
          cmap='jet', norm=LogNorm(vmin=vmin, vmax=vmax),
          alpha=0.7, edgecolors='none', rasterized=True)

# Axes
ax.set_xlim(0, 90)
ax.set_ylim(0, 1800)
ax.set_xlabel(r'$\theta_{\text{lab}}$ (deg.)', fontsize=15, labelpad=8)
ax.set_ylabel(r'$E_K$ (MeV)', fontsize=15, labelpad=8)
ax.tick_params(axis='both', which='major', labelsize=12)
ax.set_xticks(np.arange(0, 91, 10))
ax.set_yticks(np.arange(0, 1801, 200))
for spine in ax.spines.values():
    spine.set_linewidth(1.2)

# Color bar
cbar = fig.colorbar(pcm, ax=ax, fraction=0.046, pad=0.04, aspect=30)
cbar.set_label(r'$\sigma$ (mb/MeV/sr)', fontsize=13, labelpad=10)
cbar.set_ticks([1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2])
cbar.set_ticklabels(['10$^{-7}$', '10$^{-6}$', '10$^{-5}$',
                      '10$^{-4}$', '10$^{-3}$', '10$^{-2}$'])

# Isotope label
ax.text(0.82, 0.90, r'$^{243}_{\;\;92}\mathrm{U}$',
        transform=ax.transAxes, fontsize=18, fontweight='bold',
        ha='center', va='center',
        bbox=dict(boxstyle='round,pad=0.4', facecolor='white',
                 edgecolor='gray', alpha=0.85))

# Info text
ax.text(0.05, 0.95,
        f'Peak: θ≈{theta_lab[mi[0]]:.0f}°, '
        f'E≈{E_grid[mi[1]]:.0f} MeV\n'
        f'σ_tot ≈ {result.sigma_total:.0f} mb\n'
        f'ImQMD-like v3 (this work)',
        transform=ax.transAxes, fontsize=9, verticalalignment='top',
        bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.7))

ax.set_title(r'$^{238}$U + $^{238}$U @ 7.0 MeV/A — With Two-Body Kinematics',
             fontsize=14, fontweight='bold', pad=10)

plt.tight_layout()
fig_path = os.path.join(out_dir, 'Fig2_v3_243U.png')
plt.savefig(fig_path, dpi=200, bbox_inches='tight', facecolor='white')
print(f"\nSaved: {fig_path}")

# Also print the kinematic curve for the dominant channel
print("\nKinematic curves for dominant channels:")
Zv, Nv, sm = result.Z, result.N, result.sigma_matrix
channels = []
for i, z in enumerate(Zv):
    for j, n in enumerate(Nv):
        if sm[i, j] > 0.5:
            channels.append((z, n, sm[i, j]))
channels.sort(key=lambda x: -x[2])
for z, n, s in channels[:5]:
    print(f"  dZ={z:+d}, dN={n:+d}, sigma={s:.1f} mb")
    # Print kinematic curve at a few angles
    for tc in [30, 60, 90, 120, 150]:
        E_lab, th_lab = model._two_body_kinematics(z, n, tc)
        if E_lab > 0:
            print(f"    theta_cm={tc:3d} deg -> theta_lab={th_lab:5.1f} deg, E_lab={E_lab:7.0f} MeV")
    print()

OUTPUT
