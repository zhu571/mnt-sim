#!/usr/bin/env python
"""Final style: use the proven continuous model with tuned parameters."""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm
from matplotlib.cm import ScalarMappable
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from mnt_sim.cross_section.imqmd_like import ImQMDModel

print("Calculating...")
model = ImQMDModel("U", "U", 7.0, Ap=238, At=238)
result = model.calculate(dz_range=(-8, 8), dn_range=(-10, 10), n_l=300)
theta, E, d2 = model.angular_energy_dist(result, n_theta=90, nE=200)
d2d = d2[0, 0, :, :]

# Sample scatter
np.random.seed(42)
n_p = 50000
prob = d2d.flatten(); prob = prob / prob.sum()
idx = np.random.choice(len(prob), size=n_p, p=prob, replace=True)
t_i = idx // len(E); e_i = idx % len(E)
dT, dE = theta[1]-theta[0], E[1]-E[0]
t_s = theta[t_i] + (np.random.rand(n_p)-0.5)*dT
e_s = E[e_i] + (np.random.rand(n_p)-0.5)*dE
s_s = d2d[t_i, e_i]
m = (e_s>=0)&(e_s<=1800)&(t_s>=0)&(t_s<=90)
t_s, e_s, s_s = t_s[m], e_s[m], np.clip(s_s[m], 1e-7, 1e-2)

# Discrete colormap
boundaries = [1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2]
colors = ['#00008B', '#00CED1', '#228B22', '#9ACD32', '#FFD700', '#FF4500']
cmap = matplotlib.colors.ListedColormap(colors)
cmap.set_under('white')
norm = BoundaryNorm(boundaries, cmap.N, clip=False, extend='min')

fig, ax = plt.subplots(figsize=(8, 6.5))

counts, xedges, yedges, im = ax.hist2d(
    t_s, e_s, bins=[80, 160], weights=s_s,
    cmap=cmap, norm=norm)

X = (xedges[:-1]+xedges[1:])/2; Y = (yedges[:-1]+yedges[1:])/2
Xg, Yg = np.meshgrid(X, Y)
cs = ax.contour(Xg, Yg, counts.T, levels=boundaries,
                colors='black', linewidths=0.6)
ax.clabel(cs, inline=True, fontsize=6, fmt='%.0e', colors='black')

ax.scatter(t_s, e_s, c=s_s, norm=norm, cmap=cmap,
           s=6, marker='D', edgecolors='none', alpha=0.5, rasterized=True)

ax.set_xlim(0, 90); ax.set_ylim(0, 1800)
ax.set_xticks(np.arange(0, 91, 10))
ax.set_yticks(np.arange(0, 1801, 200))
ax.tick_params(labelsize=10)
ax.set_xlabel(r'$\theta_{\rm lab}$ (deg.)', fontsize=13)
ax.set_ylabel(r'$E_K$ (MeV)', fontsize=13)
ax.set_facecolor('white'); ax.grid(False)

cax = inset_axes(ax, width="15%", height="40%", loc='upper right',
                  bbox_to_anchor=(0.05, 0.05, 1, 1),
                  bbox_transform=ax.transAxes)
mappable = ScalarMappable(norm=norm, cmap=cmap); mappable.set_array([])
cb = fig.colorbar(mappable, cax=cax, ticks=boundaries)
cb.set_label(r'$\sigma$ (mb/MeV/sr)', fontsize=9, labelpad=3)
cax.tick_params(labelsize=7)
cb.set_ticks(boundaries)
cb.set_ticklabels([r'$10^{-7}$',r'$10^{-6}$',r'$10^{-5}$',r'$10^{-4}$',r'$10^{-3}$',r'$10^{-2}$'])

ax.text(55, 1550, r'$^{243}_{\;\;92}\mathrm{U}$', fontsize=14, fontweight='bold',
        bbox=dict(boxstyle='round,pad=0.2', facecolor='white', edgecolor='gray', linewidth=0.5, alpha=0.8))
ax.text(72, 1650, r'$\sigma$ (mb/MeV/sr)', fontsize=9, color='black')

od = os.path.join(os.path.dirname(__file__), '..', 'output', 'figures')
os.makedirs(od, exist_ok=True)
p = os.path.join(od, 'Fig2_continuous.png')
plt.savefig(p, dpi=300, facecolor='white', edgecolor='none')
plt.close()
print(f"Saved: {p}")

mi = np.unravel_index(np.argmax(d2d), d2d.shape)
print(f"Peak: theta={theta[mi[0]]:.1f}, E={E[mi[1]]:.0f}")
print("Distribution matrix:")
for ti in range(0, 90, 9):
    line = f'{theta[ti]:5.1f} '
    for ei in range(0, 160, 20):
        v = d2d[ti, ei]
        if v>1e-2: line+=f' {v:5.2f}'
        elif v>1e-3: line+=f' {v:5.3f}'
        elif v>1e-4: line+=f' {v:5.4f}'
        elif v>1e-5: line+=f' {v:5.5f}'
        elif v>1e-6: line+=f' {v:5.6f}'
        elif v>1e-7: line+=f' {v:5.7f}'
        else: line+='  ....'
    print(line)
