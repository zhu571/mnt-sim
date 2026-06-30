import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sys; sys.path.insert(0, '.')
from mnt_sim.cross_section.imqmd_hivap import ImQMD_HIVAP_Model

m = ImQMD_HIVAP_Model('U', 'U', 7.0, Ap=238, At=238)
t, E, d = m.calculate_d2sigma(dz_range=(-8,8), dn_range=(-10,10), n_l=300, n_theta=90, nE=200)

fig, axes = plt.subplots(1, 3, figsize=(18, 5.5))

# Plot 1: d2sigma contour
d_plot = np.clip(d.T, 1e-5, None)
im = axes[0].pcolormesh(t, E, d_plot, shading='auto', cmap='hot_r',
    norm=matplotlib.colors.LogNorm(vmin=1e-4, vmax=d.max()))
plt.colorbar(im, ax=axes[0], label='d2sigma [mb/deg/MeV]', shrink=0.8)
idx = np.unravel_index(np.argmax(d), d.shape)
axes[0].plot(t[idx[0]], E[idx[1]], 'c*', ms=15, label=f'peak {t[idx[0]]:.1f}deg, {E[idx[1]]:.0f}MeV')
axes[0].legend(fontsize=9)
axes[0].set_xlabel('theta_lab [deg]')
axes[0].set_ylabel('E_lab [MeV]')
axes[0].set_title('238U+238U @ 7.0 MeV/u')
axes[0].set_xlim(0, 90)
axes[0].set_ylim(0, 1800)

# Plot 2: angle-integrated energy spectrum per angular bin
colors = ['#e41a1c','#377eb8','#4daf4a','#984ea3']
for i, (th_lo, th_hi) in enumerate([(0,20),(20,40),(40,60),(60,90)]):
    mask = (t >= th_lo) & (t <= th_hi)
    spec = d[mask, :].sum(axis=0)
    axes[1].plot(E, spec, color=colors[i], lw=1.5, label=f'{th_lo}-{th_hi} deg')
axes[1].set_xlabel('E_lab [MeV]')
axes[1].set_ylabel('dsigma/dE [mb/MeV]')
axes[1].set_title('Energy spectra by angle bin')
axes[1].legend(fontsize=8)
axes[1].set_xlim(0, 1800)

# Plot 3: theta distribution (E-integrated)
axes[2].plot(t, d.sum(axis=1), 'b-', lw=1.5)
axes[2].axvline(t[idx[0]], color='r', ls='--', alpha=0.7, label=f'peak {t[idx[0]]:.1f}deg')
axes[2].set_xlabel('theta_lab [deg]')
axes[2].set_ylabel('dsigma/dtheta [mb/deg]')
axes[2].set_title('Angular distribution')
axes[2].legend(fontsize=8)
axes[2].set_xlim(0, 90)

plt.tight_layout()
plt.savefig('output/figures/d2sigma_new_imqmd_hivap.png', dpi=150)
print('Saved: output/figures/d2sigma_new_imqmd_hivap.png')
