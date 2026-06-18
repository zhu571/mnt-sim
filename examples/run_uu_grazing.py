#!/usr/bin/env python
"""
Run GRAZING model for 238U + 238U @ 7.0 MeV/A
and compare with Wang et al. NIM B 463 (2020) Fig. 2
"""
import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

import numpy as np
import yaml

from mnt_sim.cross_section.grazing import GrazingModel

def main():
    config_path = os.path.join(os.path.dirname(__file__), '..', 'config', 'uu.yaml')
    with open(config_path) as f:
        cfg = yaml.safe_load(f)

    rxn = cfg['reaction']
    print("=" * 70)
    print(f"GRAZING Model: {rxn['projectile_A']}{rxn['projectile']} + "
          f"{rxn['target_A']}{rxn['target']} @ {rxn['E_lab']} MeV/u")
    print("=" * 70)

    grazing = GrazingModel(
        projectile=rxn['projectile'],
        target=rxn['target'],
        E_lab=rxn['E_lab']
    )

    print(f"Projectile: Z={grazing.Zp}, A={grazing.Ap}")
    print(f"Target:      Z={grazing.Zt}, A={grazing.At}")
    print(f"E_cm = {grazing.E_cm:.1f} MeV")
    print(f"Coulomb barrier V_C = {grazing.V_C:.1f} MeV")
    print(f"E_cm - V_C = {grazing.E_cm - grazing.V_C:.1f} MeV")

    # Calculate cross sections with wider range
    result = grazing.calculate(
        delta_Z_range=(-8, 8),
        delta_N_range=(-10, 10),
        n_l=200
    )

    print(f"\nGrazing angle: {result.theta_graz:.1f}°")
    Z_vals, N_vals = result.Z, result.N
    sm = result.sigma_matrix

    # Show top channels manually
    channels = []
    for i, z in enumerate(Z_vals):
        for j, n in enumerate(N_vals):
            val = sm[i, j]
            if val > 1e-4:
                channels.append({'Z': z, 'N': n, 'A': z+n+476, 'sigma_mb': val})

    channels.sort(key=lambda x: -x['sigma_mb'])
    print("\nTop transfer channels (by cross section):")
    print(f"  {'ΔZ':>4} {'ΔN':>4} {'A_CN':>7}  {'σ (mb)':>10}")
    print(f"  {'-'*30}")
    total_cs = 0
    for ch in channels[:15]:
        print(f"  {ch['Z']:+4d} {ch['N']:+4d} {ch['A']:7d}  {ch['sigma_mb']:10.4f}")
        total_cs += ch['sigma_mb']
    print(f"\n  Total (sum of above): {total_cs:.1f} mb")

    # All channels sum
    total_all = np.sum(sm)
    print(f"  Total (all channels):  {total_all:.1f} mb")

    # Look for 243U channels
    # 243U = Z=92, N=151. For U+U symmetric: projectile A=238, target A=238
    # PLF with ΔZ=0, ΔN=+5 → Z=92+0=92, N=146+5=151 → 243U
    print(f"\n--- Channels for ²⁴³U (A=243, Z=92, N=151) ---")
    print(f"  ²⁴³U = projectile (238U) + 5n: ΔZ=0, ΔN=+5")
    print(f"  ²⁴³U = target (238U) + 5n: ΔZ=0, ΔN=-5 (from projectile perspective)")
    
    for target_dz, target_dn, label in [(0, 5, 'PLF ²⁴³U'), (0, -5, 'TLF ²⁴³U')]:
        zi = np.where(Z_vals == target_dz)[0]
        ni = np.where(N_vals == target_dn)[0]
        if len(zi) > 0 and len(ni) > 0:
            sigma_val = sm[zi[0], ni[0]]
            print(f"  {label}: σ = {sigma_val:.4f} mb (ΔZ={target_dz:+d}, ΔN={target_dn:+d})")
        else:
            print(f"  {label}: not in grid range")

    # --- Build double differential cross section ---
    nE = 50
    E_min, E_max = 200, 1400
    E_vals = np.linspace(E_min, E_max, nE)
    theta = np.linspace(0, 90, 30)  # degrees

    # Grazing angle from classical formula for near-barrier
    theta_graz = result.theta_graz  # degrees

    nZ_n, nN_n = len(Z_vals), len(N_vals)
    d2sigma = np.zeros((nZ_n, nN_n, len(theta), nE))

    for i in range(nZ_n):
        for j in range(nN_n):
            if sm[i, j] < 1e-10:
                continue
            sigma_theta_ch = 8.0 + 0.2 * abs(Z_vals[i]) + 0.2 * abs(N_vals[j])
            ang_factor = np.exp(-(theta - theta_graz - Z_vals[i] * 0.3)**2
                                / (2 * sigma_theta_ch**2))
            ang_factor /= np.sum(ang_factor) + 1e-30

            Q_approx = Z_vals[i] * 6.0 + N_vals[j] * 8.0
            E_mean = result.E_cm - Q_approx * 0.5
            E_mean = max(E_min, min(E_mean, E_max))
            sigma_E = 50 + 5 * abs(Z_vals[i]) + 5 * abs(N_vals[j])

            for k in range(len(theta)):
                d2sigma[i, j, k, :] = (sm[i, j] * ang_factor[k]
                                       * np.exp(-(E_vals - E_mean)**2 / (2 * sigma_E**2))
                                       / (np.sqrt(2 * np.pi) * sigma_E))
                total = np.sum(d2sigma[i, j, k, :]) + 1e-30
                if total > 0:
                    d2sigma[i, j, k, :] /= total

    d2sigma_total = np.sum(d2sigma, axis=(0, 1))
    max_idx = np.unravel_index(np.argmax(d2sigma_total), d2sigma_total.shape)
    print(f"\n--- Double Differential Cross Section ---")
    print(f"  Peak: θ ≈ {theta[max_idx[0]]:.1f}°, E ≈ {E_vals[max_idx[1]]:.0f} MeV")

    # Print matrix
    print(f"\n  d²σ/dE/dΩ matrix (all channels):")
    print(f"  θ↓ \\ E→  ", end="")
    for e_idx in range(0, nE, 10):
        print(f" {E_vals[e_idx]:6.0f}", end="")
    print()
    for k in range(0, len(theta), 3):
        t = theta[k]
        line = f"  {t:5.1f}°   "
        for e_idx in range(0, nE, 10):
            val = d2sigma_total[k, e_idx]
            if val > 1e-6:
                line += f" {val:6.4f}"
            else:
                line += f"       "
        print(line)

    # Save data
    out_dir = os.path.join(os.path.dirname(__file__), '..', 'output', 'uu')
    os.makedirs(out_dir, exist_ok=True)

    np.savez(os.path.join(out_dir, 'd2sigma_grazing.npz'),
             theta=theta, E=E_vals, d2sigma=d2sigma_total,
             Z_vals=Z_vals, N_vals=N_vals, sigma_matrix=sm)

    import csv
    with open(os.path.join(out_dir, 'd2sigma_grazing.csv'), 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['theta_deg'] + [f'{e:.0f}' for e in E_vals])
        for k, t in enumerate(theta):
            w.writerow([f'{t:.1f}'] + [f'{d2sigma_total[k, e_idx]:.6e}' for e_idx in range(nE)])

    # Save channel cross sections
    with open(os.path.join(out_dir, 'channels_grazing.csv'), 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['delta_Z', 'delta_N', 'sigma_mb'])
        for i, z in enumerate(Z_vals):
            for j, n in enumerate(N_vals):
                if sm[i, j] > 1e-4:
                    w.writerow([z, n, sm[i, j]])

    print(f"\n  Data saved to: {out_dir}")

    # === Comparison with paper ===
    print("\n" + "=" * 70)
    print("COMPARISON WITH Wang et al. NIM B 463 (2020)")
    print("=" * 70)
    print(f"  System: ²³⁸U + ²³⁸U @ 7.0 MeV/A")
    print(f"  {'':>5} {'Paper (ImQMD)':>18} {'Your (GRAZING)':>18}")
    print(f"  {'-'*45}")
    print(f"  {'Total σ':>12} {'(not given)':>18} {total_all:>10.1f} mb")
    print(f"  {'²⁴³U σ':>12} {'(not given)':>18} {'(see above)':>18}")
    print(f"  {'Angle peak':>12} {'35-50°':>18} {theta[max_idx[0]]:>10.1f}°")
    print(f"  {'Energy peak':>12} {'600-1000 MeV':>18} {E_vals[max_idx[1]]:>10.0f} MeV")
    print(f"\n  NOTE: Paper used ImQMD+GEMINI. GRAZING gives a simplified")
    print(f"  semi-classical estimate. For symmetric near-barrier U+U,")
    print(f"  the DNS model has insufficient driving force.")
    print(f"  The grazing angle of {theta_graz:.0f}° indicates trajectories")
    print(f"  are strongly bent at near-barrier energies.")

    print("\n✅ Done!")

if __name__ == '__main__':
    main()
