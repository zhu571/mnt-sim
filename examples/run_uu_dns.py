#!/usr/bin/env python
"""
Run DNS model for 238U + 238U @ 7.0 MeV/A
and compare with Wang et al. NIM B 463 (2020) Fig. 2
"""
import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

import numpy as np
import yaml

from mnt_sim.cross_section.dns import DNSModel

def main():
    # Load U+U config
    config_path = os.path.join(os.path.dirname(__file__), '..', 'config', 'uu.yaml')
    with open(config_path) as f:
        cfg = yaml.safe_load(f)

    rxn = cfg['reaction']
    print("=" * 70)
    print(f"DNS Model: {rxn['projectile_A']}{rxn['projectile']} + "
          f"{rxn['target_A']}{rxn['target']} @ {rxn['E_lab']} MeV/u")
    print("=" * 70)

    # Initialize DNS model
    dns = DNSModel(
        projectile=rxn['projectile'],
        target=rxn['target'],
        E_lab=rxn['E_lab']
    )

    print(f"\nProjectile: Z={dns.Zp}, A={dns.Ap}")
    print(f"Target:      Z={dns.Zt}, A={dns.At}")
    print(f"E_cm = {dns.E_cm:.1f} MeV")
    print(f"Coulomb barrier V_C = {dns.V_C:.1f} MeV")
    print(f"E_cm - V_C = {dns.E_cm_eff:.1f} MeV")

    # Calculate cross sections
    cs = cfg['cross_section']
    result = dns.calculate(
        delta_Z_range=tuple(cs['delta_Z']),
        delta_N_range=tuple(cs['delta_N']),
        n_steps=cs['n_steps']
    )

    print(f"\nTotal reaction cross section: {result.sigma_total:.1f} mb")
    df = result.to_dataframe()

    # Show top channels
    top = df.nlargest(15, 'sigma_mb')
    print("\nTop transfer channels:")
    print(f"  {'ΔZ':>4} {'ΔN':>4} {'A':>5}  {'σ (mb)':>10}  {'Nuclide':>10}")
    print(f"  {'-'*40}")
    for _, row in top.iterrows():
        z = int(dns.Zp + row['Z'])
        n = int(row['A'] - z)
        a = int(row['A'])
        nuclide = f"{a}{''}"
        print(f"  {int(row['Z']):+4d} {int(row['N']):+4d} {a:5d}  {row['sigma_mb']:10.4f}")

    # Look for 243U specifically (Z=92, N=151, A=243)
    # In the DNS convention: ΔZ = Z_transfer_to_projectile, ΔN = N_transfer_to_projectile
    # For 243U as PLF: Z=92 = 92+0 → ΔZ=0, N=151 = 146+5 → ΔN=+5
    target_243u = df[(df['Z'] == 0) & (df['A'] == 243)]
    if len(target_243u) > 0:
        row = target_243u.iloc[0]
        print(f"\n--- ²⁴³U (ΔZ=0, ΔN=+5) ---")
        print(f"  Cross section: {row['sigma_mb']:.4f} mb")
    else:
        print(f"\n  ²⁴³U not found in direct transfer channels")
        # Check nearby
        nearby = df[(df['Z'] >= -2) & (df['Z'] <= 2) & (df['A'] >= 240) & (df['A'] <= 246)]
        if len(nearby) > 0:
            print("  Nearby A=240-246, |ΔZ|<=2 channels:")
            print(f"  {'ΔZ':>4} {'ΔN':>4} {'A':>5}  {'σ (mb)':>10}")
            for _, r in nearby.iterrows():
                print(f"  {int(r['Z']):+4d} {int(r['N']):+4d} {int(r['A']):5d}  {r['sigma_mb']:10.4f}")

    # --- Generate double differential cross sections ---
    print("\n--- Double Differential Cross Section (d²σ/dE/dΩ) ---")
    # Paper: 243U from U+U at 7.0 MeV/A
    # Outgoing angles: 35-50°, Kinetic energy: 600-1000 MeV

    # Energy range covering the paper's distribution
    E_min = 200.0
    E_max = 1400.0
    nE = 50  # Must match angular_distribution's nE
    E_vals = np.linspace(E_min, E_max, nE)

    # Angular range from 0 to 90 degrees (paper shows 35-50° for U+U)
    theta = np.linspace(0, 90, 30)

    d2sigma = dns.angular_distribution(theta, result)

    # For 243U (ΔZ=0, ΔN=+5), find the indices
    z_idx = np.where(result.Z == 0)[0]
    n_idx = np.where(result.N == 5)[0]

    print(f"\nAngular distribution for ²⁴³U (ΔZ=0, ΔN=+5):")
    print(f"  {'θ (deg)':>8} {'dσ/dΩ (mb/sr)':>16}")
    print(f"  {'-'*28}")
    if len(z_idx) > 0 and len(n_idx) > 0:
        i, j = z_idx[0], n_idx[0]
        for k, t in enumerate(theta):
            ang_int = np.sum(d2sigma[i, j, k, :])
            if ang_int > 1e-6:
                print(f"  {t:8.1f}  {ang_int:16.6f}")

        # Find peak in angle-energy space
        max_idx = np.unravel_index(np.argmax(d2sigma[i, j, :, :]), d2sigma[i, j, :, :].shape)
        peak_theta = theta[max_idx[0]]
        peak_E = E_vals[max_idx[1]]
        peak_val = d2sigma[i, j, max_idx[0], max_idx[1]]
        print(f"\n  Peak d²σ/dE/dΩ at θ ≈ {peak_theta:.1f}°, E ≈ {peak_E:.0f} MeV")
        print(f"  Peak value: {peak_val:.4e} mb/sr/MeV")
    else:
        print("  ²⁴³U channel not found in angular distribution grid")

    # Also look at all channels summed distribution
    print(f"\nTotal double differential cross-section (summed over all channels):")
    d2sigma_total = np.sum(d2sigma, axis=(0, 1))  # sum over Z,N
    max_idx_tot = np.unravel_index(np.argmax(d2sigma_total), d2sigma_total.shape)
    print(f"  Peak at θ ≈ {theta[max_idx_tot[0]]:.1f}°, E ≈ {E_vals[max_idx_tot[1]]:.0f} MeV")

    # Print the full 2D distribution for all channels
    print(f"\nFull d²σ/dE/dΩ (all channels summed): [θ×E matrix]")
    print(f"  θ↓ \\ E→  ", end="")
    for e_idx in range(0, nE, 10):
        print(f" {E_vals[e_idx]:7.0f}", end="")
    print()
    for k, t in enumerate(theta):
        if k % 3 == 0:  # print every 3rd angle
            print(f"  {t:5.1f}°   ", end="")
            for e_idx in range(0, nE, 10):
                val = d2sigma_total[k, e_idx]
                if val > 0:
                    print(f" {val:7.3f}", end="")
                else:
                    print(f"        ", end="")
            print()

    # Save data for comparison
    out_dir = os.path.join(os.path.dirname(__file__), '..', 'output', 'uu')
    os.makedirs(out_dir, exist_ok=True)

    # Save 243U d2sigma
    if len(z_idx) > 0 and len(n_idx) > 0:
        i, j = z_idx[0], n_idx[0]
        np.savez(os.path.join(out_dir, 'd2sigma_243U.npz'),
                 theta=theta, E=E_vals, d2sigma=d2sigma[i, j, :, :])

    # Save total d2sigma
    np.savez(os.path.join(out_dir, 'd2sigma_total.npz'),
             theta=theta, E=E_vals, d2sigma=d2sigma_total)

    # Also save as CSV for easy inspection
    import pandas as pd
    rows = []
    for k, t in enumerate(theta):
        for e_idx, E in enumerate(E_vals):
            val = d2sigma_total[k, e_idx]
            if val > 1e-6:
                rows.append({'theta_deg': t, 'E_MeV': E, 'd2sigma': val})
    pd.DataFrame(rows).to_csv(os.path.join(out_dir, 'd2sigma_total.csv'), index=False)
    print(f"\nData saved to: {os.path.join(out_dir, '*.npz, *.csv')}")

    print("\n✅ Done!")

if __name__ == '__main__':
    main()
