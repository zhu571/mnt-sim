#!/usr/bin/env python
"""
Example: MNT cross-section calculation using the DNS model.

Run: python examples/run_cross_section.py

This script demonstrates:
1. Setting up the DNS model for a Xe + Pb system
2. Calculating transfer probabilities for various channels
3. Computing angular and energy distributions
4. Saving results to CSV and plots
"""

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

import numpy as np
from mnt_sim.cross_section.dns import DNSModel
from mnt_sim.cross_section.empirical import EmpiricalModel
from mnt_sim.plot import plot_cross_section_matrix, plot_angle_energy_map


def main():
    # === Configuration ===
    print("=" * 60)
    print("MNT Cross-Section Calculation")
    print("=" * 60)

    # --- DNS Model: Xe + Pb at 8 MeV/u ---
    print("\n[1] DNS Model: ¹³⁶Xe + ²⁰⁸Pb at 8.0 MeV/u")
    dns = DNSModel(projectile="Xe", target="Pb", E_lab=8.0)

    print(f"    Projectile: Z={dns.Zp}, A={dns.Ap}")
    print(f"    Target:      Z={dns.Zt}, A={dns.At}")
    print(f"    E_cm = {dns.E_cm:.1f} MeV")
    print(f"    Coulomb barrier V_C = {dns.V_C:.1f} MeV")
    print(f"    E_cm - V_C = {dns.E_cm_eff:.1f} MeV")

    # --- Calculate cross sections ---
    result = dns.calculate(delta_Z_range=(-4, 4), delta_N_range=(-6, 6))

    print(f"\n[2] Total MNT cross section: {result.sigma_total:.1f} mb")
    df = result.to_dataframe()

    # Top channels
    top = df.nlargest(10, 'sigma_mb')
    print("\n    Top 10 transfer channels:")
    print(f"    {'ΔZ':>3} {'ΔN':>3} {'A':>4}  {'σ (mb)':>10}")
    print(f"    {'-'*25}")
    for _, row in top.iterrows():
        print(f"    {int(row['Z']):+3d} {int(row['N']):+3d} "
              f"{int(row['A']):4d}  {row['sigma_mb']:10.3f}")

    # --- Angular distribution for the top channel ---
    print("\n[3] Angular distribution:")
    theta = np.linspace(0, 30, 15)
    d2sigma = dns.angular_distribution(theta, result)

    top_ch = top.iloc[0]
    print(f"    Top channel: ΔZ={int(top_ch['Z']):+d}, "
          f"ΔN={int(top_ch['N']):+d}")
    for i, t in enumerate(theta):
        ang_int = np.sum(d2sigma[0, 0, i, :])
        if ang_int > 0:
            print(f"    θ = {t:5.1f}°  dσ/dΩ = {ang_int:.3f} mb/sr")

    # --- Plot cross-section matrix ---
    print("\n[4] Generating plots...")
    out_dir = os.path.join(os.path.dirname(__file__), '..', 'output')
    os.makedirs(out_dir, exist_ok=True)

    plot_cross_section_matrix(
        result,
        title='136Xe + 208Pb MNT Cross Sections (DNS model)',
        save_path=os.path.join(out_dir, 'cross_section_matrix.png')
    )
    print("    Saved: output/cross_section_matrix.png")

    # --- Empirical model comparison ---
    print("\n[5] Comparison with empirical systematics:")
    emp = EmpiricalModel(Zp=54, Ap=136, Zt=82, At=208, E_lab=8.0)
    emp_result = emp.estimate()

    print(f"    Empirical total MNT: {emp_result.total_mnt:.1f} mb")
    for ch, sigma in sorted(emp_result.channels.items(),
                            key=lambda x: -x[1])[:5]:
        print(f"    {ch:>8s}: {sigma:.2f} mb")

    print(f"\n    DNS total:  {result.sigma_total:.1f} mb")
    print(f"    Empirical: {emp_result.total_mnt:.1f} mb")

    # --- Save data ---
    csv_path = os.path.join(out_dir, 'cross_sections.csv')
    df.to_csv(csv_path, index=False)
    print(f"\n[6] Data saved to: {csv_path}")

    print("\n✅ Done!")


if __name__ == '__main__':
    main()
