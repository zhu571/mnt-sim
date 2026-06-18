#!/usr/bin/env python
"""
Full pipeline: MNT cross section → double-differential → gas cell transport.

This integrates the entire workflow from reaction to deposition.
"""

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

import numpy as np
from mnt_sim.cross_section.dns import DNSModel
from mnt_sim.transport.gas_cell import GasCell
from mnt_sim.transport.monte_carlo import TransportMC
from mnt_sim.plot import (plot_cross_section_matrix,
                          plot_deposition_distribution,
                          plot_pressure_scan)


def main():
    print("=" * 60)
    print("Full MNT Simulation Pipeline")
    print("=" * 60)

    out_dir = os.path.join(os.path.dirname(__file__), '..', 'output')
    os.makedirs(out_dir, exist_ok=True)

    # Step 1: DNS cross section
    print("\n[Step 1] DNS Cross Section")
    dns = DNSModel("Xe", "Pb", E_lab=8.0)
    result = dns.calculate(delta_Z_range=(-3, 3), delta_N_range=(-5, 5))
    plot_cross_section_matrix(result,
        save_path=os.path.join(out_dir, 'pipeline_cs_matrix.png'))

    df = result.to_dataframe()
    top = df.nlargest(5, 'sigma_mb')
    print(f"  Total σ = {result.sigma_total:.1f} mb")
    print(f"  Top channels:")
    for _, row in top.iterrows():
        A_ion = 208 + int(row['N']) + int(row['Z']) - top.iloc[0]['Z'].astype(int) if False else 208 + int(row['N']) + int(row['Z']) - 82

    # Step 2: For each top channel, simulate gas cell transport
    print("\n[Step 2] Gas Cell Transport")
    cell = GasCell(gas='He', pressure_mbar=50, length_mm=200, diameter_mm=40)

    for i, (_, row) in enumerate(top.iterrows()):
        dZ = int(row['Z'])
        dN = int(row['N'])
        Z_ion = 54 + dZ
        A_ion = 136 + dZ + dN
        sigma = row['sigma_mb']

        print(f"\n  Channel ΔZ={dZ:+d} ΔN={dN:+d} (A={A_ion}, Z={Z_ion}, σ={sigma:.2f} mb)")

        if sigma < 0.01:
            print("    Skipping (cross section too small)")
            continue

        # Scale entrance energy by Q-value estimate
        E0 = dns.E_cm + 5 * (abs(dZ) + abs(dN)) * 0.5

        mc = TransportMC(cell, Z_ion, A_ion, energy=E0, n_particles=3000)
        tr = mc.run()

        print(f"    E0 = {E0:.1f} MeV")
        print(f"    Deposited: {tr.n_deposited}/{tr.n_particles} "
              f"({tr.deposition_efficiency*100:.1f}%)")

        if i == 0:
            plot_deposition_distribution(
                tr,
                title=f'Channel ΔZ={dZ:+d} ΔN={dN:+d}: Deposition',
                save_path=os.path.join(out_dir, f'pipeline_deposition.png')
            )

    # Step 3: Pressure optimization for the most important channel
    print("\n[Step 3] Pressure Optimization")
    top_row = top.iloc[0]
    dZ, dN = int(top_row['Z']), int(top_row['N'])
    Z_ion = 54 + dZ
    A_ion = 136 + dZ + dN
    E0 = dns.E_cm + 5 * (abs(dZ) + abs(dN)) * 0.5

    mc_opt = TransportMC(cell, Z_ion, A_ion, energy=E0, n_particles=2000)
    pressures = [20, 30, 40, 50, 60, 80, 100]
    scan = mc_opt.scan_pressure(pressures)

    print(f"  {'P (mbar)':>10} {'Efficiency':>12} {'Mean z (mm)':>14}")
    print(f"  {'-'*36}")
    for P in pressures:
        r = scan[P]
        mean_z = (r.deposition_positions[:, 2].mean()
                  if len(r.deposition_positions) > 0 else 0)
        print(f"  {P:10.0f} {r.deposition_efficiency*100:10.1f}% {mean_z:14.1f}")

    plot_pressure_scan(scan,
        save_path=os.path.join(out_dir, 'pipeline_pressure_scan.png'))

    print(f"\n✅ Pipeline complete. Results saved to {out_dir}/")


if __name__ == '__main__':
    main()
