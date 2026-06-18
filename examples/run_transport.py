#!/usr/bin/env python
"""
Example: Gas cell transport simulation.

Run: python examples/run_transport.py

Demonstrates:
1. Setting up a gas cell geometry (He-filled, 50 mbar)
2. Simulating transport of reaction products
3. Deposition distribution analysis
4. Pressure and length scans for optimization
"""

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

import numpy as np
from mnt_sim.transport.gas_cell import GasCell, CellGeometry
from mnt_sim.transport.monte_carlo import TransportMC
from mnt_sim.transport.stopping import StoppingPower
from mnt_sim.plot import (plot_deposition_distribution,
                          plot_range_curve, plot_pressure_scan)


def main():
    print("=" * 60)
    print("Gas Cell Transport Simulation")
    print("=" * 60)

    # === 1. Gas Cell Setup ===
    print("\n[1] Gas Cell Configuration")
    cell = GasCell(
        gas='He',
        pressure_mbar=50,
        length_mm=200,
        diameter_mm=40,
        temperature_K=293,
        geometry=CellGeometry.CYLINDRICAL
    )
    print(f"    {cell}")
    print(f"    Volume: {cell.volume_cm3:.1f} cm³")
    print(f"    Gas density: {cell.gas_density():.2e} g/cm³")

    # === 2. Stopping Power ===
    print("\n[2] Stopping Power: ²⁰⁸Pb in He")
    stopping = StoppingPower(ion_Z=82, ion_A=208,
                            gas='He', pressure_mbar=50)

    for E_test in [0.1, 1.0, 10.0, 50.0, 100.0]:
        dEdx = stopping.ziegler_stopping(E_test)
        print(f"    E={E_test:6.1f} MeV -> dE/dx = {dEdx:.4e} MeV/(mg/cm²)")

    range_cm, _ = stopping.range_energy(100.0)
    print(f"\n    Range of 100 MeV Pb in He at 50 mbar: {range_cm:.1f} cm")

    # Plot range curve
    out_dir = os.path.join(os.path.dirname(__file__), '..', 'output')
    os.makedirs(out_dir, exist_ok=True)
    plot_range_curve(stopping, 100.0,
                    save_path=os.path.join(out_dir, 'range_curve.png'))
    print("    Saved: output/range_curve.png")

    # === 3. MC Transport ===
    print("\n[3] Monte Carlo Transport: ²⁰⁸Pb @ 100 MeV")
    mc = TransportMC(
        cell,
        ion_Z=82, ion_A=208,
        energy=100.0,
        n_particles=5000,
        sigma_theta=3.0,
        seed=42
    )

    result = mc.run(step_size_mm=0.5)
    print(f"\n{result.summary()}")

    # === 4. Deposition Distribution ===
    plot_deposition_distribution(
        result,
        title='²⁰⁸Pb Deposition in He Gas Cell (50 mbar, 100 MeV)',
        save_path=os.path.join(out_dir, 'deposition_dist.png')
    )
    print("\n[4] Saved: output/deposition_dist.png")

    if len(result.deposition_positions) > 0:
        xy = result.deposition_positions
        print(f"    Deposition range:")
        print(f"    x: [{xy[:, 0].min():.1f}, {xy[:, 0].max():.1f}] mm")
        print(f"    y: [{xy[:, 1].min():.1f}, {xy[:, 1].max():.1f}] mm")
        print(f"    z: [{xy[:, 2].min():.1f}, {xy[:, 2].max():.1f}] mm")

    # === 5. Pressure Scan ===
    print("\n[5] Pressure Scan (20-100 mbar)")
    pressures = [20, 30, 50, 70, 100]
    scan_results = mc.scan_pressure(pressures)

    print(f"    {'P (mbar)':>10} {'Deposited':>10} {'Escaped':>10} {'Efficiency':>12}")
    print(f"    {'-'*42}")
    for P in pressures:
        r = scan_results[P]
        print(f"    {P:10.0f} {r.n_deposited:10d} {r.n_escaped:10d} "
              f"{r.deposition_efficiency*100:10.1f}%")

    plot_pressure_scan(scan_results,
                      save_path=os.path.join(out_dir, 'pressure_scan.png'))
    print("    Saved: output/pressure_scan.png")

    # === 6. Length Scan ===
    print("\n[6] Length Scan (100-400 mm)")
    lengths = [100, 150, 200, 300, 400]
    len_results = mc.scan_length(lengths)
    print(f"    {'L (mm)':>8} {'Deposited':>10} {'Efficiency':>12}")
    print(f"    {'-'*30}")
    for L in lengths:
        r = len_results[L]
        print(f"    {L:8.0f} {r.n_deposited:10d} {r.deposition_efficiency*100:10.1f}%")

    print("\n✅ Done!")


if __name__ == '__main__':
    main()
