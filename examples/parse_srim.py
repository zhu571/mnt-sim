#!/usr/bin/env python
"""Parse SRIM output - simpler approach."""
import numpy as np

with open("D:/software/SRIM/SR Module/243U in Uranium - Full Range") as f:
    lines = f.readlines()

E_MeV = []
dedx_total = []

for line in lines:
    line = line.strip()
    if not line:
        continue
    # Check if this is a data line (starts with a number followed by unit)
    parts = line.split()
    if len(parts) < 3:
        continue
    
    # First part should be something like "10.00" or "1.00" or "100.00"
    try:
        energy_str = parts[0].replace(',', '')
        # Check what unit is on this line
        if 'keV' in line or 'Kev' in line or 'kev' in line:
            e_mev = float(energy_str) / 1000.0
        elif 'MeV' in line or 'meV' in line or 'Mev' in line:
            e_mev = float(energy_str)
        elif 'GeV' in line or 'GeV' in line or 'gev' in line:
            e_mev = float(energy_str) * 1000.0
        else:
            continue
        
        # Second and third parts should be stopping values in sci notation
        elec = float(parts[1])
        nuc = float(parts[2])
        
        total_kev_mgcm2 = elec + nuc
        total_mev_mgcm2 = total_kev_mgcm2 / 1000.0
        
        E_MeV.append(e_mev)
        dedx_total.append(total_mev_mgcm2)
    except (ValueError, IndexError):
        continue

E_MeV = np.array(E_MeV)
dedx = np.array(dedx_total)

print(f"Parsed {len(E_MeV)} points, E range: {E_MeV[0]:.4f} - {E_MeV[-1]:.0f} MeV")

print(f"\n{'E (MeV)':>10} {'dE/dx [MeV/(mg/cm²)]':>25}")
for i in range(0, len(E_MeV), max(1, len(E_MeV)//15)):
    print(f"{E_MeV[i]:10.4f} {dedx[i]:25.4f}")

# Range
rng = np.zeros(len(E_MeV))
for i in range(1, len(E_MeV)):
    dE = E_MeV[i] - E_MeV[i-1]
    s = (dedx[i] + dedx[i-1]) / 2
    if s > 0: rng[i] = rng[i-1] + dE / s

for E_t in [200, 400, 600, 832, 1000]:
    s = np.interp(E_t, E_MeV, dedx)
    r = np.interp(E_t, E_MeV, rng)
    print(f"  E={E_t:4d}: dE/dx={s:.2f}, Range={r:.1f} mg/cm²")

def energy_loss(E_in, thick=95.0):
    E = float(E_in)
    steps = 500
    step = thick / steps
    for _ in range(steps):
        if E <= 0.01: break
        s = np.interp(E, E_MeV[::-1], dedx[::-1])
        if s <= 0: break
        E -= s * step
    return max(E, 0)

print(f"\nThrough 95 mg/cm² U target:")
for Ein in [832, 900, 1000, 1200]:
    print(f"  {Ein:4d} -> {energy_loss(Ein):5.0f} MeV")

np.savez("D:/work/agent work/mnt-sim/output/srim_data.npz",
         E_MeV=E_MeV, dedx=dedx, range_mgcm2=rng)
print("\nSaved!")
