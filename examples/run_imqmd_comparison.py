#!/usr/bin/env python
"""
Compare ImQMD-like model (v2) with Wang et al. NIM B 463 (2020) Fig. 2.
Two systems:
1. ¹³⁶Xe + ¹⁹⁸Pt @ 7.98 MeV/A → ²⁰⁰Os (θ~60°, E~150 MeV)
2. ²³⁸U + ²³⁸U @ 7.0 MeV/A → ²⁴³U (θ~35-50°, E~600-1000 MeV)
"""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
import numpy as np
from mnt_sim.cross_section.imqmd_like import ImQMDModel


def analyze(label, proj, targ, Elab, Ap, At, dzr, dnr, nl,
            target_nuclide=None, paper_angle=None, paper_energy=None):
    print(f"\n{'='*80}")
    print(f"{label}")
    print(f"{'='*80}")

    model = ImQMDModel(proj, targ, Elab, Ap=Ap, At=At)

    print(f"  {model.Ap}{proj} + {model.At}{targ} @ {Elab} MeV/u")
    print(f"  E_cm = {model.E_cm:.1f} MeV")
    print(f"  V_C  = {model.V_C:.1f} MeV  (pure Coulomb)")
    print(f"  V_nuc = {model.V_nuc_contact:.1f} MeV  (nuclear attraction)")
    print(f"  V_B_eff = {model.V_B_eff:.1f} MeV  (effective barrier)")
    print(f"  E_above = {model.E_above:.1f} MeV  (above effective barrier)")

    result = model.calculate(dz_range=dzr, dn_range=dnr, n_l=nl)

    print(f"\n  Total reaction σ = {result.sigma_total:.1f} mb")
    print(f"  Grazing angle = {result.theta_graz:.1f}°")

    # Top channels
    Zv, Nv, sm = result.Z, result.N, result.sigma_matrix
    ch = []
    for i, z in enumerate(Zv):
        for j, nv in enumerate(Nv):
            if sm[i, j] > 1e-4:
                ch.append((z, nv, sm[i, j]))
    ch.sort(key=lambda x: -x[2])
    print(f"\n  Top channels (σ > 0.1 mb):")
    for z, n, s in ch[:12]:
        print(f"    ΔZ={z:+3d} ΔN={n:+3d}  σ={s:8.4f} mb")

    # Search for specific nuclide
    if target_nuclide:
        tZ, tN = target_nuclide  # target fragment Z,N
        # For 200Os (Z=76,N=124) from 198Pt (Z=78,N=120):
        # Target loses 2p, gains 4n → dZ_proj = +2, dN_proj = -4
        # For 243U (Z=92,N=151) from 238U (Z=92,N=146):
        # Projectile gains 5n → dZ_proj=0, dN_proj=+5
        total = 0
        for z, n, s in ch:
            # PLF: projectile + (dZ, dN) → check if matches
            Z_PLF = model.Zp + z
            N_PLF = (model.Ap - model.Zp) + n
            # TLF: target - (dZ, dN) from projectile perspective
            Z_TLF = model.Zt - z
            N_TLF = (model.At - model.Zt) - n
            if (Z_PLF == tZ and N_PLF == tN) or (Z_TLF == tZ and N_TLF == tN):
                total += s
                print(f"    ✓ {tZ}{tN-match}: ΔZ={z:+d},ΔN={n:+d} σ={s:.4f} mb")
        print(f"    ²⁰⁰Os total: {total:.4f} mb")

    # Double differential cross section
    theta, E_vals, d2 = model.angular_energy_dist(result)
    d2t = np.sum(d2, axis=(0, 1))
    mi = np.unravel_index(np.argmax(d2t), d2t.shape)
    p_th, p_E = theta[mi[0]], E_vals[mi[1]]

    print(f"\n  d²σ/dE/dΩ peak:")
    print(f"    Our:   θ ≈ {p_th:.1f}°, E ≈ {p_E:.0f} MeV")
    if paper_angle:
        am = "✓" if paper_angle[0] <= p_th <= paper_angle[1] else "✗"
        em = "✓" if paper_energy[0] <= p_E <= paper_energy[1] else "✗"
        print(f"    Paper: θ ≈ {paper_angle[0]}-{paper_angle[1]}°, E ≈ {paper_energy[0]}-{paper_energy[1]} MeV  [{am}{em}]")

    # Distribution matrix
    print(f"\n  d²σ/dE/dΩ (×10³):")
    h = "  θ↓\\E→"
    for ei in range(0, len(E_vals), 10):
        h += f" {E_vals[ei]:6.0f}"
    print(h)
    for k in range(0, len(theta), 3):
        l = f"  {theta[k]:5.1f}°"
        for ei in range(0, len(E_vals), 10):
            v = d2t[k, ei] * 1000
            l += f" {v:6.3f}" if v > 0.001 else f"       "
        print(l)

    # Save
    od = os.path.join(os.path.dirname(__file__), '..', 'output')
    os.makedirs(od, exist_ok=True)
    fn = label.split('→')[1].split('@')[0].strip() if '→' in label else label[:20]
    np.savez(os.path.join(od, f'd2_{fn}.npz'), theta=theta, E=E_vals, d2=d2t)
    import csv
    with open(os.path.join(od, f'd2_{fn}.csv'), 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['theta_deg'] + [f'{e:.0f}' for e in E_vals])
        for kk in range(len(theta)):
            w.writerow([f'{theta[kk]:.1f}'] + [f'{d2t[kk, ei]:.6e}' for ei in range(len(E_vals))])

    return result, (theta, E_vals, d2t)


def main():
    rs = []

    # System 1: 136Xe + 198Pt @ 7.98 MeV/A → 200Os
    print("\n" + "=" * 80)
    print("TEST 1: ¹³⁶Xe + ¹⁹⁸Pt @ 7.98 MeV/A → ²⁰⁰Os")
    print("Paper: θ ≈ 60°, E ≈ 150 MeV")
    print("=" * 80)
    r1, d1 = analyze("Xe+Pt→200Os", "Xe", "Pt", 7.98, 136, 198,
                     (-6, 6), (-8, 8), 200,
                     target_nuclide=(76, 124),
                     paper_angle=(45, 75), paper_energy=(100, 250))
    rs.append(("Xe+Pt→²⁰⁰Os", d1, (45, 75), (100, 250)))

    # System 2: 238U + 238U @ 7.0 MeV/A → 243U
    print("\n" + "=" * 80)
    print("TEST 2: ²³⁸U + ²³⁸U @ 7.0 MeV/A → ²⁴³U")
    print("Paper: θ ≈ 35-50°, E ≈ 600-1000 MeV")
    print("=" * 80)
    r2, d2 = analyze("U+U→243U", "U", "U", 7.0, 238, 238,
                     (-8, 8), (-10, 10), 300,
                     target_nuclide=(92, 151),
                     paper_angle=(30, 60), paper_energy=(500, 1100))
    rs.append(("U+U→²⁴³U", d2, (30, 60), (500, 1100)))

    # Summary
    print("\n" + "=" * 80)
    print("FINAL COMPARISON")
    print("=" * 80)
    for name, dd, pa, pe in rs:
        th, Ev, d2t = dd
        mi = np.unravel_index(np.argmax(d2t), d2t.shape)
        a, e = th[mi[0]], Ev[mi[1]]
        am = "✓" if pa[0] <= a <= pa[1] else "✗"
        em = "✓" if pe[0] <= e <= pe[1] else "✗"
        print(f"  {name:>20}: Our(θ={a:5.1f}°, E={e:5.0f} MeV) vs "
              f"Paper(θ={pa[0]}-{pa[1]}°, E={pe[0]}-{pe[1]} MeV) [{am}{em}]")
        print(f"  {'':>20} Total σ = {r1.sigma_total if name.startswith('Xe') else r2.sigma_total:.1f} mb")

    print("\n✅ Done! Data saved to output/")

if __name__ == '__main__':
    main()
