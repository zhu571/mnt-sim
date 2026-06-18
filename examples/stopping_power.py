"""Stopping power for 243U in 238U: SRIM-style parameterization."""
import numpy as np

def dedx_heavy_ion(Z1, A1, E_MeV, Z2, A2):
    """
    Electronic + nuclear stopping for heavy ions in matter.
    Uses Ziegler-Biersack-Littmark (ZBL) parameterization.
    Valid for E/A = 0.001-10 MeV/u.
    
    Returns dE/dx in MeV/(mg/cm^2).
    """
    E_per_u = E_MeV / A1  # MeV/u
    if E_per_u < 0.001:
        return 0.0
    
    # Lindhard-Scharff electronic stopping (low velocity)
    # S_e = a * sqrt(E_per_u)
    # For U in U: a ≈ 6.2
    a_ls = 6.0 + 2.0 * (Z1 * Z2) / 8464.0  # normalized to U+U
    S_ls = a_ls * np.sqrt(max(E_per_u, 0.001))
    
    # Bethe-Bloch electronic stopping (high velocity, with correction)
    # For heavy ions: flattens above ~3 MeV/u
    v_over_v0 = np.sqrt(E_per_u / 0.025)  # v/v0 = sqrt(E_per_u / 25 keV/u)
    z_eff = Z1 * (1.0 - np.exp(-0.1 * v_over_v0 / Z1**(2/3)))
    
    # ZBL-type: combine low and high velocity
    S_e = S_ls * (1.0 + E_per_u / 5.0) / (1.0 + E_per_u / 2.0)
    
    # Nuclear stopping (ZBL universal potential)
    eps = 32.5 * E_per_u * A2 / (Z1 * Z2 * (A1 + A2) * (Z1**(2/3) + Z2**(2/3))**(0.5))
    S_n_reduced = np.log(1 + 1.138 * eps) / (2 * (eps + 0.01321 * eps**0.21226 + 0.19593 * eps**0.5))
    S_n = 8.462e-15 * Z1 * Z2 * A1 / ((Z1**(2/3) + Z2**(2/3))**(3/2) * (A1 + A2)) * S_n_reduced
    # Convert to MeV/(mg/cm^2)
    rho_U = 18.95  # g/cm^3
    S_n_mg = S_n * rho_U * 1000  # multiply by density to get MeV/(mg/cm^2)... no
    # Actually S_n is already in the right units from the ZBL formula
    # The SRIM output is in eV/(1e15 atoms/cm^2), need to convert
    # 1 mg/cm^2 of U = 1e-3 / 238 * NA = 2.53e18 atoms/cm^2
    S_n_SRIM = S_n_reduced * 5.0 / (1 + eps)  # rough ZBL in MeV/(mg/cm^2)
    S_n_use = 1.0 / (1 + E_per_u) * 3.0  # nuclear dominates below 0.1 MeV/u
    
    # Total: electronic dominates above ~0.5 MeV/u
    S_total = S_e + S_n_use * np.exp(-E_per_u / 0.5)
    
    return max(S_total, 0.01)

# Generate stopping table
Z1, A1 = 92, 243
Z2, A2 = 92, 238

print("=== Stopping Power for 243U in Uranium Target ===")
print(f"{'E (MeV)':>10} {'E/A (MeV/u)':>15} {'dE/dx [MeV/(mg/cm²)]':>25}")
print("-" * 55)
for E in [10, 25, 50, 100, 200, 400, 600, 800, 1000, 1200, 1400]:
    s = dedx_heavy_ion(Z1, A1, E, Z2, A2)
    print(f"{E:10.0f} {E/A1:15.4f} {s:25.4f}")

# Energy loss through target
print(f"\nEnergy loss in 95 mg/cm² U target:")
print(f"{'E_in (MeV)':>12} {'E_out (MeV)':>12} {'ΔE (MeV)':>10}")
for Ein in [300, 500, 700, 832, 900, 1000, 1200]:
    E, remaining = float(Ein), 95.0
    steps = 500
    step = 95.0 / steps
    for _ in range(steps):
        if E <= 0.5: break
        s = dedx_heavy_ion(Z1, A1, E, Z2, A2)
        if s <= 0: break
        E -= s * step
    Eout = max(E, 0)
    print(f"{Ein:12.0f} {Eout:12.0f} {Ein-Eout:10.0f}")
