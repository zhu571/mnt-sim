"""Debug why 243U yield is zero."""
import sys; sys.path.insert(0, '.')
import numpy as np
from mnt_sim.cross_section.imqmd_hivap import ImQMD_HIVAP_Model

model = ImQMD_HIVAP_Model('U', 'U', 7.0, Ap=238, At=238)
l = 100
T_l = model._transmission(l)
theta_cm = np.degrees(model._deflection_angle_cm(l))
sigma_l = np.pi * model._lambda_sq * (2*l+1) * T_l * 10.0

E_star = model.E_above * T_l * 0.8 * np.exp(-l/(173*0.6)) + 2.0
print(f"l={l}, E_star={E_star:.1f} MeV, n_evap_max={int(E_star/6.0)}")

# Channels that could produce 243U
test_channels = [
    (0, 5, "direct"),
    (0, 6, "direct+1n"),
    (0, 7, "direct+2n"),
    (1, 4, "1p + 4n transfer"),
    (1, 5, "1p + 5n transfer"),
    (-1, 6, "-1p + 6n transfer"),
    (2, 3, "2p + 3n"),
    (2, 4, "2p + 4n"),
]

for dZ, dN, label in test_channels:
    m_C_prim = 238 + dZ + dN
    Z_prim = 92 + dZ
    A_prim = int(m_C_prim)
    dZf = abs(Z_prim - 92)
    dAf = abs(A_prim - 243)
    
    if not (dZf <= 2 and dAf <= 8):
        print(f"  {label:>20}: dZ={dZ:+d}, dN={dN:+d}, A={A_prim} -> FILTERED (dZf={dZf}, dAf={dAf})")
        continue
    
    n_max = max(0, int(E_star / 6.0))
    prob = 0.0
    hits = []
    for n_evap in range(n_max + 1):
        for p_evap in range(min(1, dZf) + 1):
            A_f2 = A_prim - n_evap - p_evap
            Z_f2 = Z_prim - p_evap
            if Z_f2 == 92 and abs(A_f2 - 243) <= 2:
                p_n = np.exp(-n_evap * 6.0 / max(E_star/max(n_max,1), 0.5))
                p_p = 0.3 if p_evap > 0 else 1.0
                contrib = p_n * p_p * 0.1
                prob += contrib
                hits.append((n_evap, p_evap, A_f2, contrib))
    
    if prob > 0:
        print(f"  {label:>20}: dZ={dZ:+d}, dN={dN:+d}, A={A_prim}, PROB={prob:.6f}, sigma={sigma_l*prob:.4f} mb")
        for n, p, a, c in hits:
            print(f"    -> n_evap={n}, p_evap={p}, A_final={a}, contrib={c:.6f}")
    else:
        print(f"  {label:>20}: dZ={dZ:+d}, dN={dN:+d}, A={A_prim}, NO MATCH (n_max={n_max}, E_star={E_star:.1f})")
