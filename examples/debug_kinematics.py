"""Debug kinematic mapping: where does each theta_cm go in lab?"""
import numpy as np

Ap, At = 238, 238
E_lab_per_u = 7.0
Elab_total = E_lab_per_u * Ap
E_cm = E_lab_per_u * Ap * At / (Ap + At)

p_A = np.sqrt(2 * Ap * Elab_total)
v_cm = p_A / (Ap + At)
print(f"v_cm = {v_cm:.3f}")
print(f"E_cm = {E_cm:.0f} MeV")

# For 243U as PLF
for label, m_C in [("A=238 (primary)", 238), ("A=243 (243U)", 243), ("A=248 (heavy)", 248)]:
    m_D = 476 - m_C
    E_cm_C = (m_D / (m_C + m_D)) * (E_cm + 0)
    v_cm_C = np.sqrt(2 * E_cm_C / m_C)
    print(f"\n--- {label}: m_C={m_C}, v_cm_C={v_cm_C:.3f} ---")
    print(f"  theta_cm -> theta_lab -> E_lab:")
    for tc in [0, 15, 30, 45, 60, 75, 90, 105, 120, 135, 150, 165, 179]:
        tcr = np.radians(tc)
        v2 = v_cm**2 + v_cm_C**2 + 2*v_cm*v_cm_C*np.cos(tcr)
        E_lab = 0.5 * m_C * v2
        num = v_cm_C * np.sin(tcr)
        den = v_cm + v_cm_C * np.cos(tcr)
        if abs(num) < 1e-10 and abs(den) < 1e-10:
            th_lab = 0.0
        else:
            th_lab = np.degrees(np.arctan2(num, den))
        print(f"    θ_cm={tc:3d}° → θ_lab={th_lab:6.1f}°, E_lab={E_lab:7.0f} MeV")
