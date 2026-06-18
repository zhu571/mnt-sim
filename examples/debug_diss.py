"""Debug: check dissipation for dominant channel's kinematic points."""
import numpy as np, sys
sys.path.insert(0, '.')
from mnt_sim.cross_section.imqmd_like import ImQMDModel

model = ImQMDModel('U','U',7.0, Ap=238, At=238)
result = model.calculate(dz_range=(-8,8), dn_range=(-10,10), n_l=300)

Zv, Nv, sm = result.Z, result.N, result.sigma_matrix
theta_cm = np.linspace(5, 175, 120)

i_max, j_max = np.unravel_index(np.argmax(sm), sm.shape)
dZ, dN = int(Zv[i_max]), int(Nv[j_max])
sig_ch = sm[i_max, j_max]
A_tr = abs(dZ) + abs(dN)
mC = 238 + dZ + dN
mD = 238 - dZ - dN
v_cm = model.v_cm

cm_peak = 90.0
cm_sigma_val = 5.0 + 1.0 * np.sqrt(A_tr)
cm_w = np.sin(np.radians(theta_cm)) * np.exp(-(theta_cm - cm_peak)**2 / (2 * cm_sigma_val**2))
cm_w += np.flip(cm_w)
cm_w = cm_w / cm_w.sum()

Q = 6*dZ + 8*dN - 0.3*(dZ**2+dN**2) + 0.1*dZ*dN
print(f"Dominant channel: dZ={dZ:+d}, dN={dN:+d}, sigma={sig_ch:.1f} mb")
print(f"mC={mC}, mD={mD}, A_tr={A_tr}, Q={Q:.1f}")
print(f"v_cm={v_cm:.4f}, cm_sigma={cm_sigma_val:.1f}")
print(f"{'tc':>5} {'w':>10} {'E_orig':>8} {'dE_diss':>8} {'E_diss':>8} {'theta_lab':>10}")
print("-"*55)

for k in [0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110]:
    tc = theta_cm[k]
    w = cm_w[k]
    
    E_cm_C = (mD/(mC+mD))*(model.E_cm + Q)
    v_cm_C = np.sqrt(2*E_cm_C/mC)
    v2_orig = v_cm**2 + v_cm_C**2 + 2*v_cm*v_cm_C*np.cos(np.radians(tc))
    E_orig = 0.5*mC*v2_orig
    
    E_diss_val = 20.0 + 1.5*max(tc - 30.0, 0.0)
    Q_diss = -E_diss_val
    E_cm_C_d = (mD/(mC+mD))*(model.E_cm + Q + Q_diss)
    if E_cm_C_d > 0:
        v_cm_C_d = np.sqrt(2*E_cm_C_d/mC)
        v2_d = v_cm**2 + v_cm_C_d**2 + 2*v_cm*v_cm_C_d*np.cos(np.radians(tc))
        E_lab_d = 0.5*mC*max(v2_d, 0)
    else:
        E_lab_d = 0.0
    
    # Lab angle
    if E_cm_C > 0:
        v_cm_C_use = np.sqrt(2*E_cm_C/mC)
        num = v_cm_C_use * np.sin(np.radians(tc))
        den = v_cm + v_cm_C_use * np.cos(np.radians(tc))
        th_lab = np.degrees(np.arctan2(num, den)) if (abs(num)>1e-10 or abs(den)>1e-10) else 0
    else:
        th_lab = 0
    
    print(f"{tc:5.1f} {w:10.3e} {E_orig:8.0f} {E_diss_val:8.1f} {E_lab_d:8.0f} {th_lab:10.2f}")
