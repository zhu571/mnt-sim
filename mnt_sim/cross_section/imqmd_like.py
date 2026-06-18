"""ImQMD-like model with proper two-body kinematics + SRIM target energy loss."""
import numpy as np
import os
import logging
from dataclasses import dataclass

from mnt_sim.data import element_to_z, typical_mass_number

logger = logging.getLogger(__name__)

__all__ = ["ImQMDModel", "MNTResult"]


@dataclass
class MNTResult:
    Z: np.ndarray
    N: np.ndarray
    sigma_matrix: np.ndarray
    sigma_total: float
    E_cm: float


class ImQMDModel:
    """ImQMD-inspired model with kinematics + SRIM energy loss."""

    M_N = 931.494
    HBARC = 197.33
    E2 = 1.44

    def __init__(self, projectile="U", target="U", E_lab=7.0,
                 Ap=238, At=238):
        self.projectile = projectile
        self.target = target
        self.E_lab = E_lab
        self.Zp = self._Z(projectile)
        self.Zt = self._Z(target)
        self.Ap = Ap or self._A(projectile)
        self.At = At or self._A(target)
        self.mu_u = (self.Ap * self.At) / (self.Ap + self.At)
        self.mu_MeV = self.M_N * self.mu_u
        self.E_cm = self.E_lab * self.Ap * self.At / (self.Ap + self.At)
        self.E_lab_total = self.E_lab * self.Ap
        self._lambda_sq = (self.HBARC**2) / (2 * self.mu_MeV * self.E_cm)

        r0 = 1.18
        self.R_cont = r0 * (self.Ap**(1/3) + self.At**(1/3))
        self.V_C = self.E2 * self.Zp * self.Zt / self.R_cont
        self.V_nuc = -25.0
        self.V_B_eff = self.V_C + self.V_nuc
        self.E_above = max(self.E_cm - self.V_B_eff, 0.0)

        # CM velocity (correct formula)
        v_A = np.sqrt(2 * self.E_lab_total / self.Ap)
        self.v_cm = self.Ap * v_A / (self.Ap + self.At)

        # Load SRIM data
        self._has_srim = False
        self._init_srim()

    def _init_srim(self):
        srim_path = os.path.join(os.path.dirname(__file__), "..", "..", "output", "srim_data.npz")
        srim_path = os.path.normpath(srim_path)
        if not os.path.exists(srim_path):
            logger.warning(f"SRIM data not found at {srim_path}, skipping target energy loss")
            self._has_srim = False
            return
        try:
            d = np.load(srim_path)
            self._srim_E = d["E_MeV"]
            self._srim_dedx = d["dedx"]
            rng = np.zeros(len(self._srim_E))
            for i in range(1, len(self._srim_E)):
                de = self._srim_E[i] - self._srim_E[i-1]
                s = (self._srim_dedx[i] + self._srim_dedx[i-1]) / 2
                if s > 0: rng[i] = rng[i-1] + de / s
            self._srim_range = rng
            self._has_srim = True
            logger.info(f"Loaded SRIM data from {srim_path}")
        except Exception as e:
            logger.error(f"Failed to load SRIM data: {e}")
            self._has_srim = False

    def _E_out(self, E_in, remaining):
        if not self._has_srim or E_in <= 0 or remaining <= 0:
            return E_in
        E = float(E_in)
        steps = max(100, int(remaining * 3))
        step = remaining / steps
        for _ in range(steps):
            if E <= 0.01: return 0.0
            s = np.interp(E, self._srim_E[::-1], self._srim_dedx[::-1])
            if s <= 0: break
            E -= s * step
        return max(E, 0)

    def _Z(self, s):
        """Convert element symbol to atomic number."""
        return element_to_z(s)

    def _A(self, s):
        """Return typical mass number for an element."""
        return typical_mass_number(s)

    def _transmission(self, l, hw=4.0):
        E_rot = self.HBARC**2 * l * (l+1) / (2 * self.mu_MeV * self.R_cont**2)
        dE = self.V_B_eff + E_rot - self.E_cm
        if dE < -5: return 1.0
        return 1.0 / (1.0 + np.exp(2*np.pi*dE/hw))

    def _sigma_l(self, l):
        T = self._transmission(l)
        return np.pi * self._lambda_sq * (2*l+1) * T * 10.0

    def _transfer_prob(self, dZ, dN, l):
        T = self._transmission(l)
        E_star = self.E_above * T * 0.5 + 0.5
        a = (self.Ap + self.At) / 10.0
        Ttemp = np.sqrt(max(E_star / a, 0.01))
        Q = 6*dZ + 8*dN - 0.3*(dZ**2+dN**2) + 0.1*dZ*dN
        Q_opt = -2 * np.sqrt(self.E_cm / (self.Ap + self.At))
        Qf = np.exp(-(Q - Q_opt)**2 / (2 * (1.5*Ttemp)**2))
        sigma_tr = 1.5 + 0.2 * np.sqrt(max(self.E_above, 0)) + 0.3 * Ttemp
        return np.exp(-(dZ**2 + dN**2) / (2*sigma_tr**2)) * Qf

    def calculate(self, dz_range=(-8, 8), dn_range=(-10, 10), n_l=300):
        Zv = np.arange(dz_range[0], dz_range[1]+1)
        Nv = np.arange(dn_range[0], dn_range[1]+1)
        nZ, nN = len(Zv), len(Nv)
        if nZ == 0 or nN == 0:
            raise ValueError("transfer ranges must not be empty")
        sg = np.zeros((nZ, nN))
        l_max_eff = int(np.sqrt(2 * self.mu_MeV * self.R_cont**2 *
                       max(self.E_above, 0.1)) / self.HBARC)
        l_max = min(l_max_eff + 30, n_l)
        for l in range(l_max):
            Tl = self._transmission(l)
            if Tl < 1e-8: continue
            sig_l = self._sigma_l(l)
            P = np.zeros((nZ, nN))
            Pt = 0.0
            for i, dz in enumerate(Zv):
                for j, dn in enumerate(Nv):
                    p = self._transfer_prob(dz, dn, l)
                    P[i, j] = p; Pt += p
            if Pt > 0: P /= Pt
            sg += sig_l * P
        return MNTResult(Z=Zv, N=Nv, sigma_matrix=sg,
                         sigma_total=np.sum(sg), E_cm=self.E_cm)

    def _tk(self, dZ, dN, theta_cm_deg, Q=0):
        """Two-body kinematics: returns (E_lab, theta_lab)."""
        mC = self.Ap + dZ + dN
        mD = self.At - dZ - dN
        if mC <= 0 or mD <= 0: return 0.0, 0.0
        tc = np.radians(theta_cm_deg)
        E_cm_C = (mD/(mC+mD)) * (self.E_cm + Q)
        if E_cm_C <= 0: return 0.0, 0.0
        vc = np.sqrt(2 * E_cm_C / mC)
        v2 = self.v_cm**2 + vc**2 + 2*self.v_cm*vc*np.cos(tc)
        E_lab = 0.5 * mC * v2
        num = vc * np.sin(tc)
        den = self.v_cm + vc * np.cos(tc)
        if abs(num) < 1e-10 and abs(den) < 1e-10:
            th_lab = 0.0
        else:
            th_lab = np.degrees(np.arctan2(num, den))
        return max(E_lab, 0), max(0, min(90, th_lab))

    def angular_energy_dist(self, result, n_theta=90, nE=200,
                           target_mgcm2=95.0, seed=42):
        """
        Build d²σ/dE/dΩ with wider angular tail + SRIM target energy loss.
        """
        if seed is not None:
            np.random.seed(seed)
        theta_lab_grid = np.linspace(0, 90, n_theta)
        E_grid = np.linspace(0, 1800, nE)
        dE = E_grid[1] - E_grid[0]; dTh = theta_lab_grid[1] - theta_lab_grid[0]
        d2 = np.zeros((n_theta, nE))
        nZ, nN = len(result.Z), len(result.N); sm = result.sigma_matrix
        theta_cm = np.linspace(5, 175, 120)

        for i in range(nZ):
            for j in range(nN):
                sig_ch = sm[i, j]
                if sig_ch < 1e-10: continue
                dZ, dN = result.Z[i], result.N[j]; A_tr = abs(dZ)+abs(dN)

                # CM ang dist: wider tail (20%), cutoff at 10°
                core = np.exp(-(theta_cm-90)**2 / (2*(5+np.sqrt(A_tr))**2))
                raw = 1/np.sin(np.radians(theta_cm)+0.1)
                cutoff = 1/(1+np.exp((10-theta_cm)/4))
                cutoff *= 1/(1+np.exp((theta_cm-170)/4))
                cm_w = 0.80*core + 0.20*(raw*cutoff)/(raw*cutoff).sum()*core.sum()
                cm_w /= cm_w.sum()

                energies = np.zeros(len(theta_cm))
                lab_angles = np.zeros(len(theta_cm))
                for k, tc in enumerate(theta_cm):
                    E_lab, th_lab = self._tk(dZ, dN, tc)
                    energies[k] = E_lab; lab_angles[k] = th_lab
                valid = (energies>0)&(lab_angles>=0)&(lab_angles<=90)
                if not np.any(valid): continue

                sig_th = 0.5 + 0.3*np.sqrt(A_tr)
                sig_E = 10 + 5*np.sqrt(A_tr)
                for k in range(len(theta_cm)):
                    if not valid[k]: continue
                    # Initial lab energy and angle
                    E0_init = energies[k]
                    th0 = lab_angles[k]

                    # Target energy loss: random production depth
                    # Surviving fragments have remaining target in [0, R(E)]
                    # where R is SRIM range. Sample uniformly and weight.
                    if self._has_srim:
                        # Range at this energy
                        R_max = np.interp(E0_init, self._srim_E, self._srim_range
                                         ) if E0_init < self._srim_E[-1] else 100
                        R_max = min(R_max, target_mgcm2)
                        if R_max > 1:
                            remaining = np.random.uniform(0, R_max)
                            E0 = self._E_out(E0_init, remaining)
                        else:
                            E0 = E0_init
                    else:
                        E0 = E0_init

                    if E0 < 20: continue
                    g_th = np.exp(-(theta_lab_grid-th0)**2/(2*sig_th**2))
                    g_E = np.exp(-(E_grid-E0)**2/(2*sig_E**2))
                    d2 += sig_ch*cm_w[k]*np.outer(g_th, g_E)

        d2 = np.clip(d2, 0, None)
        total = np.sum(d2)*dE*dTh
        if total > 0: d2 *= result.sigma_total/total
        return theta_lab_grid, E_grid, d2[np.newaxis, np.newaxis, :, :]
