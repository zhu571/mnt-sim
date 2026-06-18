"""ImQMD-style Monte Carlo: individual events per partial wave.
Note: _theta_cm() uses a simplified linear model as documented.
Full ImQMD deflection requires numerical inversion of the potential.
"""
import numpy as np
import os
import logging
from dataclasses import dataclass

logger = logging.getLogger(__name__)


@dataclass
class MNTResult:
    Z: np.ndarray
    N: np.ndarray
    sigma_matrix: np.ndarray
    sigma_total: float
    E_cm: float


class ImQMD_MC_Model:
    """Monte Carlo version: each partial wave = individual events."""

    M_N = 931.494
    HBARC = 197.33
    E2 = 1.44

    def __init__(self, projectile="U", target="U", E_lab=7.0,
                 Ap=238, At=238):
        self._Z_map = {
            'H':1,'He':2,'Li':3,'Be':4,'B':5,'C':6,'N':7,'O':8,'F':9,'Ne':10,
            'Na':11,'Mg':12,'Al':13,'Si':14,'P':15,'S':16,'Cl':17,'Ar':18,
            'K':19,'Ca':20,'Sc':21,'Ti':22,'V':23,'Cr':24,'Mn':25,'Fe':26,
            'Co':27,'Ni':28,'Cu':29,'Zn':30,'Ga':31,'Ge':32,'As':33,'Se':34,
            'Br':35,'Kr':36,'Rb':37,'Sr':38,'Y':39,'Zr':40,'Nb':41,'Mo':42,
            'Tc':43,'Ru':44,'Rh':45,'Pd':46,'Ag':47,'Cd':48,'In':49,'Sn':50,
            'Sb':51,'Te':52,'I':53,'Xe':54,'Cs':55,'Ba':56,'La':57,'Ce':58,
            'Pr':59,'Nd':60,'Pm':61,'Sm':62,'Eu':63,'Gd':64,'Tb':65,'Dy':66,
            'Ho':67,'Er':68,'Tm':69,'Yb':70,'Lu':71,'Hf':72,'Ta':73,'W':74,
            'Re':75,'Os':76,'Ir':77,'Pt':78,'Au':79,'Hg':80,'Tl':81,'Pb':82,
            'Bi':83,'Po':84,'At':85,'Rn':86,'Fr':87,'Ra':88,'Ac':89,'Th':90,
            'Pa':91,'U':92,'Np':93,'Pu':94,'Am':95,'Cm':96,'Bk':97,'Cf':98,
        }
        self.Zp = self._Z_map.get(projectile.capitalize(), 92)
        self.Zt = self._Z_map.get(target.capitalize(), 92)
        self.Ap, self.At = Ap, At
        self.E_lab = E_lab
        self.mu_u = Ap * At / (Ap + At)
        self.mu_MeV = self.M_N * self.mu_u
        self.E_cm = E_lab * Ap * At / (Ap + At)
        self.E_lab_total = E_lab * Ap
        self._lambda_sq = self.HBARC**2 / (2 * self.mu_MeV * self.E_cm)

        r0 = 1.18
        self.R_cont = r0 * (Ap**(1/3) + At**(1/3))
        self.V_C = self.E2 * self.Zp * self.Zt / self.R_cont
        self.V_nuc = -25.0
        self.V_B_eff = self.V_C + self.V_nuc
        self.E_above = max(self.E_cm - self.V_B_eff, 0.0)

        v_A = np.sqrt(2 * self.E_lab_total / self.Ap)
        self.v_cm = self.Ap * v_A / (self.Ap + self.At)

        # Load SRIM data (cached, computed once)
        self._srim_E = np.array([0, 2000])
        self._srim_dedx = np.array([0, 50])
        self._srim_range = np.array([0, 2000])
        self._load_srim()

    def _load_srim(self):
        srim_path = os.path.join(os.path.dirname(__file__), "..", "..", "output", "srim_data.npz")
        srim_path = os.path.normpath(srim_path)
        if not os.path.exists(srim_path):
            logger.warning(f"SRIM data not found at {srim_path}")
            return
        try:
            d = np.load(srim_path)
            self._srim_E = d["E_MeV"]
            self._srim_dedx = d["dedx"]
            # Compute range table ONCE and cache
            rng = np.zeros(len(self._srim_E))
            for i in range(1, len(self._srim_E)):
                de = self._srim_E[i] - self._srim_E[i-1]
                s = (self._srim_dedx[i] + self._srim_dedx[i-1]) / 2
                if s > 0: rng[i] = rng[i-1] + de / s
            self._srim_range = rng
            logger.info(f"Loaded SRIM data from {srim_path}")
        except Exception as e:
            logger.error(f"Failed to load SRIM data: {e}")

    def _range(self, E):
        """Cached SRIM range lookup."""
        if E <= 0: return 0.0
        return max(0, np.interp(E, self._srim_E, self._srim_range))

    def _energy_at_range(self, R):
        """Cached SRIM inverse range lookup."""
        if R <= 0: return 0.0
        if R >= self._srim_range[-1]: return 0.0
        return float(np.interp(R, self._srim_range, self._srim_E))

    def _T_l(self, l, hw=4.0):
        E_rot = self.HBARC**2 * l * (l+1) / (2 * self.mu_MeV * self.R_cont**2)
        dE = self.V_B_eff + E_rot - self.E_cm
        if dE < -5: return 1.0
        return 1.0 / (1.0 + np.exp(2*np.pi*dE/hw))

    def _theta_cm(self, l):
        """Deflection: 170° at l=0 → 80° at l_max.
        NOTE: Simplified linear model. Full ImQMD requires numerical inversion."""
        l_max = np.sqrt(2 * self.mu_MeV * self.R_cont**2 *
                       max(self.E_above, 0.1)) / self.HBARC
        if l_max < 1: l_max = 1
        return np.radians(170 - 90 * (l / l_max))

    def _transfer_prob(self, dZ, dN, l):
        """ImQMD transfer probability for channel (dZ, dN) at partial wave l."""
        T_l = self._T_l(l)
        E_star = self.E_above * T_l * 0.5 + 0.5
        a = (self.Ap + self.At) / 10.0
        T_temp = np.sqrt(max(E_star / a, 0.01))
        Q = 6*dZ + 8*dN - 0.3*(dZ**2+dN**2) + 0.1*dZ*dN
        Q_opt = -2 * np.sqrt(self.E_cm / (self.Ap + self.At))
        Qf = np.exp(-(Q - Q_opt)**2 / (2 * (1.5*T_temp)**2))
        sigma_tr = 1.5 + 0.2 * np.sqrt(max(self.E_above, 0)) + 0.3 * T_temp
        return np.exp(-(dZ**2 + dN**2) / (2*sigma_tr**2)) * Qf

    def _hivap_prob(self, A_prim, E_star):
        """Probability that primary fragment A_prim → 243U."""
        if not (235 <= A_prim <= 250):
            return 0.0
        B_n = 6.0
        a = A_prim / 8.0
        T = np.sqrt(max(E_star / a, 0.01))
        n_max = max(0, int(E_star / B_n))
        prob = 0.0
        for n in range(n_max + 1):
            A_f = A_prim - n
            p = np.exp(-n * B_n / T)
            if A_f == 243:
                prob += p
            elif abs(A_f - 243) <= 2:
                prob += p * 0.3
        return min(prob, 1.0)

    def _two_body_kinematics(self, mC, mD, theta_cm_deg, Q=0):
        """Lab frame (E_lab, theta_lab)."""
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

    def _sample_transfer_channel(self, l, Zv, Nv):
        """Sample (dZ, dN) weighted by ImQMD transfer probability."""
        probs = np.array([[self._transfer_prob(dz, dn, l)
                           for dn in Nv] for dz in Zv])
        total = probs.sum()
        if total <= 0:
            return Zv[len(Zv)//2], Nv[len(Nv)//2]  # fallback to center
        probs_flat = probs.ravel() / total
        idx = np.random.choice(len(probs_flat), p=probs_flat)
        return Zv[idx // len(Nv)], Nv[idx % len(Nv)]

    def generate_events(self, n_events=50000, dz_range=(-6,6), dn_range=(-8,8),
                        target_mgcm2=95.0, n_l_max=200, seed=42):
        """
        MC event generation weighted by ImQMD transfer probability.
        Returns array of (theta_lab, E_lab, weight).
        """
        np.random.seed(seed)
        total_sigma = 0.0
        l_vals, sigma_l_vals, theta_cm_vals = [], [], []
        for l in range(n_l_max):
            T_l = self._T_l(l)
            if T_l < 1e-6: continue
            sigma_l = np.pi * self._lambda_sq * (2*l+1) * T_l * 10.0
            total_sigma += sigma_l
            l_vals.append(l)
            sigma_l_vals.append(sigma_l)
            theta_cm_vals.append(np.degrees(self._theta_cm(l)))

        l_vals = np.array(l_vals)
        sigma_l_vals = np.array(sigma_l_vals)
        theta_cm_vals = np.array(theta_cm_vals)
        prob_l = sigma_l_vals / total_sigma
        Zv = np.arange(dz_range[0], dz_range[1]+1)
        Nv = np.arange(dn_range[0], dn_range[1]+1)

        # Precompute channel prob normalization for sampled l values
        results = []
        for _ in range(n_events):
            l = np.random.choice(l_vals, p=prob_l)
            idx = np.where(l_vals == l)[0][0]
            sig_l = sigma_l_vals[idx]
            theta_cm_val = theta_cm_vals[idx]

            # Sample (dZ, dN) from ImQMD transfer probability
            dZ, dN = self._sample_transfer_channel(l, Zv, Nv)
            transfer_p = self._transfer_prob(dZ, dN, l)

            mC = self.Ap + dZ + dN
            mD = self.At - dZ - dN
            if mC <= 0 or mD <= 0: continue

            Q = 6*dZ + 8*dN - 0.3*(dZ**2+dN**2) + 0.1*dZ*dN
            E_star = self.E_above * 0.6 * np.exp(-l / (max(l_vals)*0.6)) + 2.0
            hivap_p = self._hivap_prob(int(mC), E_star)
            if hivap_p < 0.001: continue

            E_lab, th_lab = self._two_body_kinematics(mC, mD, theta_cm_val, Q)
            if E_lab <= 0 or th_lab < 0 or th_lab > 90: continue

            R_init = self._range(E_lab)
            if R_init <= 0: continue
            remaining = np.random.uniform(0, min(R_init, 95.0))
            depth_weight = np.sqrt(remaining / max(R_init, 1.0))
            E_out = self._energy_at_range(remaining)
            if E_out < 20: continue

            # Weight = sig_l × transfer_prob × hivap × depth_weight
            weight = sig_l * transfer_p * hivap_p * depth_weight
            results.append((th_lab, E_out, weight))

        return np.array(results) if results else np.zeros((0, 3))
