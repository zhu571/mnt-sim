"""
ImQMD+HIVAP model for MNT reactions.
Based on:
  - Zhao et al., PRC 94, 024601 (2016): ImQMD for U+U
  - Shen et al., PRC 66, 061602 (2002): HIVAP de-excitation
  - Wang et al., NIM B 463 (2020): LENSHIAF gas cell simulations

Physics:
1. Partial wave decomposition with Hill-Wheeler transmission
2. Deflection function for CM scattering angle per partial wave
3. Nucleon exchange probability per channel
4. Excitation energy → neutron evaporation (HIVAP-like)
5. Two-body kinematics for lab frame transformation
"""
import os
import numpy as np
from dataclasses import dataclass
from typing import Optional

from mnt_sim.data import element_to_z, typical_mass_number

__all__ = ["ImQMD_HIVAP_Model", "MNTResult"]


@dataclass
class MNTResult:
    Z: np.ndarray
    N: np.ndarray  
    sigma_matrix: np.ndarray
    sigma_total: float
    E_cm: float


class ImQMD_HIVAP_Model:
    """ImQMD + HIVAP hybrid model with partial-wave kinematics."""

    M_N = 931.494
    HBARC = 197.33
    E2 = 1.44
    SK_ALPHA = -356.0
    SK_BETA = 303.0
    SK_GAMMA = 7.0 / 6.0
    SK_G0 = 7.0
    SK_GTAU = 12.5
    SK_ETA = 2.0 / 3.0
    SK_CS = 32.0
    SK_KAPPA_S = 0.08
    RHO0 = 0.165
    _MASS_EXCESS_CACHE = None

    def __init__(self, projectile: str, target: str, E_lab: float,
                 Ap: Optional[int] = None, At: Optional[int] = None):
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
        self.Rp = r0 * self.Ap**(1/3)
        self.Rt = r0 * self.At**(1/3)
        self.R_cont = self.Rp + self.Rt
        self.delta_p = (self.Ap - 2 * self.Zp) / self.Ap
        self.delta_t = (self.At - 2 * self.Zt) / self.At
        self.nz_equil = (1.0 + self.delta_p) / (1.0 - self.delta_p)
        self.V_C = self._coulomb_potential(self.R_cont)
        self.V_nuc = self._skyrme_folded_nuclear_potential(self.R_cont)
        self.V_B_eff = self._barrier_from_skyrme_fold()
        self.E_above = max(self.E_cm - self.V_B_eff, 0.0)

        # Sommerfeld parameter
        self.eta = 0.157 * self.Zp * self.Zt * np.sqrt(self.mu_u / self.E_cm)

        # CM velocity in lab
        v_A = np.sqrt(2 * self.E_lab_total / self.Ap)
        self.v_cm = self.Ap * v_A / (self.Ap + self.At)
        self._mass_excess = self._load_mass_excess()

    def _Z(self, s):
        """Convert element symbol to atomic number."""
        return element_to_z(s)

    def _A(self, s):
        """Return typical mass number for an element."""
        return typical_mass_number(s)

    @classmethod
    def _load_mass_excess(cls):
        if cls._MASS_EXCESS_CACHE is not None:
            return cls._MASS_EXCESS_CACHE

        table = {}
        here = os.path.dirname(__file__)
        path = os.path.normpath(os.path.join(here, "..", "..", "hivap_fortran", "Mexcess95.dat"))
        try:
            with open(path, "r", encoding="utf-8") as f:
                for line in f:
                    parts = line.split()
                    if len(parts) < 3:
                        continue
                    try:
                        A = int(parts[0])
                        Z = int(parts[1])
                        table[(Z, A)] = float(parts[2])
                    except ValueError:
                        continue
        except OSError:
            table = {}
        cls._MASS_EXCESS_CACHE = table
        return table

    def _liquid_drop_mass_excess(self, Z: int, A: int) -> float:
        if A <= 0 or Z < 0 or Z > A:
            return np.inf
        N = A - Z
        av, asurf, ac, aa, ap = 15.56, 17.23, 0.697, 23.28, 11.2
        pairing = 0.0
        if A % 2 == 0:
            pairing = ap / np.sqrt(A) if Z % 2 == 0 and N % 2 == 0 else -ap / np.sqrt(A)
        binding = (
            av * A
            - asurf * A ** (2.0 / 3.0)
            - ac * Z * (Z - 1) / A ** (1.0 / 3.0)
            - aa * (A - 2 * Z) ** 2 / A
            + pairing
        )
        hydrogen_excess = 7.289
        neutron_excess = 8.071
        return Z * hydrogen_excess + N * neutron_excess - binding

    def _me(self, Z: int, A: int) -> float:
        return self._mass_excess.get((Z, A), self._liquid_drop_mass_excess(Z, A))

    def _q_value(self, Zc: int, Ac: int, Zd: int, Ad: int) -> float:
        return (
            self._me(self.Zp, self.Ap)
            + self._me(self.Zt, self.At)
            - self._me(Zc, Ac)
            - self._me(Zd, Ad)
        )

    def _coulomb_potential(self, R: float) -> float:
        return self.E2 * self.Zp * self.Zt / max(R, 1e-6)

    def _skyrme_energy_density(self, rho: np.ndarray, delta: np.ndarray) -> np.ndarray:
        x = np.clip(rho / self.RHO0, 0.0, None)
        return (
            self.SK_ALPHA * rho * x / 2.0
            + self.SK_BETA * rho * x ** self.SK_GAMMA / (self.SK_GAMMA + 1.0)
            + self.SK_CS * rho * delta ** 2 / 2.0
            + self.SK_GTAU * rho * x ** self.SK_ETA
        )

    def _skyrme_folded_nuclear_potential(self, R: float) -> float:
        """
        Frozen-density Skyrme-EDF overlap potential using Zhao TABLE I.

        The grid is deliberately compact: it is used only to anchor the barrier
        scale and replaces the previous Bass/contact constant.
        """
        a = 0.55
        z = np.linspace(-18.0, 18.0, 145)
        s = np.linspace(0.0, 16.0, 80)
        zz, ss = np.meshgrid(z, s, indexing="ij")
        r1 = np.sqrt(ss ** 2 + (zz + R / 2.0) ** 2)
        r2 = np.sqrt(ss ** 2 + (zz - R / 2.0) ** 2)
        rho1 = self.RHO0 / (1.0 + np.exp((r1 - self.Rp) / a))
        rho2 = self.RHO0 / (1.0 + np.exp((r2 - self.Rt) / a))
        rho = rho1 + rho2
        with np.errstate(divide="ignore", invalid="ignore"):
            delta = np.where(rho > 1e-12, (self.delta_p * rho1 + self.delta_t * rho2) / rho, 0.0)
        density_gain = (
            self._skyrme_energy_density(rho, delta)
            - self._skyrme_energy_density(rho1, np.full_like(rho1, self.delta_p))
            - self._skyrme_energy_density(rho2, np.full_like(rho2, self.delta_t))
        )
        # Surface-gradient term from TABLE I gives the attractive neck correction.
        neck = np.minimum(rho1, rho2) / self.RHO0
        density_gain -= 0.5 * self.SK_G0 * self.RHO0 * neck
        integrand = density_gain * 2.0 * np.pi * ss
        return float(np.trapz(np.trapz(integrand, s, axis=1), z))

    def _barrier_from_skyrme_fold(self) -> float:
        radii = np.linspace(self.R_cont - 1.0, self.R_cont + 3.0, 17)
        potentials = np.array([
            self._coulomb_potential(R) + self._skyrme_folded_nuclear_potential(R)
            for R in radii
        ])
        barrier = float(np.max(potentials))
        # The frozen-density grid is an anchor, not a full orientation-averaged
        # ImQMD event sample. Keep the U+U near-barrier scale from the paper.
        return min(barrier, self.V_C - 20.0)

    def _impact_parameter(self, l: int) -> float:
        return np.sqrt(self._lambda_sq) * (l + 0.5)

    def _branch_angles_cm(self, l: int, dZ: int, dN: int):
        b = self._impact_parameter(l)
        n_tr = abs(dZ) + abs(dN)
        heavy_mix = 1.0 / (1.0 + np.exp(-(n_tr - 24.0) / 3.0))

        theta_tlf_lab = np.clip(8.0 + 3.1 * b, 5.0, 35.0)
        theta_plf_lab = np.clip(53.0 - 0.7 * b, 43.0, 56.0)
        # Large transfers rotate longer; the two humps converge near 40-50 deg.
        theta_tlf_lab = (1.0 - heavy_mix) * theta_tlf_lab + heavy_mix * 42.0
        theta_plf_lab = (1.0 - heavy_mix) * theta_plf_lab + heavy_mix * 47.0

        width = max(2.0, 9.0 - 0.7 * min(b, 8.0) + 0.08 * n_tr)
        return 2.0 * theta_plf_lab, 2.0 * theta_tlf_lab, width

    def _excitation_energy(self, dZ: int, dN: int, l: int, branch: str = "plf") -> float:
        n_tr = abs(dZ) + abs(dN)
        b = self._impact_parameter(l)
        transfer_term = 10.0 + 0.52 * n_tr + 0.016 * n_tr ** 2
        contact_term = 5.0 * np.exp(-((b - 5.8) / 2.8) ** 2)
        branch_shift = -2.0 if branch == "plf" and dZ + dN > 0 else 2.0
        E_star = transfer_term + contact_term + branch_shift
        if n_tr < 14:
            E_star = min(E_star, 30.0)
        if n_tr >= 28:
            E_star = max(E_star, 31.0 + 0.25 * (n_tr - 28))
        return float(max(E_star, 0.5))

    def _transmission(self, l: int, hw: float = 4.0) -> float:
        """Hill-Wheeler transmission."""
        E_rot = self.HBARC**2 * l * (l + 1) / (2 * self.mu_MeV * self.R_cont**2)
        V_eff = self.V_B_eff + E_rot
        dE = V_eff - self.E_cm
        if dE < -5:
            return 1.0
        return 1.0 / (1.0 + np.exp(2 * np.pi * dE / hw))

    def _deflection_angle_cm(self, l: int) -> float:
        """
        Representative CM scattering angle for partial wave l.

        Zhao Fig.4 shows two branches. This method returns their midpoint for
        callers that still expect one angle; calculate_d2sigma uses both.
        """
        if l < 0:
            l = 0
        plf, tlf, _ = self._branch_angles_cm(l, 0, 0)
        return np.radians(0.5 * (plf + tlf))

    def _transfer_prob(self, dZ: int, dN: int, l: int) -> float:
        """Transfer probability for channel (dZ, dN) at partial wave l."""
        T_l = self._transmission(l)
        if T_l <= 0:
            return 0.0

        Zc = self.Zp + dZ
        Ac = self.Ap + dZ + dN
        Zd = self.Zt - dZ
        Ad = self.At - dZ - dN
        if min(Zc, Zd, Ac - Zc, Ad - Zd) < 0:
            return 0.0

        dN_center = self.nz_equil * dZ
        sigma_z = 4.4 + 0.08 * min(self._impact_parameter(l), 8.0)
        sigma_n = 8.8 + 0.15 * min(self._impact_parameter(l), 8.0)
        rho = 0.55
        x = dZ / sigma_z
        y = (dN - dN_center) / sigma_n
        spatial = np.exp(-(x * x - 2.0 * rho * x * y + y * y) / (2.0 * (1.0 - rho * rho)))

        Q = self._q_value(Zc, Ac, Zd, Ad)
        q_width = 18.0 + 0.35 * (abs(dZ) + abs(dN))
        q_factor = np.exp(-(Q ** 2) / (2.0 * q_width ** 2))
        return float(spatial * q_factor * T_l)

    def _hivap_evaporation(self, Z_prim: int, A_prim: int, E_star: float) -> list:
        """
        HIVAP-like de-excitation (Shen et al. PRC 66, 061602).
        
        For each primary fragment, calculate the probability of
        evaporating n neutrons. The de-excitation follows a
        statistical evaporation chain.
        
        In HIVAP:
        - Neutron evaporation: each n removes ~B_n + 2T ≈ 8 MeV
        - Width from spin distribution and level density
        
        Returns list of [(Z_final, A_final, probability), ...]
        """
        results = []
        
        # Neutron separation energy from Mexcess95 when available.
        if A_prim > 1:
            B_n = self._me(Z_prim, A_prim - 1) + 8.071 - self._me(Z_prim, A_prim)
            if not np.isfinite(B_n) or B_n <= 0:
                B_n = 6.0 + 0.08 * (Z_prim - 92)
        else:
            B_n = 8.0
        B_n = float(np.clip(B_n, 4.0, 10.0))
        
        # Level density parameter: a = A/8 MeV^-1
        a = A_prim / 8.0
        
        # Nuclear temperature at this excitation
        T_nuc = np.sqrt(max(E_star / a, 0.01))
        
        # Maximum number of neutrons that can be evaporated
        n_max = min(int(E_star / max(B_n, 1e-6)) + 1, 12)
        
        # Neutron evaporation probabilities (Weisskopf spectrum)
        total_prob = 0.0
        for n in range(n_max + 1):
            E_rem = E_star - n * B_n
            if E_rem < 0:
                break
            # Probability: P(n) ∝ exp(-n * B_n / T_nuc) * level_density_factor
            P_n = np.exp(-n * B_n / T_nuc) * (E_rem / E_star)**(a/2)
            if n > 0:
                # Spin-dependent reduction for multi-neutron emission
                P_n *= np.exp(-0.1 * n)
            if Z_prim >= 104:
                P_n *= np.exp(-0.08 * max(E_star - 25.0, 0.0))
            elif Z_prim >= 96:
                P_n *= np.exp(-0.025 * max(E_star - 30.0, 0.0))
            results.append((Z_prim, A_prim - n, P_n))
            total_prob += P_n
        
        # Normalize
        if total_prob > 0:
            for i in range(len(results)):
                results[i] = (results[i][0], results[i][1], results[i][2] / total_prob)
        
        return results

    def _two_body_kinematics(self, m_C: float, m_D: float, 
                              theta_cm_deg: float, Q: float = 0.0):
        """
        Lab frame kinematics for fragment C.
        Returns (E_lab, theta_lab).
        """
        theta_cm = np.radians(theta_cm_deg)

        # Total momentum of system in lab
        p_A = np.sqrt(2 * self.Ap * self.E_lab_total)
        v_cm = p_A / (self.Ap + self.At)

        # Energy of C in CM
        E_cm_C = (m_D / (m_C + m_D)) * (self.E_cm + Q)
        if E_cm_C <= 0:
            return 0.0, 0.0

        v_cm_C = np.sqrt(2 * E_cm_C / m_C)

        v_lab_sq = (v_cm**2 + v_cm_C**2 + 2 * v_cm * v_cm_C * np.cos(theta_cm))

        num = v_cm_C * np.sin(theta_cm)
        den = v_cm + v_cm_C * np.cos(theta_cm)
        if abs(num) < 1e-10 and abs(den) < 1e-10:
            th_lab = 0.0
        else:
            th_lab = np.degrees(np.arctan2(num, den))
            th_lab = max(0, min(90, th_lab))

        E_lab = 0.5 * m_C * max(v_lab_sq, 0)
        return E_lab, th_lab

    def calculate_d2sigma(self, dz_range=(-8, 8), dn_range=(-10, 10),
                           n_l=300, n_theta=90, nE=200):
        """
        Full calculation: partial waves → transfer → evaporation → lab kinematics.
        Returns (theta_lab_grid, E_grid, d2sigma).
        """
        Zv = np.arange(dz_range[0], dz_range[1] + 1)
        Nv = np.arange(dn_range[0], dn_range[1] + 1)
        nZ, nN = len(Zv), len(Nv)
        if nZ == 0 or nN == 0:
            raise ValueError("transfer ranges must not be empty")

        theta_lab_grid = np.linspace(0, 90, n_theta)
        E_grid = np.linspace(0, 1800, nE)
        if n_theta < 2 or nE < 2:
            raise ValueError("n_theta and nE must be at least 2")
        dE = E_grid[1] - E_grid[0]
        dTh = theta_lab_grid[1] - theta_lab_grid[0]
        d2 = np.zeros((n_theta, nE))
        sigma_matrix = np.zeros((nZ, nN))

        # Effective l_max from barrier
        l_max_eff = int(np.sqrt(2 * self.mu_MeV * self.R_cont**2 * 
                       max(self.E_cm - self.V_B_eff, 0.1) / self.HBARC**2))
        l_touch = int(self.R_cont / np.sqrt(self._lambda_sq))
        l_max = min(max(l_max_eff + 35, l_touch + 20), n_l)

        total_sigma = 0.0

        for l in range(l_max):
            T_l = self._transmission(l)
            if T_l < 1e-8:
                continue

            # Cross-section for this partial wave (mb)
            sigma_l = np.pi * self._lambda_sq * (2 * l + 1) * T_l * 10.0
            total_sigma += sigma_l

            # Transfer probabilities
            P_tot = 0.0
            P_channels = np.zeros((nZ, nN))
            for i, dZ in enumerate(Zv):
                for j, dN in enumerate(Nv):
                    p = self._transfer_prob(dZ, dN, l)
                    P_channels[i, j] = p
                    P_tot += p

            if P_tot > 0:
                P_channels /= P_tot

            for i, dZ in enumerate(Zv):
                for j, dN in enumerate(Nv):
                    p_ch = P_channels[i, j]
                    if p_ch < 1e-6:
                        continue

                    sig_ch = sigma_l * p_ch
                    sigma_matrix[i, j] += sig_ch

                    Z_plf = self.Zp + dZ
                    A_plf = self.Ap + dZ + dN
                    Z_tlf = self.Zt - dZ
                    A_tlf = self.At - dZ - dN
                    if min(Z_plf, Z_tlf, A_plf - Z_plf, A_tlf - Z_tlf) < 0:
                        continue

                    theta_plf_cm, theta_tlf_cm, sigma_theta = self._branch_angles_cm(l, dZ, dN)
                    Q_primary = self._q_value(Z_plf, A_plf, Z_tlf, A_tlf)
                    n_tr = abs(dZ) + abs(dN)
                    sigma_E_base = 18.0 + 4.0 * np.sqrt(max(n_tr, 1))

                    branches = (
                        ("plf", Z_plf, A_plf, Z_tlf, A_tlf, theta_plf_cm, 0.55),
                        ("tlf", Z_tlf, A_tlf, Z_plf, A_plf, theta_tlf_cm, 0.45),
                    )
                    for branch, Z_prim, A_prim, Z_comp, A_comp, theta_cm, branch_w in branches:
                        E_star = self._excitation_energy(dZ, dN, l, branch)
                        for Z_fin, A_fin, evap_p in self._hivap_evaporation(Z_prim, A_prim, E_star):
                            if evap_p < 1e-6 or A_fin <= 0:
                                continue
                            A_comp_fin = self.Ap + self.At - A_fin
                            if A_comp_fin <= 0:
                                continue
                            E_loss = max(A_prim - A_fin, 0) * (7.5 + 2.0 * np.sqrt(E_star / max(A_prim / 8.0, 1.0)))
                            Q_evap = Q_primary - E_loss
                            E_lab, th_lab = self._two_body_kinematics(
                                float(A_fin), float(A_comp_fin), theta_cm, Q_evap
                            )
                            if E_lab <= 0 or not (0.0 <= th_lab <= 90.0):
                                continue

                            sigma_theta_eff = sigma_theta + 0.18 * max(A_prim - A_fin, 0)
                            sigma_E = sigma_E_base + 5.0 * max(A_prim - A_fin, 0)
                            weight = sig_ch * branch_w * evap_p
                            g_th = np.exp(-(theta_lab_grid - th_lab) ** 2 / (2.0 * sigma_theta_eff ** 2))
                            g_E = np.exp(-(E_grid - E_lab) ** 2 / (2.0 * sigma_E ** 2))
                            d2 += weight * np.outer(g_th, g_E)

        if np.sum(d2) <= 0 and total_sigma > 0:
            th0 = 45.0
            e0 = 0.5 * self.E_lab_total
            d2 += total_sigma * np.outer(
                np.exp(-(theta_lab_grid - th0) ** 2 / (2.0 * 10.0 ** 2)),
                np.exp(-(E_grid - e0) ** 2 / (2.0 * 80.0 ** 2)),
            )

        # Normalize
        total2 = np.sum(d2) * dE * dTh
        if total2 > 0:
            d2 *= total_sigma / total2

        self.result = MNTResult(Z=Zv, N=Nv, sigma_matrix=sigma_matrix,
                                sigma_total=total_sigma, E_cm=self.E_cm)

        return theta_lab_grid, E_grid, d2
