"""
GRAZING-like semi-classical model for multi-nucleon transfer reactions.

Based on the semi-classical formalism of:
  - Winther, Nucl. Phys. A 572 (1994) 191
  - Corradi et al., J. Phys. G 36 (2009) 113101

Uses classical trajectories + quantum tunneling for nucleon exchange.
"""

import numpy as np
from dataclasses import dataclass
from typing import Tuple

__all__ = ["GrazingModel", "GrazingResult"]


@dataclass
class GrazingResult:
    """Container for grazing model results."""
    Z: np.ndarray
    N: np.ndarray
    sigma_matrix: np.ndarray  # mb
    theta_graz: float          # deg
    E_cm: float                # MeV


class GrazingModel:
    """Semi-classical GRAZING-like model.

    Uses classical deflection functions combined with
    quantal nucleon exchange probabilities.
    """

    def __init__(self, projectile: str, target: str, E_lab: float):
        from .dns import DNSModel
        self._dns = DNSModel(projectile, target, E_lab)
        self.Zp = self._dns.Zp
        self.Zt = self._dns.Zt
        self.Ap = self._dns.Ap
        self.At = self._dns.At
        self.E_lab = E_lab
        self.E_cm = self._dns.E_cm
        self.V_C = self._dns._compute_coulomb_barrier()

    def _interaction_time(self, l: int) -> float:
        """Interaction time [10^-22 s] for partial wave l."""
        hbar = 6.58e-22  # MeV·s
        mu = self._dns.M_N * (self.Ap * self.At) / (self.Ap + self.At)

        # Sommerfeld parameter
        eta = 1.44 * self.Zp * self.Zt / (hbar * np.sqrt(2 * self.E_cm / mu))

        # Distance of closest approach for Coulomb trajectory
        d = eta / np.sqrt(self.E_cm / mu) * (1 + 1/np.sin(np.pi/2))
        v_inf = np.sqrt(2 * self.E_cm / mu)

        # Interaction time from distance / velocity
        return 2 * d / v_inf * 1e22  # in 10^-22 s

    def _excitation_energy(self, l: int) -> float:
        """Mean excitation energy [MeV] for a given partial wave."""
        E_star = self.E_cm - self.V_C
        return max(0.1, E_star * np.exp(-l / 20.0))

    def calculate(self, delta_Z_range: Tuple[int, int] = (-6, 6),
                  delta_N_range: Tuple[int, int] = (-8, 8),
                  n_l: int = 100) -> GrazingResult:
        """Calculate transfer probabilities using GRAZING formalism.

        For each partial wave l, calculate nucleon exchange probability
        weighted by the transmission coefficient T_l.
        """
        Z_vals = np.arange(delta_Z_range[0], delta_Z_range[1] + 1)
        N_vals = np.arange(delta_N_range[0], delta_N_range[1] + 1)
        nZ, nN = len(Z_vals), len(N_vals)
        if nZ == 0 or nN == 0:
            raise ValueError("transfer ranges must not be empty")
        if n_l <= 0:
            raise ValueError("n_l must be positive")

        sigma_ch = np.zeros((nZ, nN))

        # Grazing angular momentum
        hbar = 6.58e-22  # MeV·s
        mu = self._dns.M_N * (self.Ap * self.At) / (self.Ap + self.At)
        eta = 1.44 * self.Zp * self.Zt / (hbar * np.sqrt(2 * self.E_cm / mu))

        l_graz = eta * np.sqrt(2 * self.E_cm / self.V_C - 1)

        k = np.sqrt(2 * mu * self.E_cm) / hbar
        R_int = 1.18 * (self.Ap**(1/3) + self.At**(1/3))

        for l in range(n_l):
            # Transmission coefficient (simplified Hill-Wheeler)
            E_l = hbar**2 * l * (l + 1) / (2 * mu * R_int**2)
            T_l = 1.0 / (1 + np.exp(2 * np.pi * (self.V_C + E_l - self.E_cm)
                                     / (4.0)))  # hbar_omega ≈ 4 MeV

            if T_l < 1e-6:
                continue

            # Excitation energy for this l
            E_star = self._excitation_energy(l)

            # Temperature
            T_temp = np.sqrt(8.0 * E_star / (self.Ap + self.At))

            # Transfer probability - Gaussian in Z,N space
            # Width grows with interaction time and temperature
            sigma_tr = 1.0 + 0.5 * np.sqrt(E_star)
            P_l = np.zeros((nZ, nN))

            for i, dZ in enumerate(Z_vals):
                for j, dN in enumerate(N_vals):
                    # Q-value dependence
                    Q = (dZ * 6.0 + dN * 8.0) * np.exp(-abs(dZ + dN) / 4.0)
                    P_l[i, j] = np.exp(-(dZ**2 + dN**2) / (2 * sigma_tr**2)
                                        - Q / (2 * T_temp) if Q > 0
                                        else -(dZ**2 + dN**2) / (2 * sigma_tr**2))

            # Normalize
            P_sum = np.sum(P_l)
            if P_sum > 0:
                P_l /= P_sum

            # Cross section contribution
            # dσ/dΩ = (λ²/4π) * (2l+1) * |S_l|²
            wavelength_sq = (hbar / np.sqrt(2 * mu * self.E_cm))**2
            sigma_l = (np.pi * wavelength_sq * (2 * l + 1) * T_l) * 1e3  # mb

            sigma_ch += sigma_l * P_l

        return GrazingResult(
            Z=Z_vals, N=N_vals,
            sigma_matrix=sigma_ch,
            theta_graz=np.degrees(2 * np.arctan(eta / (l_graz + 0.5))),
            E_cm=self.E_cm
        )
