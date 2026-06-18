"""
Empirical systematics for MNT cross-section estimates.

Based on published experimental data fits:
  - Zagrebaev & Greiner, J. Phys. G 34 (2007) 2265
  - Feng et al., Phys. Rev. C 108 (2023) 024609

Useful for quick estimates and as sanity checks for full DNS calculations.
"""

import numpy as np
from dataclasses import dataclass
from typing import Dict, Tuple, Optional


@dataclass
class EmpiricalResult:
    """Empirical cross section estimates."""
    channels: Dict[str, float]  # e.g. "2p2n" -> sigma [mb]
    total_mnt: float            # mb
    description: str


class EmpiricalModel:
    """Fast empirical cross-section parameterization for MNT."""

    def __init__(self, Zp: int, Ap: int, Zt: int, At: int, E_lab: float):
        self.Zp, self.Ap = Zp, Ap
        self.Zt, self.At = Zt, At
        self.E_lab = E_lab

        mu = 931.494 * (Ap * At) / (Ap + At)
        self.E_cm = E_lab * Ap * At / (Ap + At)

    def _Q_value(self, delta_p: int, delta_n: int) -> float:
        """Approximate Q-value for (delta_p, delta_n) transfer."""
        return (delta_p * self.Zp * 0.6 + delta_n * self.Ap * 0.2
                - delta_p * delta_n * 0.1)

    def _Q_opt(self) -> float:
        """Optimal Q-value for maximum transfer probability."""
        return -2.0 * np.sqrt(self.E_cm / (self.Ap + self.At))

    def estimate(self,
                 channels: list = None) -> EmpiricalResult:
        """Estimate cross sections for specified transfer channels.

        Parameters
        ----------
        channels : list of (delta_p, delta_n) tuples
            Specific channels. If None, compute all up to ±4p ±6n.

        Returns
        -------
        EmpiricalResult with channel cross sections
        """
        if channels is None:
            channels = [(dp, dn) for dp in range(-4, 5)
                        for dn in range(-6, 7)]

        results = {}
        Z_eff = np.sqrt(self.Zp * self.Zt)
        A_eff = np.sqrt(self.Ap * self.At)

        # Base cross section from systematics
        sigma_0 = 100 * np.exp(-0.5 * (self.E_cm - 1.44 * self.Zp * self.Zt
                                        / (1.18 * (self.Ap**(1/3) + self.At**(1/3))))**2 / 100)

        for dp, dn in channels:
            A_tr = dp + dn
            if A_tr == 0:
                continue

            # Q-value matching
            Q = self._Q_value(dp, dn)
            Q_opt = self._Q_opt()
            Q_factor = np.exp(-(Q - Q_opt)**2 / (2 * 5**2))

            # Mass transfer suppression
            mass_factor = np.exp(-A_tr / (2 * Z_eff / 10))

            # Isospin dependence
            N_over_Z = (self.Ap - self.Zp) / self.Zp
            iso_factor = np.exp(-(dp / (dn + 0.01) - N_over_Z)**2 / (2 * 0.5**2))

            sigma = sigma_0 * Q_factor * mass_factor * iso_factor
            label = f"{'+' if dp >= 0 else ''}{dp}p{dn:+d}n"
            results[label] = max(sigma, 1e-6)

        total = sum(results.values())
        return EmpiricalResult(
            channels=results,
            total_mnt=total,
            description=(f"Empirical estimate for Z={self.Zp}A={self.Ap} "
                        f"on Z={self.Zt}A={self.At} at E_lab={self.E_lab} MeV/u")
        )
