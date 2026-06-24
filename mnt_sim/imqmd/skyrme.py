"""Skyrme-like energy-density utilities for the Phase 1 ImQMD model.

Formula references are summarized in ``papers/imqmd_formulas.md``.  The
implemented energy is intentionally compact: it exposes the requested EDF terms
and uses finite Gaussian overlaps for static nuclei.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

HBAR_C = 197.327
M_N = 938.9
E2 = 1.43996448


@dataclass(frozen=True)
class SkyrmeParameters:
    """ImQMD Skyrme parameter set."""

    name: str
    alpha: float
    beta: float
    gamma: float
    gsur: float
    g_tau: float
    eta: float
    c_sym: float
    kappa_s: float
    rho0: float
    sigma0: float
    sigma1: float


PARAMETER_SETS = {
    "IQ1": SkyrmeParameters("IQ1", -310.0, 258.0, 7.0 / 6.0, 19.8, 9.5, 2.0 / 3.0, 32.0, 0.08, 0.165, 0.49, 0.16),
    "IQ2": SkyrmeParameters("IQ2", -356.0, 303.0, 7.0 / 6.0, 7.0, 12.5, 2.0 / 3.0, 32.0, 0.08, 0.165, 0.88, 0.09),
    "IQ3": SkyrmeParameters("IQ3", -207.0, 138.0, 7.0 / 6.0, 18.0, 14.0, 5.0 / 3.0, 32.0, 0.08, 0.165, 0.94, 0.018),
    "IQ3A": SkyrmeParameters("IQ3a", -207.0, 138.0, 7.0 / 6.0, 16.5, 14.0, 5.0 / 3.0, 34.0, 0.4, 0.165, 0.94, 0.020),
    "IQ3B": SkyrmeParameters("IQ3b", -207.0, 138.0, 7.0 / 6.0, 18.0, 14.0, 5.0 / 3.0, 34.0, 0.6, 0.165, 0.94, 0.018),
}


class SkyrmeEDF:
    """Small ImQMD EDF calculator with IQ2 defaults.

    The local terms follow Wang et al. 2014 Eq. (5).  In particle form the code
    evaluates Gaussian overlap densities at centroids, which is sufficient for
    Phase 1 ground-state acceptance and diagnostics.
    """

    def __init__(self, parameters: SkyrmeParameters | None = None, static_k: float = 5.0):
        self.parameters = parameters or PARAMETER_SETS["IQ2"]
        self.static_k = static_k

    @classmethod
    def from_name(cls, name: str, static_k: float = 5.0) -> "SkyrmeEDF":
        return cls(PARAMETER_SETS[name.upper()], static_k=static_k)

    def kinetic_energy(self, momenta: np.ndarray, mass: float = M_N) -> float:
        """Return ``sum_i p_i^2/(2m)`` in MeV for momenta in MeV/c."""

        p = np.asarray(momenta, dtype=float)
        return float(np.sum(p * p) / (2.0 * mass))

    def skyrme_potential(self, rho: np.ndarray) -> float:
        """Return local bulk Skyrme potential sampled at centroid densities."""

        p = self.parameters
        rho = np.clip(np.asarray(rho, dtype=float), 1.0e-12, None)
        u2 = p.alpha / (2.0 * p.rho0) * np.sum(rho)
        u3 = p.beta / ((p.gamma + 1.0) * p.rho0**p.gamma) * np.sum(rho**p.gamma)
        utau = p.g_tau / (p.rho0**p.eta) * np.sum(rho**p.eta)
        return float(u2 + u3 + utau)

    def coulomb_energy(self, positions: np.ndarray, is_proton: np.ndarray, sigma_r: float) -> float:
        """Finite-Gaussian proton-proton Coulomb energy in MeV."""

        pos = np.asarray(positions, dtype=float)
        proton_pos = pos[np.asarray(is_proton, dtype=bool)]
        if len(proton_pos) < 2:
            return 0.0
        diff = proton_pos[:, None, :] - proton_pos[None, :, :]
        rij = np.linalg.norm(diff, axis=-1)
        iu = np.triu_indices(len(proton_pos), 1)
        r = np.maximum(rij[iu], 1.0e-9)
        # Smooth Gaussian charge packets; erf is imported lazily to keep scipy optional.
        try:
            from scipy.special import erf
        except Exception:  # pragma: no cover
            erf = np.vectorize(np.math.erf)
        return float(np.sum(E2 * erf(r / (2.0 * sigma_r)) / r))

    def symmetry_energy(self, rho_n: np.ndarray, rho_p: np.ndarray) -> float:
        """Return a centroid-sampled symmetry term from Wang et al. 2014 Eq. (5)."""

        p = self.parameters
        rho_n = np.asarray(rho_n, dtype=float)
        rho_p = np.asarray(rho_p, dtype=float)
        rho = np.clip(rho_n + rho_p, 1.0e-12, None)
        delta = (rho_n - rho_p) / rho
        return float(p.c_sym / (2.0 * p.rho0) * np.sum(rho * delta * delta))

    def total_energy_density(self, rho: np.ndarray, rho_n: np.ndarray, rho_p: np.ndarray) -> float:
        """Return local EDF contribution excluding kinetic and Coulomb terms."""

        return self.skyrme_potential(rho) + self.symmetry_energy(rho_n, rho_p)

    def static_mean_field_energy(self, positions: np.ndarray, reference_positions: np.ndarray | None) -> float:
        """Weak harmonic static mean field used to test isolated nucleus stability."""

        if reference_positions is None or self.static_k <= 0.0:
            return 0.0
        dr = np.asarray(positions, dtype=float) - np.asarray(reference_positions, dtype=float)
        return float(0.5 * self.static_k * np.sum(dr * dr))

    def static_mean_field_force(self, positions: np.ndarray, reference_positions: np.ndarray | None) -> np.ndarray:
        """Force from the weak static mean field, ``-dU/dr``."""

        if reference_positions is None or self.static_k <= 0.0:
            return np.zeros_like(positions, dtype=float)
        return -self.static_k * (np.asarray(positions, dtype=float) - np.asarray(reference_positions, dtype=float))


__all__ = ["E2", "HBAR_C", "M_N", "PARAMETER_SETS", "SkyrmeEDF", "SkyrmeParameters"]
