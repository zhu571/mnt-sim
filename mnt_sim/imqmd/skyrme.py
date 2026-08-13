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
    w_p: float = 60.0


PARAMETER_SETS = {
    "IQ1": SkyrmeParameters("IQ1", -310.0, 258.0, 7.0 / 6.0, 19.8, 9.5, 2.0 / 3.0, 32.0, 0.08, 0.165, 0.49, 0.16),
    "IQ2": SkyrmeParameters("IQ2", -356.0, 303.0, 7.0 / 6.0, 7.0, 12.5, 2.0 / 3.0, 32.0, 0.08, 0.165, 0.88, 0.09),
    "IQ3": SkyrmeParameters("IQ3", -207.0, 138.0, 7.0 / 6.0, 18.0, 14.0, 5.0 / 3.0, 32.0, 0.08, 0.165, 0.94, 0.018),
    # IQ3a/IQ3b: variants of IQ3 with enhanced surface-symmetry for isospin-sensitive
    # observables (Li+2013 Chin. Phys. C 37, 114101 Sec. III; Yao+Wang 2017 PRC 95, 014607;
    # Li+2019 PRC 99, 034619 uses IQ3a as default for MNT).  Compared to IQ3 (κs=0.08, C_s=32):
    #   IQ3a: stronger surface-symmetry (κs=0.4) and surface tension (~8% softer, gsur=16.5)
    #         for improved neutron-skin and isospin-drift dynamics in near-barrier MNT;
    #         sigma1=0.020 gives ~1.06 fm width for 238U.
    #   IQ3b: stiffer surface-symmetry (κs=0.6) variant, used for sensitivity studies.
    "IQ3A": SkyrmeParameters("IQ3a", -207.0, 138.0, 7.0 / 6.0, 16.5, 14.0, 5.0 / 3.0, 34.0, 0.4, 0.165, 0.94, 0.020),
    "IQ3B": SkyrmeParameters("IQ3b", -207.0, 138.0, 7.0 / 6.0, 18.0, 14.0, 5.0 / 3.0, 34.0, 0.6, 0.165, 0.94, 0.018),
    "SKP*": SkyrmeParameters("SkP*", -356.0, 303.0, 7.0 / 6.0, 19.5, 13.0, 2.0 / 3.0, 35.0, 0.65, 0.162, 0.94, 0.018),
}


class SkyrmeEDF:
    """Small ImQMD EDF calculator with IQ3a defaults.

    The local terms follow Wang et al. 2014 Eq. (5).  In particle form the code
    evaluates Gaussian overlap densities at centroids, which is sufficient for
    Phase 1 ground-state acceptance and diagnostics.
    """

    def __init__(self, parameters: SkyrmeParameters | None = None, static_k: float = 5.0):
        self.parameters = parameters or PARAMETER_SETS["IQ3A"]
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
        u2 = p.alpha / (2.0 * p.rho0) * np.sum(rho**2)
        u3 = p.beta / ((p.gamma + 1.0) * p.rho0**p.gamma) * np.sum(rho**(p.gamma + 1.0))
        utau = p.g_tau / (p.rho0**p.eta) * np.sum(rho**(p.eta + 1.0))
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
            erf = np.vectorize(__import__("math").erf)
        return float(np.sum(E2 * erf(r / (2.0 * sigma_r)) / r))

    def coulomb_exchange_energy(self, positions: np.ndarray, is_proton: np.ndarray, sigma_r: float) -> float:
        """Return centroid-sampled Slater Coulomb exchange energy."""

        pos = np.asarray(positions, dtype=float)
        protons = np.asarray(is_proton, dtype=bool)
        proton_pos = pos[protons]
        if len(proton_pos) == 0:
            return 0.0
        sigma2 = float(sigma_r) ** 2
        diff = proton_pos[:, None, :] - proton_pos[None, :, :]
        weights = np.exp(-np.sum(diff * diff, axis=-1) / (2.0 * sigma2))
        norm = 1.0 / ((2.0 * np.pi * sigma2) ** 1.5)
        rho_p = np.clip(norm * np.sum(weights, axis=1), 1.0e-12, None)
        coeff = -0.75 * E2 * (3.0 / np.pi) ** (1.0 / 3.0)
        return float(coeff * np.sum(rho_p ** (4.0 / 3.0)) / self.parameters.rho0)

    def symmetry_energy(self, rho_n: np.ndarray, rho_p: np.ndarray) -> float:
        """Return a centroid-sampled symmetry term from Wang et al. 2014 Eq. (5)."""

        p = self.parameters
        rho_n = np.asarray(rho_n, dtype=float)
        rho_p = np.asarray(rho_p, dtype=float)
        return float(p.c_sym / (2.0 * p.rho0) * np.sum((rho_n - rho_p) ** 2))

    def surface_pair_energy(
        self,
        positions: np.ndarray,
        sigma_r: float,
        group_ids: np.ndarray | None = None,
    ) -> float:
        """Return a compact finite-range surface cohesion term.

        ``group_ids`` is accepted for backward compatibility, but surface
        cohesion applies between all nearby nucleons regardless of origin.
        """

        pos = np.asarray(positions, dtype=float)
        if len(pos) < 2:
            return 0.0
        diff = pos[:, None, :] - pos[None, :, :]
        weights = np.exp(-np.sum(diff * diff, axis=-1) / (4.0 * sigma_r * sigma_r))
        iu = np.triu_indices(len(pos), 1)
        strength = self.parameters.gsur / 7.0
        return float(-strength * np.sum(weights[iu]))

    def surface_symmetry_energy(self, positions: np.ndarray, is_proton: np.ndarray, sigma_r: float) -> float:
        """Finite-range approximation to the kappa_s surface-symmetry term."""

        pos = np.asarray(positions, dtype=float)
        if len(pos) < 2:
            return 0.0
        tau = np.where(np.asarray(is_proton, dtype=bool), -1.0, 1.0)
        diff = pos[:, None, :] - pos[None, :, :]
        weights = np.exp(-np.sum(diff * diff, axis=-1) / (4.0 * sigma_r * sigma_r))
        iu = np.triu_indices(len(pos), 1)
        strength = self.parameters.c_sym * self.parameters.kappa_s / 7.0
        return float(-strength * np.sum(weights[iu] * tau[iu[0]] * tau[iu[1]]))

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
