"""Nucleus containers and Gaussian density evaluation for ImQMD."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .skyrme import M_N, SkyrmeEDF


@dataclass
class GaussianPacket:
    """Centroid of a Gaussian wave packet."""

    r_i: np.ndarray
    p_i: np.ndarray
    sigma_r: float
    is_proton: bool


class ImQMDNucleus:
    """A static finite nucleus represented by Gaussian nucleon packets."""

    def __init__(
        self,
        Z: int,
        N: int,
        packets: list[GaussianPacket],
        edf: SkyrmeEDF | None = None,
        reference_positions: np.ndarray | None = None,
        energy_offset: float = 0.0,
    ):
        self.Z = int(Z)
        self.N = int(N)
        self.A = self.Z + self.N
        self.packets = packets
        self.edf = edf or SkyrmeEDF()
        self.reference_positions = None if reference_positions is None else np.asarray(reference_positions, dtype=float).copy()
        self.energy_offset = float(energy_offset)

    @property
    def positions(self) -> np.ndarray:
        return np.array([p.r_i for p in self.packets], dtype=float)

    @positions.setter
    def positions(self, value: np.ndarray) -> None:
        for packet, r_i in zip(self.packets, np.asarray(value, dtype=float)):
            packet.r_i = np.array(r_i, dtype=float)

    @property
    def momenta(self) -> np.ndarray:
        return np.array([p.p_i for p in self.packets], dtype=float)

    @momenta.setter
    def momenta(self, value: np.ndarray) -> None:
        for packet, p_i in zip(self.packets, np.asarray(value, dtype=float)):
            packet.p_i = np.array(p_i, dtype=float)

    @property
    def is_proton(self) -> np.ndarray:
        return np.array([p.is_proton for p in self.packets], dtype=bool)

    @property
    def sigma_r(self) -> float:
        return float(self.packets[0].sigma_r) if self.packets else 1.1

    def copy(self) -> "ImQMDNucleus":
        packets = [
            GaussianPacket(p.r_i.copy(), p.p_i.copy(), p.sigma_r, p.is_proton)
            for p in self.packets
        ]
        return ImQMDNucleus(self.Z, self.N, packets, self.edf, self.reference_positions, self.energy_offset)

    def density(self, r: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        """Evaluate folded neutron/proton densities at one or more points.

        Implements Wang et al. 2014 Eq. (1).
        """

        points = np.atleast_2d(np.asarray(r, dtype=float))
        pos = self.positions
        diff = points[:, None, :] - pos[None, :, :]
        norm = 1.0 / ((2.0 * np.pi * self.sigma_r**2) ** 1.5)
        weights = norm * np.exp(-np.sum(diff * diff, axis=-1) / (2.0 * self.sigma_r**2))
        protons = self.is_proton
        rho_p = np.sum(weights[:, protons], axis=1)
        rho_n = np.sum(weights[:, ~protons], axis=1)
        if np.asarray(r).ndim == 1:
            return rho_n[0], rho_p[0]
        return rho_n, rho_p

    def centroid_densities(self) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Return total, neutron, and proton densities sampled at centroids."""

        rho_n, rho_p = self.density(self.positions)
        rho = rho_n + rho_p
        return rho, rho_n, rho_p

    def energy_components(self) -> dict[str, float]:
        """Return kinetic, potential, Coulomb, symmetry, and total energies."""

        rho, rho_n, rho_p = self.centroid_densities()
        kinetic = self.edf.kinetic_energy(self.momenta, M_N)
        skyrme = self.edf.skyrme_potential(rho)
        coulomb = self.edf.coulomb_energy(self.positions, self.is_proton, self.sigma_r)
        symmetry = self.edf.symmetry_energy(rho_n, rho_p)
        static = self.edf.static_mean_field_energy(self.positions, self.reference_positions)
        potential = skyrme + coulomb + symmetry + static + self.energy_offset
        return {
            "kinetic": kinetic,
            "skyrme": skyrme,
            "coulomb": coulomb,
            "symmetry": symmetry,
            "static": static,
            "potential": potential,
            "total": kinetic + potential,
        }

    def total_energy(self) -> float:
        return self.energy_components()["total"]

    def rms_radius(self) -> tuple[float, float]:
        """Return proton and neutron RMS radii, including packet width."""

        pos = self.positions
        mask = self.is_proton

        def rms(subset: np.ndarray) -> float:
            if len(subset) == 0:
                return 0.0
            return float(np.sqrt(np.mean(np.sum(subset * subset, axis=1)) + 3.0 * self.sigma_r**2))

        return rms(pos[mask]), rms(pos[~mask])

    def max_radius(self) -> float:
        return float(np.max(np.linalg.norm(self.positions, axis=1)))


__all__ = ["GaussianPacket", "ImQMDNucleus"]
