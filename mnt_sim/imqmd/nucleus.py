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
        self.grid_nuclear_scale = 1.0
        self.use_surface_term = True
        self.use_static_stabilizer = False

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
        return float(self.packets[0].sigma_r) if self.packets else 1.3

    @property
    def packet_sigmas(self) -> np.ndarray:
        """Per-packet wave packet widths (Wang 2002 Eq. (18): σ_r depends on A).

        Projectile and target nuclei are relaxed with their own
        σ_r = σ0 + σ1·A^(1/3); the combined system must keep the per-packet
        values (essential for asymmetric systems), so density deposition,
        EDF forces and Wigner kernels all consume this array.
        """

        return np.array([p.sigma_r for p in self.packets], dtype=float)

    def copy(self) -> "ImQMDNucleus":
        packets = [
            GaussianPacket(p.r_i.copy(), p.p_i.copy(), p.sigma_r, p.is_proton)
            for p in self.packets
        ]
        copied = ImQMDNucleus(self.Z, self.N, packets, self.edf, self.reference_positions, self.energy_offset)
        copied.grid_nuclear_scale = float(getattr(self, "grid_nuclear_scale", 1.0))
        copied.use_surface_term = bool(getattr(self, "use_surface_term", True))
        copied.use_static_stabilizer = bool(getattr(self, "use_static_stabilizer", False))
        return copied

    def density(self, r: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        """Evaluate folded neutron/proton densities at one or more points.

        Implements Wang et al. 2014 Eq. (1).
        """

        points = np.atleast_2d(np.asarray(r, dtype=float))
        pos = self.positions
        sigmas = self.packet_sigmas
        diff = points[:, None, :] - pos[None, :, :]
        # Per-packet Gaussian normalization/width (system-size-dependent σ_r).
        norm = 1.0 / ((2.0 * np.pi * sigmas**2) ** 1.5)
        weights = norm[None, :] * np.exp(-np.sum(diff * diff, axis=-1) / (2.0 * sigmas[None, :] ** 2))
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

    def energy_components(
        self,
        use_surface_term: bool | None = None,
        use_static_stabilizer: bool | None = None,
    ) -> dict[str, float]:
        """Return kinetic, potential, Coulomb, symmetry, and total energies.

        DEPRECATED (legacy centroid path): this compact centroid-sampled EDF
        (including the pair-form surface terms) is kept for static tests and
        diagnostics.  Production energies/forces use the GridEDF path
        (``total_energy`` / ``propagator._snapshot``), which is the
        self-consistent implementation.
        """

        if use_surface_term is None:
            use_surface_term = bool(getattr(self, "use_surface_term", True))
        if use_static_stabilizer is None:
            use_static_stabilizer = bool(getattr(self, "use_static_stabilizer", False))
        rho, rho_n, rho_p = self.centroid_densities()
        kinetic = self.edf.kinetic_energy(self.momenta, M_N)
        skyrme = self.edf.skyrme_potential(rho)
        surface = self.edf.surface_pair_energy(self.positions, self.sigma_r) if use_surface_term else 0.0
        surface_symmetry = (
            self.edf.surface_symmetry_energy(self.positions, self.is_proton, self.sigma_r)
            if use_surface_term
            else 0.0
        )
        coulomb_direct = self.edf.coulomb_energy(self.positions, self.is_proton, self.sigma_r)
        coulomb_exchange = self.edf.coulomb_exchange_energy(self.positions, self.is_proton, self.sigma_r)
        coulomb = coulomb_direct + coulomb_exchange
        symmetry = self.edf.symmetry_energy(rho_n, rho_p)
        static = (
            self.edf.static_mean_field_energy(self.positions, self.reference_positions)
            if use_static_stabilizer
            else 0.0
        )
        potential = skyrme + surface + surface_symmetry + coulomb + symmetry + static + self.energy_offset
        return {
            "kinetic": kinetic,
            "skyrme": skyrme,
            "surface": surface,
            "surface_symmetry": surface_symmetry,
            "coulomb": coulomb,
            "coulomb_direct": coulomb_direct,
            "coulomb_exchange": coulomb_exchange,
            "symmetry": symmetry,
            "static": static,
            "potential": potential,
            "total": kinetic + potential,
        }

    def total_energy(self) -> float:
        from .grid_edf import GridEDF

        grid = GridEDF(
            self.edf.parameters,
            self.packet_sigmas,
            grid_spacing=1.0,
            nuclear_scale=float(getattr(self, "grid_nuclear_scale", 1.0)),
        )
        e_grid = grid.total_energy(
            self.positions,
            self.momenta,
            self.is_proton,
            use_surface_term=self.use_surface_term,
        )["total"]
        total = e_grid + self.energy_offset
        if self.use_static_stabilizer:
            total += self.edf.static_mean_field_energy(self.positions, self.reference_positions)
        return float(total)

    def relative_energy_drift(self, reference_energy: float) -> float:
        """Return fractional total-energy drift from ``reference_energy``."""

        return abs(self.total_energy() - float(reference_energy)) / max(abs(float(reference_energy)), 1.0)

    def rms_radius(self) -> tuple[float, float]:
        """Return proton and neutron RMS radii, including packet width."""

        pos = self.positions
        sigmas = self.packet_sigmas
        mask = self.is_proton

        def rms(subset: np.ndarray, subset_sigmas: np.ndarray) -> float:
            if len(subset) == 0:
                return 0.0
            # Packet-width correction uses each packet's own σ_r.
            return float(np.sqrt(np.mean(np.sum(subset * subset, axis=1)) + 3.0 * np.mean(subset_sigmas**2)))

        return rms(pos[mask], sigmas[mask]), rms(pos[~mask], sigmas[~mask])

    def max_radius(self) -> float:
        return float(np.max(np.linalg.norm(self.positions, axis=1)))


__all__ = ["GaussianPacket", "ImQMDNucleus"]
