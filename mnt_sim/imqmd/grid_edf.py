"""Grid-integrated Skyrme EDF for ImQMD Gaussian packets."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .skyrme import E2, M_N, PARAMETER_SETS, SkyrmeParameters


@dataclass(frozen=True)
class GridShape:
    """Cached 3D grid metadata."""

    x: np.ndarray
    y: np.ndarray
    z: np.ndarray

    @property
    def sizes(self) -> tuple[int, int, int]:
        return (len(self.x), len(self.y), len(self.z))


class GridEDF:
    """Skyrme EDF calculator using 3D spatial-grid integration.

    Densities are deposited from normalized Gaussian wave packets and truncated
    at ``n_sigma * sigma_r`` for the density sum. Local EDF terms are then
    integrated over the grid volume.
    """

    def __init__(
        self,
        parameters: SkyrmeParameters | str | None,
        sigma_r: float,
        grid_spacing: float = 1.0,
        n_sigma: float = 3,
        nuclear_scale: float = 1.0,
    ):
        """Create a grid EDF calculator.

        ``grid_spacing`` is in fm. ``n_sigma`` controls both grid padding and
        the Gaussian density cutoff radius.
        """

        if isinstance(parameters, str):
            parameters = PARAMETER_SETS[parameters.upper()]
        self.parameters = parameters or PARAMETER_SETS["IQ3A"]
        self.sigma_r = float(sigma_r)
        self.grid_spacing = float(grid_spacing)
        self.n_sigma = float(n_sigma)
        self.nuclear_scale = float(nuclear_scale)
        if self.sigma_r <= 0.0:
            raise ValueError("sigma_r must be positive")
        if self.grid_spacing <= 0.0:
            raise ValueError("grid_spacing must be positive")
        if self.n_sigma <= 0.0:
            raise ValueError("n_sigma must be positive")
        if self.nuclear_scale <= 0.0:
            raise ValueError("nuclear_scale must be positive")
        self._grid: GridShape | None = None

    def build(self, positions: np.ndarray) -> tuple[tuple[np.ndarray, np.ndarray, np.ndarray], tuple[int, int, int]]:
        """Auto-size the grid from nucleon positions and return axes and sizes."""

        pos = self._as_positions(positions)
        self._grid = self._build_grid(pos)
        return (self._grid.x, self._grid.y, self._grid.z), self._grid.sizes

    def density_on_grid(self, positions: np.ndarray, is_proton: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        """Return neutron and proton densities on the current grid."""

        pos = self._as_positions(positions)
        protons = np.asarray(is_proton, dtype=bool)
        if len(protons) != len(pos):
            raise ValueError("is_proton length must match positions")
        grid = self._require_grid(pos)
        rho_n = np.zeros(grid.sizes, dtype=float)
        rho_p = np.zeros(grid.sizes, dtype=float)

        for r_i, is_p in zip(pos, protons):
            ix, iy, iz, weights = self._density_contribution(grid, r_i)
            if weights.size == 0:
                continue
            if is_p:
                rho_p[ix, iy, iz] += weights
            else:
                rho_n[ix, iy, iz] += weights
        return rho_n, rho_p

    def gradients(self, rho: np.ndarray, grid_x: np.ndarray | None = None) -> np.ndarray:
        """Compute ``|grad rho|^2`` using finite differences on the grid."""

        del grid_x
        edge_order = 2 if min(np.shape(rho)) >= 3 else 1
        gx, gy, gz = np.gradient(np.asarray(rho, dtype=float), self.grid_spacing, edge_order=edge_order)
        return gx * gx + gy * gy + gz * gz

    def integrate(self, integrand_3d: np.ndarray) -> float:
        """Integrate a 3D integrand over the grid volume."""

        return float(np.sum(np.asarray(integrand_3d, dtype=float)) * self.grid_spacing**3)

    def skyrme_bulk(self, rho: np.ndarray, rho_n: np.ndarray, rho_p: np.ndarray) -> np.ndarray:
        """Return the bulk Skyrme energy-density integrand."""

        del rho_n, rho_p
        p = self.parameters
        rho = np.clip(np.asarray(rho, dtype=float), 0.0, None)
        return (
            p.alpha / (2.0 * p.rho0) * rho**2
            + p.beta / ((p.gamma + 1.0) * p.rho0**p.gamma) * rho ** (p.gamma + 1.0)
            + p.g_tau / (p.rho0**p.eta) * rho ** (p.eta + 1.0)
        )

    def symmetry_term(self, rho: np.ndarray, rho_n: np.ndarray, rho_p: np.ndarray) -> np.ndarray:
        """Return the volume symmetry energy-density integrand."""

        p = self.parameters
        rho = np.asarray(rho, dtype=float)
        rho_n = np.asarray(rho_n, dtype=float)
        rho_p = np.asarray(rho_p, dtype=float)
        return p.c_sym / (2.0 * p.rho0) * (rho_n - rho_p) ** 2

    def surface_term(self, grad_rho_sq: np.ndarray) -> np.ndarray:
        """Return ``g_sur/(2 rho0) * (grad rho)^2``."""

        return self.parameters.gsur / (2.0 * self.parameters.rho0) * np.asarray(grad_rho_sq, dtype=float)

    def surface_symmetry(
        self,
        grad_rho_sq: np.ndarray,
        rho: np.ndarray,
        rho_n: np.ndarray,
        rho_p: np.ndarray,
    ) -> np.ndarray:
        """Return ``C_s/(2 rho0) * [-kappa_s (grad rho)^2] * delta^2``."""

        p = self.parameters
        rho = np.asarray(rho, dtype=float)
        rho_n = np.asarray(rho_n, dtype=float)
        rho_p = np.asarray(rho_p, dtype=float)
        delta = np.zeros_like(rho)
        np.divide(rho_n - rho_p, rho, out=delta, where=rho > 1.0e-14)
        return p.c_sym / (2.0 * p.rho0) * (-p.kappa_s * np.asarray(grad_rho_sq, dtype=float)) * delta * delta

    def coulomb_exchange(self, rho_p: np.ndarray) -> np.ndarray:
        """Return Slater Coulomb exchange energy density without rho0 scaling."""

        rho_p = np.clip(np.asarray(rho_p, dtype=float), 0.0, None)
        return -0.75 * E2 * (3.0 / np.pi) ** (1.0 / 3.0) * rho_p ** (4.0 / 3.0)

    def total_energy(
        self,
        positions: np.ndarray,
        momenta: np.ndarray,
        is_proton: np.ndarray,
        use_surface_term: bool = True,
    ) -> dict[str, float]:
        """Return kinetic, grid-integrated EDF, Coulomb, and total energies."""

        pos = self._as_positions(positions)
        self.build(pos)
        rho_n, rho_p = self.density_on_grid(pos, is_proton)

        kinetic = self.kinetic_energy(momenta)
        components = self._potential_from_density(pos, is_proton, rho_n, rho_p, use_surface_term=use_surface_term)
        potential = components["potential"]
        total = kinetic + potential
        return {
            "kinetic": kinetic,
            **components,
            "potential": potential,
            "total": total,
        }

    def forces_analytical(
        self,
        positions: np.ndarray,
        is_proton: np.ndarray,
        epsilon: float = 1.0e-4,
        use_surface_term: bool = True,
    ) -> np.ndarray:
        """Return analytical grid EDF forces from one density-gradient pass."""

        del epsilon
        pos = self._as_positions(positions)
        protons = np.asarray(is_proton, dtype=bool)
        if len(protons) != len(pos):
            raise ValueError("is_proton length must match positions")

        self.build(pos)
        rho_n, rho_p = self.density_on_grid(pos, protons)
        rho = rho_n + rho_p

        p = self.parameters
        d_volume = self.grid_spacing ** 3

        x = rho / p.rho0
        bulk_prime = p.alpha * x + p.beta * x**p.gamma + p.g_tau * (p.eta + 1.0) * x**p.eta

        d_e_drho_n = self.nuclear_scale * bulk_prime
        d_e_drho_p = self.nuclear_scale * bulk_prime

        asym = rho_n - rho_p
        d_e_drho_n += self.nuclear_scale * (p.c_sym / p.rho0 * asym)
        d_e_drho_p -= self.nuclear_scale * (p.c_sym / p.rho0 * asym)

        if use_surface_term:
            lap = self._laplacian(rho)
            d_e_drho_n += self.nuclear_scale * (-(p.gsur / p.rho0) * lap)
            d_e_drho_p += self.nuclear_scale * (-(p.gsur / p.rho0) * lap)

        rho_p_safe = np.maximum(rho_p, 1.0e-12)
        d_e_drho_p += -E2 * (3.0 / np.pi) ** (1.0 / 3.0) * rho_p_safe ** (1.0 / 3.0)

        forces = np.zeros_like(pos, dtype=float)
        sigma2 = self.sigma_r**2
        grid = self._require_grid(pos)

        for i, (r_i, is_p) in enumerate(zip(pos, protons)):
            ix, iy, iz, rho_i = self._density_contribution(grid, r_i)
            if rho_i.size == 0:
                continue

            d_e_drho = d_e_drho_p if is_p else d_e_drho_n
            dx = grid.x[ix] - r_i[0]
            dy = grid.y[iy] - r_i[1]
            dz = grid.z[iz] - r_i[2]
            d_e_local = d_e_drho[ix, iy, iz]

            forces[i, 0] = -np.sum(d_e_local * (dx[:, None, None] / sigma2) * rho_i) * d_volume
            forces[i, 1] = -np.sum(d_e_local * (dy[None, :, None] / sigma2) * rho_i) * d_volume
            forces[i, 2] = -np.sum(d_e_local * (dz[None, None, :] / sigma2) * rho_i) * d_volume

        forces += self._coulomb_direct_forces(pos, protons)
        return forces

    def forces_numerical(
        self,
        positions: np.ndarray,
        is_proton: np.ndarray,
        epsilon: float = 1.0e-4,
        use_surface_term: bool = True,
    ) -> np.ndarray:
        """Return numerical EDF forces by differentiating integrated energy.

        The density grid is built once. For each coordinate shift only the
        moved nucleon's Gaussian contribution is replaced in the cached neutron
        or proton density, avoiding a full density rebuild for every force
        component.
        """

        pos = self._as_positions(positions)
        protons = np.asarray(is_proton, dtype=bool)
        if len(protons) != len(pos):
            raise ValueError("is_proton length must match positions")
        if epsilon <= 0.0:
            raise ValueError("epsilon must be positive")
        fixed_grid = self._build_grid(pos, extra_padding=epsilon)
        self._grid = fixed_grid
        rho_n, rho_p = self.density_on_grid(pos, protons)
        old_contributions = [self._density_contribution(fixed_grid, r_i) for r_i in pos]
        forces = np.zeros_like(pos, dtype=float)
        for i in range(len(pos)):
            old_ix, old_iy, old_iz, old_weights = old_contributions[i]
            for axis in range(3):
                e_plus = self._shifted_potential(
                    pos,
                    protons,
                    rho_n,
                    rho_p,
                    i,
                    axis,
                    epsilon,
                    old_ix,
                    old_iy,
                    old_iz,
                    old_weights,
                    use_surface_term=use_surface_term,
                )
                e_minus = self._shifted_potential(
                    pos,
                    protons,
                    rho_n,
                    rho_p,
                    i,
                    axis,
                    -epsilon,
                    old_ix,
                    old_iy,
                    old_iz,
                    old_weights,
                    use_surface_term=use_surface_term,
                )
                forces[i, axis] = -(e_plus - e_minus) / (2.0 * epsilon)
        self._grid = fixed_grid
        return forces

    def kinetic_energy(self, momenta: np.ndarray, mass: float = M_N) -> float:
        """Return ``sum p_i^2/(2m)`` in MeV for momenta in MeV/c."""

        p = np.asarray(momenta, dtype=float)
        return float(np.sum(p * p) / (2.0 * mass))

    def coulomb_direct(self, positions: np.ndarray, is_proton: np.ndarray) -> float:
        """Return finite-Gaussian proton-proton direct Coulomb energy."""

        pos = self._as_positions(positions)
        proton_pos = pos[np.asarray(is_proton, dtype=bool)]
        if len(proton_pos) < 2:
            return 0.0
        diff = proton_pos[:, None, :] - proton_pos[None, :, :]
        rij = np.linalg.norm(diff, axis=-1)
        iu = np.triu_indices(len(proton_pos), 1)
        r = np.maximum(rij[iu], 1.0e-12)
        try:
            from scipy.special import erf
        except Exception:  # pragma: no cover
            erf = np.vectorize(__import__("math").erf)
        return float(np.sum(E2 * erf(r / (2.0 * self.sigma_r)) / r))

    def _potential_from_density(
        self,
        positions: np.ndarray,
        is_proton: np.ndarray,
        rho_n: np.ndarray,
        rho_p: np.ndarray,
        use_surface_term: bool = True,
    ) -> dict[str, float]:
        rho = rho_n + rho_p
        skyrme_bulk = self.nuclear_scale * self.integrate(self.skyrme_bulk(rho, rho_n, rho_p))
        symmetry = self.nuclear_scale * self.integrate(self.symmetry_term(rho, rho_n, rho_p))
        if use_surface_term:
            grad_rho_sq = self.gradients(rho, None)
            surface = self.nuclear_scale * self.integrate(self.surface_term(grad_rho_sq))
            surface_sym = self.nuclear_scale * self.integrate(self.surface_symmetry(grad_rho_sq, rho, rho_n, rho_p))
        else:
            surface = 0.0
            surface_sym = 0.0
        coulomb_direct = self.coulomb_direct(positions, is_proton)
        coulomb_exchange = self.integrate(self.coulomb_exchange(rho_p))
        potential = skyrme_bulk + symmetry + surface + surface_sym + coulomb_direct + coulomb_exchange
        return {
            "skyrme_bulk": skyrme_bulk,
            "symmetry": symmetry,
            "surface": surface,
            "surface_sym": surface_sym,
            "surface_symmetry": surface_sym,
            "coulomb_direct": coulomb_direct,
            "coulomb_exchange": coulomb_exchange,
            "potential": potential,
        }

    def _density_functional_derivatives(
        self,
        rho_n: np.ndarray,
        rho_p: np.ndarray,
        use_surface_term: bool = True,
    ) -> tuple[np.ndarray, np.ndarray]:
        """Return unscaled ``dE/d rho_n`` and ``dE/d rho_p`` on the grid."""

        p = self.parameters
        rho = np.clip(rho_n + rho_p, 0.0, None)
        x = rho / p.rho0
        bulk_prime = p.alpha * x + p.beta * x**p.gamma + p.g_tau * (p.eta + 1.0) * x**p.eta

        d_n = bulk_prime.copy()
        d_p = bulk_prime.copy()

        asym = rho_n - rho_p
        d_n += p.c_sym / p.rho0 * asym
        d_p -= p.c_sym / p.rho0 * asym

        if use_surface_term:
            lap = self._laplacian(rho)
            d_n += -(p.gsur / p.rho0) * lap
            d_p += -(p.gsur / p.rho0) * lap

        d_p += -E2 * (3.0 / np.pi) ** (1.0 / 3.0) * np.maximum(rho_p, 1.0e-12) ** (1.0 / 3.0)
        return d_n, d_p

    def _gradient_components(self, rho: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        edge_order = 2 if min(np.shape(rho)) >= 3 else 1
        return np.gradient(np.asarray(rho, dtype=float), self.grid_spacing, edge_order=edge_order)

    def _laplacian(self, rho: np.ndarray) -> np.ndarray:
        """Second-order Laplacian of a 3D array on the uniform grid."""

        f = np.asarray(rho, dtype=float)
        edge_order = 2 if min(f.shape) >= 3 else 1
        gx, gy, gz = np.gradient(f, self.grid_spacing, edge_order=edge_order)
        return (
            np.gradient(gx, self.grid_spacing, axis=0, edge_order=2 if f.shape[0] >= 3 else 1)
            + np.gradient(gy, self.grid_spacing, axis=1, edge_order=2 if f.shape[1] >= 3 else 1)
            + np.gradient(gz, self.grid_spacing, axis=2, edge_order=2 if f.shape[2] >= 3 else 1)
        )

    def _divergence(self, vx: np.ndarray, vy: np.ndarray, vz: np.ndarray) -> np.ndarray:
        return (
            np.gradient(vx, self.grid_spacing, axis=0, edge_order=2 if vx.shape[0] >= 3 else 1)
            + np.gradient(vy, self.grid_spacing, axis=1, edge_order=2 if vy.shape[1] >= 3 else 1)
            + np.gradient(vz, self.grid_spacing, axis=2, edge_order=2 if vz.shape[2] >= 3 else 1)
        )

    def _coulomb_direct_forces(self, positions: np.ndarray, is_proton: np.ndarray) -> np.ndarray:
        """Return unscaled finite-Gaussian direct Coulomb forces."""

        forces = np.zeros_like(positions, dtype=float)
        proton_indices = np.flatnonzero(is_proton)
        if len(proton_indices) < 2:
            return forces
        proton_pos = positions[proton_indices]
        i_local, j_local = np.triu_indices(len(proton_pos), 1)
        diff = proton_pos[i_local] - proton_pos[j_local]
        r = np.linalg.norm(diff, axis=1)
        r_safe = np.maximum(r, 1.0e-12)
        try:
            from scipy.special import erf
        except Exception:  # pragma: no cover
            erf = np.vectorize(__import__("math").erf)
        x = r_safe / (2.0 * self.sigma_r)
        d_erf_over_r = np.exp(-x * x) / (np.sqrt(np.pi) * self.sigma_r * r_safe) - erf(x) / (r_safe * r_safe)
        pair_forces = -E2 * d_erf_over_r[:, None] * diff / r_safe[:, None]
        np.add.at(forces, proton_indices[i_local], pair_forces)
        np.add.at(forces, proton_indices[j_local], -pair_forces)
        return forces

    def _shifted_potential(
        self,
        positions: np.ndarray,
        is_proton: np.ndarray,
        rho_n: np.ndarray,
        rho_p: np.ndarray,
        index: int,
        axis: int,
        shift: float,
        old_ix: slice,
        old_iy: slice,
        old_iz: slice,
        old_weights: np.ndarray,
        use_surface_term: bool = True,
    ) -> float:
        shifted_positions = positions.copy()
        shifted_positions[index, axis] += shift
        new_ix, new_iy, new_iz, new_weights = self._density_contribution(self._require_grid(positions), shifted_positions[index])
        shifted_rho_n = rho_n.copy()
        shifted_rho_p = rho_p.copy()
        shifted_rho = shifted_rho_p if is_proton[index] else shifted_rho_n
        if old_weights.size:
            shifted_rho[old_ix, old_iy, old_iz] -= old_weights
        if new_weights.size:
            shifted_rho[new_ix, new_iy, new_iz] += new_weights
        return self._potential_from_density(
            shifted_positions,
            is_proton,
            shifted_rho_n,
            shifted_rho_p,
            use_surface_term=use_surface_term,
        )["potential"]

    def _require_grid(self, positions: np.ndarray) -> GridShape:
        if self._grid is None:
            self.build(positions)
        assert self._grid is not None
        return self._grid

    def _build_grid(self, positions: np.ndarray, extra_padding: float = 0.0) -> GridShape:
        pad = self.n_sigma * self.sigma_r + float(extra_padding)
        lower = np.min(positions, axis=0) - pad
        upper = np.max(positions, axis=0) + pad
        axes = tuple(
            np.arange(lower[axis], upper[axis] + 0.5 * self.grid_spacing, self.grid_spacing, dtype=float)
            for axis in range(3)
        )
        return GridShape(*axes)

    def _axis_window(self, axis: np.ndarray, center: float, cutoff: float) -> slice:
        start = int(np.searchsorted(axis, center - cutoff, side="left"))
        stop = int(np.searchsorted(axis, center + cutoff, side="right"))
        return slice(max(start, 0), min(stop, len(axis)))

    def _density_contribution(
        self,
        grid: GridShape,
        r_i: np.ndarray,
    ) -> tuple[slice, slice, slice, np.ndarray]:
        sigma2 = self.sigma_r * self.sigma_r
        norm = 1.0 / ((2.0 * np.pi * sigma2) ** 1.5)
        cutoff = self.n_sigma * self.sigma_r
        cutoff2 = cutoff * cutoff
        ix = self._axis_window(grid.x, float(r_i[0]), cutoff)
        iy = self._axis_window(grid.y, float(r_i[1]), cutoff)
        iz = self._axis_window(grid.z, float(r_i[2]), cutoff)
        if ix.stop <= ix.start or iy.stop <= iy.start or iz.stop <= iz.start:
            return ix, iy, iz, np.zeros((0, 0, 0), dtype=float)

        dx2 = (grid.x[ix] - r_i[0]) ** 2
        dy2 = (grid.y[iy] - r_i[1]) ** 2
        dz2 = (grid.z[iz] - r_i[2]) ** 2
        r2 = dx2[:, None, None] + dy2[None, :, None] + dz2[None, None, :]
        weights = norm * np.exp(-r2 / (2.0 * sigma2))
        return ix, iy, iz, np.where(r2 <= cutoff2, weights, 0.0)

    def _as_positions(self, positions: np.ndarray) -> np.ndarray:
        pos = np.asarray(positions, dtype=float)
        if pos.ndim != 2 or pos.shape[1] != 3:
            raise ValueError("positions must have shape (A, 3)")
        if len(pos) == 0:
            raise ValueError("positions must contain at least one nucleon")
        return pos


__all__ = ["GridEDF", "GridShape"]
