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

    ``sigma_r`` may be a scalar (all packets share one width) or a per-packet
    array.  The per-packet form implements the system-size-dependent width
    sigma_r = sigma0 + sigma1 * A^(1/3) of Wang 2002 Eq. (18) for asymmetric
    projectile-target systems, where a single width is not well defined.

    Energy/force self-consistency: the density-dependent calibration scale
    (``_effective_scale``) is applied as a LOCAL factor inside the energy
    integrals, and the analytical forces differentiate exactly that functional
    (including the ds/drho terms and the kappa_s surface-symmetry term, whose
    force was previously missing).  ``forces_numerical`` differentiates the
    same integrated energy, so the two force paths agree to finite-difference
    accuracy.
    """

    def __init__(
        self,
        parameters: SkyrmeParameters | str | None,
        sigma_r: float | np.ndarray,
        grid_spacing: float = 1.0,
        n_sigma: float = 3,
        nuclear_scale: float = 1.0,
    ):
        """Create a grid EDF calculator.

        ``grid_spacing`` is in fm. ``n_sigma`` controls both grid padding and
        the Gaussian density cutoff radius.  ``sigma_r`` accepts a scalar or a
        per-packet array (length must match the positions passed later).
        """

        if isinstance(parameters, str):
            parameters = PARAMETER_SETS[parameters.upper()]
        self.parameters = parameters or PARAMETER_SETS["IQ3A"]
        sig = np.asarray(sigma_r, dtype=float)
        if sig.ndim == 0:
            self._sigma_array: np.ndarray | None = None
            self.sigma_r = float(sig)
        else:
            # Per-packet widths (Wang 2002 Eq. (18), one value per nucleon).
            self._sigma_array = sig.ravel()
            self.sigma_r = float(np.max(self._sigma_array))
        self.grid_spacing = float(grid_spacing)
        self.n_sigma = float(n_sigma)
        self.nuclear_scale = float(nuclear_scale)
        if np.any(sig <= 0.0):
            raise ValueError("sigma_r must be positive")
        if self.grid_spacing <= 0.0:
            raise ValueError("grid_spacing must be positive")
        if self.n_sigma <= 0.0:
            raise ValueError("n_sigma must be positive")
        if self.nuclear_scale < 0.0:
            raise ValueError("nuclear_scale must be non-negative")
        self._grid: GridShape | None = None
        self._grad_matrix_cache: dict[int, np.ndarray] = {}

    def _sigmas_for(self, n: int) -> np.ndarray:
        """Return the per-packet width array for ``n`` nucleons."""

        if self._sigma_array is None:
            return np.full(n, self.sigma_r, dtype=float)
        if len(self._sigma_array) != n:
            raise ValueError("per-packet sigma_r array length must match positions")
        return self._sigma_array

    def _effective_scale(self, rho: np.ndarray) -> np.ndarray:
        """Density-dependent scale: nuclear_scale at ρ≤ρ₀, shifts toward 1.0 above.

        At collision overlap densities (ρ>ρ₀) the repulsive EOS must not be
        softened by the nuclear_scale calibration.  Only ρ>ρ₀ is affected;
        the surface and bulk at normal density keep the calibrated scale.

        The ramp is deliberately mild (slope 4, centred at 2ρ₀): the previous
        steep ramp (slope 8 at 1.5ρ₀) restored full-strength repulsion already
        in the first overlap stage.  Combined with the weak calibrated
        attraction at ρ≤ρ₀ (nuclear_scale≈0.65 for 238U), that stiff bounce
        sprayed nucleons faster than the softened field could rebind them —
        a driver of the spurious U+U multifragmentation.  Centring the
        transition at 2ρ₀ keeps the density range actually probed in
        near-barrier reactions (ρ≲2ρ₀ at E/A≈7 MeV) close to the calibrated
        EOS, so the density dependence during the collision stage is gentler.

        This LOCAL factor is part of the energy functional itself (see
        ``_potential_from_density``); the analytical forces differentiate it,
        so energy and forces stay self-consistent.
        """
        x = np.clip(np.asarray(rho, dtype=float) / self.parameters.rho0, 1.0, None)
        # Sigmoid: 0 at x=1, ~0.12 at x=1.5, 0.5 at x=2, ~0.88 at x=2.5
        ramp = 1.0 / (1.0 + np.exp(-4.0 * (x - 2.0)))
        return self.nuclear_scale + (1.0 - self.nuclear_scale) * ramp

    def _effective_scale_prime(self, rho: np.ndarray) -> np.ndarray:
        """Return ``d(_effective_scale)/d rho`` (zero below ρ₀, sigmoid tail above).

        Required for the force derived from the density-dependent scale:
        d/dρ [s(ρ)·e(ρ)] = s·e' + s'·e.  Below ρ₀ the scale is constant
        (ramp clamped by the x≥1 clip), so the derivative vanishes there.
        """

        rho = np.asarray(rho, dtype=float)
        x = rho / self.parameters.rho0
        ramp = 1.0 / (1.0 + np.exp(-4.0 * (np.clip(x, 1.0, None) - 2.0)))
        prime = (1.0 - self.nuclear_scale) * 4.0 * ramp * (1.0 - ramp) / self.parameters.rho0
        return np.where(x > 1.0, prime, 0.0)

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
        sigmas = self._sigmas_for(len(pos))
        grid = self._require_grid(pos)
        rho_n = np.zeros(grid.sizes, dtype=float)
        rho_p = np.zeros(grid.sizes, dtype=float)

        for r_i, is_p, sigma_i in zip(pos, protons, sigmas):
            ix, iy, iz, weights = self._density_contribution(grid, r_i, sigma_i)
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
        """Return analytical grid EDF forces from one density-gradient pass.

        The force on packet i is -∫ (δE/δρ_q(x)) · ∇_{r_i} ρ_i(x) dV with
        ∇_{r_i}ρ_i = (x-r_i)/σ_i² ρ_i.  δE/δρ_q now differentiates the SAME
        functional that ``_potential_from_density`` integrates, including the
        κs surface-symmetry force (previously missing) and the ds/dρ terms of
        the density-dependent calibration scale.
        """

        del epsilon
        pos = self._as_positions(positions)
        protons = np.asarray(is_proton, dtype=bool)
        if len(protons) != len(pos):
            raise ValueError("is_proton length must match positions")

        self.build(pos)
        rho_n, rho_p = self.density_on_grid(pos, protons)
        d_e_drho_n, d_e_drho_p = self._density_functional_derivatives(rho_n, rho_p, use_surface_term)

        d_volume = self.grid_spacing ** 3
        forces = np.zeros_like(pos, dtype=float)
        sigmas = self._sigmas_for(len(pos))
        grid = self._require_grid(pos)

        for i, (r_i, is_p) in enumerate(zip(pos, protons)):
            sigma2 = sigmas[i] ** 2
            ix, iy, iz, rho_i = self._density_contribution(grid, r_i, sigmas[i])
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
        component.  Differentiates exactly the functional in
        ``_potential_from_density``, so it is the reference for checking the
        analytical forces.
        """

        pos = self._as_positions(positions)
        protons = np.asarray(is_proton, dtype=bool)
        if len(protons) != len(pos):
            raise ValueError("is_proton length must match positions")
        if epsilon <= 0.0:
            raise ValueError("epsilon must be positive")
        sigmas = self._sigmas_for(len(pos))
        fixed_grid = self._build_grid(pos, extra_padding=epsilon)
        self._grid = fixed_grid
        rho_n, rho_p = self.density_on_grid(pos, protons)
        old_contributions = [
            self._density_contribution(fixed_grid, r_i, sigmas[i]) for i, r_i in enumerate(pos)
        ]
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
                    sigmas[i],
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
                    sigmas[i],
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
        """Return finite-Gaussian proton-proton direct Coulomb energy.

        For unequal packet widths the pair width entering the Gaussian charge
        overlap is σ_ij = sqrt((σ_i² + σ_j²)/2), reducing to the common-σ
        formula erf(r/(2σ))/r when all widths are equal.
        """

        pos = self._as_positions(positions)
        sigmas = self._sigmas_for(len(pos))
        proton_mask = np.asarray(is_proton, dtype=bool)
        proton_pos = pos[proton_mask]
        proton_sig = sigmas[proton_mask]
        if len(proton_pos) < 2:
            return 0.0
        diff = proton_pos[:, None, :] - proton_pos[None, :, :]
        rij = np.linalg.norm(diff, axis=-1)
        iu = np.triu_indices(len(proton_pos), 1)
        r = np.maximum(rij[iu], 1.0e-12)
        sigma_pair = np.sqrt(0.5 * (proton_sig[iu[0]] ** 2 + proton_sig[iu[1]] ** 2))
        try:
            from scipy.special import erf
        except Exception:  # pragma: no cover
            erf = np.vectorize(__import__("math").erf)
        return float(np.sum(E2 * erf(r / (2.0 * sigma_pair)) / r))

    def _potential_from_density(
        self,
        positions: np.ndarray,
        is_proton: np.ndarray,
        rho_n: np.ndarray,
        rho_p: np.ndarray,
        use_surface_term: bool = True,
    ) -> dict[str, float]:
        rho = rho_n + rho_p
        # The calibration scale is applied LOCALLY as s(ρ(x)) inside the
        # integral (previously a global constant multiplier on the integrated
        # energy, which made energy and forces inconsistent once the forces
        # used the density-dependent form).  With nuclear_scale=1.0 this
        # reduces exactly to the unscaled energy, so the raw calibration
        # baseline is unchanged.
        scale = self._effective_scale(rho)
        skyrme_bulk = self.integrate(scale * self.skyrme_bulk(rho, rho_n, rho_p))
        symmetry = self.integrate(scale * self.symmetry_term(rho, rho_n, rho_p))
        if use_surface_term:
            grad_rho_sq = self.gradients(rho, None)
            surface = self.integrate(scale * self.surface_term(grad_rho_sq))
            surface_sym = self.integrate(scale * self.surface_symmetry(grad_rho_sq, rho, rho_n, rho_p))
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
        """Return ``dE/d rho_n`` and ``dE/d rho_p`` of the integrated functional.

        This is the exact functional derivative (on the discrete grid) of the
        energy in ``_potential_from_density``:

        - bulk + volume-symmetry terms, scaled by the local s(ρ): the force
          gets s(ρ)·e'(ρ) + s'(ρ)·e(ρ) (the s' term was previously absent,
          splitting energy and forces);
        - g0 surface term: (g0/ρ₀)·Σ_ax D_axᵀ(D_ax ρ) — the EXACT transpose
          of the discrete np.gradient stencil (≈ -g0/ρ₀ ∇²ρ in the
          interior); using np.gradient twice instead leaves a boundary
          discretization mismatch between energy and force;
        - κs surface-symmetry term C_s/(2ρ₀)·(-κs)(∇ρ)²δ², whose force was
          previously missing entirely (P0 fix).  With
          δ = (ρn-ρp)/ρ, ∂δ/∂ρn = (1-δ)/ρ and ∂δ/∂ρp = -(1+δ)/ρ, giving
          dE/dρn ⊃ f[2δ(1-δ)(∇ρ)²/ρ + 2Dᵀ(δ²Dρ)] and
          dE/dρp ⊃ f[-2δ(1+δ)(∇ρ)²/ρ + 2Dᵀ(δ²Dρ)], f = -C_s·κs/(2ρ₀),
          with the same exact-transpose Dᵀ (continuum: -2∇·(δ²∇ρ));
        - Slater Coulomb exchange (unscaled, protons only).
        """

        p = self.parameters
        rho_n = np.asarray(rho_n, dtype=float)
        rho_p = np.asarray(rho_p, dtype=float)
        rho = np.clip(rho_n + rho_p, 0.0, None)
        rho_safe = np.maximum(rho, 1.0e-12)
        finite = rho > 1.0e-12
        x = rho / p.rho0
        scale = self._effective_scale(rho)
        scale_prime = self._effective_scale_prime(rho)
        asym = rho_n - rho_p

        bulk_prime = p.alpha * x + p.beta * x**p.gamma + p.g_tau * (p.eta + 1.0) * x**p.eta
        d_n = bulk_prime + p.c_sym / p.rho0 * asym
        d_p = bulk_prime - p.c_sym / p.rho0 * asym

        # Unscaled energy density of the scaled terms; enters the s'(ρ)·e(ρ)
        # part of d[s(ρ)e(ρ)]/dρ.
        e_scaled = (
            p.alpha / (2.0 * p.rho0) * rho**2
            + p.beta / ((p.gamma + 1.0) * p.rho0**p.gamma) * rho ** (p.gamma + 1.0)
            + p.g_tau / (p.rho0**p.eta) * rho ** (p.eta + 1.0)
            + p.c_sym / (2.0 * p.rho0) * asym**2
        )

        if use_surface_term:
            # NOTE: the gradient terms must be differentiated with the exact
            # TRANSPOSE of the discrete np.gradient stencil, not with a second
            # application of np.gradient (the continuum -laplacian).  The two
            # differ at one-sided boundary stencils, and the mismatch showed up
            # as an energy/force split of several MeV/fm per nucleon.
            # δ(Σ w·(Dρ)²)/δρ = 2·Dᵀ(w·Dρ); Dᵀ is built explicitly below.
            dtd = self._dtd(rho)  # Σ_ax D_axᵀ(D_ax ρ) ≈ -laplacian
            d_n = d_n + (p.gsur / p.rho0) * dtd
            d_p = d_p + (p.gsur / p.rho0) * dtd
            grad_rho_sq = self.gradients(rho, None)
            e_scaled = e_scaled + p.gsur / (2.0 * p.rho0) * grad_rho_sq

            delta = np.zeros_like(rho)
            np.divide(asym, rho_safe, out=delta, where=finite)
            e_scaled = e_scaled + p.c_sym / (2.0 * p.rho0) * (-p.kappa_s * grad_rho_sq) * delta * delta

            coeff = -p.c_sym * p.kappa_s / (2.0 * p.rho0)
            # Σ_ax D_axᵀ(δ²·D_ax ρ): exact discrete version of -∇·(δ²∇ρ).
            div_t = self._dtd(rho, weight=delta * delta)
            d_delta_n = np.zeros_like(rho)
            np.divide(2.0 * delta * (1.0 - delta), rho_safe, out=d_delta_n, where=finite)
            d_delta_p = np.zeros_like(rho)
            np.divide(-2.0 * delta * (1.0 + delta), rho_safe, out=d_delta_p, where=finite)
            d_n = d_n + coeff * (d_delta_n * grad_rho_sq + 2.0 * div_t)
            d_p = d_p + coeff * (d_delta_p * grad_rho_sq + 2.0 * div_t)

        out_n = scale * d_n + scale_prime * e_scaled
        out_p = scale * d_p + scale_prime * e_scaled
        # Coulomb exchange stays unscaled, matching _potential_from_density.
        out_p = out_p + (-E2 * (3.0 / np.pi) ** (1.0 / 3.0) * np.maximum(rho_p, 1.0e-12) ** (1.0 / 3.0))
        return out_n, out_p

    def _gradient_components(self, rho: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        edge_order = 2 if min(np.shape(rho)) >= 3 else 1
        return np.gradient(np.asarray(rho, dtype=float), self.grid_spacing, edge_order=edge_order)

    def _grad_matrix(self, n: int) -> np.ndarray:
        """Dense matrix of the 1-D np.gradient stencil (edge_order=2).

        Needed so that functional derivatives of ``∫ w·(Dρ)² dV`` use the exact
        transpose Dᵀ of the discrete gradient, which is what makes analytical
        forces exactly consistent with the gridded energy (np.gradient applied
        twice is NOT the adjoint at the one-sided boundary stencils).
        """

        cached = self._grad_matrix_cache.get(n)
        if cached is not None:
            return cached
        h = self.grid_spacing
        m = np.zeros((n, n), dtype=float)
        if n >= 3:
            idx = np.arange(1, n - 1)
            m[idx, idx + 1] = 0.5 / h
            m[idx, idx - 1] = -0.5 / h
            m[0, 0] = -1.5 / h
            m[0, 1] = 2.0 / h
            m[0, 2] = -0.5 / h
            m[n - 1, n - 1] = 1.5 / h
            m[n - 1, n - 2] = -2.0 / h
            m[n - 1, n - 3] = 0.5 / h
        elif n == 2:
            m[0, 0] = -1.0 / h
            m[0, 1] = 1.0 / h
            m[1, 0] = -1.0 / h
            m[1, 1] = 1.0 / h
        self._grad_matrix_cache[n] = m
        return m

    def _apply_grad_transpose(self, values: np.ndarray, axis: int) -> np.ndarray:
        """Apply the transpose of the 1-D gradient stencil along ``axis``."""

        m = self._grad_matrix(values.shape[axis])
        moved = np.moveaxis(np.asarray(values, dtype=float), axis, 0)
        out = np.tensordot(m.T, moved, axes=(1, 0))
        return np.moveaxis(out, 0, axis)

    def _dtd(self, rho: np.ndarray, weight: np.ndarray | None = None) -> np.ndarray:
        """Return ``Σ_axis D_axisᵀ(weight · D_axis rho)`` (exact discrete adjoint).

        Reduces to minus the Laplacian in the grid interior; used for the
        functional derivatives of the (∇ρ)² surface terms.
        """

        gx, gy, gz = self._gradient_components(rho)
        if weight is not None:
            gx, gy, gz = weight * gx, weight * gy, weight * gz
        return (
            self._apply_grad_transpose(gx, 0)
            + self._apply_grad_transpose(gy, 1)
            + self._apply_grad_transpose(gz, 2)
        )

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
        """Return unscaled finite-Gaussian direct Coulomb forces.

        Uses the same per-pair width σ_ij = sqrt((σ_i²+σ_j²)/2) as
        ``coulomb_direct`` so the two stay self-consistent.
        """

        forces = np.zeros_like(positions, dtype=float)
        sigmas = self._sigmas_for(len(positions))
        proton_indices = np.flatnonzero(is_proton)
        if len(proton_indices) < 2:
            return forces
        proton_pos = positions[proton_indices]
        proton_sig = sigmas[proton_indices]
        i_local, j_local = np.triu_indices(len(proton_pos), 1)
        diff = proton_pos[i_local] - proton_pos[j_local]
        r = np.linalg.norm(diff, axis=1)
        r_safe = np.maximum(r, 1.0e-12)
        try:
            from scipy.special import erf
        except Exception:  # pragma: no cover
            erf = np.vectorize(__import__("math").erf)
        sigma_pair = np.sqrt(0.5 * (proton_sig[i_local] ** 2 + proton_sig[j_local] ** 2))
        x = r_safe / (2.0 * sigma_pair)
        d_erf_over_r = np.exp(-x * x) / (np.sqrt(np.pi) * sigma_pair * r_safe) - erf(x) / (r_safe * r_safe)
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
        sigma_i: float,
        use_surface_term: bool = True,
    ) -> float:
        shifted_positions = positions.copy()
        shifted_positions[index, axis] += shift
        new_ix, new_iy, new_iz, new_weights = self._density_contribution(
            self._require_grid(positions), shifted_positions[index], sigma_i
        )
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
        # Pad with the largest packet width so truncated per-packet Gaussians
        # of every width stay inside the grid.
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
        sigma_r: float | None = None,
    ) -> tuple[slice, slice, slice, np.ndarray]:
        sigma = self.sigma_r if sigma_r is None else float(sigma_r)
        sigma2 = sigma * sigma
        norm = 1.0 / ((2.0 * np.pi * sigma2) ** 1.5)
        cutoff = self.n_sigma * sigma
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
