"""Time propagation for Phase 1 static ImQMD nuclei."""

from __future__ import annotations

import numpy as np

from .collisions import attempt_nn_collision, fermi_constraint_check
from .grid_edf import GridEDF
from .nucleus import ImQMDNucleus
from .skyrme import E2, M_N


def _coulomb_forces(positions: np.ndarray, is_proton: np.ndarray, sigma_r: float) -> np.ndarray:
    forces = np.zeros_like(positions, dtype=float)
    proton_indices = np.where(is_proton)[0]
    if len(proton_indices) < 2:
        return forces
    try:
        from scipy.special import erf
    except Exception:  # pragma: no cover
        erf = np.vectorize(__import__("math").erf)

    local_i, local_j = np.triu_indices(len(proton_indices), k=1)
    i_idx = proton_indices[local_i]
    j_idx = proton_indices[local_j]
    rij = positions[i_idx] - positions[j_idx]
    r = np.linalg.norm(rij, axis=1)
    r_safe = np.maximum(r, 1.0e-8)
    x = r_safe / (2.0 * sigma_r)
    d_erf_over_r = (np.exp(-x * x) / (np.sqrt(np.pi) * sigma_r * r_safe)) - (erf(x) / (r_safe * r_safe))
    # F_i = -dU/dr_i for U=e2*erf(r/2sigma)/r.
    fij = -E2 * d_erf_over_r[:, None] * rij / r_safe[:, None]
    np.add.at(forces, i_idx, fij)
    np.add.at(forces, j_idx, -fij)
    return forces


def _coulomb_exchange_force(nucleus: ImQMDNucleus, positions: np.ndarray) -> np.ndarray:
    """Force matching ``SkyrmeEDF.coulomb_exchange_energy``."""

    forces = np.zeros_like(positions, dtype=float)
    proton_indices = np.where(nucleus.is_proton)[0]
    if len(proton_indices) == 0:
        return forces
    proton_pos = positions[proton_indices]
    sigma2 = nucleus.sigma_r**2
    diff = proton_pos[:, None, :] - proton_pos[None, :, :]
    weights = np.exp(-np.sum(diff * diff, axis=-1) / (2.0 * sigma2))
    norm = 1.0 / ((2.0 * np.pi * sigma2) ** 1.5)
    rho = np.clip(norm * np.sum(weights, axis=1), 1.0e-12, None)
    coeff = -0.75 * E2 * (3.0 / np.pi) ** (1.0 / 3.0) / nucleus.edf.parameters.rho0
    d_e_d_rho = coeff * (4.0 / 3.0) * rho ** (1.0 / 3.0)

    grad_as_density_source = np.sum(
        (d_e_d_rho[:, None] * norm * weights)[:, :, None] * diff / sigma2,
        axis=0,
    )
    grad_as_sample_point = np.sum(
        (d_e_d_rho[:, None] * norm * weights)[:, :, None] * (-diff) / sigma2,
        axis=1,
    )
    forces[proton_indices] = -(grad_as_density_source + grad_as_sample_point)
    return forces


def _centroid_density_force(nucleus: ImQMDNucleus, positions: np.ndarray) -> np.ndarray:
    """Analytic force matching the centroid-sampled EDF in ``energy_components``."""

    p = nucleus.edf.parameters
    sigma2 = nucleus.sigma_r**2
    diff = positions[:, None, :] - positions[None, :, :]
    g = np.exp(-np.sum(diff * diff, axis=-1) / (2.0 * sigma2))
    norm = 1.0 / ((2.0 * np.pi * sigma2) ** 1.5)
    weights = norm * g
    is_proton = nucleus.is_proton
    rho_p = np.sum(weights[:, is_proton], axis=1)
    rho_n = np.sum(weights[:, ~is_proton], axis=1)
    rho = np.clip(rho_n + rho_p, 1.0e-12, None)

    bulk_prime = (
        p.alpha / p.rho0 * rho
        + p.beta / (p.rho0**p.gamma) * rho**p.gamma
        + p.g_tau * (p.eta + 1.0) / (p.rho0**p.eta) * rho**p.eta
    )
    delta_rho = rho_n - rho_p
    c_sym = p.c_sym / p.rho0
    sym_n = c_sym * delta_rho
    sym_p = -c_sym * delta_rho
    coeff_n = bulk_prime + sym_n
    coeff_p = bulk_prime + sym_p

    coeff_by_column_type = np.where(is_proton[None, :], coeff_p[:, None], coeff_n[:, None])
    grad_as_density_source = np.sum(
        (coeff_by_column_type * weights)[:, :, None] * diff / sigma2,
        axis=0,
    )

    coeff_by_sampled_species = np.where(is_proton[None, :], coeff_p[:, None], coeff_n[:, None])
    grad_as_sample_point = np.sum(
        (coeff_by_sampled_species * weights)[:, :, None] * (-diff) / sigma2,
        axis=1,
    )

    return -(grad_as_density_source + grad_as_sample_point)


def _surface_pair_force(nucleus: ImQMDNucleus, positions: np.ndarray) -> np.ndarray:
    """Force matching ``SkyrmeEDF.surface_pair_energy``."""

    if nucleus.A < 2:
        return np.zeros_like(positions, dtype=float)
    sigma2 = nucleus.sigma_r**2
    diff = positions[:, None, :] - positions[None, :, :]
    weights = np.exp(-np.sum(diff * diff, axis=-1) / (4.0 * sigma2))
    strength = nucleus.edf.parameters.gsur / 7.0
    return -strength * np.sum(weights[:, :, None] * diff, axis=1) / (2.0 * sigma2)


def _surface_symmetry_pair_force(nucleus: ImQMDNucleus, positions: np.ndarray) -> np.ndarray:
    """Force matching ``SkyrmeEDF.surface_symmetry_energy``."""

    if nucleus.A < 2:
        return np.zeros_like(positions, dtype=float)
    sigma2 = nucleus.sigma_r**2
    diff = positions[:, None, :] - positions[None, :, :]
    weights = np.exp(-np.sum(diff * diff, axis=-1) / (4.0 * sigma2))
    tau = np.where(nucleus.is_proton, -1.0, 1.0)
    tau_pair = tau[:, None] * tau[None, :]
    strength = nucleus.edf.parameters.c_sym * nucleus.edf.parameters.kappa_s / 7.0
    return -strength * np.sum((tau_pair * weights)[:, :, None] * diff, axis=1) / (2.0 * sigma2)


def _static_reference_force(nucleus: ImQMDNucleus, positions: np.ndarray) -> np.ndarray:
    reference_positions = nucleus.reference_positions
    group_ids = getattr(nucleus, "reference_group_ids", None)
    static_k = getattr(nucleus.edf, "static_k", 0.0)
    if reference_positions is None or static_k <= 0.0:
        return np.zeros_like(positions, dtype=float)
    if group_ids is None:
        return nucleus.edf.static_mean_field_force(positions, reference_positions)

    groups = np.asarray(group_ids)
    if groups.shape[0] != positions.shape[0]:
        return nucleus.edf.static_mean_field_force(positions, reference_positions)

    forces = np.zeros_like(positions, dtype=float)
    reference_positions = np.asarray(reference_positions, dtype=float)
    for group in np.unique(groups):
        mask = groups == group
        reference_cm = np.mean(reference_positions[mask], axis=0)
        position_cm = np.mean(positions[mask], axis=0)
        shifted_reference = reference_positions[mask] + (position_cm - reference_cm)
        forces[mask] = -static_k * (positions[mask] - shifted_reference)
    return forces


def _derivatives(
    nucleus: ImQMDNucleus,
    positions: np.ndarray,
    momenta: np.ndarray,
    use_surface_term: bool = True,
    use_static_stabilizer: bool = False,
    grid_edf: GridEDF | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    drdt = momenta / M_N
    dpdt = _static_reference_force(nucleus, positions) if use_static_stabilizer else np.zeros_like(positions, dtype=float)
    if grid_edf is None:
        dpdt += _centroid_density_force(nucleus, positions)
        if use_surface_term:
            dpdt += _surface_pair_force(nucleus, positions)
            dpdt += _surface_symmetry_pair_force(nucleus, positions)
        dpdt += _coulomb_forces(positions, nucleus.is_proton, nucleus.sigma_r)
        dpdt += _coulomb_exchange_force(nucleus, positions)
    else:
        dpdt += grid_edf.forces_analytical(
            positions,
            nucleus.is_proton,
            use_surface_term=use_surface_term,
        )
    return drdt, dpdt


def _rk4_step(
    nucleus: ImQMDNucleus,
    dt: float,
    remove_cm_drift: bool = True,
    use_surface_term: bool = True,
    use_static_stabilizer: bool = False,
    grid_edf: GridEDF | None = None,
) -> None:
    r0 = nucleus.positions
    p0 = nucleus.momenta

    k1r, k1p = _derivatives(
        nucleus,
        r0,
        p0,
        use_surface_term=use_surface_term,
        use_static_stabilizer=use_static_stabilizer,
        grid_edf=grid_edf,
    )
    k2r, k2p = _derivatives(
        nucleus,
        r0 + 0.5 * dt * k1r,
        p0 + 0.5 * dt * k1p,
        use_surface_term=use_surface_term,
        use_static_stabilizer=use_static_stabilizer,
        grid_edf=grid_edf,
    )
    k3r, k3p = _derivatives(
        nucleus,
        r0 + 0.5 * dt * k2r,
        p0 + 0.5 * dt * k2p,
        use_surface_term=use_surface_term,
        use_static_stabilizer=use_static_stabilizer,
        grid_edf=grid_edf,
    )
    k4r, k4p = _derivatives(
        nucleus,
        r0 + dt * k3r,
        p0 + dt * k3p,
        use_surface_term=use_surface_term,
        use_static_stabilizer=use_static_stabilizer,
        grid_edf=grid_edf,
    )

    nucleus.positions = r0 + (dt / 6.0) * (k1r + 2.0 * k2r + 2.0 * k3r + k4r)
    nucleus.momenta = p0 + (dt / 6.0) * (k1p + 2.0 * k2p + 2.0 * k3p + k4p)

    if remove_cm_drift:
        # Remove tiny center-of-mass drift for single-nucleus stability runs.
        positions = nucleus.positions
        momenta = nucleus.momenta
        nucleus.positions = positions - np.mean(positions, axis=0)
        nucleus.momenta = momenta - np.mean(momenta, axis=0)


def _snapshot(nucleus: ImQMDNucleus, time: float, grid_edf: GridEDF | None = None) -> dict[str, float]:
    if grid_edf is None:
        grid_edf = GridEDF(
            nucleus.edf.parameters,
            nucleus.sigma_r,
            nuclear_scale=float(getattr(nucleus, "grid_nuclear_scale", 1.0)),
        )
    components = grid_edf.total_energy(
        nucleus.positions,
        nucleus.momenta,
        nucleus.is_proton,
        use_surface_term=bool(getattr(nucleus, "use_surface_term", True)),
    )
    if bool(getattr(nucleus, "use_static_stabilizer", False)):
        static = nucleus.edf.static_mean_field_energy(nucleus.positions, nucleus.reference_positions)
    else:
        static = 0.0
    components["static"] = static
    components["potential"] += static + nucleus.energy_offset
    components["total"] += static + nucleus.energy_offset
    rp, rn = nucleus.rms_radius()
    return {
        "time": float(time),
        "total": components["total"],
        "kinetic": components["kinetic"],
        "potential": components["potential"],
        "rms_p": rp,
        "rms_n": rn,
        "rms": float(np.sqrt((nucleus.Z * rp * rp + nucleus.N * rn * rn) / nucleus.A)),
        "max_radius": nucleus.max_radius(),
    }


def propagate(
    nucleus: ImQMDNucleus,
    dt: float,
    n_steps: int,
    sample_every: int | None = None,
    with_collisions: bool = False,
    collision_dt: float | None = None,
    remove_cm_drift: bool = True,
    use_surface_term: bool = True,
    use_static_stabilizer: bool = False,
    use_grid_edf: bool = False,
    apply_fermi_constraint: bool | None = None,
) -> list[dict[str, float]]:
    """Propagate centroids with fourth-order Runge-Kutta.

    Time is in fm/c.  Momenta are MeV/c, so ``dot(r)=p/m`` has units of c.
    """

    if sample_every is None:
        sample_every = max(1, n_steps // 200)
    collision_interval = 1
    if apply_fermi_constraint is None:
        apply_fermi_constraint = True  # always apply, but relaxed during collisions
    fermi_threshold = 170.0 if with_collisions else 255.0
    fermi_time = 20.0 if use_static_stabilizer else 5.0
    fermi_interval = max(1, int(round(fermi_time / float(dt))))
    if collision_dt is not None:
        collision_interval = max(1, int(round(float(collision_dt) / float(dt))))
    nucleus.use_surface_term = bool(use_surface_term)
    nucleus.use_static_stabilizer = bool(use_static_stabilizer)
    if use_grid_edf:
        grid_edf = GridEDF(
            nucleus.edf.parameters,
            nucleus.sigma_r,
            nuclear_scale=float(getattr(nucleus, "grid_nuclear_scale", 1.0)),
        )
    else:
        grid_edf = None
    history = [_snapshot(nucleus, 0.0, grid_edf=grid_edf)]
    for step in range(1, int(n_steps) + 1):
        _rk4_step(
            nucleus,
            float(dt),
            remove_cm_drift=remove_cm_drift,
            use_surface_term=use_surface_term,
            use_static_stabilizer=use_static_stabilizer,
            grid_edf=grid_edf,
        )
        if with_collisions and step % collision_interval == 0:
            attempt_nn_collision(nucleus, dt=float(dt) * collision_interval)
        if with_collisions and apply_fermi_constraint and step % fermi_interval == 0:
            fermi_constraint_check(nucleus, threshold=fermi_threshold)
        if step % sample_every == 0 or step == n_steps:
            history.append(_snapshot(nucleus, step * dt, grid_edf=grid_edf))
    return history


__all__ = ["propagate"]
