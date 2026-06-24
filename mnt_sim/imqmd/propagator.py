"""Time propagation for Phase 1 static ImQMD nuclei."""

from __future__ import annotations

import numpy as np

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
        erf = np.vectorize(np.math.erf)

    for local_i, i in enumerate(proton_indices[:-1]):
        js = proton_indices[local_i + 1 :]
        rij = positions[i] - positions[js]
        r = np.linalg.norm(rij, axis=1)
        r_safe = np.maximum(r, 1.0e-8)
        x = r_safe / (2.0 * sigma_r)
        d_erf_over_r = (np.exp(-x * x) / (np.sqrt(np.pi) * sigma_r * r_safe)) - (erf(x) / (r_safe * r_safe))
        # F_i = -dU/dr_i for U=e2*erf(r/2sigma)/r.
        fij = -E2 * d_erf_over_r[:, None] * rij / r_safe[:, None]
        forces[i] += np.sum(fij, axis=0)
        forces[js] -= fij
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
        p.alpha / (2.0 * p.rho0)
        + p.beta * p.gamma / ((p.gamma + 1.0) * p.rho0**p.gamma) * rho ** (p.gamma - 1.0)
        + p.g_tau * p.eta / (p.rho0**p.eta) * rho ** (p.eta - 1.0)
    )
    delta_rho = rho_n - rho_p
    c_sym = p.c_sym / (2.0 * p.rho0)
    sym_n = c_sym * (2.0 * delta_rho * rho - delta_rho * delta_rho) / (rho * rho)
    sym_p = c_sym * (-2.0 * delta_rho * rho - delta_rho * delta_rho) / (rho * rho)
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


def _derivatives(nucleus: ImQMDNucleus, positions: np.ndarray, momenta: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    drdt = momenta / M_N
    dpdt = nucleus.edf.static_mean_field_force(positions, nucleus.reference_positions)
    dpdt += _centroid_density_force(nucleus, positions)
    dpdt += _coulomb_forces(positions, nucleus.is_proton, nucleus.sigma_r)
    return drdt, dpdt


def _rk4_step(nucleus: ImQMDNucleus, dt: float) -> None:
    r0 = nucleus.positions
    p0 = nucleus.momenta

    k1r, k1p = _derivatives(nucleus, r0, p0)
    k2r, k2p = _derivatives(nucleus, r0 + 0.5 * dt * k1r, p0 + 0.5 * dt * k1p)
    k3r, k3p = _derivatives(nucleus, r0 + 0.5 * dt * k2r, p0 + 0.5 * dt * k2p)
    k4r, k4p = _derivatives(nucleus, r0 + dt * k3r, p0 + dt * k3p)

    nucleus.positions = r0 + (dt / 6.0) * (k1r + 2.0 * k2r + 2.0 * k3r + k4r)
    nucleus.momenta = p0 + (dt / 6.0) * (k1p + 2.0 * k2p + 2.0 * k3p + k4p)

    # Remove tiny center-of-mass drift.
    positions = nucleus.positions
    momenta = nucleus.momenta
    nucleus.positions = positions - np.mean(positions, axis=0)
    nucleus.momenta = momenta - np.mean(momenta, axis=0)


def _snapshot(nucleus: ImQMDNucleus, time: float) -> dict[str, float]:
    components = nucleus.energy_components()
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


def propagate(nucleus: ImQMDNucleus, dt: float, n_steps: int, sample_every: int | None = None) -> list[dict[str, float]]:
    """Propagate centroids with fourth-order Runge-Kutta.

    Time is in fm/c.  Momenta are MeV/c, so ``dot(r)=p/m`` has units of c.
    """

    if sample_every is None:
        sample_every = max(1, n_steps // 200)
    history = [_snapshot(nucleus, 0.0)]
    for step in range(1, int(n_steps) + 1):
        _rk4_step(nucleus, float(dt))
        if step % sample_every == 0 or step == n_steps:
            history.append(_snapshot(nucleus, step * dt))
    return history


__all__ = ["propagate"]
