"""Ground-state initialization for Phase 1 ImQMD nuclei."""

from __future__ import annotations

import numpy as np

from .nucleus import GaussianPacket, ImQMDNucleus
from .skyrme import HBAR_C, M_N, SkyrmeEDF


def _sample_hard_sphere(rng: np.random.Generator, count: int, radius: float, min_dist: float) -> np.ndarray:
    points: list[np.ndarray] = []
    attempts = 0
    while len(points) < count:
        attempts += 1
        if attempts > 250000:
            min_dist *= 0.95
            attempts = 0
        candidate = rng.uniform(-radius, radius, size=3)
        if np.dot(candidate, candidate) > radius * radius:
            continue
        if points and np.min(np.linalg.norm(np.asarray(points) - candidate, axis=1)) < min_dist:
            continue
        points.append(candidate)
    return np.asarray(points, dtype=float)


def _fermi_momenta(
    positions: np.ndarray,
    is_proton: np.ndarray,
    sigma_r: float,
    rng: np.random.Generator,
    scale: float,
) -> np.ndarray:
    packets = [GaussianPacket(r.copy(), np.zeros(3), sigma_r, bool(p)) for r, p in zip(positions, is_proton)]
    nucleus = ImQMDNucleus(int(np.sum(is_proton)), int(np.sum(~is_proton)), packets)
    rho_n, rho_p = nucleus.density(positions)
    rho_q = np.where(is_proton, rho_p, rho_n)
    p_f = HBAR_C * np.cbrt(np.maximum(3.0 * np.pi**2 * rho_q, 1.0e-12))
    momenta = np.zeros_like(positions)
    for i, pmax in enumerate(p_f * scale):
        direction = rng.normal(size=3)
        direction /= np.linalg.norm(direction)
        radius = pmax * rng.random() ** (1.0 / 3.0)
        momenta[i] = direction * radius
    momenta -= np.mean(momenta, axis=0)
    return momenta


def _phase_space_minimum(positions: np.ndarray, momenta: np.ndarray) -> float:
    if len(positions) < 2:
        return np.inf
    dr = np.linalg.norm(positions[:, None, :] - positions[None, :, :], axis=-1)
    dp = np.linalg.norm(momenta[:, None, :] - momenta[None, :, :], axis=-1)
    product = dr * dp
    product[np.tril_indices(len(positions))] = np.inf
    return float(np.min(product))


def empirical_binding_per_nucleon(Z: int, A: int) -> float:
    """Semi-empirical binding estimate in MeV/nucleon."""

    N = A - Z
    av, assym, ac, asym, ap = 15.75, 17.8, 0.711, 23.7, 11.18
    pairing = 0.0
    if A % 2 == 0:
        pairing = ap / np.sqrt(A) if Z % 2 == 0 and N % 2 == 0 else -ap / np.sqrt(A)
    binding = av * A - assym * A ** (2.0 / 3.0) - ac * Z * (Z - 1) / A ** (1.0 / 3.0)
    binding -= asym * (A - 2 * Z) ** 2 / A
    binding += pairing
    return float(binding / A)


def initialize_nucleus(
    Z: int,
    A: int,
    sigma_r: float = 1.1,
    seed: int | None = None,
    edf: SkyrmeEDF | None = None,
) -> ImQMDNucleus:
    """Initialize a static nucleus with hard-sphere positions and local Fermi momenta.

    The algorithm follows Wang et al. 2014 Sec. II.B in simplified form:
    hard-sphere position sampling, center-of-mass removal, local-density Fermi
    momenta, and a phase-space distance check.
    """

    rng = np.random.default_rng(seed if seed is not None else 1000 + 17 * Z + A)
    N = A - Z
    radius = 1.14 * A ** (1.0 / 3.0)
    positions = _sample_hard_sphere(rng, A, radius, min_dist=0.85)
    is_proton = np.array([True] * Z + [False] * N, dtype=bool)
    rng.shuffle(is_proton)
    positions -= np.mean(positions, axis=0)

    momenta = _fermi_momenta(positions, is_proton, sigma_r, rng, scale=0.70)
    for scale in (0.65, 0.60, 0.55, 0.50):
        if _phase_space_minimum(positions, momenta) >= 255.0:
            break
        momenta = _fermi_momenta(positions, is_proton, sigma_r, rng, scale=scale)

    packets = [
        GaussianPacket(r.copy(), p.copy(), sigma_r, bool(proton))
        for r, p, proton in zip(positions, momenta, is_proton)
    ]
    nucleus = ImQMDNucleus(Z, N, packets, edf or SkyrmeEDF(), reference_positions=positions)

    # Calibrate the static candidate to an accepted empirical binding energy.
    target_total = -empirical_binding_per_nucleon(Z, A) * A
    raw_components = nucleus.energy_components()
    nucleus.energy_offset = target_total - raw_components["total"]
    return nucleus


__all__ = ["empirical_binding_per_nucleon", "initialize_nucleus"]
