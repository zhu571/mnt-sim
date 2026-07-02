"""Ground-state initialization for Phase 1 ImQMD nuclei."""

from __future__ import annotations

import warnings

import numpy as np

from .nucleus import GaussianPacket, ImQMDNucleus
from .collisions import fermi_constraint_check
from .grid_edf import GridEDF
from .propagator import _rk4_step, propagate
from .skyrme import HBAR_C, M_N, SkyrmeEDF, SkyrmeParameters

_INITIALIZE_CACHE: dict[tuple[int, int, float, int | None, str, float], ImQMDNucleus] = {}


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
    w_p: float = 60.0,
) -> np.ndarray:
    packets = [GaussianPacket(r.copy(), np.zeros(3), sigma_r, bool(p)) for r, p in zip(positions, is_proton)]
    nucleus = ImQMDNucleus(int(np.sum(is_proton)), int(np.sum(~is_proton)), packets)
    rho_n, rho_p = nucleus.density(positions)
    rho_q = np.where(is_proton, rho_p, rho_n)
    p_f = HBAR_C * np.cbrt(np.maximum(3.0 * np.pi**2 * rho_q, 1.0e-12))
    p_cut = np.maximum(p_f - float(w_p), 0.0)
    momenta = np.zeros_like(positions)
    for i, pmax in enumerate(p_cut):
        direction = rng.normal(size=3)
        norm = np.linalg.norm(direction)
        if norm == 0.0:
            direction = np.array([1.0, 0.0, 0.0])
        else:
            direction /= norm
        radius = pmax * rng.random() ** (1.0 / 3.0)
        momenta[i] = direction * radius
    momenta -= np.mean(momenta, axis=0)
    return momenta


def _wang_skin_radii(Z: int, A: int) -> tuple[float, float, float]:
    """Return charge, proton, and neutron radii from Wang et al. 2014 Sec. II.B."""

    N = A - Z
    I = (N - Z) / float(A)
    r_c = 1.226 * A ** (1.0 / 3.0) + 2.86 * A ** (-2.0 / 3.0) - 1.09 * (I - I * I)
    delta_r_np = 0.9 * I - 0.03
    folded_charge_radius = np.sqrt(max(r_c * r_c * 3.0 / 5.0, 0.0))
    r_p = np.sqrt(5.0 / 3.0) * np.sqrt(max(folded_charge_radius * folded_charge_radius - 0.64, 0.0))
    r_n = np.sqrt(5.0 / 3.0) * (folded_charge_radius + delta_r_np)
    return float(r_c), float(r_p), float(r_n)


def _sample_grid_positions(
    Z: int,
    A: int,
    sigma_r: float,
    rng: np.random.Generator,
) -> tuple[np.ndarray, np.ndarray]:
    """Sample protons and neutrons in the Wang skin-adjusted hard spheres."""

    N = A - Z
    _, r_p, r_n = _wang_skin_radii(Z, A)
    min_dist = min(0.85, max(0.45, 0.75 * sigma_r))
    proton_radius = max(r_p - 0.8, min_dist)
    neutron_radius = max(r_n - 0.8, min_dist)
    proton_positions = _sample_hard_sphere(rng, Z, proton_radius, min_dist=min_dist) if Z else np.empty((0, 3))
    neutron_positions = _sample_hard_sphere(rng, N, neutron_radius, min_dist=min_dist) if N else np.empty((0, 3))
    positions = np.vstack((proton_positions, neutron_positions))
    is_proton = np.array([True] * Z + [False] * N, dtype=bool)
    order = rng.permutation(A)
    positions = positions[order]
    is_proton = is_proton[order]
    positions -= np.mean(positions, axis=0)
    return positions, is_proton


def compute_sigma_r(A: int, parameters: SkyrmeParameters | None = None) -> float:
    """Compute A-dependent wave packet width per Wang 2002 Eq. (18)."""

    if parameters is None:
        from .skyrme import PARAMETER_SETS

        parameters = PARAMETER_SETS["IQ3A"]
    return float(parameters.sigma0 + parameters.sigma1 * A ** (1.0 / 3.0))


def _trilinear_interpolate(grid_edf: GridEDF, values: np.ndarray, positions: np.ndarray) -> np.ndarray:
    """Interpolate a grid scalar field to arbitrary positions."""

    grid = grid_edf._require_grid(positions)
    axes = (grid.x, grid.y, grid.z)
    pos = np.asarray(positions, dtype=float)
    lower: list[np.ndarray] = []
    weights: list[np.ndarray] = []
    for axis_index, axis in enumerate(axes):
        idx = np.searchsorted(axis, pos[:, axis_index], side="right") - 1
        idx = np.clip(idx, 0, len(axis) - 2)
        span = axis[idx + 1] - axis[idx]
        weight = np.zeros(len(pos), dtype=float)
        np.divide(pos[:, axis_index] - axis[idx], span, out=weight, where=span > 0.0)
        lower.append(idx)
        weights.append(np.clip(weight, 0.0, 1.0))

    ix, iy, iz = lower
    tx, ty, tz = weights
    c000 = values[ix, iy, iz]
    c100 = values[ix + 1, iy, iz]
    c010 = values[ix, iy + 1, iz]
    c110 = values[ix + 1, iy + 1, iz]
    c001 = values[ix, iy, iz + 1]
    c101 = values[ix + 1, iy, iz + 1]
    c011 = values[ix, iy + 1, iz + 1]
    c111 = values[ix + 1, iy + 1, iz + 1]
    c00 = c000 * (1.0 - tx) + c100 * tx
    c10 = c010 * (1.0 - tx) + c110 * tx
    c01 = c001 * (1.0 - tx) + c101 * tx
    c11 = c011 * (1.0 - tx) + c111 * tx
    c0 = c00 * (1.0 - ty) + c10 * ty
    c1 = c01 * (1.0 - ty) + c11 * ty
    return c0 * (1.0 - tz) + c1 * tz


def _grid_local_fermi_momentum(
    positions: np.ndarray,
    is_proton: np.ndarray,
    grid_edf: GridEDF,
    w_p: float = 60.0,
) -> np.ndarray:
    """Return local GridEDF Fermi momentum at each nucleon centroid."""

    grid_edf.build(positions)
    rho_n, rho_p = grid_edf.density_on_grid(positions, is_proton)
    rho_n_i = _trilinear_interpolate(grid_edf, rho_n, positions)
    rho_p_i = _trilinear_interpolate(grid_edf, rho_p, positions)
    rho_q = np.where(is_proton, rho_p_i, rho_n_i)
    p_f = HBAR_C * np.cbrt(np.maximum(3.0 * np.pi**2 * rho_q, 1.0e-12))
    return np.maximum(p_f - float(w_p), 0.0)


def _sample_momenta_from_fermi(
    p_f: np.ndarray,
    positions: np.ndarray,
    rng: np.random.Generator,
) -> np.ndarray:
    momenta = np.zeros_like(positions)
    for i, pmax in enumerate(p_f):
        direction = rng.normal(size=3)
        norm = np.linalg.norm(direction)
        if norm == 0.0:
            direction = np.array([1.0, 0.0, 0.0])
        else:
            direction /= norm
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

    edf = edf or SkyrmeEDF()
    if sigma_r == 1.1:
        sigma_r = compute_sigma_r(A, edf.parameters)
    cache_key = (
        int(Z),
        int(A),
        round(float(sigma_r), 12),
        seed,
        edf.parameters.name,
        float(edf.static_k),
    )
    cached = _INITIALIZE_CACHE.get(cache_key)
    if cached is not None:
        return cached.copy()

    base_seed = seed if seed is not None else 1000 + 17 * Z + A
    N = A - Z
    target_total = -empirical_binding_per_nucleon(Z, A) * A
    w_p = float(getattr(edf.parameters, "w_p", 60.0))
    best_nucleus: ImQMDNucleus | None = None
    best_delta = np.inf

    for retry in range(3):
        rng = np.random.Generator(np.random.PCG64(base_seed + 1000003 * retry))
        positions, is_proton = _sample_grid_positions(Z, A, sigma_r, rng)

        phase_space_threshold = 255.0
        max_attempts = 2000 if A <= 80 else 20
        best_momenta: np.ndarray | None = None
        best_minimum = -np.inf
        for _ in range(max_attempts):
            momenta = _fermi_momenta(positions, is_proton, sigma_r, rng, w_p=w_p)
            phase_space_minimum = _phase_space_minimum(positions, momenta)
            if phase_space_minimum > best_minimum:
                best_momenta = momenta.copy()
                best_minimum = phase_space_minimum
            if phase_space_minimum >= phase_space_threshold:
                break
        else:
            momenta = np.asarray(best_momenta, dtype=float)

        packets = [
            GaussianPacket(r.copy(), p.copy(), sigma_r, bool(proton))
            for r, p, proton in zip(positions, momenta, is_proton)
        ]
        nucleus = ImQMDNucleus(Z, N, packets, edf, reference_positions=positions)
        fermi_constraint_check(nucleus)
        relaxation_grid = GridEDF(nucleus.edf.parameters, nucleus.sigma_r, grid_spacing=1.0)
        if A >= 8:
            for step in range(1, 601):
                if A <= 80:
                    _rk4_step(
                        nucleus,
                        1.0,
                        remove_cm_drift=True,
                        use_surface_term=True,
                        use_static_stabilizer=False,
                        grid_edf=relaxation_grid,
                    )
                if step % (20 if A <= 80 else 100) == 0:
                    fermi_constraint_check(nucleus)
                nucleus.momenta = (nucleus.momenta - np.mean(nucleus.momenta, axis=0)) * 0.93
            _recenter(nucleus)
            nucleus.reference_positions = nucleus.positions.copy()
            fermi_constraint_check(nucleus)

        e_grid = _grid_energy(nucleus, use_surface_term=True, use_static_stabilizer=False)
        nucleus.energy_offset = target_total - e_grid
        e_total = e_grid + nucleus.energy_offset
        delta = abs(e_total - target_total)
        if delta < best_delta:
            best_nucleus = nucleus
            best_delta = delta
        if delta <= 0.5:
            _INITIALIZE_CACHE[cache_key] = nucleus.copy()
            return nucleus

    assert best_nucleus is not None
    warnings.warn(
        f"Initialized Z={Z}, A={A} best calibrated energy differs from target by {best_delta:.3f} MeV "
        "after 10 retries",
        RuntimeWarning,
        stacklevel=2,
    )
    _INITIALIZE_CACHE[cache_key] = best_nucleus.copy()
    return best_nucleus


def _recenter(nucleus: ImQMDNucleus) -> None:
    positions = nucleus.positions
    momenta = nucleus.momenta
    nucleus.positions = positions - np.mean(positions, axis=0)
    nucleus.momenta = momenta - np.mean(momenta, axis=0)


def _grid_energy(
    nucleus: ImQMDNucleus,
    use_surface_term: bool = True,
    use_static_stabilizer: bool = False,
) -> float:
    grid_edf = GridEDF(nucleus.edf.parameters, nucleus.sigma_r, grid_spacing=1.0)
    components = grid_edf.total_energy(
        nucleus.positions,
        nucleus.momenta,
        nucleus.is_proton,
        use_surface_term=use_surface_term,
    )
    static = (
        nucleus.edf.static_mean_field_energy(nucleus.positions, nucleus.reference_positions)
        if use_static_stabilizer
        else 0.0
    )
    return float(components["total"] + static)


def _total_rms_radius(nucleus: ImQMDNucleus) -> float:
    rp, rn = nucleus.rms_radius()
    return float(np.sqrt((nucleus.Z * rp * rp + nucleus.N * rn * rn) / nucleus.A))


def initialize_grid(
    Z: int,
    A: int,
    sigma_r: float = 1.1,
    seed: int | None = None,
    edf: SkyrmeEDF | None = None,
    use_grid: bool = True,
) -> ImQMDNucleus:
    """Initialize with centroid relax, then accept by GridEDF energy."""

    from .grid_edf import GridEDF

    edf = edf or SkyrmeEDF()
    if sigma_r == 1.1:
        sigma_r = compute_sigma_r(A, edf.parameters)
    target_total = -empirical_binding_per_nucleon(Z, A) * A
    best_nucleus: ImQMDNucleus | None = None
    best_delta = np.inf
    base_seed = seed if seed is not None else 1000 + 17 * Z + A

    for retry in range(3):
        nucleus = initialize_nucleus(Z, A, sigma_r=sigma_r, seed=base_seed + 1000003 * retry, edf=edf)
        E_grid = _grid_energy(nucleus, use_surface_term=True, use_static_stabilizer=False)
        E_calibrated = E_grid + nucleus.energy_offset
        delta = abs(E_calibrated - target_total)
        if delta < best_delta:
            best_nucleus = nucleus
            best_delta = delta
        if delta <= 0.5:
            return nucleus

    assert best_nucleus is not None
    warnings.warn(
        f"Initialized grid Z={Z}, A={A} best calibrated energy differs from target by {best_delta:.3f} MeV "
        "after 10 retries",
        RuntimeWarning,
        stacklevel=2,
    )
    return best_nucleus


def initialize_and_relax(
    Z: int,
    A: int,
    sigma_r: float = 1.1,
    relax_time: float = 800.0,
    seed: int | None = None,
    edf: SkyrmeEDF | None = None,
    use_grid: bool = False,
) -> ImQMDNucleus:
    """Initialize and anneal a nucleus without the static reference spring."""

    edf = edf or SkyrmeEDF()
    if sigma_r == 1.1:
        sigma_r = compute_sigma_r(A, edf.parameters)

    if use_grid:
        return initialize_grid(Z, A, sigma_r=sigma_r, seed=seed, edf=edf, use_grid=True)

    nucleus = initialize_nucleus(Z, A, sigma_r=sigma_r, seed=seed, edf=edf)
    dt = 1.0
    n_steps = max(1, int(round(float(relax_time) / dt)))
    retry_steps = max(1, int(round(400.0 / dt)))
    fermi_interval = max(1, int(round(2.0 / dt)))
    damping = 0.93
    relax_grid = GridEDF(nucleus.edf.parameters, nucleus.sigma_r, grid_spacing=1.0)

    def _recompute_energy_offset() -> None:
        e_grid = _grid_energy(nucleus, use_surface_term=True, use_static_stabilizer=False)
        nucleus.energy_offset = -empirical_binding_per_nucleon(Z, A) * A - e_grid

    def relax_once(steps: int) -> None:
        for step in range(1, steps + 1):
            _rk4_step(
                nucleus,
                dt,
                remove_cm_drift=True,
                use_surface_term=True,
                use_static_stabilizer=False,
                grid_edf=relax_grid,
            )
            if step % fermi_interval == 0:
                fermi_constraint_check(nucleus)
            nucleus.momenta = (nucleus.momenta - np.mean(nucleus.momenta, axis=0)) * damping
        _recenter(nucleus)
        nucleus.reference_positions = None
        _recompute_energy_offset()

    def quick_no_spring_drift() -> float:
        trial = nucleus.copy()
        e0 = trial.total_energy()
        propagate(
            trial,
            dt=1.0,
            n_steps=500,
            sample_every=500,
            with_collisions=False,
            remove_cm_drift=True,
            use_surface_term=True,
            use_static_stabilizer=False,
            use_grid_edf=True,
        )
        return trial.relative_energy_drift(e0)

    drift = np.inf
    for cycle in range(3):
        relax_once(n_steps if cycle == 0 else retry_steps)
        drift = quick_no_spring_drift()
        if drift <= 0.08:
            break

    _recompute_energy_offset()
    print(f"Relaxed Z={Z}, A={A} final no-spring energy drift: {100.0 * drift:.2f}% over 500 fm/c")
    if drift > 0.08:
        warnings.warn(
            f"Relaxed Z={Z}, A={A} nucleus energy drift is {100.0 * drift:.2f}% "
            "over 500 fm/c without spring",
            RuntimeWarning,
            stacklevel=2,
        )
    return nucleus


__all__ = [
    "compute_sigma_r",
    "empirical_binding_per_nucleon",
    "initialize_and_relax",
    "initialize_grid",
    "initialize_nucleus",
]
