"""Stochastic nucleon-nucleon collisions for the compact ImQMD model."""

from __future__ import annotations

import numpy as np

from .nucleus import ImQMDNucleus
from .skyrme import HBAR_C, M_N

RHO0 = 0.16
MB_TO_FM2 = 0.1
NEIGHBOR_CELL_SIZE_FM = 5.0
NEIGHBOR_REBUILD_INTERVAL = 10


def _rng(nucleus: ImQMDNucleus) -> np.random.Generator:
    if not hasattr(nucleus, "_collision_rng"):
        nucleus._collision_rng = np.random.Generator(np.random.PCG64(73129))
    return nucleus._collision_rng


def _pair_type(is_i_proton: bool, is_j_proton: bool) -> str:
    if is_i_proton and is_j_proton:
        return "pp"
    if (not is_i_proton) and (not is_j_proton):
        return "nn"
    return "np"


def free_nn_cross_section(E_cm: float, pair_type: str) -> float:
    """Return free elastic NN cross section in mb.

    Uses a compact Cugnon-style beta parameterization at intermediate energies
    and the low-energy ImQMD cap used by Chen-Zhang-Li 2021 below 50 MeV lab.
    """

    pair = pair_type.lower()
    if pair not in {"pp", "nn", "np", "pn"}:
        raise ValueError("pair_type must be 'pp', 'nn', or 'np'")
    e_lab = max(0.0, 2.0 * float(E_cm))
    if e_lab < 50.0:
        return 180.0 if pair in {"np", "pn"} else 60.0

    total_energy = M_N + e_lab
    beta2 = max(1.0e-8, 1.0 - (M_N / total_energy) ** 2)
    beta = float(np.sqrt(beta2))
    if pair in {"pp", "nn"}:
        sigma = 13.73 - 15.04 / beta + 8.76 / (beta * beta) + 68.67 * beta**4
    else:
        sigma = -70.67 - 18.18 / beta + 25.26 / (beta * beta) + 113.85 * beta
    return float(np.clip(sigma, 1.0, 250.0))


def in_medium_factor(rho: float, is_pp_nn_or_np: str = "np", eta: float = 0.2, rho0: float = RHO0) -> float:
    """Return density reduction factor for the NN cross section."""

    del is_pp_nn_or_np
    return float(np.clip(1.0 - eta * max(0.0, float(rho)) / rho0, 0.2, 1.0))


def in_medium_nn_cross_section(E_cm: float, pair_type: str, rho: float, eta: float = 0.2) -> float:
    """Return ``sigma_NN_med = f_med * sigma_NN_free`` in mb."""

    return in_medium_factor(rho, pair_type, eta=eta) * free_nn_cross_section(E_cm, pair_type)


def _occupation_wigner(
    r_i: np.ndarray,
    p_final: np.ndarray,
    is_proton_i: bool,
    nucleus: ImQMDNucleus,
    exclude: tuple[int, int],
    width_scale: float = 1.0,
) -> float:
    positions = nucleus.positions
    momenta = nucleus.momenta
    same_species = nucleus.is_proton == bool(is_proton_i)
    if exclude:
        same_species[list(exclude)] = False
    if not np.any(same_species):
        return 0.0

    sigma_r = nucleus.sigma_r * width_scale
    sigma_p = HBAR_C / (2.0 * sigma_r)
    dr2 = np.sum((positions[same_species] - r_i) ** 2, axis=1)
    dp2 = np.sum((momenta[same_species] - p_final) ** 2, axis=1)
    weights = np.exp(-dr2 / (2.0 * sigma_r**2) - dp2 / (2.0 * sigma_p**2))

    # Explicit Wigner normalization: f_i = 1/(pi*hbar)^3 exp(...).
    # The phase-space volume per spin-degenerate nucleon state is h^3/4,
    # so the net occupation is 2 * sum(weights), matching the ImQMD spec.
    wigner_norm = 1.0 / (np.pi * HBAR_C) ** 3
    phase_space_cell = (2.0 * np.pi * HBAR_C) ** 3 / 4.0
    occupation = phase_space_cell * wigner_norm * np.sum(weights)
    return float(np.clip(occupation, 0.0, 1.0))


def pauli_blocking_probability(
    p1: np.ndarray,
    p2: np.ndarray,
    p3: np.ndarray,
    p4: np.ndarray,
    nucleus: ImQMDNucleus,
) -> float:
    """Standard QMD final-state blocking probability.

    ``p1`` and ``p2`` are the initial pair momenta and ``p3``/``p4`` the
    proposed final momenta.  If collision indices are not supplied through
    ``nucleus._active_collision_pair``, the closest current momenta are used.
    """

    momenta = nucleus.momenta
    if hasattr(nucleus, "_active_collision_pair"):
        i, j = nucleus._active_collision_pair
    else:
        i = int(np.argmin(np.linalg.norm(momenta - p1, axis=1)))
        j = int(np.argmin(np.linalg.norm(momenta - p2, axis=1)))
    positions = nucleus.positions
    species = nucleus.is_proton
    p_i = _occupation_wigner(positions[i], np.asarray(p3, dtype=float), bool(species[i]), nucleus, (i, j))
    p_j = _occupation_wigner(positions[j], np.asarray(p4, dtype=float), bool(species[j]), nucleus, (i, j))
    return float(np.clip(1.0 - (1.0 - p_i) * (1.0 - p_j), 0.0, 1.0))


def pauli_blocking_probability_v2(
    p1: np.ndarray,
    p2: np.ndarray,
    p3: np.ndarray,
    p4: np.ndarray,
    nucleus: ImQMDNucleus,
) -> float:
    """Smoother Chen-2024-inspired blocker using a broadened Wigner kernel."""

    momenta = nucleus.momenta
    if hasattr(nucleus, "_active_collision_pair"):
        i, j = nucleus._active_collision_pair
    else:
        i = int(np.argmin(np.linalg.norm(momenta - p1, axis=1)))
        j = int(np.argmin(np.linalg.norm(momenta - p2, axis=1)))
    positions = nucleus.positions
    species = nucleus.is_proton
    p_i = _occupation_wigner(positions[i], np.asarray(p3, dtype=float), bool(species[i]), nucleus, (i, j), 1.35)
    p_j = _occupation_wigner(positions[j], np.asarray(p4, dtype=float), bool(species[j]), nucleus, (i, j), 1.35)
    return float(np.clip(1.0 - (1.0 - p_i) * (1.0 - p_j), 0.0, 1.0))


def _occupation_wigner_fast(
    r_i: np.ndarray,
    p_final: np.ndarray,
    is_proton_i: bool,
    positions: np.ndarray,
    momenta: np.ndarray,
    is_proton: np.ndarray,
    sigma_r: float,
    exclude: tuple[int, int],
) -> float:
    same_species = is_proton == bool(is_proton_i)
    same_species[list(exclude)] = False
    if not np.any(same_species):
        return 0.0

    sigma_p = HBAR_C / (2.0 * sigma_r)
    dr2 = np.sum((positions[same_species] - r_i) ** 2, axis=1)
    dp2 = np.sum((momenta[same_species] - p_final) ** 2, axis=1)
    weights = np.exp(-dr2 / (2.0 * sigma_r**2) - dp2 / (2.0 * sigma_p**2))

    wigner_norm = 1.0 / (np.pi * HBAR_C) ** 3
    phase_space_cell = (2.0 * np.pi * HBAR_C) ** 3 / 4.0
    occupation = phase_space_cell * wigner_norm * np.sum(weights)
    return float(np.clip(occupation, 0.0, 1.0))


def _pauli_blocking_probability_fast(
    i: int,
    j: int,
    p3: np.ndarray,
    p4: np.ndarray,
    positions: np.ndarray,
    momenta: np.ndarray,
    is_proton: np.ndarray,
    rho_n: np.ndarray,
    rho_p: np.ndarray,
    sigma_r: float,
) -> float:
    p_i = _occupation_wigner_fast(
        positions[i],
        np.asarray(p3, dtype=float),
        bool(is_proton[i]),
        positions,
        momenta,
        is_proton,
        sigma_r,
        (i, j),
    )
    p_j = _occupation_wigner_fast(
        positions[j],
        np.asarray(p4, dtype=float),
        bool(is_proton[j]),
        positions,
        momenta,
        is_proton,
        sigma_r,
        (i, j),
    )
    return float(np.clip(1.0 - (1.0 - p_i) * (1.0 - p_j), 0.0, 1.0))


def _scatter_isotropic(p_i: np.ndarray, p_j: np.ndarray, rng: np.random.Generator) -> tuple[np.ndarray, np.ndarray]:
    p_cm = 0.5 * (p_i + p_j)
    q = 0.5 * (p_i - p_j)
    q_mag = float(np.linalg.norm(q))
    if q_mag <= 1.0e-12:
        return p_i.copy(), p_j.copy()
    direction = rng.normal(size=3)
    direction /= np.linalg.norm(direction)
    q_new = q_mag * direction
    return p_cm + q_new, p_cm - q_new


def _grid_cell(position: np.ndarray, cell_size: float) -> tuple[int, int, int]:
    """Return a stable integer grid cell for positive and negative coordinates."""

    return tuple(np.floor(np.asarray(position, dtype=float) / cell_size).astype(int))


def _build_grid(positions: np.ndarray, cell_size: float = NEIGHBOR_CELL_SIZE_FM) -> dict[tuple[int, int, int], list[int]]:
    """Assign each nucleon to a spatial grid cell."""

    grid: dict[tuple[int, int, int], list[int]] = {}
    for i, pos in enumerate(np.asarray(positions, dtype=float)):
        cell = _grid_cell(pos, cell_size)
        grid.setdefault(cell, []).append(i)
    return grid


def _neighbor_pairs(
    grid: dict[tuple[int, int, int], list[int]],
    positions: np.ndarray | None = None,
) -> list[tuple[int, int]]:
    """Return candidate pairs from the same or adjacent spatial cells."""

    del positions
    pairs: list[tuple[int, int]] = []
    seen: set[tuple[int, int]] = set()
    for cell, indices in grid.items():
        for di in (-1, 0, 1):
            for dj in (-1, 0, 1):
                for dk in (-1, 0, 1):
                    neighbor_cell = (cell[0] + di, cell[1] + dj, cell[2] + dk)
                    for j in grid.get(neighbor_cell, []):
                        for i in indices:
                            if i >= j:
                                continue
                            pair = (int(i), int(j))
                            if pair in seen:
                                continue
                            seen.add(pair)
                            pairs.append(pair)
    return pairs


def _cached_neighbor_pairs(nucleus: ImQMDNucleus, positions: np.ndarray) -> np.ndarray:
    """Return cached spatial-neighbor collision candidates for this step."""

    cache = getattr(nucleus, "_neighbor_list", None)
    call_index = int(getattr(nucleus, "_neighbor_list_call_index", 0)) + 1
    nucleus._neighbor_list_call_index = call_index

    rebuild = cache is None
    if not rebuild:
        cached_positions = cache["positions"]
        if cache["A"] != nucleus.A:
            rebuild = True
        elif call_index - cache["call_index"] >= NEIGHBOR_REBUILD_INTERVAL:
            rebuild = True
        else:
            displacement = np.linalg.norm(positions - cached_positions, axis=1)
            rebuild = bool(np.any(displacement > 0.5 * cache["cell_size"]))

    if rebuild:
        cell_size = NEIGHBOR_CELL_SIZE_FM
        grid = _build_grid(positions, cell_size=cell_size)
        pairs = np.asarray(_neighbor_pairs(grid, positions), dtype=np.int32)
        if pairs.size == 0:
            pairs = np.empty((0, 2), dtype=np.int32)
        else:
            pairs = pairs.reshape(-1, 2)
        cache = {
            "A": nucleus.A,
            "cell_size": cell_size,
            "call_index": call_index,
            "positions": positions.copy(),
            "grid": grid,
            "pairs": pairs,
            "rebuilds": int(cache.get("rebuilds", 0) + 1) if cache else 1,
        }
        nucleus._neighbor_list = cache
    return cache["pairs"]


def attempt_nn_collision(nucleus: ImQMDNucleus, dt: float) -> dict[str, int]:
    """Find geometric NN collision candidates and apply Pauli-blocked scattering."""

    rng = _rng(nucleus)
    positions = nucleus.positions
    momenta = nucleus.momenta
    is_proton = nucleus.is_proton
    rho, rho_n, rho_p = nucleus.centroid_densities()
    pairs = _cached_neighbor_pairs(nucleus, positions)
    order = np.arange(pairs.shape[0])
    rng.shuffle(order)
    used: set[int] = set()
    stats = {"attempted": 0, "blocked": 0, "accepted": 0}

    for pair_index in order:
        i, j = pairs[int(pair_index)]
        i = int(i)
        j = int(j)
        if i in used or j in used:
            continue
        r_rel = positions[i] - positions[j]
        p_rel_vec = momenta[i] - momenta[j]
        v_rel_vec = p_rel_vec / M_N
        v2 = float(np.dot(v_rel_vec, v_rel_vec))
        if v2 <= 1.0e-12:
            continue
        t_closest = float(np.clip(-np.dot(r_rel, v_rel_vec) / v2, 0.0, dt))
        if np.dot(r_rel, v_rel_vec) > 0.0 and t_closest <= 1.0e-12:
            continue
        d_min = float(np.linalg.norm(r_rel + v_rel_vec * t_closest))
        e_cm = float(np.dot(0.5 * p_rel_vec, 0.5 * p_rel_vec) / M_N)
        pair = _pair_type(bool(is_proton[i]), bool(is_proton[j]))
        local_rho = 0.5 * (rho[i] + rho[j])
        sigma_mb = in_medium_nn_cross_section(e_cm, pair, local_rho)
        d_collision = np.sqrt(sigma_mb * MB_TO_FM2 / np.pi)
        if d_min >= d_collision:
            continue

        stats["attempted"] += 1
        p3, p4 = _scatter_isotropic(momenta[i], momenta[j], rng)
        p_block = _pauli_blocking_probability_fast(
            i,
            j,
            p3,
            p4,
            positions,
            momenta,
            is_proton,
            rho_n,
            rho_p,
            nucleus.sigma_r,
        )
        if rng.random() < p_block:
            stats["blocked"] += 1
            continue
        momenta[i] = p3
        momenta[j] = p4
        used.update((i, j))
        stats["accepted"] += 1

    nucleus.momenta = momenta - np.mean(momenta, axis=0)
    existing = getattr(nucleus, "collision_stats", {"attempted": 0, "blocked": 0, "accepted": 0})
    nucleus.collision_stats = {key: int(existing.get(key, 0) + stats[key]) for key in stats}
    return stats


def fermi_constraint_check(
    nucleus: ImQMDNucleus,
    threshold: float = 255.0,
    group_ids: np.ndarray | None = None,
) -> int:
    """Apply a CoMD-style phase-space constraint to nucleon pairs.

    When ``group_ids`` is provided, only pairs from the *same* group are
    constrained, allowing nucleon exchange across groups during collisions.
    """
    rng = _rng(nucleus)
    positions = nucleus.positions
    momenta = nucleus.momenta
    kinetic_before = float(np.sum(momenta * momenta))
    corrections = 0
    if nucleus.A < 2:
        return 0

    iu, ju = np.triu_indices(nucleus.A, 1)
    if group_ids is not None:
        same_group = group_ids[iu] == group_ids[ju]
        iu, ju = iu[same_group], ju[same_group]
        if len(iu) == 0:
            return 0
    dr_vec = positions[iu] - positions[ju]
    dp_all = momenta[iu] - momenta[ju]
    dr_all = np.linalg.norm(dr_vec, axis=1)
    dp_norm = np.linalg.norm(dp_all, axis=1)
    violating = np.flatnonzero(dr_all * dp_norm < threshold)

    for pair_index in violating:
        i = int(iu[pair_index])
        j = int(ju[pair_index])
        dr = float(np.linalg.norm(positions[i] - positions[j]))
        dp_vec = momenta[i] - momenta[j]
        dp = float(np.linalg.norm(dp_vec))
        if dr * dp >= threshold:
            continue
        p_cm = 0.5 * (momenta[i] + momenta[j])
        if dp > 1.0e-12:
            direction = dp_vec / dp
        else:
            direction = rng.normal(size=3)
            norm = np.linalg.norm(direction)
            direction = np.array([1.0, 0.0, 0.0]) if norm == 0.0 else direction / norm
        target_dp = threshold / max(dr, 0.25)
        q_new = 0.5 * target_dp * direction
        momenta[i] = p_cm + q_new
        momenta[j] = p_cm - q_new
        corrections += 1

    if corrections:
        momenta = momenta - np.mean(momenta, axis=0)
        kinetic_after = float(np.sum(momenta * momenta))
        if kinetic_before > 0.0 and kinetic_after > 0.0:
            momenta *= np.sqrt(kinetic_before / kinetic_after)
        nucleus.momenta = momenta
    return corrections


__all__ = [
    "attempt_nn_collision",
    "fermi_constraint_check",
    "free_nn_cross_section",
    "in_medium_factor",
    "in_medium_nn_cross_section",
    "pauli_blocking_probability",
    "pauli_blocking_probability_v2",
]
