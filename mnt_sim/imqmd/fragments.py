"""Fragment recognition and excitation estimates for ImQMD events."""

from __future__ import annotations

from dataclasses import dataclass
from itertools import combinations

import numpy as np

from .nucleus import GaussianPacket, ImQMDNucleus
_COLD_GROUND_STATE_CACHE: dict[tuple[int, int, float, str, float], float] = {}


@dataclass
class Fragment:
    """Recognized primary fragment."""

    Z: int
    A: int
    position: np.ndarray
    momentum: np.ndarray
    excitation_energy: float
    nucleon_indices: list[int]


def _connected_components(adjacency: np.ndarray) -> list[np.ndarray]:
    adjacency = np.asarray(adjacency, dtype=bool)
    n = adjacency.shape[0]
    if n == 0:
        return []
    try:
        from scipy.sparse import csr_matrix
        from scipy.sparse.csgraph import connected_components

        count, labels = connected_components(csr_matrix(adjacency), directed=False)
        return [np.where(labels == label)[0] for label in range(count)]
    except Exception:  # pragma: no cover - scipy is expected in normal runs
        seen = np.zeros(n, dtype=bool)
        components: list[np.ndarray] = []
        for start in range(n):
            if seen[start]:
                continue
            stack = [start]
            seen[start] = True
            component: list[int] = []
            while stack:
                i = stack.pop()
                component.append(i)
                for j in np.where(adjacency[i] & ~seen)[0]:
                    seen[j] = True
                    stack.append(int(j))
            components.append(np.asarray(component, dtype=int))
        return components


def _fragment_from_indices(nucleus: ImQMDNucleus, indices: np.ndarray | list[int]) -> Fragment:
    idx = np.asarray(indices, dtype=int)
    is_proton = nucleus.is_proton[idx]
    positions = nucleus.positions[idx]
    momenta = nucleus.momenta[idx]
    Z = int(np.sum(is_proton))
    A = int(len(idx))
    position = np.mean(positions, axis=0) if A else np.zeros(3)
    momentum = np.sum(momenta, axis=0) if A else np.zeros(3)
    fragment = Fragment(
        Z=Z,
        A=A,
        position=position,
        momentum=momentum,
        excitation_energy=0.0,
        nucleon_indices=[int(i) for i in idx],
    )
    fragment.excitation_energy = compute_fragment_excitation(nucleus, fragment)
    return fragment


def _pair_distances(values: np.ndarray) -> np.ndarray:
    diff = values[:, None, :] - values[None, :, :]
    return np.linalg.norm(diff, axis=-1)


def minimum_spanning_tree(
    nucleus: ImQMDNucleus,
    r_cut: float = 3.5,
    p_cut: float | None = 300.0,
) -> list[Fragment]:
    """Recognize fragments with the QMD minimum-spanning-tree rule.

    Despite the historical name, the algorithm returns connected components of
    the phase-space proximity graph.  p_cut default 300 MeV/c follows the
    imQMD MNT criterion (research report Sec. 4.5); the classic QMD value is
    250 MeV/c (Aichelin 1991, Zhang 2020 review).
    """

    n = nucleus.A
    if n == 0:
        return []
    dr = _pair_distances(nucleus.positions)
    adjacency = dr < float(r_cut)
    if p_cut is not None and np.isfinite(p_cut):
        dp = _pair_distances(nucleus.momenta)
        adjacency &= dp < float(p_cut)
    np.fill_diagonal(adjacency, True)
    return [_fragment_from_indices(nucleus, component) for component in _connected_components(adjacency)]


def isospin_mst(
    nucleus: ImQMDNucleus,
    r_cut_pp: float = 3.0,
    r_cut_nn: float = 6.0,
    r_cut_np: float = 6.0,
    p_cut: float | None = 300.0,
) -> list[Fragment]:
    """Recognize fragments with isospin-dependent coordinate cutoffs.

    The default cutoffs follow the iso-MST methodology used in imQMD MNT
    studies (Li+2019 PRC 99, 034619; Li+2020 PLB 808, 135697): pp=3.0 fm,
    nn=np=6.0 fm.  The factor-2 asymmetry between pp and nn/np reflects the
    physical picture in near-barrier MNT where neutrons decouple from the
    reaction dynamics later than protons (the lighter partner tends to emit
    free neutrons that the MST should NOT merge into the heavy residue).
    The 6.0 fm nn cutoff acts as a regulator: neutrons within this radius of
    a fragment are assigned to it, preventing artefactual free-neutron
    emission from the MST cut while still allowing genuine neutron evaporation
    at the statistical-decay stage.  For sensitivity studies, the Li+2019
    values of ~pp=3.0, nn=3.8, np=3.4 fm provide a tighter alternative.

    p_cut default 300 MeV/c (research report Sec. 4.5 MNT criterion).
    """

    n = nucleus.A
    if n == 0:
        return []
    is_proton = nucleus.is_proton
    dr = _pair_distances(nucleus.positions)
    pp = is_proton[:, None] & is_proton[None, :]
    nn = (~is_proton[:, None]) & (~is_proton[None, :])
    cutoffs = np.where(pp, r_cut_pp, np.where(nn, r_cut_nn, r_cut_np))
    adjacency = dr < cutoffs
    if p_cut is not None and np.isfinite(p_cut):
        dp = _pair_distances(nucleus.momenta)
        adjacency &= dp < float(p_cut)
    np.fill_diagonal(adjacency, True)
    return [_fragment_from_indices(nucleus, component) for component in _connected_components(adjacency)]


def adaptive_mst(
    nucleus: ImQMDNucleus,
    base_r_cut: float = 3.2,
    max_r_cut: float = 4.0,
    p_cut: float | None = None,
) -> list[Fragment]:
    """Recognize reaction fragments with an A-scaled coordinate cutoff.

    The first pass estimates compact seeds, then the second pass relaxes the
    coordinate cutoff for heavy, deformed residues without changing the
    standard MST defaults used by static tests.
    """

    n = nucleus.A
    if n == 0:
        return []
    seed_fragments = minimum_spanning_tree(nucleus, r_cut=base_r_cut, p_cut=p_cut)
    local_a = np.ones(n, dtype=float)
    for fragment in seed_fragments:
        local_a[fragment.nucleon_indices] = max(float(fragment.A), 1.0)
    pair_a = np.maximum(local_a[:, None], local_a[None, :])
    cutoffs = np.minimum(float(max_r_cut), float(base_r_cut) * np.cbrt(pair_a) / np.cbrt(16.0))
    cutoffs = np.maximum(cutoffs, float(base_r_cut))
    adjacency = _pair_distances(nucleus.positions) < cutoffs
    if p_cut is not None and np.isfinite(p_cut):
        adjacency &= _pair_distances(nucleus.momenta) < float(p_cut)
    np.fill_diagonal(adjacency, True)
    fragments = [_fragment_from_indices(nucleus, component) for component in _connected_components(adjacency)]
    return sorted(fragments, key=lambda fragment: fragment.A, reverse=True)


def reaction_fragments(nucleus: ImQMDNucleus, p_cut: float | None = None) -> list[Fragment]:
    """Convenience fragment finder with heavy-ion reaction defaults."""

    fragments = adaptive_mst(nucleus, base_r_cut=3.5, max_r_cut=4.0, p_cut=p_cut)
    return sorted(fragments, key=lambda fragment: fragment.A, reverse=True)


def _matches_template(is_proton: np.ndarray, indices: tuple[int, ...], z_req: int, a_req: int) -> bool:
    if len(indices) != a_req:
        return False
    return int(np.sum(is_proton[list(indices)])) == z_req


def _phase_space_compact(
    positions: np.ndarray,
    momenta: np.ndarray,
    indices: tuple[int, ...],
    r_cut: float,
    p_cut: float,
) -> bool:
    if len(indices) <= 1:
        return True
    sub_r = positions[list(indices)]
    sub_p = momenta[list(indices)]
    dr = _pair_distances(sub_r)
    dp = _pair_distances(sub_p)
    iu = np.triu_indices(len(indices), 1)
    return bool(np.all(dr[iu] <= r_cut) and np.all(dp[iu] <= p_cut))


def coalescence_light(
    nucleus: ImQMDNucleus,
    r_cut: float = 2.4,
    p_cut: float = 180.0,
    neighbor_limit: int = 14,
) -> list[Fragment]:
    """Direct phase-space coalescence for d, t, 3He, and 4He.

    The search is greedy with the usual heavy-to-light priority.  It returns
    only the light clusters found by this specialized pass.
    """

    positions = nucleus.positions
    momenta = nucleus.momenta
    is_proton = nucleus.is_proton
    unused: set[int] = set(range(nucleus.A))
    fragments: list[Fragment] = []
    templates = (
        (2, 4),  # 4He
        (2, 3),  # 3He
        (1, 3),  # t
        (1, 2),  # d
    )

    for seed in range(nucleus.A):
        if seed not in unused:
            continue
        distances = np.linalg.norm(positions - positions[seed], axis=1)
        nearby = [int(i) for i in np.argsort(distances) if i in unused and distances[i] <= r_cut]
        nearby = nearby[:neighbor_limit]
        formed: tuple[int, ...] | None = None
        for z_req, a_req in templates:
            if len(nearby) < a_req:
                continue
            candidate_iter = combinations(nearby, a_req)
            for candidate in candidate_iter:
                if seed not in candidate:
                    continue
                if not _matches_template(is_proton, candidate, z_req, a_req):
                    continue
                if not _phase_space_compact(positions, momenta, candidate, r_cut, p_cut):
                    continue
                formed = tuple(sorted(candidate))
                break
            if formed is not None:
                break
        if formed is None:
            continue
        for i in formed:
            unused.remove(i)
        fragments.append(_fragment_from_indices(nucleus, np.asarray(formed, dtype=int)))

    return fragments


def compute_fragment_excitation(nucleus: ImQMDNucleus, fragment: Fragment) -> float:
    """Estimate ``E* = E_int(fragment) - E_ground(Z,A)`` in MeV."""

    if fragment.A <= 1:
        return 0.0

    idx = np.asarray(fragment.nucleon_indices, dtype=int)
    positions = nucleus.positions[idx]
    momenta = nucleus.momenta[idx]
    is_proton = nucleus.is_proton[idx]

    p_cm = np.mean(momenta, axis=0)
    r_cm = np.mean(positions, axis=0)
    internal_positions = positions - r_cm
    internal_momenta = momenta - p_cm

    from .grid_edf import GridEDF

    fragment_scale = fragment_ground_state_scale(fragment.Z, fragment.A, nucleus)
    grid = GridEDF(
        nucleus.edf.parameters,
        nucleus.packet_sigmas[idx],
        grid_spacing=1.0,
        nuclear_scale=fragment_scale,
    )
    internal_energy = grid.total_energy(
        internal_positions,
        internal_momenta,
        is_proton,
        use_surface_term=True,
    )["total"]

    # Event records do not carry J, so rotational energy remains part of E*.
    cold_energy = cold_ground_state_energy(fragment.Z, fragment.A, nucleus)
    return float(internal_energy - cold_energy)


def cold_ground_state_energy(Z: int, A: int, reference: ImQMDNucleus) -> float:
    """Return the cached cold-nucleus energy from the same calibrated EDF."""

    if A <= 1:
        return 0.0
    sigma_r = float(reference.sigma_r)
    edf = reference.edf
    cache_key = (
        int(Z),
        int(A),
        round(sigma_r, 12),
        edf.parameters.name,
        float(edf.static_k),
    )
    cached = _COLD_GROUND_STATE_CACHE.get(cache_key)
    if cached is not None:
        return cached

    from .initializer import initialize_nucleus
    from .grid_edf import GridEDF

    cold = initialize_nucleus(Z, A, sigma_r=sigma_r, seed=1000 + 17 * int(Z) + int(A), edf=edf)
    # The event uses its inherited scale; the model ground state uses the
    # scale fitted when that cold nucleus was initialized with the same EDF.
    scale = float(getattr(cold, "grid_nuclear_scale", 1.0))
    grid = GridEDF(edf.parameters, sigma_r, grid_spacing=1.0, nuclear_scale=scale)
    energy = grid.total_energy(
        cold.positions, cold.momenta, cold.is_proton, use_surface_term=True
    )["total"]
    _COLD_GROUND_STATE_CACHE[cache_key] = float(energy)
    return float(energy)


def fragment_ground_state_scale(Z: int, A: int, reference: ImQMDNucleus) -> float:
    """Use the scale carried by the initialized/evolved system."""

    if A <= 1:
        return 1.0
    return float(getattr(reference, "grid_nuclear_scale", 1.0))


__all__ = [
    "Fragment",
    "adaptive_mst",
    "coalescence_light",
    "compute_fragment_excitation",
    "cold_ground_state_energy",
    "fragment_ground_state_scale",
    "isospin_mst",
    "minimum_spanning_tree",
    "reaction_fragments",
]
