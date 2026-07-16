"""Impact-parameter event scans and fragment-yield cross sections for ImQMD."""

from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass
from functools import lru_cache
import multiprocessing as mp

import numpy as np

from .decay import evaporate_full
from .fragments import Fragment, isospin_mst, minimum_spanning_tree
from .initializer import _wang_skin_radii, initialize_and_relax
from .nucleus import GaussianPacket, ImQMDNucleus
from .propagator import propagate
from .skyrme import M_N, SkyrmeEDF


@dataclass(frozen=True)
class EventFragmentRecord:
    Z: int
    A: int
    excitation_energy: float
    final_Z: int
    final_A: int
    e_lab: float = 0.0
    theta_lab: float = 0.0


@dataclass(frozen=True)
class ImpactParameterEvent:
    b: float
    event_index: int
    primary_fragments: tuple[EventFragmentRecord, ...]
    accepted_collisions: int
    attempted_collisions: int


@dataclass(frozen=True)
class ImpactParameterResult:
    b: float
    delta_b: float
    n_events: int
    events: tuple[ImpactParameterEvent, ...]
    sigma_by_z: dict[int, float]
    total_sigma: float


@dataclass(frozen=True)
class CrossSectionScanResult:
    projectile: tuple[int, int]
    target: tuple[int, int]
    energy_per_a: float
    b_max: float
    fragment_method: str
    impact_parameter_results: tuple[ImpactParameterResult, ...]
    dsigma_dz: dict[int, float]
    total_cross_section: float
    dsigma_da: dict = None  # ponytail: mass distribution, populated post-init
    d2sigma: tuple | None = None  # (theta_grid, e_grid, array)


@lru_cache(maxsize=32)
def _relaxed_template(
    Z: int,
    A: int,
    sigma_r: float,
    seed: int,
    relax_time: float,
) -> ImQMDNucleus:
    return initialize_and_relax(Z, A, sigma_r=sigma_r, seed=seed, relax_time=relax_time)


def _prepared_nucleus(
    Z: int,
    A: int,
    sigma_r: float,
    seed: int,
    relax_time: float,
    edf: SkyrmeEDF | None,
) -> ImQMDNucleus:
    if edf is None:
        return _relaxed_template(Z, A, sigma_r, seed, relax_time).copy()
    return initialize_and_relax(Z, A, sigma_r=sigma_r, seed=seed, relax_time=relax_time, edf=edf)


def estimate_grazing_bmax(projectile_z: int, projectile_a: int, target_z: int, target_a: int, padding: float = 1.0) -> float:
    """Estimate a grazing impact-parameter cutoff from skin-adjusted radii."""

    _, rp_p, rn_p = _wang_skin_radii(int(projectile_z), int(projectile_a))
    _, rp_t, rn_t = _wang_skin_radii(int(target_z), int(target_a))
    projectile_radius = max(rp_p, rn_p) + 0.8
    target_radius = max(rp_t, rn_t) + 0.8
    return float(projectile_radius + target_radius + float(padding))


def collision_separation(projectile_z: int, projectile_a: int, target_z: int, target_a: int) -> float:
    """Return a practical initial center-to-center separation in fm."""

    return float(max(28.0, estimate_grazing_bmax(projectile_z, projectile_a, target_z, target_a) + 16.0))


def make_collision_event(
    projectile_z: int,
    projectile_a: int,
    target_z: int,
    target_a: int,
    energy_per_a: float,
    impact_parameter: float,
    projectile_seed: int,
    target_seed: int,
    collision_seed: int,
    separation: float | None = None,
    sigma_r: float = 1.1,
    relax_time: float = 800.0,
    edf: SkyrmeEDF | None = None,
) -> ImQMDNucleus:
    """Build a two-nucleus ImQMD event for a given impact parameter."""

    separation = collision_separation(projectile_z, projectile_a, target_z, target_a) if separation is None else float(separation)
    projectile = _prepared_nucleus(projectile_z, projectile_a, float(sigma_r), int(projectile_seed), float(relax_time), edf)
    target = _prepared_nucleus(target_z, target_a, float(sigma_r), int(target_seed), float(relax_time), edf)
    edf = edf or projectile.edf
    projectile.edf = edf
    target.edf = edf
    p_beam = float(np.sqrt(2.0 * M_N * float(energy_per_a)))

    packets: list[GaussianPacket] = []
    reference_positions: list[np.ndarray] = []
    reference_group_ids: list[int] = []

    for packet in projectile.packets:
        shift = np.array([-0.5 * separation, 0.5 * impact_parameter, 0.0], dtype=float)
        packets.append(
            GaussianPacket(
                packet.r_i + shift,
                packet.p_i + np.array([p_beam, 0.0, 0.0]),
                packet.sigma_r,
                packet.is_proton,
            )
        )
        reference_positions.append(packet.r_i + shift)
        reference_group_ids.append(0)

    for packet in target.packets:
        shift = np.array([0.5 * separation, -0.5 * impact_parameter, 0.0], dtype=float)
        packets.append(
            GaussianPacket(
                packet.r_i + shift,
                packet.p_i + np.array([-p_beam, 0.0, 0.0]),
                packet.sigma_r,
                packet.is_proton,
            )
        )
        reference_positions.append(packet.r_i + shift)
        reference_group_ids.append(1)

    system = ImQMDNucleus(
        projectile_z + target_z,
        (projectile_a - projectile_z) + (target_a - target_z),
        packets,
        edf=edf,
        reference_positions=np.asarray(reference_positions, dtype=float),
    )
    system.grid_nuclear_scale = float(
        (
            projectile.A * float(getattr(projectile, "grid_nuclear_scale", 1.0))
            + target.A * float(getattr(target, "grid_nuclear_scale", 1.0))
        )
        / max(projectile.A + target.A, 1)
    )
    system.reference_group_ids = np.asarray(reference_group_ids, dtype=int)
    system._collision_rng = np.random.default_rng(int(collision_seed))
    system.collision_stats = {"attempted": 0, "blocked": 0, "accepted": 0}
    return system


def identify_fragments(
    nucleus: ImQMDNucleus,
    method: str = "iso-mst",
    p_cut: float | None = 250.0,
    iso_r_cut_pp: float = 3.0,
    iso_r_cut_nn: float = 6.0,
    iso_r_cut_np: float = 6.0,
) -> list[Fragment]:
    """Identify primary fragments with MST or iso-MST."""

    method = method.lower()
    if method == "iso-mst":
        fragments = isospin_mst(
            nucleus,
            r_cut_pp=iso_r_cut_pp,
            r_cut_nn=iso_r_cut_nn,
            r_cut_np=iso_r_cut_np,
            p_cut=p_cut,
        )
    elif method == "mst":
        fragments = minimum_spanning_tree(nucleus, p_cut=p_cut)
    else:
        raise ValueError(f"Unsupported fragment method: {method}")
    return sorted(fragments, key=lambda fragment: fragment.A, reverse=True)


def _record_fragment(
    fragment: Fragment,
    rng: np.random.Generator,
    use_hivap: bool = False,
) -> EventFragmentRecord:
    if use_hivap:
        from .hivap_wrapper import sample_residue

        final_Z, final_A = sample_residue(
            fragment.Z, fragment.A, fragment.excitation_energy, rng=rng
        )
    else:
        final_Z, final_A = evaporate_full(
            fragment.Z, fragment.A, fragment.excitation_energy, rng=rng
        )
    return EventFragmentRecord(
        Z=int(fragment.Z),
        A=int(fragment.A),
        excitation_energy=float(fragment.excitation_energy),
        final_Z=int(final_Z),
        final_A=int(final_A),
    )


def run_imqmd_event(
    projectile_z: int,
    projectile_a: int,
    target_z: int,
    target_a: int,
    energy_per_a: float,
    impact_parameter: float,
    event_index: int,
    projectile_seed: int = 8600,
    target_seed: int = 6400,
    collision_seed: int = 15000,
    decay_seed: int = 25000,
    separation: float | None = None,
    sigma_r: float = 1.1,
    relax_time: float = 800.0,
    time_fm_c: float = 1000.0,
    dt: float = 1.0,
    collision_dt: float = 1.0,
    fragment_method: str = "mst",
    p_cut: float | None = 250.0,
    iso_r_cut_pp: float = 3.0,
    iso_r_cut_nn: float = 6.0,
    iso_r_cut_np: float = 6.0,
    edf: SkyrmeEDF | None = None,
    resample_initial_nuclei: bool = False,
    use_grid_edf: bool = True,
    use_hivap: bool = False,
) -> ImpactParameterEvent:
    """Run one ImQMD event through fragment recognition and de-excitation."""

    seed_offset = int(event_index + round(10.0 * impact_parameter))
    projectile_seed_value = int(projectile_seed + 1000 * seed_offset) if resample_initial_nuclei else int(projectile_seed)
    target_seed_value = int(target_seed + 1000 * seed_offset) if resample_initial_nuclei else int(target_seed)
    system = make_collision_event(
        projectile_z=projectile_z,
        projectile_a=projectile_a,
        target_z=target_z,
        target_a=target_a,
        energy_per_a=energy_per_a,
        impact_parameter=impact_parameter,
        projectile_seed=projectile_seed_value,
        target_seed=target_seed_value,
        collision_seed=int(collision_seed + 1000 * seed_offset),
        separation=separation,
        sigma_r=sigma_r,
        relax_time=relax_time,
        edf=edf,
    )
    propagate(
        system,
        dt=float(dt),
        n_steps=max(1, int(round(float(time_fm_c) / float(dt)))),
        sample_every=max(1, int(round(100.0 / float(dt)))),
        with_collisions=True,
        collision_dt=float(collision_dt),
        remove_cm_drift=False,
        use_surface_term=True,
        use_static_stabilizer=False,
        use_grid_edf=use_grid_edf,
    )
    # ponytail: 100 fm/c cooling to reduce fragment multiplicity
    propagate(system, dt=float(dt), n_steps=50, sample_every=1000,
              with_collisions=False, remove_cm_drift=False,
              use_surface_term=True, use_static_stabilizer=False,
              use_grid_edf=use_grid_edf)
    fragments = identify_fragments(
        system,
        method=fragment_method,
        p_cut=p_cut,
        iso_r_cut_pp=iso_r_cut_pp,
        iso_r_cut_nn=iso_r_cut_nn,
        iso_r_cut_np=iso_r_cut_np,
    )
    decay_rng = np.random.default_rng(int(decay_seed + 1000 * seed_offset))
    # Two-body kinematics: fragment lab energy and angle
    v_cm = np.array([np.sqrt(2.0 * energy_per_a * M_N) * projectile_a / (projectile_a + target_a), 0.0, 0.0])
    raw_records = []
    for fragment in fragments:
        rec = _record_fragment(fragment, decay_rng, use_hivap=use_hivap)
        v_f_cm = fragment.momentum / (fragment.A * M_N)
        v_lab = v_cm + v_f_cm
        e_lab = 0.5 * fragment.A * M_N * np.dot(v_lab, v_lab)
        theta_lab = np.degrees(np.arctan2(np.linalg.norm(v_lab[1:]), v_lab[0]))
        raw_records.append(EventFragmentRecord(
            Z=rec.Z, A=rec.A, excitation_energy=rec.excitation_energy,
            final_Z=rec.final_Z, final_A=rec.final_A,
            e_lab=float(e_lab), theta_lab=float(theta_lab)))
    records = tuple(raw_records)
    collision_stats = getattr(system, "collision_stats", {})
    return ImpactParameterEvent(
        b=float(impact_parameter),
        event_index=int(event_index),
        primary_fragments=records,
        accepted_collisions=int(collision_stats.get("accepted", 0)),
        attempted_collisions=int(collision_stats.get("attempted", 0)),
    )


def _b_bin_widths(b_values: np.ndarray, b_max: float) -> np.ndarray:
    b_values = np.asarray(sorted(float(value) for value in b_values), dtype=float)
    if b_values.size == 0:
        return np.empty(0, dtype=float)
    if b_values.size == 1:
        return np.asarray([max(float(b_max), 1.0)], dtype=float)
    edges = np.empty(b_values.size + 1, dtype=float)
    edges[0] = 0.0
    edges[1:-1] = 0.5 * (b_values[:-1] + b_values[1:])
    edges[-1] = max(float(b_max), edges[-2])
    widths = np.diff(edges)
    return np.maximum(widths, 1.0e-6)


def _worker_run_event(kwargs: dict[str, object]) -> ImpactParameterEvent:
    return run_imqmd_event(**kwargs)


def impact_parameter_scan(
    projectile_z: int,
    projectile_a: int,
    target_z: int,
    target_a: int,
    energy_per_a: float,
    b_values: list[float] | tuple[float, ...] | np.ndarray | None = None,
    b_max: float | None = None,
    delta_b: float = 1.0,
    events_per_b: int = 10,
    n_workers: int = 1,
    projectile_seed: int = 8600,
    target_seed: int = 6400,
    collision_seed: int = 15000,
    decay_seed: int = 25000,
    separation: float | None = None,
    sigma_r: float = 1.1,
    relax_time: float = 800.0,
    time_fm_c: float = 1000.0,
    dt: float = 1.0,
    collision_dt: float = 1.0,
    fragment_method: str = "mst",
    p_cut: float | None = 250.0,
    edf: SkyrmeEDF | None = None,
    resample_initial_nuclei: bool = False,
    use_grid_edf: bool = True,
    use_hivap: bool = False,
) -> CrossSectionScanResult:
    """Run an impact-parameter scan and accumulate ``dσ/dZ``."""

    estimated_bmax = estimate_grazing_bmax(projectile_z, projectile_a, target_z, target_a)
    b_max = float(estimated_bmax if b_max is None else b_max)
    if b_values is None:
        b_values = np.arange(0.5 * delta_b, b_max + 0.5 * delta_b, delta_b, dtype=float)
    b_values = np.asarray(sorted(float(value) for value in b_values), dtype=float)
    widths = _b_bin_widths(b_values, b_max)

    tasks: list[dict[str, object]] = []
    for b in b_values:
        for event_index in range(int(events_per_b)):
            tasks.append(
                {
                    "projectile_z": int(projectile_z),
                    "projectile_a": int(projectile_a),
                    "target_z": int(target_z),
                    "target_a": int(target_a),
                    "energy_per_a": float(energy_per_a),
                    "impact_parameter": float(b),
                    "event_index": int(event_index),
                    "projectile_seed": int(projectile_seed),
                    "target_seed": int(target_seed),
                    "collision_seed": int(collision_seed),
                    "decay_seed": int(decay_seed),
                    "separation": separation,
                    "sigma_r": float(sigma_r),
                    "relax_time": float(relax_time),
                    "time_fm_c": float(time_fm_c),
                    "dt": float(dt),
                    "collision_dt": float(collision_dt),
                    "fragment_method": fragment_method,
                    "p_cut": p_cut,
                    "edf": edf,
                    "resample_initial_nuclei": bool(resample_initial_nuclei),
                    "use_grid_edf": bool(use_grid_edf),
                    "use_hivap": bool(use_hivap),
                }
            )

    if int(n_workers) > 1 and len(tasks) > 1:
        ctx = mp.get_context("fork")
        with ctx.Pool(processes=int(n_workers)) as pool:
            event_results = list(pool.imap_unordered(_worker_run_event, tasks, chunksize=1))
    else:
        event_results = [_worker_run_event(task) for task in tasks]

    by_b: dict[float, list[ImpactParameterEvent]] = defaultdict(list)
    for result in event_results:
        by_b[float(result.b)].append(result)
    ordered_results: list[ImpactParameterResult] = []
    sigma_total_by_z: dict[int, float] = defaultdict(float)
    sigma_total_by_a: dict[int, float] = defaultdict(float)
    total_cross_section = 0.0

    # ponytail: d2sigma binning
    theta_grid = np.linspace(0, 90, 46)
    e_grid = np.linspace(0, 1500, 51)
    d2 = np.zeros((45, 50), dtype=float)

    for b, width in zip(b_values, widths):
        events = sorted(by_b.get(float(b), []), key=lambda item: item.event_index)
        sigma_by_z: dict[int, float] = defaultdict(float)
        if events:
            weight = 2.0 * np.pi * float(b) * float(width) / float(len(events))
            for event in events:
                for fragment in event.primary_fragments:
                    sigma_by_z[int(fragment.final_Z)] += weight
                    sigma_total_by_z[int(fragment.final_Z)] += weight
                    sigma_total_by_a[int(fragment.final_A)] += weight
                    total_cross_section += weight
                    if fragment.e_lab > 0 and 0 <= fragment.theta_lab <= 90:
                        it = np.clip(np.digitize(fragment.theta_lab, theta_grid) - 1, 0, 44)
                        ie = np.clip(np.digitize(fragment.e_lab, e_grid) - 1, 0, 49)
                        d2[it, ie] += weight
        ordered_results.append(
            ImpactParameterResult(
                b=float(b),
                delta_b=float(width),
                n_events=len(events),
                events=tuple(events),
                sigma_by_z=dict(sorted(sigma_by_z.items())),
                total_sigma=float(sum(sigma_by_z.values())),
            )
        )

    d2sigma = (tuple(theta_grid.tolist()), tuple(e_grid.tolist()), tuple(d2.ravel().tolist()))
    return CrossSectionScanResult(
        projectile=(int(projectile_z), int(projectile_a)),
        target=(int(target_z), int(target_a)),
        energy_per_a=float(energy_per_a),
        b_max=float(b_max),
        fragment_method=str(fragment_method),
        impact_parameter_results=tuple(ordered_results),
        dsigma_dz=dict(sorted(sigma_total_by_z.items())),
        total_cross_section=float(total_cross_section),
        dsigma_da=dict(sorted(sigma_total_by_a.items())),
        d2sigma=d2sigma,
    )


__all__ = [
    "CrossSectionScanResult",
    "EventFragmentRecord",
    "ImpactParameterEvent",
    "ImpactParameterResult",
    "collision_separation",
    "estimate_grazing_bmax",
    "identify_fragments",
    "impact_parameter_scan",
    "make_collision_event",
    "run_imqmd_event",
]
