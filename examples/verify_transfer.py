#!/usr/bin/env python
"""Verify asymmetric ImQMD transfer cases and the heavy-system neighbor list.

Run:
    python examples/verify_transfer.py
    python examples/verify_transfer.py --include-u-benchmark
"""

from __future__ import annotations

import argparse
import os
import signal
import sys
import time
from collections import Counter
from dataclasses import dataclass

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from mnt_sim.imqmd import GaussianPacket, ImQMDNucleus, initialize_and_relax, propagate, reaction_fragments  # noqa: E402


M_N_TEST = 938.9
EVENT_TIMEOUT_S = 300


class EventTimeout(RuntimeError):
    """Raised when one verification event exceeds the interactive time budget."""


def _timeout_handler(signum, frame) -> None:
    del signum, frame
    raise EventTimeout(f"event exceeded {EVENT_TIMEOUT_S} s")


@dataclass(frozen=True)
class SystemCase:
    label: str
    projectile_z: int
    projectile_a: int
    target_z: int
    target_a: int
    energy_per_a: float
    impact_parameters: tuple[float, ...]
    n_steps: int = 100
    separation: float = 24.0
    dt: float = 1.0
    projectile_seed: int = 1001
    target_seed: int = 2001
    collision_seed: int = 3001


def make_collision_event(case: SystemCase, impact_parameter: float) -> ImQMDNucleus:
    projectile = initialize_and_relax(
        case.projectile_z,
        case.projectile_a,
        sigma_r=1.1,
        seed=case.projectile_seed + int(10 * impact_parameter),
    )
    target = initialize_and_relax(
        case.target_z,
        case.target_a,
        sigma_r=1.1,
        seed=case.target_seed + int(10 * impact_parameter),
    )
    p_beam = float(np.sqrt(2.0 * M_N_TEST * case.energy_per_a))

    packets: list[GaussianPacket] = []
    reference_positions: list[np.ndarray] = []
    reference_group_ids: list[int] = []

    for packet in projectile.packets:
        shift = np.array([-0.5 * case.separation, 0.5 * impact_parameter, 0.0])
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
        shift = np.array([0.5 * case.separation, -0.5 * impact_parameter, 0.0])
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
        case.projectile_z + case.target_z,
        (case.projectile_a - case.projectile_z) + (case.target_a - case.target_z),
        packets,
        edf=projectile.edf,
        reference_positions=np.asarray(reference_positions),
    )
    system.reference_group_ids = np.asarray(reference_group_ids, dtype=int)
    system._collision_rng = np.random.default_rng(case.collision_seed + int(10 * impact_parameter))
    system.collision_stats = {"attempted": 0, "blocked": 0, "accepted": 0}
    return system


def _fragment_record(fragment) -> dict[str, object]:
    return {
        "Z": int(fragment.Z),
        "A": int(fragment.A),
        "E_star": float(fragment.excitation_energy),
        "x": float(fragment.position[0]),
    }


def run_event(case: SystemCase, impact_parameter: float) -> dict[str, object]:
    system = make_collision_event(case, impact_parameter)
    start = time.perf_counter()
    old_handler = signal.signal(signal.SIGALRM, _timeout_handler)
    signal.setitimer(signal.ITIMER_REAL, EVENT_TIMEOUT_S)
    try:
        propagate(
            system,
            dt=case.dt,
            n_steps=case.n_steps,
            sample_every=max(1, case.n_steps // 4),
            with_collisions=True,
            collision_dt=1.0,
            remove_cm_drift=False,
            use_surface_term=True,
            use_static_stabilizer=False,
        )
    finally:
        signal.setitimer(signal.ITIMER_REAL, 0.0)
        signal.signal(signal.SIGALRM, old_handler)

    fragments = reaction_fragments(system)
    heavy = [fragment for fragment in fragments if fragment.A > 3]
    neighbor_cache = getattr(system, "_neighbor_list", {})
    return {
        "label": case.label,
        "impact_parameter": float(impact_parameter),
        "runtime_s": float(time.perf_counter() - start),
        "collision_stats": dict(getattr(system, "collision_stats", {"attempted": 0, "blocked": 0, "accepted": 0})),
        "neighbor_pairs": int(len(neighbor_cache.get("pairs", []))),
        "neighbor_rebuilds": int(neighbor_cache.get("rebuilds", 0)),
        "fragments": [_fragment_record(fragment) for fragment in fragments],
        "heavy_fragments": [_fragment_record(fragment) for fragment in heavy],
    }


def _print_result(result: dict[str, object]) -> None:
    heavy = result["heavy_fragments"]
    stats = result["collision_stats"]
    z_dist = Counter(fragment["Z"] for fragment in heavy)
    a_dist = Counter(fragment["A"] for fragment in heavy)
    print(f"\n{result['label']} b={result['impact_parameter']:.0f} fm")
    print(f"  runtime: {result['runtime_s']:.2f} s")
    print(
        "  collisions: "
        f"attempted={stats['attempted']}, blocked={stats['blocked']}, accepted={stats['accepted']}"
    )
    print(f"  neighbor list: pairs={result['neighbor_pairs']}, rebuilds={result['neighbor_rebuilds']}")
    print(f"  heavy fragments: {[(f['Z'], f['A']) for f in heavy[:6]]}")
    print(f"  Z distribution: {dict(sorted(z_dist.items()))}")
    print(f"  A distribution: {dict(sorted(a_dist.items()))}")


ASYMMETRIC_CASES = (
    SystemCase("40Ca+48Ca", 20, 40, 20, 48, 10.0, (3.0, 5.0), projectile_seed=4100, target_seed=4800, collision_seed=8900),
    SystemCase("16O+40Ca", 8, 16, 20, 40, 8.0, (2.0, 4.0), projectile_seed=1600, target_seed=4000, collision_seed=5600),
    SystemCase("86Kr+64Ni", 36, 86, 28, 64, 25.0, (3.0, 5.0), projectile_seed=8600, target_seed=6400, collision_seed=15000),
)

U_BENCHMARK = SystemCase("238U+238U", 92, 238, 92, 238, 7.0, (5.0,), n_steps=400, projectile_seed=23800, target_seed=23801, collision_seed=23802)


def run_asymmetric_transfer_suite() -> list[dict[str, object]]:
    results: list[dict[str, object]] = []
    for case in ASYMMETRIC_CASES:
        case_results = [run_event(case, b) for b in case.impact_parameters]
        results.extend(case_results)
        for result in case_results:
            _print_result(result)

        heavy = [fragment for result in case_results for fragment in result["heavy_fragments"]]
        if case.label == "40Ca+48Ca":
            assert any(fragment["A"] not in {40, 48} for fragment in heavy), (
                "Expected neutron transfer in 40Ca+48Ca: at least one heavy fragment with A not in {40,48}"
            )
        elif case.label == "16O+40Ca":
            assert any(fragment["Z"] not in {8, 20} for fragment in heavy), (
                "Expected charge transfer in 16O+40Ca: at least one heavy fragment with Z outside {8,20}"
            )
        elif case.label == "86Kr+64Ni":
            assert any(fragment["Z"] not in {36, 28} for fragment in heavy), (
                "Expected charge spread in 86Kr+64Ni: at least one heavy fragment with Z outside {36,28}"
            )
    return results


def run_u_benchmark() -> dict[str, object]:
    result = run_event(U_BENCHMARK, 5.0)
    _print_result(result)
    assert result["runtime_s"] < 300.0, "Expected 238U+238U b=5 fm runtime below 5 min"
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--include-u-benchmark", action="store_true", help="also run 238U+238U b=5 fm timing check")
    args = parser.parse_args()

    run_asymmetric_transfer_suite()
    if args.include_u_benchmark:
        run_u_benchmark()
    print("\nAll transfer verification assertions passed.")


if __name__ == "__main__":
    main()
