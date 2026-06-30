#!/usr/bin/env python
"""Run the requested U+U benchmark and Kr+Ni impact-parameter scan.

Run:
    python examples/run_uu_krni_benchmark.py
"""

from __future__ import annotations

import json
import os
import signal
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from mnt_sim.imqmd import (  # noqa: E402
    GaussianPacket,
    ImQMDNucleus,
    attempt_nn_collision,
    fermi_constraint_check,
    initialize_and_relax,
    reaction_fragments,
    weisskopf_evaporation,
)
from mnt_sim.imqmd.propagator import _rk4_step  # noqa: E402


M_N_TEST = 938.9
UU_TIMEOUT_S = 8 * 60


class EventTimeout(RuntimeError):
    """Raised when an event exceeds its wall-clock budget."""


def _timeout_handler(signum, frame) -> None:
    del signum, frame
    raise EventTimeout("event exceeded wall-clock budget")


def make_collision_event(
    projectile_z: int,
    projectile_a: int,
    target_z: int,
    target_a: int,
    energy_per_a: float,
    separation: float,
    impact_parameter: float,
    projectile_seed: int,
    target_seed: int,
    collision_seed: int,
) -> ImQMDNucleus:
    projectile = initialize_and_relax(projectile_z, projectile_a, sigma_r=1.1, seed=projectile_seed)
    target = initialize_and_relax(target_z, target_a, sigma_r=1.1, seed=target_seed)
    p_beam = float(np.sqrt(2.0 * M_N_TEST * energy_per_a))

    packets: list[GaussianPacket] = []
    reference_positions: list[np.ndarray] = []
    reference_group_ids: list[int] = []

    for packet in projectile.packets:
        shift = np.array([-0.5 * separation, 0.5 * impact_parameter, 0.0])
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
        shift = np.array([0.5 * separation, -0.5 * impact_parameter, 0.0])
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
        edf=projectile.edf,
        reference_positions=np.asarray(reference_positions),
    )
    system.reference_group_ids = np.asarray(reference_group_ids, dtype=int)
    system._collision_rng = np.random.default_rng(collision_seed)
    system.collision_stats = {"attempted": 0, "blocked": 0, "accepted": 0}
    return system


def run_timed_event(
    system: ImQMDNucleus,
    dt: float,
    n_steps: int,
    collision_dt: float,
    timeout_s: int | None = None,
) -> dict[str, object]:
    collision_interval = max(1, int(round(float(collision_dt) / float(dt))))
    fermi_interval = max(1, int(round(5.0 / float(dt))))
    neighbor_pairs_by_step: list[int] = []
    collision_step_stats: list[dict[str, int]] = []
    state = {"step": 0, "phase": "initializing"}
    system.use_surface_term = True
    system.use_static_stabilizer = False

    old_handler = None
    if timeout_s is not None:
        old_handler = signal.signal(signal.SIGALRM, _timeout_handler)
        signal.setitimer(signal.ITIMER_REAL, float(timeout_s))

    start = time.perf_counter()
    try:
        for step in range(1, int(n_steps) + 1):
            state["step"] = step
            state["phase"] = "rk4"
            _rk4_step(
                system,
                float(dt),
                remove_cm_drift=False,
                use_surface_term=True,
                use_static_stabilizer=False,
            )
            if step % collision_interval == 0:
                state["phase"] = "collisions"
                step_stats = attempt_nn_collision(system, dt=float(dt) * collision_interval)
                collision_step_stats.append({key: int(value) for key, value in step_stats.items()})
                neighbor_cache = getattr(system, "_neighbor_list", {})
                neighbor_pairs_by_step.append(int(len(neighbor_cache.get("pairs", []))))
            if step % fermi_interval == 0:
                state["phase"] = "fermi_constraint"
                fermi_constraint_check(system)
        timed_out = False
        timeout_where = None
    except EventTimeout:
        timed_out = True
        timeout_where = f"step={state['step']} phase={state['phase']}"
    finally:
        if timeout_s is not None:
            signal.setitimer(signal.ITIMER_REAL, 0.0)
            signal.signal(signal.SIGALRM, old_handler)

    runtime_s = time.perf_counter() - start
    neighbor_cache = getattr(system, "_neighbor_list", {})
    return {
        "runtime_s": float(runtime_s),
        "completed_steps": int(state["step"]),
        "timed_out": bool(timed_out),
        "timeout_where": timeout_where,
        "collision_stats": dict(getattr(system, "collision_stats", {"attempted": 0, "blocked": 0, "accepted": 0})),
        "neighbor_pairs_by_step": neighbor_pairs_by_step,
        "neighbor_rebuilds": int(neighbor_cache.get("rebuilds", 0)),
        "collision_step_stats": collision_step_stats,
    }


def fragment_record(fragment) -> dict[str, object]:
    return {
        "Z": int(fragment.Z),
        "A": int(fragment.A),
        "E_star": float(fragment.excitation_energy),
        "position": [float(x) for x in fragment.position],
        "momentum": [float(x) for x in fragment.momentum],
    }


def evaporation_record(fragment) -> dict[str, object] | None:
    if fragment.excitation_energy <= 5.0:
        return None
    channels = weisskopf_evaporation(fragment.Z, fragment.A, fragment.excitation_energy, n_max=10)
    lower_channels = [item for item in channels if item[1] < fragment.A and item[2] > 0.0]
    if lower_channels:
        final = max(lower_channels, key=lambda item: item[2])
    else:
        final = max(channels, key=lambda item: item[2])
    return {
        "before_Z": int(fragment.Z),
        "before_A": int(fragment.A),
        "E_star": float(fragment.excitation_energy),
        "final_Z": int(final[0]),
        "final_A": int(final[1]),
        "final_probability": float(final[2]),
        "dominant_Z": int(max(channels, key=lambda item: item[2])[0]),
        "dominant_A": int(max(channels, key=lambda item: item[2])[1]),
        "has_lower_A_channel": bool(lower_channels),
    }


def save_json_npz(path: Path, result: dict[str, object]) -> None:
    path.parent.mkdir(exist_ok=True)
    np.savez(path, summary_json=np.array(json.dumps(result, indent=2)), result=np.array(result, dtype=object))


def print_neighbor_summary(prefix: str, pairs: list[int]) -> None:
    if not pairs:
        print(f"{prefix} neighbor_list pairs per step: none")
        return
    arr = np.asarray(pairs, dtype=int)
    print(
        f"{prefix} neighbor_list pairs per step: "
        f"mean={arr.mean():.1f}, min={arr.min()}, max={arr.max()}, last={arr[-1]}, n={len(arr)}"
    )


def run_uu() -> dict[str, object]:
    print("U+U: E/A=7 MeV, b=5 fm, 1200 fm/c, dt=1.0")
    system = make_collision_event(
        projectile_z=92,
        projectile_a=238,
        target_z=92,
        target_a=238,
        energy_per_a=7.0,
        separation=28.0,
        impact_parameter=5.0,
        projectile_seed=23850,
        target_seed=23851,
        collision_seed=23852,
    )
    run = run_timed_event(system, dt=1.0, n_steps=1200, collision_dt=1.0, timeout_s=UU_TIMEOUT_S)
    fragments = reaction_fragments(system)
    largest = sorted(fragments, key=lambda item: item.A, reverse=True)[:3]
    result = {
        "system": "238U+238U",
        "energy_per_a": 7.0,
        "impact_parameter": 5.0,
        "n_steps": 1200,
        "dt": 1.0,
        "separation": 28.0,
        "remove_cm_drift": False,
        **run,
        "fragments": [fragment_record(fragment) for fragment in fragments],
        "largest_3_Z": [int(fragment.Z) for fragment in largest],
    }
    save_json_npz(Path("output/uu_benchmark.npz"), result)

    stats = result["collision_stats"]
    print(f"U+U runtime: {result['runtime_s']:.2f} s ({result['runtime_s'] / 60.0:.2f} min)")
    print_neighbor_summary("U+U", result["neighbor_pairs_by_step"])
    print(
        "U+U collisions: "
        f"attempted={stats['attempted']}, blocked={stats['blocked']}, accepted={stats['accepted']}"
    )
    print(f"U+U fragments found: {len(fragments)}")
    print(f"U+U Z of 3 largest fragments: {result['largest_3_Z']}")
    if result["timed_out"]:
        print(f"U+U killed after timeout at {result['timeout_where']}")
    else:
        assert not result["timed_out"], "U+U should complete without timeout"
    print("Saved output/uu_benchmark.npz")
    return result


def run_krni_event(b: float) -> dict[str, object]:
    system = make_collision_event(
        projectile_z=36,
        projectile_a=86,
        target_z=28,
        target_a=64,
        energy_per_a=25.0,
        separation=24.0,
        impact_parameter=b,
        projectile_seed=8600 + int(10 * b),
        target_seed=6400 + int(10 * b),
        collision_seed=15000 + int(10 * b),
    )
    run = run_timed_event(system, dt=1.0, n_steps=600, collision_dt=1.0, timeout_s=None)
    fragments = reaction_fragments(system)
    heavy = [fragment for fragment in fragments if fragment.A > 10]
    evaporated = [evaporation_record(fragment) for fragment in heavy]
    evaporated = [item for item in evaporated if item is not None]
    return {
        "system": "86Kr+64Ni",
        "energy_per_a": 25.0,
        "impact_parameter": float(b),
        "n_steps": 600,
        "dt": 1.0,
        "remove_cm_drift": False,
        **run,
        "fragments": [fragment_record(fragment) for fragment in fragments],
        "heavy_fragments": [fragment_record(fragment) for fragment in heavy],
        "evaporation": evaporated,
    }


def run_krni_scan() -> list[dict[str, object]]:
    print("\nKr+Ni scan: E/A=25 MeV, b=[2, 3, 4, 5, 6, 7] fm")
    results = []
    for b in [2.0, 3.0, 4.0, 5.0, 6.0, 7.0]:
        result = run_krni_event(b)
        results.append(result)
        stats = result["collision_stats"]
        print(
            f"b={b:.0f} runtime={result['runtime_s']:.2f}s "
            f"collisions att/block/acc={stats['attempted']}/{stats['blocked']}/{stats['accepted']}"
        )
        print_neighbor_summary(f"b={b:.0f}", result["neighbor_pairs_by_step"])

    save_json_npz(Path("output/krni_scan.npz"), {"results": results})
    print("\nKr+Ni summary")
    print(f"{'b':>3} {'n_fragments':>11} {'Z_largest_2':>14} {'A_largest_2':>14} {'E*_largest_2':>20}")
    for result in results:
        heavy = result["heavy_fragments"]
        z2 = [fragment["Z"] for fragment in heavy[:2]]
        a2 = [fragment["A"] for fragment in heavy[:2]]
        e2 = [round(fragment["E_star"], 2) for fragment in heavy[:2]]
        print(f"{result['impact_parameter']:3.0f} {len(heavy):11d} {str(z2):>14} {str(a2):>14} {str(e2):>20}")
    print("Saved output/krni_scan.npz")

    b3 = next(result for result in results if result["impact_parameter"] == 3.0)
    assert any(
        fragment["Z"] not in {36, 28} for fragment in b3["heavy_fragments"]
    ), "Kr+Ni at b=3 should show transfer"
    for result in results:
        if result["impact_parameter"] >= 6.0:
            z_values = {fragment["Z"] for fragment in result["heavy_fragments"][:2]}
            assert {36, 28}.issubset(z_values), "Kr+Ni at b>=6 should show elastic original-Z fragments"
    assert any(
        evap["final_A"] < evap["before_A"]
        for result in results
        for evap in result["evaporation"]
    ), "Evaporation should reduce A for excited fragments"
    print("Kr+Ni in-code assertions passed.")
    return results


def main() -> None:
    uu_result = run_uu()
    if uu_result["timed_out"]:
        print("Continuing with Kr+Ni after U+U timeout.")
    run_krni_scan()


if __name__ == "__main__":
    main()
