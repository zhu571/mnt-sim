#!/usr/bin/env python
"""Fast end-to-end ImQMD MNT pipeline verification.

Run: python examples/verify_mnt_pipeline.py
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
    initialize_and_relax,
    propagate,
    reaction_fragments,
    weisskopf_evaporation,
)
from verify_transfer import run_asymmetric_transfer_suite  # noqa: E402


M_N_TEST = 938.9
EVENT_TIMEOUT_S = 180


class EventTimeout(RuntimeError):
    """Raised when one verification event exceeds the interactive time budget."""


def _timeout_handler(signum, frame) -> None:
    del signum, frame
    raise EventTimeout(f"event exceeded {EVENT_TIMEOUT_S} s")


def make_collision_event(
    Z: int,
    A: int,
    energy_per_a: float,
    separation: float,
    impact_parameter: float,
    projectile_seed: int,
    target_seed: int,
    collision_seed: int,
) -> ImQMDNucleus:
    projectile = initialize_and_relax(Z, A, sigma_r=1.1, seed=projectile_seed)
    target = initialize_and_relax(Z, A, sigma_r=1.1, seed=target_seed)
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
        2 * Z,
        2 * (A - Z),
        packets,
        edf=projectile.edf,
        reference_positions=np.asarray(reference_positions),
    )
    system.reference_group_ids = np.asarray(reference_group_ids, dtype=int)
    system._collision_rng = np.random.default_rng(collision_seed)
    system.collision_stats = {"attempted": 0, "blocked": 0, "accepted": 0}
    return system


def _fragment_record(fragment) -> dict[str, object]:
    return {
        "Z": int(fragment.Z),
        "A": int(fragment.A),
        "E_star": float(fragment.excitation_energy),
        "position": [float(x) for x in fragment.position],
        "momentum": [float(x) for x in fragment.momentum],
    }


def _evaporation_records(heavy_fragments) -> list[dict[str, object]]:
    records: list[dict[str, object]] = []
    for fragment in heavy_fragments:
        if fragment.excitation_energy <= 5.0:
            continue
        channels = weisskopf_evaporation(fragment.Z, fragment.A, fragment.excitation_energy, n_max=8)
        dominant = max(channels, key=lambda item: item[2])
        mean_a = float(sum(a_final * p for _, a_final, p in channels))
        records.append(
            {
                "before_Z": int(fragment.Z),
                "before_A": int(fragment.A),
                "E_star": float(fragment.excitation_energy),
                "dominant_after_Z": int(dominant[0]),
                "dominant_after_A": int(dominant[1]),
                "dominant_probability": float(dominant[2]),
                "mean_after_A": mean_a,
                "has_lower_A_channel": bool(any(a_final < fragment.A and p > 0.0 for _, a_final, p in channels)),
            }
        )
    return records


def _print_event_summary(label: str, result: dict[str, object], reference_z: int) -> None:
    fragments = result["fragments"]
    heavy = result["heavy_fragments"]
    largest = heavy[:2]
    stats = result["collision_stats"]
    transfer = [fragment for fragment in fragments if fragment["Z"] != reference_z]

    print(f"\n{label}")
    print(f"  runtime: {result['runtime_s']:.2f} s")
    print(f"  fragments: total={len(fragments)}, heavy(A>10)={len(heavy)}")
    for i, fragment in enumerate(largest, start=1):
        print(
            f"  largest {i}: Z={fragment['Z']}, A={fragment['A']}, "
            f"E*={fragment['E_star']:.2f} MeV"
        )
    for i in range(len(largest) + 1, 3):
        print(f"  largest {i}: none")
    print(
        "  collisions: "
        f"attempted={stats['attempted']}, blocked={stats['blocked']}, accepted={stats['accepted']}"
    )
    if transfer:
        transfer_text = ", ".join(f"Z={fragment['Z']} A={fragment['A']}" for fragment in transfer)
    else:
        transfer_text = "none"
    print(f"  fragments with Z != {reference_z}: {transfer_text}")

    evaporation = result["evaporation"]
    if evaporation:
        print("  Weisskopf evaporation for heavy E*>5 MeV fragments:")
        for item in evaporation:
            print(
                "    "
                f"before Z={item['before_Z']} A={item['before_A']} E*={item['E_star']:.2f} -> "
                f"dominant Z={item['dominant_after_Z']} A={item['dominant_after_A']} "
                f"p={item['dominant_probability']:.3f}, mean A={item['mean_after_A']:.2f}"
            )
    else:
        print("  Weisskopf evaporation: no heavy fragments with E*>5 MeV")


def run_event(
    label: str,
    Z: int,
    A: int,
    energy_per_a: float,
    separation: float,
    impact_parameter: float,
    dt: float,
    n_steps: int,
    projectile_seed: int,
    target_seed: int,
    collision_seed: int,
) -> dict[str, object]:
    system = make_collision_event(
        Z=Z,
        A=A,
        energy_per_a=energy_per_a,
        separation=separation,
        impact_parameter=impact_parameter,
        projectile_seed=projectile_seed,
        target_seed=target_seed,
        collision_seed=collision_seed,
    )
    start = time.perf_counter()
    old_handler = signal.signal(signal.SIGALRM, _timeout_handler)
    signal.setitimer(signal.ITIMER_REAL, EVENT_TIMEOUT_S)
    try:
        propagate(
            system,
            dt=dt,
            n_steps=n_steps,
            sample_every=max(1, n_steps // 4),
            with_collisions=True,
            collision_dt=1.0,
            remove_cm_drift=False,
            use_surface_term=True,
            use_static_stabilizer=False,
        )
    finally:
        signal.setitimer(signal.ITIMER_REAL, 0.0)
        signal.signal(signal.SIGALRM, old_handler)
    runtime_s = time.perf_counter() - start

    fragments = reaction_fragments(system)
    heavy = [fragment for fragment in fragments if fragment.A > 10]
    return {
        "label": label,
        "Z_reference": int(Z),
        "A_reference": int(A),
        "impact_parameter": float(impact_parameter),
        "energy_per_a": float(energy_per_a),
        "runtime_s": float(runtime_s),
        "collision_stats": dict(getattr(system, "collision_stats", {"attempted": 0, "blocked": 0, "accepted": 0})),
        "fragments": [_fragment_record(fragment) for fragment in fragments],
        "heavy_fragments": [_fragment_record(fragment) for fragment in heavy],
        "evaporation": _evaporation_records(heavy),
    }


def main() -> None:
    results: list[dict[str, object]] = []
    print("40Ca+40Ca verification: E/A=10 MeV, separation=20 fm, p_beam=%.2f MeV/c" % np.sqrt(2.0 * M_N_TEST * 10.0))
    for i, b in enumerate([2.0, 4.0, 6.0], start=1):
        result = run_event(
            label=f"40Ca+40Ca b={b:.0f} fm",
            Z=20,
            A=40,
            energy_per_a=10.0,
            separation=20.0,
            impact_parameter=b,
            dt=1.0,
            n_steps=400,
            projectile_seed={1: 6000, 2: 7010, 3: 4020}[i],
            target_seed={1: 6001, 2: 7011, 3: 4021}[i],
            collision_seed={1: 6002, 2: 7012, 3: 4022}[i],
        )
        results.append(result)
        _print_event_summary(result["label"], result, reference_z=20)
        if result["runtime_s"] > 180.0:
            raise TimeoutError(f"{result['label']} took {result['runtime_s']:.1f} s; stopping per task limit")

    print("\n16O+16O quick test: E/A=8 MeV, b=3 fm, separation=20 fm")
    oxygen = run_event(
        label="16O+16O b=3 fm",
        Z=8,
        A=16,
        energy_per_a=8.0,
        separation=20.0,
        impact_parameter=3.0,
        dt=1.0,
        n_steps=300,
        projectile_seed=1601,
        target_seed=1602,
        collision_seed=1603,
    )
    results.append(oxygen)
    _print_event_summary(oxygen["label"], oxygen, reference_z=8)

    oxygen_heavy = oxygen["heavy_fragments"]
    if len(oxygen_heavy) >= 2:
        x_positions = [fragment["position"][0] for fragment in oxygen_heavy[:2]]
        print(f"  two largest O fragments x positions: {x_positions[0]:.2f}, {x_positions[1]:.2f} fm")
    oxygen_transfer = any(fragment["Z"] != 8 for fragment in oxygen_heavy)
    print(f"  O transfer observed at this low energy: {oxygen_transfer}")

    print("\nAsymmetric transfer verification")
    asymmetric_results = run_asymmetric_transfer_suite()
    results.extend(asymmetric_results)

    ca_results = [result for result in results if result.get("Z_reference") == 20]
    assert all(
        result["collision_stats"]["accepted"] > 0
        for result in ca_results
        if result["impact_parameter"] <= 4.0
    ), "Expected accepted NN collisions for symmetric Ca+Ca b <= 4 fm"
    assert any(
        evap["has_lower_A_channel"]
        for result in ca_results
        for evap in result["evaporation"]
    ), "Expected Weisskopf evaporation to provide lower-A channels for some excited heavy fragment"
    assert len(oxygen_heavy) >= 2, "Expected at least two heavy O+O fragments after collision"
    assert oxygen_heavy[0]["position"][0] * oxygen_heavy[1]["position"][0] < 0.0, (
        "Expected the two largest O+O fragments to separate to opposite x sides"
    )

    output_dir = Path("output")
    output_dir.mkdir(exist_ok=True)
    output_path = output_dir / "ca40_verify.npz"
    np.savez(
        output_path,
        summary_json=np.array(json.dumps(results, indent=2)),
        labels=np.array([result["label"] for result in results]),
        impact_parameters=np.array([result["impact_parameter"] for result in results], dtype=float),
        accepted_collisions=np.array([result["collision_stats"]["accepted"] for result in results], dtype=int),
        attempted_collisions=np.array([result["collision_stats"]["attempted"] for result in results], dtype=int),
        blocked_collisions=np.array([result["collision_stats"]["blocked"] for result in results], dtype=int),
    )
    print(f"\nSaved summary to {output_path}")
    print("All in-code assertions passed.")


if __name__ == "__main__":
    main()
