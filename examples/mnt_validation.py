#!/usr/bin/env python
"""Validate near-barrier U+U MNT events with ImQMD fragments.

Default run:
    python examples/mnt_validation.py

Use --full to run all requested impact parameters.
"""

from __future__ import annotations

import argparse
import os
import pathlib
import sys
import time

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

M_N = 938.9
DEFAULT_B_VALUES = (5.0, 7.0, 10.0)
FULL_B_VALUES = (0.0, 2.0, 4.0, 5.0, 6.0, 7.0, 8.0, 10.0, 12.0)


def make_uu_event(impact_parameter: float, energy_per_a: float, seed: int) -> ImQMDNucleus:
    """Create a fixed-target U+U event at projectile lab energy per nucleon."""

    projectile = initialize_and_relax(92, 238, sigma_r=1.1, seed=seed)
    target = initialize_and_relax(92, 238, sigma_r=1.1, seed=seed + 1)
    p_beam = np.sqrt(2.0 * M_N * energy_per_a)
    separation = 28.0
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
        packets.append(GaussianPacket(packet.r_i + shift, packet.p_i.copy(), packet.sigma_r, packet.is_proton))
        reference_positions.append(packet.r_i + shift)
        reference_group_ids.append(1)

    nucleus = ImQMDNucleus(184, 292, packets, edf=projectile.edf, reference_positions=np.asarray(reference_positions))
    nucleus.reference_group_ids = np.asarray(reference_group_ids, dtype=int)
    nucleus._collision_rng = np.random.default_rng(seed + 100_000)
    return nucleus


def fragment_observables(fragment, b_value: float, event_index: int, stage: str) -> dict[str, float | int | str]:
    p = np.asarray(fragment.momentum, dtype=float)
    p_norm = float(np.linalg.norm(p))
    theta_lab = 0.0 if p_norm <= 1.0e-12 else float(np.degrees(np.arccos(np.clip(p[0] / p_norm, -1.0, 1.0))))
    e_lab = float(np.dot(p, p) / (2.0 * M_N * max(fragment.A, 1)))
    return {
        "event": int(event_index),
        "b": float(b_value),
        "stage": stage,
        "Z": int(fragment.Z),
        "A": int(fragment.A),
        "theta_lab": theta_lab,
        "E_lab": e_lab,
        "E_star": float(fragment.excitation_energy),
        "px": float(p[0]),
        "py": float(p[1]),
        "pz": float(p[2]),
    }


def evaporate_fragment(fragment, rng: np.random.Generator) -> tuple[int, int]:
    if fragment.A <= 100 or fragment.excitation_energy <= 5.0:
        return int(fragment.Z), int(fragment.A)
    channels = weisskopf_evaporation(fragment.Z, fragment.A, fragment.excitation_energy, n_max=12)
    probabilities = np.asarray([channel[2] for channel in channels], dtype=float)
    probabilities /= np.sum(probabilities)
    choice = int(rng.choice(len(channels), p=probabilities))
    z_final, a_final, _ = channels[choice]
    return int(z_final), int(a_final)


def save_results(path: pathlib.Path, primary_rows: list[dict], final_rows: list[dict], summaries: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    np.savez(
        path,
        primary=np.asarray(primary_rows, dtype=object),
        final=np.asarray(final_rows, dtype=object),
        summaries=np.asarray(summaries, dtype=object),
    )


def run_event(args: argparse.Namespace, b_value: float, event_index: int) -> tuple[list[dict], list[dict], dict]:
    seed = int(args.seed + 1000 * event_index + round(10.0 * b_value))
    rng = np.random.default_rng(seed + 200_000)
    nucleus = make_uu_event(b_value, args.energy_per_a, seed)

    started = time.time()
    print(f"b={b_value:g} fm: propagating {args.time:g} fm/c, dt={args.dt:g}, seed={seed}", flush=True)
    propagate(
        nucleus,
        dt=args.dt,
        n_steps=int(round(args.time / args.dt)),
        sample_every=max(1, int(round(100.0 / args.dt))),
        with_collisions=True,
        collision_dt=args.collision_dt,
        remove_cm_drift=False,
        use_surface_term=True,
        use_static_stabilizer=False,
    )

    fragments = reaction_fragments(nucleus)
    primary_rows = [fragment_observables(fragment, b_value, event_index, "primary") for fragment in fragments]
    final_rows: list[dict] = []
    for fragment in fragments:
        row = fragment_observables(fragment, b_value, event_index, "evaporated")
        z_final, a_final = evaporate_fragment(fragment, rng)
        row["Z"] = z_final
        row["A"] = a_final
        final_rows.append(row)

    heavy = [fragment for fragment in fragments if fragment.A > 100]
    summary = {
        "event": int(event_index),
        "b": float(b_value),
        "seed": int(seed),
        "fragments": int(len(fragments)),
        "heavy_fragments": int(len(heavy)),
        "mean_E_star_heavy": float(np.mean([fragment.excitation_energy for fragment in heavy])) if heavy else 0.0,
        "accepted_collisions": int(getattr(nucleus, "collision_stats", {}).get("accepted", 0)),
        "attempted_collisions": int(getattr(nucleus, "collision_stats", {}).get("attempted", 0)),
        "elapsed_s": float(time.time() - started),
    }
    return primary_rows, final_rows, summary


def print_summary(final_rows: list[dict], summaries: list[dict]) -> None:
    z_values = np.asarray([row["Z"] for row in final_rows], dtype=int)
    heavy_e = [summary["mean_E_star_heavy"] for summary in summaries if summary["heavy_fragments"]]
    if z_values.size:
        unique_z, counts = np.unique(z_values, return_counts=True)
        z_peak = int(unique_z[int(np.argmax(counts))])
    else:
        z_peak = -1
    print("\nMNT validation summary")
    print(f"  events: {len(summaries)}")
    print(f"  total final fragments: {len(final_rows)}")
    print(f"  Z distribution peak: {z_peak}")
    print(f"  mean heavy-fragment E*: {np.mean(heavy_e) if heavy_e else 0.0:.2f} MeV")
    for summary in summaries:
        print(
            "  b={b:g} fm: fragments={fragments}, heavy={heavy_fragments}, "
            "accepted collisions={accepted_collisions}, mean heavy E*={mean_E_star_heavy:.2f} MeV, "
            "elapsed={elapsed_s:.1f}s".format(**summary)
        )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--full", action="store_true", help="run all requested b-values instead of the 3-event smoke set")
    parser.add_argument("--b-values", type=float, nargs="+", default=None, help="override impact parameters in fm")
    parser.add_argument("--energy-per-a", type=float, default=7.0, help="projectile lab energy in MeV/A")
    parser.add_argument("--time", type=float, default=2000.0, help="propagation time in fm/c")
    parser.add_argument("--dt", type=float, default=1.0, help="time step in fm/c")
    parser.add_argument("--collision-dt", type=float, default=1.0, help="collision search interval in fm/c")
    parser.add_argument("--seed", type=int, default=7000, help="base random seed")
    parser.add_argument("--output", default="output/mnt_validation.npz", help="output npz path")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    b_values = args.b_values if args.b_values is not None else (FULL_B_VALUES if args.full else DEFAULT_B_VALUES)
    output_path = pathlib.Path(args.output)
    primary_rows: list[dict] = []
    final_rows: list[dict] = []
    summaries: list[dict] = []

    for event_index, b_value in enumerate(b_values):
        event_primary, event_final, event_summary = run_event(args, float(b_value), event_index)
        primary_rows.extend(event_primary)
        final_rows.extend(event_final)
        summaries.append(event_summary)
        save_results(output_path, primary_rows, final_rows, summaries)
        print(f"b={b_value:g} fm: saved intermediate results to {output_path}", flush=True)

    print_summary(final_rows, summaries)
    print(f"Saved {output_path}")


if __name__ == "__main__":
    main()
