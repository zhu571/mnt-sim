#!/usr/bin/env python
"""End-to-end ImQMD event with fragment recognition.

Run: python examples/reaction_event.py
"""

from __future__ import annotations

import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from mnt_sim.imqmd import (
    GaussianPacket,
    ImQMDNucleus,
    initialize_and_relax,
    isospin_mst,
    minimum_spanning_tree,
    propagate,
    reaction_fragments,
)


def make_uu_event(impact_parameter: float = 7.0, energy_per_a: float = 7.0) -> ImQMDNucleus:
    projectile = initialize_and_relax(92, 238, sigma_r=1.1, seed=2381)
    target = initialize_and_relax(92, 238, sigma_r=1.1, seed=2382)
    p_beam = np.sqrt(2.0 * 938.9 * energy_per_a)
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

    nucleus = ImQMDNucleus(184, 292, packets, edf=projectile.edf, reference_positions=np.asarray(reference_positions))
    nucleus.reference_group_ids = np.asarray(reference_group_ids, dtype=int)
    return nucleus


def _print_fragments(label: str, fragments) -> None:
    fragments = sorted(fragments, key=lambda frag: frag.A, reverse=True)
    print(f"{label}: recognized {len(fragments)} fragments")
    print(f"{'#':>3} {'Z':>4} {'A':>4} {'E* (MeV)':>10} {'x':>9} {'y':>9} {'z':>9} {'px':>10} {'py':>10} {'pz':>10}")
    for i, fragment in enumerate(fragments[:5], start=1):
        x, y, z = fragment.position
        px, py, pz = fragment.momentum
        print(
            f"{i:3d} {fragment.Z:4d} {fragment.A:4d} {fragment.excitation_energy:10.2f} "
            f"{x:9.2f} {y:9.2f} {z:9.2f} {px:10.2f} {py:10.2f} {pz:10.2f}"
        )


def main() -> None:
    nucleus = make_uu_event(impact_parameter=7.0, energy_per_a=7.0)
    print("U+U event: b=7 fm, E/A=7 MeV, A=476")
    print("Propagating 2000 fm/c with dt=2 fm/c...")
    propagate(
        nucleus,
        dt=2.0,
        n_steps=1000,
        sample_every=100,
        with_collisions=False,
        remove_cm_drift=False,
        use_surface_term=True,
        use_static_stabilizer=False,
    )

    standard = reaction_fragments(nucleus)
    isotope_aware = isospin_mst(nucleus, r_cut_pp=3.5, r_cut_nn=3.5, r_cut_np=4.0, p_cut=None)
    tight_standard = minimum_spanning_tree(nucleus, r_cut=3.0, p_cut=None)
    _print_fragments("Reaction MST", standard)
    _print_fragments("Isospin MST", isotope_aware)
    print(f"Tight r_cut=3.0 MST count: {len(tight_standard)}")


if __name__ == "__main__":
    main()
