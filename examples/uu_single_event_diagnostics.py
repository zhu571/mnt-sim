"""Run one U+U event and dump energy, transfer, and fragment diagnostics."""

from __future__ import annotations

import argparse
import pathlib
import sys

import numpy as np

ROOT = pathlib.Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from mnt_sim.imqmd import (
    SkyrmeEDF,
    grid_energy_diagnostics,
    identify_fragments,
    make_collision_event,
    propagate,
)
from mnt_sim.imqmd.skyrme import PARAMETER_SETS


def _fragment_summary(label: str, fragments: list) -> None:
    heavy = sorted((fragment for fragment in fragments if fragment.A >= 20), key=lambda fragment: fragment.A, reverse=True)
    print(f"\n[{label}] {len(fragments)} fragments")
    for fragment in heavy[:10]:
        print(
            f"  Z={fragment.Z:3d} A={fragment.A:3d} "
            f"E*={fragment.excitation_energy:8.1f} MeV "
            f"x={fragment.position[0]:7.2f} fm"
        )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--parameter-set", default="IQ3A", choices=sorted(PARAMETER_SETS), help="EDF parameter set")
    parser.add_argument("--impact-parameter", type=float, default=4.0)
    parser.add_argument("--energy-per-a", type=float, default=7.0)
    parser.add_argument("--time", type=float, default=600.0, help="Propagation time in fm/c")
    parser.add_argument("--dt", type=float, default=1.0)
    parser.add_argument("--collision-dt", type=float, default=1.0)
    parser.add_argument("--relax-time", type=float, default=200.0)
    parser.add_argument("--sigma-r", type=float, default=1.1)
    args = parser.parse_args()

    edf = SkyrmeEDF(PARAMETER_SETS[args.parameter_set])
    diag = grid_energy_diagnostics(92, 238, sigma_r=args.sigma_r, seed=23801, edf=edf)
    print("[U-238 grid diagnostics]")
    for key in (
        "sigma_r",
        "target_total",
        "raw_total",
        "raw_skyrme_bulk",
        "raw_surface",
        "raw_surface_symmetry",
        "raw_symmetry",
        "raw_coulomb_direct",
        "raw_coulomb_exchange",
        "nuclear_scale",
        "scaled_total",
        "energy_offset",
    ):
        print(f"  {key:22s} {diag[key]:12.3f}")

    system = make_collision_event(
        92,
        238,
        92,
        238,
        args.energy_per_a,
        args.impact_parameter,
        projectile_seed=23801,
        target_seed=23802,
        collision_seed=23804,
        sigma_r=args.sigma_r,
        relax_time=args.relax_time,
        edf=edf,
    )
    propagate(
        system,
        dt=args.dt,
        n_steps=max(1, int(round(args.time / args.dt))),
        sample_every=max(1, int(round(100.0 / args.dt))),
        with_collisions=True,
        collision_dt=args.collision_dt,
        remove_cm_drift=False,
        use_surface_term=True,
        use_static_stabilizer=False,
        use_grid_edf=True,
    )
    print("\n[Collision stats]")
    for key, value in system.collision_stats.items():
        print(f"  {key:12s} {value}")

    left = system.positions[:, 0] < 0.0
    projectile_mask = system.reference_group_ids == 0
    projectile_to_right = int(np.sum(projectile_mask & ~left))
    target_to_left = int(np.sum(~projectile_mask & left))
    print("\n[Transfer by side crossing]")
    print(f"  projectile nucleons on right: {projectile_to_right}")
    print(f"  target nucleons on left:     {target_to_left}")
    print(f"  projectile protons on right: {int(np.sum(projectile_mask & system.is_proton & ~left))}")
    print(f"  target protons on left:      {int(np.sum(~projectile_mask & system.is_proton & left))}")

    default_fragments = identify_fragments(system, method="iso-mst", p_cut=250.0)
    tight_fragments = identify_fragments(
        system,
        method="iso-mst",
        p_cut=250.0,
        iso_r_cut_pp=2.5,
        iso_r_cut_nn=4.0,
        iso_r_cut_np=3.5,
    )
    _fragment_summary("iso-mst default", default_fragments)
    _fragment_summary("iso-mst tight", tight_fragments)


if __name__ == "__main__":
    main()
