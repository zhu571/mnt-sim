"""Inspect GridEDF raw binding and the energy_offset calibration for benchmark nuclei."""

from __future__ import annotations

import argparse
import pathlib
import sys
from dataclasses import dataclass

import numpy as np

ROOT = pathlib.Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from mnt_sim.imqmd import SkyrmeEDF, initialize_nucleus, propagate
from mnt_sim.imqmd.grid_edf import GridEDF
from mnt_sim.imqmd.skyrme import PARAMETER_SETS


@dataclass(frozen=True)
class Benchmark:
    label: str
    z: int
    a: int
    binding_mev: float
    rms_fm: float
    binding_tolerance_mev: float


BENCHMARKS = (
    Benchmark("40Ca", 20, 40, 342.0, 3.48, 20.0),
    Benchmark("208Pb", 82, 208, 1636.0, 5.50, 100.0),
)


def iq_width(a: int, parameter_set: str) -> float:
    params = PARAMETER_SETS[parameter_set.upper()]
    return float(params.sigma0 + params.sigma1 * a ** (1.0 / 3.0))


def mass_rms(nucleus) -> float:
    rp, rn = nucleus.rms_radius()
    return float(np.sqrt((nucleus.Z * rp * rp + nucleus.N * rn * rn) / nucleus.A))


def raw_grid_components(nucleus, grid_spacing: float, n_sigma: float) -> dict[str, float]:
    grid = GridEDF(
        nucleus.edf.parameters,
        nucleus.sigma_r,
        grid_spacing=grid_spacing,
        n_sigma=n_sigma,
    )
    components = grid.total_energy(nucleus.positions, nucleus.momenta, nucleus.is_proton)
    components["rms"] = mass_rms(nucleus)
    return components


def run_stability(nucleus, time_fm_c: float, dt: float) -> tuple[float, dict[str, float]]:
    trial = nucleus.copy()
    history = propagate(
        trial,
        dt=dt,
        n_steps=max(1, int(round(time_fm_c / dt))),
        sample_every=max(1, int(round(time_fm_c / dt))),
        with_collisions=False,
        remove_cm_drift=True,
        use_surface_term=True,
        use_static_stabilizer=False,
        use_grid_edf=True,
    )
    e0 = history[0]["total"]
    final = history[-1]
    drift = abs(final["total"] - e0) / max(abs(e0), 1.0)
    return float(drift), final


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--parameter-set", default="IQ3A", help="Skyrme parameter set, default: IQ3A")
    parser.add_argument("--grid-spacing", type=float, default=1.0, help="Grid spacing in fm")
    parser.add_argument("--n-sigma", type=float, default=3.0, help="Gaussian cutoff in packet widths")
    parser.add_argument("--seed", type=int, default=40, help="Base initialization seed")
    parser.add_argument(
        "--sigma-mode",
        choices=("code-default", "parameter-set"),
        default="code-default",
        help="Use initialize_nucleus default sigma_r=1.1 or sigma0+sigma1*A^(1/3)",
    )
    parser.add_argument(
        "--stability",
        action="store_true",
        help="Run a 500 fm/c grid-EDF propagation stability check for 40Ca",
    )
    parser.add_argument("--stability-time", type=float, default=500.0, help="Stability propagation time in fm/c")
    parser.add_argument("--dt", type=float, default=1.0, help="Stability time step in fm/c")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    edf = SkyrmeEDF.from_name(args.parameter_set)
    rows = []
    nuclei = []

    for index, bench in enumerate(BENCHMARKS):
        sigma_r = 1.1 if args.sigma_mode == "code-default" else iq_width(bench.a, args.parameter_set)
        nucleus = initialize_nucleus(
            bench.z,
            bench.a,
            sigma_r=sigma_r,
            seed=args.seed + index,
            edf=edf,
        )
        # Show the raw GridEDF result before the built-in empirical calibration.
        nucleus.energy_offset = 0.0
        components = raw_grid_components(nucleus, args.grid_spacing, args.n_sigma)
        rows.append((bench, components))
        nuclei.append((bench, nucleus))

    print("Raw GridEDF (energy_offset = 0)")
    print("nucleus  binding(MeV)  bind/A  rms(fm)  target_bind  target_rms  kinetic  potential")
    for bench, comp in rows:
        binding = -comp["total"]
        print(
            f"{bench.label:>7}  {binding:12.2f}  {binding / bench.a:6.2f}  "
            f"{comp['rms']:7.3f}  {bench.binding_mev:11.1f}  {bench.rms_fm:10.2f}  "
            f"{comp['kinetic']:7.2f}  {comp['potential']:8.2f}"
        )

    print("\nEnergy_offset calibration (built into initialize_nucleus)")
    for bench, nucleus in nuclei:
        # Re-apply the calibration computed by initialize_nucleus.
        nucleus.energy_offset = -bench.binding_mev - raw_grid_components(
            nucleus, args.grid_spacing, args.n_sigma
        )["total"]
        calibrated_binding = -nucleus.total_energy()
        print(
            f"{bench.label}: offset={nucleus.energy_offset:8.1f} MeV, "
            f"calibrated binding={calibrated_binding:8.1f} MeV "
            f"({calibrated_binding / bench.a:5.2f} MeV/A)"
        )

    if args.stability:
        bench, nucleus = nuclei[0]
        drift, final = run_stability(nucleus, args.stability_time, args.dt)
        print(
            f"\nStability {bench.label}: drift={100.0 * drift:.2f}% over {final['time']:.0f} fm/c, "
            f"final rms={final['rms']:.3f} fm, final binding={-final['total']:.1f} MeV"
        )


if __name__ == "__main__":
    main()
