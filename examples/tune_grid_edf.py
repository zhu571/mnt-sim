"""Tune the global GridEDF potential scale against benchmark nuclei."""

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
        eta=1.0,
    )
    components = grid.total_energy(nucleus.positions, nucleus.momenta, nucleus.is_proton)
    components["rms"] = mass_rms(nucleus)
    return components


def binding_at_eta(components: dict[str, float], eta: float) -> float:
    total = components["kinetic"] + eta * components["potential_unscaled"]
    return -float(total)


def eta_for_target(components: dict[str, float], binding_mev: float) -> float:
    target_total = -float(binding_mev)
    return float((target_total - components["kinetic"]) / components["potential_unscaled"])


def eta_interval_for_binding(components: dict[str, float], bench: Benchmark) -> tuple[float, float]:
    raw_attraction = -float(components["potential_unscaled"])
    lower_binding = bench.binding_mev - bench.binding_tolerance_mev
    upper_binding = bench.binding_mev + bench.binding_tolerance_mev
    return (
        max(0.0, (lower_binding + components["kinetic"]) / raw_attraction),
        max(0.0, (upper_binding + components["kinetic"]) / raw_attraction),
    )


def print_baseline(rows: list[tuple[Benchmark, dict[str, float]]]) -> None:
    print("Baseline GridEDF at eta=1.0")
    print("nucleus  binding(MeV)  bind/A  rms(fm)  target_bind  target_rms  kinetic  pot_raw")
    for bench, comp in rows:
        binding = -comp["total"]
        print(
            f"{bench.label:>7}  {binding:12.2f}  {binding / bench.a:6.2f}  "
            f"{comp['rms']:7.3f}  {bench.binding_mev:11.1f}  {bench.rms_fm:10.2f}  "
            f"{comp['kinetic']:7.2f}  {comp['potential_unscaled']:8.2f}"
        )


def print_sweep(rows: list[tuple[Benchmark, dict[str, float]]], etas: list[float]) -> None:
    print()
    print("Eta sweep, using E_total = kinetic + eta * potential_raw")
    header = "eta    " + "  ".join(f"{bench.label:>18}" for bench, _ in rows)
    print(header)
    print("       " + "  ".join("bind(MeV)  bind/A".rjust(18) for _ in rows))
    for eta in etas:
        values = []
        for bench, comp in rows:
            binding = binding_at_eta(comp, eta)
            values.append(f"{binding:9.1f}  {binding / bench.a:6.2f}")
        print(f"{eta:5.3f}  " + "  ".join(values))


def run_stability(nucleus, eta: float, time_fm_c: float, dt: float) -> tuple[float, dict[str, float]]:
    trial = nucleus.copy()
    trial.energy_offset = 0.0
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
        grid_eta=eta,
    )
    e0 = history[0]["total"]
    final = history[-1]
    drift = abs(final["total"] - e0) / max(abs(e0), 1.0)
    return float(drift), final


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--parameter-set", default="IQ2", help="Skyrme parameter set, default: IQ2")
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
        "--eta",
        type=float,
        default=None,
        help="Report this eta in addition to the Ca-fitted value",
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
        nucleus.energy_offset = 0.0
        components = raw_grid_components(nucleus, args.grid_spacing, args.n_sigma)
        rows.append((bench, components))
        nuclei.append((bench, nucleus))

    print_baseline(rows)

    ca_bench, ca_components = rows[0]
    fitted_eta = eta_for_target(ca_components, ca_bench.binding_mev)
    intervals = [(bench, eta_interval_for_binding(comp, bench)) for bench, comp in rows]
    overlap_low = max(interval[0] for _, interval in intervals)
    overlap_high = min(interval[1] for _, interval in intervals)
    recommended_eta = 0.5 * (overlap_low + overlap_high) if overlap_low <= overlap_high else fitted_eta

    etas = [0.30, 0.35, 0.40, 0.45, 0.50, 0.55, recommended_eta, 0.60, fitted_eta]
    if args.eta is not None:
        etas.append(args.eta)
    etas = sorted(set(round(float(eta), 6) for eta in etas if eta >= 0.0))
    print_sweep(rows, etas)

    print()
    print(f"Ca-fitted eta: {fitted_eta:.6f}")
    for bench, interval in intervals:
        print(
            f"{bench.label} acceptable eta for binding tolerance: "
            f"{interval[0]:.6f} .. {interval[1]:.6f}"
        )
    if overlap_low <= overlap_high:
        print(f"Recommended eta from tolerance overlap: {recommended_eta:.6f}")
    else:
        print("No shared eta interval satisfies all binding tolerances; using the Ca-fitted eta.")

    for bench, comp in rows:
        binding = binding_at_eta(comp, recommended_eta)
        print(
            f"{bench.label}: binding={binding:.1f} MeV ({binding / bench.a:.2f} MeV/A), "
            f"rms={comp['rms']:.3f} fm; targets {bench.binding_mev:.1f} MeV, {bench.rms_fm:.2f} fm"
        )

    if args.stability:
        bench, nucleus = nuclei[0]
        drift, final = run_stability(nucleus, recommended_eta, args.stability_time, args.dt)
        print()
        print(
            f"Stability {bench.label}: drift={100.0 * drift:.2f}% over {final['time']:.0f} fm/c, "
            f"final rms={final['rms']:.3f} fm, final binding={-final['total']:.1f} MeV"
        )


if __name__ == "__main__":
    main()
