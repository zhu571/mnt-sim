"""Compact statistical de-excitation helpers for ImQMD fragments."""

from __future__ import annotations

import os
import warnings
from dataclasses import dataclass

import numpy as np

_MASS_EXCESS_CACHE: dict[tuple[int, int], float] | None = None

NEUTRON_MASS_EXCESS = 8.071
HYDROGEN_MASS_EXCESS = 7.289
ALPHA_MASS_EXCESS = 2 * HYDROGEN_MASS_EXCESS + 2 * NEUTRON_MASS_EXCESS - 28.3


def _load_mass_excess() -> dict[tuple[int, int], float]:
    global _MASS_EXCESS_CACHE
    if _MASS_EXCESS_CACHE is not None:
        return _MASS_EXCESS_CACHE

    table: dict[tuple[int, int], float] = {}
    here = os.path.dirname(__file__)
    path = os.path.normpath(os.path.join(here, "..", "..", "hivap_fortran", "Mexcess95.dat"))
    try:
        with open(path, "r", encoding="utf-8") as f:
            for line in f:
                parts = line.split()
                if len(parts) < 3:
                    continue
                try:
                    A = int(parts[0])
                    Z = int(parts[1])
                    table[(Z, A)] = float(parts[2])
                except ValueError:
                    continue
    except OSError:
        table = {}
    _MASS_EXCESS_CACHE = table
    return table


def liquid_drop_binding_energy(Z: int, A: int) -> float:
    """Return a semi-empirical nuclear binding energy in MeV."""

    Z = int(Z)
    A = int(A)
    if A <= 1 or Z < 0 or Z > A:
        return 0.0
    N = A - Z
    av, asurf, ac, aa, ap = 15.56, 17.23, 0.697, 23.28, 11.2
    pairing = 0.0
    if A % 2 == 0:
        pairing = ap / np.sqrt(A) if Z % 2 == 0 and N % 2 == 0 else -ap / np.sqrt(A)
    binding = (
        av * A
        - asurf * A ** (2.0 / 3.0)
        - ac * Z * (Z - 1) / A ** (1.0 / 3.0)
        - aa * (A - 2 * Z) ** 2 / A
        + pairing
    )
    return float(max(binding, 0.0))


def mass_excess(Z: int, A: int) -> float:
    """Return mass excess in MeV from Mexcess95 or a liquid-drop fallback."""

    Z = int(Z)
    A = int(A)
    if A <= 0 or Z < 0 or Z > A:
        return float("inf")
    table = _load_mass_excess()
    if (Z, A) in table:
        return table[(Z, A)]
    N = A - Z
    return float(Z * HYDROGEN_MASS_EXCESS + N * NEUTRON_MASS_EXCESS - liquid_drop_binding_energy(Z, A))


def binding_energy(Z: int, A: int) -> float:
    """Return ground-state binding energy in MeV."""

    Z = int(Z)
    A = int(A)
    if A <= 1 or Z < 0 or Z > A:
        return 0.0
    me = mass_excess(Z, A)
    if np.isfinite(me):
        N = A - Z
        return float(max(Z * HYDROGEN_MASS_EXCESS + N * NEUTRON_MASS_EXCESS - me, 0.0))
    return liquid_drop_binding_energy(Z, A)


def _raw_neutron_separation_energy(Z: int, A: int) -> float:
    if A <= 1 or A - 1 < Z:
        return float("nan")
    return mass_excess(Z, A - 1) + NEUTRON_MASS_EXCESS - mass_excess(Z, A)


def _raw_proton_separation_energy(Z: int, A: int) -> float:
    if A <= 1 or Z <= 0 or Z > A:
        return float("nan")
    return mass_excess(Z - 1, A - 1) + HYDROGEN_MASS_EXCESS - mass_excess(Z, A)


def _raw_alpha_separation_energy(Z: int, A: int) -> float:
    if A <= 4 or Z < 2 or Z > A:
        return float("nan")
    return mass_excess(Z - 2, A - 4) + ALPHA_MASS_EXCESS - mass_excess(Z, A)


def neutron_separation_energy(Z: int, A: int) -> float:
    """Return one-neutron separation energy in MeV."""

    b_n = _raw_neutron_separation_energy(Z, A)
    if np.isfinite(b_n) and b_n > 0.0:
        return float(np.clip(b_n, 0.1, 20.0))
    if A <= 1 or A - 1 < Z:
        return 8.0
    b_n = binding_energy(Z, A) - binding_energy(Z, A - 1)
    return float(np.clip(b_n, 0.1, 15.0))


def _coulomb_barrier(z_emitted: int, z_daughter: int, a_daughter: int, a_emitted: int = 1) -> float:
    """Return a touching-spheres Coulomb barrier estimate in MeV."""

    if z_emitted <= 0 or z_daughter <= 0 or a_daughter <= 0:
        return 0.0
    r0 = 1.2
    radius = r0 * (a_daughter ** (1.0 / 3.0) + max(a_emitted, 1) ** (1.0 / 3.0))
    return float(1.44 * z_emitted * z_daughter / max(radius, 1.0e-9))


def proton_separation_energy(Z: int, A: int) -> float:
    """Return proton separation plus Coulomb barrier in MeV."""

    Z = int(Z)
    A = int(A)
    if A <= 1 or Z <= 0 or Z > A:
        return 12.0
    s_p = _raw_proton_separation_energy(Z, A)
    if not np.isfinite(s_p) or s_p <= 0.0:
        s_p = binding_energy(Z, A) - binding_energy(Z - 1, A - 1)
    s_p = float(np.clip(s_p, 0.1, 20.0))
    return s_p + _coulomb_barrier(1, Z - 1, A - 1, 1)


def alpha_separation_energy(Z: int, A: int) -> float:
    """Return alpha separation plus Coulomb barrier in MeV."""

    Z = int(Z)
    A = int(A)
    if A <= 4 or Z < 2 or Z > A:
        return 30.0
    s_alpha = _raw_alpha_separation_energy(Z, A)
    if not np.isfinite(s_alpha) or s_alpha <= 0.0:
        s_alpha = binding_energy(Z, A) - binding_energy(Z - 2, A - 4) - 28.3
    s_alpha = float(np.clip(s_alpha, 0.0, 35.0))
    return s_alpha + _coulomb_barrier(2, Z - 2, A - 4, 4)


def fission_barrier(Z: int, A: int) -> float:
    """Liquid-drop fission barrier with a compact shell-stabilization term."""

    Z = int(Z)
    A = int(A)
    if A <= 0 or Z <= 0:
        return float("inf")
    N = A - Z
    fissility = (Z * Z) / max(47.0 * A, 1.0)
    liquid_drop = 98.0 * max(1.0 - fissility, 0.0) ** 2
    shell = (
        2.0 * np.exp(-((N - 126.0) / 18.0) ** 2)
        + 1.0 * np.exp(-((Z - 82.0) / 12.0) ** 2)
    )
    return float(max(liquid_drop + shell, 0.25))


def _level_density_parameter(A: int, channel: str = "compound") -> float:
    scale = 1.0
    if channel == "fission":
        scale = 1.04
    return float(max(scale * A / 8.0, 1.0))


def _fermi_gas_level_density(U: np.ndarray | float, a: float) -> np.ndarray | float:
    energy = np.asarray(U, dtype=float)
    rho = np.zeros_like(energy, dtype=float)
    valid = energy > 1.0e-8
    if np.any(valid):
        e = energy[valid]
        exponent = np.clip(2.0 * np.sqrt(np.maximum(a * e, 0.0)), -80.0, 80.0)
        denominator = np.maximum(e, 1.0e-6) ** 1.25
        rho[valid] = np.exp(exponent) / denominator
    if np.isscalar(U):
        return float(rho)
    return rho


def _inverse_cross_section(A_daughter: int, particle: str, kinetic_energy: np.ndarray | float, barrier: float) -> np.ndarray:
    energy = np.asarray(kinetic_energy, dtype=float)
    geom = np.pi * (1.2 * max(A_daughter, 1) ** (1.0 / 3.0)) ** 2
    if particle == "n":
        sigma = geom * np.ones_like(energy, dtype=float)
    else:
        sigma = geom * np.maximum(energy, 0.0) / np.maximum(energy + max(barrier, 0.0), 1.0e-9)
    return sigma


def _weisskopf_width(A_daughter: int, E_available: float, a_daughter: float, particle: str, barrier: float) -> float:
    if E_available <= 0.0 or A_daughter <= 0:
        return 0.0
    n_grid = int(np.clip(48 + 6 * E_available, 48, 256))
    epsilon = np.linspace(0.0, E_available, n_grid)
    daughter_excitation = np.maximum(E_available - epsilon, 0.0)
    sigma_inv = _inverse_cross_section(A_daughter, particle, epsilon, barrier)
    rho_d = _fermi_gas_level_density(daughter_excitation, a_daughter)
    integral = np.trapezoid(sigma_inv * rho_d * epsilon, epsilon)
    return float(max(integral, 0.0))


def _gamma_width(A: int, E_star: float) -> float:
    if E_star <= 0.0:
        return 0.0
    a = _level_density_parameter(A)
    temperature = np.sqrt(max(E_star / a, 1.0e-9))
    return float(max(0.02 * max(A, 1) ** (2.0 / 3.0) * temperature**5, 1.0e-12))


def _bohr_wheeler_width(Z: int, A: int, E_star: float) -> float:
    if E_star <= 0.0 or Z < 70:
        return 0.0
    b_f = fission_barrier(Z, A)
    e_f = E_star - b_f
    if e_f <= 0.0:
        return 0.0
    a_p = _level_density_parameter(A)
    a_f = _level_density_parameter(A, channel="fission")
    rho_parent = _fermi_gas_level_density(E_star, a_p)
    if rho_parent <= 0.0:
        return 0.0
    n_grid = int(np.clip(48 + 6 * e_f, 48, 256))
    epsilon = np.linspace(0.0, e_f, n_grid)
    rho_saddle = _fermi_gas_level_density(np.maximum(e_f - epsilon, 0.0), a_f)
    integral = np.trapezoid(rho_saddle, epsilon)
    return float(max(integral / (2.0 * np.pi * rho_parent), 0.0))


def _fission_to_neutron_ratio(Z: int, A: int, E_star: float) -> float:
    if E_star <= 0.0 or Z < 70:
        return 0.0
    b_f = fission_barrier(Z, A)
    b_n = neutron_separation_energy(Z, A)
    e_f = E_star - b_f
    e_n = E_star - b_n
    if e_f <= 0.0:
        return 0.0
    if e_n <= 0.0:
        return float("inf")
    a_n = _level_density_parameter(A)
    a_f = _level_density_parameter(A, channel="fission")
    exponent = 2.0 * np.sqrt(max(a_f * e_f, 0.0)) - 2.0 * np.sqrt(max(a_n * e_n, 0.0))
    prefactor = 0.5 + 0.02 * max(Z - 90, 0)
    return float(prefactor * np.exp(np.clip(exponent, -80.0, 80.0)))


def _channel_width_dict(Z: int, A: int, E_star: float) -> dict[str, tuple[int, int, float]]:
    widths: dict[str, tuple[int, int, float]] = {}
    width_n = 0.0

    if A - 1 >= Z:
        b_n = neutron_separation_energy(Z, A)
        e_n = E_star - b_n
        width_n = _weisskopf_width(A - 1, e_n, _level_density_parameter(A - 1), "n", 0.0)
        if width_n > 0.0:
            widths["n"] = (Z, A - 1, width_n)

    if Z >= 1 and A - 1 >= Z - 1:
        b_p = proton_separation_energy(Z, A)
        e_p = E_star - b_p
        width_p = _weisskopf_width(
            A - 1,
            e_p,
            _level_density_parameter(A - 1),
            "p",
            _coulomb_barrier(1, Z - 1, A - 1, 1),
        )
        if width_p > 0.0:
            widths["p"] = (Z - 1, A - 1, width_p)

    if Z >= 2 and A >= 5 and A - 4 >= Z - 2:
        b_alpha = alpha_separation_energy(Z, A)
        e_alpha = E_star - b_alpha
        width_alpha = _weisskopf_width(
            A - 4,
            e_alpha,
            _level_density_parameter(A - 4),
            "alpha",
            _coulomb_barrier(2, Z - 2, A - 4, 4),
        )
        if width_alpha > 0.0:
            widths["alpha"] = (Z - 2, A - 4, width_alpha)

    widths["gamma"] = (Z, A, _gamma_width(A, E_star))

    if Z > 90:
        width_f = _bohr_wheeler_width(Z, A, E_star)
        ratio_floor = _fission_to_neutron_ratio(Z, A, E_star)
        if np.isfinite(ratio_floor) and width_n > 0.0:
            width_f = max(width_f, width_n * ratio_floor)
        elif not np.isfinite(ratio_floor):
            width_f = max(width_f, 1.0e12)
        if width_f > 0.0:
            widths["fission"] = (Z, A, width_f)

    return widths


def fission_competition(Z: int, A: int, E_star: float) -> float:
    """Return Bohr-Wheeler fission probability relative to neutron emission."""

    E_star = float(max(E_star, 0.0))
    widths = _channel_width_dict(int(Z), int(A), E_star)
    width_f = widths.get("fission", (0, 0, 0.0))[2]
    if width_f <= 0.0:
        return 0.0
    width_n = widths.get("n", (0, 0, 0.0))[2]
    if width_n <= 0.0:
        return 1.0
    return float(np.clip(width_f / max(width_f + width_n, 1.0e-12), 0.0, 1.0))


def weisskopf_evaporation_multi(Z: int, A: int, E_star: float) -> list[tuple[int, int, str, float]]:
    """Return one-step n, p, alpha, gamma, and fission branching ratios."""

    Z = int(Z)
    A = int(A)
    E_star = float(max(E_star, 0.0))
    if A <= 0 or Z < 0 or Z > A:
        return []
    if E_star <= 0.0:
        return [(Z, A, "gamma", 1.0)]

    channel_widths = _channel_width_dict(Z, A, E_star)
    widths = [(z_final, a_final, channel, width) for channel, (z_final, a_final, width) in channel_widths.items()]
    total = sum(width for *_, width in widths)
    if total <= 0.0 or not np.isfinite(total):
        return [(Z, A, "gamma", 1.0)]
    return [(z_final, a_final, channel, float(width / total)) for z_final, a_final, channel, width in widths]


def _channel_energy_cost(Z: int, A: int, channel: str, E_star: float) -> float:
    a = _level_density_parameter(A)
    temperature = np.sqrt(max(E_star / a, 1.0e-9))
    kinetic = 3.0 * temperature
    if channel == "n":
        return neutron_separation_energy(Z, A) + kinetic
    if channel == "p":
        return proton_separation_energy(Z, A) + kinetic
    if channel == "alpha":
        return alpha_separation_energy(Z, A) + kinetic
    if channel == "gamma":
        return min(E_star, max(1.0, 2.0 * temperature))
    return E_star


@dataclass(frozen=True)
class DecayStep:
    channel: str
    Z_before: int
    A_before: int
    E_star_before: float
    Z_after: int
    A_after: int
    E_star_after: float


@dataclass(frozen=True)
class EvaporationChainResult:
    residue_Z: int
    residue_A: int
    residue_E_star: float
    steps: tuple[DecayStep, ...]
    terminated_by: str


def evaporate_chain(
    Z: int,
    A: int,
    E_star: float,
    max_steps: int = 64,
    rng: np.random.Generator | None = None,
) -> EvaporationChainResult:
    """Sample a full evaporation chain and return the full history."""

    Z = int(Z)
    A = int(A)
    E_star = float(max(E_star, 0.0))
    rng = rng if rng is not None else np.random.default_rng()
    steps: list[DecayStep] = []
    terminated_by = "energy_floor"

    for _ in range(max(0, int(max_steps))):
        if A <= 1 or E_star < 1.0:
            break
        channels = weisskopf_evaporation_multi(Z, A, E_star)
        if not channels:
            terminated_by = "no_channels"
            break
        probabilities = np.asarray([prob for *_, prob in channels], dtype=float)
        probabilities /= probabilities.sum()
        idx = int(rng.choice(len(channels), p=probabilities))
        z_final, a_final, channel, _ = channels[idx]
        if channel == "fission":
            terminated_by = "fission"
            steps.append(
                DecayStep(
                    channel=channel,
                    Z_before=Z,
                    A_before=A,
                    E_star_before=E_star,
                    Z_after=Z,
                    A_after=A,
                    E_star_after=0.0,
                )
            )
            E_star = 0.0
            break
        cost = _channel_energy_cost(Z, A, channel, E_star)
        next_E_star = max(E_star - cost, 0.0)
        steps.append(
            DecayStep(
                channel=channel,
                Z_before=Z,
                A_before=A,
                E_star_before=E_star,
                Z_after=int(z_final),
                A_after=int(a_final),
                E_star_after=next_E_star,
            )
        )
        Z, A = int(z_final), int(a_final)
        E_star = next_E_star
    else:
        terminated_by = "max_steps"

    return EvaporationChainResult(
        residue_Z=Z,
        residue_A=A,
        residue_E_star=E_star,
        steps=tuple(steps),
        terminated_by=terminated_by,
    )


def evaporate_full(
    Z: int,
    A: int,
    E_star: float,
    max_steps: int = 64,
    rng: np.random.Generator | None = None,
) -> tuple[int, int]:
    """Sample a multi-step evaporation chain and return the surviving residue."""

    result = evaporate_chain(Z, A, E_star, max_steps=max_steps, rng=rng)
    return result.residue_Z, result.residue_A


def weisskopf_evaporation(Z: int, A: int, E_star: float, n_max: int = 10) -> list[tuple[int, int, float]]:
    """Return neutron evaporation channels ``(Z, A_final, probability)``.

    Only neutron evaporation is included.  The no-evaporation channel is kept;
    callers can evaluate fission separately with ``fission_competition``.
    """

    Z = int(Z)
    A = int(A)
    E_star = float(max(E_star, 0.0))
    if A <= 0:
        return []
    if E_star <= 1.0 or n_max <= 0:
        return [(Z, A, 1.0)]

    a = max(A / 8.0, 1.0)
    temperature = np.sqrt(max(E_star / a, 1.0e-6))
    channels: list[tuple[int, int, float]] = []
    total = 0.0
    max_neutrons = min(int(n_max), A - Z, A - 1)
    cumulative_cost = 0.0

    for n in range(max_neutrons + 1):
        if n > 0:
            cumulative_cost += neutron_separation_energy(Z, A - n + 1)
        E_rem = E_star - cumulative_cost
        if E_rem < 0.0:
            break
        if n == 0:
            weight = np.exp(-E_star / max(temperature, 1.0e-6))
        else:
            phase = (max(E_rem, 1.0e-9) / max(E_star, 1.0e-9)) ** (0.5 * a)
            # The leading E*/B factor represents the available sequential
            # emission window and keeps moderate 20 MeV cases from becoming
            # artificially dominated by the unevaporated channel.
            window = (E_star / max(cumulative_cost, 1.0)) ** n
            weight = np.exp(-cumulative_cost / max(temperature, 1.0e-6)) * phase * window
        channels.append((Z, A - n, float(weight)))
        total += float(weight)

    if total <= 0.0:
        return [(Z, A, 1.0)]

    return [(z, a_final, p / total) for z, a_final, p in channels]


def neutron_separation_discrepancies(
    nuclei: list[tuple[int, int]] | tuple[tuple[int, int], ...] | None = None,
    tolerance_mev: float = 2.0,
) -> list[dict[str, float | int]]:
    """Return neutron-separation discrepancies against raw Mexcess95 values."""

    table = _load_mass_excess()
    if nuclei is None:
        nuclei = sorted((Z, A) for (Z, A) in table if (Z, A - 1) in table)
    discrepancies: list[dict[str, float | int]] = []
    for Z, A in nuclei:
        raw = _raw_neutron_separation_energy(int(Z), int(A))
        model = neutron_separation_energy(int(Z), int(A))
        if not np.isfinite(raw):
            continue
        delta = abs(model - raw)
        if delta > float(tolerance_mev):
            discrepancies.append(
                {
                    "Z": int(Z),
                    "A": int(A),
                    "raw_mev": float(raw),
                    "model_mev": float(model),
                    "delta_mev": float(delta),
                }
            )
    return discrepancies


def validate_neutron_separation_energies(
    nuclei: list[tuple[int, int]] | tuple[tuple[int, int], ...] | None = None,
    tolerance_mev: float = 2.0,
    warn: bool = True,
) -> list[dict[str, float | int]]:
    """Validate neutron separation energies against the mass table."""

    discrepancies = neutron_separation_discrepancies(nuclei=nuclei, tolerance_mev=tolerance_mev)
    if warn and discrepancies:
        sample = ", ".join(
            f"(Z={item['Z']},A={item['A']},Δ={item['delta_mev']:.2f} MeV)"
            for item in discrepancies[:5]
        )
        warnings.warn(
            f"Found {len(discrepancies)} neutron-separation discrepancies above {tolerance_mev:.1f} MeV: {sample}",
            RuntimeWarning,
            stacklevel=2,
        )
    return discrepancies


def benchmark_evaporation_chain(
    test_points: list[tuple[int, int, float]] | None = None,
    n_trials: int = 256,
    seed: int = 12345,
) -> list[dict[str, float | int | str]]:
    """Run a simple statistical benchmark for representative evaporation chains."""

    if test_points is None:
        test_points = [
            (28, 64, 20.0),
            (50, 132, 40.0),
            (92, 238, 60.0),
        ]
    rng = np.random.default_rng(seed)
    summaries: list[dict[str, float | int | str]] = []
    for Z, A, E_star in test_points:
        chains = [evaporate_chain(Z, A, E_star, rng=rng) for _ in range(int(n_trials))]
        residue_a = np.asarray([chain.residue_A for chain in chains], dtype=float)
        residue_z = np.asarray([chain.residue_Z for chain in chains], dtype=float)
        step_counts = np.asarray([len(chain.steps) for chain in chains], dtype=float)
        fission_fraction = float(np.mean([chain.terminated_by == "fission" for chain in chains]))
        gamma_terminal_fraction = float(np.mean([chain.steps[-1].channel == "gamma" if chain.steps else True for chain in chains]))
        summaries.append(
            {
                "Z": int(Z),
                "A": int(A),
                "E_star": float(E_star),
                "n_trials": int(n_trials),
                "mean_residue_Z": float(np.mean(residue_z)),
                "mean_residue_A": float(np.mean(residue_a)),
                "mean_steps": float(np.mean(step_counts)),
                "fission_fraction": fission_fraction,
                "gamma_terminal_fraction": gamma_terminal_fraction,
            }
        )
    return summaries


__all__ = [
    "DecayStep",
    "EvaporationChainResult",
    "alpha_separation_energy",
    "benchmark_evaporation_chain",
    "binding_energy",
    "evaporate_chain",
    "evaporate_full",
    "fission_barrier",
    "fission_competition",
    "liquid_drop_binding_energy",
    "mass_excess",
    "neutron_separation_energy",
    "neutron_separation_discrepancies",
    "proton_separation_energy",
    "validate_neutron_separation_energies",
    "weisskopf_evaporation",
    "weisskopf_evaporation_multi",
]
