"""Compact statistical de-excitation helpers for ImQMD fragments."""

from __future__ import annotations

import os

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


def neutron_separation_energy(Z: int, A: int) -> float:
    """Return one-neutron separation energy in MeV."""

    if A <= 1 or A - 1 < Z:
        return 8.0
    b_n = mass_excess(Z, A - 1) + NEUTRON_MASS_EXCESS - mass_excess(Z, A)
    if not np.isfinite(b_n) or b_n <= 0.0:
        b_n = binding_energy(Z, A) - binding_energy(Z, A - 1)
    return float(np.clip(b_n, 3.0, 12.0))


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
    s_p = mass_excess(Z - 1, A - 1) + HYDROGEN_MASS_EXCESS - mass_excess(Z, A)
    if not np.isfinite(s_p) or s_p <= 0.0:
        s_p = binding_energy(Z, A) - binding_energy(Z - 1, A - 1)
    s_p = float(np.clip(s_p, 1.0, 18.0))
    return s_p + _coulomb_barrier(1, Z - 1, A - 1, 1)


def alpha_separation_energy(Z: int, A: int) -> float:
    """Return alpha separation plus Coulomb barrier in MeV."""

    Z = int(Z)
    A = int(A)
    if A <= 4 or Z < 2 or Z > A:
        return 30.0
    s_alpha = mass_excess(Z - 2, A - 4) + ALPHA_MASS_EXCESS - mass_excess(Z, A)
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


def fission_competition(Z: int, A: int, E_star: float) -> float:
    """Return Bohr-Wheeler fission probability relative to neutron emission."""

    E_star = float(max(E_star, 0.0))
    if E_star <= 0.0 or Z < 70:
        return 0.0
    a_n = max(A / 8.0, 1.0)
    a_f = 1.02 * a_n
    b_f = fission_barrier(Z, A)
    b_n = neutron_separation_energy(Z, A)
    e_f = E_star - b_f
    e_n = E_star - b_n
    if e_f <= 0.0:
        return 0.0
    if e_n <= 0.0:
        return 1.0
    exponent = 2.0 * np.sqrt(max(a_f * e_f, 0.0)) - 2.0 * np.sqrt(max(a_n * e_n, 0.0))
    ratio = 0.5 * np.exp(float(np.clip(exponent, -80.0, 80.0)))
    return float(np.clip(ratio / (1.0 + ratio), 0.0, 1.0))


def _evaporation_width(A: int, e_available: float, level_density_a: float, barrier_scale: float = 1.0) -> float:
    """Relative Weisskopf-Ewing evaporation width."""

    if e_available <= 0.0 or A <= 1:
        return 0.0
    phase_space = max(e_available, 1.0e-9)
    radius_factor = max(A, 1) ** (2.0 / 3.0)
    exponent = 2.0 * np.sqrt(max(level_density_a * e_available, 0.0))
    return float(radius_factor * phase_space * barrier_scale * np.exp(np.clip(exponent, -80.0, 80.0)))


def weisskopf_evaporation_multi(Z: int, A: int, E_star: float) -> list[tuple[int, int, str, float]]:
    """Return one-step n, p, alpha, gamma, and fission branching ratios."""

    Z = int(Z)
    A = int(A)
    E_star = float(max(E_star, 0.0))
    if A <= 0 or Z < 0 or Z > A:
        return []
    if E_star <= 0.0:
        return [(Z, A, "gamma", 1.0)]

    widths: list[tuple[int, int, str, float]] = []

    if A - 1 >= Z:
        b_n = neutron_separation_energy(Z, A)
        e_n = E_star - b_n
        width_n = _evaporation_width(A - 1, e_n, max((A - 1) / 8.0, 1.0))
        if width_n > 0.0:
            widths.append((Z, A - 1, "n", width_n))
    else:
        width_n = 0.0

    if Z >= 1 and A - 1 >= Z - 1:
        b_p = proton_separation_energy(Z, A)
        e_p = E_star - b_p
        width_p = _evaporation_width(A - 1, e_p, max((A - 1) / 8.0, 1.0), barrier_scale=0.6)
        if width_p > 0.0:
            widths.append((Z - 1, A - 1, "p", width_p))

    if Z >= 2 and A >= 5 and A - 4 >= Z - 2:
        b_alpha = alpha_separation_energy(Z, A)
        e_alpha = E_star - b_alpha
        width_alpha = _evaporation_width(A - 4, e_alpha, max((A - 4) / 8.0, 1.0), barrier_scale=0.35)
        if width_alpha > 0.0:
            widths.append((Z - 2, A - 4, "alpha", width_alpha))

    gamma_reference = width_n if width_n > 0.0 else _evaporation_width(A, E_star, max(A / 8.0, 1.0))
    gamma_width = max(0.1 * gamma_reference, 1.0e-12)
    widths.append((Z, A, "gamma", gamma_width))

    if Z > 90:
        fission_probability = fission_competition(Z, A, E_star)
        if fission_probability >= 1.0:
            widths.append((Z, A, "fission", 1.0e80))
        elif fission_probability > 0.0:
            non_fission_width = sum(width for *_, width in widths)
            fission_width = non_fission_width * fission_probability / max(1.0 - fission_probability, 1.0e-12)
            widths.append((Z, A, "fission", fission_width))

    total = sum(width for *_, width in widths)
    if total <= 0.0 or not np.isfinite(total):
        return [(Z, A, "gamma", 1.0)]
    return [(z_final, a_final, channel, float(width / total)) for z_final, a_final, channel, width in widths]


def _channel_energy_cost(Z: int, A: int, channel: str, E_star: float) -> float:
    a = max(A / 8.0, 1.0)
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


def evaporate_full(
    Z: int,
    A: int,
    E_star: float,
    max_steps: int = 10,
    rng: np.random.Generator | None = None,
) -> tuple[int, int]:
    """Sample a multi-step evaporation chain and return the surviving residue."""

    Z = int(Z)
    A = int(A)
    E_star = float(max(E_star, 0.0))
    rng = rng if rng is not None else np.random.default_rng()

    for _ in range(max(0, int(max_steps))):
        if A <= 1 or E_star <= 1.0:
            break
        channels = weisskopf_evaporation_multi(Z, A, E_star)
        if not channels:
            break
        probabilities = np.asarray([prob for *_, prob in channels], dtype=float)
        probabilities /= probabilities.sum()
        idx = int(rng.choice(len(channels), p=probabilities))
        z_final, a_final, channel, _ = channels[idx]
        if channel == "fission":
            break
        cost = _channel_energy_cost(Z, A, channel, E_star)
        Z, A = int(z_final), int(a_final)
        E_star = max(E_star - cost, 0.0)
        if channel == "gamma" and E_star <= 1.0:
            break

    return Z, A


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


__all__ = [
    "alpha_separation_energy",
    "binding_energy",
    "evaporate_full",
    "fission_barrier",
    "fission_competition",
    "liquid_drop_binding_energy",
    "mass_excess",
    "neutron_separation_energy",
    "proton_separation_energy",
    "weisskopf_evaporation",
    "weisskopf_evaporation_multi",
]
