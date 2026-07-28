"""Phase 2 collision, Pauli blocking, and Fermi-constraint checks."""

from __future__ import annotations

import pathlib
import sys

import numpy as np

ROOT = pathlib.Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from mnt_sim.imqmd import (
    attempt_nn_collision,
    fermi_constraint_check,
    free_nn_cross_section,
    in_medium_factor,
    initialize_nucleus,
    pauli_blocking_probability,
    propagate,
)


def test_free_cross_sections():
    assert free_nn_cross_section(25.0, "np") > free_nn_cross_section(25.0, "pp")
    assert free_nn_cross_section(25.0, "nn") == free_nn_cross_section(25.0, "pp")


def test_in_medium_scaling():
    # Enhancement form sigma_med = (1 + eta*sqrt(s)*rho/rho0)*sigma_free
    # (research report Sec. 4.4 / Chen 2024 Eq. (5)): factor > 1 at rho0 and
    # growing with sqrt(s).  At threshold sqrt(s) = 2*M_N ~ 1.878 GeV.
    f_thr = in_medium_factor(0.165, "np")
    assert f_thr > 1.0
    assert np.isclose(f_thr, 1.0 + 0.2 * (2.0 * 938.9 / 1000.0), rtol=1e-3)
    # sqrt(s) dependence: higher CM energy -> larger enhancement
    assert in_medium_factor(0.165, "np", e_cm=100.0) > f_thr
    # density dependence: vacuum limit is the free cross section
    assert np.isclose(in_medium_factor(0.0, "np"), 1.0)


def test_pauli_blocking_in_nucleus():
    nucleus = initialize_nucleus(20, 40, sigma_r=1.1)
    nucleus._active_collision_pair = (0, 1)
    p1, p2 = nucleus.momenta[0].copy(), nucleus.momenta[1].copy()
    p_block = pauli_blocking_probability(p1, p2, p1, p2, nucleus)
    del nucleus._active_collision_pair
    assert p_block > 0.5

    attempted = blocked = accepted = 0
    for _ in range(20):
        stats = attempt_nn_collision(nucleus, dt=1.0)
        attempted += stats["attempted"]
        blocked += stats["blocked"]
        accepted += stats["accepted"]
    if attempted:
        assert blocked >= accepted


def test_fermi_constraint_stability():
    nucleus = initialize_nucleus(20, 40, sigma_r=1.1)
    momenta = nucleus.momenta
    same_species = np.where(nucleus.is_proton == nucleus.is_proton[0])[0]
    momenta[same_species[1]] = momenta[same_species[0]]
    nucleus.momenta = momenta

    corrections = fermi_constraint_check(nucleus)
    assert corrections > 0

    history = propagate(
        nucleus,
        dt=1.0,
        n_steps=100,
        sample_every=50,
        with_collisions=True,
        use_surface_term=True,
        use_static_stabilizer=True,
    )
    assert history[-1]["max_radius"] < 15.0


if __name__ == "__main__":
    for test in (
        test_free_cross_sections,
        test_in_medium_scaling,
        test_pauli_blocking_in_nucleus,
        test_fermi_constraint_stability,
    ):
        test()
    print("Collision checks passed")
