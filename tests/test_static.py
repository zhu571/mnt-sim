"""Phase 1 static-nucleus stability checks."""

from __future__ import annotations

import pathlib
import sys

import numpy as np

ROOT = pathlib.Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from mnt_sim.imqmd import GaussianPacket, ImQMDNucleus, SkyrmeEDF, initialize_and_relax, initialize_nucleus, propagate
from mnt_sim.imqmd.initializer import empirical_binding_per_nucleon
from mnt_sim.imqmd.propagator import _derivatives


def _run_stability(Z: int, A: int, dt: float = 2.0, time: float = 2000.0):
    nucleus = initialize_nucleus(Z, A, sigma_r=1.1)
    initial_rms = sum(nucleus.rms_radius()) / 2.0
    history = propagate(
        nucleus,
        dt=dt,
        n_steps=int(time / dt),
        sample_every=100,
        use_surface_term=True,
        use_static_stabilizer=True,
        use_grid_edf=True,
    )
    final = history[-1]
    final_rms = 0.5 * (final["rms_p"] + final["rms_n"])
    binding_per_a = -final["total"] / A
    expected = empirical_binding_per_nucleon(Z, A)

    assert abs(final_rms - initial_rms) / initial_rms < 0.10
    assert abs(binding_per_a - expected) < 2.0
    assert final["max_radius"] < 15.0
    return history, binding_per_a, initial_rms, final_rms


def test_ca40_stability():
    _run_stability(20, 40)


def test_pb208_stability():
    _run_stability(82, 208)


def test_single_nucleus_energy_conservation():
    nucleus = initialize_nucleus(20, 40, sigma_r=1.1)
    history = propagate(
        nucleus,
        dt=1.0,
        n_steps=1000,
        sample_every=100,
        use_surface_term=True,
        use_static_stabilizer=True,
        use_grid_edf=True,
    )
    e0 = history[0]["total"]
    e1 = history[-1]["total"]
    assert abs(e1 - e0) / max(abs(e0), 1.0) < 0.01


def test_static_no_spring():
    nucleus = initialize_and_relax(20, 40, sigma_r=1.1, seed=40)
    initial_rms = sum(nucleus.rms_radius()) / 2.0
    history = propagate(
        nucleus,
        dt=1.0,
        n_steps=500,
        sample_every=100,
        with_collisions=False,
        remove_cm_drift=True,
        use_surface_term=True,
        use_static_stabilizer=False,
        use_grid_edf=True,
    )
    final = history[-1]
    final_rms = 0.5 * (final["rms_p"] + final["rms_n"])
    e0 = history[0]["total"]
    e1 = final["total"]

    assert abs(e1 - e0) / max(abs(e0), 1.0) < 0.15


def test_initialization_keeps_physical_grid_scale():
    nucleus = initialize_nucleus(20, 40, sigma_r=1.1, seed=40)

    assert nucleus.grid_nuclear_scale == 1.0


def test_collision_preserves_offsets_without_scaling_mean_field(monkeypatch):
    from mnt_sim.imqmd import reaction

    def prepared(Z, A, sigma_r, seed, relax_time, edf):
        packet = GaussianPacket(np.zeros(3), np.zeros(3), sigma_r, bool(Z))
        return ImQMDNucleus(Z, A - Z, [packet], energy_offset=10.0 + Z)

    monkeypatch.setattr(reaction, "_prepared_nucleus", prepared)
    system = reaction.make_collision_event(1, 1, 0, 1, 1.0, 0.0, 1, 2, 3, separation=10.0)

    assert system.grid_nuclear_scale == 1.0
    assert system.energy_offset == 21.0


def test_physical_energy_force_consistency():
    packets = [
        GaussianPacket(np.array([-1.0, 0.0, 0.0]), np.zeros(3), 1.1, True),
        GaussianPacket(np.array([1.0, 0.1, 0.0]), np.zeros(3), 1.1, False),
        GaussianPacket(np.array([0.0, 1.0, 0.2]), np.zeros(3), 1.1, True),
    ]
    nucleus = ImQMDNucleus(2, 1, packets, edf=SkyrmeEDF(static_k=0.0))
    positions = nucleus.positions
    _, force = _derivatives(
        nucleus,
        positions,
        nucleus.momenta,
        use_surface_term=True,
        use_static_stabilizer=False,
    )

    eps = 1.0e-5
    for index, axis in ((0, 0), (1, 1), (2, 2)):
        plus = positions.copy()
        minus = positions.copy()
        plus[index, axis] += eps
        minus[index, axis] -= eps
        nucleus.positions = plus
        e_plus = nucleus.energy_components(use_surface_term=True, use_static_stabilizer=False)["total"]
        nucleus.positions = minus
        e_minus = nucleus.energy_components(use_surface_term=True, use_static_stabilizer=False)["total"]
        nucleus.positions = positions
        numerical_gradient = (e_plus - e_minus) / (2.0 * eps)
        assert abs(numerical_gradient + force[index, axis]) < 1.0e-5


def test_final_state_pauli():
    nucleus = initialize_nucleus(20, 40, sigma_r=1.1)
    initial_rms = sum(nucleus.rms_radius()) / 2.0
    history = propagate(
        nucleus,
        dt=1.0,
        n_steps=500,
        sample_every=100,
        with_collisions=True,
        collision_dt=1.0,
        use_surface_term=True,
        use_static_stabilizer=True,
        use_grid_edf=True,
    )
    final = history[-1]
    final_rms = 0.5 * (final["rms_p"] + final["rms_n"])
    stats = getattr(nucleus, "collision_stats", {"attempted": 0, "blocked": 0, "accepted": 0})

    assert abs(final_rms - initial_rms) / initial_rms < 0.12
    assert final["max_radius"] < 15.0
    assert stats["blocked"] >= stats["accepted"]


def test_collision_rate_sanity():
    nucleus = initialize_nucleus(10, 20, sigma_r=1.1)
    positions = nucleus.positions
    momenta = nucleus.momenta
    half = nucleus.A // 2
    positions[:half, 0] -= 2.0
    positions[half:, 0] += 2.0
    momenta[:half, 0] += 220.0
    momenta[half:, 0] -= 220.0
    nucleus.positions = positions
    nucleus.momenta = momenta - momenta.mean(axis=0)

    propagate(
        nucleus,
        dt=1.0,
        n_steps=20,
        sample_every=20,
        with_collisions=True,
        collision_dt=1.0,
        use_surface_term=True,
        use_static_stabilizer=True,
        use_grid_edf=True,
    )
    stats = getattr(nucleus, "collision_stats", {"attempted": 0, "blocked": 0, "accepted": 0})
    assert stats["attempted"] > 0


if __name__ == "__main__":
    for test in (
        test_ca40_stability,
        test_pb208_stability,
        test_single_nucleus_energy_conservation,
        test_static_no_spring,
        test_physical_energy_force_consistency,
        test_final_state_pauli,
        test_collision_rate_sanity,
    ):
        test()
    ca_hist, ca_bind, ca_r0, ca_r1 = _run_stability(20, 40)
    pb_hist, pb_bind, pb_r0, pb_r1 = _run_stability(82, 208)
    drift = abs(ca_hist[-1]["total"] - ca_hist[0]["total"]) / abs(ca_hist[0]["total"]) * 100.0
    print("Static ImQMD Phase 1 checks passed")
    print(f"40Ca: rms {ca_r0:.3f} -> {ca_r1:.3f} fm, binding {ca_bind:.3f} MeV/A")
    print(f"208Pb: rms {pb_r0:.3f} -> {pb_r1:.3f} fm, binding {pb_bind:.3f} MeV/A")
    print(f"40Ca total-energy drift over 2000 fm/c: {drift:.3f}%")
