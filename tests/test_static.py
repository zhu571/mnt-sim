"""Phase 1 static-nucleus stability checks."""

from __future__ import annotations

import pathlib
import sys

ROOT = pathlib.Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from mnt_sim.imqmd import initialize_nucleus, propagate
from mnt_sim.imqmd.initializer import empirical_binding_per_nucleon


def _run_stability(Z: int, A: int, dt: float = 2.0, time: float = 2000.0):
    nucleus = initialize_nucleus(Z, A, sigma_r=1.1)
    initial_rms = sum(nucleus.rms_radius()) / 2.0
    history = propagate(nucleus, dt=dt, n_steps=int(time / dt), sample_every=100)
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
    history = propagate(nucleus, dt=1.0, n_steps=1000, sample_every=100)
    e0 = history[0]["total"]
    e1 = history[-1]["total"]
    assert abs(e1 - e0) / max(abs(e0), 1.0) < 0.01


if __name__ == "__main__":
    for test in (test_ca40_stability, test_pb208_stability, test_single_nucleus_energy_conservation):
        test()
    ca_hist, ca_bind, ca_r0, ca_r1 = _run_stability(20, 40)
    pb_hist, pb_bind, pb_r0, pb_r1 = _run_stability(82, 208)
    drift = abs(ca_hist[-1]["total"] - ca_hist[0]["total"]) / abs(ca_hist[0]["total"]) * 100.0
    print("Static ImQMD Phase 1 checks passed")
    print(f"40Ca: rms {ca_r0:.3f} -> {ca_r1:.3f} fm, binding {ca_bind:.3f} MeV/A")
    print(f"208Pb: rms {pb_r0:.3f} -> {pb_r1:.3f} fm, binding {pb_bind:.3f} MeV/A")
    print(f"40Ca total-energy drift over 2000 fm/c: {drift:.3f}%")
