"""Initialize and briefly evolve a static ImQMD nucleus."""

from __future__ import annotations

import pathlib
import sys

ROOT = pathlib.Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from mnt_sim.imqmd import initialize_nucleus, propagate


def main() -> None:
    nucleus = initialize_nucleus(20, 40, sigma_r=1.1)
    rp0, rn0 = nucleus.rms_radius()
    e0 = nucleus.total_energy()
    history = propagate(nucleus, dt=1.0, n_steps=200, sample_every=20)
    final = history[-1]
    drift = (final["total"] - e0) / abs(e0) * 100.0

    print("40Ca static ImQMD example")
    print(f"Initial binding: {-e0 / nucleus.A:.3f} MeV/A")
    print(f"Initial rms radii: proton={rp0:.3f} fm neutron={rn0:.3f} fm")
    print(f"Final rms radii: proton={final['rms_p']:.3f} fm neutron={final['rms_n']:.3f} fm")
    print(f"Energy drift over {final['time']:.0f} fm/c: {drift:.4f}%")
    print(f"Max centroid radius: {final['max_radius']:.3f} fm")


if __name__ == "__main__":
    main()
