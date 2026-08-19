"""Sensitivity scan of the dynamics->decay switch time (research report Sec. 4.6).

The imQMD+GEMINI/HIVAP MNT studies find an optimal switch time of about
500 fm/c: primary fragments are formed but still excited.  This script runs
a light system at several switch times and compares the heavy-fragment
charge/mass statistics so the choice can be checked for this code base.

Usage:
    python3 examples/switch_time_scan.py [n_events] [b_fm]

Default: 2 events at b = 3 fm for 40Ca+40Ca at 8 MeV/u, switch times
300 / 500 / 800 fm/c.
"""

from __future__ import annotations

import pathlib
import sys
import time

import numpy as np

ROOT = pathlib.Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from mnt_sim.imqmd import run_imqmd_event

SWITCH_TIMES_FM_C = (300.0, 500.0, 800.0)


def main() -> None:
    n_events = int(sys.argv[1]) if len(sys.argv) > 1 else 2
    impact = float(sys.argv[2]) if len(sys.argv) > 2 else 3.0

    print(f"40Ca+40Ca @ 8 MeV/u, b={impact} fm, {n_events} events per switch time")
    for switch in SWITCH_TIMES_FM_C:
        t0 = time.time()
        heavy_a: list[int] = []
        heavy_z: list[int] = []
        for event_index in range(n_events):
            event = run_imqmd_event(
                projectile_z=20,
                projectile_a=40,
                target_z=20,
                target_a=40,
                energy_per_a=8.0,
                impact_parameter=impact,
                event_index=event_index,
                time_fm_c=switch,
                fragment_method="mst",
            )
            for frag in event.primary_fragments:
                if frag.final_Z > 2:
                    heavy_a.append(frag.final_A)
                    heavy_z.append(frag.final_Z)
        if heavy_a:
            print(
                f"  t_switch={switch:5.0f} fm/c: <A>={np.mean(heavy_a):6.1f} "
                f"<Z>={np.mean(heavy_z):5.1f} n_heavy={len(heavy_a)} "
                f"({time.time() - t0:.0f}s)"
            )
        else:
            print(f"  t_switch={switch:5.0f} fm/c: no heavy fragments ({time.time() - t0:.0f}s)")


if __name__ == "__main__":
    main()
