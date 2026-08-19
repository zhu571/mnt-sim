"""Single U+U @ 7 MeV/A, b=4 fm sanity check after Fermi/EOS/cooling fixes."""

import pathlib
import sys
import time

import numpy as np

ROOT = pathlib.Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from mnt_sim.imqmd.reaction import run_imqmd_event

t0 = time.time()
time_fm_c = float(sys.argv[1]) if len(sys.argv) > 1 else 600.0
event = run_imqmd_event(
    projectile_z=92,
    projectile_a=238,
    target_z=92,
    target_a=238,
    energy_per_a=7.0,
    impact_parameter=float(sys.argv[2]) if len(sys.argv) > 2 else 4.0,
    event_index=0,
    time_fm_c=time_fm_c,
    dt=1.0,
    fragment_method="iso-mst",
)
frags = event.primary_fragments
heavy = [f for f in frags if f.final_Z > 2]
light = [f for f in frags if f.final_Z <= 2]
print(f"runtime: {time.time() - t0:.1f} s")
print(f"total fragments: {len(frags)}  heavy(Z>2): {len(heavy)}  light(Z<=2): {len(light)}")
print(f"collisions attempted/accepted: {event.attempted_collisions}/{event.accepted_collisions}")
heavy_sorted = sorted(heavy, key=lambda f: f.A, reverse=True)
for f in heavy_sorted[:8]:
    print(
        f"  Z={f.Z:3d} A={f.A:3d} E*={f.excitation_energy:8.1f} MeV -> "
        f"final Z={f.final_Z} A={f.final_A}  E_lab={f.e_lab:7.1f} theta={f.theta_lab:6.1f}"
    )
zs = sorted(f.Z for f in heavy)
print("heavy Z range:", zs[:5], "...", zs[-5:] if len(zs) > 5 else "")
amax = max((f.A for f in heavy), default=0)
print(f"heaviest heavy fragment A={amax}")
