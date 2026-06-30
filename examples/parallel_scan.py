#!/usr/bin/env python3
"""Parallel ImQMD scan on nucimp — Kr+Ni + U+U."""
import sys, os, time, json, multiprocessing as mp
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from mnt_sim.imqmd import (
    GaussianPacket, ImQMDNucleus, initialize_and_relax,
    propagate, reaction_fragments, weisskopf_evaporation,
)

M_N = 938.9

PROJECTILE_TARGET_Z = {
    "KrNi": (36, 28),
    "UU": (92, 92),
}
EXCITATION_ENERGY_NOTE = (
    "E* is the raw ImQMD fragment estimate and is known to run high for "
    "hot/deformed residues; use Z/A transfer observables for robust scan triage."
)

def make_event(Zp, Ap, Zt, At, E_per_A, b, seed, separation=28.0):
    """Build two-nucleus collision event."""
    proj = initialize_and_relax(Zp, Ap, sigma_r=1.1, seed=seed)
    targ = initialize_and_relax(Zt, At, sigma_r=1.1, seed=seed+1)
    p_beam = np.sqrt(2.0 * M_N * E_per_A)
    packets = []
    ref_pos = []
    ref_grp = []

    for pkt in proj.packets:
        shift = np.array([-0.5*separation, 0.5*b, 0.0])
        packets.append(GaussianPacket(pkt.r_i+shift, pkt.p_i+np.array([p_beam,0,0]), pkt.sigma_r, pkt.is_proton))
        ref_pos.append(pkt.r_i+shift)
        ref_grp.append(0)
    for pkt in targ.packets:
        shift = np.array([0.5*separation, -0.5*b, 0.0])
        packets.append(GaussianPacket(pkt.r_i+shift, pkt.p_i.copy(), pkt.sigma_r, pkt.is_proton))
        ref_pos.append(pkt.r_i+shift)
        ref_grp.append(1)

    nuc = ImQMDNucleus(Zp+Zt, Ap+At-Zp-Zt, packets, edf=proj.edf, reference_positions=np.asarray(ref_pos))
    nuc.reference_group_ids = np.asarray(ref_grp, dtype=int)
    nuc._collision_rng = np.random.default_rng(seed + 100_000)
    return nuc

def _print_progress(label, b, step, n_steps, t0):
    print(
        f"{label} b={b:g} step {step}/{n_steps} elapsed={time.time()-t0:.0f}s",
        flush=True,
    )


def _propagate_with_progress(nuc, label, b, n_steps, dt, progress_every):
    completed = 0
    t0 = time.time()
    while completed < n_steps:
        chunk = min(progress_every, n_steps - completed)
        propagate(nuc, dt=dt, n_steps=chunk, sample_every=chunk,
                  with_collisions=True, collision_dt=dt, remove_cm_drift=False,
                  use_surface_term=True, use_static_stabilizer=False)
        completed += chunk
        _print_progress(label, b, completed, n_steps, t0)
    return completed


def _transfer_happened(label, heavy):
    projectile_z, target_z = PROJECTILE_TARGET_Z.get(label.split("_", 1)[0], (None, None))
    if projectile_z is None:
        return False
    return any(fragment.Z not in (projectile_z, target_z) for fragment in heavy)


def _fragment_summary(fragments, limit=3):
    return ', '.join(
        f"Z={fragment.Z} A={fragment.A} E*={fragment.excitation_energy:.0f}MeV"
        for fragment in fragments[:limit]
    ) or "none"


def run_one(args):
    """Single event runner for multiprocessing."""
    label, Zp, Ap, Zt, At, E_per_A, b, n_steps, dt, progress_every = args
    t0 = time.time()
    seed = int(7000 + b*10 + Zp*100 + hash(label) % 10000)
    print(
        f"{label} b={b:g} starting: E/A={E_per_A:g} MeV, steps={n_steps}, dt={dt:g} fm/c",
        flush=True,
    )
    nuc = make_event(Zp, Ap, Zt, At, E_per_A, b, seed)
    completed = _propagate_with_progress(nuc, label, b, n_steps, dt, progress_every)
    frags = reaction_fragments(nuc)
    heavy = [f for f in frags if f.A > 10]
    cs = getattr(nuc, 'collision_stats', {})
    runtime_s = time.time() - t0
    transfer = _transfer_happened(label, heavy)
    print(
        f"{label} b={b:g} complete: steps={completed}, runtime={runtime_s:.1f}s, "
        f"collisions attempted={cs.get('attempted',0)} blocked={cs.get('blocked',0)} "
        f"accepted={cs.get('accepted',0)}, transfer={transfer}",
        flush=True,
    )
    print(f"{label} b={b:g} top heavy fragments: {_fragment_summary(heavy)}", flush=True)
    return {
        'label': label, 'b': b,
        'steps': completed, 'dt': dt,
        'runtime_s': runtime_s,
        'heavy': [(f.Z, f.A, float(f.excitation_energy), float(f.position[0])) for f in heavy[:4]],
        'transfer': transfer,
        'coll_att': cs.get('attempted',0), 'coll_acc': cs.get('accepted',0),
        'coll_blk': cs.get('blocked',0),
    }

def main():
    try:
        sys.stdout.reconfigure(line_buffering=True)
    except AttributeError:
        pass

    tasks = []
    # Kr+Ni scan
    for b in [2,3,4,5,6,7]:
        tasks.append((f'KrNi_b{b:.0f}', 36,86, 28,64, 25.0, b, 1500, 1.0, 100))
    # U+U at b=5,7
    for b in [5.0, 7.0]:
        tasks.append((f'UU_b{b:.0f}', 92,238, 92,238, 7.0, b, 2500, 2.0, 100))

    n_workers = min(16, len(tasks))
    print(f"Running {len(tasks)} tasks with {n_workers} workers", flush=True)
    print(EXCITATION_ENERGY_NOTE, flush=True)
    t0 = time.time()

    with mp.Pool(n_workers) as pool:
        results = pool.map(run_one, tasks)

    print(f"\nTotal wall time: {time.time()-t0:.0f}s", flush=True)
    print(f"{'Task':<12} {'b':>4} {'steps':>6} {'dt':>4} {'time':>7} {'coll_att':>9} {'coll_blk':>9} {'coll_acc':>9} {'transfer':>8} {'heavy_frags'}", flush=True)
    for r in sorted(results, key=lambda x: x['label']):
        h = r['heavy']
        h_str = ', '.join(f"Z={z}A={a}E*={e:.0f}" for z,a,e,_ in h[:3])
        print(
            f"{r['label']:<12} {r['b']:4.0f} {r['steps']:6d} {r['dt']:4.1f} "
            f"{r['runtime_s']:6.1f}s {r['coll_att']:9d} {r['coll_blk']:9d} "
            f"{r['coll_acc']:9d} {str(r['transfer']):>8} {h_str}",
            flush=True,
        )

    # Save
    out = {
        'results': results,
        'n_workers': n_workers,
        'excitation_energy_note': EXCITATION_ENERGY_NOTE,
    }
    with open('output/nucimp_scan.json', 'w') as f:
        json.dump(out, f, indent=2, default=str)
    print("\nSaved output/nucimp_scan.json", flush=True)

if __name__ == '__main__':
    main()
