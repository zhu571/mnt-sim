#!/usr/bin/env python3
"""P1-fix verification on nucimp — 8 events, check for nucleon transfer."""
import sys, os, time, json, multiprocessing as mp
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from mnt_sim.imqmd import (
    GaussianPacket, ImQMDNucleus, initialize_and_relax,
    propagate, reaction_fragments,
)

M_N = 938.9

def make_event(Zp, Ap, Zt, At, E_per_A, b, seed, separation=28.0):
    proj = initialize_and_relax(Zp, Ap, sigma_r=1.1, seed=seed)
    targ = initialize_and_relax(Zt, At, sigma_r=1.1, seed=seed+1)
    p_beam = np.sqrt(2.0 * M_N * E_per_A)
    packets, ref_pos, ref_grp = [], [], []
    for pkt in proj.packets:
        shift = np.array([-0.5*separation, 0.5*b, 0.0])
        packets.append(GaussianPacket(pkt.r_i+shift, pkt.p_i+np.array([p_beam,0,0]), pkt.sigma_r, pkt.is_proton))
        ref_pos.append(pkt.r_i+shift); ref_grp.append(0)
    for pkt in targ.packets:
        shift = np.array([0.5*separation, -0.5*b, 0.0])
        packets.append(GaussianPacket(pkt.r_i+shift, pkt.p_i.copy(), pkt.sigma_r, pkt.is_proton))
        ref_pos.append(pkt.r_i+shift); ref_grp.append(1)
    nuc = ImQMDNucleus(Zp+Zt, Ap+At-Zp-Zt, packets, edf=proj.edf, reference_positions=np.asarray(ref_pos))
    nuc.reference_group_ids = np.asarray(ref_grp, dtype=int)
    nuc._collision_rng = np.random.default_rng(seed + 100_000)
    return nuc

def run_one(args):
    label, Zp, Ap, Zt, At, E_per_A, b, n_steps, dt = args
    t0 = time.time()
    seed = int(7000 + b*10 + Zp*100 + hash(label) % 10000)
    nuc = make_event(Zp, Ap, Zt, At, E_per_A, b, seed)
    propagate(nuc, dt=dt, n_steps=n_steps, sample_every=max(1,n_steps//5),
              with_collisions=True, collision_dt=1.0,
              remove_cm_drift=False, use_surface_term=True, use_static_stabilizer=False)
    frags = reaction_fragments(nuc)
    heavy = [(f.Z, f.A, float(f.excitation_energy)) for f in frags if f.A > 4]
    cs = getattr(nuc, 'collision_stats', {})
    z_vals = sorted(set(f[0] for f in heavy))
    has_transfer = any(z not in {Zp, Zt} for z in z_vals)
    print(f"  {label}: {len(heavy)} heavy frags, Z={z_vals}, transfer={has_transfer}, {time.time()-t0:.1f}s", flush=True)
    return {'label': label, 'b': b, 'runtime_s': time.time()-t0,
            'coll_acc': cs.get('accepted',0), 'coll_att': cs.get('attempted',0),
            'coll_blk': cs.get('blocked',0), 'transfer': has_transfer,
            'heavy_Z': z_vals, 'heavy': heavy[:5]}

def main():
    tasks = [
        ('O+Ca b4', 8,16, 20,40, 8.0, 4.0, 400, 1.0),
        ('Ca+Ca b5', 20,40, 20,48, 10.0, 5.0, 600, 1.0),
        ('Kr+Ni b3', 36,86, 28,64, 25.0, 3.0, 800, 1.0),
        ('Kr+Ni b5', 36,86, 28,64, 25.0, 5.0, 800, 1.0),
        ('UU b5', 92,238, 92,238, 7.0, 5.0, 1500, 2.0),
        ('UU b7', 92,238, 92,238, 7.0, 7.0, 1500, 2.0),
    ]
    n_workers = min(len(tasks), 12)
    print(f"P1-fix verification: {len(tasks)} tasks, {n_workers} workers", flush=True)
    t0 = time.time()
    with mp.Pool(n_workers) as pool:
        results = pool.map(run_one, tasks)
    results.sort(key=lambda r: r['label'])
    print(f"\nTotal wall: {time.time()-t0:.0f}s", flush=True)
    for r in results:
        print(f"{r['label']:<12} b={r['b']:.0f} {r['runtime_s']:.0f}s acc={r['coll_acc']} transfer={r['transfer']} Z={r['heavy_Z']}", flush=True)
    with open('output/p1_fix_test.json', 'w') as f:
        json.dump(results, f, indent=2, default=str)
    print("Saved output/p1_fix_test.json", flush=True)

if __name__ == '__main__':
    main()
