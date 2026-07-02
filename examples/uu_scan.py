#!/usr/bin/env python3
"""U+U MNT scan — nucimp parallel."""
import sys, os, time, json, multiprocessing as mp
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from mnt_sim.imqmd import (
    initialize_and_relax, GaussianPacket, ImQMDNucleus,
    propagate, reaction_fragments,
)

M_N = 938.9
B_VALUES = [5, 7, 9, 11, 13]
N_STEPS = 800
DT = 2.0

def make_uu_event(b, seed):
    proj = initialize_and_relax(92, 238, sigma_r=1.1, seed=seed)
    targ = initialize_and_relax(92, 238, sigma_r=1.1, seed=seed+1)
    p_beam = np.sqrt(2.0 * M_N * 7.0)
    sep = 40.0  # Zhao 2016 uses d0=40 fm
    packets, ref_pos, ref_grp = [], [], []
    for pkt in proj.packets:
        shift = np.array([-0.5*sep, 0.5*b, 0.0])
        packets.append(GaussianPacket(pkt.r_i+shift, pkt.p_i+np.array([p_beam,0,0]), pkt.sigma_r, pkt.is_proton))
        ref_pos.append(pkt.r_i+shift); ref_grp.append(0)
    for pkt in targ.packets:
        shift = np.array([0.5*sep, -0.5*b, 0.0])
        packets.append(GaussianPacket(pkt.r_i+shift, pkt.p_i.copy(), pkt.sigma_r, pkt.is_proton))
        ref_pos.append(pkt.r_i+shift); ref_grp.append(1)
    nuc = ImQMDNucleus(184, 292, packets, edf=proj.edf, reference_positions=np.asarray(ref_pos))
    nuc.reference_group_ids = np.asarray(ref_grp, dtype=int)
    nuc._collision_rng = np.random.Generator(np.random.PCG64(seed + 100_000))
    return nuc

def run_one(b):
    seed = int(9200 + b*10)
    t0 = time.time()
    nuc = make_uu_event(b, seed)
    propagate(nuc, dt=DT, n_steps=N_STEPS, sample_every=500,
              with_collisions=True, collision_dt=1.0,
              remove_cm_drift=False, use_surface_term=True, use_static_stabilizer=False)
    frags = reaction_fragments(nuc)
    cs = getattr(nuc, 'collision_stats', {})
    heavy = [(f.Z, f.A, f.excitation_energy, float(np.linalg.norm(f.position))) for f in frags if f.A > 10]
    zs = sorted(set(z for z, _, _, _ in heavy))
    transfer = any(z != 92 for z in zs)
    elapsed = time.time() - t0
    print(f"  U+U b={b:.0f}: {len(heavy)} heavy, Z={zs}, transfer={transfer}, "
          f"acc={cs.get('accepted',0)}, {elapsed:.0f}s", flush=True)
    return {'b': b, 'heavy': heavy[:10], 'transfer': transfer, 'elapsed': elapsed,
            'coll_acc': cs.get('accepted',0), 'coll_att': cs.get('attempted',0),
            'coll_blk': cs.get('blocked',0)}

def main():
    print(f"U+U scan: E/A=7 MeV, {N_STEPS} steps, dt={DT}, b={B_VALUES}\n")
    n_workers = min(len(B_VALUES), 16)
    t0 = time.time()
    with mp.Pool(n_workers) as pool:
        results = pool.map(run_one, B_VALUES)
    results.sort(key=lambda r: r['b'])
    print(f"\nTotal wall: {time.time()-t0:.0f}s")
    for r in results:
        h = r['heavy']
        h_str = ', '.join(f"Z={z}A={a}" for z,a,_,_ in h[:3])
        print(f"b={r['b']:3.0f}: {len(h)} heavy, transfer={r['transfer']}, "
              f"acc={r['coll_acc']}, {r['elapsed']:.0f}s, {h_str}")
    os.makedirs('output', exist_ok=True)
    with open('output/uu_scan.json', 'w') as f:
        json.dump(results, f, indent=2, default=str)
    print("Saved output/uu_scan.json")

if __name__ == '__main__':
    main()
