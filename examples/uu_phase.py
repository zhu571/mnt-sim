#!/usr/bin/env python3
"""U+U MNT scan — two-phase: spring ON approach, spring OFF collision."""
import sys, os, time, json, multiprocessing as mp
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from mnt_sim.imqmd import (
    initialize_and_relax, GaussianPacket, ImQMDNucleus,
    propagate, reaction_fragments,
)

M_N, B_VALUES = 938.9, [5, 9, 13, 15]

def run_one(b):
    seed = int(9200 + b*10)
    t0 = time.time()
    proj = initialize_and_relax(92, 238, sigma_r=1.1, seed=seed)
    targ = initialize_and_relax(92, 238, sigma_r=1.1, seed=seed+1)
    p_beam = np.sqrt(2.0*M_N*7.0); sep = 40.0
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
    nuc._collision_rng = np.random.Generator(np.random.PCG64(seed+100000))
    nuc._grid_eta = 0.5

    # Phase 1: approach with spring stabilizer (200 steps, dt=2.0)
    propagate(nuc, dt=2.0, n_steps=200, sample_every=200,
              with_collisions=False,
              remove_cm_drift=False, use_surface_term=True, use_static_stabilizer=True)
    
    # Phase 2: collision without spring (400 steps, dt=2.0)
    propagate(nuc, dt=2.0, n_steps=400, sample_every=200,
              with_collisions=True, collision_dt=1.0,
              remove_cm_drift=False, use_surface_term=True, use_static_stabilizer=False)
    
    frags = reaction_fragments(nuc)
    cs = getattr(nuc, 'collision_stats', {})
    heavy = [(f.Z, f.A, f.excitation_energy) for f in frags if f.A > 10]
    zs = sorted(set(z for z,_,_ in heavy))
    transfer = any(z != 92 for z in zs)
    elapsed = time.time() - t0
    print(f"  b={b:.0f}: {len(heavy)} heavy, Z={zs}, transfer={transfer}, acc={cs.get('accepted',0)}, {elapsed:.0f}s", flush=True)
    return {'b': b, 'heavy': heavy, 'transfer': transfer, 'elapsed': elapsed}

def main():
    print(f"U+U two-phase: 200(spring)+400(free), b={B_VALUES}\n")
    t0 = time.time()
    with mp.Pool(min(len(B_VALUES), 8)) as pool:
        results = pool.map(run_one, B_VALUES)
    results.sort(key=lambda r: r['b'])
    print(f"\nTotal: {time.time()-t0:.0f}s")
    for r in results:
        h = r['heavy']
        h_str = ', '.join(f"Z={z}A={a}" for z,a,_ in h[:3]) if h else 'NONE'
        print(f"b={r['b']:3.0f}: {len(h)} heavy, transfer={r['transfer']}, {r['elapsed']:.0f}s, {h_str}")
    os.makedirs('output', exist_ok=True)
    with open('output/uu_phase.json', 'w') as f: json.dump(results, f, indent=2)

if __name__ == '__main__':
    main()
