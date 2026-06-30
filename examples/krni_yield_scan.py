#!/usr/bin/env python3
"""Kr+Ni b-scan for yield data — runs on laptop."""
import sys, os, time, json
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from mnt_sim.imqmd import (
    GaussianPacket, ImQMDNucleus, initialize_and_relax,
    propagate, reaction_fragments, weisskopf_evaporation,
)

M_N = 938.9
B_VALUES = [2, 3, 4, 5, 6, 7]
N_STEPS = 800
DT = 1.0

def make_event(Zp, Ap, Zt, At, E_per_A, b, seed):
    proj = initialize_and_relax(Zp, Ap, sigma_r=1.1, seed=seed)
    targ = initialize_and_relax(Zt, At, sigma_r=1.1, seed=seed+1)
    p_beam = np.sqrt(2.0 * M_N * E_per_A)
    sep = 28.0
    packets, ref_pos, ref_grp = [], [], []
    for pkt in proj.packets:
        shift = np.array([-0.5*sep, 0.5*b, 0.0])
        packets.append(GaussianPacket(pkt.r_i+shift, pkt.p_i+np.array([p_beam,0,0]), pkt.sigma_r, pkt.is_proton))
        ref_pos.append(pkt.r_i+shift); ref_grp.append(0)
    for pkt in targ.packets:
        shift = np.array([0.5*sep, -0.5*b, 0.0])
        packets.append(GaussianPacket(pkt.r_i+shift, pkt.p_i.copy(), pkt.sigma_r, pkt.is_proton))
        ref_pos.append(pkt.r_i+shift); ref_grp.append(1)
    nuc = ImQMDNucleus(Zp+Zt, Ap+At-Zp-Zt, packets, edf=proj.edf, reference_positions=np.asarray(ref_pos))
    nuc.reference_group_ids = np.asarray(ref_grp, dtype=int)
    nuc._collision_rng = np.random.default_rng(seed + 100_000)
    return nuc

def run_one(b):
    seed = int(7000 + b*10 + 36*100 + 28*1000)
    t0 = time.time()
    nuc = make_event(36, 86, 28, 64, 25.0, b, seed)
    propagate(nuc, dt=DT, n_steps=400, sample_every=200,
              with_collisions=True, collision_dt=1.0,
              remove_cm_drift=False, use_surface_term=True, use_static_stabilizer=False)
    frags = reaction_fragments(nuc)
    cs = getattr(nuc, 'collision_stats', {})
    
    primary = []
    final = []
    for f in frags:
        if f.A <= 1:
            continue
        theta = np.degrees(np.arccos(np.clip(f.momentum[0] / max(np.linalg.norm(f.momentum), 1e-12), -1, 1)))
        e_lab = np.dot(f.momentum, f.momentum) / (2*M_N*max(f.A,1))
        primary.append({'Z': int(f.Z), 'A': int(f.A), 'E*': float(f.excitation_energy), 
                        'theta': float(theta), 'E_lab': float(e_lab)})
        # Evaporation
        if f.A > 10 and f.excitation_energy > 5.0:
            ch = weisskopf_evaporation(f.Z, f.A, f.excitation_energy, n_max=8)
            probs = [c[2] for c in ch]
            idx = np.random.default_rng(seed+len(final)).choice(len(ch), p=np.array(probs)/sum(probs))
            z_f, a_f, _ = ch[idx]
        else:
            z_f, a_f = f.Z, f.A
        final.append({'Z': int(z_f), 'A': int(a_f), 'Z_primary': int(f.Z), 'A_primary': int(f.A)})
    
    runtime = time.time() - t0
    print(f"  b={b:.0f}: {len(primary)} frags, {cs.get('accepted',0)} coll, {runtime:.0f}s, Z={sorted(set(p['Z'] for p in primary))}", flush=True)
    return {'b': b, 'runtime': runtime, 'primary': primary, 'final': final,
            'coll_acc': cs.get('accepted',0), 'coll_att': cs.get('attempted',0)}

def main():
    print(f"Kr+Ni b-scan: E/A=25 MeV, {N_STEPS} steps, dt={DT}")
    print(f"b = {B_VALUES}\n")
    results = []
    t0 = time.time()
    for b in B_VALUES:
        results.append(run_one(b))
    
    # Summary
    print(f"\nTotal: {time.time()-t0:.0f}s")
    print(f"{'b':>3} {'frags':>5} {'acc_col':>7} {'Z_range'}")
    for r in results:
        zs = sorted(set(p['Z'] for p in r['primary']))
        print(f"{r['b']:3.0f} {len(r['primary']):5d} {r['coll_acc']:7d} {zs}")
    
    # Yield data
    print(f"\n=== Primary Z yields (heavy only, A>10) ===")
    for r in results:
        heavy = [p for p in r['primary'] if p['A'] > 10]
        if heavy:
            zs = [h['Z'] for h in heavy]
            As = [h['A'] for h in heavy]
            print(f"  b={r['b']:.0f}: Z={zs}, A={As}")
    
    print(f"\n=== After evaporation (A>10) ===")
    for r in results:
        heavy_f = [f for f in r['final'] if f['A'] > 10]
        if heavy_f:
            print(f"  b={r['b']:.0f}: {[(f['Z'], f['A']) for f in heavy_f]}")
    
    os.makedirs('output', exist_ok=True)
    with open('output/krni_yield.json', 'w') as f:
        json.dump(results, f, indent=2, default=str)
    print("Saved output/krni_yield.json")

if __name__ == '__main__':
    main()
