# ImQMD Diagnosis and Solutions

This audit covers the seven core files in `mnt_sim/imqmd/`: `skyrme.py`, `nucleus.py`, `initializer.py`, `propagator.py`, `collisions.py`, `fragments.py`, and `decay.py`.

## P1: No nucleon transfer observed

### Root cause

The reaction dynamics still uses the artificial static ground-state stabilizer. `initialize_nucleus()` stores `reference_positions` for every initialized nucleus (`initializer.py:106`), and reaction builders copy those references into the combined system (`examples/run_uu_krni_benchmark.py:63-99`, `examples/verify_transfer.py:72-108`, `examples/reaction_event.py:33-62`). During every RK4 substep, `_derivatives()` adds `_static_reference_force()` (`propagator.py:102-106`). With `reference_group_ids` present, `_static_reference_force()` computes a group-centered harmonic restoring force to each nucleon's original projectile/target configuration (`propagator.py:78-99`).

That means a nucleon can scatter, but it remains tethered to its original nucleus. Accepted NN collisions in `attempt_nn_collision()` only rotate pair momenta (`collisions.py:352-370`); they do not change ownership or remove the restoring field. Fragment recognition then reads final spatial membership from MST (`fragments.py:161-165`), so if the tether prevents a centroid from migrating across the neck, the fragments keep the original `Z,A`.

There is a secondary setup issue for physical reaction kinematics: benchmark builders give projectile and target equal and opposite per-nucleon beam momenta (`examples/run_uu_krni_benchmark.py:66-88`, `examples/verify_transfer.py:75-97`). For asymmetric systems this is not the physical center-of-mass assignment unless momenta are mass-weighted; for fixed-target comparisons it is also not the lab setup.

### Evidence

The static term is not part of the ImQMD Hamiltonian described in `papers/imqmd_formulas.md:26-40`. It was added for static test stability (`skyrme.py:111-124`, `nucleus.py:108-109`) and calibrated by an energy offset (`initializer.py:108-111`). In a reaction, physics demands that nucleons move under the mutual EDF, Coulomb, and collision terms, not under a private spring to their initial sampled coordinates.

The current code path is:

1. `make_collision_event()` creates a combined `ImQMDNucleus` with `reference_positions`.
2. `propagate()` calls `_rk4_step()` (`propagator.py:166-169`).
3. `_rk4_step()` calls `_derivatives()` four times (`propagator.py:114-120`).
4. `_derivatives()` always includes `_static_reference_force()` (`propagator.py:102-106`).
5. `_static_reference_force()` keeps each original group internally near its initial shape (`propagator.py:91-99`).

This directly explains "collisions happen but nucleons do not move from one nucleus to the other."

### Fix

Separate static-nucleus stabilization from reaction propagation.

Specific changes:

1. Add a reaction mode or flag to `propagate()` and `_rk4_step()` such as `use_static_reference_force: bool = True`, defaulting to `True` only for existing static tests.
2. Pass that flag into `_derivatives()` and skip `_static_reference_force()` when propagating a two-body reaction.
3. In reaction builders, create the combined system with `reference_positions=None`, or set `system.disable_static_reference_force = True`, instead of passing projectile/target reference positions.
4. Keep `reference_positions` only for accepted isolated nuclei and static stability tests.
5. Fix asymmetric reaction kinematics in event builders: compute projectile/target momenta in the desired frame, mass-weight them so total center-of-mass momentum is correct, and avoid relying on `attempt_nn_collision()`'s global mean subtraction (`collisions.py:374`) as a hidden correction.

### Risk

High for existing static tests if the default behavior changes globally. Low if the new flag defaults to the current static behavior and reaction scripts explicitly disable it. Tests that currently assume a static spring, especially `tests/test_static.py:16-57`, may fail if run in pure ImQMD mode without improving initialization/EDF stability first.

## P2: Excitation energy `E*` unreliable

### Root cause

`compute_fragment_excitation()` compares two incompatible energy scales. It rebuilds a sub-nucleus with `energy_offset=0.0` and `reference_positions=None` (`fragments.py:258-269`), computes the raw compact-QMD energy (`fragments.py:270`), then subtracts a mass-table or liquid-drop ground-state energy (`fragments.py:287-288`, `decay.py:77-88`).

The raw compact EDF is not calibrated to physical nuclear masses. Calibration is stored as `nucleus.energy_offset` during initialization (`initializer.py:108-111`) and included in `energy_components()` (`nucleus.py:108-117`), but that offset is deliberately omitted in the fragment sub-nucleus (`fragments.py:262-269`). The result is a subtraction of "raw QMD Hamiltonian" minus "physical mass-table ground state."

The function also adds static deformation energy when the parent still has a static reference field (`fragments.py:271-286`). That folds the artificial stabilizer from P1 into `E*`.

### Evidence

For cold initialized nuclei, the sub-fragment raw QMD energy is often hundreds of MeV more bound than `-binding_energy(Z,A)`, so `max(..., 0.0)` (`fragments.py:288`) hides the scale error by clamping `E*` to zero. For hot or fused systems, deformation and kinetic changes can push the same uncalibrated difference positive, producing arbitrary hundreds-to-thousands of MeV values.

Physics demands the Zhao-style rest-frame subtraction stated in `papers/imqmd_formulas.md:326-341`, but both terms must be on the same Hamiltonian/mass scale. The current code uses the same formula shape but not a calibrated ground-state reference.

### Fix

Make `E*` a calibrated relative quantity.

Specific changes:

1. Add a ground-state energy provider for fragments, e.g. `ground_state_qmd_energy(Z, A, edf, sigma_r)` that initializes or loads a stable nucleus with the same EDF and returns its calibrated `total_energy()`.
2. In `compute_fragment_excitation()`, compute `E_int` with reaction-only physical terms: kinetic in the fragment rest frame plus EDF and Coulomb, excluding static reference/deformation energy.
3. Recenter fragment positions before building the sub-nucleus. Pair energies are translationally invariant, but recentering avoids future non-invariant terms and keeps diagnostics sane.
4. Subtract a calibrated ground-state energy on the same scale, not raw compact EDF against `decay.binding_energy()`.
5. Keep `decay.binding_energy()` for statistical evaporation thresholds and mass tables, not for raw QMD calibration unless an explicit conversion is introduced.
6. Remove or gate the static deformation addition (`fragments.py:271-286`) behind a static-test-only diagnostic path.

### Risk

Medium. `tests/test_fragments.py:64-68` only asserts non-negative `E*`, so it may still pass while values change. Any downstream evaporation tests or benchmark assertions that expect current high `E*` values will need updated tolerances. This fix is likely to expose rather than create existing physics issues.

## P3: U+U too slow

### Root cause

The neighbor list only reduces candidate NN collision pairs. The dominant propagation and blocking paths remain quadratic or worse for `A=476`.

Hot spots by code path:

1. `_rk4_step()` calls `_derivatives()` four times per step (`propagator.py:114-120`).
2. `_centroid_density_force()` builds full `A x A x 3` differences and full density weights every derivative call (`propagator.py:40-75`).
3. `_coulomb_forces()` evaluates all proton-proton pairs (`propagator.py:12-34`).
4. `attempt_nn_collision()` computes all centroid densities once per collision step through `nucleus.centroid_densities()` (`collisions.py:317-322`, `nucleus.py:93-98`), which itself calls the full-density evaluator (`nucleus.py:75-91`).
5. `_pauli_blocking_probability_fast()` still scans all same-species nucleons for each attempted pair (`collisions.py:162-179`, `collisions.py:194-218`).
6. `fermi_constraint_check()` loops over all identical pairs every 20 fm/c (`collisions.py:380-408`).

For U+U, `output/uu_benchmark.npz` shows neighbor lists around tens of thousands of candidate pairs per collision step and a timeout at `step=876 phase=collisions`, so the neighbor list is not enough.

### Evidence

A spatial neighbor list exists in `collisions.py:239-311`, but it is only used in `attempt_nn_collision()` (`collisions.py:322`). The force and density routines ignore it and recompute dense all-pairs arrays. Since RK4 uses four force evaluations per step, 2500 steps means 10,000 dense mean-field evaluations for 476 nucleons before counting collisions.

Physics does not demand all-pairs work for Gaussian terms at arbitrarily long range. The Gaussian nuclear density terms have finite practical support; Coulomb can use optimized pair kernels or approximate far fields.

### Fix

Optimize the shared pair infrastructure, not only collision candidates.

Specific changes:

1. Refactor `_centroid_density_force()` to use a cutoff neighbor list for Gaussian nuclear terms. A cutoff of several `sigma_r` should be validated against energy conservation and static stability.
2. Cache pair displacements, distances, Gaussian weights, and densities once per RK4 derivative call and reuse them for force components.
3. Replace `_coulomb_forces()` with a vectorized neighbor/far-field split, a Barnes-Hut-style approximation for heavy systems, or a compiled kernel. At minimum, avoid repeated Python/scipy overhead inside every derivative.
4. Refactor `_pauli_blocking_probability_fast()` to inspect local phase-space neighbors only, using the same spatial cells and a momentum cutoff from `sigma_p = hbar/(2 sigma_r)` (`collisions.py:167`).
5. Make `fermi_constraint_check()` cell-based rather than all identical pairs (`collisions.py:388-389`).
6. Consider switching reaction propagation from RK4 to a symplectic velocity-Verlet/leapfrog integrator after validating energy conservation; it would reduce force evaluations from four to one or two per step.

### Risk

Medium to high. Performance changes can alter numerical trajectories and collision ordering. Existing tests should still pass if tolerances are physical, but static energy conservation (`tests/test_static.py:39-44`) and fragment separation (`tests/test_fragments.py:89-131`) need to be rerun with stricter diagnostic logging. Cutoff choices can bias transfer yields if too aggressive.

## P4: Kr+Ni `b=2-4` all fuse or lacks a deep-inelastic/transfer window

### Root cause

The code has no reliable reaction-stage classification or reseparation stopping condition. `run_krni_event()` propagates a fixed 600 steps (`examples/run_uu_krni_benchmark.py:271`) and immediately calls `reaction_fragments()` (`examples/run_uu_krni_benchmark.py:272`). `reaction_fragments()` uses `adaptive_mst()` with coordinate-only membership (`fragments.py:161-165`), and the default call sets `p_cut=None` (`fragments.py:161-165`). `adaptive_mst()` can connect a whole necked system into one component through any chain of links below the adaptive cutoff (`fragments.py:146-158`).

Combined with P1, the dynamics is biased toward two extremes: nucleons are held to original residues when the system separates, or MST sees a connected composite and labels it fused. There is no explicit deep-inelastic category based on contact time, reseparation, two-heavy-fragment recognition, total kinetic energy loss, or transferred nucleon counts.

### Evidence

Zhao-style production runs terminate after reseparation and then recognize primary fragments, not after a fixed arbitrary step count. The paper text in the local PDF states that simulations are terminated at `1000 fm/c` after reseparation, with fragments recognized then. The current Kr+Ni scan stops at `600 fm/c` regardless of whether the composite has reseparated.

The fragment code also defaults to coordinate-only reaction fragments (`fragments.py:161-165`). Momentum information exists in `minimum_spanning_tree()` and `adaptive_mst()` (`fragments.py:85`, `fragments.py:134`, `fragments.py:154-155`), but the reaction convenience path does not use it.

### Fix

Add reaction outcome tracking before changing physics parameters.

Specific changes:

1. Implement a reseparation monitor in propagation or event scripts: track the two largest clusters or projectile-like/target-like centroids, contact start/end, and require a post-reseparation wait before final MST.
2. Replace fixed `n_steps=600` in `run_krni_event()` with `run_until_reseparated(max_time, settle_time)` for transfer/deep-inelastic events.
3. In `reaction_fragments()`, use a momentum cutoff for final primary fragments or expose it in benchmark scripts, e.g. `reaction_fragments(system, p_cut=250.0)`, after validating against static nuclei.
4. Disable the static reference force during reactions as described in P1.
5. Add event classification fields: `fused`, `deep_inelastic`, `quasi_elastic`, `transfer_counts`, `contact_time`, and `final_TKE`.
6. Use physically correct Kr+Ni center-of-mass or lab kinematics in event construction.

### Risk

Medium. Existing fragment tests include coordinate-only assumptions (`tests/test_fragments.py:43-50`, `tests/test_fragments.py:80-87`). Keep `minimum_spanning_tree()` defaults unchanged and update only the reaction convenience path or benchmark scripts to limit test fallout.

## P5: Not compared with Zhao 2016

### Root cause

The current scripts produce a few event summaries, not Zhao-comparable production cross sections. `examples/run_uu_krni_benchmark.py` runs one U+U event at `b=5 fm` (`examples/run_uu_krni_benchmark.py:210-255`) and six Kr+Ni impact parameters (`examples/run_uu_krni_benchmark.py:290-328`). It saves fragment lists, but it does not accumulate `N_frag(Z,A,b,E*) / N_tot(b)` with `2*pi*b*db` weights.

The local Zhao paper text specifies `bmax = 15 fm`, `Delta b = 0.15 fm`, many events per impact parameter, primary fragments at reseparation plus de-excitation, and angle-binned yields. The current decay module is also not Zhao/HIVAP-equivalent: it implements neutron evaporation and a simple fission probability (`decay.py:118-179`), but not gamma, proton, alpha, and full fission-channel competition.

### Evidence

Zhao's cross-section formula is also summarized in the local PDF text: primary fragment cross sections require summing over impact parameters with event counts and `2*pi*b*Delta b` weights. The code has no accumulator for that observable. `fragment_record()` stores individual fragment fields (`examples/run_uu_krni_benchmark.py:162-169`), and `save_json_npz()` stores summaries (`examples/run_uu_krni_benchmark.py:194-197`), but no yield grid.

### Fix

Create a dedicated Zhao comparison pipeline.

Specific changes:

1. Add a scan driver for `238U+238U` at `7.0 MeV/A` with `b` from `0` to `15 fm` and configurable `db` (`0.15 fm` for final validation, coarser for development).
2. For each `b`, run many seeds/orientations and store primary fragments after reseparation.
3. Accumulate `sigma_primary[Z,A,E_bin,theta_bin] += 2*pi*b*db / N_tot(b)` for every recognized primary fragment.
4. Add residual-fragment accumulation after decay. Either couple to HIVAP/GEMINI or clearly label the compact `decay.py` output as not Zhao-equivalent.
5. Produce comparison plots matching Zhao Fig. 1/2/4/5 shapes: primary/residual `Z,A` maps, angle-binned primary yields, uranium isotope angular distributions, and average `E*`/contact-time maps.
6. Add a small regression fixture that runs a coarse low-stat scan and verifies nonzero neutron-transfer uranium-like channels such as `Z=92, A>=244` before attempting high-stat comparison.

### Risk

Low for existing unit tests because this is mostly new analysis code. High compute cost. Physics results will remain non-comparable until P1, P2, and P4 are fixed.

## Priority order

1. **P1: Disable static reference force in reactions.** This is the main blocker for nucleon exchange and contaminates fragment excitation through static deformation energy.
2. **P2: Recalibrate `E*`.** Transfer yields are not useful if the decay stage receives arbitrary excitation energies.
3. **P4: Add reseparation/outcome classification.** This turns "fused vs elastic" into measurable fusion, deep-inelastic, transfer, and quasi-elastic categories.
4. **P3: Optimize heavy-system propagation/collisions.** Do this after the physics path is corrected so optimizations preserve the right behavior.
5. **P5: Build Zhao comparison pipeline.** Run serious yield comparisons only after the event physics and `E*` are credible.

## Expected outcome after fixes

For `86Kr+64Ni` at `25 MeV/A`, a successful scan should show a transition with impact parameter: central events dominated by fusion or long-contact composites, intermediate events with two heavy fragments and nonzero transfer (`Z` and `A` outside `{36,86}` and `{28,64}`), and peripheral events dominated by quasi-elastic original-like fragments. Excitation energies should be tens of MeV for ordinary transfer residues, not hundreds or thousands of MeV by construction.

For `238U+238U` at `7 MeV/A`, successful events should produce uranium-like fragments with neutron pickup/loss around the original `Z=92,A=238`, plus lower-probability transuranium primary fragments. The final validation target should be a weighted `Z,A` yield landscape and angle-binned distributions comparable in structure to Zhao 2016, not single-event fragment lists.
