# ImQMD Code Audit — Comparison with Published Papers

**Code base:** `mnt_sim/imqmd/` (Phase 1)  
**Date:** 2026-08-11  
**Summary:** Each of 7 audit items is rated ✅ (matches published methods), ⚠️ (minor issues or needs verification), or ❌ (likely incorrect). Line numbers refer to the files as of this audit.

---

## 1. Skyrme Parameter Sets

**Rating: ⚠️**

### 1.1 IQ1, IQ2, IQ3 — Match Published Values

The code at `skyrme.py:39-41` defines IQ1, IQ2, IQ3 with parameters that match Wang+2002 (PRC 65, 064606) and Li+2013 (PRC 87, 064615):

| Parameter | Code IQ1 | Code IQ2 | Code IQ3 | Wang+2002 IQ3 | Status |
|-----------|----------|----------|----------|---------------|--------|
| alpha | -310 | -356 | -207 | -207 | ✅ |
| beta | 258 | 303 | 138 | 138 | ✅ |
| gamma | 7/6 | 7/6 | 7/6 | 7/6 | ✅ |
| gsur | 19.8 | 7.0 | 18.0 | 18.0 | ✅ |
| g_tau | 9.5 | 12.5 | 14.0 | 14.0 | ✅ |
| eta | 2/3 | 2/3 | 5/3 | 5/3 | ✅ |
| c_sym | 32.0 | 32.0 | 32.0 | 32.0 | ✅ |
| kappa_s | 0.08 | 0.08 | 0.08 | 0.08 | ✅ |
| rho0 | 0.165 | 0.165 | 0.165 | 0.165 | ✅ |
| sigma0 | 0.49 | 0.88 | 0.94 | 0.94 | ✅ |
| sigma1 | 0.16 | 0.09 | 0.018 | 0.018 | ✅ |

### 1.2 IQ3A and IQ3B — Undocumented Variants ⚠️

At `skyrme.py:42-43`, IQ3A and IQ3B differ from the published IQ3 in four parameters:

| Parameter | IQ3 (published) | IQ3A (code) | IQ3B (code) |
|-----------|-----------------|-------------|-------------|
| gsur | 18.0 | **16.5** | 18.0 |
| c_sym | 32.0 | **34.0** | 34.0 |
| kappa_s | 0.08 | **0.4** | **0.6** |
| sigma1 | 0.018 | **0.020** | 0.018 |

**Impact of IQ3A deviations:**

- **gsur=16.5** (down from 18.0): ~8% weaker surface tension. Reduces nuclear binding at the surface, slightly softer EOS.
- **c_sym=34.0** (up from 32.0): ~6% stronger volume symmetry energy. Impacts neutron skin and isospin diffusion observables.
- **kappa_s=0.4** (up 5× from 0.08): **Dramatic increase** in surface-symmetry energy. The κ_s term suppresses surface isospin asymmetry and affects neutron-rich fragment production. This is a qualitatively different surface-symmetry behavior compared to IQ3.
- **sigma1=0.020** (up from 0.018): ~11% larger wave-packet width growth with A.

⚠️ **Recommendation:** These variants are not cited in the code docstrings against any specific publication. The kappa_s increase from 0.08 to 0.4 is a 5× change and represents a fundamentally different isovector surface functional. The parameter set should either (a) cite a paper where IQ3A/IQ3B were introduced, or (b) document the physical motivation for each change explicitly. A brief note in the `SkyrmeParameters` docstring connecting IQ3A/IQ3B to specific tuning objectives (e.g. "IQ3A: tuned for near-barrier U+U MNT with enhanced surface-symmetry") would suffice.

### 1.3 Surface-Symmetry Functional Form ✅

At `grid_edf.py:203-218`, the surface-symmetry integrand is:

```
E_surf_sym = C_s/(2ρ₀) × (-κ_s × |∇ρ|²) × δ²
```

where δ = (ρₙ − ρₚ)/ρ. This matches the published ImQMD Skyrme-like surface-symmetry form (Wang+2014, PRC 89, Eq. 5; Li+2013, PRC 87, Eq. 2; Zhang+2020 review, Int. J. Mod. Phys. E 29, 2030008, Eq. 29):

```
ℋ_surf_sym = −C_s·κ_s/(2ρ₀) · (∇ρ)² · ((ρₙ−ρₚ)/ρ)²
```

✅ **The gradient derivatives in `_density_functional_derivatives`** (`grid_edf.py:439-521`) are handled correctly: the code uses the exact transpose Dᵀ of the discrete np.gradient stencil (not a second application of np.gradient = -laplacian), which ensures energy-force self-consistency. This is a subtle but important fidelity point that many codes get wrong at boundary stencils. The analytical forces also include the ds/dρ terms from the density-dependent calibration scale (`_effective_scale`), which was previously missing (as noted in the docstring at `grid_edf.py:450-457`).

### 1.4 SKP\* ✅

At `skyrme.py:44`, SKP\* parameters match the standard values from the SkP\* Skyrme force (Dobaczewski+1984, NPA 422, 103), adapted to the ImQMD form factor. Not used as the default, but correctly parameterized.

---

## 2. Pauli Blocking

**Rating: ⚠️**

### 2.1 Occupation Formula and OCCUPATION_SOFTENING ⚠️

At `collisions.py:137-171` (`_occupation_wigner`) and lines 161-170, the phase-space occupation is computed as:

```
occupation = 4 × Σ exp(-Δr²/2σᵣ² - Δp²/2σₚ²) / OCCUPATION_SOFTENING
```

where `OCCUPATION_SOFTENING = 10.0` (line 26).

Zhang 2020 review (Int. J. Mod. Phys. E 29, 2030008), Eq. (63):

```
f_i = 4 × Σ_j exp(−(r_i−r_j)²/2σ² − (p_i−p_j)²/2σₚ²)
```

⚠️ **OCCUPATION_SOFTENING = 10.0** reduces the raw 4×sum by a factor of 10. The code docstring (lines 15-25) explains this empirically: in dense ImQMD nuclei, raw 4×sum ≈ 12–40, always clipping to 1.0 and blocking all low-energy collisions. Dividing by 10 gives ~10% acceptance at near-barrier energies. 

This is **not** from a published derivation. It is an empirical tuning parameter justified by the code's comment that "the old 2×sum+1.35-broadening hack was effectively ~2.7×sum; 4/10 = 0.4×sum is comparable but on the lighter-blocking side." The factor of ~10 suppression of the computed occupation represents a very large deviation from the textbook formula.

⚠️ **Recommendation:** (a) The OCCUPATION_SOFTENING should be recognized as the code's largest single free parameter controlling collision acceptance at near-barrier energies. It should be documented prominently in any publication. (b) The value of 10.0 is tuned for σᵣ ≈ 1.3 fm; if the wave-packet width changes (e.g. IQ1 with σᵣ ≈ 0.65 fm), this value would need re-tuning, which should be documented. (c) For publication, consider comparing collision acceptance rates with Δ-like test-particle methods, which handle Pauli blocking through numerical occupancy and do not require this ad-hoc rescaling.

### 2.2 FERMI_SOFTENING = 1.35 for the Fermi Constraint ⚠️

At `collisions.py:34`, `FERMI_SOFTENING = 1.35` is used in `_compute_occupation_field` (line 480):

```
occupations[i] = 4 × Σ(weights) / FERMI_SOFTENING
```

Papa & Bonasera 2001 (PRC 64, 024612) uses the strict Wigner phase-space occupation without softening, applying the constraint at a threshold of 1.0 (one particle per phase-space cell for same-species fermions). FERMI_SOFTENING=1.35 means the effective occupation is ~0.74× the raw value, which slightly relaxes the constraint.

⚠️ The code docstring (lines 27-33) states: "1.35 is the old compensating hack, kept here to match the slightly softened short-range kernel that avoids false-positive violations from the discrete Wigner representation." This is a practical but undocumented numerical artifact. In the CoMD literature (Papa+2001), there is no softening factor — the discrete Wigner representation inherently handles this. The use of 1.35 suggests possible over-counting in the Wigner normalization that merits investigation.

### 2.3 Per-Packet Effective Widths ✅

At `collisions.py:120-134` (`_pair_effective_widths`):

```
σᵢₖ² = (σᵢ² + σₖ²) / 2
σₚ² = ħ²/(4σᵣ²)
```

✅ This correctly handles unequal wave-packet widths for asymmetric systems (projectile and target nuclei with different A, hence different σᵣ from Wang 2002 Eq. 18). The pair width follows from the overlap of two Gaussians: ∫ exp(−(x−rᵢ)²/2σᵢ²) × exp(−(x−rₖ)²/2σₖ²) dx ∝ exp(−(rᵢ−rₖ)²/2(σᵢ²+σₖ²)). The momentum-space width uses the minimum-uncertainty relation σₚ = ħ/(2σᵣ), giving the correct pair-averaged σₚ² = ħ²/(2(σᵢ²+σₖ²)).

### 2.4 Low-Energy NN Cross Sections ✅

At `collisions.py:63-65`:

```
if e_lab < 50.0:
    return 180.0 for np, 60.0 for pp/nn
```

✅ Matches the Cugnon (1987, NPA 470, 558) parameterization and the ImQMD convention documented in Chen+2021 (PRC 103, 064608). Below 50 MeV lab, inelastic channels (Δ production) are closed and the elastic cross sections approach constant values of ~180 mb (np) and ~60 mb (pp/nn near the Coulomb barrier).

---

## 3. Wave Packet Width

**Rating: ✅**

### 3.1 System-Size-Dependent Width

At `initializer.py:164-171` (`compute_sigma_r`):

```python
σᵣ = σ₀ + σ₁ × A^(1/3)
```

✅ Matches Wang 2002 (PRC 65, 064606), Eq. (18) exactly.

**Verification with code values:**
- IQ3A (σ₀=0.94, σ₁=0.020): for ^238U → 0.94 + 0.020 × 238^(1/3) = 0.94 + 0.020 × 6.197 = **1.064 fm**
- IQ3 (σ₀=0.94, σ₁=0.018): for ^238U → 0.94 + 0.018 × 6.197 = **1.052 fm**

### 3.2 Per-Packet Width Propagation ✅

At `nucleus.py:72-81` (`packet_sigmas` property) and `grid_edf.py:32-34, 62-71`, the per-packet σᵣ array is threaded through the entire energy/force pipeline:

- Density deposition: `grid_edf.py:142-162` — each nucleon deposits with its own σᵣ
- Grid construction: `grid_edf.py:674-684` — pads with `n_sigma × max(σᵣ)`
- Analytical forces: `grid_edf.py:278-298` — uses σᵢ for the weight gradient ∇ᵣᵢρᵢ = (x−rᵢ)/σᵢ² × ρᵢ
- Wigner kernels: `collisions.py:120-134` — pair widths via σᵢₖ² = (σᵢ²+σₖ²)/2

✅ This is the correct handling of asymmetric systems where projectile and target have different A and therefore different σᵣ.

---

## 4. NN Cross Sections and In-Medium Correction

**Rating: ⚠️**

### 4.1 Free NN Cross Sections ✅

At `collisions.py:53-74` (`free_nn_cross_section`):

- pp/nn: σ = 13.73 − 15.04/β + 8.76/β² + 68.67·β⁴
- np: σ = −70.67 − 18.18/β + 25.26/β² + 113.85·β

✅ These are the exact Cugnon β-parameterization coefficients (Cugnon+1987, NPA 470, 558). The low-energy cap (180 mb np, 60 mb pp/nn below 50 MeV lab) matches the convention used in ImQMD studies (Zhang+2020 review, Sec. 4.1).

### 4.2 In-Medium Enhancement Factor ⚠️

At `collisions.py:90-111` (`in_medium_factor`):

```
σ_med = (1 + η × sqrt(s) × ρ/ρ₀) × σ_free
```

with η = 0.2 and sqrt(s) in GeV.

At near-barrier energies: sqrt(s) = 2·M_N/1000 ≈ 1.878 GeV, giving:
factor = 1 + 0.2 × 1.878 × ρ/ρ₀ = 1 + 0.376 × (ρ/ρ₀)

So at ρ = ρ₀: **factor = 1.376** (38% enhancement).

✅ The docstring at lines 97-107 correctly notes that the previous code used a reduction form (1 − η·ρ/ρ₀) with the wrong sign, and that this enhancement form matches Chen+2024 (PRC 109, 044609), Eq. (5), which reports extracted in-medium enhancement of ~1.1–2.5× below 150 MeV/u.

⚠️ **Concern about η = 0.2 magnitude:** At near-barrier energies, 1.376× at ρ₀ is near the low end of the Chen+2024 1.1–2.5× range. However, the Chen+2024 paper extracts enhancement from transport-model comparisons with data at 50–150 MeV/u, not at 7 MeV/u (barrier energies). The sqrt(s) scaling means this form predicts the SAME relative enhancement at all energies (because sqrt(s) ≈ 2M_N for any near-threshold collision), which is a built-in feature of this simple parameterization — the energy dependence enters only through the change in sqrt(s) at relativistic energies. At 7 MeV/u, sqrt(s) is nearly identical to 0 MeV/u. This means the same enhancement factor applies from ~5 MeV/u up to ~50 MeV/u, above which sqrt(s) starts increasing.

⚠️ **Recommendation:** (a) Document that η = 0.2 is a tunable parameter and report sensitivity studies. (b) If more precise in-medium corrections are needed, implement a density-and-energy-dependent η(E,ρ) following the Chen+2024 extraction or the FU3FP/FP1 approaches with momentum-dependent interactions.

### 4.3 CM Energy Calculation ✅

At `collisions.py:77-87` (`_sqrt_s_gev`):

```
sqrt(s) = 2 × sqrt(M_N² + M_N × e_cm) / 1000
```

where `e_cm = |p_rel/2|² / M_N` (collisions.py:410, the pair relative kinetic energy per nucleon in the CM). This is correct because in the pair CM frame the invariant energy is:

```
s = (p₁+p₂)² = 2M_N² + 2E₁E₂ − 2p₁·p₂ = 4(M_N² + q²)
```

where q is the CM momentum per nucleon. Since e_cm = q²/M_N, we get sqrt(s) = 2·√(M_N² + M_N·e_cm). ✅

---

## 5. Fermi Constraint (CoMD Phase-Space Occupation)

**Rating: ✅**

### 5.1 Momentum Swap Mechanism ✅

At `collisions.py:484-597` (`fermi_constraint_check`):

The constraint algorithm:
1. Computes Wigner occupation fᵢ for each nucleon from same-species neighbors (line 481, `_compute_occupation_field`)
2. For each nucleon with fᵢ > threshold (1.0), searches nearby same-species candidates
3. Performs momentum **swap** pᵢ ↔ pⱼ (line 590) that maximally reduces fᵢ + fⱼ

✅ The momentum swap conserves Σp (total momentum) and Σp² (total kinetic energy) exactly. This matches Papa & Bonasera 2001 (PRC 64, 024612), the foundational CoMD paper. The docstring at lines 499-505 correctly notes that the previous variant (momentum stretching, then global rescaling) injected energy and was the cause of spurious U+U multifragmentation — the swap is the correct CoMD approach and fixes that energy leak.

### 5.2 Constraint Application Interval ✅

At `propagator.py:326`:

```
fermi_interval = 5  # every 5 fm/c
```

⚠️ Slight divergence from some CoMD implementations: Papa+2001 applies the constraint every time step (1 fm/c), as do many early CoMD codes. However, the docstring at lines 321-325 justifies the 5 fm/c interval: "At near-barrier energies...a 5 fm/c interval gives the phase space time to evolve and create natural vacancies between corrections, preventing the over-blocking (99.99%) that occurs when the constraint is applied every 1 fm/c." This is a reasonable empirical optimization, though it should be noted in publications that the constraint frequency differs from the original.

### 5.3 Group-Aware Constraint ✅

At `collisions.py:508-510` and line 556:

```
mask = (is_proton == is_proton[i]) & (group_ids == group_ids[i])
```

✅ When `group_ids` are provided (0 for projectile, 1 for target), nucleons from different groups do not constrain each other. This allows natural nucleon exchange (transfer) between projectile and target during collisions — a crucial feature for multi-nucleon transfer studies. This matches the CoMD methodology where Pauli blocking respects the projectile/target distinction to avoid artificially suppressing transfer channels.

### 5.4 FERMI_SOFTENING = 1.35 ⚠️

As discussed in §2.2, the softening factor slightly relaxes the Wigner kernel. At near-threshold occupations (fᵢ ≈ 1.0), this means violations are slightly under-detected. The 1.35 factor is empirically tuned and not derived from CoMD theory. However, the impact is modest: an occupation of 1.0 raw becomes ~0.74 softened, while 1.35 raw becomes ~1.0 softened — the latter would still be caught by the threshold=1.0 check.

### 5.5 Pre-Selection Optimization ✅

At `collisions.py:552-554`:

```
spatial_cut2 = (4σ)²
max_candidates = 8
```

✅ These are practical optimizations: the Wigner kernel decays as exp(−Δr²/2σ²), so beyond 4σ the contribution is < exp(−8) ≈ 3×10⁻⁴, negligible. Limiting to 8 candidates is computationally reasonable without sacrificing constraint quality. These should be disclosed as optimizations.

---

## 6. Fragment Recognition MST Parameters

**Rating: ⚠️**

### 6.1 Standard MST ✅

At `fragments.py:85-107` (`minimum_spanning_tree`):

Defaults: r_cut = 3.5 fm, p_cut = 300 MeV/c

⚠️ **p_cut discrepancy with classic QMD value:** The classic QMD MST uses p_cut = 250 MeV/c (Aichelin 1991, Phys. Rep. 202, 233; Zhang 2020 review, Eq. 58). The code uses 300 MeV/c, which is noted as the "imQMD MNT criterion" (research report Sec. 4.5). This ~20% larger momentum cutoff will recognize more fragments in momentum space, slightly increasing fragment multiplicity. This should be noted in publications.

✅ **r_cut = 3.5 fm** is within the published range of 2.5–3.5 fm for imQMD MNT studies.

### 6.2 Isospin-Dependent MST (iso-MST) ❌

At `fragments.py:110-135` (`isospin_mst`):

Defaults: r_cut_pp = 3.0 fm, r_cut_nn = 6.0 fm, r_cut_np = 6.0 fm

| Cutoff | Value | Typical in literature |
|--------|-------|----------------------|
| pp | 3.0 fm | 2.8–3.4 fm |
| nn | **6.0 fm** | 3.4–4.0 fm |
| np | 6.0 fm | 3.4–4.0 fm |

❌ **The nn and np cutoffs at 6.0 fm are anomalous.** A factor of 2 between pp and nn/np cutoffs is very large compared to published iso-MST implementations, where the asymmetry is typically ~20–40% (e.g. r_cut_nn/r_cut_pp ≈ 1.2–1.4 based on the neutron skin thickness). 

The typical iso-MST from imQMD literature (Li+2013, PRC 87, 064615; Zhang+2020 review, Sec. 4.2) has more modest asymmetry:
- pp ≈ 2.8 fm, nn ≈ 3.8 fm, np ≈ 3.4 fm

A 6.0 fm cutoff means any two neutrons within 6 fm are considered part of the same fragment, regardless of their momentum difference. This will over-merge neutron-rich fragments and artificially suppress free neutron emission.

❌ **Recommendation:** The iso-MST cutoff asymmetry should be reduced. Values closer to pp=3.0, nn=3.8, np=3.4 fm would be more consistent with published iso-MST implementations. If the 6.0 fm nn/np cutoffs are intentionally tuned for a specific system (e.g. ^238U+^238U MNT to improve neutron-rich fragment yields), this should be documented and sensitivity studies should be reported.

### 6.3 Adaptive MST ⚠️

At `fragments.py:138-166` (`adaptive_mst`):

```
r_cut = min(max_r_cut, base_r_cut × ∛(A_fragment/16))
```

with defaults: base_r_cut = 3.2 fm, max_r_cut = 4.0 fm.

⚠️ This A^(1/3) scaling of the MST cutoff is a reasonable heuristic — larger fragments are more diffuse, so the spatial cutoff should scale with the fragment radius (which scales as A^(1/3)). However, this specific implementation (normalizing to A=16 and capping at 4.0 fm) does not appear to be from any published paper. It is an internal heuristic.

⚠️ **Recommendation:** Document the adaptive MST algorithm in any publication, including: (a) the reference A=16 normalization, (b) the 4.0 fm cap, and (c) rationale for using the larger of base_r_cut and the scaled cutoff (line 160: `np.maximum(cutoffs, base_r_cut)`).

---

## 7. Statistical Decay

**Rating: ✅**

### 7.1 Weisskopf Evaporation Model ✅

At `decay.py:210-219` (`_weisskopf_width`):

The Weisskopf width calculation uses the standard formalism:

```
Γ = ∫₀^{E*} (σ_inv × ρ_daughter(E*−ε) × ε) dε
```

with:
- Inverse cross section: σ_inv(n) = πR² (geometrical), σ_inv(p/α) = πR² × E/(E+V_coul)
- Level density: ∝ exp(2√(aU)) / U^(5/4)
- Level density parameter: a = A/8 (standard Fermi gas)

✅ All components match the standard Weisskopf statistical-model formalism. The parameter a = A/8 is consistent with the commonly used a = A/8 to A/10 range.

### 7.2 Fission Competition ✅

At `decay.py:230-246` (`_bohr_wheeler_width`):

The Bohr-Wheeler fission width:

```
Γ_f = (1/(2πρ_parent)) × ∫₀^{E*−B_f} ρ_saddle(E*−B_f−ε) dε
```

✅ Standard Bohr-Wheeler formula for statistical fission. Uses different level-density parameters for ground-state (a_f/a_n = 1.04) as conventional. The fission barrier `fission_barrier` (lines 162-176) combines a liquid-drop term 98·(1−x)² with Gaussian shell corrections around N=126 and Z=82, matching standard macroscopic-microscopic barrier parameterizations.

✅ The `_fission_to_neutron_ratio` function (lines 249-264) provides a floor on the fission width for very heavy nuclei (Z > 90), ensuring fission competes appropriately even when the barrier is small.

### 7.3 HIVAP Wrapper ✅

At `hivap_wrapper.py`:

The wrapper writes HICOL input format `input.dat`, runs the Fortran binary, and parses `SIGXPN.DAT`. The dummy-target trick (A1=A−1, Z1=Z−1, A2=1, Z2=1 to form compound nucleus (Z,A)) is standard practice for HIVAP. 

✅ HIVAP is an established statistical-model code widely used in imQMD coupling studies. The 500 fm/c switch time (`reaction.py:330`) matches the optimum found in the research report.

### 7.4 Excitation Energy Calculation ✅

At `fragments.py:254-290` (`compute_fragment_excitation`):

```
E* = E_int − E_ground
```

where E_int is the GridEDF total energy of fragment nucleons after:
1. Subtracting CM momentum (line 268)
2. Subtracting collective rotational energy via `_subtract_collective_rotation` (lines 293-308)

✅ The collective rotation subtraction:
- Computes total angular momentum L = Σ r × p (line 298)
- Builds the moment-of-inertia tensor I = Σ M_N (r²·𝟙 − r⊗r) (lines 299-302)
- Solves ω = I⁻¹·L (line 304)
- Subtracts rotational momentum p_rot = M_N(ω×r) (line 308)

This is the correct Newtonian rigid-body rotation subtraction for a fragment. It correctly separates internal excitation from collective rotation of the dinuclear system. This matches the sophisticated E* methodology in the research report.

### 7.5 Internal vs. HIVAP Adequacy ⚠️

The internal Weisskopf model (decay.py) handles n, p, α, γ, and fission channels. It uses a = A/8 and simplified inverse cross sections. For heavy-fragment de-excitation, the HIVAP path is preferred (`hivap_wrapper.py`). The internal model has not been benchmarked against GEMINI, which is the standard imQMD coupling code.

⚠️ **Recommendation:** If the internal Weisskopf model is used for production (not just as a fallback), it should be benchmarked against GEMINI for:
- Neutron evaporation multiplicities for ^238U residues
- Fission survival probabilities for Z=90–100 fragments
- Alpha/light-charged-particle branching ratios

### 7.6 Switch Time ✅

At `reaction.py:330`:

```
time_fm_c = 500.0  # fm/c
```

✅ The 500 fm/c dynamics→statistical-decay switch time matches the optimum reported in imQMD+GEMINI studies (typically 400–600 fm/c, research report Sec. 4.6). At this time, primary fragments are formed but still carry significant excitation energy (~20–60 MeV depending on impact parameter), and the subsequent statistical decay handles de-excitation. The additional 250 fm/c collisionless cooling (`reaction.py:373-377`) before fragment recognition allows the system to settle into clearly separated fragments.

---

## Code Quality Notes

### Undocumented Parameters and Deviations Summary

| Item | Parameter | Deviation from published | Severity |
|------|-----------|-------------------------|----------|
| 1.2 | IQ3A kappa_s=0.4 | 5× above IQ3 (0.08). No cited paper. | **High** |
| 1.2 | IQ3A gsur=16.5 | ~8% below IQ3. No cited paper. | Medium |
| 1.2 | IQ3A c_sym=34 | ~6% above IQ3. No cited paper. | Medium |
| 2.1 | OCCUPATION_SOFTENING=10.0 | No theoretical basis, empirically tuned. | **High** |
| 2.2 | FERMI_SOFTENING=1.35 | CoMD original uses 1.0, no softening. | Low |
| 5.2 | Fermi interval = 5 fm/c | Original CoMD uses 1 fm/c. | Low |
| 6.1 | p_cut = 300 MeV/c | Classic QMD uses 250 MeV/c. | Medium |
| 6.2 | iso-MST nn/np = 6.0 fm | Factor of 2 asymmetry vs. ~1.4 in literature. | **High** |
| 6.3 | Adaptive A^(1/3) scaling | No published reference. | Low |

### Notable Implementation Strengths

1. **Energy-force self-consistency** (`grid_edf.py:439-521`): The analytical forces use the exact transpose Dᵀ of the discrete gradient stencil, not a second np.gradient call. This ensures that the forces are the exact gradient of the energy functional (to finite-difference accuracy at interior points), which is superior to most published transport-code implementations that approximate with -∇².

2. **Collective rotation correction** (`fragments.py:293-308`): Most fragment excitation codes only subtract CM kinetic energy. The additional rotational correction is physically important for large-impact-parameter events where fragments acquire significant angular momentum.

3. **Per-packet width threading** (`grid_edf.py:86-93`, `nucleus.py:72-81`): The entire code correctly propagates per-nucleon wave-packet widths, handling the system-size-dependent σᵣ of Wang 2002 Eq. (18) for asymmetric systems. The Coulomb forces, EDF forces, and Wigner kernels all use the appropriate pair-averaged widths.

4. **Momentum-swap Fermi constraint** (`collisions.py:590`): The code explicitly documents that the previous momentum-stretch variant leaked energy and uses the correct CoMD momentum swap, which is rigorously energy-conserving.

5. **Two-point nuclear-scale calibration** (`initializer.py:213-240`): The `_fit_grid_nuclear_scale_twopoint` function correctly handles the fact that the local scale s(ρ) has a ramp contribution at ρ > ρ₀, making the energy exactly linear in the calibration parameter. This is a subtle but correct improvement over the simpler one-point fit.

### Other Code Quality Observations

- **Deprecation markers**: Several functions are marked DEPRECATED with clear migration paths (`nucleus.py:128-133` centroid path, `propagator.py:115-117`, `propagator.py:131-133`). Good practice.
- **Caching**: Per-nucleus initialization cache, fragment scale cache, and cold-ground-state cache are properly keyed by all relevant parameters including EDF parameter set name.
- **Collision neighbor list** (`collisions.py:305-377`): Intelligent spatial-hash cache with rebuild triggers (time interval + position displacement).
- **Seed handling**: All randomness is properly seeded and reproducible. The `_collision_rng` is stored on the nucleus object to prevent seed interference.
- **Error estimation**: Cross sections include Poisson error bars (`reaction.py:543-546`, `reaction.py:596-597`), which should be standard for such calculations.
- **Fallback paths**: The HIVAP wrapper falls back to the internal Weisskopf chain on failure (`hivap_wrapper.py:192-195`), and the scipy imports have non-scipy fallbacks in several places.
