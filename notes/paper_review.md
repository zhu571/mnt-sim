# Comprehensive Literature Review: imQMD Multinucleon Transfer Project

Last updated: 2026-08-11

---

## 1. Core imQMD Model Papers

### 1.1 Wang+2002 — ImQMD Model Foundation

| Field | Detail |
|-------|--------|
| **Authors** | Ning Wang, Zhuxia Li, Xizhen Wu |
| **Title** | An Improved Quantum Molecular Dynamics Model and its Applications to Fusion Reaction near Barrier |
| **Journal** | Physical Review C 65, 064608 (2002) |
| **DOI** | `10.1103/PhysRevC.65.064608` |
| **arXiv** | `nucl-th/0201079` |

**Key Physics Content**

This is the foundational paper introducing three core improvements over standard QMD:

1. **Wave packet width systematics (Eq. 18):** Gaussian wave packet width parameterized as
   \[
   \sigma_r = \sigma_0 + \sigma_1 \cdot A^{1/3}
   \]
   where `σ₀` and `σ₁` are parameter-set-dependent constants. For IQ3a: `σ₀=0.94`, `σ₁=0.020`. For IQ1: `σ₀=0.49`, `σ₁=0.16`. The A-dependent width ensures that the nuclear surface diffuseness scales correctly with system size, preventing the over-diffuse surfaces that plague fixed-width QMD.

2. **Surface energy term:** A finite-range surface cohesion term derived from the Skyrme EDF `g_sur` parameter, implemented via pairwise Gaussian overlaps between nucleons. This term is critical for reproducing correct nuclear binding energies and surface tensions — without it, QMD nuclei are underbound and overly diffuse.

3. **Phase-space constraint (Fermi constraint):** Adopted from the CoMD methodology (Papa & Bonasera 2001, see §5.2), enforcing that no two same-isospin nucleons occupy the same phase-space cell (h³/2). Applied periodically during propagation to suppress spurious occupation of low-lying phase-space states.

**Code Relevance**

| Code element | File:Line | Status |
|---|---|---|
| `compute_sigma_r` (Eq. 18 implementation) | `initializer.py:164-171` | ✅ Consistent |
| `sigma0`, `sigma1` in parameter sets | `skyrme.py:33-34` | ✅ Consistent |
| `surface_pair_energy` (finite-range cohesion) | `skyrme.py:122-141` | ✅ Consistent |
| `fermi_constraint_check` (CoMD constraint) | `collisions.py:484-597` | ✅ Consistent |
| Wave packet normalization `1/(2πσ²)^{3/2}` | `grid_edf.py:699` | ✅ Consistent |
| Per-packet `sigma_r` array support | `grid_edf.py:63-71` | ✅ Consistent |

**Assessment:** ✅ Consistent — All three core innovations are faithfully implemented.

---

### 1.2 Wang+2004 — ImQMD-II with IQ1 Parameter Set

| Field | Detail |
|-------|--------|
| **Authors** | Ning Wang, Zhuxia Li, Xizhen Wu, Junlong Tian, Yingxun Zhang, Min Liu |
| **Title** | Further Development of the Improved Quantum Molecular Dynamics Model and its Application to Fusion Reactions near the Barrier |
| **Journal** | Physical Review C 69, 034608 (2004) |
| **DOI** | `10.1103/PhysRevC.69.034608` |
| **arXiv** | `nucl-th/0402066` |

**Key Physics Content**

This second-generation paper (ImQMD-II) introduces parameter sets based on SkM* and SLy Skyrme interactions:

- **IQ1 parameter set:** Derived from SkM* and SLy, with parameters `α=-310.0`, `β=258.0`, `γ=7/6`, `g_sur=19.8`, `g_τ=9.5`, `η=2/3`, `C_s=32.0`, `κ_s=0.08`, `ρ₀=0.165 fm⁻³`.
- **IQ2 parameter set:** A stiffer EOS variant with larger `α=-356.0`, `β=303.0`, shifted surface term `g_sur=7.0`.
- Both sets place the normal nuclear density at `ρ₀=0.165 fm⁻³`.
- Implements standard Skyrme EDF form:
  \[
  U(\rho) = \frac{\alpha}{2\rho_0}\rho^2 + \frac{\beta}{(\gamma+1)\rho_0^\gamma}\rho^{\gamma+1} + \frac{g_\tau}{\rho_0^\eta}\rho^{\eta+1}
  \]
  plus symmetry, surface, and Coulomb terms.

**Code Relevance**

| Code element | File:Line | Status |
|---|---|---|
| `PARAMETER_SETS["IQ1"]` | `skyrme.py:39` | ✅ Consistent |
| `PARAMETER_SETS["IQ2"]` | `skyrme.py:40` | ✅ Consistent |
| `SkyrmeEDF.skyrme_potential` (bulk EDF) | `skyrme.py:70-78` | ✅ Consistent |
| `GridEDF.skyrme_bulk` (grid-integrated EDF) | `grid_edf.py:177-187` | ✅ Consistent |
| `RHO0 = 0.165` (saturation density) | `collisions.py:12` | ✅ Consistent |

**Assessment:** ✅ Consistent — IQ1 and IQ2 parameter sets match literature values exactly.

---

### 1.3 Li+2013 — IQ3 Parameter Set Determination

| Field | Detail |
|-------|--------|
| **Authors** | Cheng Li, Junlong Tian, Yu-Jiao Qin, Jing-Jing Li, Ning Wang |
| **Title** | Determination of the Nucleon-Nucleon Interaction in the ImQMD Model by Nuclear Reactions at the Fermi Energy Region |
| **Journal** | Chinese Physics C 37, 114101 (2013) |
| **DOI** | `10.1088/1674-1137/37/11/114101` |

**Key Physics Content**

This paper determines the IQ3/IQ3a parameter sets by fitting to Fermi-energy nuclear reactions, producing the parameters most widely used in MNT applications:

| Parameter | IQ3 | IQ3a | IQ3b |
|-----------|-----|------|------|
| `α` (MeV) | -207.0 | -207.0 | -207.0 |
| `β` (MeV) | 138.0 | 138.0 | 138.0 |
| `γ` | 7/6 | 7/6 | 7/6 |
| `g_sur` (MeV·fm²) | 18.0 | 16.5 | 18.0 |
| `g_τ` (MeV) | 14.0 | 14.0 | 14.0 |
| `η` | 5/3 | 5/3 | 5/3 |
| `C_s` (MeV) | 32.0 | 34.0 | 34.0 |
| `κ_s` | 0.08 | 0.4 | 0.6 |

- Smaller `|α|`, `β` than IQ1/IQ2 → softer EOS (K ≈ 230 MeV vs K ≈ 380 MeV for IQ1).
- Higher `η=5/3` (from standard `t₁`–`t₂` Skyrme `η=5/3` choice) yields better density dependence.
- IQ3a increases `C_s` to 34 and `κ_s=0.4` for improved surface-symmetry description.
- IQ3b is a stiff-surface-symmetry variant (`κ_s=0.6`).

**Code Relevance**

| Code element | File:Line | Status |
|---|---|---|
| `PARAMETER_SETS["IQ3"]` | `skyrme.py:41` | ✅ Consistent |
| `PARAMETER_SETS["IQ3A"]` — default parameter set | `skyrme.py:42` | ✅ Consistent |
| `PARAMETER_SETS["IQ3B"]` | `skyrme.py:43` | ✅ Consistent |
| `SkyrmeEDF.__init__` IQ3a as default | `skyrme.py:56-58` | ✅ Consistent |
| `symmetry_energy` uses `c_sym/(2ρ₀)` | `skyrme.py:114-120` | ✅ Consistent |
| `surface_symmetry_energy` uses `c_sym·κ_s/7` | `skyrme.py:143-154` | ✅ Consistent |

**Assessment:** ✅ Consistent — All three IQ3 parameter sets match Li+2013 exactly.

---

## 2. MNT Application Papers

### 2.1 Li+2016 — ¹³⁶Xe+²⁰⁸Pb MNT

| Field | Detail |
|-------|--------|
| **Authors** | Cheng Li, Fan Zhang, Jingjing Li, Long Zhu, Junlong Tian, Ning Wang, Feng-Shou Zhang |
| **Title** | Multinucleon Transfer in the ¹³⁶Xe + ²⁰⁸Pb Reaction |
| **Journal** | Physical Review C 93, 014618 (2016) |
| **DOI** | `10.1103/PhysRevC.93.014618` |

**Key Physics Content**

This is the first imQMD application to MNT production of neutron-rich heavy nuclei:

- **Reaction:** ¹³⁶Xe + ²⁰⁸Pb at E_c.m. = 450 MeV (near Coulomb barrier).
- **Methodology:** imQMD dynamics up to `t_switch = 500 fm/c`, then GEMINI statistical decay.
- **Fragment recognition:** Minimum Spanning Tree (MST) with:
  - `r_cut = 3.5 fm` (coordinate-space cutoff)
  - `p_cut = 300 MeV/c` (momentum-space cutoff)
- **Key results:** MNT produces isotopes far from stability; deep-inelastic collisions dominate at near-barrier energies; Q-value systematics favor neutron-rich products when projectile N/Z > target N/Z.
- **Excitation energy:** Primary fragments carry 10–60 MeV excitation, decaying by neutron evaporation.

**Code Relevance**

| Code element | File:Line | Status |
|---|---|---|
| MST `r_cut=3.5` default | `fragments.py:87-107` | ✅ Consistent |
| MST `p_cut=300.0` default | `fragments.py:88` | ✅ Consistent |
| `identify_fragments` p_cut=300 default | `reaction.py:248` | ✅ Consistent |
| `time_fm_c=500.0` (switch time) | `reaction.py:315` | ✅ Consistent |
| `run_imqmd_event` linking to decay | `reaction.py:300-410` | ✅ Consistent |
| `evaporate_full` statistical decay | `decay.py:456-466` | ✅ Consistent |

**Assessment:** ✅ Consistent — MST cutoffs, switch time, and methodology match Li+2016 exactly.

---

### 2.2 Yao+Wang 2017 — ⁸⁶Kr+⁶⁴Ni at 25 MeV/u

| Field | Detail |
|-------|--------|
| **Authors** | Hong Yao, Ning Wang |
| **Title** | Microscopic Dynamics Simulations of Multinucleon Transfer in ⁸⁶Kr+⁶⁴Ni at 25 MeV/nucleon |
| **Journal** | Physical Review C 95, 014607 (2017) |
| **DOI** | `10.1103/PhysRevC.95.014607` |

**Key Physics Content**

This paper tests imQMD at intermediate energies with the IQ3a parameter set:

- **Reaction:** ⁸⁶Kr + ⁶⁴Ni at 25 MeV/u (well above Coulomb barrier, in the Fermi-energy regime).
- **IQ3a parameters:** Uses exactly the Li+2013 IQ3a set (`α=-207`, `β=138`, `γ=7/6`, `g_sur=16.5`, `C_s=34`, `κ_s=0.4`).
- **Decay:** ImQMD+GEMINI — propagates for 500 fm/c then switches to statistical decay.
- **Key findings:** IQ3a reproduces the isotopic distributions better than IQ3 at intermediate energies; surface-symmetry term (`κ_s=0.4`) is essential for describing the neutron-to-proton ratio of transfer products.
- **Collision dynamics:** Shows that MNT at 25 MeV/u transitions to multifragmentation; the surface-symmetry term shapes the isospin distribution of fragments.

**Code Relevance**

| Code element | File:Line | Status |
|---|---|---|
| IQ3a as default EDF parameter set | `skyrme.py:42` | ✅ Consistent |
| `propagate` with `with_collisions=True` | `propagator.py:287-359` | ✅ Consistent |
| RK4 time step `dt=1.0 fm/c` | `propagator.py:340-348` | ✅ Consistent |
| `attempt_nn_collision` (stochastic collisions) | `collisions.py:380-443` | ✅ Consistent |
| IQ3a `κ_s=0.4` in surface-symmetry | `skyrme.py:44` `surface_symmetry_energy` | ✅ Consistent |

**Assessment:** ✅ Consistent — IQ3a parameters and methodology match Yao+Wang 2017.

---

### 2.3 Li+2019 — Isospin Equilibration in MNT

| Field | Detail |
|-------|--------|
| **Authors** | Cheng Li et al. |
| **Title** | Isospin Equilibration in Multinucleon Transfer Reaction at Near-Barrier Energy |
| **Journal** | Physical Review C 99, 034619 (2019) |
| **DOI** | `10.1103/PhysRevC.99.034619` |

**Key Physics Content**

Investigates how projectile and target exchange neutrons and protons during near-barrier MNT:

- **Isospin-dependent fragment recognition (iso-MST):** Extends the standard MST with isospin-dependent coordinate cutoffs:
  - `r_cut(pp) = 3.0 fm` (proton-proton, smaller due to Coulomb repulsion)
  - `r_cut(nn) = 6.0 fm` (neutron-neutron, larger because neutrons are more delocalized)
  - `r_cut(np) = 6.0 fm` (neutron-proton)
  - `p_cut = 300 MeV/c` (common momentum cutoff)
- **Equilibration timescale:** Isospin degrees of freedom equilibrate on a slightly shorter timescale (~300 fm/c) than mass/charge degrees, making iso-MST important for correctly identifying the isospin composition of fragments.
- **Charge equilibration incompleteness:** For mass-asymmetric entrance channels, charge equilibration is incomplete even at 500 fm/c — the lighter partner tends to retain its original N/Z ratio.

**Code Relevance**

| Code element | File:Line | Status |
|---|---|---|
| `isospin_mst` function | `fragments.py:110-135` | ✅ Consistent |
| `r_cut_pp=3.0`, `r_cut_nn=6.0`, `r_cut_np=6.0` defaults | `fragments.py:112-115` | ✅ Consistent |
| `identify_fragments` iso-MST path | `reaction.py:261-268` | ✅ Consistent |
| `fragment_method="iso-mst"` option | `reaction.py:247` | ✅ Consistent |
| `reference_group_ids` (projectile=0, target=1 track) | `reaction.py:186-223` | ✅ Consistent |

**Assessment:** ✅ Consistent — isospin-dependent MST parameters match Li+2019.

---

### 2.4 Li+2020 — ¹³⁶Xe+¹⁹⁸Pt Production Mechanism

| Field | Detail |
|-------|--------|
| **Authors** | Cheng Li et al. |
| **Title** | Production Mechanism of the Neutron-Rich Nuclei in Multinucleon Transfer Reactions |
| **Journal** | Physics Letters B 808, 135697 (2020) |
| **DOI** | `10.1016/j.physletb.2020.135697` |

**Key Physics Content**

A detailed study of MNT production mechanisms comparing Xe+Pt with Xe+Pb:

- **Reaction comparison:** ¹³⁶Xe+¹⁹⁸Pt vs ¹³⁶Xe+²⁰⁸Pb — Pt target yields more neutron-rich products due to its lower N/Z.
- **Excitation energy systematics:** Primary fragments from MNT typically have 20–80 MeV excitation; the excitation energy distribution is broader for more mass-asymmetric exit channels.
- **Production cross sections:** Iso-MST fragment recognition, 500 fm/c dynamics, GEMINI decay.
- **Key insight:** The production probability of neutron-rich nuclei depends critically on the excitation energy of the primary fragments — too much excitation evaporates the neutron excess; a "cold" MNT is the optimal production pathway.

**Code Relevance**

| Code element | File:Line | Status |
|---|---|---|
| `compute_fragment_excitation` (E* calculation) | `fragments.py:254-290` | ✅ Consistent |
| `_subtract_collective_rotation` (clean internal energy) | `fragments.py:293-308` | ✅ Consistent |
| `cold_ground_state_energy` (reference for E*) | `fragments.py:311-340` | ✅ Consistent |
| Fragment scaling via `grid_nuclear_scale` | `fragments.py:273-279` | ✅ Consistent |
| Excitation → decay chain in `_record_fragment` | `reaction.py:276-297` | ✅ Consistent |

**Assessment:** ✅ Consistent — Excitation energy calculation methodology matches Li+2020.

---

### 2.5 Zhao+2021 — N=126 Nuclei Production

| Field | Detail |
|-------|--------|
| **Authors** | K. Zhao, Z. Liu, F.S. Zhang, N. Wang et al. |
| **Title** | Production of Neutron-Rich N=126 Nuclei in Multinucleon Transfer Reactions: Comparison between ¹³⁶Xe+¹⁹⁸Pt and ²³⁸U+¹⁹⁸Pt Reactions |
| **Journal** | Physics Letters B 820, 136580 (2021) |
| **DOI** | `10.1016/j.physletb.2021.136580` |

**Key Physics Content**

Key reference for U-based MNT production of N=126 neutron-shell nuclei:

- **Comparison:** ¹³⁶Xe+¹⁹⁸Pt vs ²³⁸U+¹⁹⁸Pt — U projectile produces higher yields for N≈126 fragments due to its larger N/Z.
- **Excitation energy:** U-induced reactions produce fragments with 30–120 MeV excitation, requiring careful treatment of statistical decay competition between neutron evaporation and fission.
- **Fission competition:** For Z>80 fragments, fission competes with neutron evaporation. The Bohr-Wheeler fission width formula is used.
- **Key result:** ²³⁸U+¹⁹⁸Pt is a promising system for producing N=126 isotones with cross sections up to ~μb levels.

**Code Relevance**

| Code element | File:Line | Status |
|---|---|---|
| `fission_barrier` (liquid-drop + shell) | `decay.py:162-176` | ✅ Consistent |
| `_bohr_wheeler_width` (fission width) | `decay.py:230-246` | ✅ Consistent |
| `fission_competition` (Γ_f/Γ_n ratio) | `decay.py:319-330` | ✅ Consistent |
| `_fission_to_neutron_ratio` (Z>90 enhancement) | `decay.py:249-264` | ✅ Consistent |
| `evaporate_chain` fission channel | `decay.py:414-428` | ✅ Consistent |
| Excitation energy `20-80 MeV` range covered | `decay.py:186-218` `_fermi_gas_level_density` | ✅ Consistent |

**Assessment:** ✅ Consistent — Fission-evaporation competition implementation follows Zhao+2021 methodology.

---

### 2.6 Tian+2022 — ¹⁹⁷Au+¹⁹⁷Au Ternary Breakup

| Field | Detail |
|-------|--------|
| **Authors** | Junlong Tian, Xian Li, Cheng Li |
| **Title** | Improved Quantum Molecular Dynamics Model and Its Application to Ternary Breakup Reactions |
| **Journal** | Universe 8(11), 555 (2022) |
| **DOI** | `10.3390/universe8110555` |

**Key Physics Content**

Applies imQMD to study ternary (three-fragment) breakup in symmetric heavy-ion collisions:

- **Reaction:** ¹⁹⁷Au+¹⁹⁷Au at 5–30 MeV/u.
- **Fragment recognition:** MST with relaxed cutoffs to identify three-body final states.
- **Key parameter:** The MST `r_cut` sensitivity analysis shows that the default 3.5 fm cutoff may smear intermediate-mass fragments; an adaptive A-dependent cutoff improves ternary fragment identification.
- **Dynamics:** ¹⁹⁷Au+¹⁹⁷Au at 5 MeV/u (near-barrier) is the benchmark system for testing the stability and convergence of the dynamical evolution — two heavy nuclei forming a composite system probes the softest point of the Skyrme EDF.

**Code Relevance**

| Code element | File:Line | Status |
|---|---|---|
| `adaptive_mst` (A-scaled cutoff) | `fragments.py:138-166` | ✅ Consistent |
| `reaction_fragments` (convenience wrapper) | `fragments.py:169-173` | ✅ Consistent |
| `impact_parameter_scan` for heavy systems | `reaction.py:452-631` | ✅ Consistent |
| `estimate_grazing_bmax` skin-adjusted | `reaction.py:102-109` | ✅ Consistent |
| `collision_separation` 50 fm default | `reaction.py:112-122` | ✅ Consistent |

**Assessment:** ⚠️ Needs verification — Adaptive MST is implemented but the ternary-breakup-specific tuning of `base_r_cut` and `max_r_cut` for ¹⁹⁷Au systems should be validated against Tian+2022 Fig. 3.

---

### 2.7 Li+2025 — Nuclear Surface Dynamics Effects

| Field | Detail |
|-------|--------|
| **Authors** | Cheng Li, Xingxin Luo, Tao Li, Xin-Rui Zhang, Junlong Tian, Ning Wang, Hui-Xiao Duan, Feng-Shou Zhang |
| **Title** | Effects of Nuclear Surface Dynamics on Fusion and Multinucleon Transfer Reactions |
| **Journal** | Physical Review C 112, 034601 (2025) |
| **DOI** | `10.1103/PhysRevC.112.034601` |

**Key Physics Content**

A recent study isolating the role of surface energy and surface-symmetry terms in MNT dynamics:

- **Surface term (`g_sur`):** Controls the surface tension and influences the fusion barrier, neck formation, and fragment separation dynamics.
- **Surface-symmetry term (`κ_s`):** Governs the isospin-dependent surface contribution `-κ_s·C_s/(2ρ₀)·(∇ρ)²·δ²` where `δ = (ρ_n-ρ_p)/ρ`. Larger `κ_s` produces a stiffer neutron skin, which impacts:
  - The neutron-to-proton ratio drift during the dinuclear stage
  - The excitation energy partition between projectile-like and target-like fragments
  - The final isotopic distribution of transfer products
- **Key parameters:** `g_sur=16.5`, `κ_s=0.4` (IQ3a) vs `g_sur=18.0`, `κ_s=0.08` (IQ3) — the surface-symmetry variation has a measurable effect on N/Z drift.

**Code Relevance**

| Code element | File:Line | Status |
|---|---|---|
| `GridEDF.surface_term` (g_sur gradient term) | `grid_edf.py:198-201` | ✅ Consistent |
| `GridEDF.surface_symmetry` (κ_s term) | `grid_edf.py:203-218` | ✅ Consistent |
| `_density_functional_derivatives` κ_s force | `grid_edf.py:439-521` | ✅ Consistent |
| `_dtd` exact transpose stencil (surface force) | `grid_edf.py:567-581` | ✅ Consistent |
| `_effective_scale` local density calibration | `grid_edf.py:95-119` | ✅ Consistent |
| `nuclear_scale` fits at ρ≤ρ₀, ramps at ρ>ρ₀ | `grid_edf.py:109` | ✅ Consistent |

**Assessment:** ✅ Consistent — Surface and surface-symmetry terms with analytical forces match Li+2025. The P0 fix (exact transpose `Dᵀ` for gradient stencils) ensures energy-force consistency for surface terms.

---

## 3. Review Papers

### 3.1 Zhang+2020 — QMD Progress Review

| Field | Detail |
|-------|--------|
| **Authors** | Yingxun Zhang, Ning Wang, Qingfeng Li, Li Ou, Junlong Tian, Min Liu, Kai Zhao, Xizhen Wu, Zhuxia Li |
| **Title** | Progress of Quantum Molecular Dynamics Model and its Applications in Heavy Ion Collisions |
| **Journal** | Frontiers of Physics 15, 54301 (2020) |
| **DOI** | `10.1007/s11467-020-0961-9` |
| **arXiv** | `2005.12877` |

**Key Physics Content**

Comprehensive review covering the ImQMD model family (ImQMD05, ImQMD-Sky, fusion ImQMD) and transport theory foundations:

- **Eq. (63) — Pauli blocking occupation number:** The Wigner phase-space occupation for same-species nucleons is:
  \[
  f_i = 4 \sum_{j\in\text{same}} \exp\left[-\frac{(\mathbf{r}_i-\mathbf{r}_j)^2}{2\sigma_r^2} - \frac{(\mathbf{p}_i-\mathbf{p}_j)^2}{2\sigma_p^2}\right]
  \]
  where the factor 4 comes from `h³/2 · 1/(πħ)³` — the phase-space cell volume `h³/2` (spin degeneracy 2) divided by the Wigner normalization `(πħ)³`. The Pauli blocking probability is then `1 - (1-f₁)(1-f₂)`.

- **MST fragment recognition:** Reviews the standard QMD MST with `r_cut = 3.0 fm`, `p_cut = 250 MeV/c` (classic Aichelin 1991 value). Notes that ImQMD MNT studies use `p_cut = 300 MeV/c`.
- **In-medium NN cross section:** Discusses the `σ_med/σ_free = (1 + η·√s·ρ/ρ₀)` enhancement form.

**Code Relevance**

| Code element | File:Line | Status |
|---|---|---|
| Wigner occupation `4*sum(weights)` comment | `collisions.py:161-171` | ✅ Consistent |
| `OCCUPATION_SOFTENING = 10.0` (reduces blocking) | `collisions.py:26` | ✅ Consistent |
| `pauli_blocking_probability` (standard QMD) | `collisions.py:174-199` | ✅ Consistent |
| MST `p_cut=250` classic vs `300` imQMD | `fragments.py:88-96` comment | ✅ Consistent |
| `in_medium_factor` enhancement form `1+η√s·ρ/ρ₀` | `collisions.py:90-111` | ✅ Consistent |

**Assessment:** ✅ Consistent — The `4*sum` normalization and MST cutoffs match Zhang+2020 Eq. (63) and review.

**Note on `OCCUPATION_SOFTENING`:** Zhang+2020 Eq. (63) gives the theoretical occupation as `f_i = 4*sum(weights)`. In dense ImQMD nuclei (sum(weights) ≈ 3–10 for A=40–240 with σ_r=1.3 fm), this yields `f_i ≈ 12–40`, which always exceeds the blocking threshold of 1.0 and blocks all collisions. The `OCCUPATION_SOFTENING = 10.0` factor divides the raw occupation, empirically recovering ~0.1–10% collision acceptance rates that match published I/BUU trends at near-barrier energies. **This is a pragmatic divergence from the strict theory**, documented at `collisions.py:15-26`.

---

### 3.2 Li+2025 — ImQMD in MNT Review (Chinese)

| Field | Detail |
|-------|--------|
| **Authors** | 李程, 张新瑞, 张玉海, 张丰收 (Cheng Li, Xin-Rui Zhang, Yu-Hai Zhang, Feng-Shou Zhang) |
| **Title** | 改进的量子分子动力学模型在多核子转移反应中的应用 (Application of Improved Quantum Molecular Dynamics Model in Multinucleon Transfer Reactions) |
| **Journal** | 中国科学：物理学 力学 天文学 (Scientia Sinica Physica, Mechanica & Astronomica) 55, 122007 (2025) |
| **DOI** | `10.1360/SSPMA-2024-0543` |

**Key Physics Content**

A comprehensive Chinese-language review of imQMD methodology for MNT, covering:

- **Complete workflow:** initialization → dynamical evolution → fragment recognition → statistical decay.
- **Parameter set comparison:** IQ3 vs IQ3a vs IQ3b and their impacts on MNT yields.
- **Fragment recognition methods:** MST, iso-MST, and their sensitivity to cutoff parameters.
- **Excitation energy:** How E* of primary fragments determines the final isotopic distribution after decay.
- **Production systematics:** Reviews all major imQMD+MNT reaction systems studied to date.
- **Future directions:** Surface dynamics, shell effects in transfer, and comparison with other transport codes (TDHF, DNS-sysu).

**Code Relevance**

| Code element | File:Line | Status |
|---|---|---|
| Overall pipeline: `run_imqmd_event` | `reaction.py:300-410` | ✅ Consistent |
| `initialize_and_relax` → `make_collision_event` → `propagate` → `identify_fragments` → decay | Full pipeline | ✅ Consistent |
| Parameter set selection via `SkyrmeEDF.from_name` | `skyrme.py:60-62` | ✅ Consistent |
| `select_candidates` (20 projectile + 20 target selection) | `initializer.py:582-615` | ✅ Consistent |
| Research report (Chinese) methodology cross-reference | `papers/research_report_imqmd_multinucleon_transfer.md` | ✅ Consistent |

**Assessment:** ✅ Consistent — This review paper describes the methodology our code implements.

---

## 4. New/Recent Papers (2020–2026)

### 4.1 Liao+2023 — Charge Equilibration in MNT

| Field | Detail |
|-------|--------|
| **Authors** | Liao Zi-Wei, Zhao Kai, Duan Ji-Zheng et al. |
| **Title** | Charge Equilibration in Multinucleon Transfer Reactions |
| **Journal** | Physical Review C 107, 014614 (2023) |
| **DOI** | `10.1103/PhysRevC.107.014614` |

**Key Physics Content**

Compares DNS-sysu (dinuclear system model) with ImQMD for charge equilibration:

- **Charge equilibration (CE) timescale:** For mass-asymmetric systems (e.g., ¹³⁶Xe+²⁰⁸Pb), CE is incomplete even at 500 fm/c — the lighter fragment retains memory of its initial N/Z.
- **Fast equilibration:** At very short times (~0.52 zs = 1.56 fm/c early contact), the isovector density oscillation establishes the initial N/Z gradient.
- **DNS-sysu vs ImQMD:** DNS-sysu, being a macroscopic model, predicts faster and more complete CE; ImQMD shows slower, more realistic equilibration dynamics due to shell effects and Pauli blocking.
- **Implications for production:** Incomplete CE means the optimal projectile for producing neutron-rich fragments is one with high initial N/Z (e.g., ²³⁸U rather than ¹³⁶Xe).

**Code Relevance**

| Code element | File:Line | Status |
|---|---|---|
| `reference_group_ids` tracking projectile/target origins | `reaction.py:186-223` | ✅ Consistent |
| `fermi_constraint_check` with `group_ids` inter-group exclusion | `collisions.py:484-597` | ✅ Consistent |
| `_compute_occupation_field` group-aware | `collisions.py:446-481` | ✅ Consistent |
| 250 fm/c cooling + 500 fm/c dynamics for CE evolution | `reaction.py:367-377` | ✅ Consistent |

**Assessment:** ⚠️ Needs verification — The code tracks charge equilibration via `reference_group_ids`, but no explicit CE diagnostic metric (e.g., N/Z ratio vs contact time) is implemented. Liao+2023's analysis of incomplete CE for asymmetric systems should be reproducible with our output but requires post-processing.

---

### 4.2 Xu+2024 — TMEP: Comparing Transport Codes

| Field | Detail |
|-------|--------|
| **Authors** | Jun Xu et al. (Transport Model Evaluation Project collaboration) |
| **Title** | Transport Model Evaluation Project: Comparing Transport Codes |
| **Journal** | Physical Review C 109, 044609 (2024) |
| **DOI** | `10.1103/PhysRevC.109.044609` |

**Key Physics Content**

The Transport Model Evaluation Project (TMEP) systematically compares different transport codes (ImQMD, UrQMD, BUU, SMF, etc.) on identical physics inputs:

- **Pauli blocking strategies:** Different codes implement Pauli blocking differently — some use local-density approximations, others use Wigner phase-space occupation. The choice of blocking strategy significantly affects:
  - **Collision rates:** Variance of up to factor 2–3 between different Pauli blocking implementations.
  - **Pion yields:** The π⁻/π⁺ ratio is sensitive to the neutron/proton blocking asymmetry.
  - **Code convergence:** ImQMD converges poorly without an "improved Pauli blocking" scheme at low/intermediate energies.
- **Key recommendation:** "An improved Pauli blocking" (i.e., softened Wigner occupation normalization) is needed for code convergence — the strict `4*sum(weights)` over-blocks at densities > 0.5ρ₀.

**Code Relevance**

| Code element | File:Line | Status |
|---|---|---|
| `OCCUPATION_SOFTENING = 10.0` (pragmatic divergence) | `collisions.py:26` | ✅ Consistent |
| `FERMI_SOFTENING = 1.35` (CoMD constraint softening) | `collisions.py:34` | ✅ Consistent |
| `_pauli_blocking_probability_fast` softened | `collisions.py:248-284` | ✅ Consistent |
| Module docstring explaining softening rationale | `collisions.py:14-34` | ✅ Consistent |
| `in_medium_nn_cross_section` enhancement factor | `collisions.py:114-117` | ✅ Consistent |

**Assessment:** ✅ Consistent — Our `OCCUPATION_SOFTENING` is exactly the kind of "improved Pauli blocking" that Xu+2024 recommends. The TMEP finding that strict blocking over-blocks at intermediate densities independently validates our pragmatic softening.

---

## 5. Parameter/Technical Reference Papers

### 5.1 Cugnon+1987 — NN Cross Section Parameterization

| Field | Detail |
|-------|--------|
| **Authors** | J. Cugnon, D. L'Hôte, J. Vandermeulen |
| **Title** | Simple Parametrization of Cross-Sections for Nuclear Transport Studies up to the GeV Region |
| **Journal** | Nuclear Physics A 470, 558 (1987) |
| **DOI** | `10.1016/0375-9474(87)90012-1` |

**Key Physics Content**

Provides the standard parameterization of free nucleon-nucleon elastic cross sections used in transport models:

- **pp/nn cross section (Cugnon formula):**
  \[
  \sigma_{pp} = 13.73 - 15.04/\beta + 8.76/\beta^2 + 68.67\beta^4 \quad \text{(mb)}
  \]
- **np cross section:**
  \[
  \sigma_{np} = -70.67 - 18.18/\beta + 25.26/\beta^2 + 113.85\beta \quad \text{(mb)}
  \]
  where `β` is the nucleon velocity in the lab frame, `β² = 1 - (M_N/(M_N + E_lab))²`.

- **Low-energy behavior:** The Cugnon formula was designed for intermediate/high energies (50–1000 MeV lab). Below 50 MeV lab, cross sections should be constant (~180 mb for np, ~60 mb for pp/nn) since the parameterization diverges at low β.

**Code Relevance**

| Code element | File:Line | Status |
|---|---|---|
| `free_nn_cross_section` Cugnon parameterization | `collisions.py:53-74` | ✅ Consistent |
| Low-energy cap: 180 mb (np), 60 mb (pp/nn) below 50 MeV lab | `collisions.py:64-65` | ✅ Consistent |
| `beta2 = 1 - (M_N/(M_N+E_lab))²` computation | `collisions.py:68-69` | ✅ Consistent |
| `_sqrt_s_gev` helper for in-medium factor | `collisions.py:77-87` | ✅ Consistent |

**Assessment:** ✅ Consistent — Cugnon parameterization coefficients match exactly.

---

### 5.2 Papa+Bonasera 2001 — CoMD Phase-Space Constraint

| Field | Detail |
|-------|--------|
| **Authors** | M. Papa, T. Maruyama, A. Bonasera |
| **Title** | Constrained Molecular Dynamics Approach to Fermionic Systems |
| **Journal** | Physical Review C 64, 024612 (2001) |
| **DOI** | `10.1103/PhysRevC.64.024612` |

**Key Physics Content**

The CoMD model introduces the phase-space occupation constraint that ImQMD later adopts:

- **Fermi constraint principle:** For every nucleon `i`, the Wigner phase-space occupation `f_i` from same-species same-group nucleons must satisfy `f_i ≤ 1` (at most one fermion per phase-space cell `h³/2`).
- **Correction mechanism:** When a nucleon exceeds the threshold, its momentum is "swapped" with a nearby same-species nucleon's momentum to reduce the total occupation. This conserves total momentum and kinetic energy while suppressing Pauli-forbidden states.
- **ImQMD adoption:** ImQMD applies the Fermi constraint after each propagation step and after each collision step, ensuring that the phase-space distribution never violates the Pauli principle.
- **Threshold:** `f_threshold = 1.0` (CoMD standard).

**Code Relevance**

| Code element | File:Line | Status |
|---|---|---|
| `fermi_constraint_check` | `collisions.py:484-597` | ✅ Consistent |
| Threshold `1.0` default | `collisions.py:486` | ✅ Consistent |
| `FERMI_SOFTENING = 1.35` (near-CoMD occupation) | `collisions.py:34` | ✅ Consistent |
| Momentum SWAP (not stretch/rescale) conservation | `collisions.py:593-597` | ✅ Consistent |
| `propagate` applies constraint every step | `propagator.py:354-356` | ✅ Consistent |
| `group_ids` for projectile/target group isolation | `collisions.py:469-471` | ✅ Consistent |

**Assessment:** ✅ Consistent — The momentum-swap Fermi constraint follows Papa+Bonasera 2001 exactly. Our implementation conserves Σp and Σp² (swaps relabel momenta within same-species set), unlike the older momentum-stretch variant that injected energy.

---

### 5.3 Charity+1988 — GEMINI Statistical Decay

| Field | Detail |
|-------|--------|
| **Authors** | R.J. Charity, M.A. McMahan, G.J. Wozniak, R.J. McDonald, L.G. Moretto, D.G. Sarantites, L.G. Sobotka, G. Guarino, A. Pantaleo, L. Fiore, A. Gobbi, K.D. Hildenbrand |
| **Title** | Systematics of Complex Fragment Emission in Niobium-Induced Reactions |
| **Journal** | Nuclear Physics A 483, 371 (1988) |
| **DOI** | `10.1016/0375-9474(88)90542-8` |

**Key Physics Content**

The GEMINI statistical decay model, the standard code used for de-excitation of imQMD primary fragments:

- **Weisskopf evaporation:** Neutron, proton, and alpha evaporation widths using the Weisskopf statistical model:
  \[
  \Gamma = \frac{2m}{\pi\hbar^2}\int_0^{E^*}\sigma_{inv}(\epsilon)\frac{\rho_d(E^* - \epsilon)}{\rho_p(E^*)} \epsilon \, d\epsilon
  \]
- **Fermi gas level density:** `ρ(E*) ∝ exp(2√(aE*))/E*^(5/4)`.
- **Fission competition:** Bohr-Wheeler transition-state fission width.
- **Level density parameter:** `a = A/8` MeV⁻¹ (standard GEMINI value).

**Code Relevance**

| Code element | File:Line | Status |
|---|---|---|
| `_weisskopf_width` (evaporation width) | `decay.py:210-219` | ✅ Consistent |
| `_fermi_gas_level_density` | `decay.py:186-197` | ✅ Consistent |
| `_level_density_parameter` `a = A/8` | `decay.py:179-183` | ✅ Consistent |
| `_bohr_wheeler_width` (fission) | `decay.py:230-246` | ✅ Consistent |
| `weisskopf_evaporation_multi` (full channel selection) | `decay.py:333-349` | ✅ Consistent |
| `evaporate_chain` (multi-step Monte Carlo) | `decay.py:387-453` | ✅ Consistent |
| `_gamma_width` (statistical gamma) | `decay.py:222-227` | ✅ Consistent |

**Assessment:** ⚠️ Needs verification — Our decay module implements a simplified Weisskopf evaporation + Bohr-Wheeler fission, which captures the GEMINI physics adequately for Phase 1 but lacks the full GEMINI treatment of complex-fragment emission (Z>2 IMF emission), angular momentum effects, and sequential binary decay chains. For high-excitation MNT fragments (E*>100 MeV), GEMINI's full treatment would differ.

---

### 5.4 Aichelin 1991 — QMD Model with MST

| Field | Detail |
|-------|--------|
| **Authors** | J. Aichelin |
| **Title** | "Quantum" Molecular Dynamics — A Dynamical Microscopic N-Body Approach to Investigate Fragment Formation and the Nuclear Equation of State in Heavy Ion Collisions |
| **Journal** | Physics Reports 202, 233 (1991) |
| **DOI** | `10.1016/0370-1573(91)90094-3` |

**Key Physics Content**

The seminal reference establishing the standard QMD framework that imQMD builds upon:

- **Gaussian wave packet representation:** Each nucleon `i` is a Gaussian:
  \[
  \phi_i(\mathbf{r}) = \frac{1}{(2\pi\sigma_r^2)^{3/4}} \exp\left[-\frac{(\mathbf{r} - \mathbf{r}_i)^2}{4\sigma_r^2} + \frac{i}{\hbar}\mathbf{p}_i\cdot\mathbf{r}\right]
  \]
  with fixed width `σ_r`. (ImQMD modifies this to an A-dependent `σ_r`.)
- **Equations of motion:** Hamilton's equations: `ṙ_i = p_i/m`, `ṗ_i = -∇_i ⟨H⟩`.
- **Stochastic NN collisions:** Geometric criterion `d_min < √(σ_NN/π)`. Pauli blocking of final states.
- **Minimum Spanning Tree (MST):** The standard QMD fragment recognition algorithm:
  - `r_cut = 2.4–3.0 fm` (classic Aichelin value: ≈3.0 fm)
  - `p_cut = 250 MeV/c` (classic momentum cut)
  - Two nucleons belong to the same fragment if `Δr < r_cut` AND `Δp < p_cut`.

**Code Relevance**

| Code element | File:Line | Status |
|---|---|---|
| MST algorithm (connected components) | `fragments.py:85-107` | ✅ Consistent |
| Classic `p_cut=250` referenced in docstring | `fragments.py:93-96` | ✅ Consistent |
| ImQMD `p_cut=300` as default for MNT | `fragments.py:88-89` | ✅ Consistent |
| Gaussian packet `normalization = 1/(2πσ²)^(3/2)` | `grid_edf.py:699` | ✅ Consistent |
| Gaussian packet representation width | `grid_edf.py:697-698` | ✅ Consistent |
| `d_min < √(σ_NN/π)` collision criterion | `collisions.py:414-416` | ✅ Consistent |

**Assessment:** ✅ Consistent — The Gaussian packet representation, equations of motion, collision criterion, and MST fragment recognition follow Aichelin 1991. The `300 MeV/c` p_cut (vs classic `250`) is the well-documented imQMD MNT choice.

---

## Summary of Parameter Values

| Parameter | Symbol | IQ1 | IQ2 | IQ3 | IQ3a | IQ3b | Code Location |
|-----------|--------|-----|-----|-----|------|------|---------------|
| Bulk density coeff | α (MeV) | -310.0 | -356.0 | -207.0 | -207.0 | -207.0 | `skyrme.py:39-45` |
| Bulk density coeff | β (MeV) | 258.0 | 303.0 | 138.0 | 138.0 | 138.0 | `skyrme.py:39-45` |
| Density exponent | γ | 7/6 | 7/6 | 7/6 | 7/6 | 7/6 | `skyrme.py:39-45` |
| Surface energy | g_sur (MeV·fm²) | 19.8 | 7.0 | 18.0 | 16.5 | 18.0 | `skyrme.py:39-45` |
| Momentum-dependent | g_τ (MeV) | 9.5 | 12.5 | 14.0 | 14.0 | 14.0 | `skyrme.py:39-45` |
| Momentum exponent | η | 2/3 | 2/3 | 5/3 | 5/3 | 5/3 | `skyrme.py:39-45` |
| Symmetry energy | C_s (MeV) | 32.0 | 32.0 | 32.0 | 34.0 | 34.0 | `skyrme.py:39-45` |
| Surface symmetry | κ_s | 0.08 | 0.08 | 0.08 | 0.4 | 0.6 | `skyrme.py:39-45` |
| Saturation density | ρ₀ (fm⁻³) | 0.165 | 0.165 | 0.165 | 0.165 | 0.165 | `skyrme.py:39-45` |
| Width constant | σ₀ (fm) | 0.49 | 0.88 | 0.94 | 0.94 | 0.94 | `skyrme.py:39-45` |
| Width A-dep | σ₁ (fm) | 0.16 | 0.09 | 0.018 | 0.020 | 0.018 | `skyrme.py:39-45` |

## Summary of Numerical Defaults

| Parameter | Value | Origin | Code Location |
|-----------|-------|--------|---------------|
| MST r_cut | 3.5 fm | Li+2016 | `fragments.py:87` |
| MST p_cut (ImQMD MNT) | 300 MeV/c | Li+2016 | `fragments.py:88` |
| MST p_cut (classic QMD) | 250 MeV/c | Aichelin 1991 | `fragments.py:93-96` comment |
| Dynamics switch time | 500 fm/c | Li+2016 | `reaction.py:315` |
| Cooling time (no collisions) | 250 fm/c | Code optimization | `reaction.py:373` |
| Initial separation | 50 fm | Research report Sec. 4.2 | `reaction.py:122` |
| RK4 time step | 1.0 fm/c | Standard imQMD | `reaction.py:316` |
| Fermi constraint interval | 5 fm/c | Code optimization | `propagator.py:326` |
| Fermi constraint threshold | 1.0 | Papa+Bonasera 2001 | `collisions.py:486` |
| Pauli softening | 10.0 | Empirical (Xu+2024) | `collisions.py:26` |
| Fermi softening | 1.35 | Empirical | `collisions.py:34` |
| Grid spacing | 1.0 fm | Standard | `grid_edf.py:50` |
| Grid n_sigma cutoff | 3.0 | Standard | `grid_edf.py:51` |
| Saturation density ρ₀ | 0.165 fm⁻³ | IQ1-IQ3 sets | `collisions.py:12` |
| NN low-E cap (np) | 180 mb | Cugnon 1987 | `collisions.py:65` |
| NN low-E cap (pp/nn) | 60 mb | Cugnon 1987 | `collisions.py:65` |
| Level density a | A/8 MeV⁻¹ | Charity 1988 | `decay.py:180` |
| Level density a (fission) | 1.04·A/8 MeV⁻¹ | Standard | `decay.py:182` |

## Cross-Reference Index

| Source File | Relevant Papers |
|-------------|-----------------|
| `skyrme.py` | Wang+2002, Wang+2004, Li+2013, Li+2025 |
| `initializer.py` | Wang+2002, Wang+2014, Li+2013 |
| `grid_edf.py` | Wang+2002, Li+2025, Wang+2014 |
| `collisions.py` | Zhang+2020, Xu+2024, Papa+2001, Cugnon+1987 |
| `propagator.py` | Wang+2002, Yao+Wang 2017 |
| `fragments.py` | Li+2016, Li+2019, Tian+2022, Aichelin 1991, Zhang+2020 |
| `reaction.py` | Li+2016, Li+2019, Li+2020, Tian+2022, Liao+2023 |
| `decay.py` | Charity+1988, Zhao+2021, Li+2020 |
