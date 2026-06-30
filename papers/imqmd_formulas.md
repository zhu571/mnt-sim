# ImQMD Phase 1 Formula Notes

Sources:

- Wang, Li, Wu, "An improved quantum molecular dynamics model and its applications to fusion reaction near barrier", Phys. Rev. C 65, 064608 (2002), arXiv:nucl-th/0201079.
- Wang, Ou, Zhang, Li, "Microscopic dynamics simulations of heavy-ion fusion reactions induced by neutron-rich nuclei", Phys. Rev. C 89, 064601 (2014), arXiv:1405.5271.
- Li, Tian, Qin, Li, Wang, "Determination of the nucleon-nucleon interaction in the ImQMD model by nuclear reactions at the Fermi energy region", Chinese Phys. C 37, 114101 (2013), for IQ1/IQ2/IQ3 parameter table.

## 1. Gaussian wave packet and density

Wang et al. 2014 Eq. (1) writes the one-body density as

```text
rho(r) = sum_i [1 / (2*pi*sigma_r^2)^(3/2)]
         exp[-(r - r_i)^2 / (2*sigma_r^2)].
```

The nucleon centroids are `r_i` and `p_i`.  `sigma_r` is the coordinate-space wave-packet width.  Later ImQMD work often uses a system-size dependent width

```text
sigma_r^n = sigma0 + sigma1 * A_n^(1/3),  n = projectile,target
```

from the Chinese Phys. C 2013 Eq. (3).  For simple static nuclei in this Phase 1 implementation the default fixed value is `sigma_r = 1.1 fm`, in the typical 1.0-1.2 fm range requested by the task.

## 2. Hamiltonian

Wang et al. 2014 Eqs. (2)-(4):

```text
dot(r_i) =  dH / d p_i
dot(p_i) = -dH / d r_i

H = T + U
T = sum_i p_i^2 / (2m)
U = U_loc + U_Coul
U_loc = integral V_loc(r) dr.
```

The Coulomb term contains direct and exchange contributions in the paper.  The code uses a finite-Gaussian direct Coulomb pair term for protons.

## 3. Local Skyrme energy-density functional

Wang et al. 2014 Eq. (5), without spin-orbit:

```text
V_loc =
  alpha/2 * rho^2/rho0
+ beta/(gamma+1) * rho^(gamma+1)/rho0^gamma
+ g_sur/(2*rho0) * (grad rho)^2
+ g_tau * rho^(eta+1)/rho0^eta
+ C_s/(2*rho0) * [rho^2 - k_s*(grad rho)^2] * delta^2

delta = (rho_n - rho_p) / (rho_n + rho_p).
```

Symmetry energy in cold matter, Wang et al. 2014 Eq. (6):

```text
E_sym(rho) = 13 * (rho/rho0)^(2/3) + C_s/2 * (rho/rho0).
```

The slope parameter is Wang et al. 2014 Eq. (7):

```text
L = 3*rho0 * [d E_sym / d rho]_(rho=rho0).
```

The original 2002 paper and later ImQMD variants may include momentum-dependent terms for intermediate-energy heavy-ion collisions.  Wang et al. 2014 notes this extension, but for near-barrier fusion and this static-nucleus phase the collision term and explicit momentum-dependent interaction are not included.

## 4. Ground-state initialization

Original 2002 preparation:

1. Obtain neutron/proton density distributions, then sample nucleon positions.
2. Compute local Fermi momentum from the local density.
3. Sample momenta inside a reduced local Fermi sphere to account for packet momentum width, with the uncertainty relation `sigma_r*sigma_p = hbar/2`.
4. Evolve candidate nuclei and retain only stable nuclei with acceptable binding energy, radius, density, momentum, and no spurious emission.

Wang et al. 2014 Sec. II.B details the hard-sphere version:

```text
R_c = 1.226*A^(1/3) + 2.86*A^(-2/3) - 1.09*(I - I^2) + 0.99*DeltaE/A       (8)
Delta R_np = <r_n^2>^(1/2) - <r_p^2>^(1/2) = 0.9*I - 0.03                 (9)
```

The paper samples protons inside a hard sphere of radius `R_p - w_r` and neutrons inside `R_n - w_r`, with `w_r = 0.8 fm`.  This Phase 1 implementation uses the requested simpler radius

```text
R0 = 1.14*A^(1/3)
```

and then subtracts the center-of-mass coordinate.

For momenta, Wang et al. 2014 samples the `i`-th nucleon from a local Fermi sphere:

```text
p_F,q(r_i) = hbar * [3*pi^2*rho_q(r_i)]^(1/3) - w_p,  q = n,p.
```

The sampled nucleus is accepted in the paper if:

```text
E_ground in BE +/- 0.05 MeV
|r_i - r_j| * |p_i - p_j| >= 255 fm MeV/c for all pairs.
```

For the requested Phase 1 tests the acceptance tolerance is relaxed to about `2 MeV/nucleon`; the initialized object stores a calibrated energy offset so the static total energy matches the accepted empirical binding.

## 5. Equations of motion

Wang et al. 2014 Eq. (2):

```text
dot(r_i) =  dH / d p_i
dot(p_i) = -dH / d r_i.
```

The implementation uses fourth-order Runge-Kutta and analytic gradients for the finite-Gaussian Coulomb plus a weak calibrated static mean field.

## 6. IQ2 and IQ3 parameter values

Chinese Phys. C 2013 Table 1 gives:

```text
set   alpha  beta  gamma  g0/gsur  g_tau  eta  C_s  kappa_s  rho0   sigma0  sigma1
      MeV    MeV          MeV fm2   MeV         MeV  fm2      fm^-3  fm      fm
IQ1   -310   258   7/6    19.8      9.5    2/3  32   0.08     0.165  0.49    0.16
IQ2   -356   303   7/6    7.0       12.5   2/3  32   0.08     0.165  0.88    0.09
IQ3   -207   138   7/6    18.0      14.0   5/3  32   0.08     0.165  0.94    0.018
```

Wang et al. 2014 Table I gives updated neutron-rich fusion sets:

```text
set    alpha  beta  gamma  gsur  g_tau  eta  C_s  kappa_s  rho0   sigma0  sigma1
SkP*   -356   303   7/6    19.5  13     2/3  35   0.65     0.162  0.94    0.018
IQ3a   -207   138   7/6    16.5  14     5/3  34   0.4      0.165  0.94    0.020
IQ3b   -207   138   7/6    18.0  14     5/3  34   0.6      0.165  0.94    0.018
```

This code defaults to IQ2 because the task requested IQ2, and exposes IQ3/IQ3a/IQ3b through `SkyrmeEDF.from_name`.

## 7. Fermi constraint

The 2014 paper adopts a modified Fermi constraint following CoMD and checks total energy after two-body elastic scattering.  Its practical acceptance condition is stated as:

```text
|r_i - r_j| * |p_i - p_j| >= 255 fm MeV/c.
```

This Phase 1 initializer applies this as a candidate quality check and resamples momenta when the minimum phase-space distance is too small.

## Phase 2: Collisions

Sources:

- Zhang et al., Phys. Rev. C 85, 024602 (2012), arXiv:1009.1928.
- Chen, Zhang, Li, Chinese Phys. C 45, 074106 (2021), arXiv:2103.13218.
- Chen et al., Phys. Rev. C 109, 034611 (2024), arXiv:2403.00343.
- Papa et al., CoMD Fermi constraint idea as used by later ImQMD initialization/constraint checks.

### 1. Geometric NN collision criterion

For each pair over one collision step `dt`, estimate the closest approach of the two packet centroids using their relative coordinate and velocity:

```text
r_ij(t) = r_i - r_j + (v_i - v_j) t
t_min = clamp[-r_ij(0) dot v_ij / |v_ij|^2, 0, dt]
d_min = |r_ij(t_min)|
```

An attempted collision is made when

```text
d_min < sqrt(sigma_NN_med / pi),
```

with `sigma` converted from mb to fm^2 by `1 mb = 0.1 fm^2`.  The elastic final state conserves the pair center-of-mass momentum and relative momentum magnitude, while the outgoing relative direction is sampled isotropically.

The pair center-of-mass kinetic energy used by the compact implementation is

```text
E_cm = q^2 / m_N,    q = |p_i - p_j| / 2.
```

### 2. Free NN cross sections

Zhang et al. state that ImQMD05 uses isospin-dependent free NN cross sections from Cugnon and applies an in-medium factor.  Chen-Zhang-Li 2021 also notes the low-energy ImQMD caps

```text
sigma_nn/pp_free = 60 mb,
sigma_np_free    = 180 mb,       for p_lab < 0.3 GeV/c, roughly E_lab < 50 MeV.
```

Above that low-energy cap, the code uses a compact Cugnon-style beta parameterization:

```text
sigma_pp/nn = 13.73 - 15.04/beta + 8.76/beta^2 + 68.67 beta^4
sigma_np    = -70.67 - 18.18/beta + 25.26/beta^2 + 113.85 beta
beta^2      = 1 - [m_N / (m_N + E_lab)]^2
E_lab       = 2 E_cm
```

The parameterization is clipped to a positive finite range for numerical robustness.

### 3. In-medium scaling

Zhang et al. 2012 Sec. II uses

```text
sigma*_nn/np = (1 - xi(E_beam) rho/rho0) sigma_free_nn/np
xi(E_beam = 50 AMeV) = 0.2.
```

The Phase 2 implementation follows the requested common-channel form

```text
sigma_NN_med = f_med sigma_NN_free
f_med        = 1 - eta rho/rho0,   eta = 0.2,
rho0         = 0.16 fm^-3.
```

### 4. Pauli blocking

For an attempted collision `i + j -> i' + j'`, positions are unchanged and final momenta are tested.  Chen et al. 2024 gives the Uehling-Uhlenbeck blocking form

```text
P_block = 1 - (1 - P_i)(1 - P_j),
```

where `P_i` and `P_j` are the occupation probabilities of the two final states.

The standard ImQMD/QMD Wigner occupation method in Chen-Zhang-Li 2021 and Chen 2024 evaluates same-isospin neighbors around the outgoing state:

```text
P_i = P(r_i, p'_i)
    proportional to sum_{k != i, tau_k = tau_i}
      exp[-(r_i - R_k)^2 / (2 sigma_r^2)]
      exp[-(p'_i - P_k)^2 / (2 sigma_p^2)]

sigma_r sigma_p = hbar/2,
P_i -> min(P_i, 1).
```

Chen-Zhang-Li 2021 compares Wigner, Husimi, and hard-sphere overlap algorithms and defines the blocking ratio

```text
R_block = 1 - (dN_coll^suc/dt) / (dN_coll^att/dt).
```

Chen et al. 2024 proposes a smoother PB(W*) occupation by replacing the single momentum state of each neighbor with an average over sampled momentum states:

```text
c_j(p'_i) = (1/N) sum_lambda exp[-(p'_i - P_{j,lambda})^2 / (2 sigma_p^2)]
P_i       proportional to sum_{j != i} exp[-(r_i - R_j)^2/(2 sigma_r^2)] c_j(p'_i).
```

The code exposes this as an optional broadened-kernel blocker while using the standard Wigner blocker by default.

### 5. Fermi / CoMD constraint

The CoMD-style constraint enforces a maximum allowed phase-space occupation for identical nucleons.  In the simplified implementation this is represented by the same practical phase-space separation used during Phase 1 initialization:

```text
|r_i - r_j| |p_i - p_j| >= 255 fm MeV/c,     tau_i = tau_j.
```

When an identical pair violates the threshold, the pair center-of-mass momentum is preserved and the relative momentum is reshuffled.  To avoid unphysical heating of an already accepted cold nucleus, the implementation preserves the pair kinetic energy whenever the relative momentum is nonzero, and applies only a small separating kick for exactly degenerate momenta.  This is applied periodically, about every `20 fm/c`, during propagation with collisions enabled.

### 6. Collision scheduling and rate

Chen-Zhang-Li 2021 gives the analytical attempted collision rate in uniform matter as

```text
<dN_coll^att/dt> = (1/2) A rho <v_rel sigma_NN_med>.
```

The code schedules geometric stochastic attempts every `collision_dt` (default: every propagation step), so its pair attempts are controlled by the same density, relative velocity, and in-medium cross-section dependence.  At low energies in a cold finite nucleus, most attempted collisions should be Pauli-blocked and the successful collision rate should be near zero.

## Phase 3: Fragment Recognition and De-excitation

Sources:

- Zhang, Li, Zhou, Tsang, "Effect of isospin dependent cluster recognition on the observables in heavy ion collisions", arXiv:1205.1605.
- Wei, Wang, Ou, "Mechanism of production of light complex particles in nucleon-induced reactions", arXiv:1309.7534.
- Zhao et al., "The production of unknown neutron-rich isotopes in 238U+238U collisions at near-barrier energy", arXiv:1605.07393.
- Wang et al., "Further Development of the Improved QMD Model and its Applications to Fusion Reaction near Barrier", arXiv:nucl-th/0402066.

### 1. Minimum spanning tree fragment recognition

QMD calculations identify primary fragments after the dynamical stage with a minimum spanning tree (MST) coalescence rule.  Nucleons `i` and `j` belong to the same fragment if they are connected, directly or through neighbors, by

```text
|r_i - r_j| <= R_cut
|p_i - p_j| <= P_cut.
```

Zhang 2012 summarizes common QMD values as `R_cut` of order the nucleon interaction range, about `3-3.5 fm`, and `P_cut` around `250 MeV/c`.  Wei 2014 uses `R_c = 4.5 fm` and `P_c = 250 MeV/c` for a spallation application.  For near-barrier ImQMD fragment tagging the coordinate cut alone is often the most robust default; in code this corresponds to `p_cut = None`.

### 2. Isospin-dependent MST

The ordinary MST treats all nucleon pairs with one common distance threshold:

```text
R0_nn = R0_np = R0_pp = R0.
```

Zhang 2012 introduces an isospin-dependent MST by using different coordinate cutoffs for neutron-neutron, neutron-proton, and proton-proton links:

```text
|r_i - r_j| <= R0_tau_i_tau_j
|p_i - p_j| <= P0.
```

Their sensitivity study used a deliberately strong example,

```text
R0_nn = R0_np = 6 fm
R0_pp = 3 fm
P0    = 250 MeV/c,
```

motivated by neutron skins/halos and proton Coulomb repulsion.  The implementation keeps conservative defaults near the requested Phase 3 scale, e.g. `R_pp = R_nn = 2.8 fm`, `R_np = 3.2 fm`, with the stronger Zhang values available through arguments.

### 3. Fragment excitation energy

After recognizing primary fragments, Zhao 2016 and Wei 2014 compute each fragment's total energy in its own rest frame and subtract the corresponding ground-state energy:

```text
E*(Z,A) = E_int(fragment) - E_ground(Z,A).
```

Operationally,

```text
E_int = sum_{i in F} (p_i - P_F/A)^2 / (2 m_N) + U_F
E_ground(Z,A) = -B_ground(Z,A),
```

where `P_F` is the fragment center-of-mass momentum and `U_F` is the fragment potential energy from the same compact EDF/Coulomb model used by the ImQMD nucleus object.  The code estimates `B_ground` from the HIVAP `Mexcess95.dat` mass table when available and otherwise from a liquid-drop binding formula.

### 4. Statistical decay and fission competition

Zhao 2016 couples ImQMD primary fragments to HIVAP.  HIVAP evaluates survival with branching ratios from relative decay widths,

```text
Gamma_i(Z,A,E*) / Gamma_tot(Z,A,E*),
Gamma_tot = sum_i Gamma_i,   i = gamma, n, p, alpha, fission.
```

The Phase 3 compact decay module implements the neutron-evaporation part with a Weisskopf-like statistical weight,

```text
P(n) proportional to exp[-n B_n / T] * (E_rem / E*)^(a/2)
E_rem = E* - n B_n
a     = A/8 MeV^-1
T     = sqrt(E*/a),
```

with `B_n` from the mass table or liquid-drop fallback.  A simple Bohr-Wheeler fission competition factor is used for heavy systems,

```text
P_fiss = Gamma_f / (Gamma_f + Gamma_evap)
       approx logistic[exp(-B_fiss/T)].
```

The fission barrier is a liquid-drop-scale estimate that decreases with fissility, and fission is forced to dominate for superheavy fragments (`Z >= 100`) at sizable excitation.

### 5. Light-cluster coalescence

Wei 2014 introduces direct light complex particle formation for

```text
d  = 1p + 1n
t  = 1p + 2n
3He = 2p + 1n
4He = 2p + 2n.
```

Candidate nucleons are required to be close in phase space, expressed there as a Jacobian-coordinate product cut,

```text
R_im * P_im <= h0,     R_im >= 1 fm,
```

with priority `4He > 3He > t > d`.  The compact implementation uses equivalent tight pairwise cuts in coordinate and momentum space for final-state recognition, also applying the same heavy-to-light priority so a nucleon is not assigned to multiple light clusters.
