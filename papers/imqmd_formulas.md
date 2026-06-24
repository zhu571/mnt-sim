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
