# Complete ImQMD Model Specification from the Literature

This document specifies the ImQMD model family as described in the supplied literature, with emphasis on the original ImQMD model, the 2014 neutron-rich fusion initialization/EDF paper, the ImQMD05 collision and fragment-recognition papers, and the 2016 U+U MNT plus HIVAP application. It is intended as an audit target for a Python implementation.

Primary sources:

- Wang, Li, Wu, PRC 65, 064608 (2002), `nucl-th_0201079.pdf`: original ImQMD.
- Wang, Ou, Zhang, Li, PRC 89, 064601 (2014), `1405.5271.pdf`: neutron-rich initialization and IQ3a/IQ3b fusion parameters.
- Zhang et al., PRC 85, 024602 (2012), `1009.1928.pdf`: ImQMD05 collision term, in-medium NN cross sections, MST settings.
- Zhang, Li, Zhou, Tsang, PRC 85, 051602(R) (2012), `1205.1605.pdf`: isospin-dependent cluster recognition.
- Zhao et al., PRC 94, 024601 (2016), `1605.07393.pdf`: U+U MNT application, primary-fragment cross sections, HIVAP coupling.
- Zhang et al., Front. Phys. 15, 54301 (2020), `2005.12877.pdf`: review and model-family overview.
- Chen et al., PRC 109, L021604 (2024), `2403.00343.pdf`: improved Pauli blocking.
- Wei, Wang, Ou, JPG 41, 035104 (2014), `1309.7534.pdf`: optional light-complex-particle surface coalescence.
- Li, Tian, Qin, Li, Wang, Chinese Phys. C 37, 114101 (2013): IQ1/IQ2/IQ3 parameter table, as recorded in `papers/imqmd_formulas.md`.

## 1. Gaussian Wave Packet Representation

The original ImQMD model represents each nucleon by a coherent Gaussian wave packet. Wang 2002 Eq. (1) gives

$$
\phi_i(\mathbf r)=\frac{1}{(2\pi\sigma_r^2)^{3/4}}
\exp\left[-\frac{(\mathbf r-\mathbf r_i)^2}{4\sigma_r^2}
+\frac{i}{\hbar}\mathbf r\cdot \mathbf p_i\right].
\tag{W2002-1}
$$

The many-body trial state is the direct product

$$
\Phi=\prod_i\phi_i.
$$

Antisymmetrization is not explicit; fermionic behavior is approximated by Pauli blocking in collisions and by a Fermi/phase-space constraint during initialization and propagation.

The single-particle Wigner distribution is Wang 2002 Eq. (2):

$$
f_i(\mathbf r,\mathbf p)=\frac{1}{(\pi\hbar)^3}
\exp\left[-\frac{(\mathbf r-\mathbf r_i)^2}{2\sigma_r^2}
-\frac{2\sigma_r^2}{\hbar^2}(\mathbf p-\mathbf p_i)^2\right].
\tag{W2002-2}
$$

The coordinate and momentum densities are Wang 2002 Eqs. (3)-(6):

$$
\rho(\mathbf r)=\sum_i \rho_i(\mathbf r),\qquad
g(\mathbf p)=\sum_i g_i(\mathbf p),
\tag{W2002-3,4}
$$

$$
\rho_i(\mathbf r)=\frac{1}{(2\pi\sigma_r^2)^{3/2}}
\exp\left[-\frac{(\mathbf r-\mathbf r_i)^2}{2\sigma_r^2}\right],
\tag{W2002-5}
$$

$$
g_i(\mathbf p)=\frac{1}{(2\pi\sigma_p^2)^{3/2}}
\exp\left[-\frac{(\mathbf p-\mathbf p_i)^2}{2\sigma_p^2}\right].
\tag{W2002-6}
$$

The widths satisfy the minimum uncertainty relation, Wang 2002 Eq. (7):

$$
\sigma_r\sigma_p=\frac{\hbar}{2}.
\tag{W2002-7}
$$

The dynamical variables are the centroid coordinates and momenta,

$$
\{(\mathbf r_i,\mathbf p_i)\}_{i=1}^A.
$$

### Width prescription

The original ImQMD paper found that a fixed width cannot describe light and heavy nuclei equally well. It proposed Wang 2002 Eq. (18):

$$
\sigma_r = 0.16 A_N^{1/3}+0.49\ {\rm fm},
\tag{W2002-18}
$$

where \(A_N\) is the number of nucleons bound in the current nucleus, cluster, or compound system. This is the IQ1 width prescription.

Later ImQMD parameter sets use the general form

$$
\sigma_r(A)=\sigma_0+\sigma_1 A^{1/3},
\tag{width}
$$

with \(\sigma_0,\sigma_1\) supplied by the parameter set. For reactions, the projectile and target can each be initialized with their own \(A^{1/3}\)-dependent width. When a single system-dependent width is used after contact, the literature treats \(A\) as the bound system size; an implementation must define this transition consistently.

## 2. Hamiltonian

The centroids propagate under a self-consistent Hamiltonian. Wang 2002 Eqs. (8)-(12) and Wang 2014 Eqs. (2)-(4) write

$$
\dot{\mathbf r}_i=\frac{\partial H}{\partial \mathbf p_i},
\qquad
\dot{\mathbf p}_i=-\frac{\partial H}{\partial \mathbf r_i},
\tag{W2014-2}
$$

$$
H=T+U,\qquad
T=\sum_i\frac{\mathbf p_i^2}{2m},
\tag{W2002-9,10}
$$

$$
U=U_{\rm loc}+U_{\rm Coul},\qquad
U_{\rm loc}=\int V_{\rm loc}(\mathbf r)\,d^3r.
\tag{W2014-3,4}
$$

### Modern local Skyrme EDF

For fusion-oriented ImQMD, Wang 2014 Eq. (5) gives the local energy-density functional without spin-orbit:

$$
V_{\rm loc}
=\frac{\alpha}{2}\frac{\rho^2}{\rho_0}
+\frac{\beta}{\gamma+1}\frac{\rho^{\gamma+1}}{\rho_0^\gamma}
+\frac{g_{\rm sur}}{2\rho_0}(\nabla\rho)^2
+g_\tau\frac{\rho^{\eta+1}}{\rho_0^\eta}
+\frac{C_s}{2\rho_0}\left[\rho^2-\kappa_s(\nabla\rho)^2\right]\delta^2,
\tag{W2014-5}
$$

where

$$
\rho=\rho_n+\rho_p,\qquad
\delta=\frac{\rho_n-\rho_p}{\rho_n+\rho_p}.
$$

Term-by-term audit expectations:

- Two-body bulk term:

$$
V_{t0}=\frac{\alpha}{2}\frac{\rho^2}{\rho_0}.
$$

- Density-dependent bulk term:

$$
V_{t3}=\frac{\beta}{\gamma+1}\frac{\rho^{\gamma+1}}{\rho_0^\gamma}.
$$

- Isoscalar surface term:

$$
V_{\rm sur}=\frac{g_{\rm sur}}{2\rho_0}(\nabla\rho)^2.
$$

- Momentum/Thomas-Fermi \(rho tau\)-like term in this fusion EDF:

$$
V_\tau=g_\tau\frac{\rho^{\eta+1}}{\rho_0^\eta}.
$$

- Symmetry and surface-symmetry term:

$$
V_{\rm sym}=
\frac{C_s}{2\rho_0}\left[\rho^2-\kappa_s(\nabla\rho)^2\right]\delta^2.
$$

The cold nuclear matter symmetry energy corresponding to the Wang 2014 EDF is Eq. (6):

$$
E_{\rm sym}(\rho)=13\left(\frac{\rho}{\rho_0}\right)^{2/3}
+\frac{C_s}{2}\left(\frac{\rho}{\rho_0}\right).
\tag{W2014-6}
$$

The slope parameter is Wang 2014 Eq. (7):

$$
L=3\rho_0\left(\frac{\partial E_{\rm sym}}{\partial\rho}\right)_{\rho=\rho_0}.
\tag{W2014-7}
$$

### Original 2002 EDF

The original model used the simpler Skyrme-like local functional, Wang 2002 Eq. (13):

$$
V_{\rm loc}
=\frac{\alpha}{2}\frac{\rho^2}{\rho_0}
+\frac{\beta}{3}\frac{\rho^3}{\rho_0^2}
+\frac{C_s}{2}\frac{(\rho_p-\rho_n)^2}{\rho_0}
+\frac{g_1}{2}(\nabla\rho)^2.
\tag{W2002-13}
$$

After the \(N^3\) density-dependent term is approximated by an \(N^2\) expression, Wang 2002 Eq. (16), the extra gradient contribution is absorbed into

$$
g_0=g_1+g_2.
\tag{W2002 after 16}
$$

The original parameter table is

| parameter | value |
|---|---:|
| \(\alpha\) | \(-0.124\ {\rm GeV}\) |
| \(\beta\) | \(0.071\ {\rm GeV}\) |
| \(\rho_0\) | \(0.165\ {\rm fm^{-3}}\) |
| \(g_0\) | \(0.96\ {\rm GeV\,fm^{-5}}\) |
| \(C_s\) | \(0.032\ {\rm GeV}\) |

### ImQMD05 EDF and explicit momentum-dependent interaction

For intermediate-energy ImQMD05, Zhang 2012 Eq. (2) writes

$$
U=U_\rho+U_{\rm md}+U_{\rm Coul}.
\tag{Z2012-2}
$$

The local term, Zhang 2012 Eq. (4), is

$$
u_\rho=
\frac{\alpha}{2}\frac{\rho^2}{\rho_0}
+\frac{\beta}{\eta+1}\frac{\rho^{\eta+1}}{\rho_0^\eta}
+\frac{g_{\rm sur}}{2\rho_0}(\nabla\rho)^2
+\frac{g_{\rm sur,iso}}{\rho_0}\left[\nabla(\rho_n-\rho_p)\right]^2
+\frac{C_s}{2}\left(\frac{\rho}{\rho_0}\right)^{\gamma_i}\delta^2\rho
+g_{\rho\tau}\frac{\rho^{8/3}}{\rho_0^{5/3}}.
\tag{Z2012-4}
$$

The explicit momentum-dependent energy density is Zhang 2012 Eq. (5):

$$
u_{\rm md}=
\frac{1}{2\rho_0}\sum_{N_1,N_2}\frac{1}{16\pi^6}
\int d^3p_1\,d^3p_2\,
f_{N_1}(\mathbf p_1)f_{N_2}(\mathbf p_2)
\,1.57\left[\ln\left(1+5\times10^{-4}(\Delta p)^2\right)\right]^2,
\tag{Z2012-5}
$$

with \(\Delta p=|\mathbf p_1-\mathbf p_2|\), momenta in MeV/c, and energy in MeV. In this paper the local coefficients are

$$
\alpha=-356\ {\rm MeV},\quad
\beta=303\ {\rm MeV},\quad
\eta=7/6,\quad
g_{\rm sur}=19.47\ {\rm MeV\,fm^2},
$$

$$
g_{\rm sur,iso}=-11.35\ {\rm MeV\,fm^2},\quad
C_s=35.19\ {\rm MeV},\quad
g_{\rho\tau}=0.
$$

The ImQMD05 symmetry energy is Zhang 2012 Eq. (6):

$$
S(\rho)=\frac{1}{3}\frac{\hbar^2}{2m}\rho_0^{2/3}
\left(\frac{3\pi^2}{2}\frac{\rho}{\rho_0}\right)^{2/3}
+\frac{C_s}{2}\left(\frac{\rho}{\rho_0}\right)^{\gamma_i}.
\tag{Z2012-6}
$$

### Coulomb energy

The original Coulomb energy is Wang 2002 Eq. (17):

$$
U_{\rm Coul}^{\rm direct}
=\frac{1}{2}\sum_{i\ne j}
\int \rho_i(\mathbf r)\frac{e^2}{|\mathbf r-\mathbf r'|}
\rho_j(\mathbf r')\,d^3r\,d^3r',
\tag{W2002-17}
$$

where the sum is over charged particles. Later MNT/fusion applications include direct plus exchange Coulomb. Zhao 2016 Eq. (2) gives

$$
U_{\rm Coul}=
\frac{1}{2}\int\!\!\int \rho_p(\mathbf r)
\frac{e^2}{|\mathbf r-\mathbf r'|}
\rho_p(\mathbf r')\,d^3r\,d^3r'
-\frac{3}{4}e^2\left(\frac{3}{\pi}\right)^{1/3}
\int \rho_p(\mathbf R)^{4/3}\,d^3R.
\tag{Zhao2016-2}
$$

The second term is the Slater exchange approximation.

### Density and gradients from wave packets

At any point,

$$
\rho_q(\mathbf r)=\sum_{i\in q}\frac{1}{(2\pi\sigma_{r,i}^2)^{3/2}}
\exp\left[-\frac{(\mathbf r-\mathbf r_i)^2}{2\sigma_{r,i}^2}\right],
\qquad q=n,p.
\tag{density}
$$

For a common width \(\sigma_r\),

$$
\nabla\rho_i(\mathbf r)=
-\frac{\mathbf r-\mathbf r_i}{\sigma_r^2}\rho_i(\mathbf r).
\tag{grad-density}
$$

The implementation may evaluate EDF integrals on a mesh or use analytic Gaussian-folded sums where available. The audit target is equivalence to the EDF above, not the particular quadrature method.

## 3. Ground State Initialization

### Original 2002 initialization

Wang 2002 Sec. 2.4 initializes nuclei as follows:

1. Generate neutron and proton density distributions, originally from RMF calculations.
2. Sample nucleon positions according to those densities.
3. Use a local density approximation to determine the Fermi momentum.
4. Sample centroid momenta from a reduced Fermi sphere because each centroid momentum is itself smeared by \(\sigma_p\), with \(\sigma_r\sigma_p=\hbar/2\).
5. Evolve the prepared nucleus for at least \(600\ {\rm fm}/c\).
6. Retain only nuclei with stable binding energy, rms radius, density, momentum distribution, phase-space distribution, and no spurious particle emission.

The empirical rms-radius comparison used in Wang 2002 Eq. (19) is

$$
\langle r^2\rangle^{1/2}=0.82A^{1/3}+0.58\ {\rm fm}.
\tag{W2002-19}
$$

### 2014 hard-sphere initialization with neutron skin

Wang 2014 Sec. II.B uses a hard-sphere initialization corrected for neutron skin and Gaussian packet width. The charge-radius formula is Eq. (8):

$$
R_c=1.226A^{1/3}+2.86A^{-2/3}
-1.09(I-I^2)+0.99\frac{\Delta E}{A},
\tag{W2014-8}
$$

where \(I=(N-Z)/A\). The neutron skin formula is Eq. (9):

$$
\Delta R_{np}=\langle r_n^2\rangle^{1/2}
-\langle r_p^2\rangle^{1/2}
=0.9I-0.03.
\tag{W2014-9}
$$

The paper converts charge radius to proton hard-sphere radius using

$$
\langle r_c^2\rangle^{1/2}=\sqrt{\frac{3}{5}}R_c,
\qquad
R_p=\sqrt{\frac{5}{3}}\sqrt{\langle r_c^2\rangle-0.64},
\tag{init-rp}
$$

and then

$$
R_n=\sqrt{\frac{5}{3}}\left(\langle r_p^2\rangle^{1/2}+\Delta R_{np}\right).
\tag{init-rn}
$$

Positions are sampled uniformly inside

$$
R_p-w_r \quad\text{for protons},\qquad
R_n-w_r \quad\text{for neutrons},
\tag{init-pos}
$$

with

$$
w_r=0.8\ {\rm fm}.
$$

For each sampled position, momenta are sampled inside a local Fermi sphere:

$$
p_{F,q}(\mathbf r_i)=
\hbar\left[3\pi^2\rho_q(\mathbf r_i)\right]^{1/3}-w_p,
\qquad q=n,p.
\tag{init-pf}
$$

The value of \(w_p\) is tuned so that the sampled nucleus has the experimental ground-state binding energy. A sampled nucleus is accepted only when

$$
E_{\rm ground}\in BE\pm0.05\ {\rm MeV}
\tag{init-energy}
$$

and every pair satisfies the phase-space distance condition

$$
|\mathbf r_i-\mathbf r_j|\,|\mathbf p_i-\mathbf p_j|
\ge 255\ {\rm fm\,MeV}/c.
\tag{init-phase}
$$

The accepted nucleus must remain stable for thousands of fm/c in fusion applications. Wang 2014 reports, for IQ3a, average spurious emissions at \(t=2000\ {\rm fm}/c\) of about 1.1 nucleons for \(^{92}{\rm Zr}\) and 2.6 nucleons for \(^{132}{\rm Sn}\).

### Phase-space constraint during propagation

The original model uses a phase-space-density constraint following CoMD/Papa. If the one-body occupation \(\bar f_i>1\), the code performs many-body elastic scattering to reduce occupation and applies Pauli blocking as in binary collisions. Later papers describe this as a phase-space occupation constraint applied each time step, with an energy check after the corrective elastic scattering.

### Center-of-mass correction

A complete implementation should remove initial center-of-mass drift:

$$
\mathbf R_{\rm cm}=\frac{1}{A}\sum_i\mathbf r_i,\qquad
\mathbf P_{\rm cm}=\sum_i\mathbf p_i,
$$

$$
\mathbf r_i\leftarrow \mathbf r_i-\mathbf R_{\rm cm},\qquad
\mathbf p_i\leftarrow \mathbf p_i-\frac{\mathbf P_{\rm cm}}{A}.
\tag{cm-correction}
$$

For fragment excitation energies, the same rest-frame correction is required before computing internal kinetic energy.

## 4. Equations of Motion

The fundamental equations are Hamilton's equations for centroids:

$$
\dot{\mathbf r}_i=\frac{\partial H}{\partial \mathbf p_i},
\qquad
\dot{\mathbf p}_i=-\frac{\partial H}{\partial \mathbf r_i}.
\tag{eom}
$$

The force is the EDF gradient with respect to the Gaussian centroid:

$$
\mathbf F_i=-\frac{\partial U}{\partial \mathbf r_i}
=-\frac{\partial}{\partial \mathbf r_i}
\left[
\int V_{\rm loc}(\rho_n,\rho_p,\nabla\rho_n,\nabla\rho_p)\,d^3r
+U_{\rm Coul}
\right].
\tag{force}
$$

For explicit momentum-dependent interactions,

$$
\dot{\mathbf r}_i=\frac{\mathbf p_i}{m}
+\frac{\partial U_{\rm md}}{\partial \mathbf p_i},
\qquad
\dot{\mathbf p}_i=-
\frac{\partial (U_{\rho}+U_{\rm Coul}+U_{\rm md})}{\partial \mathbf r_i}.
\tag{mdi-eom}
$$

The literature specifies Hamiltonian propagation and uses event time steps; for example, Yao and Wang 2017 state \(\Delta t=1\ {\rm fm}/c\) for an ImQMD MNT calculation. The supplied priority papers do not explicitly mandate RK4. If a Python implementation uses RK4, the RK4 step must integrate the Hamiltonian right-hand side above and conserve total energy to the accuracy expected by the paper-level stability tests. A valid RK4 step is

$$
y_{n+1}=y_n+\frac{\Delta t}{6}(k_1+2k_2+2k_3+k_4),
\tag{rk4}
$$

with \(y=(\mathbf r_1,\ldots,\mathbf r_A,\mathbf p_1,\ldots,\mathbf p_A)\) and \(k_a=F_{\rm Ham}(y_a)\).

## 5. Nucleon-Nucleon Collisions

### Geometric collision criterion

The QMD-family collision prescription schedules an attempted two-body collision when two centroids approach within the geometric cross-section radius:

$$
d_{\rm min}\le \sqrt{\frac{\sigma_{NN}^{\rm med}(\sqrt{s},\rho)}{\pi}}.
\tag{collision-geom}
$$

The 2020 review writes this criterion in the same form, \(d_{12}\le\sqrt{\sigma_{\rm tot}/\pi}\). In a time-step implementation, compute closest approach over \(t\in[0,\Delta t]\):

$$
\mathbf r_{ij}(t)=\mathbf r_i-\mathbf r_j+(\mathbf v_i-\mathbf v_j)t,
\qquad
t_*=\mathrm{clip}\left[
-\frac{\mathbf r_{ij}(0)\cdot\mathbf v_{ij}}{|\mathbf v_{ij}|^2},
0,\Delta t\right],
\tag{closest}
$$

$$
d_{\rm min}=|\mathbf r_{ij}(t_*)|.
\tag{dmin}
$$

Final relative momentum directions are sampled from the adopted differential NN cross section, while conserving total pair momentum and pair energy. Collision attempts must avoid double-counting the same pair within one collision step.

### Free NN cross sections

Zhang 2012 states that the isospin-dependent free NN cross sections are taken from Cugnon. A commonly used Cugnon-style elastic parameterization in transport codes is

$$
\sigma_{pp}^{\rm free}=\sigma_{nn}^{\rm free}
=13.73-\frac{15.04}{\beta}+\frac{8.76}{\beta^2}
+68.67\beta^4\quad {\rm mb},
\tag{cugnon-pp}
$$

$$
\sigma_{np}^{\rm free}
=-70.67-\frac{18.18}{\beta}+\frac{25.26}{\beta^2}
+113.85\beta\quad {\rm mb},
\tag{cugnon-np}
$$

with

$$
\beta^2=1-\left(\frac{m_N}{m_N+E_{\rm lab}}\right)^2.
\tag{cugnon-beta}
$$

At low energy, later ImQMD discussions often cap free cross sections approximately at

$$
\sigma_{nn/pp}^{\rm free}=60\ {\rm mb},\qquad
\sigma_{np}^{\rm free}=180\ {\rm mb}
\quad (p_{\rm lab}<0.3\ {\rm GeV}/c).
\tag{low-energy-cap}
$$

Because the exact Cugnon implementation is code-level rather than fully restated in Zhang 2012, the audit should check that the selected free-cross-section routine is explicitly documented and consistently converted to fm\(^2\) using

$$
1\ {\rm mb}=0.1\ {\rm fm^2}.
$$

### In-medium scaling

Zhang 2012 uses the phenomenological reduction

$$
\sigma_{nn/np}^*
=\left(1-\xi(E_{\rm beam})\frac{\rho}{\rho_0}\right)
\sigma_{nn/np}^{\rm free},
\tag{Z2012-medium}
$$

with

$$
\xi(E_{\rm beam}=50A{\rm MeV})=0.2.
\tag{Z2012-xi}
$$

The paper compares three cases:

$$
\sigma_{nn/np}^*=\sigma_{nn/np}^{\rm free},
\tag{case1}
$$

$$
\sigma_{nn/np}^*=
\left(1-0.2\frac{\rho}{\rho_0}\right)\sigma_{nn/np}^{\rm free},
\tag{case2}
$$

and an isospin-independent reduced case

$$
\sigma_{nn}^*=\sigma_{pp}^*=\sigma_{np}^*
=\left(1-0.2\frac{\rho}{\rho_0}\right)\sigma',
\tag{case3}
$$

where

$$
\sigma'=
\frac{
2N_{np}\sigma_{np}^{\rm free}
+(N_{nn}+N_{pp})\sigma_{nn/pp}^{\rm free}
}{N_{NN}}.
\tag{sigma-prime}
$$

Chen 2024 uses a more general extraction form:

$$
\sigma_{NN}^{\rm med}
=\left(1+\eta(\sqrt{s})\frac{\rho}{\rho_0}\right)
\sigma_{NN}^{\rm free}.
\tag{Chen2024-5}
$$

This form allows fitted \(\eta\) values that may be negative for suppression.

### Pauli blocking

For an attempted collision \(i,j\), compute final-state occupations \(P_i,P_j\). The collision is blocked with Uehling-Uhlenbeck probability

$$
P_{\rm block}=1-(1-P_i)(1-P_j).
\tag{pauli-block}
$$

The conventional QMD occupation from phase-space density is Chen 2024 Eq. (1):

$$
P_i=P(\mathbf r_i,\mathbf p'_i)
=\frac{1}{4/h^3}\sum_{k\ne i}\frac{1}{(\pi\hbar)^3}
\exp\left[
-\frac{(\mathbf r_i-\mathbf R_k)^2}{2\sigma_r^2}
-\frac{(\mathbf p'_i-\mathbf P_k)^2}{2\sigma_p^2}
\right].
\tag{Chen2024-1}
$$

The hard-sphere overlap alternative is Chen 2024 Eq. (2):

$$
P_i=P(\mathbf r_i,\mathbf p'_i)
=\frac{1}{4/h^3}\sum_{k\ne i}O_{ik}^{(x)}O_{ik}^{(p)},
\tag{Chen2024-2}
$$

where \(O_{ik}^{(x)}\) and \(O_{ik}^{(p)}\) are overlap volumes of coordinate and momentum spheres.

### Chen 2024 improved Pauli blocking PB(W*)

Chen 2024 reduces occupation-probability fluctuations by treating each nucleon as occupying a series of momentum states. The improved occupation is Eq. (6):

$$
P_i=P(\mathbf r_i,\mathbf p'_i)
=\frac{1}{4/h^3}\frac{1}{(\pi\hbar)^3}
\sum_{j\ne i}
\exp\left[-\frac{(\mathbf r_i-\mathbf R_j)^2}{2\sigma_r^2}\right]
c_j(\mathbf p'_i),
\tag{Chen2024-6}
$$

where Eq. (7) defines

$$
c_j(\mathbf p'_i)
=\frac{1}{N}\sum_{\lambda=1}^{N}
\exp\left[-\frac{(\mathbf p'_i-\mathbf P_{j,\lambda})^2}{2\sigma_p^2}\right].
\tag{Chen2024-7}
$$

The conventional method corresponds to Eq. (8):

$$
c_j(\mathbf p'_i)
=\exp\left[-\frac{(\mathbf p'_i-\mathbf P_j)^2}{2\sigma_p^2}\right].
\tag{Chen2024-8}
$$

Chen 2024 reports that PB(W*) reduces the variance of \(P(\mathbf p')\) by about 70% in nuclear matter at \(T=5\) MeV and suppresses spurious first-step surface collisions in \(^{124}{\rm Sn}\) by about 60%. A complete modern implementation should support both PB(W) and PB(W*) as selectable algorithms because older ImQMD papers used conventional blocking.

## 6. Fragment Recognition

### Minimum Spanning Tree

The standard QMD/ImQMD recognition is minimum spanning tree (MST). Zhang 2012 and Zhang-Li-Zhou-Tsang 2012 state that nucleons belong to the same fragment if connected by neighbor links satisfying

$$
|\mathbf r_i-\mathbf r_j|\le R_0,
\qquad
|\mathbf p_i-\mathbf p_j|\le P_0.
\tag{MST}
$$

The relation is transitive: if \(i\) is linked to \(j\), and \(j\) to \(k\), all three are in one cluster even if \(i\) and \(k\) are not directly linked.

Typical ImQMD05 values in Zhang 2012 are

$$
R_0=3.5\ {\rm fm},\qquad
P_0=250\ {\rm MeV}/c.
\tag{Z2012-MST}
$$

The nucleon-induced ImQMD05 paper uses

$$
R_c=4.5\ {\rm fm},\qquad
P_c=250\ {\rm MeV}/c,
\tag{Wei2014-MST}
$$

for its specific spallation/light-cluster study.

### Isospin-dependent MST

Zhang, Li, Zhou, Tsang 2012 introduces iso-MST. The coordinate cutoffs become pair-isospin dependent:

$$
R_{0,nn}=6\ {\rm fm},\qquad
R_{0,np}=6\ {\rm fm},\qquad
R_{0,pp}=3\ {\rm fm},
\tag{iso-MST-R}
$$

while

$$
P_0=250\ {\rm MeV}/c
\tag{iso-MST-P}
$$

is unchanged. The motivation is to account phenomenologically for neutron skins/halos and Coulomb repulsion among protons. The paper reports that iso-MST suppresses \(Z=1\) yields, enhances fragments especially for \(Z\ge 12\), and increases neutron-rich fragment production at mid-rapidity.

### Light complex particle surface coalescence

For nucleon-induced reactions, Wei-Wang-Ou 2014 adds surface coalescence for

$$
d,\quad t,\quad ^3{\rm He},\quad ^4{\rm He}.
$$

The model defines a core radius \(R_0\) and a surface thickness \(D_0\). A fast nucleon leaving the compound system through the surface is the leading nucleon. Candidate cluster nucleons are added using the phase-space condition, Wei 2014 Eq. (12):

$$
R_{im}P_{im}\le h_0,
\qquad
R_{im}\ge 1\ {\rm fm}.
\tag{Wei2014-12}
$$

Here \(R_{im}\) and \(P_{im}\) are Jacobian relative coordinate and momentum of candidate nucleon \(i\) with respect to subgroup \(m\). Candidate emission is checked in priority order

$$
^4{\rm He} > ^3{\rm He} > t > d.
\tag{LCP-priority}
$$

The candidate LCP kinetic energy is Wei 2014 Eq. (13):

$$
E_{\rm lcp}=\sum_{i=1}^{A_{\rm lcp}}(E_i+V_i)+B_{\rm lcp}.
\tag{Wei2014-13}
$$

The cluster is emitted only if it can overcome/tunnel through the Coulomb barrier; otherwise constituents remain in the ImQMD system and the leading nucleon may be emitted as a free nucleon.

Fitted surface coalescence parameters:

$$
R_0\simeq 1.4 A^{1/3}\ {\rm fm},
\qquad
D_0=2.3\ {\rm fm},
\tag{Wei-RD}
$$

and

$$
h_0=
\begin{cases}
200\ {\rm MeV\,fm}/c, & E_{\rm lab}\le 300\ {\rm MeV},\\
260\ {\rm MeV\,fm}/c, & 300<E_{\rm lab}\le 500\ {\rm MeV},\\
330\ {\rm MeV\,fm}/c, & E_{\rm lab}>500\ {\rm MeV}.
\end{cases}
\tag{Wei-h0}
$$

This mechanism is optional for near-barrier fusion/MNT but belongs to the broader complete ImQMD-family behavior for LCP production.

## 7. De-excitation

The dynamical ImQMD calculation produces primary fragments. The fragment rest-frame internal energy is computed from constituent kinetic energies, mutual potential energies, and Coulomb terms with the fragment center-of-mass motion removed. The excitation energy is

$$
E^*(Z,A)=E_{\rm int}(Z,A)-E_{\rm ground}(Z,A).
\tag{excitation}
$$

Zhao 2016 states this as subtracting the corresponding ground-state energy from the total energy of the excited fragment in its rest frame. Yao and Wang 2017 use the same definition for GEMINI coupling.

### HIVAP coupling

For \(^{238}{\rm U}+^{238}{\rm U}\) at \(7.0\ {\rm MeV}/A\), Zhao 2016 terminates ImQMD at

$$
t=1000\ {\rm fm}/c
$$

after reseparation of the composite system, recognizes primary fragments, and passes each \((Z,A,E^*)\) fragment to HIVAP. HIVAP treats

$$
\gamma,\quad n,\quad p,\quad \alpha,\quad {\rm fission}
\tag{HIVAP-channels}
$$

channels. Survival/decay probabilities are calculated from branching ratios

$$
P_i(Z,A,E^*)=
\frac{\Gamma_i(Z,A,E^*)}{\Gamma_{\rm tot}(Z,A,E^*)},
\qquad
\Gamma_{\rm tot}=\sum_i\Gamma_i,
\tag{HIVAP-branch}
$$

with \(i=\gamma,n,p,\alpha,\) fission.

### Statistical evaporation and fission model expectations

A complete de-excitation module should at minimum reproduce the HIVAP/GEMINI coupling pattern:

- Weisskopf-Ewing evaporation for \(n,p,\alpha\) and optionally heavier LCP channels:

$$
\Gamma_j(E^*)\propto
\int_0^{E^*-B_j}
\sigma_j^{\rm inv}(\epsilon)\,
\rho_d(E^*-B_j-\epsilon)\,
\epsilon\,d\epsilon,
\tag{Weisskopf}
$$

where \(B_j\) is separation plus Coulomb barrier energy and \(\rho_d\) is daughter level density.

- Gamma cooling, if included, as a competing width \(\Gamma_\gamma\).

- Bohr-Wheeler fission competition:

$$
\Gamma_f(E^*)\propto
\frac{1}{2\pi\rho_p(E^*)}
\int_0^{E^*-B_f}\rho_s(E^*-B_f-\epsilon)\,d\epsilon.
\tag{Bohr-Wheeler}
$$

The supplied priority ImQMD papers do not rederive the Weisskopf or Bohr-Wheeler formulas; they use established external codes. Therefore, for auditing, an internal de-excitation implementation should be checked against external-code behavior, not merely against the schematic formulas above.

## 8. Reaction Observables

### Impact-parameter integration

For fusion/capture, Wang 2014 Eq. (10) gives

$$
\sigma_{\rm fus}(E_{\rm c.m.})
=2\pi\int b\,g_{\rm fus}(E_{\rm c.m.},b)\,db
\simeq 2\pi\sum_b b\,g_{\rm fus}(E_{\rm c.m.},b)\Delta b.
\tag{W2014-10}
$$

For primary fragments in MNT, Zhao 2016 Eq. (3) gives

$$
\sigma(Z,A,E^*)=
\int_0^{b_{\rm max}}2\pi b\,db\,
\frac{N_{\rm frag}(Z,A,b,E^*)}{N_{\rm tot}(b)}
\simeq
\sum_{b=0}^{b_{\rm max}}2\pi b\Delta b\,
\frac{N_{\rm frag}(Z,A,b,E^*)}{N_{\rm tot}(b)}.
\tag{Zhao2016-3}
$$

For energy spectra and double-differential cross sections in nucleon-induced reactions, Wei 2014 Eqs. (14)-(15) give

$$
\frac{d\sigma}{dE}
=\sum_i^{i_{\rm max}}2\pi b_i\Delta b\,f(E,b_i),
\tag{Wei2014-14}
$$

$$
\frac{d^2\sigma}{d\Omega\,dE}
=\sum_i^{i_{\rm max}}2\pi b_i\Delta b\,f(E,\Omega,b_i).
\tag{Wei2014-15}
$$

Discrete isotope and mass distributions are obtained by summing the primary or residual cross section over unobserved variables:

$$
\frac{d\sigma}{dZ}(Z)=\sum_A\int dE^*\,\sigma(Z,A,E^*),
\qquad
\frac{d\sigma}{dA}(A)=\sum_Z\int dE^*\,\sigma(Z,A,E^*).
\tag{dZdA}
$$

Angular distributions follow from fragment momenta:

$$
\theta_{\rm lab}=\arccos\left(\frac{p_z}{|\mathbf p|}\right),
\qquad
\frac{d\sigma}{d\Omega}\simeq
\sum_b 2\pi b\Delta b\,
\frac{N_{\rm frag}(\theta\in\Delta\theta;b)}
{N_{\rm tot}(b)\Delta\Omega}.
\tag{angular}
$$

### Fusion event setup

Wang 2014 boosts sampled nuclei with

$$
E_{\rm kin}=E_{\rm c.m.}-\frac{Z_1Z_2e^2}{R_0}.
\tag{boost}
$$

The initial separation along the beam direction is

$$
d_0=30\ {\rm fm}
$$

for intermediate fusion systems and

$$
d_0=40\ {\rm fm}
$$

for stronger Coulomb systems such as \(^{132}{\rm Sn}+^{40}{\rm Ca}\). A simulated event is counted as fusion/capture when the center-to-center distance becomes smaller than the compound-nucleus radius. The paper uses about 100-200 events per \((E_{\rm c.m.},b)\).

### Zhao 2016 U+U benchmark expectations

Zhao 2016 applies ImQMD+HIVAP to

$$
^{238}{\rm U}+^{238}{\rm U}
\quad \text{at}\quad
7.0\ {\rm MeV}/A.
$$

Initialization and event settings:

$$
E_{\rm gs}=7.37\ {\rm MeV/nucleon},\qquad
\beta_2=0.215,\qquad
\beta_4=0.093,
\tag{U-init}
$$

random uranium orientations with equal probability,

$$
b_{\rm max}=15\ {\rm fm},\qquad
\Delta b=0.15\ {\rm fm},\qquad
d_0=40\ {\rm fm},
\tag{U-b}
$$

and

$$
N_{\rm events}(b)=100000.
\tag{U-events}
$$

The paper's qualitative benchmark expectations are:

- Primary fragments with \(Z=70\) to 120 are produced by proton and neutron transfer; most primary products lie near the \(^{238}{\rm U}\) isospin asymmetry \(I=0.227\).
- Residual transuranium survival is strongly reduced by fission; the plotted lower bound is \(10^{-8}\) mb.
- About sixty unknown neutron-rich isotopes from Ra (\(Z=88\)) to Db (\(Z=105\)) are predicted above \(10^{-8}\) mb.
- Almost all unknown neutron-rich isotopes are emitted at \(\theta_{\rm lab}\le 60^\circ\).
- Unknown uranium residues with \(A\ge244\) mainly arise from \(b=4\)-8 fm, contact times about \(200\)-\(400\ {\rm fm}/c\), primary \(E^*\le30\) MeV, and outgoing angles \(30^\circ\)-\(60^\circ\).
- Unknown rutherfordium residues with \(Z=104,A\ge269\) require longer contact times, have higher primary excitation energies, much smaller residual cross sections, and narrower outgoing angles around \(40^\circ\)-\(50^\circ\).

## 9. Complete Parameter Table

The following table consolidates the requested IQ parameter sets. Values are in MeV, MeV fm\(^2\), fm\(^2\), fm\(^{-3}\), and fm as indicated. Widths use \(\sigma_r=\sigma_0+\sigma_1A^{1/3}\).

| set | \(\alpha\) (MeV) | \(\beta\) (MeV) | \(\gamma\) | \(g_{\rm sur}\) or \(g_0\) (MeV fm\(^2\)) | \(g_\tau\) (MeV) | \(\eta\) | \(C_s\) (MeV) | \(\kappa_s\) (fm\(^2\)) | \(\rho_0\) (fm\(^{-3}\)) | \(\sigma_0\) (fm) | \(\sigma_1\) (fm) | source |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|
| IQ1 | -310 | 258 | 7/6 | 19.8 | 9.5 | 2/3 | 32 | 0.08 | 0.165 | 0.49 | 0.160 | Li et al. Chinese Phys. C 37, 114101 (2013); width equals Wang 2002 Eq. (18) |
| IQ2 | -356 | 303 | 7/6 | 7.0 | 12.5 | 2/3 | 32 | 0.08 | 0.165 | 0.88 | 0.090 | Li et al. Chinese Phys. C 37, 114101 (2013); Zhao 2016 Table I gives same non-width EDF values |
| IQ3 | -207 | 138 | 7/6 | 18.0 | 14.0 | 5/3 | 32 | 0.08 | 0.165 | 0.94 | 0.018 | Li et al. Chinese Phys. C 37, 114101 (2013); discussed as prior IQ3 in Wang 2014 |
| IQ3a | -207 | 138 | 7/6 | 16.5 | 14.0 | 5/3 | 34 | 0.40 | 0.165 | 0.94 | 0.020 | Wang 2014 Table I |
| IQ3b | -207 | 138 | 7/6 | 18.0 | 14.0 | 5/3 | 34 | 0.60 | 0.165 | 0.94 | 0.018 | Wang 2014 Table I |

Wang 2014 also defines SkP*:

| set | \(\alpha\) | \(\beta\) | \(\gamma\) | \(g_{\rm sur}\) | \(g_\tau\) | \(\eta\) | \(C_s\) | \(\kappa_s\) | \(\rho_0\) | \(\sigma_0\) | \(\sigma_1\) | source |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|
| SkP* | -356 | 303 | 7/6 | 19.5 | 13 | 2/3 | 35 | 0.65 | 0.162 | 0.94 | 0.018 | Wang 2014 Table I |

Original Wang 2002 parameters, in MeV units for comparison:

| parameter | value |
|---|---:|
| \(\alpha\) | -124 MeV |
| \(\beta\) | 71 MeV |
| \(\rho_0\) | 0.165 fm\(^{-3}\) |
| \(g_0\) | 960 MeV fm\(^{-5}\) in the paper table notation |
| \(C_s\) | 32 MeV |
| \(\sigma_r\) | \(0.16A^{1/3}+0.49\) fm |

Note on \(g_0/g_{\rm sur}\): the original 2002 paper's \(g_0\) table uses the dimensions and normalization of its analytic surface-energy expression. Later tables use \(g_{\rm sur}\) in Eq. (W2014-5), with units MeV fm\(^2\). An implementation must not mix these without checking the EDF normalization.

## 10. Expected Complete-Model Behavior Checklist

A complete ImQMD implementation should satisfy the following audit checks:

- Wave packets: Gaussian density and Wigner functions match Eqs. (W2002-1)-(W2002-7), including \(\sigma_r\sigma_p=\hbar/2\).
- EDF: selectable original ImQMD, fusion-oriented IQ EDF, and ImQMD05/Skyrme-style EDF terms are implemented with correct normalizations and units.
- Coulomb: direct finite-density Coulomb is included; Slater exchange is included for modern fusion/MNT calculations.
- Initialization: hard-sphere/neutron-skin sampling, local Fermi momentum, \(w_r=0.8\) fm correction, \(BE\pm0.05\) MeV energy acceptance, and \(255\ {\rm fm\,MeV}/c\) pair phase-space cut are available.
- Stability: isolated nuclei remain bound for at least 600 fm/c in original checks and about 2000 fm/c in modern fusion checks, with minimal spurious emission.
- Propagation: Hamiltonian equations are integrated with energy conservation adequate for the stability tests.
- Collisions: geometric NN collision scheduling, Cugnon/free NN cross sections, in-medium scaling, and Pauli blocking are consistently applied.
- Pauli blocking: conventional PB(W) and improved Chen 2024 PB(W*) are both available for comparison.
- Fragments: MST with \(R_0=3.5\) fm and \(P_0=250\) MeV/c is implemented; iso-MST with \(R_{nn}=R_{np}=6\) fm, \(R_{pp}=3\) fm is available.
- De-excitation: primary fragment \(E^*\) is computed in the rest frame and passed to HIVAP/GEMINI or an equivalent evaporation/fission module.
- Cross sections: impact-parameter integration follows Eqs. (W2014-10), (Zhao2016-3), and, for spectra, (Wei2014-14)-(Wei2014-15).
- U+U benchmark: the implementation can reproduce the Zhao 2016 workflow settings and qualitative angular/excitation/contact-time trends before being trusted for MNT predictions.
