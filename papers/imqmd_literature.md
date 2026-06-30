# ImQMD Model Literature for a From-Scratch Implementation

Scope: papers found from the requested 2002-2024 search phrases:

- `"ImQMD model" Skyrme initialization equations`
- `"improved quantum molecular dynamics" nucleon-nucleon collision Pauli blocking`
- `"ImQMD" fragment recognition coalescence algorithm`
- `"ImQMD" "multi-nucleon transfer"`
- `"ImQMD" fusion reaction heavy-ion`
- `Wang Ning "quantum molecular dynamics" review`
- `"ImQMD" symmetry energy equation of state`

The list is sorted in the order I would read for implementation: model foundations first, then initialization and EDF details, then collisions/Pauli blocking, then cluster recognition, then fusion/MNT applications.

## Highest-Priority Core Set

If time is limited, start with these 12 papers:

1. Wang, Li, Wu, PRC 65 (2002) 064608.
2. Wang et al., PRC 69 (2004) 034608.
3. Zhang et al., Front. Phys. 15 (2020) 54301.
4. Wang, Ou, Zhang, Li, PRC 89 (2014) 064601.
5. Zhang et al., PRC 85 (2012) 024602.
6. Zhang et al., PRC 85 (2012) 051602(R).
7. Zhang et al., PLB 732 (2014) 186.
8. Chen, Zhang, Li, arXiv:2103.13218 (2021).
9. Chen et al., PRC 109 (2024) L021604.
10. Zhao et al., PRC 94 (2016) 024601.
11. Li et al., PLB 776 (2018) 278.
12. Yao and Wang, PRC 95 (2017) 014607.

## Reading Order and Implementation Notes

### 1. General QMD foundation

#### J. Aichelin, "Quantum molecular dynamics--a dynamical microscopic n-body approach to investigate fragment formation and the nuclear equation of state in heavy ion collisions", Physics Reports 202 (1991) 233.

- Link: https://doi.org/10.1016/0370-1573(91)90094-3
- Why it matters: This is the base QMD reference for Gaussian wave packets, centroid equations of motion, stochastic two-body collisions, Pauli blocking, and fragment formation. Read this before ImQMD-specific modifications so the assumptions in ImQMD are clear.
- Implementation use: baseline Gaussian packet representation, Hamiltonian propagation, collision scheduling, and minimum spanning tree style clustering.

#### Yingxun Zhang, Ning Wang, Qingfeng Li, Li Ou, Junlong Tian, Min Liu, Kai Zhao, Xizhen Wu, Zhuxia Li, "Progress of Quantum Molecular Dynamics model and its applications in Heavy Ion Collisions", Frontiers of Physics 15 (2020) 54301.

- Links: https://doi.org/10.1007/s11467-020-0961-9 and https://arxiv.org/abs/2005.12877
- Why it matters: This is the most useful review for locating the ImQMD family in the broader QMD ecosystem. It reviews ImQMD and UrQMD developments, fusion, MNT, fragmentation, collective flow, particle production, EOS constraints, and in-medium NN cross sections.
- Implementation use: map of model variants and references; sanity check for which mechanisms belong in ImQMD05, ImQMD-Sky, fusion-oriented ImQMD, and application-specific versions.

### 2. Original ImQMD model and low-energy fusion implementation

#### Ning Wang, Zhuxia Li, Xizhen Wu, "Improved quantum molecular dynamics model and its applications to fusion reaction near barrier", Physical Review C 65 (2002) 064608.

- Links: https://doi.org/10.1103/PhysRevC.65.064608 and https://arxiv.org/abs/nucl-th/0201079
- Why it matters: This is the original ImQMD paper. It introduces the improved QMD model, demonstrates stable ground-state nuclei from Li to Pb with one parameter set, and applies it to near-barrier fusion of Ca+Zr systems.
- Implementation use: first source for ImQMD wave packets, total Hamiltonian, effective interaction choices, initialization/stability criteria, fusion event definition, and capture/fusion cross section extraction.

#### Ning Wang, Zhuxia Li, Xizhen Wu, Junlong Tian, Yingxun Zhang, Min Liu, "Further development of the improved quantum molecular dynamics model and its application to fusion reactions near the barrier", Physical Review C 69 (2004) 034608.

- Link: https://doi.org/10.1103/PhysRevC.69.034608
- Why it matters: This is the main follow-up for the "ImQMD-II" development. It is essential because later papers cite it for the mature low-energy fusion implementation.
- Implementation use: improved initialization, interaction terms, stability controls, and fusion application details not fully settled in the 2002 paper.

#### M. Papa, T. Maruyama, A. Bonasera, "Constrained molecular dynamics approach to fermionic systems", Physical Review C 64 (2001) 024612.

- Link: https://doi.org/10.1103/PhysRevC.64.024612
- Why it matters: Later ImQMD papers explicitly adopt the Fermi constraint idea from CoMD to enforce fermionic phase-space occupation and stabilize individual nuclei.
- Implementation use: phase-space occupation constraint; handling spurious over-occupation by momentum reshuffling or elastic scattering-like corrections.

#### M. Papa, G. Giuliani, A. Bonasera, "Constrained molecular dynamics II: An N-body approach to nuclear systems", Journal of Computational Physics 208 (2005) 403.

- Link: https://doi.org/10.1016/j.jcp.2005.02.014
- Why it matters: Computational companion for CoMD constraints. Useful for writing robust code rather than only reproducing formulas.
- Implementation use: numerical implementation of fermionic constraints and event stability diagnostics.

### 3. Skyrme energy-density functional and initialization

#### Ning Wang, Li Ou, Yingxun Zhang, Zhuxia Li, "Microscopic dynamics simulations of heavy-ion fusion reactions induced by neutron-rich nuclei", Physical Review C 89 (2014) 064601.

- Links: https://doi.org/10.1103/PhysRevC.89.064601 and https://arxiv.org/abs/1405.5271
- Why it matters: This is one of the best implementation papers for modern fusion-oriented ImQMD. It gives neutron-skin-aware initialization, local Fermi sphere momentum sampling, phase-space rejection criteria, IQ parameter sets, surface/surface-symmetry terms, and stability checks against Skyrme-Hartree-Fock densities.
- Implementation use: hard-sphere proton/neutron position sampling with wave-packet correction, charge-radius and neutron-skin formulas, local density-dependent Fermi momentum sampling, binding-energy acceptance, phase-space distance cut, IQ2/IQ3/IQ3a-style EDF parameters, and fusion-event criteria.

#### D. Vautherin, D. M. Brink, "Hartree-Fock calculations with Skyrme's interaction. I. Spherical nuclei", Physical Review C 5 (1972) 626.

- Link: https://doi.org/10.1103/PhysRevC.5.626
- Why it matters: Foundational Skyrme-Hartree-Fock paper. ImQMD uses Skyrme-like energy-density functionals, usually omitting spin-orbit in transport.
- Implementation use: reference for Skyrme EDF structure and how nuclear ground-state densities are generated for comparison/initialization validation.

#### J. Bartel, P. Quentin, M. Brack, C. Guet, H.-B. Hakansson, "Towards a better parametrisation of Skyrme-like effective forces: A critical study of the SkM force", Nuclear Physics A 386 (1982) 79.

- Link: https://doi.org/10.1016/0375-9474(82)90403-1
- Why it matters: Source of SkM* and related Skyrme-force comparisons used as density benchmarks in ImQMD initialization studies.
- Implementation use: Skyrme-force parameters and density profiles for validating static nuclei.

#### M. Dutra et al., "Skyrme interaction and nuclear matter constraints", Physical Review C 85 (2012) 035201.

- Link: https://doi.org/10.1103/PhysRevC.85.035201
- Why it matters: Catalogs Skyrme parameterizations and nuclear matter constraints. Useful when choosing parameter sets beyond IQ2/IQ3 or implementing ImQMD-Sky.
- Implementation use: EDF parameter database and constraints on incompressibility, symmetry energy, effective masses, and gradients.

#### Y. Zhang, M. B. Tsang, Z. Li, H. Liu, "Constraints on nucleon effective mass splitting with heavy ion collisions", Physics Letters B 732 (2014) 186-190.

- Link: https://doi.org/10.1016/j.physletb.2014.03.030
- Why it matters: Introduces/uses the ImQMD-Sky code where the mean field is taken directly from Skyrme interactions, including momentum-dependent terms and effective mass splitting.
- Implementation use: mapping standard Skyrme parameters to ImQMD-Sky mean field, Lane potential, symmetry energy, and effective mass terms.

### 4. Equations of motion, NN collisions, Pauli blocking, and EOS observables

#### Yingxun Zhang, D. D. S. Coupland, P. Danielewicz, Zhuxia Li, Hang Liu, Fei Lu, W. G. Lynch, M. B. Tsang, "The Influence of in-medium NN cross-sections, symmetry potential and impact parameter on the isospin observables", Physical Review C 85 (2012) 024602.

- Links: https://doi.org/10.1103/PhysRevC.85.024602 and https://arxiv.org/abs/1009.1928
- Why it matters: This is the key ImQMD05 paper for in-medium NN cross sections, symmetry potential sensitivity, impact parameter dependence, and isospin observables.
- Implementation use: in-medium NN cross-section handling, collision-rate constraints, impact-parameter sampling, symmetry-potential parameterization, neutron/proton double ratios, and isospin transport ratios.

#### M. B. Tsang et al., "Constraints on the Density Dependence of the Symmetry Energy", Physical Review Letters 102 (2009) 122701.

- Links: https://doi.org/10.1103/PhysRevLett.102.122701 and https://arxiv.org/abs/0811.3107
- Why it matters: Uses ImQMD calculations to constrain symmetry energy from Sn+Sn isospin diffusion and neutron/proton yield ratios.
- Implementation use: symmetry energy parameter scan, gamma_i parameterization, connection from ImQMD fragments to EOS constraints.

#### Yingxun Zhang, P. Danielewicz, M. Famiano, Zhuxia Li, W. G. Lynch, M. B. Tsang, "The influence of cluster emission and the symmetry energy on neutron-proton spectral double ratios", Physics Letters B 664 (2008) 145-148.

- Link: https://doi.org/10.1016/j.physletb.2008.05.009
- Why it matters: Shows that cluster emission affects neutron/proton spectral observables used for symmetry-energy extraction.
- Implementation use: validation target for fragment emission effects and free nucleon spectra.

#### Xiang Chen, Yingxun Zhang, Zhuxia Li, "Effects of Pauli blocking and in-medium nucleon-nucleon cross sections on the stopping power at low-intermediate energy heavy ion collisions", arXiv:2103.13218 (2021).

- Link: https://arxiv.org/abs/2103.13218
- Why it matters: Compares three Pauli-blocking algorithms in QMD-type models and quantifies underblocking in nuclear matter and spurious collisions at nuclear surfaces.
- Implementation use: benchmark suite for Pauli blocking: nuclear matter blocking ratio, finite-nucleus spurious collision radial profile, stopping-power sensitivity.

#### Xiang Chen, Junping Yang, Ying Cui, Kai Zhao, Zhuxia Li, Yingxun Zhang, "Novel Pauli blocking method in quantum molecular dynamics type models", Physical Review C 109 (2024) L021604.

- Links: https://doi.org/10.1103/PhysRevC.109.L021604 and https://arxiv.org/abs/2403.00343
- Why it matters: Proposes a smoother occupation-probability calculation for QMD Pauli blocking and evaluates it in ImQMD-like simulations.
- Implementation use: candidate modern Pauli-blocking upgrade, finite-nucleus stability improvement, and revised extraction of in-medium NN cross sections.

#### J. Cugnon, D. L'Hote, J. Vandermeulen, "Simple parametrization of cross-sections for nuclear transport studies up to the GeV range", Nuclear Instruments and Methods in Physics Research B 111 (1996) 215.

- Link: https://doi.org/10.1016/0168-583X(95)01384-9
- Why it matters: Frequently cited for free NN elastic cross-section parameterizations used inside transport models.
- Implementation use: free pp, nn, np cross sections before in-medium scaling.

### 5. Fragment recognition and coalescence

#### Yingxun Zhang, Zhuxia Li, Chengshuang Zhou, M. B. Tsang, "Effect of isospin-dependent cluster recognition on the observables in heavy ion collisions", Physical Review C 85 (2012) 051602(R).

- Link: https://doi.org/10.1103/PhysRevC.85.051602
- Why it matters: Key ImQMD paper for cluster recognition. It introduces isospin dependence into fragment recognition and shows that fragment yields and isospin observables change significantly.
- Implementation use: cluster-recognition criteria, isospin-dependent coalescence/MST tuning, fragment yields, isotope distributions, isoscaling, and equilibration observables.

#### Dexian Wei, Ning Wang, Li Ou, "Mechanism of production of light complex particles in nucleon-induced reactions", Journal of Physics G: Nuclear and Particle Physics 41 (2014) 035104.

- Links: https://doi.org/10.1088/0954-3899/41/3/035104 and https://arxiv.org/abs/1309.7534
- Why it matters: Adds a phenomenological surface coalescence and emission mechanism to ImQMD for light complex particles d, t, 3He, and 4He.
- Implementation use: optional light-cluster afterburner/coalescence module, especially for nucleon-induced or spallation-like reactions.

#### Li Mao, Ning Wang, Li Ou, "Dynamical and statistical description of multifragmentation in heavy-ion collisions", Physical Review C 91 (2015) 044604.

- Link: https://doi.org/10.1103/PhysRevC.91.044604
- Why it matters: Uses ImQMD with statistical decay to describe multifragmentation, and is useful for recognizing the handoff point from dynamical fragments to de-excitation.
- Implementation use: fragment excitation, primary-to-final fragment treatment, GEMINI/statistical de-excitation coupling.

### 6. Fusion applications and capture/fusion validation

#### Y. Y. Jiang, Ning Wang, Zhuxia Li, W. Scheid, "Dynamical nucleus-nucleus potential at short distances", Physical Review C 81 (2010) 044602.

- Link: https://doi.org/10.1103/PhysRevC.81.044602
- Why it matters: Extracts dynamical nucleus-nucleus potentials from ImQMD simulations.
- Implementation use: potential extraction from event trajectories, barrier evolution, and comparison to ETF/frozen-density barriers.

#### Ning Wang, Kai Zhao, Zhuxia Li, "Systematic study of 16O-induced fusion with the improved quantum molecular dynamics model", Physical Review C 90 (2014) 054610.

- Link: https://doi.org/10.1103/PhysRevC.90.054610
- Why it matters: Systematic fusion validation for 16O-induced reactions.
- Implementation use: capture/fusion cross-section workflow, impact parameter integration, barrier/fusion systematics.

#### Kai Wen, Fumihiko Sakata, Zhuxia Li, Xizhen Wu, Yingxun Zhang, Shan-Gui Zhou, "Non-Gaussian fluctuation-dissipation dynamics in heavy-ion fusion", Physical Review Letters 111 (2013) 012501.

- Link: https://doi.org/10.1103/PhysRevLett.111.012501
- Why it matters: Uses ImQMD trajectories to study dissipation and random-force distributions during fusion.
- Implementation use: extracting collective coordinates, friction, random force, and memory effects from microscopic trajectories.

#### Kai Wen, Fumihiko Sakata, Zhuxia Li, Xizhen Wu, Yingxun Zhang, Shan-Gui Zhou, "Energy dependence of the nucleus-nucleus potential and the friction parameter in fusion reactions", Physical Review C 90 (2014) 054613.

- Link: https://doi.org/10.1103/PhysRevC.90.054613
- Why it matters: Extends ImQMD-based extraction of energy-dependent potential and friction parameters.
- Implementation use: validation and analysis layer for dissipative fusion dynamics.

### 7. MNT, heavy fragments, and HIVAP/GEMINI coupling

#### Kai Zhao, Xizhen Wu, Zhuxia Li, "Quantum molecular dynamics study of the mass distribution of products in 7.0 A MeV 238U+238U collisions", Physical Review C 80 (2009) 054607.

- Link: https://doi.org/10.1103/PhysRevC.80.054607
- Why it matters: Early ImQMD application to very heavy U+U collisions, before the later production-mechanism papers.
- Implementation use: event classification, mass distributions, long contact-time dynamics in very heavy systems.

#### Kai Zhao, Zhuxia Li, Ning Wang, Yingxun Zhang, Qingfeng Li, Yongjia Wang, Xizhen Wu, "Production mechanism of neutron-rich transuranium nuclei in 238U+238U collisions at near-barrier energies", Physical Review C 92 (2015) 024613.

- Link: https://doi.org/10.1103/PhysRevC.92.024613
- Why it matters: Direct precursor to the 2016 paper. It connects primary fragments from ImQMD to evaporation residues for transuranium production.
- Implementation use: production mechanism analysis, excitation energy and angle cuts, de-excitation coupling, near-barrier U+U event selection.

#### Kai Zhao, Zhuxia Li, Yingxun Zhang, Ning Wang, Qingfeng Li, Caiwan Shen, Yongjia Wang, Xizhen Wu, "Production of unknown neutron-rich isotopes in 238U+238U collisions at near-barrier energy", Physical Review C 94 (2016) 024601.

- Links: https://doi.org/10.1103/PhysRevC.94.024601 and https://arxiv.org/abs/1605.07393
- Why it matters: You already have this paper, but it remains the central MNT+HIVAP application paper for your current code context.
- Implementation use: ImQMD primary fragment production, HIVAP de-excitation, residual fragment cross sections, emission angles, contact time, excitation energy, and transfer-product interpretation.

#### Cheng Li, Peng Wen, Jiajun Li, Guoqiang Zhang, Bao-An Li, Xinxin Xu, Zhongzhou Ren, "Production mechanism of new neutron-rich heavy nuclei in the 136Xe+198Pt reaction", Physics Letters B 776 (2018) 278-283.

- Link: https://doi.org/10.1016/j.physletb.2017.11.060
- Why it matters: Clear MNT application with model equations summarized in the text. It uses ImQMD for primary fragments and GEMINI for de-excitation, and explicitly lists the IQ2 parameter set.
- Implementation use: density function, Hamilton equations, Coulomb direct/exchange term, local Skyrme EDF, IQ2 parameters, TKE-mass distributions, angular cuts, and GEMINI coupling.

#### Hong Yao, Ning Wang, "Microscopic dynamics simulations of multinucleon transfer in 86Kr+64Ni at 25 MeV/nucleon", Physical Review C 95 (2017) 014607.

- Link: https://doi.org/10.1103/PhysRevC.95.014607
- Why it matters: Focused MNT dynamics paper for a measured system at higher energy than near-barrier actinide cases.
- Implementation use: MNT observables, isotope distributions, impact-parameter dependence, and validation against experimental data.

#### Xiang Jiang, Nan Wang, "Probing the production mechanism of neutron-rich nuclei in multinucleon transfer reactions", Physical Review C 101 (2020) 014604.

- Link: https://doi.org/10.1103/PhysRevC.101.014604
- Why it matters: Later MNT mechanism study comparing production routes for neutron-rich nuclei.
- Implementation use: mechanism-level diagnostics after event generation: transfer paths, excitation energy, survival after de-excitation.

#### Walter Loveland, "The Synthesis of New Neutron-Rich Heavy Nuclei", Frontiers in Physics 7 (2019) 23.

- Link: https://doi.org/10.3389/fphy.2019.00023
- Why it matters: Experimental and model-comparison review for neutron-rich heavy nuclei production. It notes ImQMD predictions can be closer to data but may underestimate yields by factors of 10-100, possibly because shell effects are absent.
- Implementation use: external validation expectations and known limitations of ImQMD for shell-sensitive MNT yields.

#### W. Loveland et al. / V. V. Desai, W. Loveland, K. McCaleb, R. Yanez et al., "The 136Xe+198Pt reaction: A test of models of multi-nucleon transfer reactions", Physical Review C 99 (2019) 044604.

- Link: https://doi.org/10.1103/PhysRevC.99.044604
- Why it matters: Important benchmark for MNT models including ImQMD, DNS, and GRAZING-type approaches.
- Implementation use: validation dataset for target-like/projectile-like fragments and model discrepancy checks.

#### L. Zhu, J. Su, W. J. Xie, F. S. Zhang, "Production of neutron-rich nuclei around N = 126 in multinucleon transfer reactions", Physics Letters B 767 (2017) 437-441.

- Link: https://doi.org/10.1016/j.physletb.2017.01.082
- Why it matters: Not an ImQMD foundation paper, but useful as a competing MNT transport-model application around N=126.
- Implementation use: comparison target for reaction selection and production cross-section expectations.

### 8. Recent and optional extensions through 2024

#### Xiang Chen, Guo-Fang Dai, "Influence of the treatment of initialization and mean-field potential on the neutron-to-proton yield ratios", Physical Review C 104 (2021) 024605.

- Link: https://doi.org/10.1103/PhysRevC.104.024605
- Why it matters: Directly relevant to initialization sensitivity in ImQMD-like calculations.
- Implementation use: uncertainty budget for initialization choices, mean-field potential implementation, and neutron/proton observables.

#### C. Y. Tsang, M. Kurata-Nishimura, M. B. Tsang, W. G. Lynch, Yingxun Zhang et al., "Constraining nucleon effective masses with flow and stopping observables from the SRIT experiment", Physics Letters B 853 (2024) 138661.

- Link: https://doi.org/10.1016/j.physletb.2024.138661
- Why it matters: Recent ImQMD-Sky application constraining effective masses and in-medium cross-section parameters using flow and stopping data.
- Implementation use: modern parameter-range guidance for ImQMD-Sky and validation observables beyond MNT/fusion.

#### Hong Yao, Hui Yang, Ning Wang, "Systematic study of capture thresholds with time dependent Hartree-Fock theory", Physical Review C 110 (2024) 014602.

- Link: https://doi.org/10.1103/PhysRevC.110.014602
- Why it matters: Not ImQMD, but useful as a 2024 microscopic comparison for capture thresholds and fusion barriers.
- Implementation use: benchmark for capture threshold predictions when validating ImQMD fusion barriers.

## Implementation Dependency Map

- Wave packets and propagation: Aichelin 1991; Wang, Li, Wu 2002; Wang et al. 2004.
- Ground-state initialization: Wang et al. 2002; Wang et al. 2004; Wang, Ou, Zhang, Li 2014; Chen and Dai 2021.
- Skyrme EDF and parameter mapping: Vautherin and Brink 1972; Bartel et al. 1982; Dutra et al. 2012; Zhang et al. 2014; Wang, Ou, Zhang, Li 2014.
- NN collisions and Pauli blocking: Zhang et al. 2012 PRC 85 024602; Chen, Zhang, Li 2021; Chen et al. 2024; Cugnon et al. 1996.
- Fermionic phase-space constraint: Papa, Maruyama, Bonasera 2001; Papa, Giuliani, Bonasera 2005.
- Fragment recognition: Zhang, Li, Zhou, Tsang 2012; Wei, Wang, Ou 2014; Mao, Wang, Ou 2015.
- Symmetry energy/EOS validation: Tsang et al. 2009; Zhang et al. 2008; Zhang et al. 2012; Zhang et al. 2014; Tsang et al. 2024.
- Fusion validation: Wang et al. 2002; Wang et al. 2004; Jiang et al. 2010; Wang, Ou, Zhang, Li 2014; Wang, Zhao, Li 2014; Wen et al. 2013/2014.
- MNT/HIVAP/GEMINI applications: Zhao, Wu, Li 2009; Zhao et al. 2015; Zhao et al. 2016; Li et al. 2018; Yao and Wang 2017; Jiang and Wang 2020; Loveland 2019; Desai/Loveland et al. 2019.

## Chinese Summary

共整理出约 27 篇核心/关键论文，其中最核心的实现论文约 12 篇。建议阅读顺序如下：

1. 先读 Wang, Li, Wu 2002 和 Wang et al. 2004：这是 ImQMD 模型本体，包含高斯波包、哈密顿量、初始化、稳定性和近垒融合计算。
2. 再读 2020 年 Zhang/Wang/Li 等综述：快速建立 ImQMD、ImQMD05、ImQMD-Sky、UrQMD 以及融合/MNT/碎裂应用的全局图景。
3. 接着读 Wang, Ou, Zhang, Li 2014：这是实现初始化最重要的论文，包含中子皮、局域费米球采样、相空间约束、Skyrme-Hartree-Fock 密度对比和 IQ 参数组。
4. 碰撞项和 Pauli blocking 读 Zhang et al. 2012、Chen et al. 2021、Chen et al. 2024：分别覆盖介质中 NN 截面、传统 Pauli blocking 的问题和新的改进算法。
5. 碎片识别读 Zhang, Li, Zhou, Tsang 2012：这是 ImQMD 中同位旋相关 cluster recognition 的关键论文；轻复合粒子再看 Wei, Wang, Ou 2014。
6. MNT 应用从 Zhao et al. 2016 开始，然后读 Li et al. 2018 和 Yao/Wang 2017：它们给出 ImQMD + HIVAP/GEMINI 的实际产额计算、角分布、激发能和反应机制分析。

如果要从零实现，第一阶段应先复现静态核稳定性和能量守恒；第二阶段加入 NN collision、Pauli blocking 和 Fermi constraint；第三阶段加入 fragment recognition 与 de-excitation 接口；最后再做 fusion/MNT 反应截面验证。
