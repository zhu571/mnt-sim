# imQMD计算多核子转移反应的代码实现方法研究报告

## 摘要

本报告基于对十余篇使用改进量子分子动力学模型（imQMD）计算多核子转移反应的学术论文的系统梳理，总结了imQMD模型的代码实现方法和计算过程。imQMD模型由王宁、李祝霞、吴锡真等中国核物理学者提出，在标准QMD模型基础上引入了表面能量项、系统尺寸依赖的波包宽度以及相空间占据约束三大改进。计算多核子转移反应的完整流程包括基态核初始化、运动方程演化、二体碰撞处理、碎片识别和统计衰变五个核心模块，通常采用imQMD+GEMINI的两阶段方法实现从动力学到统计衰变的全过程描述。

---

## 一、研究背景

多核子转移反应（Multinucleon Transfer Reaction, MNT）是近库仑势垒能量下重离子碰撞中产生丰中子核素的重要途径，尤其对于合成超重元素和探索N=126中子闭壳附近的丰中子核具有关键意义。量子分子动力学（QMD）模型作为一种非平衡的微观输运模型，能够自洽描述重离子碰撞中从准弹性到深度非弹性的多核子转移过程。然而，标准QMD模型存在核表面过于弥散、动量分布偏离费米分布、波包宽度参数不统一等问题。imQMD模型针对这些问题进行了系统性改进，成为中国学者在重离子反应动力学领域的重要贡献。

---

## 二、论文清单

以下整理了10篇使用imQMD计算多核子转移反应的核心论文：

### 论文1：imQMD模型奠基论文
- 标题：An Improved Quantum Molecular Dynamics Model and its Applications to Fusion Reaction near Barrier
- 作者：Ning Wang, Zhuxia Li, Xizhen Wu
- 期刊：Physical Review C 65, 064608 (2002)
- 链接：https://arxiv.org/abs/nucl-th/0201079
- 反应体系：⁴⁰Ca+⁹⁰Zr, ⁴⁰Ca+⁹⁶Zr, ⁴⁸Ca+⁹⁰Zr等
- 贡献：提出imQMD模型的三个核心改进

### 论文2：imQMD-II模型进一步发展
- 标题：Further Development of the Improved QMD Model and its Applications to Fusion Reaction near Barrier
- 作者：Ning Wang, Zhuxia Li, Xizhen Wu, Junlong Tian, Yingxun Zhang, Min Liu
- 期刊：Physical Review C 69, 034608 (2004)
- 链接：https://arxiv.org/abs/nucl-th/0402066
- 贡献：引入基于SkM*和SLy系列Skyrme相互作用的新参数集IQ1

### 论文3：¹³⁶Xe+²⁰⁸Pb多核子转移
- 标题：Multinucleon transfer in the ¹³⁶Xe + ²⁰⁸Pb reaction
- 作者：Cheng Li, Fan Zhang, Jingjing Li, Long Zhu, Junlong Tian, Ning Wang, Feng-Shou Zhang
- 期刊：Physical Review C 93, 014618 (2016)
- 链接：https://doi.org/10.1103/PhysRevC.93.014618
- 反应体系：¹³⁶Xe+²⁰⁸Pb, E_c.m.=450 MeV

### 论文4：⁸⁶Kr+⁶⁴Ni多核子转移微观动力学模拟
- 标题：Microscopic dynamics simulations of multinucleon transfer in ⁸⁶Kr+⁶⁴Ni at 25 MeV/nucleon
- 作者：Hong Yao, Ning Wang
- 期刊：Physical Review C 95, 014607 (2017)
- 链接：https://doi.org/10.1103/PhysRevC.95.014607
- 反应体系：⁸⁶Kr+⁶⁴Ni, 25 MeV/nucleon

### 论文5：¹²⁴Xe+²⁰⁸Pb同位旋平衡
- 标题：Isospin equilibration in multinucleon transfer reaction at near-barrier energy
- 作者：Cheng Li等
- 期刊：Physical Review C 99, 034619 (2019)
- 链接：https://doi.org/10.1103/PhysRevC.99.034619

### 论文6：¹³⁶Xe+¹⁹⁸Pt丰中子核产生机制
- 标题：Production mechanism of the neutron-rich nuclei in multinucleon transfer reactions
- 作者：Cheng Li等
- 期刊：Physics Letters B 808, 135697 (2020)
- 链接：https://doi.org/10.1016/j.physletb.2020.135697

### 论文7：N=126丰中子核产生
- 标题：Production of neutron-rich N=126 nuclei in multinucleon transfer reactions: Comparison between ¹³⁶Xe+¹⁹⁸Pt and ²³⁸U+¹⁹⁸Pt reactions
- 作者：K. Zhao, Z. Liu, F.S. Zhang, N. Wang等
- 期刊：Physics Letters B 820, 136580 (2021)
- 链接：https://www.sciencedirect.com/science/article/pii/S0370269321000411

### 论文8：Q值效应对¹³⁶Xe+²⁰⁸Pb多核子转移的影响
- 标题：Q值效应对多核子转移反应¹³⁶Xe+²⁰⁸Pb的影响
- 作者：赵凯、夏政通、段济正
- 期刊：原子核物理评论 37, 2 (2020)
- 链接：http://npr.xml-journal.net/article/doi/10.11804/NuclPhysRev.37.2020022
- 反应体系：¹³⁶Xe+²⁰⁸Pb, E_c.m.=617 MeV

### 论文9：¹⁹⁷Au+¹⁹⁷Au三体破裂
- 标题：Improved Quantum Molecular Dynamics Model and Its Application to Ternary Breakup Reactions
- 作者：Junlong Tian, Xian Li, Cheng Li
- 期刊：Universe 8(11), 555 (2022)
- 链接：https://www.mdpi.com/2218-1997/8/11/555
- 反应体系：¹⁹⁷Au+¹⁹⁷Au, 5-30 A MeV

### 论文10：核表面动力学对融合和MNT的影响
- 标题：Effects of nuclear surface dynamics on fusion and multinucleon transfer reactions
- 作者：Cheng Li, Xingxin Luo, Tao Li, Xin-Rui Zhang, Junlong Tian, Ning Wang, Hui-Xiao Duan, Feng-Shou Zhang
- 期刊：Physical Review C 112, 034601 (2025)
- 链接：https://doi.org/10.1103/PhysRevC.112.034601

### 补充论文：imQMD参数集确定
- 标题：Determination of the nucleon-nucleon interaction in the ImQMD model by nuclear reactions at the Fermi energy region
- 作者：Cheng Li, Junlong Tian, Yu-Jiao Qin, Jing-Jing Li, Ning Wang
- 期刊：Chinese Physics C 37, 114101 (2013)
- 链接：https://iopscience.iop.org/article/10.1088/1674-1137/37/11/114101

### 补充综述：imQMD在多核子转移反应中的应用
- 标题：改进的量子分子动力学模型在多核子转移反应中的应用
- 作者：李程、张新瑞、张玉海、张丰收
- 期刊：中国科学：物理学 力学 天文学 55, 122007 (2025)
- 链接：https://doi.org/10.1360/SSPMA-2024-0543

### 补充综述：QMD模型进展
- 标题：Progress of Quantum Molecular Dynamics model and its applications in Heavy Ion Collisions
- 作者：Yingxun Zhang, Ning Wang等
- 期刊：Frontiers of Physics 15, 54301 (2020)
- 链接：https://arxiv.org/abs/2005.12877

---

## 三、imQMD模型的理论基础

### 3.1 高斯波包表示

在imQMD模型中，每个核子用一个高斯相干态波包表示：

$$\phi_i(\vec{r}) = \frac{1}{(2\pi\sigma_r^2)^{3/4}} \exp\left[-\frac{(\vec{r}-\vec{r}_i)^2}{4\sigma_r^2} + \frac{i\vec{p}_i \cdot \vec{r}}{\hbar}\right]$$

其中 $\vec{r}_i$ 和 $\vec{p}_i$ 分别为第 $i$ 个波包在坐标空间和动量空间的中心，$\sigma_r$ 为波包在坐标空间中的展宽。整个 $N$ 体波函数为所有核子波包的直积。通过Wigner变换，单粒子相空间分布函数为：

$$f_i(\vec{r},\vec{p}) = \frac{1}{(\pi\hbar)^3} \exp\left[-\frac{(\vec{r}-\vec{r}_i)^2}{2\sigma_r^2} - \frac{2\sigma_r^2(\vec{p}-\vec{p}_i)^2}{\hbar^2}\right]$$

核子密度分布和动量分布分别为所有波包贡献的求和：

$$\rho(\vec{r}) = \sum_i \frac{1}{(2\pi\sigma_r^2)^{3/2}} \exp\left[-\frac{(\vec{r}-\vec{r}_i)^2}{2\sigma_r^2}\right]$$

坐标空间和动量空间宽度满足最小不确定关系 $\sigma_r \cdot \sigma_p = \hbar/2$。

### 3.2 imQMD相比标准QMD的三大改进

imQMD模型相对于标准QMD模型有三项关键改进，这些改进使其能够更好地描述核基态性质和重离子反应动力学：

第一，在势能密度泛函中引入表面能量项和表面对称能量项。标准QMD模型由于高斯波包的长尾部，导致核表面过于弥散、中心密度不合理升高，产生虚假核子发射。imQMD通过引入密度梯度平方项 $g_0(\nabla\rho)^2/(2\rho_0)$ 来描述表面能量，使中心区域粒子受到排斥（抑制过高中心密度），表面区域粒子受到吸引（防止表面过度延展）。

第二，引入系统尺寸依赖的波包宽度。标准QMD模型对不同质量数的核使用不同的固定波包宽度参数（如Ca+Ca用 $L=4.33$ fm²，Au+Au用 $L=8.66$ fm²），缺乏统一性。imQMD采用公式 $\sigma_r = \sigma_0 + \sigma_1 A^{1/3}$（单位fm），使得一套参数可以描述从轻核到重核的基态性质。

第三，采用相空间占据约束（Fermi约束）模拟费米子特性。标准QMD模型中动量分布随时间演化趋向经典Boltzmann分布，低动量粒子数量不合理增加。imQMD遵循Papa等人CoMD模型的方法，要求单体相空间占据 $\bar{f}_i \leq 1$，当 $\bar{f}_i > 1$ 时执行多体弹性散射以降低相空间占据。

### 3.3 哈密顿量与势能

系统总哈密顿量为 $H = T + U_{\text{loc}} + U_{\text{Coul}}$，其中动能 $T = \sum_i \vec{p}_i^2/(2m)$。

imQMD的核心改进在于采用Skyrme型能量密度泛函，局部相互作用势能密度为：

$$V_{\text{loc}} = \frac{\alpha}{2}\frac{\rho^2}{\rho_0} + \frac{\beta}{\gamma+1}\frac{\rho^{\gamma+1}}{\rho_0^\gamma} + \frac{g_0}{2\rho_0}(\nabla\rho)^2 + g_\tau\frac{\rho^{\eta+1}}{\rho_0^\eta} + \frac{C_s}{2\rho_0}\left(\rho^2 - \kappa_s(\nabla\rho)^2\right)\delta^2$$

各项物理含义分别为：第一项为体积结合能（来自Skyrme $t_0$），第二项为密度依赖体积能（来自Skyrme $t_3$），第三项为表面能项（imQMD关键改进），第四项为动量依赖有效质量相关项，第五项为对称能项（含体对称能和表面对称能，imQMD关键改进）。其中 $\delta = (\rho_n - \rho_p)/(\rho_n + \rho_p)$ 为同位旋不对称度。

对 $V_{\text{loc}}$ 积分后得到局部相互作用势能：

$$U_{\text{loc}} = \frac{\alpha}{2}\sum_{i,j\neq i}\frac{\rho_{ij}}{\rho_0} + \frac{\beta}{\gamma+1}\sum_i\left(\sum_{j\neq i}\frac{\rho_{ij}}{\rho_0}\right)^\gamma + \frac{g_0}{2}\sum_{i,j\neq i} f_s^{(ij)}\frac{\rho_{ij}}{\rho_0} + g_\tau\sum_i\left(\sum_{j\neq i}\frac{\rho_{ij}}{\rho_0}\right)^\eta + \frac{C_s}{2}\sum_{i,j\neq i} t_{iz}t_{jz}\frac{\rho_{ij}}{\rho_0}(1-\kappa_s f_s^{(ij)})$$

其中核子间重叠函数和表面因子分别为：

$$\rho_{ij} = \frac{1}{(4\pi\sigma_r^2)^{3/2}} \exp\left[-\frac{(\vec{r}_i-\vec{r}_j)^2}{4\sigma_r^2}\right], \quad f_s^{(ij)} = \frac{3}{2\sigma_r^2} - \frac{(\vec{r}_i-\vec{r}_j)^2}{4\sigma_r^4}$$

$t_{iz} = +1$（质子）或 $-1$（中子）。库仑能包含直接项和Slater交换近似项：

$$U_{\text{Coul}} = \frac{1}{2}\int\frac{\rho_P(\vec{r})e^2}{|\vec{r}-\vec{r}'|}\rho_P(\vec{r}')d\vec{r}d\vec{r}' - \frac{e^2}{4}\left(\frac{3}{\pi}\right)^{1/3}\int\rho_P^{4/3}dR$$

由于高斯密度分布的解析形式，所有积分均可解析完成，这是imQMD计算效率的重要保证。

---

## 四、代码实现方法与计算过程

### 4.1 参数集

imQMD模型发展了三套参数集IQ1、IQ2、IQ3，其完整对比见下表：

| 参数 | IQ1 | IQ2 | IQ3 | 物理含义 |
|------|-----|-----|-----|----------|
| α (MeV) | -310 | -356 | -207 | 体积结合系数 |
| β (MeV) | 258 | 303 | 138 | 高阶密度系数 |
| γ | 7/6 | 7/6 | 7/6 | 密度幂指数 |
| g₀ (MeV·fm²) | 19.8 | 7.0 | 18.0 | 表面能系数 |
| gτ (MeV) | 9.5 | 12.5 | 14.0 | 动量依赖项系数 |
| η | 2/3 | 2/3 | 5/3 | 动量依赖幂指数 |
| Cs (MeV) | 32.0 | 32.0 | 32.0 | 对称能系数 |
| κs (fm²) | 0.08 | 0.08 | 0.08 | 表面对称能参数 |
| ρ₀ (fm⁻³) | 0.165 | 0.165 | 0.165 | 饱和核密度 |
| σ₀ (fm) | 0.49 | 0.88 | 0.94 | 波包宽度参数1 |
| σ₁ (fm) | 0.16 | 0.09 | 0.018 | 波包宽度参数2 |
| K∞ (MeV) | 165 | 195 | 226 | 不可压缩系数 |

IQ1参数集基于SkM*和SLy系列Skyrme力标定，对应最软的状态方程（K∞=165 MeV）。IQ3参数集的不可压缩系数K∞=226 MeV最接近经验值约230 MeV，被证明是最适合描述重离子碰撞的参数集。在实际多核子转移反应计算中，⁸⁶Kr+⁶⁴Ni使用了IQ3a参数集（IQ3的变种），¹⁹⁷Au+¹⁹⁷Au使用了IQ2参数集。

### 4.2 基态核初始化

基态核的制备是imQMD计算的关键起点，直接影响后续反应模拟的可靠性。初始化过程包括以下步骤：

核子位置根据从相对论平均场（RMF）理论计算获取的中子和质子密度分布进行抽样。核子动量基于局部密度近似计算费米动量 $P_F$，并考虑高斯波包动量宽度修正 $P_F' = P_F - \Delta P_F$。

制备后的核系统需要进行严格的稳定性检验。预制备核系统需演化至少600 fm/c（重核聚变研究要求约3000 fm/c，部分研究要求6000 fm/c），检查均方根半径、结合能、密度分布、动量分布和相空间分布是否保持稳定，且无虚假核子发射。以¹⁹⁷Au为例，其精选标准为结合能 $7.92 \pm 0.05$ MeV/核子，均方根电荷半径 $5.44 \pm 0.2$ fm。

通过稳定性检验后，从数千个采样的核中精选20个弹核和20个靶核作为初始构型。每个初始核围绕质心以随机选择的欧拉角旋转，以消除空间取向的影响。弹核与靶核的初始间距设为50 fm（部分早期工作使用20 fm），确保初始时两核间无相互作用。碰撞参数 $b$ 从0到 $b_{\text{max}}$（典型12-14 fm）取值，步长 $\Delta b = 1.0$ fm。对于每个碰撞参数，通常模拟超过100,000个碰撞事件以保证统计精度。

### 4.3 运动方程演化

核子在自洽生成的平均场中的传播由哈密顿运动方程控制：

$$\dot{\vec{r}}_i = \frac{\partial H}{\partial \vec{p}_i}, \quad \dot{\vec{p}}_i = -\frac{\partial H}{\partial \vec{r}_i}$$

数值积分采用Euler方法，时间步长 $\Delta t = 1$ fm/c。在⁸⁶Kr+⁶⁴Ni的计算中，总模拟时间为1000 fm/c；在¹⁹⁷Au+¹⁹⁷Au的计算中，每100 fm/c记录一次粒子位置和动量用于分析直接过程和级联过程。

在每一步演化中，除了哈密顿方程传播核子位置和动量外，还需施加相空间占据约束。当单体相空间占据 $\bar{f}_i > 1$ 时，执行多体弹性散射以降低相空间占据，同时检查泡利阻塞概率。这一步骤有效维持了费米子特性，防止动量分布向经典Boltzmann分布退化。

### 4.4 二体碰撞处理

在平均场演化之外，imQMD还处理核子-核子二体碰撞。碰撞处理中核子位置保持不变，仅动量发生改变（$\vec{p}_i + \vec{p}_j \to \vec{p}_i' + \vec{p}_j'$），碰撞顺序按碰撞时间排序处理。

核子-核子碰撞截面采用Cugnon参数化的自由核子-核子截面：在低动量区（$p_{\text{lab}} < 0.3$ GeV/c），$\sigma_{nn/pp} = 60$ mb，$\sigma_{np} = 180$ mb；高能量时截面随束能减小。介质内NN截面通过引入密度依赖修正因子 $\sigma_{NN}^{\text{med}} = (1 + \eta\sqrt{s}\cdot\rho/\rho_0)\sigma_{NN}^{\text{free}}$ 来考虑核介质效应。

泡利阻塞采用Wigner分布函数方法计算。对于碰撞后动量为 $\vec{p}_i'$ 的核子，其泡利阻塞概率为：

$$P_\tau(\vec{r}_i, \vec{p}_i') = 4 \sum_{j \in \tau, j\neq i} \exp\left[-\frac{(\vec{r}_i - \vec{r}_j)^2}{2\sigma_r^2} - \frac{(\vec{p}_i' - \vec{p}_j)^2}{2\sigma_p^2}\right]$$

总碰撞被阻塞的概率为 $P = 1 - (1-P_i)(1-P_j)$，若 $P_\tau > 1$ 则取 $\min(P_\tau, 1)$。2024年有研究提出了改进的PB-W*方法，对动量空间因子进行N个态的平均以减小涨落。

### 4.5 碎片识别算法

碎片识别是多核子转移反应计算中的关键环节，直接决定产物截面的计算结果。imQMD在模拟结束时采用最小生成树（Minimum Spanning Tree, MST）算法或coalescence模型识别碎片。

MST方法的基本原理是将粒子视为图的顶点，边表示邻近关系，算法找到连接所有顶点的无环子图使总边权重最小。在粒子对的质心系中，两个核子若同时满足相对距离 $r_{\text{clust}} \leq r_{\text{max}}$（典型值2.5-3.5 fm）和相对动量 $p_{\text{rel}} \leq P_{\text{max}}$（典型值300 MeV/c）的条件，则归入同一团簇。实现上使用并查集数据结构，通过Find和Union操作高效构建团簇。

coalescence模型采用类似的双判据：空间判据 $\Delta r^2 \leq R^2$（$R \approx 4.0$ fm）和动量判据 $\Delta p^2 \leq P^2$（$P \approx 300$ MeV/c）。对形成后，下一个粒子通过满足相同的 $\Delta r$ 和 $\Delta p$ 条件加入，对动量和对位置分别更新为 $\vec{p}_{\text{pair}} = \vec{p}_1 + \vec{p}_2$ 和 $\vec{r}_{\text{pair}} = (\vec{r}_1 + \vec{r}_2)/2$。

对于三体事件（如¹⁹⁷Au+¹⁹⁷Au三元裂变），需要额外的事件筛选判据：质量平衡要求 $A_P + A_T - 70 \leq A_1 + A_2 + A_3 \leq A_P + A_T$，纵向动量和要求 $\sum_{i=1}^{3}\vec{P}_{\text{long}}^{(i)} > 0.8 P_0$，横向动量和要求 $\sum_{i=1}^{3}\vec{P}_{\text{trans}}^{(i)} < 0.04 P_0$。

产生截面通过碰撞参数积分计算：$\sigma(A_i) = 2\pi\int_0^{b_{\text{max}}} b \cdot P(A_i, b) \, db$，其中 $P(A_i, b) = N(A_i, b)/N_0(b)$ 为碰撞参数 $b$ 处产生碎片 $A_i$ 的概率。

### 4.6 与统计衰变模型的耦合

imQMD动力学演化产生的是激发态初级碎片，需要与统计衰变模型耦合以获得最终产物分布。常用方案为imQMD+GEMINI两阶段方法。

GEMINI统计衰变模型（Charity等，Nucl. Phys. A 483, 371, 1988）用于描述初级碎片的退激发过程，包括轻粒子蒸发、不对称裂变等模式。碎片激发能计算为碎片在体坐标系中的总能量减去对应基态结合能：$E_{\text{frag}} = \sum_i(\vec{p}_i - \vec{p}_c)^2/(2m) + U - E_{\text{bind}}^{\text{gs}}$。

GEMINI的关键参数设置包括：能级密度类型aden_type = -23（Grimes case B修改形式），中间质量碎片裂变模式imf_option = 2（考虑所有不对称分裂），最小中间质量碎片电荷 $Z_{\text{imf,min}} = 5$（Z≤5的粒子按轻粒子蒸发处理）。imQMD与GEMINI的最优切换时间约为500 fm/c，此时初级碎片基本形成但仍保持较高激发能。未测量核的质量取自Weizsäcker-Skyrme (WS4)质量模型预测值。

部分研究也采用imQMD+HIVAP方案，HIVAP统计蒸发模型用于计算剩余碎片截面，在N=126丰中子核产生的研究中有所应用。

### 4.7 Q值效应处理

在¹³⁶Xe+²⁰⁸Pb多核子转移反应研究中，赵凯等人提出了imQMD+Q方法，将不同核子数转移前后形成的类弹与类靶碎片按Q值进行修正。考虑Q值效应后（imQMD+Q），能更好地再现实验质量分布和总动能损失（TKEL）分布，且Q值效应明显增加了比²⁰⁸Pb更丰中子的N=126初级碎块的生成截面。

### 4.8 代码获取与版本

imQMD代码可从王宁教授的个人网站 http://www.imqmd.com/code/ 下载。模型经历了多个版本演进：imQMD（2002年原始版本）、imQMD-II（2004年引入SkM*和SLy系列参数）、imQMD-v2.2（2014-2016年，参数集IQ3a，为目前多核子转移反应计算中使用的最新版本）。WS4质量表可从 http://www.imqmd.com/mass/WS4.txt 获取。此外，Geant4中也包含了QMD的实现（G4QMD），其源代码在GitHub上公开，包含G4QMDCollision、G4QMDGroundStateNucleus、G4QMDMeanField等8个核心模块，可作为理解QMD代码结构的参考。

---

## 五、计算流程总结

imQMD计算多核子转移反应的完整流程可概括为以下步骤：

第一步，基态核制备。根据RMF密度分布抽样核子位置，根据局部密度近似费米动量抽样核子动量，通过能量最小化调整使结合能和电荷半径满足实验值。预演化6000 fm/c检验稳定性，从数千个采样中精选20个弹核和20个靶核。

第二步，碰撞构型设置。弹核与靶核初始间距50 fm，各自绕质心旋转随机欧拉角，碰撞参数从0到12-14 fm取值，步长1 fm，每碰撞参数模拟100,000+事件。

第三步，动力学演化。采用Euler方法积分哈密顿运动方程，时间步长1 fm/c，每步施加相空间占据约束。同时处理核子-核子二体碰撞，采用Cugnon参数化截面和Wigner分布泡利阻塞。

第四步，碎片识别。模拟结束时采用MST算法或coalescence模型识别碎片，基于核子间相对距离（≤2.5-3.5 fm）和相对动量（≤300 MeV/c）判据构建团簇。

第五步，统计衰变。将激发态初级碎片输入GEMINI或HIVAP统计衰变模型，模拟轻粒子蒸发和裂变退激发，获得最终产物分布。最优切换时间约500 fm/c。

第六步，截面计算。通过碰撞参数积分计算产生截面 $\sigma(A_i) = 2\pi\int_0^{b_{\text{max}}} b \cdot P(A_i, b) \, db$，与实验数据进行对比。

---

## 六、模型对比与研究脉络

imQMD与标准QMD和BUU模型的主要区别在于：BUU模型基于单粒子分布函数，无法产生碎片涨落；标准QMD模型使用固定波包宽度且无表面能量项，核基态描述较差；imQMD通过三大改进实现了从⁶Li到²⁰⁸Pb基态性质的良好描述，并能自洽处理动力学效应、同位旋效应和弹靶质量不对称效应。

从研究脉络来看，imQMD在多核子转移反应中的应用经历了以下发展：2002-2004年王宁等提出imQMD模型并引入IQ1参数集；2013年李程、田俊龙等确定IQ3为最优参数集；2016年李程等首次将imQMD应用于¹³⁶Xe+²⁰⁸Pb多核子转移反应；2017年姚宏、王宁将imQMD+GEMINI应用于⁸⁶Kr+⁶⁴Ni；2020年赵凯等研究Q值效应；2021年利用反应时间分析丰中子核产生机制；2025年SSPMA发表系统综述，涵盖同位旋输运、能量耗散和反转移三大机制。

---

## 七、局限性

本研究存在以下局限性。第一，碎片识别算法的具体实现细节在多数应用论文中未被详细描述，本文中的MST和coalescence判据参数主要来自QMD类模型的通用文献和部分可获取的论文，imQMD原始代码中的具体实现可能有所不同，需参考王宁2002和2004年的原始论文全文确认。第二，部分Physical Review C论文的全文由于PDF访问限制未能完全获取，相关技术细节主要来自arXiv预印本HTML版本和开放获取论文。第三，不同研究团队使用的imQMD版本和参数集存在差异，实际应用时需参考具体论文的设置。第四，数值积分方法在不同文献中记载不一，⁸⁶Kr+⁶⁴Ni论文明确使用Euler方法，但其他版本可能采用蛙跳法或Velocity-Verlet法，需进一步核实。

---

## References

1. [An Improved Quantum Molecular Dynamics Model (arXiv)](https://arxiv.org/abs/nucl-th/0201079)
2. [Further Development of the Improved QMD Model (arXiv)](https://arxiv.org/abs/nucl-th/0402066)
3. [Multinucleon transfer in 136Xe+208Pb (Physical Review C)](https://doi.org/10.1103/PhysRevC.93.014618)
4. [Microscopic dynamics simulations of MNT in 86Kr+64Ni (Physical Review C)](https://doi.org/10.1103/PhysRevC.95.014607)
5. [Isospin equilibration in MNT (Physical Review C)](https://doi.org/10.1103/PhysRevC.99.034619)
6. [Production mechanism of neutron-rich nuclei in 136Xe+198Pt (Physics Letters B)](https://doi.org/10.1016/j.physletb.2020.135697)
7. [Production of neutron-rich N=126 nuclei (Physics Letters B)](https://www.sciencedirect.com/science/article/pii/S0370269321000411)
8. [Q值效应对多核子转移反应136Xe+208Pb的影响 (原子核物理评论)](http://npr.xml-journal.net/article/doi/10.11804/NuclPhysRev.37.2020022)
9. [ImQMD Application to Ternary Breakup 197Au+197Au (Universe)](https://www.mdpi.com/2218-1997/8/11/555)
10. [Effects of nuclear surface dynamics on fusion and MNT (Physical Review C)](https://doi.org/10.1103/PhysRevC.112.034601)
11. [Determination of NN interaction in ImQMD (Chinese Physics C)](https://iopscience.iop.org/article/10.1088/1674-1137/37/11/114101)
12. [改进的量子分子动力学模型在多核子转移反应中的应用 (中国科学)](https://doi.org/10.1360/SSPMA-2024-0543)
13. [Progress of QMD model and its applications (arXiv)](https://arxiv.org/abs/2005.12877)
14. [利用反应时间分析多核子转移反应中丰中子核的产生机制 (原子核物理评论)](http://npr.xml-journal.net/article/doi/10.11804/NuclPhysRev.38.2020045)
15. [多核子转移反应的微观动力学模拟 (学位论文)](https://www.zhangqiaokeyan.com/academic-degree-domestic_mphd_thesis/020311804309.html)
16. [ImQMD代码官方网站](http://www.imqmd.com/code/)
17. [Inverse transfer mechanisms in MNT (Physical Review C)](https://doi.org/10.1103/PhysRevC.110.044615)
18. [QMD模型的改进及其在低能重离子反应中的应用 (原子核物理评论)](http://npr.xml-journal.net/article/doi/10.11804/NuclPhysRev.21.04.374)
