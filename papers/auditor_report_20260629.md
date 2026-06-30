# ImQMD 代码审计报告 — vs 论文模型规格

**审计任务**: `t_d2f4132e`  
**审计对象**: `mnt_sim/imqmd/` 全部模块  
**规格文档**: `papers/imqmd_model_spec.md`, `papers/imqmd_formulas.md`, `papers/imqmd_literature.md`  
**审计日期**: 2026-06-29  
**审计人**: Code Auditor (`auditor`)

---

## 综合结论

**❌ 需修改**

当前 `mnt_sim/imqmd/` 实现了 ImQMD 家族的核心框架，参数表、高斯波包密度、RK4 传播、几何 NN 碰撞、MST/iso-MST 碎片识别、质量表读取与衰变通道等基础设施均已就位，且 `tests/test_static.py`、`tests/test_collisions.py` 可通过。但存在 **1 项 CRITICAL** 与 **3 项 HIGH** 问题，导致代码并未真正使用论文规格的 EDF 势能面，同时 Pauli 阻塞的实现加入了非物理的硬地板，会系统性压低有效碰撞率。在这些问题修复之前，代码产出的物理结果不能视为 ImQMD 论文模型的一致实现。

---

## 1. 按维度逐项审查

### 1. 波包表示与宽度参数化

| 检查项 | 标记 | 文件:行号 | 说明 |
|---|---|---|---|
| 高斯密度 Eq. (W2002-5) | ✓ | `nucleus.py:91`, `grid_edf.py:506` | 归一化高斯权重、指数形式正确。 |
| Wigner 相空间密度 | ⚠ | `collisions.py:85-87` | 仅用于 Pauli 阻塞，缺少 `1/(πħ)^3` 归一化（见维度 5）。 |
| `σ_r σ_p = ħ/2` | ✓ | `collisions.py:84` | 用 `HBAR_C/(2σ_r)` 得到 `σ_p`，单位自洽。 |
| 宽度公式 `σ_r = σ_0 + σ_1 A^{1/3}` | ✓ | `initializer.py:99-106` | `compute_sigma_r` 实现与 IQ1/IQ2/IQ3 表一致。 |
| 反应系统宽度过渡 | ✗ | 未实现 | 论文要求弹靶接触后使用统一的系统相关宽度；代码未定义该过渡。 |

**结论**: 高斯表示正确，但反应动力学中的系统宽度过渡缺失。

---

### 2. Skyrme EDF 能量密度泛函

| 检查项 | 标记 | 文件:行号 | 说明 |
|---|---|---|---|
| IQ1/IQ2/IQ3/IQ3a/IQ3b/SkP* 参数表 | ✓ | `skyrme.py:38-45` | 数值与 `imqmd_model_spec.md` 表 9 完全一致。 |
| 现代 EDF Eq. (W2014-5) — Grid 路径 | ✓ | `grid_edf.py:106-147` | `skyrme_bulk`、`symmetry_term`、`surface_term`、`surface_symmetry`、`coulomb_exchange` 公式正确。 |
| 体项 `α/2·ρ²/ρ₀` | ✗ | `skyrme.py:75` | `skyrme_potential` 使用 `np.sum(rho)` 而非 `np.sum(rho**2)`。 |
| 体项 `β/(γ+1)·ρ^{γ+1}/ρ₀^γ` | ✗ | `skyrme.py:76` | 使用 `np.sum(rho**gamma)` 而非 `np.sum(rho**(gamma+1))`。 |
| `g_τ` 项 `g_τ·ρ^{η+1}/ρ₀^η` | ✗ | `skyrme.py:77` | 使用 `np.sum(rho**eta)` 而非 `np.sum(rho**(eta+1))`。 |
| 对称能体项 | ✗ | `skyrme.py:121-122` | 使用 `np.sum(rho*delta**2)` = `sum((ρ_n-ρ_p)²/ρ)`，应为 `sum((ρ_n-ρ_p)²)`。 |
| 表面/表面-对称项 | ⚠ | `skyrme.py:124-156` | 用经验高斯对势近似 `g_sur/7.0`、`-C_s κ_s/7.0`，不是论文的 `(∇ρ)²` 形式。 |
| Coulomb 直接项 + Slater 交换 | ✓ | `skyrme.py:80-112`, `grid_edf.py:149-153` | 有限高斯直接 Coulomb 与 `ρ_p^{4/3}` 交换项均实现。 |
| 能量-力一致性 | ⚠ | `propagator.py:66-104` | `_centroid_density_force` 与错误的 `skyrme_potential` 自洽，导致 `test_static.py::test_physical_energy_force_consistency` 无法发现错误。 |

**结论**: GridEDF 路径的 EDF 公式正确；但 `SkyrmeEDF` 的非 Grid 路径（能量诊断、`energy_components`、无 Grid 传播）存在幂次错误，会输出与论文不符的能量分量。

---

### 3. 基态初始化

| 检查项 | 标记 | 文件:行号 | 说明 |
|---|---|---|---|
| 电荷半径 Eq. (W2014-8) | ✓ | `initializer.py:67` | `R_c` 公式实现正确。 |
| 中子皮厚度 Eq. (W2014-9) | ✓ | `initializer.py:68` | `ΔR_np = 0.9I - 0.03` 正确。 |
| 质子/中子硬球半径转换 | ✓ | `initializer.py:69-71` | 按论文 `⟨r_c²⟩`、`R_p`、`R_n` 步骤实现。 |
| 位置采样 `R_p-w_r`、`R_n-w_r`，`w_r=0.8 fm` | ✓ | `initializer.py:86-87` | 实现并减去 0.8 fm。 |
| 局部 Fermi 动量 Eq. (init-pf) | ✓ | `initializer.py:46-47`, `grid_edf.py:158-159` | `p_F = ħ(3π²ρ_q)^{1/3} - w_p` 正确。 |
| 相空间距离 `|r_ij||p_ij| ≥ 255` | ✓ | `initializer.py:245`, `collisions.py:396` | 初始化与 Fermi 约束均使用 255 fm·MeV/c。 |
| 质心修正 | ✓ | `initializer.py:95`, `propagator.py:230-235` | 位置与动量均去质心。 |
| 能量验收 `BE ± 0.05 MeV` | ⚠ | `initializer.py:290` | 实际使用 `0.5 MeV`，比论文宽松一个量级。 |
| 能量偏移 `energy_offset` | ⚠ | `nucleus.py:136` | 用经验结合能校准总能量，掩盖 EDF 自身精度。 |
| 弛豫时间 | ✓ | `initializer.py:268-283` | 轻核 600 fm/c 阻尼弛豫，与论文思路一致。 |

**结论**: 初始化流程框架与论文一致，但能量验收过宽且依赖 `energy_offset` 校准。

---

### 4. 运动方程 / 传播算法

| 检查项 | 标记 | 文件:行号 | 说明 |
|---|---|---|---|
| Hamilton 方程 `ṙ = ∂H/∂p`, `ṗ = -∂H/∂r` | ✓ | `propagator.py:165-180` | `drdt = p/M_N`，`dpdt` 由各势能梯度组成。 |
| RK4 四阶积分 | ✓ | `propagator.py:183-235` | 标准 RK4 实现正确。 |
| 能量守恒测试 | ✓ | `tests/test_static.py:49-61` | 带 stabilizer 2000 fm/c 漂移约 0.28%。 |
| 无弹簧稳定性 | ✓ | `tests/test_static.py:64-82` | 500 fm/c 无弹簧漂移 < 8%，满足当前测试门限。 |
| 静态 stabilizer 弹簧 | ⚠ | `skyrme.py:163-176`, `propagator.py:133-154` | 论文未要求；用于稳定性测试，非物理外力。 |

**结论**: RK4 与 Hamilton 结构正确；stabilizer 弹簧属于测试辅助项，不影响论文模型主路径。

---

### 5. NN 碰撞

| 检查项 | 标记 | 文件:行号 | 说明 |
|---|---|---|---|
| 几何碰撞判据 / 最近距离 | ✓ | `collisions.py:340-343` | `t_*` 裁剪到 `[0, dt]`，`d_min` 计算正确。 |
| Cugnon 自由截面公式 | ✓ | `collisions.py:47-50` | pp/nn、np 公式与论文一致。 |
| 低能截面截断 | ✓ | `collisions.py:41-42` | 60 mb / 180 mb 截断，与 Chen-Zhang-Li 2021 一致。 |
| 介质修正 `(1 - ξ ρ/ρ₀)` | ✓ | `collisions.py:58` | 实现 `1 - 0.2 ρ/ρ₀`。 |
| 介质因子下限截断 | ⚠ | `collisions.py:58` | 强制下限 0.2，论文未要求。 |
| Zhang 2012 三种截面情形 | ✗ | 未完整实现 | 论文要求 case1/2/3 与等效截面 `σ'`；代码仅提供单一通道因子。 |
| Pauli 阻塞概率公式 `P_block = 1-(1-P_i)(1-P_j)` | ✓ | `collisions.py:126`, `218` | 公式结构正确。 |
| Wigner 占据 `P_i` 归一化 | ✗ | `collisions.py:89` | 缺少论文 `1/(πħ)^3` 归一化，导致数值占据远小于 1。 |
| 任意 Fermi  floor = 0.75 | ✗ | `collisions.py:95-98` | 对费米面内/附近状态硬编码最低 0.75 占据，不是论文公式，会过度阻塞碰撞。 |
| PB(W*) 宽化核 | ✓ | `collisions.py:129-148` | 提供 `pauli_blocking_probability_v2`，对应 Chen 2024 思路。 |
| CoMD Fermi 约束 | ✓ | `collisions.py:380-425` | `|r_ij||p_ij| ≥ 255` 并尝试保持动能。 |
| 碰撞避免重复计数 | ✓ | `collisions.py:325`, `332`, `371` | 每步内已参与碰撞的核子不再参与。 |

**结论**: 几何截面与介质缩放正确；Pauli 阻塞实现偏离论文，存在缺失归一化与非物理 floor。

---

### 6. 碎片识别

| 检查项 | 标记 | 文件:行号 | 说明 |
|---|---|---|---|
| MST `R_cut=3.5 fm`, `P_cut=250 MeV/c` | ✓ | `fragments.py:82-86` | 默认值与 Zhang 2012 一致。 |
| iso-MST `R_nn=R_np=6 fm`, `R_pp=3 fm` | ✓ | `fragments.py:105-111` | 函数签名支持论文值，默认保守。 |
| 传递性连通分量 | ✓ | `fragments.py:26-53` | 使用 SciPy/回退 DFS 实现。 |
| 轻复合核表面聚合 | ⚠ | `fragments.py:191-243` | 用两体 `r_cut/p_cut` 搜索，未实现论文 Jacobian `R_im·P_im ≤ h_0`。 |
| 碎片激发能 `E* = E_int - E_ground` | ✓ | `fragments.py:246-288` | 去质心系，用 GridEDF 或子系统能量减基态能量。 |

**结论**: MST/iso-MST 与激发能计算基本合规；轻核聚合仅为近似实现。

---

### 7. 退激发

| 检查项 | 标记 | 文件:行号 | 说明 |
|---|---|---|---|
| Mexcess95 质量表读取 | ✓ | `decay.py:16-39` | 从 `hivap_fortran/Mexcess95.dat` 读取。 |
| 液滴 fallback 结合能 | ✓ | `decay.py:42-61` | 提供质量表缺失时的 fallback。 |
| 中子/质子/α 分离能 + Coulomb 位垒 | ✓ | `decay.py:92-138` | 用质量差与液滴 fallback 估算。 |
| 多通道分支比 n, p, α, γ, fission | ✓ | `decay.py:190-242` | 实现 HIVAP 通道模式。 |
| Weisskopf-Ewing 积分公式 | ⚠ | `decay.py:179-188` | 用简化相位空间 `exp(2√(aE))` 近似，未执行论文所示的能量积分。 |
| Bohr-Wheeler 裂变竞争 | ⚠ | `decay.py:158-176` | 用简化指数比，未按论文积分鞍点能级密度。 |

**结论**: 退激发框架与 HIVAP 耦合模式一致，但核心公式为简化版，不能替代 HIVAP/GEMINI 验证。

---

### 8. 能量-力一致性

| 检查项 | 标记 | 文件:行号 | 说明 |
|---|---|---|---|
| GridEDF 能量与力的 `η` 缩放一致性 | ✓ | `grid_edf.py:238`, `347` | 力与能同时乘以 `self.eta`，内部自洽。 |
| 能量-力对应论文 EDF | ✗ | `grid_edf.py:347`, `initializer.py:266` | `eta=0.5` 使势能与力均为论文 EDF 的 50%，不符合论文参数。 |
| `SkyrmeEDF` 能量-力自洽 | ⚠ | `skyrme.py:70-176`, `propagator.py:66-131` | 两者自洽，但对应的不是论文 EDF。 |

**结论**: 主路径的能量-力自洽，但 Hamiltonian 被全局缩放 0.5，偏离论文规格。

---

## 2. 问题汇总表

| 编号 | 严重程度 | 维度 | 问题简述 | 位置 | 修复建议 |
|---|---|---|---|---|---|
| C1 | **CRITICAL** | EDF / 能量-力 | `GridEDF.eta=0.5` 全局缩放势能与力，实际 EDF 强度为论文一半 | `initializer.py:266,291`; `grid_edf.py:54,347,238` | 移除 `eta` 缩放（默认 1.0），用正确参数或修正 GridEDF 实现使 `eta=1.0` 复现结合能；若必须保留，将其显式写入参数表并说明不是论文 IQ3a。 |
| H1 | **HIGH** | EDF | `SkyrmeEDF.skyrme_potential` 幂次错误：`sum(rho)` 应为 `sum(rho**2)`，`sum(rho**gamma)` 应为 `sum(rho**(gamma+1))` 等 | `skyrme.py:75-77` | 改为与 `grid_edf.py:106-116` 一致的体项；或删除/标记为内部调试用途。 |
| H2 | **HIGH** | EDF | `SkyrmeEDF.symmetry_energy` 公式错误：`sum(rho*delta**2)` 应为 `sum((rho_n-rho_p)**2)` | `skyrme.py:121-122` | 改为 `p.c_sym/(2*p.rho0) * np.sum((rho_n - rho_p)**2)`。 |
| H3 | **HIGH** | 碰撞 / Pauli | Pauli 占据缺少 `1/(πħ)^3` 归一化，并用 `fermi_floor=0.75` 人为抬高阻塞 | `collisions.py:89,95-98` | 加入 Wigner 核归一化；移除或重新论证 `fermi_floor`，按论文纯 Wigner/Husimi 实现。 |
| M1 | **MEDIUM** | EDF | `SkyrmeEDF` 表面项与表面-对称项用高斯对势近似，非 `(∇ρ)²` 形式 | `skyrme.py:124-156` | 若保留 centroid 路径，改用有限差分计算 `∇ρ` 并匹配 GridEDF 公式。 |
| M2 | **MEDIUM** | 初始化 | 能量验收阈值 `0.5 MeV` 比论文 `0.05 MeV` 宽 10 倍 | `initializer.py:290` | 收紧到 `0.05 MeV`；若因此导致初始化失败，先修复 C1/H1/H2。 |
| M3 | **MEDIUM** | 初始化 | `energy_offset` 用经验结合能校准，掩盖 EDF 真实精度 | `nucleus.py:136` | 仅在调试模式保留；正式发布时置 0，让 EDF 自身给出结合能。 |
| M4 | **MEDIUM** | 碎片 | 轻复合核聚合用两体截断，未实现 Jacobian `R_im·P_im ≤ h_0` | `fragments.py:191-243` | 按 Wei 2014 Eq. (12) 实现 Jacobian 相空间积。 |
| M5 | **MEDIUM** | 退激发 | Weisskopf/Bohr-Wheeler 为简化近似，未按论文积分 | `decay.py:158-188` | 文档化“示意实现”；与 HIVAP 输出做系统对比标定。 |
| L1 | **LOW** | 碰撞 | 介质截面因子强制下限 0.2 | `collisions.py:58` | 移除 clip 下限或说明高密度行为依据。 |
| L2 | **LOW** | 代码 / 初始化 | 警告文本写“after 10 retries”，实际循环仅 3 次 | `initializer.py:241,297-298` | 统一文本与实际循环次数。 |
| L3 | **LOW** | 波包 | 反应系统宽度从弹靶独立到系统统一的过渡未定义 | 未实现 | 在反应事件设置中补充宽度过渡规则。 |
| L4 | **LOW** | 碰撞 | Zhang 2012 的 case1/2/3 与等效 `σ'` 未完整暴露 | `collisions.py:54-65` | 增加可选的三种截面模式供比较。 |

---

## 3. 统计表

| 类别 | 计数 | 编号 |
|---|---|---|
| ✓ 通过 | 24 项 | 见各维度表格 |
| ⚠ 偏差 / 近似 | 6 项 | EDF 表面项、能量验收、energy_offset、轻核聚合、退激发近似、stabilizer 弹簧 |
| ✗ 未达标 | 6 项 | C1, H1, H2, H3, M2(已计入), L3, L4 |
| **CRITICAL** | 1 | C1 |
| **HIGH** | 3 | H1, H2, H3 |
| **MEDIUM** | 5 | M1, M2, M3, M4, M5 |
| **LOW** | 4 | L1, L2, L3, L4 |

---

## 4. P0 / P1 / P2 / P3 优先级路线图

### P0 — 发布前必须修复（阻塞级）

1. **C1**: 决定 `GridEDF.eta` 的物理意义。若它只是调试参数，应默认 `eta=1.0`；若发现 `eta=1.0` 时结合能偏差大，则意味着 EDF 参数、网格积分或密度归一化存在更底层错误，需一并修复。不能让生产运行使用 50% 的势能面。

### P1 — 应立即修复（高影响）

2. **H1/H2**: 修正 `SkyrmeEDF.skyrme_potential` 与 `symmetry_energy` 的幂次/公式，使其与 `GridEDF` 一致。否则 `energy_components`、无 Grid 传播、用户诊断都会给出错误能量分解。
3. **H3**: 重写 Pauli 阻塞占据计算：加入 `1/(πħ)^3` 归一化，移除 `fermi_floor=0.75`，按论文提供标准 Wigner 阻塞；将当前 floor 版本保留为可选开关并文档化。

### P2 — 重要但非阻塞

4. **M1**: 将 `SkyrmeEDF` 的表面项改为 `(∇ρ)²` 形式，或明确标注为经验近似。
5. **M2**: 将初始化能量验收收紧到 `BE ± 0.05 MeV`；若初始化失败率高，回到 P0 先修正 EDF。
6. **M3**: 移除生产代码中的 `energy_offset`，仅在调试/校准模式保留。
7. **M4/M5**: 轻核聚合与退激发的简化公式需文档化，并制定与 HIVAP 的对比验证计划。

### P3 — 建议项

8. **L1/L4**: 完善介质截面选项（下限 clip、case1/2/3）。
9. **L2**: 修正 retry 警告文本。
10. **L3**: 补充反应宽度过渡规则。

---

## 5. 验证命令

```bash
cd "/home/zhuhaofan/work/agent work/mnt-sim"
python tests/test_static.py   # 当前通过，但验证的是缩放后的 Hamiltonian
python tests/test_collisions.py
```

建议在修复 C1/H1/H2/H3 后新增以下回归测试：

- `GridEDF.eta=1.0` 时 `^{40}Ca` 结合能应落在 `empirical_binding_per_nucleon ± 0.5 MeV/A` 内。
- `SkyrmeEDF.skyrme_potential(rho)` 与 `GridEDF.skyrme_bulk` 在相同 centroid 密度向量上给出可比的体项能量。
- Pauli 阻塞：在冷核中随机弹性散射后，`accepted_collisions / attempted_collisions` 接近 0；在高激发热核中接受率上升。

---

## 6. 最终结论

代码结构完整、可运行、测试通过，但**核心物理势能面与 Pauli 阻塞实现偏离论文规格**。在修复 **C1（GridEDF 势能面被缩放 0.5）**、**H1/H2（SkyrmeEDF 公式错误）**、**H3（Pauli 阻塞非物理 floor）** 之前，不能把该实现宣称为符合 Wang 2014 / Zhang 2012 / Zhao 2016 等论文的 ImQMD 模型。
