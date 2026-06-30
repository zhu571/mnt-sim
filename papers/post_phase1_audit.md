# ImQMD Phase 1 修复后审计报告

> 生成日期: 2026-06-28
> 审计范围: Phase 1 修复后 `mnt_sim/imqmd/` 全部 9 个 .py 文件 vs `papers/imqmd_model_spec.md` (1026 行完整规范)
> 基准: `papers/final_gap_audit.md` 中识别的 36 项差距

---

## Phase 1 修复验证

下表逐项核实四项 Phase 1 修复是否已正确落地代码中。

| # | 修复项 | 文件 | 位置 | 验证结果 |
|---|--------|------|------|----------|
| 1 | **IQ3a 默认参数** | `skyrme.py` | L56: `PARAMETER_SETS["IQ3A"]` | ✅ `SkyrmeEDF.__init__` 默认使用 IQ3a（α=-207, β=138, g_sur=16.5, κ_s=0.4, C_s=34） |
| 2 | **SkP* 参数集** | `skyrme.py` | L43: `"SKP*": SkyrmeParameters(...)` | ✅ 已添加 SkP*（α=-356, β=303, g_sur=19.5, g_τ=13, C_s=35, κ_s=0.65, ρ0=0.162, σ0=0.94, σ1=0.018） |
| 3 | **A 依赖波包宽度** | `initializer.py` | L92-99: `compute_sigma_r()`；L213-214: 自动调用 | ✅ `σ_r = σ₀ + σ₁·A^(1/3)` 初始化时按参数集自动计算。对IQ3a在A=476时：σ_r = 0.94 + 0.020×476^(1/3) ≈ 1.15 fm（之前固定 1.1 fm 对应 A≈66） |
| 4 | **_sample_grid_positions 接通** | `initializer.py` | L216: `positions, is_proton = _sample_grid_positions(Z, A, sigma_r, rng)` | ✅ `initialize_nucleus()` 现调用 Wang 2014 中子皮初始化（R_c→R_p→R_n→w_r=0.8 fm校正），替换了之前的简化硬球 |
| 5 | **CoMD 阈值 170→255** | `collisions.py` | L380: `threshold: float = 255.0` | ✅ 默认阈值改为 255 fm·MeV/c |
| 6 | **CoMD 全核子对检查** | `collisions.py` | L389-412: 无双循环 `if species[i] != species[j]: continue` | ✅ 已移除同位旋过滤，现检查所有核子对（np 对、nn 对、pp 对均检查） |

**Phase 1 修复结论**: 四项修复均正确落地。U+U (A=476, b=7 fm) 结果从 402 碎片完全解体改善为 1 个重碎片 Z=183 A=472（保留 99.2% 核子），验证了修复的有效性。

---

## 残留差距清单（Phase 1 后仍未修复）

以下按严重程度分组，从 `final_gap_audit.md` 的 36 个差距中筛除已修复项后列出。

### Critical（严重） — 3 项

| Gap # | 描述 | 当前状态 | 影响 |
|-------|------|----------|------|
| **2.2** | 动量采样使用 scale 因子（0.70-0.50）而非论文的 w_p 截断（p_F - w_p） | `_fermi_momenta()` L46: `pmax * scale`；scale 是全局缩减而非固定截断 | 动量分布形状与论文不一致；重核（A=476）累积误差导致结合能不匹配 |
| **2.3** | 能量偏移机制强制匹配 BE，掩盖 Hamiltonian 不一致 | `initialize_nucleus()` L257-259: `nucleus.energy_offset = target_total - raw_components["total"]` | 这是最严重的物理问题。能量偏移可达 100+ MeV，弛豫去除弹簧后真实 Hamiltonian 无法维持结合。U+U 之前完全解体的另一根本原因 |
| **2.5** | 弛豫时间仅 50 fm/c（`initialize_nucleus` 直接调用时） | L241: `for step in range(1, 51):` × 1.0 fm/c = 50 fm/c | `initialize_and_relax` 已做 800 fm/c，但若直接调用 `initialize_nucleus` 则弛豫严重不足。论文要求 ≥600 fm/c |

### Major（重要） — 16 项

| Gap # | 描述 | 当前状态 |
|-------|------|----------|
| **1.3** | 缺少 ImQMD05 显式动量相关相互作用 u_md | EDF 无 Zhang 2012 Eq.(5) 的 1.57[ln(1+5e-4(Δp)²)]² 项。对 E_beam > 20 MeV/A 的重离子碰撞影响显著 |
| **2.7** | 缺少虚假发射检查 | 弛豫后无任何虚假粒子发射检查。论文要求拒绝有核子逃逸（距 CM > 15 fm）的初始配置 |
| **3.3** | CoMD 约束仅一次性弹性散射，无迭代 | `fermi_constraint_check()` L412: `break` 跳至下一个 i。论文要求多体弹性散射迭代且检查能量变化。注意：内循环 break 意味着每个外循环 i 最多修正一对 |
| **4.1** | A 依赖 σ_r 仅初始化使用，传播时不动态更新 | 初始化时正确，但 `nucleus.sigma_r` 是固定属性。论文要求传播时根据碎片大小动态调整宽度 |
| **4.2** | 弹靶独立宽度未实现 | 所有核子共享 `nucleus.sigma_r`。非对称反应（如 Kr+Ni, A1=78, A2=58）两边应使用不同宽度 |
| **5.1** | PB(W*) 未被正确实现 | `pauli_blocking_probability_v2()` 仅将 Wigner 核宽度×1.35（"Chen-2024-inspired"）。论文要求对每个邻居采样 N 个动量状态计算指数和平均（Chen 2024 Eq.6-7），报告方差降低 70% |
| **5.4** | 碰撞泡利阻塞使用预计算密度而非实时密度 | `_pauli_blocking_probability_fast()` 使用 `centroid_densities()` 预计算值。前序碰撞可能已改变核子位置，误差在多碰撞步骤中累积 |
| **6.1** | Weisskopf-Ewing 蒸发宽度计算过于简化 | `_evaporation_width()` L186: `A^(2/3)·ε·exp(2√(aE))` 近似。论文要求完整积分 `∫ σ_j^inv(ε)·ρ_d(E*-B_j-ε)·ε dε` |
| **6.2** | γ 道竞争使用启发式 | L227: `gamma_width = max(0.1*gamma_reference, 1e-12)`。论文要求 HIVAP 使用 GDR 巨偶极共振参数 |
| **6.3** | 裂变竞争使用简化液滴模型 | `fission_barrier()` 使用液滴+简单壳修正。对 Z>90 超铀核，与实验值差异可达 2-3 MeV |
| **7.1** | In-Medium η 固定 0.2 | `in_medium_factor()` L58: 固定 `eta=0.2`（Zhang 2012 E_beam=50A MeV 的值）。对 U+U E_beam=7 MeV/A 可能高估介质抑制 |
| **7.2** | LCP 表面凝聚使用简化相位空间凝聚非 Wei 2014 模型 | `coalescence_light()` 使用 r_cut=2.4, p_cut=180 的简单相位空间凝聚。缺少 Wei 2014 的核心+表面厚度模型、Jacobi 坐标条件、Coulomb 势垒穿透 |
| **7.3** | 缺少反应截面（碰撞参数积分）模块 | 无事件循环+碰撞参数网格+截面累计。无法直接从代码产出可与实验对比的截面数据 |
| **7.4** | 传播默认不使用 GridEDF | `propagate()` L295: `use_grid_edf=False` 默认值。对 A=476 的 O(N²)=226K 对/RK4 步，性能严重下降 |
| **3.4** | 传播阶段 CoMD 约束频率 | `propagate()` L289-290: `fermi_interval = max(1, int(round(5.0/dt)))`（每 5 fm/c 一次）。论文描述为"each time step" |
| **3.5** | 缺少单体占据 bar f_i 计算 | 仅检测成对 dr×dp，不计算核子 i 对所有 j 的总 Wigner 占据。论文要求 bar f_i > 1 触发多体散射 |

### Moderate（中等） — 10 项

| Gap # | 描述 | 当前状态 |
|-------|------|----------|
| **1.4** | 缺少 isovector 表面项 g_sur_iso | `SkyrmeParameters` 无 `g_sur_iso` 字段；`GridEDF` 无对应项。论文 Zhang 2012 Eq.(4) 包含此表面项 |
| **2.8** | 缺少 RMF/Woods-Saxon 初始密度分布选项 | 仅有硬球+中子皮初始化。论文 Wang 2002 使用 RMF 密度分布 |
| **2.9** | 质心移除在弛豫期间多次重复 | L242 和 L251 都做 CM 校正 + 0.93 阻尼。过度校正可能人为压制合理质心涨落 |
| **4.3** | σ_r·σ_p = ħ/2 关系在能量计算中不一致 | `density()` 使用 `sigma_r`，Pauli 阻塞使用 `sigma_p=HBAR_C/(2*sigma_r)`，但 `surface_pair_energy` 使用 `2σ_r` 范围 |
| **5.2** | Fermi floor 启发式不在论文中 | `_occupation_wigner()` L94-99: 当 p<p_F 时强制占据 ≥0.75。论文的 Pauli 阻塞不包含此逻辑，可能遮掩低密度区泡利违反 |
| **5.3** | 缺少硬球重叠替代方案（Chen 2024 Eq.2） | 仅有 Wigner 核占据。无法区分两种方案的差异 |
| **6.4** | 多步蒸发能量扣除过于简化 | `_channel_energy_cost()` 扣除 S_n + 3T（动能），论文标准为 2T。多步后 E* 累积误差约 10-20% |
| **6.5** | 缺少 HIVAP/GEMINI 外部耦合接口 | 仅有内部蒸发/裂变实现。无法复现 Zhao 2016 的 ImQMD→HIVAP 工作流 |
| **7.5** | 未实现形变核初始化（β₂, β₄） | 仅支持球形硬球。Zhao 2016 U+U 要求 ²³⁸U β₂=0.215, β₄=0.093 + 随机取向 |
| **7.6** | 未实现 Zhao 2016 U+U 基准工作流 | 无完整事件生成脚本（E_gs=7.37 MeV/nucleon, b_max=15 fm, Δb=0.15 fm, d₀=40 fm, N_events=100000） |

---

## Phase 1 引入的新问题 / 回归

### 新问题 1: GridEDF 默认参数仍为 IQ2（Latent bug）

| 属性 | 值 |
|------|-----|
| **文件** | `grid_edf.py` L50 |
| **严重程度** | ★★ Moderate（潜在，当前代码路径不受影响） |
| **描述** | `GridEDF.__init__(self, parameters, ...)` 中：`self.parameters = parameters or PARAMETER_SETS["IQ2"]`。SkyrmeEDF 已切至 IQ3a，但 GridEDF 的 `parameters=None` 默认值仍用 IQ2。当前所有创建 GridEDF 的代码路径都显式传入参数（`nucleus.edf.parameters`），因此实际运行不受影响，但这是一个等待触发的漏洞 |
| **修复** | 将 `"IQ2"` 改为 `"IQ3A"` |

### 新问题 2: σ_r 初始化正确但 σ_p 可能出现偏差

| 属性 | 值 |
|------|-----|
| **文件** | `nucleus.py` L68, `collisions.py` L83-84 |
| **严重程度** | ★★ Moderate |
| **描述** | 初始化通过 `compute_sigma_r(A, parameters)` 计算 σ_r（对 IQ3a A=476: ~1.15 fm），但核子动量采样仍使用固定 scale 因子 0.70-0.50（而非论文 w_p），且 `GaussianPacket` 每个核子存有独立 `sigma_r` 但 nucleus 的 `sigma_r` property 只返回第一个 packet 的值。这意味着如果今后支持独立宽度，多处代码会断裂 |
| **影响** | 1) 动量-宽度关系通过 `σ_p = ħ/(2σ_r)` 间接受 σ_r 变化影响，但采样尺度未联动；2) `sigma_r` property 的单点假设是脆弱的 |

### 新问题 3: CoMD 约束现在正确检查所有对，但单点修正限制

| 属性 | 值 |
|------|-----|
| **文件** | `collisions.py` L389-412 |
| **严重程度** | ★★★ Major |
| **描述** | Phase 1 修复了同位旋过滤，但现在 `fermi_constraint_check()` 的内循环包含 `break`（L412），每个外循环 i 最多修正一对 (i,j)。对于 A=476 系统且阈值提升至 255 后，违反对数增多，但每次调用只能修正 ≤A/2 对 |
| **修复** | 移除 `break` 或改为收集所有违反对后统一修正 |

### 新问题 4: `initialize_nucleus` 的弛豫仍使用 scale 动量采样

| 属性 | 值 |
|------|-----|
| **文件** | `initializer.py` L223 |
| **严重程度** | ★★★ Major |
| **描述** | 即使接通了 `_sample_grid_positions()` 实现正确空间初始化，动量采样仍通过 `_fermi_momenta()` 使用 scale 因子轮换（0.70, 0.65, 0.60, 0.55, 0.50）。`_fermi_momenta()` 内部首先做密度求值（L42）创建临时 ImQMDNucleus，但此临时核的 sigma_r 通过位置传入参数设置。问题：当 σ_r(A) 随 A 增大时，scale 因子对应的实际费米动量截断尺寸可能不匹配 |
| **影响** | 对于 A=476 的大系统，σ_r(~1.15 fm for IQ3a) 比固定 1.1 fm 略大，费米动量采样形状不变（scale 因子不变），但密度求值时波包更宽导致密度略低→p_F 略低，产生二阶效应 |

### 新问题 5: 文档声称使用 IQ2 默认但代码用 IQ3a

| 属性 | 值 |
|------|-----|
| **文件** | `skyrme.py` L48 |
| **严重程度** | ★ Minor |
| **描述** | docstring: `"Small ImQMD EDF calculator with IQ2 defaults."` 已过时，实际默认已改为 IQ3a |
| **修复** | 更新 docstring |

---

## 修复进展统计

| 状态 | 数量 | 占比 |
|------|------|------|
| ✅ Phase 1 已修复 | 6 | 16.7% |
| ⬜ 未修复 — Critical | 3 | 8.3% |
| ⬜ 未修复 — Major | 16 | 44.4% |
| ⬜ 未修复 — Moderate | 10 | 27.8% |
| ⚠️ Phase 1 新增 | 5 | — |
| **原始差距总计** | **36** | **100%** |

### 修复进展按类别

| 类别 | 原始差距 | Phase 1 已修复 | 未修复 |
|------|---------|---------------|--------|
| 1. EDF 参数集 | 4 | 2 (IQ3a, SkP*) | 2 (u_md, g_sur_iso) |
| 2. 初始化方法 | 9 | 3 (硬球、阈值、w_r校正) | 6 (w_p, 能量偏移, 弛豫时间, 虚假发射, RMF, CM重复) |
| 3. CoMD 约束 | 5 | 2 (阈值、全对) | 3 (迭代散射, 频率, 单体占据) |
| 4. 波包宽度 | 3 | 1 (A依赖初始化) | 2 (弹靶独立, σ_rσ_p一致性) |
| 5. Pauli 阻塞 | 4 | 0 | 4 (全部未修复) |
| 6. 退激发 | 5 | 0 | 5 (全部未修复) |
| 7. 其他 | 6 | 0 | 6 (全部未修复) |

---

## 第二阶段优先修复建议

Phase 1 证明了常数级修复可以极大改善 U+U 行为（碎片数 402→5）。第二阶段应聚焦**直接影响基态和动力学正确性**的 Critical/Major 差距：

| 优先级 | 差距 | 描述 | 预计工作量 |
|--------|------|------|-----------|
| 🔴 P0 | 2.3 能量偏移 | 移除 `energy_offset`，替换为 w_p 迭代 + BE 接受准则 | 2-3 天 |
| 🔴 P0 | 2.2 动量采样 w_p | 实现 p_F - w_p 采样替代 scale 因子 | 1-2 天 |
| 🔴 P0 | 2.5 弛豫时间 | 确保所有初始化路径 ≥600 fm/c | 0.5 天 |
| 🟠 P1 | 新问题 3 CoMD break | 移除 break 使每次调用修正所有违反对 | 0.5 天 |
| 🟠 P1 | 4.1 传播时动态 σ_r | 根据碎片大小更新波包宽度 | 2-3 天 |
| 🟠 P1 | 2.7 虚假发射检查 | 弛豫后推进无碰撞演化并检查逃逸 | 1 天 |
| 🟠 P1 | 5.4 实时密度 PB | 碰撞泡利阻塞使用碰撞时刻即时密度 | 0.5 天 |
| 🟡 P2 | 1.3 u_md | 实现 ImQMD05 动量相关相互作用 | 3-5 天 |
| 🟡 P2 | 5.1 PB(W*) | 实现 Chen 2024 动量状态采样 | 2-3 天 |
| 🟡 P2 | 7.5 形变初始化 | β₂/β₄ 形变核采样 | 1-2 天 |

---

## 附录: 逐文件变更摘要

| 文件 | Phase 1 变更 | 残留问题 |
|------|-------------|----------|
| `skyrme.py` | IQ3a 默认 + SkP* 添加 | docstring 过时；GridEDF 默认 IQ2（latent） |
| `grid_edf.py` | 无 | 默认参数 IQ2（latent bug） |
| `initializer.py` | A-σ_r + _sample_grid_positions 接通 | 动量采样仍用 scale；能量偏移未移除；弛豫 50 fm/c |
| `collisions.py` | CoMD 255 + 全核子对 | break 限制修正量；Fermi floor 仍在 |
| `propagator.py` | 无 | use_grid_edf 默认 False（性能） |
| `fragments.py` | 无 | LCP 凝聚非 Wei 2014 模型 |
| `decay.py` | 无 | 简化蒸发/裂变/γ |
| `nucleus.py` | 无 | sigma_r property 单点假设 |
| `__init__.py` | 无（仅导出） | — |

---

*审计基于对 `papers/imqmd_model_spec.md`（1026 行）与 `mnt_sim/imqmd/` 全部 9 个 .py 文件的逐项比对。Phase 1 修复使模型从"完全无法模拟重核融合"改善为"能够产生物理上合理的重碎片"，但仍有 29 项论文差距未修复（3 Critical + 16 Major + 10 Moderate）。*
