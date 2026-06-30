# ImQMD 全面审计报告

> **审计日期**: 2026-06-27  
> **审计范围**: `mnt_sim/imqmd/` 全部 7 个核心模块 + 3 个测试文件 + `verify_transfer.py`  
> **模型版本**: Phase 1-4, IQ2 默认参数集  
> **审计目标**: 代码质量、物理正确性、稳定性（为什么 no-spring 在 nucimp 上失败）、性能热点、架构评估、优先级建议

---

## 1. 执行摘要

ImQMD 模型代码结构良好，测试覆盖率合理（16 项测试全部通过），但在 **物理正确性** 上存在一个关键缺陷：**表面力（surface pair force）被错误地与人工静态弹簧（static reference spring）捆绑在一起**。当 `use_static_reference=False`（反应模式）时，哈密顿量中仍包含表面能量项，但运动方程中缺少对应的表面力。这导致：

- **笔记本上**：轻核（⁴⁰Ca）短期（500 fm/c）测试勉强通过（8% 能量漂移容差很宽松）
- **nucimp 上**：重核（²³⁸U+²³⁸U, 476 核子）长期（2500 步）演算中，缺乏表面张力 + 碰撞加热 → 核完全解体

其他主要发现包括：激发能 `E*` 使用了不兼容的能量标度、Fermi 约束在 no-spring 模式下频率过高（每 5 fm/c）、O(N²) 力计算是主要性能瓶颈。

---

## 2. 逐个模块审计

### 2.1 `skyrme.py` — Skyrme EDF 参数和能量密度（148 行）

**代码质量**: ★★★★☆  
**物理正确性**: ★★★★☆（surface_pair_energy 的物理角色需澄清）

| 方面 | 评估 |
|------|------|
| 结构 | `SkyrmeParameters` 是不可变 dataclass，5 个参数集清晰枚举，设计良好 |
| 数值稳定性 | `rho = np.clip(rho, 1e-12, None)` 在所有地方防止除零 |
| 公式映射 | `skyrme_potential()` 正确实现了 Wang 2014 Eq.(5) 的 α, β, γ, g_τ, η 项 |
| 对称能 | `symmetry_energy()` 使用 centroid 密度采样，符合 Phase 1 简化范围 |

**问题**:

1. **`surface_pair_energy()` 的物理角色模糊** (skyrme.py:106-125):
   - 文档注释说"surface cohesion applies between all nearby nucleons regardless of origin"
   - 但实际上这个函数实现的是 Skyrme EDF 中 `g_sur/(2*ρ₀) (∇ρ)²` 表面梯度项的粒子形式（等价于有限程高斯对相互作用）
   - `gsur/7.0` 归一化因子需要文献验证——在 IQ2 参数集 `gsur=7.0` 时该因子退化为 1.0，巧合但合理
   - **关键问题**: `group_ids` 参数标注为"backward compatibility"，但注释又说对所有核子一视同仁——这自相矛盾

2. **`static_k` 污染了 EDF 类** (skyrme.py:54-56):
   - 人工谐波弹簧常数存储在 `SkyrmeEDF` 对象上，但这是一个数值稳定器，不是物理 EDF 参数
   - 应移到 `ImQMDNucleus` 或 propagator 上下文

3. **`static_mean_field_force` 和 `static_mean_field_energy` 命名误导** (skyrme.py:132-145):
   - 这不是"平均场力"（mean field force），而是纯数值人工弹簧
   - 应重命名为 `harmonic_stabilizer_force/energy` 或类似名称

### 2.2 `nucleus.py` — 核容器和密度评估（151 行）

**代码质量**: ★★★★☆  
**物理正确性**: ★★★☆☆（能量分量存在不一致）

| 方面 | 评估 |
|------|------|
| 结构 | `GaussianPacket` + `ImQMDNucleus` 分离清晰 |
| 密度计算 | `density()` 正确实现 Wang 2014 Eq.(1)，支持单点和多点评估 |
| 属性设计 | `positions`/`momenta` 作为 getter/setter 正确同步到 `GaussianPacket` |
| `copy()` | 深拷贝实现正确 |

**问题**:

1. **`energy_components()` 无条件包含 `surface_pair_energy`** (nucleus.py:106-109):
   ```python
   surface = self.edf.surface_pair_energy(
       self.positions, self.sigma_r,
       getattr(self, "reference_group_ids", None),
   )
   ```
   - 无论 `use_static_reference` 标志如何，表面能量总是计算
   - 但在 `propagator.py` 中，`_surface_pair_force()` 只在 `use_static_reference=True` 时执行
   - **这是哈密顿量与力的不一致**: 能量包含表面项，但力不包含 → 能量不守恒
   - 这是 no-spring 模式下核解体的**直接原因**（详见第 3 节）

2. **`reference_group_ids` 副作用**:
   - 当存在 `reference_group_ids` 时，`surface_pair_energy` 只考虑同组内的对
   - 在反应中这意味着表面力只在弹核/靶核内部作用，不跨颈部 → 阻止融合和转移
   - 这是 P1 的第二个机制（除了静态弹簧之外）

3. **`rms_radius()` 包含 packet width 修正** (nucleus.py:142):  
   `+ 3.0 * self.sigma_r**2` — 对于 σ_r = 1.1 fm 这是 3.63 fm² 的修正。这是标准 QMD 做法，但需要确认与实验比较时的一致性。

### 2.3 `initializer.py` — 基态初始化（204 行）

**代码质量**: ★★★★☆  
**物理正确性**: ★★★☆☆（能量偏移方法值得商榷）

| 方面 | 评估 |
|------|------|
| 硬球采样 | `_sample_hard_sphere()` 正确，有自适应最小距离回退 |
| Fermi 动量 | `_fermi_momenta()` 基于局部密度计算 p_F，采样在缩小的 Fermi 球内 |
| 相空间约束 | `_phase_space_minimum()` 与 CoMD 标准一致 |

**问题**:

1. **能量标度偏移方法** (initializer.py:132-135):
   ```python
   target_total = -empirical_binding_per_nucleon(Z, A) * A
   raw_components = nucleus.energy_components()
   nucleus.energy_offset = target_total - raw_components["total"]
   ```
   - 将 raw QMD 哈密顿量强行偏移到经验结合能
   - 这不是物理校准——它是在 QMD 能量尺度上添加一个任意常数
   - 后果1: 当在 `fragments.py` 中计算 `E*` 时，片段子核不带此偏移 → 能量尺度不匹配
   - 后果2: 两个不同初始化的核（不同 seed）有不同偏移，合并时总偏移无定义

2. **`initialize_nucleus()` 的松弛阶段使用 `use_static_reference=False`** (initializer.py:124):
   - 与方案描述矛盾——文档说松弛使用静态弹簧
   - 实际代码在 `_rk4_step(nucleus, 1.0, remove_cm_drift=True, use_static_reference=False)` 中跳过了弹簧和表面力
   - 依赖手动动量阻尼 (`*0.93`) 和 Fermi 约束来稳定核

3. **`initialize_and_relax()` 的容差** (initializer.py:190):
   - `drift <= 0.08` 即 8% 能量漂移被认为可接受
   - 对于反应模拟（需要追踪几百 MeV 的能量差），8% 的标度误差可能淹没物理信号

4. **相空间阈值** (initializer.py:100):
   - `phase_space_threshold = 170.0`（代码内）vs 论文中的 255 fm·MeV/c
   - 这是有意的放宽（为了初始化成功率），但应该在文档中明确标注为"开发简化"

### 2.4 `propagator.py` — 时间传播（229 行）

**代码质量**: ★★★☆☆  
**物理正确性**: ★★☆☆☆（**力不匹配哈密顿量，核心 Bug**）

| 方面 | 评估 |
|------|------|
| RK4 实现 | 标准四阶龙格-库塔，实现正确 |
| 库仑力 | 解析梯度正确，与 `coulomb_energy` 一致 |
| CM 漂移修正 | 每步移除，适合单核稳定性测试 |

**严重问题**:

1. **`_derivatives()` 中 `use_static_reference` 标志捆绑了两个独立概念** (propagator.py:119-131):
   ```python
   def _derivatives(nucleus, positions, momenta, use_static_reference=True):
       drdt = momenta / M_N
       dpdt = (_static_reference_force(nucleus, positions) 
               if use_static_reference else np.zeros_like(positions))
       dpdt += _centroid_density_force(nucleus, positions)
       if use_static_reference:                          # ← BUG: 表面力也受此门控
           dpdt += _surface_pair_force(nucleus, positions)
       dpdt += _coulomb_forces(positions, nucleus.is_proton, nucleus.sigma_r)
       return drdt, dpdt
   ```
   - `_static_reference_force` 是人工弹簧（应只在稳定性测试中使用）
   - `_surface_pair_force` 是物理 Skyrme 表面梯度项的力（应始终存在）
   - **这两者被错误地捆绑在同一个 `use_static_reference` 标志下**
   - 当 `use_static_reference=False` 时，表面力也被移除，但 `energy_components()` 仍然包含表面能量

2. **力的能量不一致**:
   | 能量项 | energy_components() | _derivatives() with use_static_reference=False |
   |--------|---------------------|-----------------------------------------------|
   | 体 Skyrme | ✓ | ✓ (via _centroid_density_force) |
   | 对称能 | ✓ | ✓ (via _centroid_density_force) |
   | 库仑 | ✓ | ✓ (via _coulomb_forces) |
   | 表面对能量 | ✓ | ✗ (_surface_pair_force 被跳过) |
   | 静态弹簧 | ✓ | ✗ (正确跳过) |

   这意味着 **no-spring 模式的哈密顿力学不是保守的**。已经缺失的力的总能量效应取决于系统大小和构型。

3. **`_surface_pair_force()` 使用 `reference_group_ids` 边界** (propagator.py:86-91):
   - 当 `group_ids` 存在时，权重矩阵被同组掩码截断
   - 这意味着在反应中，即使表面力被启用，它也只作用于弹核/靶核内部的核子对
   - 颈部核子之间没有表面张力 → 进一步阻止融合/转移

4. **Fermi 约束时间间隔** (propagator.py:207):
   ```python
   fermi_time = 20.0 if use_static_reference else 5.0
   ```
   - no-spring 模式下约束检查频率是 with-spring 的 4 倍
   - 逻辑合理（没有弹簧需要更频繁的相空间管制），但 5 fm/c 对于重核 O(N²) 约束检查非常昂贵

5. **`_centroid_density_force()` 的内存复杂度** (propagator.py:37-75):
   - 构建 `A × A × 3` 的 diff 张量 → 476² × 3 × 8 bytes ≈ 5.4 MB per call
   - 每步 4 次导数调用 → 每步约 22 MB 临时分配
   - 对称能梯度的分析推导看起来正确，但未经验证（没有针对数值梯度的单元测试）

### 2.5 `collisions.py` — NN 碰撞和 Pauli 阻塞（429 行）

**代码质量**: ★★★☆☆  
**物理正确性**: ★★★☆☆

| 方面 | 评估 |
|------|------|
| 截面模型 | free_nn_cross_section 使用了正确的 Cugnon 参数化和低能 cap |
| 介质修正 | in_medium_factor 使用 η=0.2 的标准形式 |
| 邻居列表 | 网格化空间搜索减少了候选对，有自适应重建 |
| Pauli 阻塞 | 标准 Wigner 核方法，包含 Fermi 球下限保护 |

**问题**:

1. **`fermi_constraint_check()` 的 O(N²) 循环** (collisions.py:380-414):
   - 对相同同位旋的所有对进行双重循环
   - 对 U+U (476 核子，~238 质子)，这需要检查约 28k 对
   - no-spring 模式每 5 fm/c 调用一次 → 2500 步需约 500 次完整扫描
   - 有 `break` 优化（找到一对违规即跳出内层），但外层仍扫描所有未使用的 i
   - **建议**: 使用空间格点限制到邻近对

2. **`attempt_nn_collision()` 中的质心动量归零** (collisions.py:374):
   ```python
   nucleus.momenta = momenta - np.mean(momenta, axis=0)
   ```
   - 每次碰撞调用后移除总质心动量
   - 对于反应事件（两个核相向运动），这相当于每碰撞间隔施加一个微小的反冲校正
   - 在 `remove_cm_drift=False` 的反应传播中，这是一个隐藏的动量修正

3. **Pauli 阻塞性能** (collisions.py:151-218):
   - `_pauli_blocking_probability_fast` 对每个尝试的对扫描所有相同同位旋的核子
   - 对 U+U 有 ~17000 次接受碰撞（每个需要 1 次 Pauli 检查），每次检查扫描 ~238 个核子
   - 总计约 4M 次距离计算——可优化为空间邻居

4. **`_occupation_wigner_fast` 的密度参数冗余** (collisions.py:151-179):
   - 接受 `rho_q` 参数用于 Fermi 动量计算，但同时调用 `nucleus.density(r_i)` 会更一致
   - 当前使用 centroid 密度，可能与实际局部密度有偏差

5. **模板重复**: `pauli_blocking_probability` 和 `pauli_blocking_probability_v2` 几乎相同（仅 `width_scale` 参数不同）。应合并为一个带参数的方法。

### 2.6 `fragments.py` — 碎片识别（303 行）

**代码质量**: ★★★★☆  
**物理正确性**: ★★★☆☆（E* 标度不匹配）

| 方面 | 评估 |
|------|------|
| MST 实现 | 正确使用 scipy connected_components（有回退实现） |
| 同位旋 MST | 实现正确，阈值来自文献 |
| 自适应 MST | A 标度截断设计合理 |
| 聚结 | 贪婪算法，模板优先级正确（d,t,³He,⁴He） |

**问题**:

1. **`compute_fragment_excitation()` 能量标度不匹配** (fragments.py:246-292):
   - 片段子核以 `energy_offset=0.0` 构建 → raw QMD 能量（未校准）
   - 减去 `binding_energy(Z,A)` → 经验/质量表结合能
   - **这两个值在不同的物理标度上**
   - Parent nucleus 的 `energy_offset`（~几百 MeV）不传给子核
   - 结果: 冷碎片 E* 被 clamp 到 0（掩盖误差），热碎片 E* 可达数百 MeV（夸大）

2. **静态变形能量折入 E\*** (fragments.py:271-286):
   - 当 parent 仍有 `reference_positions` 和 `static_k > 0` 时，添加人工变形能
   - 此项不是物理激发能的一部分

3. **`reaction_fragments()` 的 `p_cut=None` 默认值** (fragments.py:161-165):
   - 反应模式默认不使用动量截断
   - 可能导致两个已分离的碎片因动量接近而被 MST 连接
   - `verify_transfer.py` 和测试正确使用了 `reaction_fragments(system)` 的默认值
   - 但生产级分析可能需要动量截断（~250 MeV/c）

4. **`adaptive_mst` A-标度公式** (fragments.py:151):  
   `cutoffs = base_r_cut * cbrt(pair_A) / cbrt(16.0)` — 以 A=16 为参考，物理上合理但除数为硬编码。

### 2.7 `decay.py` — 统计退激（190 行）

**代码质量**: ★★★★☆  
**物理正确性**: ★★★★☆（简化但合理）

| 方面 | 评估 |
|------|------|
| 质量表 | 正确加载 Mexcess95.dat，有液滴回退 |
| 结合能 | `binding_energy()` 正确使用质量过剩公式 |
| 中子蒸发 | Weisskopf 统计权重实现正确 |
| 裂变竞争 | Bohr-Wheeler 简化形式，重核增强合理 |

**问题**:

1. **`weisskopf_evaporation()` 的窗口因子** (decay.py:171):
   ```python
   window = (E_star / max(cumulative_cost, 1.0)) ** n
   ```
   - 此因子使中等激发能（20 MeV）的情况不受未蒸发通道支配
   - 物理动机合理但需要文献引用
   - 对小 E* 可能导致蒸发通道低概率——但实现意图正确

2. **蒸发通道不含质子/α/γ**: 只有中子蒸发。对接近滴线的富质子碎片会有偏差。

3. **`fission_competition()` Z≥104 增强因子** (decay.py:132): `width_f *= 20.0` 是为超重核设计的简单增强，但数值缺乏论文引用。

---

## 3. 稳定性分析：为什么 no-spring 在 nucimp 上失败

### 3.1 症状

| 环境 | 核系统 | 行为 |
|------|--------|------|
| 笔记本 | ⁴⁰Ca, 500 fm/c, no-spring, 无碰撞 | 通过 (8% 能量漂移) |
| 笔记本 | ⁴⁰Ca+⁴⁰Ca, 150 fm/c, no-spring, 无碰撞 | 通过 (两碎片分离) |
| nucimp | ⁸⁶Kr+⁶⁴Ni, 600 fm/c, no-spring, 有碰撞 | 熔合或弹性（无转移） |
| nucimp | ²³⁸U+²³⁸U, 2500 fm/c, no-spring, 有碰撞 | **核完全解体** |

### 3.2 根因分析

no-spring 模式（`use_static_reference=False`）有**两个独立的问题**：

#### 根因 1: 表面力缺失（哈密顿-力不一致）

```python
# propagator.py:119-131 — BUG 位置
dpdt = _static_reference_force(nucleus, positions) if use_static_reference else 0
dpdt += _centroid_density_force(nucleus, positions)   # 只有体和对称能
if use_static_reference:                               # ← 表面力也受此门控!
    dpdt += _surface_pair_force(nucleus, positions)    # ← 缺失
dpdt += _coulomb_forces(...)                           # 库仑力
```

当 `use_static_reference=False`:
- **体 Skyrme 力**存在（提供体引力）
- **表面张力**缺失（没有有限程表面对力）
- **库仑排斥力**存在

净效应：**短程核引力不足以对抗长程库仑排斥和热运动**。对于小核（⁴⁰Ca），体/面比大，体引力占主导，短期稳定；对于大核（²³⁸U），库仑排斥与 Z² 成正比增长但表面张力不增长，**核逐渐膨胀直至解体**。

#### 根因 2: `reference_group_ids` 破坏跨颈部相互作用

`_surface_pair_force()` 和 `_static_reference_force()` 在检测到 `reference_group_ids` 时对组内执行弹簧/表面力。在反应中（`verify_transfer.py:108`: `system.reference_group_ids = ...`），弹核=组0，靶核=组1。这意味着：

- 即使启用表面力，颈部核子之间也没有吸引配对（属于不同组）
- 静态弹簧也将每个核子拴在原始组内

这是 P1"无转移"的**双重机制**：弹簧拴住 + 表面力不跨组。

#### 根因 3: 阻尼在传播中不存在

`initialize_and_relax()` 使用 `momenta * 0.93` 阻尼来稳定 no-spring 松弛。但 `propagate()` 没有此阻尼。在 nucimp 上运行的反应（2500 步）纯粹依赖保守（或不守恒）动力学加上碰撞加热，**没有任何冷却机制**。结果是热失控。

### 3.3 笔记本 vs nucimp 差异解释

| 因素 | 笔记本测试 | nucimp 运行 |
|------|-----------|------------|
| 核子数 A | 40-56 | 150-476 |
| 时间步数 | 500 | 600-2500 |
| 碰撞 | 否（或仅测试中） | 是（~5000-17000 次接受） |
| 库仑力 | ~190 对 | ~113k 对（U+U） |
| 体/面比 | 高（表面张力相对不重要） | 低（表面张力关键） |
| Fermi 约束频率 | N/A | 每 5 fm/c（4× 更频繁） |

笔记本测试在 500 fm/c 内能通过，是因为：小核 + 无碰撞 + 短期。nucimp 失败是因为：大核 + 碰撞加热 + 长期 + 缺失表面张力 = 解体。

---

## 4. 性能热点

### 4.1 热点分布

对于 ²³⁸U+²³⁸U（A=476），2500 步，单核：

| 热点 | 位置 | 每步调用 | 每次复杂度 | 估计占比 |
|------|------|---------|-----------|---------|
| **密度力** | `_centroid_density_force` | 4 (RK4) | O(A²) = ~227k 对 | **~45%** |
| **库仑力** | `_coulomb_forces` | 4 (RK4) | O(Z²) = ~113k 对 | **~20%** |
| **NN 碰撞** | `attempt_nn_collision` | 1 per 碰撞间隔 | O(pairs × A) | **~15%** |
| **Fermi 约束** | `fermi_constraint_check` | 1 per 5 fm/c | O(A²_samespecies) | **~10%** |
| **表面力** | `_surface_pair_force` | 4 (RK4, when enabled) | O(A²) | **~5%** |
| **能量诊断** | `centroid_densities` | 1 per sample | O(A²) | **~2%** |
| 其他 | — | — | — | ~3% |

### 4.2 可优化点

1. **密度力 + 库仑力缓存**: `_centroid_density_force` 和 `_coulomb_forces` 各自独立计算对距离。合并为单次对距离计算并可复用。
2. **密度力的空间截断**: 高斯核在 3σ_r ≈ 3.3 fm 处已衰减到 <1%。使用邻居列表限制计算范围可将 O(A²) → 有效 O(A)。
3. **Fermi 约束邻居化**: 当前对所有相同同位旋对进行全扫描。使用空间格点限制到接近对可将大核成本降低 10-100×。
4. **RK4 → 速度-Verlet**: 将力评估从每步 4 次减到 1-2 次。需验证能量守恒。
5. **Numba/Cython 编译**: 浮点密集型循环（密度力、库仑力、邻居列表）是 JIT 编译的理想候选。

### 4.3 邻居列表评估

`_cached_neighbor_pairs()` (collisions.py:275-311) 实现正确：
- 5 fm 格点大小
- 自适应重建（位移 > 2.5 fm 或每 10 步）
- 仅用于 NN 碰撞候选，**不用于力计算或密度计算**

扩展邻居列表到力计算可将最昂贵的 O(A²) 内核转化为接近 O(A)。

---

## 5. 架构评估

### 5.1 优势

| 方面 | 说明 |
|------|------|
| 模块化 | 7 个模块，每个有明确的物理职责 |
| 类型注释 | 全面使用 `from __future__ import annotations` 和类型提示 |
| 不可变参数 | `SkyrmeParameters` 使用 `frozen=True` dataclass |
| 测试驱动 | 3 个测试文件覆盖静态、碰撞、碎片阶段 |
| 文档 | `papers/imqmd_formulas.md` 详细记录公式来源 |
| 可配置性 | 5 个 Skyrme 参数集可通过名称切换 |

### 5.2 劣势

| 方面 | 问题 |
|------|------|
| **力/能量不一致** | `_derivatives()` 产生的力不对应 `energy_components()` 的能量（见 §2.4） |
| **混用物理与数值概念** | `static_k`、`use_static_reference`、`energy_offset`、`reference_group_ids` 跨越模块边界，缺少明确的 API 契约 |
| **Nucleus 上的隐式状态** | `_collision_rng`, `_neighbor_list`, `_active_collision_pair`, `collision_stats`, `reference_group_ids` 不是构造函数参数而是直接设置在 `ImQMDNucleus` 实例上 |
| **缺少能量守恒验证** | 没有自动化的数值梯度检查（对比解析力和有限差分） |
| **无中间表示** | 传播的 `history` 是 dict 列表，碎片是 `Fragment` 列表，两者之间缺少结构化的"事件状态"对象 |

### 5.3 建议的接口清理

```python
# 当前（有问题）
nucleus._collision_rng = np.random.default_rng(seed)
nucleus.reference_group_ids = np.array([0, 0, 1, 1])

# 建议
class ReactionState:
    collision_rng: np.random.Generator
    neighbor_cache: NeighborCache
    group_assignments: np.ndarray  # 0=projectile, 1=target
    collision_stats: dict[str, int]

# 传播
def propagate(nucleus, dt, n_steps, 
              use_physics_surface=True,      # 物理表面力
              use_numerical_stabilizer=False, # 人工弹簧
              reaction_state=None):
    ...
```

**关键修复**: 将 `use_static_reference` 分解为两个独立标志：
- `use_surface_term: bool = True` — 控制物理 Skyrme 表面力
- `use_static_stabilizer: bool = False` — 控制人工谐波弹簧（仅用于静态稳定性测试）

---

## 6. 测试审计

### 6.1 `tests/test_static.py`（6 项测试）

| 测试 | 评估 |
|------|------|
| `test_ca40_stability` | ✓ 合理：2000 fm/c, RMS/能量/半径约束 |
| `test_pb208_stability` | ✓ 合理：重核实测 |
| `test_single_nucleus_energy_conservation` | ⚠ 1% 容差——在有弹簧(use_static_reference=True)的情况下易通过 |
| `test_static_no_spring` | ⚠ 仅 500 fm/c, 8% 容差——太宽松，可能掩盖解体趋势 |
| `test_final_state_pauli` | ✓ 验证 Pauli 阻塞率 |
| `test_collision_rate_sanity` | ✓ 验证碰撞发生 |

**缺失的测试**:
- No-spring 长期稳定性（2000+ fm/c, 重核）
- 力的数值梯度一致性（解析 `dpdt` vs 有限差分）
- No-spring 模式的能量守恒（期望高漂移，需要量化）
- 多核系统 no-spring 传播（如 2×⁴⁰Ca 分离）

### 6.2 `tests/test_collisions.py`（4 项测试）
所有测试合理且通过。截面值、介质缩放、Pauli 阻塞率、Fermi 约束均有覆盖。

### 6.3 `tests/test_fragments.py`（7 项测试）
- MST 分离、聚结、蒸发测试均合理
- `test_mnt_event_has_transfer` 仅在 `RUN_SLOW_MNT_TESTS=1` 时运行——这是合理的慢测试门控

### 6.4 `examples/verify_transfer.py`
良好的端到端验证。`use_static_reference=False` 使用正确。断言合理（检查 Z/A 转移）。

---

## 7. 优先建议

### P0: 紧急（阻塞物理结果）

| # | 建议 | 文件 | 预计工作量 |
|---|------|------|-----------|
| **1** | **拆分 `use_static_reference` 为两个独立标志** | propagator.py | 1h |
| | - `use_surface_term: bool = True` 控制物理表面力 | | |
| | - `use_static_stabilizer: bool = False` 控制人工弹簧 | | |
| | - 默认: surface=True, stabilizer=False（物理模式） | | |
| **2** | **使 `_surface_pair_force` 在物理模式下始终启用** | propagator.py:128-129 | 紧跟 #1 |
| | - 当 `use_surface_term=True` 时无条件启用 | | |
| | - 当存在 `reference_group_ids` 时不按组掩码截断 | | |
| **3** | **移除 `energy_components()` 中 `surface_pair_energy` 的 `group_ids` 参数** | nucleus.py:106-109, skyrme.py:106-125 | 紧跟 #2 |
| | - 表面能量应始终在所有对之间计算 | | |

### P1: 高优先（提高物理准确性）

| # | 建议 | 文件 | 预计工作量 |
|---|------|------|-----------|
| **4** | **新增力/能量一致性测试** | tests/ | 2h |
| | - 对随机构型比较解析力与有限差分 | | |
| | - 验证 no-spring 传播的瞬时 `dE/dt ≈ 0` | | |
| **5** | **校准 `compute_fragment_excitation` 能量标度** | fragments.py | 4h |
| | - 使用相同基态的 calibrated QMD 能量作为地面参考 | | |
| | - 移除静态变形能添加 | | |
| | - 重正片段位置后再计算内能 | | |
| **6** | **重评估 `initialize_and_relax` 的 no-spring 松弛** | initializer.py | 3h |
| | - 使用 `use_surface_term=True` 重新测试 | | |
| | - 收紧能量漂移容差（8% → 3%） | | |
| | - 评估是否需要表面项进行初始松弛 | | |

### P2: 中优先（性能与可维护性）

| # | 建议 | 文件 | 预计工作量 |
|---|------|------|-----------|
| **7** | **将邻居列表扩展到力计算** | propagator.py | 4h |
| | - 在 `_centroid_density_force` 中使用空间截断 | | |
| | - 在 `_coulomb_forces` 中使用邻居列表 | | |
| **8** | **邻居化 `fermi_constraint_check`** | collisions.py | 2h |
| | - 使用空间格点限制对检查范围 | | |
| **9** | **清理隐式状态 API** | nucleus.py, collisions.py | 3h |
| | - 将 `_collision_rng`, `_neighbor_list`, `collision_stats` 移入显式状态对象 | | |
| | - 将 `static_k` 移出 `SkyrmeEDF` | | |
| **10** | **添加数值梯度单元测试** | tests/ | 2h |
| | - 验证所有力分量与能量的解析梯度匹配 | | |
| | - 为未来重构提供安全网 | | |

### P3: 低优先（增强功能）

| # | 建议 | 预计工作量 |
|---|------|-----------|
| **11** | 为 `weisskopf_evaporation` 添加质子/α 蒸发通道 | 4h |
| **12** | 实现反应终止条件（重新分离检测）替代固定步数 | 6h |
| **13** | 速度-Verlet 积分器作为 RK4 的替代方案 | 4h |
| **14** | 构建 Zhao 2016 产额比较管线 | 8h |
| **15** | 添加 `_centroid_density_force` 中对称能梯度的独立推导文件 | 2h |

---

## 8. 代码统计

| 文件 | 行数 | 函数 | 类 | 注释/总行比 |
|------|------|------|-----|------------|
| skyrme.py | 148 | 8 | 2 | 30% |
| nucleus.py | 151 | 10 | 2 | 25% |
| initializer.py | 204 | 6 | 0 | 20% |
| propagator.py | 229 | 8 | 0 | 25% |
| collisions.py | 429 | 18 | 0 | 20% |
| fragments.py | 303 | 13 | 1 | 20% |
| decay.py | 190 | 9 | 0 | 25% |
| **合计** | **1654** | **72** | **5** | **~23%** |

测试代码: ~390 行（3 个测试文件），覆盖 16 项功能/集成测试。

---

## 9. 结论

ImQMD 代码库架构良好，公式映射仔细，但存在一个影响所有反应结果的**根本性物理缺陷**：

> **表面力被错误地捆绑在人工弹簧门控下。no-spring 模式缺少 Skyrme EDF 的表面梯度贡献，导致哈密顿量与运动方程不一致。这是 nucimp 上核解体的直接原因，也是转移/融合被抑制的促成因素。**

修复方案清晰且风险有限：拆分 `use_static_reference` 为 `use_surface_term` 和 `use_static_stabilizer`，默认启用物理表面力。这可以在不影响现有静态测试（继续使用 `stabilizer=True`）的情况下修复反应动力学。

修复后，应优先解决激发能标度问题（P1-#5）和力/能量一致性验证（P1-#4），然后进行性能优化（P2-#7,#8）。

---

*审计人: Hermes Agent (Nous Research), 基于 deepseek-v4-pro 模型*  
*参考文献: Wang et al. PRC 65, 064608 (2002); Wang et al. PRC 89, 064601 (2014); Chen et al. PRC 109, L021604 (2024); Zhao et al. PRC 94, 024601 (2016); 及 papers/imqmd_formulas.md 中引用的其他文献*
