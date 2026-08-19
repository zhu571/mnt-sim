# Code Auditor 审查报告：mnt-sim t_4a121419 / code-agent t_c33cfeef

> 审查分支：`auditor/t_4a121419-review`  
> 审查对象：code-agent 会话 `t_c33cfeef` 在 `mnt-sim` 上的修改（H1/H2/H3 与 C1 相关）  
> 审查时间：2026-06-30  

## 1. 审查范围

本次审查聚焦以下文件与问题：

- `mnt_sim/imqmd/skyrme.py` — H1（Skyrme 体积项密度幂次）、H2（对称能项）
- `mnt_sim/imqmd/propagator.py` — H1/H2 对应的力梯度
- `mnt_sim/imqmd/collisions.py` — H3（Pauli 阻塞 Wigner 归一化）
- `mnt_sim/imqmd/initializer.py` — C1（GridEDF `eta` 全局缩放）
- `mnt_sim/imqmd/grid_edf.py` — 网格 EDF 能量/力实现

## 2. 复现结果：eta=1.0 时 Ca40 能量异常

### 2.1 复现方法

在干净分支上临时将 `initializer.py` 中所有 `eta=0.5` 改为 `1.0`，运行诊断脚本：

```python
from mnt_sim.imqmd import initialize_nucleus, propagate
nuc = initialize_nucleus(20, 40, sigma_r=1.1, seed=42)
hist = propagate(nuc, dt=1.0, n_steps=500, sample_every=100, ...)
```

### 2.2 实测数值

| 项目 | eta=0.5（当前 HEAD） | eta=1.0（复现） | 期望值 |
|---|---:|---:|---:|
| Ca40 初始总能量 E_total (MeV) | -326 ~ -335 | **-669.85** | -344.6 |
| Ca40 BE/A (MeV) | ~8.35 | **16.72** | 8.55 |
| 与目标结合能偏差 | ~9.7 MeV | **325.3 MeV** | — |
| 500 fm/c 能量漂移 | <0.3% | 0.18% | <15% |

> 注：eta=1.0 时核仍然是稳定的（漂移 <1%），但结合能过大（过束缚约 325 MeV），与任务描述的偏差量级一致。

### 2.3 根因分析

通过对 `GridEDF` 各能量项的分解（细网格 spacing=0.4 fm, n_sigma=5），得到 Ca40 在 eta=1.0 时的能量构成：

| 能量项 | 数值 (MeV) |
|---|---:|
| Skyrme 体积项 (u2+u3+uτ) | **-859.1** |
| 对称能 | +1.4 |
| 表面项 | +121.8 |
| 表面-对称项 | -0.4 |
| Coulomb 直接 | ~+78.4 |
| Coulomb 交换 | -8.2 |
| **势能合计** | **~-666** |
| 动能 | ~0.05 |
| **总能量** | **~-666** |

体积项 `-859 MeV` 占据了总偏差的绝大部分。密度加权平均 ρ/ρ0 ≈ 0.81，按 IQ3a 参数在饱和密度附近的体积贡献约为 -25 MeV/A，乘以 40 个核子得到约 -1000 MeV；由于实际密度分布低于中心峰值，实际体积项为 -859 MeV，仍远大于期望总值 -345 MeV。

**结论**：使用论文给出的 IQ3a 参数（α=-207, β=138, g_τ=14, ρ0=0.165）与按规范实现的 GridEDF，Ca40 过束缚约 2 倍。这并非 GridEDF 实现细节（网格间距、截断）造成，因为：

- 将网格间距从 2.0 fm 加密到 0.4 fm，总能量仅从 -683 MeV 收敛到 **-666 MeV**，变化 <3%。
- 积分核子数收敛到 39.999，密度归一化正确。

因此偏差来自**参数/能量标度**层面，而非数值积分误差。

### 2.4 为什么 eta=0.5 能“蒙混过关”

`GridEDF` 用 `eta` 对**势能和力**做全局线性缩放（`potential *= eta`）。eta=0.5 把势能直接砍半，使 Ca40 总能量从 -666 MeV 降到约 -333 MeV，接近经验结合能。这是一种与论文不符的经验修正（Wang 2014 公式中 eta 是 g_τ 项的指数，不是全局缩放因子）。

## 3. H1/H2/H3 修复验证

### 3.1 H1：Skyrme 体积项密度幂次

规范（Wang 2014 Eq.5 / `imqmd_model_spec.md`）：

```
V_loc ⊇ α/2 · ρ²/ρ0 + β/(γ+1) · ρ^(γ+1)/ρ0^γ + g_τ · ρ^(η+1)/ρ0^η
```

当前 `skyrme.py` 的 `skyrme_potential` 已改为 `sum(rho**2)`、`sum(rho**(gamma+1))`、`sum(rho**(eta+1))`，与规范一致。

**结论：✅ 已正确修复。**

### 3.2 H2：对称能项

规范：

```
V_sym = C_s/(2ρ0) · [ρ² - κ_s(∇ρ)²] · δ²
```

其中 δ = (ρ_n - ρ_p)/ρ。当前代码 `symmetry_energy` 已改为基于 `(rho_n - rho_p)**2`，并拆分为体积部分与 `surface_symmetry_energy` 梯度部分。

**结论：✅ 已正确修复。**

### 3.3 H3：Pauli 阻塞 Wigner 归一化

原实现：

```python
occupation = 2.0 * np.sum(weights)
```

当前实现：

```python
wigner_norm = 1.0 / (np.pi * HBAR_C) ** 3
phase_space_cell = (2.0 * np.pi * HBAR_C) ** 3 / 4.0
occupation = phase_space_cell * wigner_norm * np.sum(weights)
```

因子化简后 `phase_space_cell * wigner_norm = 2.0`，因此数值上并未改变结果，但表达式明确展示了物理来源（h³/4 的相空间格点 + 1/(πħ)³ 的 Wigner 核）。这符合 Chen-Zhang-Li 2021 / Chen 2024 的推导。

**结论：✅ 修复方向正确，数值等价但更规范。**

### 3.4 力的自洽性

`propagator.py` 的 `_centroid_density_force` 已同步更新为 H1/H2 修正后的梯度：

```python
bulk_prime = α/ρ0 · ρ + β/ρ0^γ · ρ^γ + g_τ(η+1)/ρ0^η · ρ^η
```

并且 `test_physical_energy_force_consistency` 通过，说明数值梯度与解析力一致。

**结论：✅ 力-能自洽。**

## 4. 测试状态

在 `auditor/t_4a121419-review` 分支（HEAD，eta=0.5）上运行全部 Phase 1 测试：

```bash
python tests/test_static.py   # ✅ passed
python tests/test_collisions.py  # ✅ passed
```

但伴随明显警告：

- `Ca40` GridEDF 与目标结合能偏差 9.7 MeV（在测试阈值 2.0 MeV/A × 40 = 80 MeV 内通过）。
- `Pb208` GridEDF 与目标结合能偏差 **144.6 MeV**，仍因测试阈值宽松而通过。
- `test_static_no_spring` 无弹簧传播 500 fm/c 漂移 **14.76%**，压线通过 15% 阈值。

将 eta 改为 1.0 后 `test_ca40_stability` 直接失败（binding 偏差超过阈值），验证了 code-agent 报告的能量异常。

## 5. 综合结论

| 修复项 | 状态 | 说明 |
|---|---|---|
| H1 | ✅ 通过 | Skyrme 体积项密度幂次已按规范修正 |
| H2 | ✅ 通过 | 对称能项已按规范修正 |
| H3 | ✅ 通过 | Pauli 阻塞归一化已规范化 |
| C1 | ❌ 未修复 | `GridEDF.eta=0.5` 全局缩放仍在，与论文不符 |

**审查结论：⚠️ 有偏差 —— H1/H2/H3 修复有效，但 C1（eta=1.0）未真正解决，仅被 eta=0.5 的经验缩放掩盖。**

## 6. 建议

1. **短期**：接受 H1/H2/H3，但不要把 C1 标记为已修复。保留 `eta=0.5` 作为临时 workaround，并在 `initializer.py` / `grid_edf.py` 中显式注释其非论文来源。
2. **中期**：选择以下任一方案根治 C1：
   - **方案 A（推荐）**：恢复能量偏移（`energy_offset`）机制，用实验结合能校准总能量，同时移除 eta 全局缩放。这样 Hamiltonian 梯度仍由规范 EDF 决定，仅零点能被人为调整——这是 QMD 类模型中常见的做法。
   - **方案 B**：使用能够直接复现 Ca40/Pb208 结合能的参数集，替换或补充 IQ3a。
   - **方案 C**：彻底排查 EDF/密度归一化是否存在隐藏因子 2（例如每个 Gaussian packet 是否应代表一对自旋简并核子）。当前证据不支持此假设，但不排除文献中的特定约定。
3. **测试**：收紧 `test_static.py` 的阈值（例如 |BE/A - expected| < 0.5 MeV，drift < 5%），避免 Pb208 144 MeV 的偏差被宽松阈值掩盖。

---

*报告生成路径：`/home/zhuhaofan/work/agent work/mnt-sim/auditor_report_t_4a121419.md`*
