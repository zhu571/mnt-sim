# 审计报告 t_3a5c03e2 — GridEDF 碰撞动力学: U+U 碎裂根因分析

**审计日期**: 2026-07-07  
**审计范围**: `grid_edf.py`, `propagator.py`, `fragments.py`, `initializer.py`, `reaction.py`  
**关键参数集**: IQ3A (α=-207, β=138, γ=7/6, g_tau=14, η=5/3, ρ0=0.165)

---

## 审计结论: ❌ 需修改

发现 1 项严重问题（nuclear_scale 扭曲高密度 EOS）和 2 项偏差，综合判断碰撞动力学中的 GridEDF 力计算存在根本性物理缺陷。

---

## 🔴 严重: nuclear_scale 全局缩放力导致高密度 EOS 软化

**位置**: `grid_edf.py:forces_analytical` L198-213, `initializer.py:_fit_grid_nuclear_scale` L143-152

**问题描述**:

`forces_analytical` 对所有核力项（体项、对称项、表面项）施加统一的 `self.nuclear_scale` 因子：

```python
# grid_edf.py L203-213
d_e_drho_n = self.nuclear_scale * bulk_prime
d_e_drho_p = self.nuclear_scale * bulk_prime
d_e_drho_n += self.nuclear_scale * (p.c_sym / p.rho0 * asym)
...
d_e_drho_n += self.nuclear_scale * (-(p.gsur / p.rho0) * lap)
```

`nuclear_scale` 通过 `_fit_grid_nuclear_scale` 在 **正常核密度 (ρ≈ρ0)** 下校准：
```python
scale = (target_total - kinetic - coulomb) / nuclear
return np.clip(scale, 0.35, 1.25)
```

这个单一标量在 **碰撞重叠区高密度 (ρ≫ρ0)** 下被不加区分地施加到所有密度依赖项上。

**定量证据** (IQ3A 参数):

| ρ/ρ0 | bulk_prime (ns=1.0) | ns=0.95 | ns=0.90 |
|-------|---------------------|---------|---------|
| 1.0   | -31.7 (吸引)        | -30.1   | -28.5   |
| 1.5   | -15.6 (吸引)        | -14.9   | -14.1   |
| 1.8   | ≈0   (**零交叉**)    | —       | —       |
| 2.0   | +14.3 (排斥)        | +13.6   | +12.9   |
| 2.5   | +56.3 (排斥)        | +53.5   | +50.7   |
| 3.0   | +109.2 (排斥)       | +103.7  | +98.2   |

- 零交叉点: ρ/ρ0 = 1.80  
- nuclear_scale < 1 会均匀减弱所有密度下的力，但高密度排斥力的绝对减弱量远大于低密度吸引力  
- 对于重核 (如 U-238, scale 典型 ~0.90-0.95): 排斥力减弱 5-10%，允许更深穿透 → 更高压缩 → 爆炸性碎裂

**物理后果**:

U+U 碰撞中 nuclear_scale ≈ 0.92-0.95:
1. 重叠区排斥力比应有值低 5-8%
2. 原子核比物理上正确的情况穿透更深（ρ 可达 3ρ0 而非 2.5ρ0）
3. 额外储存的压缩能 → 爆炸性解压 → 50+ 碎片（观测到 52-58）
4. 这正是静态测试通过（单核 ρ≈ρ0 下校准准确）但碰撞动力学失败的原因

**根因**: 单一全局标量 scale 无法捕捉能量密度泛函的密度依赖性修正。不同 Skyrme 项 (α·ρ², β·ρ^(γ+1), g_tau·ρ^(η+1)) 有不同的密度幂律，用一个 scale 因子统一缩放扭曲了状态方程的形状。

---

## 🟡 偏差: 初始化松弛与碰撞动力学使用不同的 nuclear_scale

**位置**: `initializer.py:initialize_nucleus` L351, `initializer.py:initialize_and_relax` L491, `propagator.py:propagate` L302-306

**问题描述**:

初始化松弛阶段使用的 GridEDF 具有 `nuclear_scale=1.0`（默认值）:
```python
# initializer.py L351
relaxation_grid = GridEDF(nucleus.edf.parameters, nucleus.sigma_r, grid_spacing=1.0)
# → nuclear_scale=1.0

# initializer.py L491
relax_grid = GridEDF(nucleus.edf.parameters, nucleus.sigma_r, grid_spacing=1.0)
# → nuclear_scale=1.0
```

但碰撞动力学使用 `nuclear_scale=grid_nuclear_scale`（通常 ≠ 1.0）:
```python
# propagator.py L301-306
grid_edf = GridEDF(nucleus.edf.parameters, nucleus.sigma_r,
                   nuclear_scale=float(getattr(nucleus, "grid_nuclear_scale", 1.0)))
```

`quick_no_spring_drift` 检查使用动力学 scale 验证稳定性，但容忍度仅为 8% 能量漂移（500 fm/c）。这不是零—意味着原子核在碰撞前就已微弱不稳定。

---

## 🟡 偏差: 合并系统的 energy_offset 未被继承

**位置**: `reaction.py:make_collision_event` L155-161

**问题描述**:

```python
system = ImQMDNucleus(
    projectile_z + target_z,
    (projectile_a - projectile_z) + (target_a - target_z),
    packets, edf=edf,
    reference_positions=...,
    # energy_offset 未传递，默认 = 0.0
)
```

每个核有自己的 `energy_offset`（补偿 GridEDF 能量与经验结合能的差异），但合并系统的 `energy_offset=0.0`。这仅影响 `_snapshot` 中的能量报告（不影响力），但使总能量快照偏离了各核 energy_offset 之和（典型值各自可达 10-30 MeV）。

---

## ✅ 通过项

- ✅ **Coulomb 力未受 nuclear_scale 影响**: `forces_analytical` L216, L237 中 Coulomb 交换和直接力不包含 nuclear_scale — 物理正确
- ✅ **力的数学推导正确**: bulk_prime、对称项导数、表面项 Laplacian 公式均与能量泛函一致
- ✅ **对称项力符号正确**: dE/dρ_n = +Csym/ρ0*(ρ_n-ρ_p), dE/dρ_p = -Csym/ρ0*(ρ_n-ρ_p)
- ✅ **`energy_offset` 不参与力**: 常数偏移对力的贡献为零 (F = -dU/dr, d(常数)/dr = 0) — 物理正确
- ✅ **fragments.py**: `compute_fragment_excitation` 对每个碎片独立计算 GridEDF 能量，使用该 (Z,A) 对应的 `fragment_ground_state_scale`

---

## 🟢 建议

1. **`_density_functional_derivatives` 未使用** (grid_edf.py:356): 定义了但从未调用。`forces_analytical` 内联计算了相同的导数。建议删除死代码或统一使用该方法以减少重复。

2. **建议在松弛中也使用 nuclear_scale**: 修改 `initialize_nucleus` 和 `initialize_and_relax` 中的松弛 GridEDF 使用 `nuclear_scale=grid_nuclear_scale`。这需要先估算 scale 值（可能需两遍松弛），或先做一次带 scale=1.0 的预松弛再切换到实际 scale。

---

## 修复建议

### 方案 A (推荐): 密度依赖的 nuclear_scale

将全局标量 `nuclear_scale` 替换为在 ρ≫ρ0 时趋于 1.0 的密度依赖函数：

```python
# grid_edf.py
def _effective_scale(self, rho: np.ndarray) -> np.ndarray:
    x = np.clip(rho / self.parameters.rho0, 0.0, None)
    # 在正常密度使用校准的 scale，在高密度平滑过渡到 1.0
    w = np.exp(-(x - 1.0)**2 / 0.5)  # 密度权重
    return self.nuclear_scale * w + 1.0 * (1.0 - w)
```

在 `forces_analytical` 中将 `self.nuclear_scale` 替换为 `self._effective_scale(rho)`（逐格点）。

### 方案 B: 重新拟合 Skyrme 参数

直接调整 Skyrme 参数 (α, β, γ, g_tau, η) 以匹配经验结合能，避免使用后验 scale 因子。这更基础但工作量更大。

### 方案 C (应急): 限制 nuclear_scale 下限

将 `_fit_grid_nuclear_scale` 中的 `np.clip(scale, 0.35, 1.25)` 下限提高到 0.90，减小 scale 偏离 1.0 的程度。这会降低静态结合能精度但改善碰撞动力学。

---

## 验证建议

1. 对 U+U b=8 fm 事件分别使用 `nuclear_scale=1.0` 和 `nuclear_scale=grid_nuclear_scale` 运行，对比碎片多重数
2. 在固定 snapshot 时刻输出碰撞重叠区的最大密度 ρ_max，确认是否到达非物理值 (>3ρ0)
3. 用 centroid force 路径 (`use_grid_edf=False`) 运行同一碰撞作为对照
