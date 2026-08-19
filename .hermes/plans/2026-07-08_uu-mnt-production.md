# U+U MNT 反应生产级改进计划

> **For Hermes:** Use plan mode — no execution until reviewed.

**Goal:** 将当前代码从"Kr+Ni 可跑"提升到"U+U 多核子转移反应可产出截面"

## 核心诊断 (7 次 U+U 扫描总结)

| # | 配置 | 碎片 | 转移 | 结论 |
|:--:|------|:--:|:--:|:--:|
| 1 | **Fermi 全关 + MST** | **15-21** | Z=80-83 ✅ | **唯一成功** |
| 2-7 | Fermi=170/CoMD/Group/Sigmoid | 50-182 | 碎裂/爆碎 | 全失败 |

**根因:** Fermi 约束在重离子碰撞中必须关闭。论文 Fortran 版使用 σ_r(A) 自适应宽度 + 未公开调参使 CoMD 可用，我们的 Python 实现用固定 σ_r=1.1 fm 做不到。7 次尝试全失败。**接受 Fermi-off 作为工作假设。**

**次因:** HIVAP Fortran 过度蒸发碎片的 E*，因为论文给 HIVAP 的是 E*=30 MeV 级碎片，我们的是 80-1200 MeV。**不可接论文原版 HIVAP, 用 Python decay 过渡。**

## 其他已知差距 (接受, 不修)

- U-238 形变 β₂≈0.28: 论文随机取向但 Fermi-off 已有 Z=80-83 转移。跳过。
- GridEDF 高密度 EOS 软化 (nuclear_scale): sigmoid 尝试无效。已有 E* 不高于论文。跳过。
- CoMD Wigner 占有数约束: σ_r(A) 重参化需数月。跳过。

**策略:** 不动 Fermi/EDF/碰撞。补冷却 + 确认 Python decay 可用 → 跑 U+U 生产。

---

## Task 1: 制作当前版本基线

**目标:** 记录当前 commit 的 Kr+Ni 和 U+U 完整结果

**Step 1:** 在 nucimp 服务器上跑 Kr+Ni 生产扫描, 保存 dσ/dZ

```python
impact_parameter_scan(36,86,28,64,7.0,
    b_max=14,delta_b=0.5,events_per_b=20,
    time_fm_c=800,n_workers=8,fragment_method='mst',use_hivap=False)
# 过滤 Z>=3, 保存 JSON
```

**Step 2:** 在 nucimp 上跑 U+U 基准 (b=[5,8,11], 3事件/b, 1000 fm/c)

**验证:** 两个 JSON 文件存在, U+U 有 Z=80-83 转移产物

---

## Task 2: 传播后冷却

**目标:** 传播完成后, 碎片识别前, 进行 100 fm/c 无碰撞冷却

**文件:** `mnt_sim/imqmd/reaction.py:run_imqmd_event`

**修改:** `propagate()` 后加一行:
```python
propagate(system, dt=1.0, n_steps=100, with_collisions=False,
          use_surface_term=True, use_static_stabilizer=False, use_grid_edf=True)
```

**验证:** U+U b=8 碎片数从 15-21 降到 4-8

---

## Task 3: 确认 Python decay 可用

**目标:** 确认 `evaporate_full()` 对 U+U 热碎片不产生 Z=0/1 主导

**文件:** 不改代码, 运行验证

**验证:** 单事件 U+U b=8, 最终产物 Z≥3 的截面 > 总截面的 50%

---

## Task 4: U+U 生产扫描

**目标:** 完整 U+U b-扫描 → dσ/dZ

```python
impact_parameter_scan(92,238,92,238,7.0,
    b_max=18,delta_b=1.0,events_per_b=10,
    time_fm_c=1000,n_workers=8,
    fragment_method='mst',use_hivap=False)
# 过滤 Z>=3
```

**验证:** Z=80-90 区域有峰, σ_total ~1000-2000 mb, 与 Zhao 2016 Fig.1 定性对比
