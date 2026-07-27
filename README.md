# MNT-SIM: ImQMD 多核子转移反应模拟

基于 ImQMD (Improved Quantum Molecular Dynamics) 的 MNT 反应模拟框架。笔记本开发 + nucimp 批量计算。

## 当前状态 (2026-06-28)

### 统一平台

| 平台 | Python | numpy | scipy | CPU |
|------|--------|-------|-------|-----|
| 笔记本 | 3.12.3 | 2.2.6 | 1.16.3 | Intel |
| nucimp | 3.12.3 | 2.2.6 | 1.16.3 | AMD EPYC 48核 |

### 已验证物理

| 系统 | 平台 | 结果 |
|------|------|------|
| ⁸⁶Kr+⁶⁴Ni b=5 (800步) | 笔记本 | Z=32(Ge) + Z=26(Fe) ✅ |
| ⁸⁶Kr+⁶⁴Ni b=5 (800步) | nucimp | Z=37(Rb) + Z=24(Cr) ✅ |
| ⁴⁰Ca 静态 (2000 fm/c) | 两台 | 0.00% 漂移 ✅ |

### 已完成工作（4天）

**实现（Phase 1-4）**：
| 模块 | 文件 | 内容 |
|------|------|------|
| EDF | `skyrme.py`, `grid_edf.py` | IQ1-3b 参数集 + 质心采样 + 网格积分 |
| 核 | `nucleus.py` | 高斯波包 + 密度 + 能量 |
| 初始化 | `initializer.py` | 硬球 + 弛豫 `initialize_and_relax` |
| 传播 | `propagator.py` | RK4 + Surface/Stabilizer 分拆 + GridEDF |
| 碰撞 | `collisions.py` | NN 碰撞 + Pauli + Fermi + 邻居列表 |
| 碎片 | `fragments.py` | MST/iso-MST + 激发能 |
| 退激 | `decay.py` | Weisskopf 蒸发 + 裂变 |

**关键修复**：
| 修复 | 说明 |
|------|------|
| `use_static_reference` 拆分 | → `use_surface_term` + `use_static_stabilizer` |
| 参数对标论文 | MST 3.0→3.5, iso-MST nn/np 2.8→6.0, Fermi 140→170 |
| RNG 锁定 PCG64 | 跨 numpy 1.x/2.x 一致 |
| nucimp 部署 | numpy 1.25→2.2.6, Python 3.10→3.12（`--user` 安装，不影响他人） |
| GridEDF | 271 行网格积分替代质心采样（待参数调优后启用） |

**审计文档**（`papers/` 目录）：
- `imqmd_model_spec.md` — 1026 行完整论文规格
- `spec_vs_code_audit.md` — 88 项差距分析
- `shortcut_audit.md` — 37 项实现捷径
- `diagnosis_and_solutions.md` + `auditor_report.md`
- `numpy_fix_proposal.md` — RNG 修复方案

### 核心参数

| 参数 | 值 | 来源 |
|------|-----|------|
| EDF | IQ2 | Wang 2014 |
| σ_r | 1.1 fm | — |
| MST | 3.5 fm / 250 MeV/c | Z2012 |
| iso-MST nn/np/pp | 6.0/6.0/3.0 fm | iso-MST-R |
| Fermi | 1.0（Wigner 占据数阈值，CoMD 动量交换，能量守恒） | CoMD |
| 弛豫 | 800 fm/c | — |
| 表面力 | ON | 物理 |
| 静态弹簧 | OFF | 反应模式 |

### 待完成（按优先级）

| 优先级 | 任务 | 说明 |
|:--:|------|------|
| P0 | GridEDF 参数调优 | IQ2 在网格上需重校 |
| P1 | 初始化重写 | 中子皮 + w_r/w_p + BE±0.05 |
| P1 | 激发能 E* 校准 | 同尺度假基态参考 |
| P2 | 退激发完整 | 完整 Weisskopf + p/α + 裂变 |
| P2 | 反应截面框架 | b 扫描 + dσ/dZ/dA |
| P3 | Zhao 2016 对比 | U+U 基准产额 |

### 环境

**笔记本**：`~/work/agent work/mnt-sim/`
**nucimp**：`zhuhf@210.77.75.5:1800` → `~/work/mnt-sim/`

nucimp Python：`~/.local/python312/bin/python3`（不影响系统 Python 3.10）

### 文献

14 篇 ImQMD 论文 PDF + Zotero (合集 `ImQMD Model`, userId=9566388)

### 作者

朱浩钒 — 中国科学院近代物理研究所 / 惠州 HIAF
