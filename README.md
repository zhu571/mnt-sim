# MNT-SIM: Multi-Nucleon Transfer Reaction Simulation

多核子转移反应截面计算 + 气体单元传输沉积模拟工作区

## 概述

本项目提供一套完整的 Python 框架，用于：

1. **MNT 反应截面计算** — DNS（双核系统）模型与经验参数化方法，计算多核子转移反应的微分截面
2. **双微分截面生成** — 生成 d²σ/dE/dΩ 格式的双微分截面，作为气体单元传输的输入
3. **气体单元传输模拟** — 基于 Monte Carlo 方法，模拟反冲核在气体单元（如 HIAF/IMP 低温气体单元）内的能量损失、射程分布与沉积位置

## 项目结构

```
mnt-sim/
├── mnt_sim/
│   ├── cross_section/        # 反应截面计算
│   │   ├── dns.py            # DNS (双核系统) 模型
│   │   ├── grazing.py        # GRAZING-like 半经典方法
│   │   └── empirical.py      # 经验系统学参数化
│   ├── transport/            # 气体单元传输
│   │   ├── stopping.py       # 阻止本领 (Bethe-Bloch, Ziegler)
│   │   ├── monte_carlo.py    # MC 传输模拟
│   │   └── gas_cell.py       # 气体单元几何与边界条件
│   ├── data/                 # 参考数据与核素表
│   │   └── nuclides.py       # 核素质量、电荷等数据
│   └── plot/                 # 可视化工具
│       └── plotting.py       # 截面图、能谱图、空间分布图
├── examples/
│   ├── run_cross_section.py  # 截面计算示例
│   └── run_transport.py      # 传输模拟示例
├── config/
│   └── default.yaml          # 默认配置文件
├── notebooks/                # Jupyter 交互式分析
└── README.md
```

## 安装

```bash
cd "D:/work/agent work/mnt-sim"
pip install -r requirements.txt
```

依赖：numpy, scipy, matplotlib, pandas, pyyaml

## 快速入门

### 1. 计算反应截面

```python
from mnt_sim.cross_section.dns import DNSModel

dns = DNSModel(projectile="Xe", target="Pb", E_lab=8.0)  # MeV/u
xs = dns.calculate(delta_Z=2, delta_N=2)  # 2p2n 转移道
xs.to_dataframe()
```

### 2. 模拟气体单元传输

```python
from mnt_sim.transport.gas_cell import GasCell
from mnt_sim.transport.monte_carlo import TransportMC

cell = GasCell(gas="He", pressure=50, length=200)  # mbar, mm
sim = TransportMC(cell, ion="Pb-208", energy=10.0)  # MeV
result = sim.run(n_particles=10000)
result.plot_deposition()
```

## 物理模型

### DNS 模型

双核系统 (Di-Nuclear System) 模型基于在相互作用过程中 projectile 和 target 保持各自独立核子体系的假设，核子通过势垒扩散转移。采用主方程 (Master Equation) 描述核子转移概率的时间演化：

```
dP(Z,N,t)/dt = Σ ΔZ,ΔN [Λ(Z-ΔZ,N-ΔN→Z,N) P(Z-ΔZ,N-ΔN,t)
                         - Λ(Z,N→Z+ΔZ,N+ΔN) P(Z,N,t)]
```

### 气体单元能量损失

反冲核在 He 气中的能量损失采用 Bethe-Bloch 公式（高能区）和 Ziegler 参数化（低能区）计算，结合 Monte Carlo 模拟能散（straggling）和多次散射。

## 引用

如使用本代码进行研究，请引用：

> 朱浩钒 et al., 多核子转移反应在低温气体单元中的传输模拟, *Nuclear Science and Techniques* (2026)

## 作者

朱浩钒 (Hao-Fan Zhu) — 中国科学院上海应用物理研究所 / 上海同步辐射光源
