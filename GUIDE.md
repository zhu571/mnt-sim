# MNT-SIM 代码阅读指南

> 写给会基础 Python 但不熟悉核物理的读者。
> 按顺序从头读到尾，读完你就能理解整个项目在做什么、每个文件每个函数的作用。

---

## 目录

1. [这个项目是干什么的](#1-这个项目是干什么的)
2. [你需要知道的背景知识](#2-你需要知道的背景知识)
3. [安装与运行](#3-安装与运行)
4. [项目文件地图](#4-项目文件地图)
5. [核心思路：三步流水线](#5-核心思路三步流水线)
6. [第一层：核素数据表 (data/)](#6-第一层核素数据表)
7. [第二层：反应截面计算 (cross_section/)](#7-第二层反应截面计算)
8. [第三层：气体单元传输 (transport/)](#8-第三层气体单元传输)
9. [可视化工具 (plot/)](#9-可视化工具)
10. [配置文件详解](#10-配置文件详解)
11. [怎么跑示例](#11-怎么跑示例)
12. [如何修改和扩展](#12-如何修改和扩展)
13. [术语表](#13-术语表)

---

## 1. 这个项目是干什么的

### 一句话概括

用 Python 模拟"两个原子核碰撞后，产物在气体中减速停下来"的全过程。

### 稍微展开

在 HIAF（惠州加速器设施）上，科学家会把一束高速原子核（比如铀-238）打到另一个原子核靶（也是铀-238）上。
碰撞时，两个核之间会"交换"一些质子和中子，产生新的核素——这叫**多核子转移反应**（Multi-Nucleon Transfer, MNT）。

产生的新核素会高速飞出，进入一个充了氦气的容器（**气体单元**），在氦气中逐渐减速、最终停下来，被收集和研究。

这个项目做的就是模拟这两步：
1. **碰撞**：计算会产生哪些新核素，产生多少（截面），以什么角度和能量飞出
2. **减速**：模拟新核素在氦气中飞行、减速、停下来的过程

### 为什么有用

- 设计实验前，先用模拟预测"能不能看到想要的核素"
- 优化气体单元的参数（气压、长度、温度），让更多产物停下来
- 和实验数据对比，验证物理模型

---

## 2. 你需要知道的背景知识

### 2.1 原子核基础

```
原子核 = 质子（Z 个）+ 中子（N 个）
质量数 A = Z + N

例子：
  铀-238 (U-238)：Z=92 个质子, N=146 个中子, A=238
  铅-208 (Pb-208)：Z=82, N=126, A=208
  氙-136 (Xe-136)：Z=54, N=82, A=136
```

不同的 (Z, N) 组合 = 不同的核素（nuclide）。改变 Z 就变成了不同元素。

### 2.2 什么是"截面"（cross section）

截面 σ 是衡量"反应有多容易发生"的量，单位是 **mb**（毫靶恩，1 mb = 10^-27 cm^2）。

- σ 越大 --> 反应越容易发生 --> 实验中能看到更多产物
- σ 越小 --> 越稀有

打个比方：截面就像一个"靶子面积"。靶子越大，子弹越容易打中。

### 2.3 什么是"双微分截面"

普通截面 σ 只告诉你"总共产生多少"。
**双微分截面** d^2σ/dΩdE 告诉你"在某个角度 θ、某个能量 E 处，产生多少"。

这就像不只知道"总共卖了多少杯咖啡"，还知道"每个时间段、每种口味各卖了多少"。

### 2.4 什么是"阻止本领"（stopping power）

高速离子飞进气体后，会和气体原子碰撞，逐渐损失能量。
**阻止本领** dE/dx = 每飞过单位距离损失多少能量。

- 能量高时：主要和电子碰撞（电子阻止），用 **Bethe-Bloch 公式**
- 能量低时：和整个原子碰撞（核阻止），用 **Ziegler/SRIM 参数化**

### 2.5 什么是 Monte Carlo 模拟

Monte Carlo = 用随机数模拟大量粒子的行为，统计结果。

比如模拟 10000 个离子在气体中飞行：
- 每个离子的初始方向稍有不同（随机）
- 每一步的能量损失有随机涨落（散射）
- 最终统计：多少停下来了？停在哪里？

### 2.6 MeV/u 是什么

能量单位。1 MeV = 百万电子伏特。"/u" 表示"每核子"。

所以"7 MeV/u 的铀-238"意味着这个核的总动能 = 7 x 238 = 1666 MeV。

---

## 3. 安装与运行

```bash
# 进入项目目录
cd "~/work/agent work/mnt-sim"

# 安装依赖（系统 Python 3.12）
pip install numpy scipy matplotlib pandas pyyaml

# 验证安装
python -c "from mnt_sim.cross_section.dns import DNSModel; print('OK')"
```

依赖很轻量，只需要 5 个常见科学计算库。

---

## 4. 项目文件地图

```
mnt-sim/
|
+-- mnt_sim/                      <-- Python 包（所有核心代码在这里）
|   +-- __init__.py               <-- 包入口，定义版本号
|   |
|   +-- data/                     <-- [最底层] 核素数据表
|   |   +-- __init__.py           <-- 元素符号表、核素质量、结合能
|   |
|   +-- cross_section/            <-- [第一步] 反应截面计算
|   |   +-- __init__.py           <-- 模块说明
|   |   +-- dns.py                <-- DNS 模型（最基础的截面计算）
|   |   +-- grazing.py            <-- GRAZING 模型（半经典方法）
|   |   +-- empirical.py          <-- 经验公式（快速估算）
|   |   +-- imqmd_hivap.py        <-- ImQMD+HIVAP（最完整的模型）
|   |   +-- imqmd_like.py         <-- ImQMD 变种（含 SRIM 能损）
|   |   +-- imqmd_mc.py           <-- ImQMD Monte Carlo 版
|   |
|   +-- transport/                <-- [第二步] 气体传输模拟
|   |   +-- __init__.py
|   |   +-- stopping.py           <-- 阻止本领计算
|   |   +-- gas_cell.py           <-- 气体单元的几何形状定义
|   |   +-- monte_carlo.py        <-- MC 传输模拟主程序
|   |
|   +-- plot/                     <-- 画图工具
|       +-- __init__.py           <-- 各种可视化函数
|
+-- config/                       <-- 配置文件
|   +-- default.yaml              <-- 默认配置（Xe+Pb 反应）
|   +-- uu.yaml                   <-- U+U 反应配置
|
+-- examples/                     <-- 示例脚本（约 20 个）
|   +-- run_cross_section.py      <-- 跑截面计算
|   +-- run_transport.py          <-- 跑传输模拟
|   +-- run_pipeline.py           <-- 跑完整流水线
|   +-- run_uu_dns.py             <-- U+U 反应 DNS 计算
|   +-- plot_*.py                 <-- 各种绘图脚本
|   +-- debug_*.py                <-- 调试脚本
|
+-- output/                       <-- 输出结果
|   +-- figures/                  <-- 生成的图片
|   +-- *.npz / *.csv             <-- 数值数据
|   +-- uu/                       <-- U+U 反应专用输出
|
+-- wang2020.pdf                  <-- 参考论文
+-- requirements.txt              <-- Python 依赖列表
+-- README.md                     <-- 项目简介
```

**阅读顺序建议**：data/ --> cross_section/dns.py --> transport/stopping.py --> transport/gas_cell.py --> transport/monte_carlo.py

**代码规范**（适用于所有模块）：
- 每个模块顶部都有 `__all__ = [...]`，列出该模块对外暴露的类和函数
- 所有公开函数/类的构造器都有**输入验证**：传入非法参数（负数、空值等）会立即抛 `ValueError`，不会静默算出错误结果
- 元素/核素数据统一从 `mnt_sim.data` 导入，不在各模块重复定义

---

## 5. 核心思路：三步流水线

整个项目的数据流是一条**三步流水线**：

```
+----------------+      +--------------------+      +--------------------+
| 第一步：截面    | ---> | 第二步：双微分截面   | ---> | 第三步：传输模拟    |
| "产生多少"      |      | "什么角度什么能量"   |      | "飞多远停在哪"      |
+----------------+      +--------------------+      +--------------------+

输入：                   中间产物：                    输出：
- 弹核（如 U-238）       - sigma(Z,N) 矩阵            - 沉积位置分布
- 靶核（如 U-238）       - d2sigma/dOmega/dE          - 收集效率
- 束流能量（7 MeV/u）    - 能量-角度分布               - 逃逸比例
- 转移核子范围                                        - 能量损失曲线
```

用代码表示：

```python
# 第一步：计算截面
from mnt_sim.cross_section.dns import DNSModel
dns = DNSModel("U", "U", E_lab=7.0)
result = dns.calculate(delta_Z_range=(-4, 4), delta_N_range=(-6, 6))
# result.sigma_matrix = 一个矩阵，每个格子是一个转移道的截面

# 第二步：角度-能量分布
theta = np.arange(0, 30, 2)
d2sigma = dns.angular_distribution(theta, result)
# d2sigma[iZ, iN, iTheta, iE] = 在某角度某能量的截面

# 第三步：传输模拟
from mnt_sim.transport.gas_cell import GasCell
from mnt_sim.transport.monte_carlo import TransportMC
cell = GasCell(gas='He', pressure_mbar=50, length_mm=200)
mc = TransportMC(cell, ion_Z=92, ion_A=243, energy=800, n_particles=1000)
result = mc.run()
print(f"收集效率: {result.deposition_efficiency*100:.1f}%")
```

---

## 6. 第一层：核素数据表

**文件**：`mnt_sim/data/__init__.py`

这是整个项目最基础的一层——提供核素的基本物理数据。

### 6.1 元素符号表

项目维护了两张互为反向的映射表：

```python
# Z -> 符号（用于把数字转成人能读的名字）
ELEMENT_SYMBOLS = {
    0: 'n', 1: 'H', 2: 'He', 3: 'Li', ..., 92: 'U', ..., 118: 'Og'
}

# 符号 -> Z（所有截面模块统一查这张表）
ELEMENTS = {'n': 0, 'H': 1, 'He': 2, ..., 'U': 92, ...}
```

还有一张**典型质量数表**，给简化模型用（不指定 A 时取默认值）：

```python
TYPICAL_MASS_NUMBERS = {
    'Xe': 136, 'Pb': 208, 'U': 238, 'Ca': 40, 'Pt': 198, ...
}
```

### 6.2 质量过剩表 MASS_EXCESS

```python
MASS_EXCESS = {
    (92, 146): -58.164,  # U-238 的质量过剩 = -58.164 MeV
    (82, 126): -73.175,  # Pb-208
    ...
}
```

**质量过剩**是什么？原子核的实际质量 ≠ Z x 质子质量 + N x 中子质量。差值就是质量过剩。
它和结合能直接相关，在计算反应 Q 值（能量释放/吸收）时会用到。

### 6.3 关键函数

| 函数 | 输入 | 输出 | 用途 |
|------|------|------|------|
| `element_to_z(symbol)` | 元素符号如 'U' | 原子序数 92 | **所有截面模块统一调用这个**，不再各自维护映射表 |
| `typical_mass_number(symbol)` | 元素符号 | 默认质量数 | 不指定 A 时的默认值（如 'U' -> 238） |
| `get_nuclide(Z, A)` | 原子序数, 质量数 | Nuclide 对象 | 创建一个核素对象 |
| `mass_excess(Z, N)` | 质子数, 中子数 | MeV | 查表或用半经验公式估算 |
| `target_projectile_pairs()` | 无 | 列表 | 常用的弹靶组合 |

> 注：传入无效符号时会抛 `ValueError`，不会静默返回 0。

### 6.4 Nuclide 数据类

```python
@dataclass
class Nuclide:
    Z: int       # 质子数
    N: int       # 中子数
    A: int       # 质量数 = Z + N
    symbol: str  # 元素符号，如 'U'
    mass: float  # 原子质量 [MeV/c^2]
    abundance: float  # 天然丰度
```

就是一个"核素身份证"，把一个核素的所有基本信息打包在一起。

---

## 7. 第二层：反应截面计算

**目录**：`mnt_sim/cross_section/`

这是项目的核心计算部分。提供了 4 种不同的物理模型来计算截面，从简单到复杂：

```
简单 <---------------------------------> 复杂

经验模型          GRAZING          DNS         ImQMD+HIVAP
(快速估算)      (半经典)       (主方程)      (最完整)
```

### 7.1 DNS 模型 -- dns.py（最核心，最先读）

**DNS = Di-Nuclear System（双核系统）**

#### 物理图像

想象两个核碰撞时，就像两个水滴轻轻碰到一起：
- 它们在接触面保持各自形状（不融合）
- 质子和中子可以通过接触面从一边"扩散"到另一边
- 过一段时间后分开，各自带着不同数量的核子

```
碰撞前：  O --> <-- @         碰撞中：  O@       分开后：  o  @
         弹核    靶核              核子交换           新核素
         (Z1,N1) (Z2,N2)                         (Z1+dZ, N1+dN)
```

#### 代码结构

```python
class DNSModel:
    def __init__(self, projectile, target, E_lab):
        # 1. 查元素的 Z 和 A（调用共享数据模块）
        self.Zp = self._element_to_Z(projectile)  # -> element_to_z('U') = 92
        self.Ap = self._element_to_A(projectile)  # -> typical_mass_number('U') = 238
        # ...靶核同理

        # 2. 计算运动学参数
        self.E_cm = self._compute_E_cm()           # 质心系能量
        self.V_C = self._compute_coulomb_barrier()  # 库仑势垒
```

> **新代码变化**：`_element_to_Z()` 和 `_element_to_A()` 不再各自维护一大段硬编码字典，
> 而是只有一行 `return element_to_z(symbol)` / `return typical_mass_number(symbol)`，
> 数据统一在 `mnt_sim/data/__init__.py` 管理。其他模型（GRAZING、ImQMD 等）也一样。

**关键概念解释**：

| 变量 | 含义 | 比喻 |
|------|------|------|
| `E_lab` | 实验室系束流能量 (MeV/u) | 子弹飞行速度 |
| `E_cm` | 质心系能量 (MeV) | "两人都在走，相对速度才是真正的碰撞速度" |
| `V_C` | 库仑势垒 (MeV) | "两个正电荷核之间的排斥壁垒"——能量必须超过它才能碰到 |
| `E_cm_eff` | E_cm - V_C | 超过势垒后还剩多少能量可以"做事" |

#### calculate() 方法 -- 核心计算

这是整个 DNS 模型的核心，做的事情是：

```
输入：dZ 范围 (-4 到 +4)、dN 范围 (-6 到 +6)
输出：P(dZ, dN) 概率矩阵 --> sigma(dZ, dN) 截面矩阵
```

算法步骤：

```
1. 创建一个 (dZ x dN) 的概率矩阵 P，初始值全为 0，只有 P[0,0] = 1
   （一开始没有核子转移，概率 100% 在"不转移"状态）

2. 循环 n_steps 次（时间演化）：
   对矩阵中每个格子 (i, j)：
     a. 计算转移率 Lambda（核子跳过去的速率）
     b. 从相邻格子 "流入" 概率（gain 项）
     c. 从当前格子 "流出" 概率（loss 项）
   归一化

3. 用 Wong 公式算总反应截面 sigma_total

4. 各转移道的截面 = sigma_total x P(dZ, dN)
```

**通俗理解**：就像一个棋盘上的墨水扩散。一开始墨水点在中心 (0,0)，然后每一步都往四周扩散一点。最终的墨水分布就是各个转移道的概率。

#### _binding_energy() -- 半经验质量公式

```python
def _binding_energy(self, Z, N):
    B = (a_v * A           # 体积项：核子越多，结合越强
         - a_s * A^(2/3)   # 表面项：表面核子结合弱一些
         - a_c * Z*(Z-1)/A^(1/3)  # 库仑项：质子排斥
         - a_a * (A-2*Z)^2/A      # 对称项：Z约等于N时最稳定
         + 配对项)                  # 偶偶核更稳定
```

这是核物理中的"万能近似公式"，不需要查表就能估算任何核素的结合能。

#### angular_distribution() -- 角度分布

把截面矩阵"展开"成 d2sigma(theta, E) 的形式：
- 角度分布：以掠射角（grazing angle）为中心的高斯分布
- 能量分布：以质心能量为中心的高斯分布
- 宽度随转移核子数增加而变宽

#### DNSResult 数据类

```python
@dataclass
class DNSResult:
    Z: np.ndarray          # dZ 的值，如 [-4,-3,...,+4]
    N: np.ndarray          # dN 的值，如 [-6,-5,...,+6]
    P: np.ndarray          # 概率矩阵 P[iZ, iN]
    E_cm: float            # 质心能量
    sigma_total: float     # 总截面 [mb]
    sigma_matrix: np.ndarray  # 各道截面 [mb]

    def to_dataframe(self):  # 转成 pandas 表格，方便查看
```

---

### 7.2 GRAZING 模型 -- grazing.py

**特点**：比 DNS 更"量子化"——按部分波展开（partial wave）计算。

#### 什么是部分波？

把入射束流按角动量 l 分解。每个 l 对应不同的"碰撞距离"：
- l = 0：迎头碰撞
- l 很大：擦边而过
- l_grazing（掠射角动量）：恰好碰到核表面

```
l = 0:    --> O@ <--    迎头碰
l 中等:   --> O @       擦边
l 大:     -->   O  @    飞过去了
```

#### 代码流程

```python
def calculate(self):
    for l in range(n_l):          # 对每个部分波 l
        T_l = Hill_Wheeler(l)     # 1. 计算透射系数（能不能碰到）
        P_l = transfer_prob(l)    # 2. 核子转移概率
        sigma_ch += sigma_l * P_l # 3. 贡献到各转移道
```

和 DNS 的区别：DNS 用时间演化求转移概率，GRAZING 用量子隧穿概率。

---

### 7.3 经验模型 -- empirical.py

**特点**：最简单，不做动力学计算，直接用经验公式估算截面。

```python
class EmpiricalModel:
    def estimate(self):
        for dp, dn in channels:
            sigma = sigma_0 * Q_factor * mass_factor * iso_factor
```

三个因子：
- **Q_factor**：反应 Q 值越接近最优值，截面越大
- **mass_factor**：转移的核子越多，截面越小（指数衰减）
- **iso_factor**：同位旋匹配因子

**用途**：快速估算，5 秒出结果。适合初步判断"这个反应有没有戏"。

---

### 7.4 ImQMD+HIVAP 模型 -- imqmd_hivap.py（最复杂）

**特点**：最完整的模型，包含了全部物理过程。

#### 物理过程链

```
部分波展开 --> 核子转移 --> 激发态退激（HIVAP）--> 实验室系运动学变换
    |            |            |                    |
  T_l(l)    P(dZ,dN,l)   中子蒸发              E_lab, theta_lab
                          (冷却到基态)
```

#### 什么是 HIVAP？

碰撞后产生的核素是"激发态"的（内部很热），需要通过蒸发中子来冷却：

```
U-245* (E*=40 MeV)  -->  U-245* - n --> U-244* - n --> U-243 (基态)
                          蒸发第1个n      蒸发第2个n     冷下来了
```

每蒸发一个中子，大约带走 ~6-8 MeV 的能量。

代码中：

```python
def _hivap_evaporation(self, Z_prim, A_prim, E_star):
    """计算初级碎片蒸发中子后能到达哪些最终核素"""
    B_n = 6.0   # 中子分离能 约 6 MeV
    n_max = int(E_star / B_n)  # 最多能蒸发几个
    for n in range(n_max + 1):
        A_final = A_prim - n
        P_n = exp(-n * B_n / T_nuc)  # 蒸发概率
```

#### 两体运动学 _two_body_kinematics()

碰撞在质心系（CM）发生，但实验在实验室系（Lab）测量。需要做坐标变换：

```python
def _two_body_kinematics(self, m_C, m_D, theta_cm_deg, Q):
    """质心系 --> 实验室系 的速度/能量/角度变换"""
    # 1. 质心系中碎片 C 的能量
    E_cm_C = (m_D / (m_C + m_D)) * (E_cm + Q)
    # 2. 质心系中碎片 C 的速度
    v_cm_C = sqrt(2 * E_cm_C / m_C)
    # 3. 叠加质心运动速度，得到实验室系速度
    v_lab^2 = v_cm^2 + v_cm_C^2 + 2 * v_cm * v_cm_C * cos(theta_cm)
    # 4. 计算实验室角度
    theta_lab = arctan2(v_cm_C * sin(theta_cm), v_cm + v_cm_C * cos(theta_cm))
```

**比喻**：你在火车上扔球（质心系），地面上的人看到球的速度 = 火车速度 + 球相对火车的速度。

---

### 7.5 ImQMD 变种 -- imqmd_like.py 和 imqmd_mc.py

| 文件 | 特点 |
|------|------|
| `imqmd_like.py` | ImQMD + SRIM 靶能损（连续模式），考虑碎片在靶中的能量损失 |
| `imqmd_mc.py` | Monte Carlo 事件生成器，逐事件模拟（更像真实实验） |

`imqmd_mc.py` 的特色是用随机采样代替遍历：

```python
def generate_events(self, n_events=50000):
    for _ in range(n_events):
        l = 随机选一个部分波（按截面权重）
        dZ, dN = 随机选一个转移道（按转移概率权重）
        E_lab, theta_lab = 运动学变换
        E_out = SRIM 靶能损
        记录这个事件 (theta_lab, E_out, weight)
```

---

## 8. 第三层：气体单元传输

**目录**：`mnt_sim/transport/`

### 8.1 阻止本领 -- stopping.py

**核心问题**：一个高速离子在氦气中飞行时，每飞 1 cm 损失多少能量？

#### 三个能区，三种公式

```
能量 E/A:   <-- 低 ------------------- 中 ------------------- 高 -->

            Lindhard-Scharff     Ziegler/SRIM      Bethe-Bloch
            (< 25 keV/u)        (25 keV/u ~ 1 MeV/u) (> 1 MeV/u)
            核散射为主           过渡区              电子碰撞为主
```

代码中 `ziegler_stopping()` 方法自动选择：

```python
def ziegler_stopping(self, E_MeV):
    E_per_u = E_MeV / self.ion_A
    if E_per_u > 1.0:
        return self.bethe_bloch(E_MeV)       # 高能
    elif E_per_u > 0.025:
        return Ziegler参数化(E_MeV)            # 中能
    else:
        return self._lindhard_scharff(E_MeV)  # 低能
```

#### Bethe-Bloch 公式（高能区）

这是最经典的能量损失公式：

```
dE/dx ~ z^2 * Z_gas / (beta^2 * A_gas) * [ln(2*me*beta^2*gamma^2*W_max/I^2) - beta^2 - delta/2 - C/Z]
```

其中：
- z = 离子电荷
- Z_gas, A_gas = 气体的原子序数和质量数
- beta = v/c（速度/光速）
- I = 气体的平均激发能（He: 41.8 eV）

#### StoppingPower 类

```python
sp = StoppingPower(ion_Z=92, ion_A=243, gas='He', pressure_mbar=50)
dEdx = sp(E_MeV=100)  # 调用返回 dE/dx [MeV/(mg/cm^2)]
```

#### 射程计算 range_energy()

积分 1/dEdx 得到射程：

```python
range_cm, trajectory = sp.range_energy(E0_MeV=800)
# trajectory[0] = 距离数组 [cm]
# trajectory[1] = 能量数组 [MeV]
```

#### 能量歧离 energy_straggling()

能量损失不是精确值，有统计涨落（像骰子有随机性）。Bohr 歧离给出涨落的方差。

---

### 8.2 气体单元 -- gas_cell.py

定义气体单元的"容器"。

#### GasCell 数据类

```python
cell = GasCell(
    gas='He',              # 气体种类
    pressure_mbar=50,      # 气压 [mbar]
    length_mm=200,         # 长度 [mm]（沿束流方向）
    diameter_mm=40,        # 直径 [mm]
    temperature_K=293,     # 温度 [K]（293K = 室温）
)
```

#### 三种几何形状

```python
class CellGeometry(Enum):
    CYLINDRICAL = "cylindrical"  # 圆柱形（最常见）
    CONICAL = "conical"          # 锥形入口
    TAPERED = "tapered"          # 渐缩形（出口收窄）
```

#### 边界检测 is_within_cell()

```python
def is_within_cell(self, x_mm, y_mm, z_mm) -> bool:
    """粒子还在容器里吗？"""
    # z 方向：0 <= z <= length
    # 径向：r = sqrt(x^2 + y^2) <= diameter/2
```

这个函数在 MC 模拟中每一步都会调用——如果粒子飞出了边界，就标记为"逃逸"。

#### 窗口材料 Foil

```python
@dataclass
class Foil:
    material: WindowMaterial  # Ti, Havar, Mylar, Ni
    thickness_um: float       # 厚度 [um]
    diameter_mm: float        # 直径 [mm]
```

真实的气体单元有入射窗（防止气体泄漏但让束流穿过）。目前代码中窗口能损还没有完全集成到传输模拟中。

---

### 8.3 Monte Carlo 传输 -- monte_carlo.py（最值得仔细读）

这是传输模拟的主程序。

#### TransportMC 类初始化

```python
mc = TransportMC(
    gas_cell=cell,        # 用哪个气体单元
    ion_Z=92, ion_A=243,  # 模拟哪种离子
    energy=800,           # 初始能量 [MeV]
    n_particles=1000,     # 模拟多少个粒子
    sigma_theta=3.0,      # 初始角展度 [度]
)
```

#### run() 方法 -- 核心模拟循环

```
对每个粒子 (共 N 个):
    1. 给定初始位置 (入口)、初始能量、随机初始角度
    2. while 粒子还活着:
        a. 计算这一步的能量损失 dE
        b. 加上随机的能量歧离
        c. 计算多次散射带来的角度偏转
        d. 更新粒子位置
        e. 检查是否出界 --> 标记"逃逸"
        f. 检查能量是否 <= 0 --> 标记"停下来了"
    3. 记录最终状态（位置、能量）
```

每一步的细节在 `_step_particle()` 方法中：

```python
def _step_particle(self, state, step_size_mm=0.5):
    # 1. 能量损失
    dEdx = self.stopping(E)       # 查阻止本领
    dE = dEdx * dx_mgcm2          # 损失 = 阻止本领 x 路径长度

    # 2. 能量歧离（随机涨落）
    omega2 = self.stopping.energy_straggling(E, dx_mgcm2)
    dE_strag = np.random.normal(0, sqrt(omega2))

    # 3. 多次散射（角度偏转）
    theta_scat, phi_scat = self._multiple_scattering(E, dx_cm)

    # 4. 更新位置
    state.x += step_size * sin(theta) * cos(phi)
    state.y += step_size * sin(theta) * sin(phi)
    state.z += step_size * cos(theta)  # z = 束流方向

    # 5. 边界检查
    if not self.cell.is_within_cell(x, y, z):
        state.escaped = True
```

#### 多次散射 _multiple_scattering()

用 **Highland 公式**计算散射角的均方根值：

```
theta_rms = (13.6 / (beta * p)) * z * sqrt(x/X0) * (1 + 0.038 * ln(x/X0))
```

其中 X0 是辐射长度（He: 94.3 g/cm^2）。

散射角服从正态分布，所以每一步都用 `np.random.normal(0, theta_rms)` 采样。

#### TransportResult 数据类

```python
@dataclass
class TransportResult:
    n_particles: int          # 总粒子数
    n_deposited: int          # 停在容器内的粒子数
    n_escaped: int            # 飞出容器的粒子数
    deposition_positions: np.ndarray  # 停下位置 (N,3)
    deposition_energies: np.ndarray   # 停下时的剩余能量
    escape_energies: np.ndarray       # 逃逸时的剩余能量

    @property
    def deposition_efficiency(self):
        """收集效率 = 停下的/总数"""
        return self.n_deposited / self.n_particles
```

#### 参数扫描 scan_pressure() / scan_length()

```python
# 扫描不同气压下的收集效率
results = mc.scan_pressure(
    pressures=[20, 40, 60, 80, 100],  # 扫描 5 个气压点
    n_particles=500
)
# results[20] = TransportResult at 20 mbar
# results[40] = TransportResult at 40 mbar
# ...
```

---

## 9. 可视化工具

**文件**：`mnt_sim/plot/__init__.py`

提供 5 个绘图函数，全部基于 matplotlib：

| 函数 | 画什么 | 输入 |
|------|--------|------|
| `plot_cross_section_matrix()` | Z-N 平面截面热力图 | DNSResult |
| `plot_angle_energy_map()` | theta-E 双微分截面等高线 | d2sigma 数组 |
| `plot_deposition_distribution()` | 沉积位置分布（4 个子图） | TransportResult |
| `plot_range_curve()` | 能量-距离曲线 + 阻止本领曲线 | StoppingPower |
| `plot_pressure_scan()` | 效率 vs 气压 | scan_pressure 结果 |

每个函数都支持 `save_path` 参数直接保存图片。

---

## 10. 配置文件详解

**文件**：`config/default.yaml` 和 `config/uu.yaml`

```yaml
# === 反应体系 ===
reaction:
  projectile: "U"       # 弹核元素
  projectile_A: 238     # 弹核质量数
  target: "U"           # 靶核元素
  target_A: 238         # 靶核质量数
  E_lab: 7.0            # 束流能量 [MeV/u]

# === 截面计算设置 ===
cross_section:
  model: "dns"          # 用哪个模型：dns / grazing / empirical
  delta_Z: [-6, 6]      # 质子转移范围
  delta_N: [-10, 10]    # 中子转移范围
  n_steps: 500          # DNS 时间步数（越大越精确，越慢）

# === 气体单元参数 ===
gas_cell:
  gas: "He"             # 气体种类
  pressure_mbar: 80     # 气压 [mbar]
  length_mm: 1600       # 长度 [mm]
  diameter_mm: 1600     # 直径 [mm]
  temperature_K: 90     # 温度 [K]（90K = 低温）

# === MC 传输设置 ===
transport:
  n_particles: 10000    # MC 粒子数（越大统计越好，越慢）
  step_size_mm: 0.5     # 步长 [mm]
  sigma_theta_deg: 3.0  # 初始角展度 [度]
```

---

## 11. 怎么跑示例

### 11.1 最简单的截面计算

```bash
cd "~/work/agent work/mnt-sim"
python examples/run_cross_section.py
```

或者在 Python 中交互运行：

```python
import sys
sys.path.insert(0, '.')
from mnt_sim.cross_section.dns import DNSModel

# 创建模型：Xe-136 + Pb-208 @ 8 MeV/u
dns = DNSModel("Xe", "Pb", E_lab=8.0)
print(f"质心能量: {dns.E_cm:.1f} MeV")
print(f"库仑势垒: {dns.V_C:.1f} MeV")

# 计算截面
result = dns.calculate()
print(f"总截面: {result.sigma_total:.1f} mb")

# 查看最大截面的通道
df = result.to_dataframe()
top5 = df.nlargest(5, 'sigma_mb')
print(top5[['Z', 'N', 'A', 'sigma_mb']])
```

### 11.2 U+U 反应（论文重现）

```python
from mnt_sim.cross_section.imqmd_hivap import ImQMD_HIVAP_Model

model = ImQMD_HIVAP_Model("U", "U", E_lab=7.0, Ap=238, At=238)
theta, E, d2sigma = model.calculate_d2sigma()
# d2sigma[i_theta, i_E] = 双微分截面
# 可以画 contour 图
```

### 11.3 完整流水线

```python
from mnt_sim.transport.monte_carlo import run_doublediff_example
result_xs, result_transport = run_doublediff_example()
```

这个函数会依次执行截面计算 --> 角度分布 --> 传输模拟，打印所有中间结果。

---

## 12. 如何修改和扩展

### 12.1 换一个反应体系

修改 `config/default.yaml` 或创建新的配置文件：

```yaml
reaction:
  projectile: "Ca"      # 改成 Ca-40
  projectile_A: 40
  target: "Pb"          # 打 Pb-208
  target_A: 208
  E_lab: 7.0
```

### 12.2 换一种气体

目前支持的气体在 `stopping.py` 的 `GAS_PROPERTIES` 字典中：

```python
GAS_PROPERTIES = {
    'He':  (2,  4.0026, 1.786e-4),   # Z, A, 密度@STP
    'H2':  (1,  1.0079, 8.988e-5),
    'N2':  (7,  14.0067, 1.251e-3),
    'Ar':  (18, 39.948,  1.784e-3),
    'Ne':  (10, 20.1797, 9.002e-4),
    'CH4': (6,  16.043,  7.168e-4),
}
```

要添加新气体，在这个字典里加一行即可。

### 12.3 调整 MC 精度

- **粒子数** (`n_particles`)：1000 = 快速测试，10000 = 正式计算，100000 = 高精度
- **步长** (`step_size_mm`)：0.5 mm = 默认，0.1 mm = 更精确但慢 25 倍
- **DNS 步数** (`n_steps`)：200 = 默认，500 = 重系统，50 = 快速测试

### 12.4 添加新的截面模型

在 `cross_section/` 目录下创建新文件，实现 `calculate()` 方法，返回一个带有 `Z, N, sigma_matrix, sigma_total, E_cm` 属性的结果对象。

---

## 13. 术语表

| 术语 | 英文 | 含义 |
|------|------|------|
| 弹核 | Projectile | 加速器打出的束流粒子 |
| 靶核 | Target | 被打的固定靶 |
| 核子 | Nucleon | 质子或中子的统称 |
| 转移道 | Transfer channel | 一种特定的 (dZ, dN) 转移组合 |
| 截面 | Cross section (sigma) | 反应发生概率的度量，单位 mb |
| 双微分截面 | d2sigma/dOmega/dE | 按角度和能量展开的截面 |
| 库仑势垒 | Coulomb barrier | 两个正电核之间的静电排斥壁垒 |
| 质心系 | Center of mass (CM) | 质心静止的参考系 |
| 实验室系 | Lab frame | 靶核静止的参考系 |
| 部分波 | Partial wave (l) | 按角动量量子数分解的入射波 |
| 透射系数 | Transmission T_l | 部分波穿过势垒的概率 |
| 阻止本领 | Stopping power dE/dx | 单位路径长度的能量损失 |
| 能量歧离 | Energy straggling | 能量损失的统计涨落 |
| 掠射角 | Grazing angle | 恰好碰到核表面的散射角 |
| 退激 | De-excitation | 激发态核释放能量回到基态 |
| 中子蒸发 | Neutron evaporation | 通过发射中子释放激发能 |
| SRIM | -- | 离子阻止与射程计算程序（经验数据库） |
| MC | Monte Carlo | 用随机数做统计模拟 |
| CSDA | Continuous Slowing Down Approximation | 连续减速近似 |

---

> 最后更新：2026-06-19
> 作者：Hermes Agent（为浩钒生成）
