# NumPy 版本依赖根因分析与 RNG 修复方案

## 一、问题现象

ImQMD 代码（`mnt_sim/imqmd/` 下全部 7 个文件）在两组环境中产生**系统性不同**的物理结果：

| 环境 | NumPy 版本 | Python | CPU | 结果 |
|------|-----------|--------|-----|------|
| 笔记本电脑 | 2.4.6 | 3.12 | Intel | 正确核物理碎片分布 |
| nucimp 服务器 | 1.25.2 | 3.10 | AMD EPYC | 完全核解体（disintegration） |

注意：**相同代码、相同参数**（费米阈值 170.0、弛豫时间 800.0 fm/c、表面力 gsur 等）在不同 NumPy 版本上产生的不只是噪声级别的差异，而是**系统性**的完全不同——核碎片模式从合理的重离子碰撞碎片分布变成所有核子散开的完全解体。

## 二、根因分析

### 2.1 所有 RNG 使用点

代码中共有 **2 个独立的 `np.random.Generator` 实例**，以及它们的 **7 处调用点**：

#### RNG 实例 A：`collisions.py` — `_rng()` 函数（第 16–19 行）

```python
def _rng(nucleus: ImQMDNucleus) -> np.random.Generator:
    if not hasattr(nucleus, "_collision_rng"):
        nucleus._collision_rng = np.random.default_rng(73129)
    return nucleus._collision_rng
```

**硬编码种子 `73129`**，所有核碰撞共享同一个 RNG（通过 nucleus 对象上的延迟初始化属性）。该 RNG 在以下位置被调用：

| 位置 | 行号 | 调用 | 作用 |
|------|------|------|------|
| `attempt_nn_collision()` | 324 | `rng.shuffle(order)` | 随机打乱碰撞候选对顺序 |
| `_scatter_isotropic()` | 227 | `rng.normal(size=3)` | 生成随机散射方向 |
| `attempt_nn_collision()` | 366 | `rng.random()` | 泡利阻塞概率判定 |
| `fermi_constraint_check()` | 406 | `rng.normal(size=3)` | 动量修正方向 |

#### RNG 实例 B：`initializer.py` — `initialize_nucleus()` 函数（第 92 行）

```python
rng = np.random.default_rng(seed if seed is not None else 1000 + 17 * Z + A)
```

**按 (Z, A) 计算种子**（无显式传入时），每个核有不同的初始化种子。该 RNG 在以下位置被调用：

| 位置 | 行号 | 调用 | 作用 |
|------|------|------|------|
| `_sample_hard_sphere()` | 23 | `rng.uniform(-radius, radius)` | 硬球位置采样 |
| `initialize_nucleus()` | 97 | `rng.shuffle(is_proton)` | 随机打乱质子/中子标签 |
| `_fermi_momenta()` | 46, 48 | `rng.normal()`, `rng.random()` | 局域费米动量采样 |

### 2.2 NumPy 1.x vs 2.x 的 RNG 差异

**核心根因：`np.random.default_rng(seed)` 在 NumPy 1.x 和 2.x 中使用不同的默认 BitGenerator。**

| NumPy 版本 | `default_rng()` 底层算法 | 相同种子的输出序列 |
|-----------|------------------------|-------------------|
| 1.17 – 1.26 | **PCG64** | 序列 A |
| 2.0+ | **PCG64DXSM** | 序列 B（与序列 A 完全不同） |

虽然两个算法都是高质量、可复现的 PRNG，但**相同的种子在两种算法中产生完全不同的随机数序列**。这不是 bug——NumPy 2.0 更改默认 BitGenerator 是经过文档公告的有意行为。但对于像 ImQMD 这样混沌的物理模拟，这导致**相同的代码和参数产生系统性的不同物理结果**。

### 2.3 为什么差异是"系统性的"而非噪声

ImQMD 的 NN 碰撞动力学是**确定性混沌系统**：

1. `rng.shuffle(order)` 打乱碰撞对的遍历顺序。碰撞顺序改变意味着**完全不同的碰撞对**在 `used` 集合中被排除。
2. `rng.random() < p_block` 的泡利阻塞判定是二元的——每个碰撞要么被接受，要么被拒绝。一个判定不同，后续所有核子的动量分布就不同。
3. `rng.normal(size=3)` 生成散射的最终态动量方向——这是三维连续随机量，改变后核子的轨迹完全不同。

**一次泡利阻塞判定的差异 → 单次碰撞接受/拒绝不同 → 动量分布演化分叉 → 经过数千步传播 → 完全不同的碎片模式**。这正是我们在服务器上观察到的"完全核解体"——RNG 序列不同导致碰撞统计偏差，参数组合（是在笔记本 numpy 2.x 序列上调出来的）碰巧不适用于 numpy 1.x 序列。

### 2.4 参数调优依赖 RNG 序列的恶性循环

当前参数（`fermi_constraint_check` 阈值 170.0、`initialize_and_relax` 的 `relax_time=800.0`、`SkyrmeEDF` 表面力 `gsur` 等）是在**笔记本电脑（numpy 2.4.6）**上调优的：

- **初始化阶段**（`initializer.py:100-113`）：用 numpy 2.x 的 RNG 序列生成初态位置和动量，经过 50 步 RK4 + 阻尼弛豫（`initializer.py:123-134`）。
- **弛豫阶段**（`initializer.py:200-201`）：`relax_once()` 在每一步调用 `fermi_constraint_check(nucleus)`（`initializer.py:179`），该函数使用 `_collision_rng` 的 `rng.normal()`。
- **传播阶段**（`propagator.py:273-274`）：`attempt_nn_collision()` 使用 `_collision_rng` 的 `shuffle()`, `normal()`, `random()`。

整个过程从初始化到传播全部依赖 RNG 序列。参数调优结果是**隐式地适配了 numpy 2.x PCG64DXSM 序列**，在 numpy 1.x PCG64 序列上当然失效。

## 三、修复方案

### 方案 A（推荐）：显式锁定 BitGenerator + 参数重调优

**核心思路**：使用 `np.random.PCG64(seed)` 显式指定 BitGenerator，确保 numpy 1.17+ 和 numpy 2.x 产生相同序列。然后在此固定 RNG 上重调参数。

**修改点**：

1. **`collisions.py` 第 18 行**：
   ```python
   # 修改前
   nucleus._collision_rng = np.random.default_rng(73129)
   # 修改后
   nucleus._collision_rng = np.random.Generator(np.random.PCG64(73129))
   ```

2. **`initializer.py` 第 92 行**：
   ```python
   # 修改前
   rng = np.random.default_rng(seed if seed is not None else 1000 + 17 * Z + A)
   # 修改后
   rng = np.random.Generator(np.random.PCG64(seed if seed is not None else 1000 + 17 * Z + A))
   ```

**优点**：
- 最小改动（2 行代码）
- `PCG64` 在 numpy 1.17+（服务器 1.25.2 满足）和 numpy 2.x 上均可用且行为一致
- 利用 numpy 的 `Generator` API，保留所有现有调用语法（`rng.shuffle()`, `rng.random()`, `rng.normal()`）

**缺点**：
- 需要**重调优物理参数**：锁定的 PCG64 序列与当前笔记本序列（PCG64DXSM）不同，也与服务器序列（PCG64）不同。所有参数（费米阈值、弛豫时间、表面力等）需在 3 套 RNG 中选一个重新校准
- 依赖 numpy 的 PCG64 实现细节在未来版本中不被"改进"

**验证方法**：
```python
import numpy as np
# 在两台机器上分别运行
rng = np.random.Generator(np.random.PCG64(73129))
assert rng.random() == expected_value  # 应输出相同值
```

### 方案 B（最稳健）：内嵌独立确定性 RNG

**核心思路**：实现一个不依赖任何外部库的确定性 PRNG，直接从 Python 内置模块构建。

**实现**：利用 Python `hashlib` + `struct` 构建基于 SHAKE256 的 CSPRNG（或简化为 splitmix64 等轻量算法），确保跨平台跨版本 100% 一致。

**伪代码骨架**：
```python
import hashlib
import struct

class StableRNG:
    """Deterministic RNG using SHAKE-256, version-independent."""
    def __init__(self, seed: int):
        self._state = hashlib.shake_256(seed.to_bytes(8, 'little'))
    
    def _next_u64(self) -> int:
        block = self._state.digest(8)
        self._state = hashlib.shake_256(block)  # 反馈模式
        return struct.unpack('<Q', block)[0]
    
    def random(self) -> float:
        return self._next_u64() / float(2**64)
    
    def normal(self, size=None) -> np.ndarray:
        # Box-Muller transform
        ...
    
    def shuffle(self, arr) -> None:
        # Fisher-Yates
        ...
```

**优点**：
- **完全独立于 numpy 版本**——任何 Python 3.8+ 环境都产生相同序列
- 不受 numpy 未来版本更改的影响
- 可移植到任何语言/平台

**缺点**：
- 实现量大（需要包装为 `np.random.Generator` 兼容接口或替换所有调用点）
- 性能：纯 Python 的 RNG 比 C 实现的 PCG64 慢约 10-100×，但碰撞采样不是主要性能瓶颈

### 方案 C（快速但不彻底）：legacy RandomState

**核心思路**：使用 `np.random.RandomState(seed)`（Mersenne Twister），该 API 自 numpy 1.x 以来行为稳定。

**修改点**：
- `collisions.py` 和 `initializer.py` 将 `np.random.default_rng(seed)` 改为 `np.random.RandomState(seed)`

**优点**：
- 一行代码改动
- Mersenne Twister 在 numpy 全版本中行为一致

**缺点**：
- `RandomState` 是 legacy API，未来可能被废弃
- Mersenne Twister 统计质量不如 PCG 系列（但不影响物理正确性）
- 同样需要参数重调优

### 方案 D（接受现状）：文档化版本依赖

记录需求的 numpy 版本，在服务器上升级或降级 numpy。

**优点**：零代码改动

**缺点**：
- 没有真正解决问题——如果笔记本坏掉或升级 numpy 到 3.x，代码又坏了
- 服务器管理员可能不允许升级 numpy

## 四、推荐实施路径

**分三步走**：

### 第一步：确认 RNG 差异（立即）

在服务器和笔记本上分别运行验证脚本，确认 `default_rng(73129)` 输出不同：

```python
import numpy as np
rng = np.random.default_rng(73129)
print(f"NumPy {np.__version__}: first 5 randoms = {[rng.random() for _ in range(5)]}")
```

### 第二步：实施方案 A（本次迭代）

将两处 `default_rng()` 改为显式 `np.random.Generator(np.random.PCG64(seed))`。

**需要同步调整的代码位置汇总**：

| 文件 | 行号 | 当前代码 | 目标代码 |
|------|------|---------|---------|
| `collisions.py` | 18 | `np.random.default_rng(73129)` | `np.random.Generator(np.random.PCG64(73129))` |
| `initializer.py` | 92 | `np.random.default_rng(seed if seed is not None else 1000 + 17 * Z + A)` | `np.random.Generator(np.random.PCG64(seed if seed is not None else 1000 + 17 * Z + A))` |

### 第三步：参数重调优（后续）

在锁定 RNG 后，重新标定物理参数：

1. **费米约束阈值**（`collisions.py:380`，当前 `threshold=170.0`）：调整 `phase_space_threshold` 使初态费米相空间距离达标。
2. **弛豫时间**（`initializer.py:158`，当前 `relax_time=800.0`）：调整使 `relative_energy_drift` < 0.08。
3. **表面力参数**：如需要，微调 `SkyrmeParameters` 中的 `gsur`。
4. **阻尼系数**（`initializer.py:133,167,180`，当前 `0.93`）：可能需要小幅调整。

## 五、附录：RNG 调用链全景图

```
initialize_nucleus(seed=None → 1000+17*Z+A)
  │
  ├── np.random.default_rng(seed)          ← [修改点 2]
  │   ├── rng.uniform() × N           ← _sample_hard_sphere()
  │   ├── rng.shuffle()               ← 质子打乱
  │   ├── rng.normal() × N            ← _fermi_momenta() 方向
  │   └── rng.random() × N            ← _fermi_momenta() 大小
  │
  ├── fermi_constraint_check()
  │   └── _rng() → _collision_rng
  │       └── rng.normal()            ← 动量修正方向
  │
  └── _rk4_step() × 50 (弛豫循环)
      └── (纯确定性, 不用 RNG)

propagate()
  │
  ├── _rk4_step() × n_steps (纯确定性)
  │
  ├── attempt_nn_collision()
  │   └── _rng() → _collision_rng       ← [修改点 1]
  │       ├── rng.shuffle(order)     ← 碰撞对顺序
  │       ├── rng.normal(size=3)     ← _scatter_isotropic()
  │       └── rng.random()           ← Pauli 阻塞判定
  │
  └── fermi_constraint_check()
      └── _rng() → _collision_rng
          └── rng.normal()            ← 动量修正方向
```

## 六、总结

| 项目 | 内容 |
|------|------|
| **根因** | `np.random.default_rng()` 在 NumPy 1.x 使用 PCG64、在 NumPy 2.x 使用 PCG64DXSM，相同种子产生不同随机序列，导致混沌碰撞动力学的系统性分叉 |
| **影响范围** | `collisions.py` 1 处 RNG 创建 + `initializer.py` 1 处 RNG 创建，共 7 处调用点 |
| **推荐方案** | 方案 A：显式 `np.random.Generator(np.random.PCG64(seed))` + 参数重调优 |
| **最小改动** | 2 行代码（`collisions.py:18` 和 `initializer.py:92`） |
| **后续工作** | 费米阈值、弛豫时间、表面力等参数需在锁定 RNG 后重新标定 |
