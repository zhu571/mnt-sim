# 审查报告 — t_a4414f67: U+U MST/iso-MST 碎片识别 + E* 异常 + 碎片过多

**审查对象**: commit `0add78a` (`feat(imqmd): add reaction cross-section framework, fix U+U transfer freeze`)
**审查文件**: `mnt_sim/imqmd/fragments.py`, `propagator.py`, `reaction.py`, `grid_edf.py`, `initializer.py`, `collisions.py`
**审查时间**: 2026-07-04

---

## 1.  empirical 验证

使用仓库自带脚本与默认 IQ3A 参数集，在 b=4.0 fm、E/A=7.0 MeV、t=600 fm/c、relax_time=200 fm/c 条件下运行单个 U+U 事件：

```text
[iso-mst default] 47 fragments
  Z= 43 A=107 E*=   283.2 MeV
  Z= 24 A= 60 E*=   175.4 MeV
  ...
```

近垒 U+U 碰撞产生了 **47 个碎片**，最重碎片仅 A=107。这属于完全多重碎裂（multifragmentation），与近垒反应应以准弹性/深非弹性/少量核子转移为主的物理图像严重不符，也直接 contradict 提交信息所称的 "transfer products Z=93(Np), Z=91(Pa), Z=80-81"。

对冷 U-238 核的静态测试显示，mst/iso-mst/adaptive_mst 都能正确识别为 1 个碎片；因此问题不在碎片识别算法本身，而在 **碰撞传播后的位形已经被打得过碎**。

---

## 2.  🔴 未达标项

### 2.1 碰撞传播使用未充分调参的 GridEDF 力

`reaction.py:268` 在 `run_imqmd_event` 中硬编码：

```python
use_grid_edf=True,
```

README 明确写明：

> GridEDF | 271 行网格积分替代质心采样（**待参数调优后启用**）

在未完成参数调优前把 GridEDF 力用于重离子碰撞动力学，是导致 U+U 过碎的最可能根因之一。GridEDF 力与质心采样力还存在已知的 Coulomb exchange 归一化不一致（见 3.2）。

### 2.2 multiprocessing 在 Windows 上会崩溃

`reaction.py:373`：

```python
ctx = mp.get_context("fork")
```

Windows 平台不支持 `fork` 启动方式（仅支持 `spawn`）。用户环境为 Windows，设置 `n_workers>1` 会直接抛出 `ValueError`。

### 2.3 Fermi 约束在碰撞期间被完全关闭

`propagator.py:292-294`：

```python
if apply_fermi_constraint is None:
    apply_fermi_constraint = not with_collisions
```

虽然这确实解除了原先阻止核子穿过弹靶边界的 Pauli 阻塞，但**完全关闭** Fermi 约束会导致碰撞期间允许违反 Pauli 原理的相空间态，进而产生非物理的高激发/高动量核子，加剧碎裂。更合理的修复应是按 `reference_group_ids` 分组施加约束或放宽阈值，而不是整段碰撞期间禁用。

---

## 3.  🟡 偏差项

### 3.1 默认 `fragment_method` 不统一

- `reaction.py:329` / `run_imqmd_event:230` / `impact_parameter_scan:329` 默认 `"iso-mst"`
- `fragments.py:167` 的 `reaction_fragments()` 使用 `adaptive_mst`
- `reaction.py:identify_fragments()` 只支持 `"iso-mst"` / `"mst"`，不支持 `"adaptive-mst"`

README 列出的碎片方法是 MST/iso-MST，而实际反应截面扫描入口却默认 iso-MST，与 `reaction_fragments` 不一致，易造成 API 混淆。

### 3.2 GridEDF 与 SkyrmeEDF 的 Coulomb exchange 归一化不一致

`skyrme.py:112`（质心采样）：

```python
coeff * np.sum(rho_p ** (4.0 / 3.0)) / self.parameters.rho0
```

`grid_edf.py:152`（网格积分）：

```python
return -0.75 * E2 * (3.0 / np.pi) ** (1.0 / 3.0) * rho_p ** (4.0 / 3.0)
```

网格版缺少 `/rho0` 因子（约 6 倍差异）。虽然总能量可通过 `nuclear_scale` 与 `energy_offset` 整体校准，但力场细节不一致会改变碰撞动力学。

### 3.3 `cold_ground_state_energy` 的标度因子偏离物理 1.0

对 U-238 的 `grid_energy_diagnostics` 显示：

```text
raw_total=-3258 MeV, target_total=-1815 MeV, nuclear_scale=0.659
```

虽然 `cold_ground_state_energy` 通过 `energy_offset` 把总能校准到经验值，但核力标度只有 0.66，说明 GridEDF 尚未真正复现 IQ3A 的饱和性质。以此为基础计算碎片激发能 E* 时，对重核（Z~92）的可靠性显著低于中等质量核（Z~80-83）。

### 3.4 Fermi 阈值与 README 不一致

README 写 "Fermi = 170"，但 `collisions.py:367` 默认 `threshold=255.0`，且传播器没有把 README 中的 170 传入。参数文档与实现不一致。

---

## 4.  🟢 建议项

1. **在 GridEDF 完成参数调优前**，把 `run_imqmd_event` 的 `use_grid_edf` 默认改为 `False`，或至少暴露为可调参数。
2. **统一碎片方法**：在 `identify_fragments` 中增加 `"adaptive-mst"`，并明确 `impact_parameter_scan` 的默认方法；若研究转移则默认 `"mst"` 更合理。
3. **Windows 兼容性**：把 `mp.get_context("fork")` 改为根据平台自动选择 `spawn`/`fork`。
4. **修复 Coulomb exchange 归一化**，使 GridEDF 与 SkyrmeEDF 一致。
5. **Fermi 约束分段/分组启用**：仅在碰撞剧烈阶段放宽，或按 `reference_group_ids` 在组内维持 Pauli 抑制。
6. 增加物理合理性检查：例如 U+U 近垒事件的最重碎片 A 应 >150，否则报错或提示。

---

## 5.  综合结论

**❌ 需修改**

commit `0add78a` 虽然搭建了反应截面扫描框架并解决了 "transfer freeze" 的表象问题，但当前默认参数下 U+U 近垒事件产生完全多重碎裂，与物理预期严重不符。根因集中在：

- 过早启用未调优的 GridEDF 碰撞力；
- Fermi 约束整段关闭导致 Pauli 相空间失控；
- 部分实现细节（Coulomb exchange、Windows 多进程、默认碎片方法）尚未对齐。

建议在修复上述问题并重新通过 U+U 物理合理性检查前，不要将该版本用于截面扫描生产。
