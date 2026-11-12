> **2026-10-06 最终状态**：本文是第一轮复查结论，部分状态已过时（生产源码已安装；数据库重启问题已由 Codex 修复；
> 第二轮又修复一处三维迁移错误并完成三维全路径分类、13 面墙验算）。当前权威说明见
> [`PROGRAM_DESCRIPTION.md`](PROGRAM_DESCRIPTION.md)，Codex 的后续开发见 `DEVELOPMENT_REPORT.txt`。

# ADINA 混凝土迁移复查与实用化修改（2026-10-06）

详细程序说明见 `PROGRAM_DESCRIPTION.md`，源码差异见 `CHANGES.diff`。

本目录是对 `LLM_HANDOFF.txt` 所述迁移的独立复查结果、修改后的源码和可重复的验证工具。
结论只覆盖下面实际运行过的内容。

## 1. 现状与交付方式

- 生产源码 `SRC/material/nD/AdinaConcrete*`、`adinaConcrete2D/`、`adinaConcrete3D/`
  属于另一个 Windows 账户（CodexSandbox），本账户无写权限。修改后的完整源码在 `work/`，
  修改前快照及哈希在 `src-before/`。生产源码**未被改动**。
- 安装：放开写权限后运行 `python install_sources.py`（只检查）和
  `python install_sources.py --apply`（备份、复制、在 CMakeLists 中登记新头文件），再重新构建 bm。
  若生产文件在此期间被他人改过，脚本拒绝覆盖。
- 未安装前可用 `build_variant.py --out <目录>` 以 bm 的编译参数单独链接出一个
  OpenSees.exe（不触碰 bm）。本报告的最终可执行文件：
  `tmp/claude-review/v10/OpenSees.exe`，SHA256 `93d6747d4380040e3711bc3e4ae65c34efacd154c6a7942290d5df7ba5b6344f`
  （与 v9 源码只差行尾符，sw1-1 结果逐位相同），
  引擎库 libOpenSeesLIB.a SHA256 `b864d1be…312ed6`。

## 2. 迁移是否正确：与原 CDMOD 逐步整路径对比

工具：`difftest/`。直接从只读 `C:/Work/Adina/all_adina.for`（SHA256 `9a0e2731…ffb1b`）
抽取 CDMOD/PRNCPL/CRAKID/DCRACK 原文编译，与 C++ 内核在 44 条完整加载路径 × 2 套参数上逐步比较
（开裂、闭合、再张开、剪切、双向开裂、压碎、软化、压碎后反向、随机非比例路径）。

与此前审计的两点不同：
1. 原驱动 `tests/materials/adina_2d/oracle_driver.f90` 设 `ieqref=1`，此时 CDMOD 只返回弹性矩阵（标号 85）；
   真正的迭代矩阵在 `IEQREF=0`（标号 74→5→100→110），此前从未检验。本工具用 `IEQREF=0`。
2. 此前整步对比在开裂前（第 49 步）截断；本工具比较全路径。

参考程序有两份：`raw`（原文），`corrected`（仅三处已记录的一行修正：CRAKID 工程剪应变少乘 2、
PRNCPL 将新局部剪应变置零、`57.296` 度换算）。

结果（`difftest/report_before.json` 为修改前内核，`report_v7.json` 为修改后）：
- 大多数路径与 `corrected` 吻合到 1e-14·|fc|，说明迁移整体正确。
- **发现并修复的迁移错误（二维、三维都有）**：压碎后原程序在最小主应力不再受压
  （`SIGP ≥ 0.001·SIGCP`）时取 E=0（应力冻结），C++ 缺此判断，继续以负模量把应力推向受拉。
  `biaxial_compress_crush` 误差由 6%·fc 降到 1e-6·fc。
- **二维修复**：压碎后应力全释放的判断，原程序用“更新后应力”的主轴应变并检查厚度应变 TEP(4)，
  C++ 用旧主轴且漏检厚度应变。修复后压碎路径与原程序逐步一致。
- 其余差异均已定位，不是迁移错误：
  - 恰好等于 ft 的平局（路径应变恰为开裂应变，舍入决定先后一步）；
  - 原 CRAKID 缺陷：两条裂缝中一条闭合时另一条也被丢失；闭合裂缝在垂直方向重新开裂时释放了错误平面的应力。
    C++ 的处理在物理上正确，差异 0.2–0.3·fc 属有意保留；
  - ≤1e-7·fc 的舍入累积。
- 三维的 CMOD3D 尚未做同样的全路径对比（只用单元/结构算例做了运行检查）。

## 3. 为什么构件算不动：原因与修改

| 问题 | 证据 | 修改 |
|---|---|---|
| 每次 Newton 试算都从上一收敛态重放，脆性开裂（ft 处应力突降）在迭代间出现又消失，平衡点落在突降区间内时无解 | sw1-1 第一个侧移步 100 次迭代，残差在 110–280 N（≈ft×一个积分点层面积）间循环 | 新增 `-update iterative\|committed`，**默认 iterative**：与 ADINA 本身一致。ADINA 的单元循环在每次平衡迭代中原位更新 WA 历史（ECHECK 置 IUPDT=0，CDMOD 覆写 SIG/EPS/EVMAX/ANGLE/CRKSTR/PGRAV） |
| 有限差分“算法切线”每列一次完整重放，开销大，且在脆性跳变处无意义 | 单步耗时 2.6 s → 0.4 s | 新增 `-tangent secant\|algorithmic`，**默认 secant**：即 CDMOD/CMOD3D 在 IEQREF=0 时重建的材料矩阵（加载曲线切线、卸载 E、开裂 STIFAC/SHEFAC、压碎 E·STIFAC）。与原 CDMOD 返回矩阵对比吻合到 4–5 位有效数字。断裂能模式下开裂面用黏聚割线 |
| 原 ADINA 压碎后卸载沿软化曲线“爬回”，越过峰值应变后模量为 0：**应力锁定在峰值压应力，应变变成受拉也不变** | 原 CDMOD 实测（`crush_reverse` 路径）；墙 sw1-1：+10 mm 时 5 kN、−10 mm 时 168 kN | 新增可选扩展 `-crushedUnloading elastic`（**默认 original**，不改原模式）：低于历史最大压应变时按 E 卸载/再加载，主应力不受拉、压应力不超过压碎曲线在历史应变处的值；超过历史压应变及极限应变释放仍按原律 |

其他：
- 序列化升级为 v3（2D 118 项、3D 198 项），保存上述选项、压碎历史应变和二维释放标志；
  v1/v2 旧检查点按其产生时的公式恢复（committed + algorithmic + original）。
- 迭代模式下，同一应变的重复请求不再施加零增量（2D 原有、3D 新增）。
- `-DADINA_TRACE` 编译时可输出试算间状态变化，仅用于诊断构建。
- 试过“锁定本步已出现裂缝”的折中方案（latched），sw1-1 在裂缝闭合/再张开上仍两步循环，已删除。
- 试过在压碎下降段让迭代矩阵取负切线，墙 sw1-3 首步即失败，已撤回。

旧行为仍可用：`-tangent algorithmic -update committed`。

## 4. 构件结果（与实验骨架点对比）

驱动 `wall_driver_adina.tcl`：与原墙驱动相同，只在 `newton` 档先用完整 Newton；力容差、位移协议、
加载数据均未改。步长 dt=0.02（0.2 mm）。实验点取自冻结的 `reference_summary_csv`。

此前 22 次构件运行均未走完；修改后：

| 案例 | 原模式（默认） | `-crushedUnloading elastic` | PlainConcrete（原模型，同驱动） |
|---|---|---|---|
| sw1-1 | 到 t=48.7 停滞；平均误差 30.1% | **走完**，平均误差 17.8%；7.17 mm：167.9 vs 160.3 kN | 到 t=45.1 停滞，4.4% |
| sw1-3 | **走完**，46.2% | **走完**，40.6%（≤8 mm 时 ≤5%） | t=12.2 停滞 |
| sw2-2 | **走完**，但负位移处力的符号错误（锁定应力所致） | **走完**，48.4% | t=45.8 停滞，43.9% |
| 2H16-Ref 梁 | **走完**（此前停滞/超时） | 同左 | — |

结论：
- 原 ADINA 模式在往复荷载下受“压碎应力锁定”缺陷影响，结果不可用于循环分析；
  `-crushedUnloading elastic` 消除了该问题，峰值前和峰值附近与实验相差约 5–10%。
- 15 mm 以后普遍低估 40–60%。原因是 ADINA 模型在 εu 处释放全部应力，墙脚压碎集中在单排
  0.2 m 单元内，且约束区混凝土仍用未约束参数（εu=−0.006）。这是模型参数和单元尺寸问题，
  不是程序错误；如需改善，应提供约束混凝土参数，或另行研究受压能量正则化（未实现）。
- 迭代原位更新的代价是结果依赖迭代路径。sw1-1 在 dt=0.05 时早期循环与 dt=0.02/0.01 差约 20–50%；
  dt≤0.02 时 12 mm 以内大致收敛（约 10%）。建议使用小步长并做步长敏感性检查。
- 梁 2H16-Ref：ShellMITC4 与 ShellNLDKGQ 的初始刚度都约为实验的 3 倍，而弹性梁理论
  （跨度 1.5 m、截面 0.16×0.34 m、E=30 GPa）给出与有限元相同的初始刚度，因此差异来自试验装置
  （支座/垫层柔度、位移测点）或假定的 E，不是单元问题。开裂后 NLDKGQ 比 MITC4 好
  （超估约 30% 对 46%），建议在面内弯曲的粗网格模型中优先用 ShellNLDKGQ。

## 5. 三维

- 单个 stdBrick 单轴（侧向自由）加载：默认设置在压碎起始处不收敛。原因是压碎后峰值应变前
  原律模量为 0，第一次迭代预测应力中的残余侧向应力无法再被平衡；`committed + algorithmic` 可以走完。
- 81 个 stdBrick 的素混凝土柱，轴压加往复侧移：默认设置推进到 t=4.1（旧默认只到 t=1.85），
  单步约快一倍；两者都未走完全程（素混凝土柱在约 0.9% 侧移时可能本已接近破坏，未作物理判别）。
- 三维尚无实验对比。

## 6. 其他发现（未修改）

- 数据库重启（`database File` / `save` / `restore`）对 ShellMITC4 + LayeredShell + 非线性材料
  不能复现不中断的结果。这一点与本次修改无关：修改前的可执行文件、OpenSees 自带的
  PlasticDamageConcretePlaneStress 也一样；弹性材料可以复现。恢复后材料点的应力、应变和裂缝状态
  逐位相同，但单元抗力不同，问题在单元/截面层（ShellMITC4 含其他工作的修改）。
  见 `restart/restart_test.tcl`。此前报告中“FileDatastore 生命周期通过”的结论不成立。
- 迭代失败回滚后，断裂能模式的应力与回滚前差 1 ulp（提交后同应变再算一次所致），无实际影响。
- `tests/concreteValidationADINA2D` 中 sw2-2 实验骨架点与 PlainConcrete 结果相差 36–60%，
  sw2-2 实验资料本身需要核实。

## 7. 使用建议

```tcl
# 墙、梁等往复加载（推荐）
nDMaterial AdinaConcrete2DPlateFiber 1 $E $nu $ft $fc $epsc $fu $epsu -crushedUnloading elastic
# 完全按原 ADINA（单调加载、与原程序对照）
nDMaterial AdinaConcrete2DPlateFiber 1 $E $nu $ft $fc $epsc $fu $epsu
```
求解：先用完整 `Newton`（KrylovNewton 的外推在原位更新下可能锁定虚假损伤）；
位移增量取约 0.2 mm 或更小，并做一次步长减半检查。

## 8. 重复验证的命令

```
python difftest/difftest.py --out difftest/out                               # 与原 CDMOD 全路径对比
python build_variant.py --out v-check                                        # 链接试验可执行文件
python launch_wall.py r-sw1-1 v-check/OpenSees.exe --dt 0.02 --knot 0.02 --material-extra "-crushedUnloading elastic"
python compare_exp.py sw1-1 r-sw1-1
OpenSees restart/restart_test.tcl straight|save|restore|rollback [选项]
OpenSees brick3d/column_test.tcl [选项]
```
本次的全部运行目录（含 stdout、status、curve）保留在 `<合并前工作区>/tmp/claude-review/`
（r14–r18、run3–run13、r15/r16 梁）。
