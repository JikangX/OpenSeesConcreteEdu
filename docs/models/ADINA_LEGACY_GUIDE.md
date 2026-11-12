# OpenSees 中的 ADINA 混凝土本构（模式 5）——使用、开发与验证说明

> 目录迁移说明：本文保留原开发阶段的描述，其中未随主线提供的 `results/`、Audit、Walls等路径属于原版本。当前使用入口见 [ADINA二维与三维说明](ADINA_CONCRETE_CN.md)；示例位于 `EXAMPLES/AdinaConcrete/`，当前测试位于 `tests/materials/adina_2d/`、`adina_3d/`。

> 2026-10-08 文档索引：本目录原始说明保留。当前参数、版本及验证边界汇总于
> [统一 ADINA 说明](ADINA_CONCRETE_CN.md)，全部材料与壳接口见
> 模型总览（OpenSeesCUDA 仓库 `docs/models/README.md`）。该句属于合并前文档核查；本次整合后的验证见 整合记录（合并前记录：`../docs/integration/VALIDATION_CN.md`）。

本目录是 ADINA 混凝土模型（ADINA 材料模式 5，Bathe–Ramaswamy 型固定裂缝/压碎模型）迁移到 OpenSees 的
**入口文档**。面向两类读者：

- 要在自己的项目中**使用**该本构的人：看第 1–4 节（命令、参数、推荐设置、例子、限制）。
- 要**检查或继续开发**的人：看第 5–9 节（迁移过程、修改清单、验证方法与结果、复现命令、文件索引）。

最详细的技术说明（逐条修改依据、存盘格式、全部数值）在
[`tests/materials/adina_review/PROGRAM_DESCRIPTION.md`](../../tests/materials/adina_review/PROGRAM_DESCRIPTION.md)，
其中 **R4 节为最新优选案例复核，R3 节描述当前材料修订与默认值；与旧章节冲突时以 R4/R3 为准**。

---

## 0. 一页摘要

| 项目 | 内容 |
|---|---|
| 材料命令 | `AdinaConcrete2D`（平面应力单元，如 `quad`）、`AdinaConcrete2DPlateFiber`（分层壳 `LayeredShell` 的混凝土层）、`AdinaConcrete3D`（实体单元，如 `stdBrick`） |
| 对应 ADINA 例程 | 2D：A06.FOR 的 CDMOD / PRNCPL / CRAKID / DCRACK；3D：A16.FOR 的 CMOD3D / PRNCP3 / CRAK3D / DCRAK3 |
| 默认设置 | **与原 ADINA 程序一致**：每次平衡迭代从上一收敛态积分（`-update committed`）、每个增量调用一次本构（`-maxStrainIncrement 0`）、迭代矩阵为 ADINA 材料矩阵（`-tangent secant`）、压碎后按原律（`-crushedUnloading original`） |
| 正确性证据 | ① 逐子程序与原 Fortran 例程逐步比对（2D 88 组路径、3D 92 组路径，全部差异均归因于原程序已识别的缺陷或阈值平局）；② **把原 ADINA 程序完整编译运行**，在单元与构件层面对同一有限元模型逐步比对：单元层差 < 30 N（峰值 3e7 N），2H20-Ref 梁全程差 ≤ 0.022 kN（对修正已知缺陷后的 ADINA），W28 墙 3 mm 前差 ≤ 2.3% |
| 主要局限（属 ADINA 模型本身） | 剪切控制构件强度高估（无受拉软化、固定裂缝常数剪力保留 0.5 → 应力锁定）；脆性开裂使 Newton 在开裂/闭合间循环、往复加载反向后难以收敛；压碎后单点集中力无法继续承载。原 ADINA 程序在相同模型上表现相同 |
| 可选的非 ADINA 扩展 | `-update iterative|hybrid`、`-crushedUnloading elastic`、`-fractureEnergy Gf Lch`（裂缝带受拉软化）。默认关闭；使用时结果**不再是 ADINA 结果** |
| 代码位置 | `SRC/material/nD/adinaConcrete2D/`、`SRC/material/nD/adinaConcrete3D/`、`SRC/material/nD/AdinaConcrete*.{h,cpp}`、`AdinaConcreteOptions.h` |

---

## 1. 命令接口

```tcl
nDMaterial AdinaConcrete2D           $tag $E $nu $ft $fc $epsc $fu $epsu  <options>
nDMaterial AdinaConcrete2DPlateFiber $tag $E $nu $ft $fc $epsc $fu $epsu  <options>   ;# LayeredShell 层
nDMaterial AdinaConcrete3D           $tag $E $nu $ft $fc $epsc $fu $epsu  <options>
```

OpenSeesPy 中同名：`ops.nDMaterial('AdinaConcrete2D', tag, E, nu, ft, fc, epsc, fu, epsu, *options)`。
**材料名区分大小写**（`AdinaConcrete2D`，首字母大写）。

### 1.1 必需参数（单位自洽；受拉为正，**受压量输入为负**）

| 参数 | ADINA 卡片名 | 含义 |
|---|---|---|
| `E` | E0 | 零应变处初始切线模量 |
| `nu` | ν | 泊松比（全过程常数） |
| `ft` | SIGMAT | 单轴抗拉强度（开裂应力，> 0）。注意：试验报告中的抗折强度（modulus of rupture）不是单轴抗拉强度 |
| `fc` | SIGMAC | 单轴抗压强度（< 0） |
| `epsc` | EPSC | 峰值压应变（< 0） |
| `fu` | SIGMAU | 极限压应力（< 0，`fc ≤ fu < 0`） |
| `epsu` | EPSU | 极限压应变（< 0，`epsu < epsc`），超过后压碎点应力全部释放 |

单轴受压曲线为 ADINA 有理式 σ = E ε / (1 + A x + B x² + C x³)，x = ε/εc，A、B、C 由 E、fc、εc、fu、εu 确定（手册 p.105）。
要求 `fc/epsc ≤ E`。

### 1.2 选项

| 选项 | 默认 | 是否 ADINA | 说明 |
|---|---|---|---|
| `-retention STIFAC SHEFAC` | `1e-4 0.5` | 是 | 开裂面法向/剪切刚度折减系数（手册卡片 5） |
| `-controls BETA GAMA KAPA ALFA` | `0.5 1 0.7 0` | 是 | 破坏包络 β、临界应变放大 γ、正交异性切换 κ、加载函数 α |
| `-envelope 24 个数` | 手册默认 | 是 | 三轴破坏包络 SP1(6)、SP3(6,1..3)（均除以 fc，正值） |
| `-G G` | E/(2(1+ν)) | （仅 PlateFiber）| 壳层横向剪切模量（横向剪切弹性） |
| `-rho rho` | 0 | — | 质量密度 |
| `-tangent secant\|algorithmic` | secant | secant = ADINA | secant：ADINA 在 IEQREF=0 时重建的材料矩阵；algorithmic：对整个试算的有限差分 Jacobian |
| `-update committed\|iterative\|hybrid` | **committed** | committed = ADINA | committed：每次试算从上一收敛态积分（原程序 UNBLD 的做法）；iterative：从上一次迭代状态积分；hybrid：iterative 但以收敛态加载函数最大值作为加卸载参照 |
| `-maxStrainIncrement x` | **0** | 0 = ADINA | 0：一次更新（原程序每增量一次 CDMOD/CMOD3D）；> 0：按应变张量范数 x 细分子步 |
| `-crushedUnloading original\|elastic` | original | original = ADINA | elastic：压碎点在历史最大压应变以内按 E 卸载/再加载（非 ADINA 扩展） |
| `-fractureEnergy Gf Lch` | 关闭 | 否 | 裂缝带线性受拉软化（非 ADINA 扩展）；开启后默认 algorithmic + committed |

### 1.3 输出（`eleResponse ... material i [fiber j] <name>`）

| 名称 | 2D（10 个值） | 3D（9 个值） |
|---|---|---|
| `stress` / `strain` / `tangent` | 标准 NDMaterial 输出 | 同 |
| `crackState` | open1, open2, crushed, maximumLoading, hasCrackFrame, angleRadians, crackStrain1, crackStrain2, crushingStress, crushingStrain | openCracks, crushed, maximumLoading, hasCrackFrame, crackStrain1..3, crushingStress, crushingStrain |
| `crackFrame` | — | 裂缝坐标 3×3（q11..q33） |
| `fractureState` | 断裂能历史（开启 `-fractureEnergy` 时） | 同 |

存盘/重启（`database File`）支持，格式 v3（2D 118 项、3D 198 项），兼容旧记录。

---

## 2. 推荐用法

1. **要得到 ADINA 的结果：用默认选项。** 不要加 `-update iterative`、`-maxStrainIncrement`。
2. **求解器**：完整 Newton（`algorithm Newton`），力平衡判据 `test NormUnbalance`；失败时依次尝试
   `NewtonLineSearch`、`KrylovNewton`，最后 `KrylovNewton -iterate initial -increment initial`（初始刚度 Krylov 加速，
   允许较多迭代，如 2000），再失败才二分步长。这些都只改迭代策略，不放宽容差。
3. **位移加载请用 `DisplacementControl`**（在控制节点施加单位参考荷载，多个加载点用 `equalDOF` 绑定）。
   `LoadControl` + `sp` 时 OpenSees 在首次求解前把 sp 值写入节点，单元先在“只有加载节点移动”的非平衡状态下被计算；
   对 committed 默认这只影响迭代路径，但对非 ADINA 的 iterative/hybrid 会留下永久损伤。
   与原 ADINA 逐步对照时，OpenSees 侧用 `constraints Lagrange`（与 ADINA 首次迭代相同）。
4. **分层壳**：混凝土层厚度不宜过大（例如 ≤ 15 mm）。脆性开裂一次释放整层应力，2 层混凝土的板在首次开裂时无法收敛；
   细分为 10 层后可继续计算。
5. **集中荷载**：压碎后 ADINA 压碎律在峰值后不再承载，单节点集中轴压/集中力下的单元会先压碎。请用加载板/加载梁分布荷载。
6. **单个均匀单元**：侧向自由的单个单元压碎后问题病态（ADINA 同样），材料检验请用多单元或约束侧向。
7. 断裂能、iterative、hybrid、elastic 压碎卸载只作研究用途，报告结果时须注明“非 ADINA”。

---

## 3. 最小例子（`EXAMPLES/AdinaConcrete/`，已用本分支程序运行）

| 文件 | 内容 | 本机结果 |
|---|---|---|
| `ex1_quad_beam.tcl` | 平面应力悬臂条（`quad` + 桁杆钢筋），DisplacementControl 推覆至 8 mm | 走完，8.05 mm 时 35.6 kN |
| `ex2_layered_shell.tcl` | `ShellNLDKGQ` + `LayeredShell`（`AdinaConcrete2DPlateFiber` + `PlateRebar`）墙推覆至 3 mm | 走完，3.05 mm 时 149.6 kN；示范 `crackState` 输出 |
| `ex3_brick.tcl` | 两个 `stdBrick` 素混凝土柱单轴受压至 ε = 0.002 | 峰值 30.17 MPa（输入 fc = 30 MPa） |

运行：`OpenSees ex1_quad_beam.tcl`（Windows 需 PATH 中有 MinGW 运行库）。

---

## 4. 已知局限（不是迁移错误）

经原 ADINA 程序在同一模型上确认（第 7.3 节）：

1. **剪切控制构件强度偏高**：2H20-Ref 梁 ADINA 105 kN / OpenSees 105 kN / 试验 66 kN；W28 墙 3 mm 时 ADINA 892 kN /
   OpenSees 912 kN / 试验峰值 290 kN。原因：无受拉软化、固定裂缝常数剪力保留（SHEFAC 0.5）造成应力锁定。
2. **脆性开裂导致不收敛**：B07 梁原 ADINA 在 1.69 mm 停止（试验峰值在 7.66 mm）；W28 往复加载时原 ADINA 在第一次反向、
   裂缝闭合后即不收敛。OpenSees 默认设置下也常在中途停止。
3. **集中力下压碎**：sw6-3 墙单节点 1155 kN 轴压，原 ADINA 只能加到 0.85–0.884，OpenSees 0.85–0.891。
4. 原程序自身的缺陷（PRNCPL 180° 角、CRAKID/CRAK3D 闭合后丢失或错放裂缝平面、DCA 方向余弦符号丢失、PRNCP3 近重根近似）
   **在 C++ 中按物理正确的方式实现**，因此在这些特定路径上与原程序不同；差异均已分类说明。

---

## 5. 迁移与开发过程

| 轮次 | 执行者 | 内容 |
|---|---|---|
| 初次迁移 | 前一 LLM | 把 CDMOD 等 Fortran 例程改写为 C++ 内核与 OpenSees 适配器；当时测试效果不好（`tests/adinaConcreteAudit/LLM_HANDOFF.txt`） |
| 第一轮复查 | Claude | 与原例程逐子程序对比；修复压碎后模量、二维压碎释放等迁移错误；增加 secant 迭代矩阵、存盘格式 v3 |
| 中间开发 | Codex | 安装进生产源码；7 项修复（三维压碎释放主轴、getCopy 的 dbTag、LayeredShellFiberSection 初始切线（影响重启复现）、开裂比较舍入邻域、断裂能默认等） |
| 第二轮 | Claude | 三维 CRAK3D 单裂缝包络（标号 111）迁移错误修复；三维 92 组全路径分类；13 面剪力墙验算；说明文档 |
| 第三轮 | Claude | ① 用离线试验库（25 个 HEOM 有限元模型）与 4 根三维实体梁调试；② `ShellNLDKGQ` 失败传递；③ **完整编译并运行原 ADINA 程序**，发现原程序每次迭代从收敛态积分，据此更正默认算法；④ 单元与构件层 ADINA–OpenSees 对比；⑤ 文档 |

### 5.1 第三轮的程序修改（相对第二轮，`tests/materials/adina_review/CHANGES_ROUND3.diff`）

| 编号 | 文件 | 内容 |
|---|---|---|
| C1 | `SRC/element/shell/ShellNLDKGQ.{h,cpp}` | 新增 `update()`；截面拒绝试算时返回失败，分析回退并可二分（此前失败被忽略，继续用旧应力组装，并打印上百万条警告） |
| C2 | `AdinaConcrete2D.cpp`、`AdinaConcrete3D.cpp` | 拒绝信息限频（前 20 条及每 100000 条） |
| C3 | `Concrete3D.h`、`Concrete2D/3D.cpp`、适配器、`AdinaConcreteOptions.h`、工厂 | 默认 `-update committed`、`-maxStrainIncrement 0`；新增 `-update hybrid`；加载判别写成 `f ≥ fmax`（原算法中与相等判别等价）；存盘枚举兼容 |
| C4 | 测试 | 新增默认值与 hybrid 检查；依赖子步的旧检查显式设置 `-maxStrainIncrement 2e-5` |

此前各轮的全部修改见 `PROGRAM_DESCRIPTION.md` 第 5 节（A1–A4、B1–B10）。本分支还包含同一工作区中与本构配套构建的其他修改
（LayeredShellFiberSection、ShellMITC4、Truss、PlateFromPlaneStress、PlateRebar、Domain/LoadPattern 重启修复、
NewtonLineSearch 线搜索失败传递、AbaqusSmearedCracking 与 PlaneStressUserMaterial 的 forumat.f90 等），
以保证本分支按原样可以构建出经过测试的程序。

---

## 6. 验证方法

1. **逐子程序对比**（`tests/adinaConcreteAudit/run_audit.py`）：7,877 条记录、26 组，0 失败。
2. **整路径逐步对比**（`tests/materials/adina_review/difftest/`）：编译原 CDMOD/CMOD3D 与“只修正已知缺陷”的对照版本，
   与 C++ 逐步比较；2D 88 组、3D 92 组（79,552 步）路径全部分类，无未解释差异；C++ 在所有路径上满足固定裂缝不变量。
3. **回归测试**（`tests/materials/adina_2d/run_tests.py`、`tests/materials/adina_3d/run_tests.py`）：内核、适配器、Tcl 壳/砖、
   失败传递、RC 膜片生命周期（3,961,660 项检查）、存盘重启（`tests/materials/adina_review/restart/`，8 组选项）。
4. **原 ADINA 程序整体对比**（`tests/materials/adina_review/adina84/`）：从只读源码编译完整程序，用生成的卡片运行同一模型。
5. **试验对比**：离线试验库（`cases2/`）、三维实体梁（`cases3d/`）、13 面墙（`run_all_walls.py`）。

---

## 7. 结果（第三轮生产版本）

### 7.1 回归

| 测试 | 结果 |
|---|---|
| 2D 回归 | 234,848 内核、6,755 适配器/截面、失败传递、3,961,660 RC 膜片检查，全部通过 |
| 3D 回归 | 37,358 内核、1,980 适配器、954 Tcl，全部通过 |
| 审计 | 7,877 条，26 组 0 失败 |
| 断裂能 | 62,736 内核、69,094 分层壳 |
| 重启 | 8 组选项全部通过（3 组逐位一致，其余 ≤ 3e-12） |
| difftest | 2D/3D 与第二轮逐项相同，0 未解释 |

### 7.2 单元层：原 ADINA 程序 vs OpenSees（同一单元、同一步、完全 Newton，力容差 0.03 N）

| 算例 | 最大差 | 说明 |
|---|---|---|
| 2D 单轴压缩至峰值 | < 30 N（峰值 2.99e7 N） | 压碎后侧向自由单元病态，两程序侧向位移均 17.4 m |
| 2D 单轴拉伸开裂 | 9e-10 N | |
| 2D 压后剪 | 1.5%（对原程序）；21 N（对修正已知缺陷的 ADINA） | 原程序 PRNCPL/CRAKID 缺陷 |
| 桁杆往复（ADINA 模型 5 随动强化 vs Steel01） | 4e-9 N | |
| 3D 单轴压缩 / 拉伸 / 不等侧压三轴 | 0.065 N / 5e-10 N / 0.05 N | |
| 3D 等侧压三轴 | 0.74% | 重复主应力，原 PRNCP3 近重根近似 |

### 7.3 构件层：原 ADINA 程序 vs OpenSees vs 试验

| 模型 | 原 ADINA | OpenSees 默认 | 试验 |
|---|---|---|---|
| 2H20-Ref 梁（剪力 V） | 完成，峰值 105.1 kN | 完成，105.0 kN（对修正缺陷 ADINA 全程差 ≤ 0.022 kN） | 66.0 kN |
| B07 梁（总荷载 P） | 1.69 mm、171 kN 处不收敛停止 | 同类停止 | 333.6 kN @ 7.66 mm |
| W28 墙（平面应力等效，推覆） | 3 mm 892 kN、6 mm 974 kN | 3 mm 912 kN、峰值 975 kN | 峰值 290 kN |
| W28 往复 | 第一次反向后不收敛 | 走完两圈 | — |
| sw6-3 墙重力阶段（单节点 1155 kN） | 0.85–0.884 倍荷载达倒塌荷载 | 0.85–0.891 | — |

### 7.4 离线试验库（`OpenSees-Cases/offline-downloads/concrete-validation-cases-offline-*`，默认 = ADINA 设置）

只替换混凝土本构（HEOM 属性映射 E=p8、ν=p14、ft=p2、fc=−p1、εc=p4、fu=p3、εu=p5；SHEFAC 取 ADINA 默认 0.5），
严格驱动：力容差 1e-4×试验峰值，不放宽、不跳步。覆盖率 = 计算到的挠度 / 试验最大挠度。

| 算例 | 状态 | 覆盖率 | NMAE % | 峰值误差 % |
|---|---|---|---|---|
| 2H16-Ref | 完成 | 1.00 | 9.3 | +29.9 |
| 2H20-Ref | 完成 | 1.00 | 22.1 | +56.2 |
| B07 | 试验范围后停止 | 1.00 | 14.4 | +4.7 |
| B08 / B13 / B14 | 中途停止 | 0.86 / 0.60 / 0.92 | 9.3 / 6.7 / 14.8 | — |
| VS-C3 / VS-OA1 | 中途停止 | 0.38 / 0.48 | 6.8 / — | — |
| deep W1 | 时限 | 0.46 | 7.4 | — |
| 2H16-Void / 2H20-Void | 中途停止（加载点下压碎） | 0.08 / 0.15 | 45.6 / 19.9 | — |
| 2H16-Truss / 2H20-Truss | 停止于模型自带“粘结链接”2 mm 激活之后 | 0.11 / 0.18 | 1.1 / 0.9 | — |
| W28–W32 层壳墙 | 早期停止 | — | — | —（W28 见 7.3） |
| 6 块 GFRP 板 | 本轮未用新默认重跑（每步耗时长，原 2 层混凝土无法收敛） | — | — | — |

### 7.5 三维实体梁（Bresler–Scordelis；数据包冻结主模型，只换本构；同一评分工具）

| 试件 | 状态 | 计算到挠度 mm | 重合段 NRMSE % | 末点误差 % | 数据包 ASDConcrete3D 主试算 NRMSE % |
|---|---|---|---|---|---|
| BS-OA1 | 完成执行器目标 | 8.83 | 4.6 | +4.1 | 20.6 |
| BS-C1 | 最小步长停止 | 12.90 | 4.2 | — | 11.0 |
| BS-A1 | 时限 5400 s | 11.96 | 4.3 | — | 7.2 |
| BS-B1 | 时限 5400 s | 11.24 | 7.8 | — | 10.9 |

所有保存步残差 ≤ 10 N。A1/B1/C1 未到试验末点，末点误差为空（不是 0）。

### 7.6 13 面剪力墙（`run_all_walls.py --driver dc`，默认 = ADINA 设置）

7 面（sw1-4、sw2-3、sw4-1、sw4-2、sw5-1、sw6-1、sw6-3）在重力阶段 78–98% 轴压处停止（单节点集中轴压，见 7.3 sw6-3）；
其余在横向阶段停止：sw1-1 ±10 mm（达到点平均误差 14.5%）、sw1-2 7.65 mm（82.4%）、sw2-1 4 mm（15.0%）、
sw2-2 ±4 mm（36.4%）、sw1-3 0.95 mm、sw5-3 1.66 mm。

### 7.7 优选案例（2026-10-07，19 例，第四轮复核）

用 `github.com/luxinzheng/concrete-validation-cases` 的 `preferred/2026-10-07` 新目标曲线重新运行与评价，
并对效果不好的算例用原 ADINA 程序跑同一模型。工具、全部数字与结论：
[`tests/concretePreferred/RESULTS.md`](https://github.com/luxinzheng/OpenSeesCUDA/blob/33c20525d5c73cbc5661cd782e434068275ed12e/tests/concretePreferred/RESULTS.md)，程序说明 R4 节。要点：

- 原 ADINA 对照：VS-OA1 原 ADINA 1.0 mm 停止（OpenSees 4.7 mm）；VS-C3 4.4–4.6 mm（OpenSees 30.6 mm）；2H16-Void 1.70 mm
  （OpenSees 2.73 mm，二者初始刚度相同，均约为试验 2.7 倍）；7 面墙（平面应力等效 + 加载梁，往复）原 ADINA 均在早期循环停止。
  **效果不好的算例原 ADINA 同样停止且更早，属 ADINA 模式 5 的局限；本轮未发现迁移错误，未修改材料程序。**
- 表现良好：BS-OA1 三维 NMAE 7.3%（末点 −9.2%）；BS-A1/C1 已算段 NMAE 5.5%/3.7%；2H16-Ref NMAE 9.2%（峰值 +29.9%）。
- 墙的加载：存档模型单节点集中轴压使 sw1-3/1-4/2-3/4-1 在重力阶段失败；`run_preferred.py --wall-loading beam`
  用 `rigidLink beam` 表示试件加载梁后全部通过重力（往复中仍在 2–7 mm 停止）。两种加载的结果都保留。

---

## 8. 复现

构建（MSYS2 MinGW；必须给编译器可写临时目录）：
```
export PATH=/mingw64/bin:/usr/bin:$PATH TMP=<可写目录> TEMP=<可写目录>
cmake -S . -B bm -DCMAKE_BUILD_TYPE=Release && cmake --build bm --target OpenSees -j8
```
测试（输出目录须为新目录；原 ADINA 源码 `all_adina.for` 为只读输入，**不在本仓库中**，由使用者自备）：
```
python tests/materials/adina_2d/run_tests.py --adina <all_adina.for> --opensees-build bm --opensees bm/OpenSees.exe --build tmp/x-2d
python tests/materials/adina_3d/run_tests.py --adina <all_adina.for> --opensees-build bm --opensees bm/OpenSees.exe --build tmp/x-3d
python tests/adinaConcreteAudit/run_audit.py --output tmp/x-audit
python tests/adinaConcreteAudit/run_fracture_tests.py --output tmp/x-fr --opensees bm/OpenSees.exe --default-evidence tmp/x-audit
python tests/materials/adina_review/restart/verify_restart.py --opensees bm/OpenSees.exe --out tmp/x-restart
python tests/materials/adina_review/difftest/difftest.py --kernel SRC/material/nD --out tmp/x-d2
python tests/materials/adina_review/difftest/difftest3d.py --out tmp/x-d3
# 原 ADINA 程序对比（先按 adina84/README.md 从 all_adina.for 编译 adina84.exe）
python tests/materials/adina_review/adina84/validation_a.py <tag>
python tests/materials/adina_review/adina84/validation_3d.py <tag>
# 试验算例（需要 OpenSees-Cases 数据包）
python tests/materials/adina_review/cases2/run_offline_cases.py --out tmp/x-off --variants original --jobs 3 --timeout 1800
python tests/materials/adina_review/cases3d/run_solid_cases.py --out tmp/x-solid --variants original --seconds 5400 --jobs 3
python tests/materials/adina_review/run_all_walls.py --bin <含 manifest.json 的程序目录> --out tmp/x-walls --materials original --driver dc
```
各脚本中的数据包路径、Python 路径为本机路径，迁移到其他机器时需要修改。

---

## 9. 文件索引

| 路径 | 内容 |
|---|---|
| `SRC/material/nD/adinaConcrete2D/Concrete2D.{h,cpp}` | 2D 内核（CDMOD/PRNCPL/CRAKID/DCRACK） |
| `SRC/material/nD/adinaConcrete3D/Concrete3D.{h,cpp}`、`BranchTangent.h`、`TensionSoftening.h` | 3D 内核与共用参数、有限差分切线、可选断裂能 |
| `SRC/material/nD/AdinaConcrete2D.*`、`AdinaConcrete3D.*`、`*Factory.cpp`、`AdinaConcreteOptions.h` | OpenSees 适配器与命令解析 |
| `tests/materials/adina_2d/`、`tests/materials/adina_3d/` | 回归测试 |
| `tests/adinaConcreteAudit/` | 逐子程序审计与断裂能测试；`LLM_HANDOFF.txt` 为最初的问题交接说明 |
| `tests/materials/adina_review/PROGRAM_DESCRIPTION.md` | **详细技术说明（R3 节最新）** |
| `tests/materials/adina_review/adina84/` | 原 ADINA 程序卡片生成、运行、解析与对比脚本；`README.md` 记录映射方法和该版本程序的输出缺陷；`ADINA84_INPUT_CARDS.md` 为手册输入卡片格式整理 |
| `tests/materials/adina_review/difftest/` | 整路径逐步对比 |
| `tests/materials/adina_review/cases2/`、`cases3d/`、`run_all_walls.py`、`wall_driver_dc.tcl` | 试验算例运行器与严格驱动 |
| `tests/materials/adina_review/production-*`、`src-before/`、`history/`、`CHANGES_*.diff` | 各阶段源码快照与差异 |
| `EXAMPLES/AdinaConcrete/` | 使用例子 |
| `AdinaConcrete/results/` | 第三轮结果：试验算例的生成模型、曲线、状态与汇总；ADINA 对比表；回归日志 |

## 10. 未包含的内容

- 原 ADINA 源码 `all_adina.for`、手册 PDF、编译出的 ADINA 程序（第三方专有内容）；
- 试验数据包本体（`OpenSees-Cases/...`，体积大，由数据提供方单独分发）及文献 PDF；
- 编译好的 OpenSees 可执行文件（本机生产版本 SHA256 `3ef7ff6f0b19d7cd2f4c07952d1163caa7685d9cb30172757d701d6f839a1512`）。
