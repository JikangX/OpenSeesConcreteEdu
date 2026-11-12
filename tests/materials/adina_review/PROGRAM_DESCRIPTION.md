# OpenSees 中 ADINA 混凝土本构（模式 5）程序说明（第四轮修订版）

供独立项目检查。本文描述 `<合并前工作区>/SRC/material/nD/` 中**当前生产源码**：
组成、相对原迁移版本的全部修改及依据、命令接口、状态与存盘格式、验证范围与结果、
已知问题，以及可直接运行的复现命令。

- 日期：2026-10-07（第三轮）；第二轮 2026-10-06（已归档 `history/PROGRAM_DESCRIPTION_round2.md`）
- 参与者：前一 LLM（初次迁移）、Claude 第一轮复查、Codex 后续开发（见 `DEVELOPMENT_REPORT.txt`）、
  Claude 第二轮（3D 迁移错误修复、3D 全路径分类、13 面墙）、Claude 第三轮（**原 ADINA 程序整体对比**、
  默认算法更正、离线试验库与三维实体梁、ShellNLDKGQ 失败传递）。
- **第四轮（2026-10-07，优选案例复核）见 R4 节；第三轮结论仍有效**：R3 节与第 0–11 节冲突时以 R3 为准（第 0–11 节中已在相应处标注更正）。

---

## R4. 第四轮：优选案例（2026-10-07）复核（以本节为准）

用户在 `github.com/luxinzheng/concrete-validation-cases` 发布了 19 例“最新优选案例（2026-10-07）”（原图重描的目标曲线；
VS-C3 纠正了历史档案误取分析线的问题，VS-OA1 补齐上升段，7 面 SW 墙改为正、负骨架）。本轮用这批目标曲线
重新运行和评价，并对效果不好的算例用原 ADINA 程序跑同一模型。工具与全部数字：`tests/concretePreferred/`
（`README.md`、`RESULTS.md`）。

**ADINA 材料源码本轮没有修改**（默认仍为第三轮的原 ADINA 设置）。在所有可以与原程序对照的算例中，没有发现迁移错误：
效果不好的算例原 ADINA 同样（或更早）停止。

### R4.1 结果概要（默认 = 原 ADINA 设置；覆盖率 = 计算到的位移 / 试验最大位移）

| 案例 | OpenSees（AdinaConcrete） | 原 ADINA（同一模型） | 判断 |
|---|---|---|---|
| 2H16-Ref | 完成；NMAE 9.2%，峰值 +29.9% | — | 剪切梁强度高估（第三轮 2H20-Ref 原 ADINA +59%） |
| 2H20-Ref | 完成；NMAE 22.2%，峰值 +56.2% | 第三轮：完成，+59% | ADINA 模型局限 |
| VS-OA1 | 4.7 mm 停止（0.48） | 1.0 mm 停止（步长 /8，ITEMAX 200） | ADINA 模型局限（脆性开裂不收敛） |
| VS-C3 | 30.6 mm 停止（0.41） | 4.4 mm（原步长）/ 4.6 mm（/4）停止 | 同上 |
| 2H16-Void | 2.7 mm 停止（0.08），初始刚度为试验 2.7 倍 | 1.70 mm 停止（/8），刚度与 OpenSees 相同 | ADINA 与有限元模型局限（存档 HEOM 模拟也约 2 倍刚） |
| 2H20-Void | 3.6 mm 停止（0.15） | 未运行（模型含 ElasticPPGap，转换器不支持） | 同类 |
| 2H16/2H20-Truss | 1.65 / 1.53 mm 停止（0.11/0.18），在存档模型的“tie-bond 链接”激活之后 | 未运行（零长单元 Parallel/ElasticPPGap） | 第三轮：去掉该钩子仍在 3.7 mm 停止 |
| 0F-60S（GFRP 板，分层壳） | 0.11 mm、90 kN 首次开裂即停止 | 无（ADINA84 转换器无壳单元） | 与 B07/VS 同一机理（脆性开裂），未能直接对照 |
| BS-OA1（三维） | 完成；NMAE 7.3%，末点 −9.2% | 无（三维梁含嵌入筋/接触） | 良好 |
| BS-A1 / BS-C1（三维） | 11.96 mm（0.86）/ 12.90 mm（0.72），均以最小位移增量停止（A1 时限 4 h） | 无 | 已算段 NMAE 5.5% / 3.7% |
| 7 面 SW 墙（存档加载） | sw1-3/1-4/2-3/4-1 重力阶段停止；sw1-2/2-1/2-2 横向 4–7.7 mm 停止 | 第三轮 sw6-3：单节点轴压原 ADINA 同样停止 | 有限元模型的加载方式问题 |
| 7 面 SW 墙（刚性加载梁） | 全部通过重力；横向 2–7.3 mm 停止（覆盖 0.24–0.56） | 平面应力等效 + 加载梁：原 ADINA 在 +3.8～+5.0 / −0.2～−2.8 mm 停止（ITEMAX 200） | 往复中的不收敛是 ADINA 模型局限；OpenSees 不早于原 ADINA 停止 |

### R4.2 本轮的程序/工具修改

| 编号 | 内容 |
|---|---|
| D1 | 不修改 ADINA 材料。 |
| D2 | 墙模型的加载梁（`run_preferred.py --wall-loading beam`）：`rigidLink beam` 把顶边节点与控制节点刚接，总轴力与水平协议不变；驱动的提交后残差检查对刚接组按合力（绕控制节点）计算，否则会把刚接约束力误计为不平衡力（未修正时 887 kN 假“残差”）。这是有限元模型对试件加载梁的表示，不是材料修改；存档加载方式的结果同时保留。 |
| D3 | `adina84/beam_run.py`：新增 VS-OA1、VS-C3、2H16-Void、2H20-Void 的定义（ft 按 0.33√fc 规则，与 OpenSees 运行一致）；VS-C3 使用 W28 已记录的桁杆输出缺陷绕过方法（桁杆组不存盘、IREFOR=0），否则原程序在第 1 步后以 “INPUT ERROR IN STRAIN-RATE EFFECT” 停止。 |
| D4 | 新增 `adina84/sw_cyclic.py`：墙的平面应力等效（同 `sw_run.build_2d`）+ 刚性弹性加载梁（一排 E×100 的弹性单元）+ 重力 + 往复位移，原 ADINA 与配对 OpenSees 二维模型。sw1-4：原 ADINA 在 +4.6 mm 停止；配对 OpenSees 越过 +5 mm，在 −4.8 mm 停止。 |

### R4.3 结论

按“若原 ADINA 效果也不好，即为 ADINA 模型自身局限”的原则：本批优选案例中效果不好的梁（VS、Void）和墙（往复）
在原 ADINA 中同样停止，且原 ADINA 停止得更早；剪切梁强度高估（2H 系列）与原 ADINA 相同。OpenSees 中的 ADINA
迁移在这些算例上没有发现新的错误，因此本轮不修改材料程序；唯一的模型层面改进是墙的加载梁（D2）。

---

## R3. 第三轮修订（2026-10-07，以本节为准）

### R3.1 用原 ADINA 程序运行同一有限元模型（新增验证手段）

按用户要求“模型不应违反 ADINA 理论；某算例效果不好时与 ADINA 的计算结果对比，若 ADINA 也不好即为 ADINA
模型自身局限”，将只读源码 `C:/Work/Adina/all_adina.for` 原样编译为完整程序
（`tmp/adina84-build/adina84.exe`，`gfortran -std=legacy -fno-automatic -ffixed-line-length-72 -O1`，
源码未改），并编写输入卡片生成器与结果解析器（`tests/materials/adina_review/ADINA84/`：`deck.py`、`parse.py`、`run.py`、
`README.md` 记录全部映射与该程序版本的输出缺陷及绕过方法）。所有运行在 `tmp/adina84-runs/`。

### R3.2 重要更正：原 ADINA 每次平衡迭代都从收敛态积分（第 5.2 节 B2 原结论错误）

UNBLD 中每次平衡迭代从单元组文件 IT2 读出单元数据（应力、应变、历史数组 WA），调用 ELEMNT 后**不写回**；
历史只在下一步形成刚度时由收敛增量写入。插桩副本（只在 CDMOD 中加 WRITE）运行 B07 卡片：同一步内全部
迭代（例如第 6 步 7 次迭代）的 CDMOD 输入旧应力与 EVMAX 完全相同，均为上一收敛态。且每个增量 CDMOD 只调用
一次（不分子步）。因此与原程序一致的设置是 `-update committed` 且不细分子步。

**生产默认值已改为**：`-update committed`、`-maxStrainIncrement 0`（0 = 不细分；>0 仍可细分），
`-tangent secant`、`-crushedUnloading original` 不变。`iterative`、`hybrid` 保留为**非 ADINA** 的可选算法。
（`hybrid`：迭代原位更新，但在本步无新开裂/压碎事件的点上采用收敛态的加载函数最大值作为加/卸载参照，
消除 iterative 中“超调—回退”迭代把压应力棘轮式推高到曲线以上的问题，单点测试中最多 +36%。）
加载判别由“f 与最大值相等”改写为“f ≥ 最大值（相对 1e-8）”：在原算法中 f 不会超过已记录最大值，两者等价，
回归与 difftest 结果逐位不变；hybrid 原点下二者才有区别。

### R3.3 同一模型 ADINA 与 OpenSees 的对比结果

| 对比 | 内容 | 结果 |
|---|---|---|
| 二维单元 A1 单轴压缩 | 4 节点平面应力单元，E=30e9，ft=3e6，fc=-30e6，完全 Newton，力容差 0.03 N | 至峰值（2.99e7 N）逐步差 < 30 N；压碎步之后单个侧向自由单元病态（两程序侧向位移均 17.4 m） |
| 二维 A2 单轴拉伸开裂 | 同上 | 差 9e-10 N |
| 二维 A3 压后剪 | 同上 | 与原程序差 1.5%；与修正三处已记录缺陷的 ADINA 差 21 N |
| 桁杆 A4 往复 | ADINA 桁杆模型 5（随动强化）对 Steel01 | 差 4e-9 N |
| 三维单砖 B1 单轴压缩 / B2 拉伸 / B4 侧压不等的三轴 | `validation_3d.py` | 0.065 N / 5e-10 N / 0.05 N（峰值 3e7–6.8e7 N） |
| 三维单砖 B3 等侧压三轴 | 同上 | 0.74%：重复主应力，原 PRNCP3 近重根近似（第 8.2 节已记录） |
| B07（Bresler–Scordelis OA-1，平面应力梁） | 完全 Newton，力容差 1e-4×试验峰值 | **原 ADINA 在 1.69 mm、171 kN 处不收敛停止**（步长 /8、/32，ITEMAX 200/999 均同一步）；试验峰值 333.6 kN 于 7.66 mm |
| 2H20-Ref（Mak–Lees 梁） | 同上，ITEMAX 200，步长 /4 | **原 ADINA 完成，V 峰值 105.1 kN，试验 66.0 kN（高估 59%）**；OpenSees 默认 105.0 kN，与修正缺陷后的 ADINA 全程差 ≤ 0.022 kN |
| W28 墙（层壳模型转为同网格平面应力 + 桁杆，单调推覆至 6 mm） | 0.01 mm 步 | **原 ADINA 3 mm 时 892 kN，6 mm 时 974 kN；试验峰值 290 kN**。OpenSees 默认 912 kN（3 mm 前与原 ADINA 差 ≤ 2.3%；峰后差异来自原程序三处已记录缺陷，修正缺陷后的 ADINA 与 OpenSees 全程差 ≤ 1.6 kN）。往复：原 ADINA 在第一次反向、裂缝闭合后即不收敛 |

结论：OpenSees 默认设置与原 ADINA 在单元层和构件层一致（差异只来自原程序已记录的缺陷）。
**剪切控制构件的强度高估与往复/卸载中的不收敛是 ADINA 模式 5 本身的局限**（无受拉软化、固定裂缝常数剪力保留
SHEFAC=0.5 导致应力锁定、脆性开裂使 Newton 循环），不是迁移错误。按用户要求，不再以非 ADINA 方法追求这些算例的精度：
第三轮试过的可变剪力保留、受压软化、统一 ft 规则、裂缝带在梁中失稳等均只作诊断记录，未进入生产。

### R3.4 其他程序修改

| 编号 | 内容 |
|---|---|
| C1 | `ShellNLDKGQ`：新增 `update()`，截面（材料）拒绝试算时把失败返回给 domain::update，分析回退并可二分步长；此前失败步仍以旧应力参与组装并打印上百万条警告（ShellMITC4 已有同样处理） |
| C2 | 适配器的拒绝信息限频：前 20 条及每 100000 条打印一次 |
| C3 | `-update hybrid`、`-maxStrainIncrement 0`、默认 committed（R3.2）；存盘格式枚举 0/1/2 兼容旧记录 |
| C4 | 测试：内核回归新增 hybrid 与默认值检查；原依赖 2e-5 子步的检查显式设置 `-maxStrainIncrement 2e-5`；失败传递测试同样显式设置；`-maxStrainIncrement 0` 改为合法值、负值拒绝 |

### R3.5 回归（第三轮生产版本，bm/OpenSees.exe）

| 测试 | 结果 |
|---|---|
| 二维 `adinaConcrete2D/run_tests.py`（`reg2d-adef5`） | 234,848 内核、6,755 适配器/截面、壳失败传递、3,961,660 RC 膜片检查，全部通过 |
| 三维 `adinaConcrete3D/run_tests.py`（`reg3d-adef2`） | 37,358 内核、1,980 适配器、954 Tcl，全部通过 |
| 审计 `run_audit.py`（`audit-adef`） | 7,877 条，26 组 0 失败 |
| 断裂能（`fracture-adef`） | 62,736 内核、69,094 分层壳 |
| 重启（`restart-adef`） | 8 组选项（含 hybrid）全部通过 |
| difftest 2D（`d2-r3`）/ 3D（`d3-r3`） | 与第二轮逐项相同；3D 92 组 0 未解释 |

### R3.6 离线试验库（`OpenSees-Cases/offline-downloads/concrete-validation-cases-offline-20261006-7e7ae3874d8b`）与三维实体梁

运行器：`cases2/run_offline_cases.py`（只替换混凝土；HEOM 属性映射 E=p8, ν=p14, ft=p2, fc=−p1, εc=p4, fu=p3, εu=p5，
SHEFAC 取 ADINA 默认 0.5；严格驱动 `cases2/offline_driver.tcl`：NormUnbalance = 1e-4×试验峰值，不放宽、不跳步，
失败二分；原 sp+LoadControl 位移协议改为等价的 DisplacementControl——同一位移历程，避免非平衡预测步）。
三维：`cases3d/run_solid_cases.py`（`OpenSees-Cases/selected-solid-3d-benchmark-20261006` 冻结主模型，仅把 ASDConcrete3D 换为 adinaConcrete3D，
保留其驱动、10 N 残差门槛与评分工具）。结果目录 `tmp/r5-offline`、`tmp/r5-solid`。

二维（默认 = ADINA 设置；覆盖率 = 计算达到的挠度/试验最大挠度；峰值误差只在覆盖到试验峰值挠度时给出）：

| 算例 | 状态 | 覆盖率 | NMAE % | 峰值误差 % | ADINA 对照 |
|---|---|---|---|---|---|
| 2H16-Ref | 完成 | 1.00 | 9.3 | +29.9 | — |
| 2H20-Ref | 完成 | 1.00 | 22.1 | +56.2 | 原 ADINA +59%（R3.3） |
| B07 | 试验范围后停止 | 1.00 | 14.4 | +4.7 | 原 ADINA 1.69 mm 即停止 |
| B08 / B13 / B14 | 停止 | 0.86 / 0.60 / 0.92 | 9.3 / 6.7 / 14.8 | — | 未运行 |
| VS-C3 / VS-OA1 | 停止 | 0.38 / 0.48 | 6.8 / — | — | 未运行 |
| deep W1 | 时限 | 0.46 | 7.4 | — | 未运行 |
| 2H16-Void / 2H20-Void | 停止 | 0.08 / 0.15 | 45.6 / 19.9 | — | 未运行 |
| 2H16-Truss / 2H20-Truss | 停止 | 0.11 / 0.18 | 1.1 / 0.9 | — | 停止发生在模型自带的“tie-bond 链接”在 2 mm 激活之后；去掉该钩子可继续（诊断 `tmp/truss-nohook`）。属历史模型附加功能，不是混凝土问题 |
| W28–W32（层壳墙） | 早期停止 | — | — | — | W28 原 ADINA 同样高估约 3 倍且反向即不收敛（R3.3） |

三维实体梁（Bresler–Scordelis，总荷载 P，与数据包同一评分）：

| 试件 | 状态 | 计算到挠度 mm | 重合段 NRMSE % | 末点误差 % | 数据包 ASDConcrete3D 主试算 NRMSE % |
|---|---|---|---|---|---|
| BS-OA1 | 完成执行器目标 | 8.83 | 4.6 | +4.1 | 20.6 |
| BS-C1 | 最小步长停止 | 12.90 | 4.2 | — | 11.0 |
| BS-A1 | 时限 5400 s | 11.96 | 4.3 | — | 7.2 |
| BS-B1 | 时限 5400 s | 11.24 | 7.8 | — | 10.9 |

所有保存步残差 ≤ 10 N、数值有限。OA1 计算最大力 394 kN 出现在试验末点之后，没有试验峰后数据可比；
A1/B1/C1 未到试验末点，末点误差为空（不是 0）。

### R3.7 13 面墙（`run_all_walls.py --driver dc`，`tmp/r5-walls`，默认 = ADINA 设置）

13 面墙中 7 面（sw1-4、sw2-3、sw4-1、sw4-2、sw5-1、sw6-1、sw6-3，两种卸载选项相同）在**重力阶段**（约 78–98% 轴压）
即停止：轴压以单个顶部节点集中力施加（如 sw6-3 节点 53 的 1155 kN），节点下单元压碎后 ADINA 压碎律无刚度。
其余在横向阶段停止：sw1-1 到 ±10 mm（误差 14.5%），sw1-2 到 7.65 mm（82.4%），sw2-1 到 4 mm（15.0%），
sw2-2 到 ±4 mm（36.4%），sw1-3 0.95 mm，sw5-3 1.66 mm。第二轮的 iterative 默认下走得更远，但那不是 ADINA 的算法。
（重力阶段的原 ADINA 对照见 R3.8。）

### R3.8 重力阶段的原 ADINA 对照（sw6-3，`ADINA84/sw_run.py`，`adina84-runs/sw6-3-*`）

把 sw6-3 层壳墙按同网格转为平面应力 + 桁杆，只施加顶部节点 53 的 1155 kN 集中力：
原 ADINA 在 0.85–0.884 倍荷载处以 “STIFFNESS MATRIX NOT POSITIVE DEFINITE”（达到倒塌荷载）停止，与步数无关；
OpenSees 默认 0.85–0.891（与修正缺陷后的 ADINA 同一步停止），`-update iterative` 0.90–0.908。
位移控制的承载力检查：原 ADINA 峰值 1022–1024 kN，OpenSees 默认 1033–1035 kN，层壳原模型 1122 kN（0.971，
与其在 t≈0.979 停止一致）。荷载点正下方两个积分点压碎（约 −30.2 MPa）。
**结论：单节点集中轴压超过该模型在 ADINA 压碎律下的局部承载力，原 ADINA 同样无法施加全部轴压**，
属 ADINA 模式 5 与模型加载方式的组合局限，不是迁移或求解器缺陷。其余 6 面重力阶段停止的墙未逐一用 ADINA 运行。

### R3.10 复现命令（在 `<合并前工作区>` 下，输出目录须为新目录）

```
python tests/materials/adina_review/ADINA84/validation_a.py <tag>          # 二维单元：原 ADINA 对 OpenSees
python tests/materials/adina_review/ADINA84/validation_3d.py <tag>         # 三维单砖
python tests/materials/adina_review/ADINA84/beam_run.py ...                # B07 / 2H20-Ref（见该目录 README）
python tests/materials/adina_review/ADINA84/wall_run.py ...                # W28 平面应力等效墙（见 README）
python tests/materials/adina_review/cases2/run_offline_cases.py --out tmp/<new> --variants original --jobs 3 --timeout 1800
python tests/materials/adina_review/cases3d/run_solid_cases.py --out tmp/<new> --variants original --seconds 5400 --jobs 3
python tests/materials/adina_review/run_all_walls.py --bin <pair dir> --out tmp/<new> --materials original elastic --driver dc
```
`cases2/table.py <dir> ...` 汇总离线结果。源码快照：`production-r3/`；相对第二轮最终版的差异：`CHANGES_ROUND3.diff`。

### R3.9 使用建议

1. 与 ADINA 一致：使用默认选项。需要与原程序逐步对照时，OpenSees 侧用 `constraints Lagrange`（Transformation 会在首次求解前写入 sp 值）。
2. 结构算例推进不了时，先用 `tests/materials/adina_review/ADINA84` 在原程序中跑同一模型再判断。
3. iterative/hybrid 只作为非 ADINA 的数值选项，结果与原程序不同，不得称为 ADINA 结果。
4. 集中轴压请用加载梁分布；W28 存档模型的截面只用了竖向钢筋层（0° 水平分布筋层定义了但未引用），且无轴压，需对照试件图纸核实。


---

## 0. 交付物与身份

| 项目 | 位置 / 哈希 |
|---|---|
| 生产源码 | `SRC/material/nD/` 中 13 个文件（第 2 节）；快照 `production-final/` |
| 截面修改 | `SRC/material/section/LayeredShellFiberSection.cpp/.h`；快照 `history/section-final/`，相对 git HEAD 的差异 `history/section-final/LayeredShellFiberSection_vs_git_HEAD.diff`（含此前 LLM 的修改） |
| 相对最初迁移版本的全部差异 | `CHANGES_FINAL.diff`（`src-before` → `production-final`，10 个文件，+493/−69） |
| 第二轮差异 | `CHANGES_SECOND_ROUND.diff`（Codex 报告时的生产版本 → 最终，仅 `Concrete3D.cpp` +12/−1） |
| 生产可执行文件 | `bm/OpenSees.exe` SHA256 `2e4ae86f1799dad52c3fcaac0f754140d6482c49d2909596107351b96cac3a90` |
| 配对可执行文件 | `tmp/claude-review/pair-final/`：`OpenSees-ADINA2D.exe`（同上）、`OpenSees-PlainConcrete.exe` `250816602bea21a056ba82fb0ffb697f3346f57817a1e4bc0adeda2144889474`（原 forumat.for + PlainConcrete.f90，与 ADINA 同一引擎对象链接） |
| 原 ADINA 源码（只读） | `C:/Work/Adina/all_adina.for` SHA256 `9a0e2731af05ad4eb76b8fabd890e9492896a296aa3770a560fd7fde180ffb1b` |

## 1. 范围

ADINA 模式 5 混凝土：
- 二维平面应力 `AdinaConcrete2D` / `AdinaConcrete2DPlateFiber`（原 A06.FOR：CDMOD、PRNCPL、CRAKID、DCRACK），
  用于 `LayeredShell` + `ShellMITC4` / `ShellNLDKGQ`；
- 三维 `AdinaConcrete3D`（原 A16.FOR：CMOD3D、PRNCP3、CRAK3D、DCRAK3）。

不在范围：温度、有限变形、模式 4、平面应变/轴对称、率相关。可选断裂能（前一 LLM 的扩展）保留并回归通过，
其多轴/混合模式的能量客观性未证明。

## 2. 源文件与职责

| 文件 | 职责 |
|---|---|
| `adinaConcrete3D/Concrete3D.h/.cpp` | 两维共用的参数 `Parameters`；三维状态、谱分解、包络、有理曲线、裂缝状态机、压碎、迭代矩阵 |
| `adinaConcrete3D/BranchTangent.h` | 有限差分算法切线（`-tangent algorithmic`） |
| `adinaConcrete3D/TensionSoftening.h` | 可选断裂能裂缝带模块 |
| `adinaConcrete2D/Concrete2D.h/.cpp` | 二维平面应力积分、DCRACK 刚度、迭代矩阵 |
| `AdinaConcrete2D.h/.cpp`、`AdinaConcrete3D.h/.cpp` | NDMaterial 适配：试算/提交/回滚/复制/存盘/响应 |
| `AdinaConcrete2DFactory.cpp`、`AdinaConcrete3DFactory.cpp` | Tcl/Python 参数解析 |
| `AdinaConcreteOptions.h` | 关键字选项解析；断裂能模式的求解选项默认值 |
| `section/LayeredShellFiberSection.cpp/.h` | 分层壳截面（Codex 修复 `getInitialTangent`，见 5.2-B7） |

内核不依赖 OpenSees 类型、无全局可写状态。

## 3. 约定

受拉为正；`fc, epsc, fu, epsu` 输入为负。二维 `[εxx, εyy, γxy]`，PlateFiber `[xx, yy, xy, yz, xz]`
（横向剪切独立弹性 `τ = Gγ`）。三维 Voigt `[xx, yy, zz, xy, yz, xz]`，工程剪应变。单位自洽。

## 4. 与原 ADINA 例程的对应

| 原例程（all_adina.for 行号） | 迁移位置 |
|---|---|
| CDMOD（25035） | `Concrete2D.cpp` `step()`、`stiffness()`、`secantStiffness()`、`integrate()` |
| PRNCPL（25581） | `principalAngle()`、`strengthEnvelope()` |
| CRAKID（25709） | `classify()`、`step()` 首次开裂段 |
| DCRACK（25930） | `stiffness()`、`step()` 释放段 |
| CMOD3D（35351）/PRNCP3（35866）/CRAK3D（36171）/DCRAK3（36410） | `Concrete3D.cpp` `step()`、`principalStress()`、`updateCracks()`、`localStiffness()` |
| 单元循环 + ECHECK（10017，IUPDT=0） | 适配器 `setTrialStrain()` 的更新起点（5.2-B2） |
| 平衡迭代（约 4440–4600，IEQREF=0 用材料矩阵） | 迭代矩阵（5.2-B1） |

## 5. 修改明细（相对最初迁移版本）

“轮次”：R1 = Claude 第一轮，C = Codex，R2 = Claude 第二轮。

### 5.1 迁移错误修复（与原程序逐步对比发现）

| 编号 | 轮次 | 内容与依据 | 位置 |
|---|---|---|---|
| A1 | R1 | 压碎后，最受压主方向应力不再受压（`≥0.001·SIGCP`）时模量取 0（应力冻结），不再以负模量把应力推向受拉。依据 CDMOD 标号 9 循环、CMOD3D 第 35480–35500 行 | 2D `step()` 压碎分支；3D 同 |
| A2 | R1 | 二维压碎后全释放：用**更新后应力**的主轴（CDMOD 标号 60 再调用 PRNCPL）并检查厚度应变 TEP(4)（DCRACK 第 26111–26116 行） | 2D `step()` |
| A3 | C | 三维压碎后全释放同样用更新后应力主轴（CMOD3D 标号 60 + DCRAK3 IK=4） | 3D `step()` 压碎分支 |
| A4 | R2 | **三维单裂缝时增量前的强度包络**：CRAK3D（KKK=1）标号 111 在恰有一条张开裂缝时，用裂缝法向应力和**横向 2×2 块的两个主值**求 SIGCP/EPSCP；C++ 原用横向的法向分量。有横向剪应力时曲线系数、YP 和 MOD45 阈值不同，造成随机路径中逐步漂移。修复后 `random7` 由偏差变为逐步一致，`random0` 首次分歧由 178 步推迟到 372 步（剩余为原程序缺陷，见 8.2） | 3D `step()` `envelopeStress` |

### 5.2 数值方法与 OpenSees 接口

| 编号 | 轮次 | 内容 |
|---|---|---|
| B1 | R1 | `-tangent secant\|algorithmic`，默认 secant：返回 CDMOD/CMOD3D 在 IEQREF=0 时重建的材料矩阵（加载曲线切线、卸载 E、开裂 STIFAC/SHEFAC、压碎 E·STIFAC；断裂能模式开裂面用黏聚割线）；对称半正定，只用于迭代。与原 CDMOD 返回矩阵抽查吻合到 4–5 位有效数字 |
| B2 | R1，**R3 更正** | `-update committed\|iterative\|hybrid`。R1 称“ADINA 原位更新 WA，故默认 iterative”，**此结论错误**：原程序每次平衡迭代都从上一收敛态积分（见 R3 节），故 R3 起默认 committed；iterative、hybrid 为非 ADINA 的可选算法 |
| B3 | R1/C | 迭代模式下，对与当前试算完全相同的应变请求不再施加零增量（3D R1，2D C），使提交、`Domain::update`、恢复、回滚不改变历史 |
| B4 | R1 | 存盘格式 v3（2D 118 项、3D 198 项），新增选项、`crushedStrain`、二维 `released`；v1/v2 旧记录按其产生时的公式恢复（algorithmic + committed + original），坏枚举拒绝 |
| B5 | C | 3D `getCopy()` 将副本 dbTag 置 0，避免已存盘原型的积分点在数据库中互相覆盖 |
| B6 | C | 2D/3D `recvSelf`：压碎后拉伸历史冻结，不再用旧拉伸最大应变约束后续压碎应变（此前会拒绝合法检查点） |
| B7 | C | `LayeredShellFiberSection::getInitialTangent` 原返回当前损伤切线；改为由各层 `getInitialTangent` 组装。原问题使恢复后 ShellMITC4 钻转惩罚 Ktt 改变，**数据库重启无法复现**（对任何非线性材料，含 OpenSees 自带的 PlasticDamageConcretePlaneStress）。修复后 6 组选项逐步一致 |
| B8 | C | 三维开裂比较统一用 64ε 应力尺度舍入邻域（`exceedsTensileCutoff`），保留 CMOD3D 首次严格比较、CRAK3D 指定分支包含等号的语义 |
| B9 | C | 断裂能模式（`-fractureEnergy`）在未显式指定时默认 `algorithmic + committed`，与选项顺序无关；显式选项优先 |
| B10 | C | 三维断裂能扩展：交互抗拉阈值非正时（已超出压缩破坏面）不建立负强度的拉伸面，交由压碎检查处理 |

### 5.3 可选扩展（非 ADINA 方程）

`-crushedUnloading elastic`（默认 original）：原 CDMOD/CMOD3D 中压碎点卸载时沿软化曲线“爬回”，越过峰值应变后模量为 0，
**应力锁定在峰值压应力，应变变成受拉也不变**（用原 CDMOD 复现：应变 +0.00026 时应力仍为 −32.2 MPa）。
扩展：在历史最大压应变以内按初始 E 卸载/再加载，主应力截断在 `[σ曲线(历史应变), 0]`；超过历史压应变及
极限应变释放仍按原律。往复荷载下推荐使用。

### 5.4 试过并撤销的方案（未进入生产）

| 方案 | 结果 |
|---|---|
| 压碎下降段迭代矩阵取负切线 | sw1-3 首步失败（R1） |
| latched（本步已开裂面保持张开） | 裂缝闭合/再张开两步循环（R1） |
| staggered（试算中延迟脆性事件、提交时施加） | RC 膜片 3967400 项检查含对称性全部通过，但 sw1-1 提交后的释放力在下一步无法平衡，dt→0 时残差仍为 300–550 N（R2） |
| 压碎起始所在步改为从收敛态重放 | 单个 3D 砖多走 3 步，但 sw1-3、sw2-2 回退（R2） |

## 6. 命令接口

```
nDMaterial AdinaConcrete2D           tag E nu ft fc epsc fu epsu [选项]
nDMaterial AdinaConcrete2DPlateFiber tag E nu ft fc epsc fu epsu [选项]
nDMaterial AdinaConcrete3D           tag E nu ft fc epsc fu epsu [选项]
```

| 选项 | 默认 | 说明 |
|---|---|---|
| `-tangent secant\|algorithmic` | secant（断裂能模式：algorithmic） | 迭代矩阵 |
| `-update committed\|iterative\|hybrid` | committed（与原程序相同，R3） | 试算起点；iterative/hybrid 为非 ADINA 可选算法 |
| `-crushedUnloading original\|elastic` | original | 5.3 |
| `-retention STIFAC SHEFAC` | 1e-4 0.5 | |
| `-controls beta gamma kappa alpha` | 0.5 1 0.7 0 | |
| `-envelope 24 个数` | 手册默认 | SP1、SP31、SP32、SP33 |
| `-maxStrainIncrement x` | 0（R3；不细分，与原程序每增量调用一次 CDMOD/CMOD3D 相同） | >0 时按此应变范数细分子步（单次试算 >10000 子步即失败） |
| `-fractureEnergy Gf Lch` | 关闭 | 可选断裂能扩展 |
| `-G G`（仅二维） | E/(2(1+ν)) | 板横向剪切 |
| `-rho rho` | 0 | |

推荐：墙、梁等往复加载 `-crushedUnloading elastic`；求解先用完整 `Newton`；每步位移约 0.2 mm 或更小，并做一次步长减半检查。
与原程序对照时不加选项。

## 7. 状态、存盘与响应

- 2D 状态：应变、应力、加载函数最大值、裂缝角、两个裂缝面的开闭与形成应变、压碎强度/应变、
  `crushedStrain`、`released`、断裂能历史。3D：完整 3×3 裂缝坐标、开裂数（张开面在前）、三个面的历史等。
- 存盘 v3 见 5.2-B4。
- 响应：2D `crackState`（10 项）、3D `crackState`/`crackFrame`（各 9 项）、`fractureState`。

## 8. 验证（最终生产版本上重新运行）

所有输出目录在 `<合并前工作区>/tmp/claude-review/` 下。

### 8.1 二维与原 CDMOD 整路径逐步对比（`difftest/difftest.py`，输出 `d2-final`）

88 组路径（44 种 × 两套参数）、130,544 步，原程序按 IEQREF=0 驱动，另编一份只含三处一行修正的 corrected 参考：
68 组逐步一致（<1e-9·|fc|），6 组舍入级（≤1e-6），6 组开裂阈值平局，8 组原 CRAKID 多裂缝缺陷
（一条闭合另一条丢失；闭合裂缝在垂直方向再开裂时释放错误平面）。无未解释项。

### 8.2 三维与原 CMOD3D 整路径逐步对比（`difftest/difftest3d.py`，输出 `d3k`）

92 组、79,552 步（含第二轮新增的纯剪往复、开裂后剪切、两裂缝闭合再张开、往复压缩及 8 条随机非比例路径）。
同时编译 corrected 参考（P1 第二/三裂缝应变按应力角张量旋转；P2 保留完整坐标系而不用 DCA 的正平方根重建；
P4 闭合比较加 64ε 舍入下限），并对 C++ 每一步检查与两程序都无关的固定裂缝不变量
（张开面法向应力为 0；未张开面法向应力 ≤ ft；无裂缝时最大主应力 ≤ ft）：

| 分类 | 组数 |
|---|---|
| 与原程序逐步一致 | 40 |
| 两者在同一步因主应力超出 SP1 包络停止（ADINA STOP），之前一致到 ≤1.4e-13·|fc| | 16 |
| 阈值平局（应变恰为 0 或恰为存储的开裂应变） | 5 |
| 与 corrected 参考一致（P1/P2/P4） | 7 |
| 原 CRAK3D 三裂缝闭合缺陷（独立参考校核到 1e-17·|fc|） | 4 |
| 原程序在首次分歧处违反固定裂缝不变量（释放/保留了错误的物理面） | 20 |
| 未解释 | **0** |

C++ 在全部 92 组路径的每一步上不变量偏差 ≤1.5e-10·|fc|。

### 8.3 回归测试

| 测试 | 结果 |
|---|---|
| 逐子程序审计 `run_audit.py`（`audit-final`） | 7,877 条，26 组全部 0 失败 |
| 二维回归 `adinaConcrete2D/run_tests.py`（`reg2d-b`） | 225,475 内核、6,755 适配器/截面、7,136 Tcl 壳、失败传播、3,967,400 RC 膜片（committed 公式，见 9.3） |
| 三维回归 `adinaConcrete3D/run_tests.py`（`reg3d`） | 35,193 内核、1,980 适配器、954 Tcl |
| 断裂能扩展 `run_fracture_tests.py`（`fracture-final`） | 62,736 内核、69,252 分层壳 |
| 跨进程数据库重启 `restart/verify_restart.py`（`restart-final`） | 6 组选项各 41 个续算状态，逐步比较位移/抗力/截面/各层；3 组逐位一致，其余 3 组缩放差 ≤2.9e-12 |

### 8.4 13 面剪力墙（`run_all_walls.py`，输出 `walls13b`）

本地全部 13 面墙（OpenSees-Cases/fagelongxia），模型由原案例文件生成（`run_walls.prepare`），三种材料设置使用
相同几何、钢筋、边界、轴压、位移协议和驱动（`wall_driver_adina.tcl`，完整 Newton，力容差不变），每步 0.2 mm，
未拟合任何参数；所有已接受步的提交后残差都在求解容差内。实验指标：同号、位移在实验骨架点 [0.95, 1.02] 区间内的最大反力。

| 墙 | ADINA original | ADINA elastic | PlainConcrete |
|---|---|---|---|
| sw1-1 | 完成，30.0% | 完成，17.8% | t=45.1 停滞，4.4% |
| sw1-2 | 完成，85.0% | 完成，46.1% | t=19.4，5.0% |
| sw1-3 | t=3.4，5.9% | 完成，40.9% | t=12.2，15.1% |
| sw1-4 | 重力阶段失败 | 完成，19.1% | t=9.9，30.9% |
| sw2-1 | t=11.3，17.6% | t=25.0，18.1% | t=23.0，21.3% |
| sw2-2 | 完成，123.9%（负向力符号错误） | 完成，49.0% | t=45.8，43.9% |
| sw2-3 | t=2.1，17.3% | t=16.7，15.5% | t=8.2，56.7% |
| sw4-1 | 重力阶段失败 | 重力阶段失败 | t=6.5，68.3% |
| sw4-2 | 完成，22.8% | 完成，57.8% | t=12.1，46.5% |
| sw5-1 | 重力阶段失败 | t=7.1，64.2% | t=6.5，26.7% |
| sw5-3 | t=8.8，16.4% | t=2.5，— | t=8.2，19.7% |
| sw6-1 | t=19.5，7.9% | 完成，56.2% | t=8.2，18.1% |
| sw6-3 | t=8.0，7.2% | t=36.5，37.6% | t=12.0，50.3% |

完整走完：ADINA elastic 7/13，original 4/13，PlainConcrete 0/13。
两者都到达的 30 个实验点上平均绝对误差：ADINA elastic 32.5%，PlainConcrete 27.3%。
重力阶段失败的三例全部轴压加在顶部单个节点（模型无加载梁）：ADINA 压碎律在压碎后不能再承受应力增量，
节点下方单元在力控制下局部压溃（PlainConcrete 保留残余压应力）。属模型边界与 ADINA 压碎律的组合，不是程序错误。

### 8.5 梁 2H16-Ref

走到实验最大挠度（最初迁移版本停滞或超时）。初始刚度约为实验的 3 倍，ShellMITC4 与 ShellNLDKGQ 相同，
弹性梁理论（跨度 1.5 m、截面 0.16×0.34 m、E=30 GPa）给出与有限元一致的刚度，差异来自试验装置或假定的 E（Codex
指出成因仍未独立确认）；开裂后 NLDKGQ 超估约 30%，MITC4 约 46%，粗网格面内弯曲建议 NLDKGQ。

### 8.6 三维构件

81 个 stdBrick 素混凝土柱（`brick3d/column_test.tcl`）：默认推进到 t=4.1，旧设置（algorithmic + committed）t=1.85，
均未走完。单个 stdBrick 侧向自由单轴加载见 9.4。三维无实验对比。

## 9. 已知问题与使用限制

1. **迭代路径/步长依赖**（iterative 模式，ADINA 本身的算法特性）：sw1-1 elastic，dt 0.02 与 0.01 在 0–2 mm 的反力 RMS 差约 9%，
   1–7 mm 约 15%，峰后最大差约峰值的 55%。PlainConcrete 约 2%。建议小步长与步长减半检查；
   committed 模式路径无关，但构件中因脆性开裂突降不收敛。
2. **大变形低估**：15 mm 以后低估 40–60%（εu 处全释放、墙脚单排单元局部化、约束区用未约束参数）。
3. **均匀状态对称性**：迭代模式下均匀拉伸膜片出现 7e-8 m 横向位移（舍入经迭代选择不同裂缝事件）；
   committed 模式严格对称。回归测试以 committed 公式检验该严格门槛，门槛未放宽。
4. **三维单个砖压碎**：压碎起始时保留未平衡迭代的预测应力并冻结（峰值应变前模量 0），侧向自由时 10.6 kN 残差永久保留，
   之后侧向刚度为零、问题病态。ADINA 本身同样。单个均匀单元请用 `-update committed -tangent algorithmic`。
5. **顶部集中轴压**：见 8.4，建议模型中设置加载梁。
6. 断裂能扩展：仅 Mode-I 单轴能量经解析验证。

## 10. 构建与复现

环境：MSYS2 MinGW；构建时必须设置可写临时目录，否则 `cc1plus` 报 “Cannot create temporary file in C:\windows\”：
```
C:/msys64/usr/bin/bash.exe --noprofile --norc -c "export PATH=/mingw64/bin:/usr/bin:$PATH TMP=<合并前工作区>/tmp/build-tmp TEMP=<合并前工作区>/tmp/build-tmp; cd /c/Work/OpenSees; cmake --build bm --target OpenSees -j8"
```
运行时 PATH 加入 `C:/msys64/mingw64/bin`。Python：`C:/Users/Administrator/.cache/codex-runtimes/codex-primary-runtime/dependencies/python/python.exe`。
以下命令在 `<合并前工作区>` 下，输出目录须为新目录：
```
python tests/materials/adina_review/difftest/difftest.py --kernel SRC/material/nD --out tmp/x-d2
python tests/materials/adina_review/difftest/difftest3d.py --out tmp/x-d3
python tests/adinaConcreteAudit/run_audit.py --output tmp/x-audit
python tests/materials/adina_2d/run_tests.py --adina C:/Work/Adina/all_adina.for --opensees-build bm --opensees bm/OpenSees.exe --build tmp/x-2d
python tests/materials/adina_3d/run_tests.py --adina C:/Work/Adina/all_adina.for --opensees-build bm --opensees bm/OpenSees.exe --build tmp/x-3d
python tests/adinaConcreteAudit/run_fracture_tests.py --output tmp/x-fr --opensees bm/OpenSees.exe --default-evidence tmp/x-audit
python tests/materials/adina_review/restart/verify_restart.py --opensees bm/OpenSees.exe --out tmp/x-restart
python tests/adinaConcrete2DWalls/build_baseline.py --output tmp/x-pair
python tests/materials/adina_review/run_all_walls.py --bin tmp/x-pair --out tmp/x-walls --jobs 4
```

## 11. 建议的检查要点

1. 对照 `CHANGES_FINAL.diff` 与第 5 节，逐条核对所引原程序语句（尤其 A4：CRAK3D 标号 111 在 KKK=1 时的 SIG(2)/SIG(3)）。
2. 重跑 `difftest3d.py`，抽查 20 组“违反不变量”路径的首次分歧（`report.json` 中
   `raw_invariant_violation_at_first_divergence_over_fc`），判断是否同意原程序缺陷的结论。
3. 检查 corrected 参考的 P1/P2/P4 补丁（`CORRECTIONS`）是否只修复所述缺陷。注意 P2 依赖单点测试中 CMOD3D
   局部变量的静态存储，个别路径（如 `cyclic_compression_rotated`）中该补丁本身会偏离原程序，已在分类中说明。
4. 评估 iterative 更新的路径依赖是否可接受（第 9.1 条），可用 committed 对照。
5. 评估 `-crushedUnloading elastic` 作为非 ADINA 扩展的合理性。
6. 第 8.4 节墙算例的模型边界（单点加载、约束区参数）与实验资料仍需核实，结果不构成实验精度认证。
