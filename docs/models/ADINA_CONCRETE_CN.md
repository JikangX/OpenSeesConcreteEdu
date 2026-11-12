# ADINA 混凝土：当前实现、修订与使用

核查日期：2026-10-08。本文描述整合后的仓库中的材料代码。详细历次修订以 [第四轮程序说明的 R4/R3 节](../../tests/materials/adina_review/PROGRAM_DESCRIPTION.md) 为主要资料；较早段落与之冲突时按 R3/R4 解释。

## 1. 实现组成与理论范围

- [二维适配器](../../SRC/material/nD/AdinaConcrete2D.cpp)、[二维工厂](../../SRC/material/nD/AdinaConcrete2DFactory.cpp) 和 [二维内核](../../SRC/material/nD/adinaConcrete2D/Concrete2D.cpp)。
- [三维适配器](../../SRC/material/nD/AdinaConcrete3D.cpp)、[三维工厂](../../SRC/material/nD/AdinaConcrete3DFactory.cpp) 和 [三维内核](../../SRC/material/nD/adinaConcrete3D/Concrete3D.cpp)。
- [共享参数及默认值](../../SRC/material/nD/adinaConcrete3D/Concrete3D.h)、[工厂选项](../../SRC/material/nD/AdinaConcreteOptions.h)、[分支切线](../../SRC/material/nD/adinaConcrete3D/BranchTangent.h)、[受拉正则化扩展](../../SRC/material/nD/adinaConcrete3D/TensionSoftening.h)。

基础模型迁移 ADINA 模式 5 的小应变混凝土算法：有理压缩曲线、多轴破坏包络、加载/卸载判别、固定正交裂缝、裂缝闭合/张开与压碎释放。基础受拉开裂是脆性的。可选 Gf 受拉软化、压碎后弹性卸载和不同迭代历史更新属于扩展，应与原模型分开描述。

## 2. 当前命令与默认值

```tcl
# N-mm-MPa；负值为压缩。
nDMaterial AdinaConcrete2D 1 30000 0.2 3 -30 -0.002 -6 -0.006

# 五分量壳材料，可直接供 LayeredShell 使用。
nDMaterial AdinaConcrete2DPlateFiber 2 30000 0.2 3 -30 -0.002 -6 -0.006

# 六分量实体材料。
nDMaterial AdinaConcrete3D 3 30000 0.2 3 -30 -0.002 -6 -0.006
```

以上是接口示例，不是试件标定。`ft>0`，`fc<=fu<0`，`epsu<epsc<0`，`E>0`，`-1<nu<0.5`；其他输入限制以工厂和 Parameters 校验为准。

| 选项 | 当前基础模式默认 | 作用与解释 |
|---|---|---|
| `-update committed\|iterative\|hybrid` | **committed** | committed 每次试算从上一收敛态重算；iterative/hybrid 改变迭代期间历史更新，不只是改变切线 |
| `-tangent secant\|algorithmic` | **secant** | 选择迭代矩阵；algorithmic 含各分支数值算法的导数处理 |
| `-maxStrainIncrement value` | **0** | 0 不细分；正数允许子步。这是数值设置，不是材料常数 |
| `-crushedUnloading original\|elastic` | **original** | elastic 为可选压碎后卸载扩展 |
| `-fractureEnergy Gf Lch` | 不启用 | 启用连续受拉软化；Gf、Lch 必须使用一致单位 |
| `-retention normal shear` | `0.0001 0.5` | 开裂后的增量刚度保留，不等价于具有已验证骨料咬合/滑移能力的剪切本构 |
| `-controls beta gamma kappa alpha` | `0.5 1 0.7 0` | 原算法控制常数 |
| `-envelope 24_values` | 内置包络 | 按 SP1/SP31/SP32/SP33 数组顺序输入，不能任意排列 |
| `-rho density` | 0 | 密度 |
| `-G transverseG` | `E/(2*(1+nu))` | 二维工厂支持；用于 PlateFiber 的弹性横向剪切 |

**启用 Gf 时，Tcl/Python 工厂对未显式指定的选项采用 `committed + algorithmic`。** 用户显式选项优先，且不依赖选项先后顺序。直接在 C++ 中创建 Parameters 的调用方应自行设置这些控制，不能假定工厂行为也会自动发生。

```tcl
# Gf=0.1 N/mm，Lch=100 mm；示意参数。
nDMaterial AdinaConcrete2DPlateFiber 10 30000 0.2 3 -30 -0.002 -6 -0.006 \
    -fractureEnergy 0.1 100 -update committed -tangent algorithmic
section LayeredShell 10 3 10 100 10 100 10 100
```

Lch 代表软化方向的特征长度，不能自动取壳层厚度。三维分量为 `[xx,yy,zz,xy,yz,xz]`，二维 PlaneStress 为 `[xx,yy,xy]`，剪应变为工程剪应变。PlateFiber 接口含面内及横向剪切分量；原生二维 PlateFiber 的横向剪切为弹性补充。

## 3. 最重要的历史更正

早期复查曾认为 ADINA 每次平衡迭代在上一迭代态上更新，因此把 iterative 设为默认。第三轮读取原程序 UNBLD/ELEMNT/CDMOD 调用链并对原程序插桩后，确认每次平衡迭代从 IT2 读取同一已收敛历史，试算历史不在该迭代中写回，且每增量没有材料子步。因此当前默认恢复为 committed、无子步、secant。

`iterative` 和 `hybrid` 仍保留供研究。它们可能把 Newton 超调带入历史；hybrid 对无新开裂/压碎事件的点调整加载函数参照，减轻已记录的压应力棘轮。结果可能依赖平衡迭代路径，应单独验证。早期“iterative 更忠实于 ADINA”的文字不再适用。

后续材料修订还涉及二维/三维压碎释放、极限应变、裂缝方向、三维压碎后的活动状态、有限性检查、试算失败回退和历史复制/存盘。原 Fortran 的工程剪切、主方向及固定裂缝边界缺陷有单独说明和修正 oracle；因此不应把“原源码逐点相同”作为所有路径的唯一正确性标准。具体改动编号、条件及补丁见 PROGRAM_DESCRIPTION 和同目录 CHANGES 文件。

## 4. 存盘与相关壳修订

当前材料存盘格式为 v3：二维数据包 118 项、三维 198 项，保存试算/收敛历史和相关矩阵。旧 v1/v2 按兼容路径恢复，旧算法默认值与当前新建材料默认值不能混为一谈。数据库可用性还取决于 Domain、LoadPattern、包装器、截面和 broker；单一材料 sendSelf 通过不代表完整结构重启通过。

LayeredShell 当前分离初始刚度与已损伤试算刚度，影响 ShellMITC4 的 drilling 刚度及重启一致性。ShellNLDKGQ 增加材料/截面失败向分析层传递。这些属于 接口调用链修订（合并前记录：`EXISTING_MATERIALS_AND_SHELL_CN.md`；参见 历史核查（OpenSeesCUDA 仓库 `docs/models/history-local-review/`，未随本仓库迁移）），不是 ADINA 本构理论变更。

## 5. 已有验证记录应如何解读

下表来自已保存记录，本次没有重跑。

| 资料 | 记录与适用范围 |
|---|---|
| 分组审计（合并前记录：`../../tmp/audit-adef/group_summary.json`；参见 历史核查（OpenSeesCUDA 仓库 `docs/models/history-local-review/`，未随本仓库迁移）） | 26 组、7877 个采样、0 个记录失败；有限测试集合 |
| 重启报告（合并前记录：`../../tmp/restart-adef/results.json`；参见 历史核查（OpenSeesCUDA 仓库 `docs/models/history-local-review/`，未随本仓库迁移）） | 8 个算例全部通过；41 步，记录最大缩放差约 8.13e-13；应匹配报告中的 executable SHA |
| 二维 R3 路径报告（合并前记录：`../../tmp/d2-r3/report.json`；参见 历史核查（OpenSeesCUDA 仓库 `docs/models/history-local-review/`，未随本仓库迁移）） | 88 组；包含紧容差一致、较宽容差、阈值平局和已知原 Fortran 多裂缝差异，不能合称全部逐位一致 |
| 三维 R3 路径报告（合并前记录：`../../tmp/d3-r3/report.json`；参见 历史核查（OpenSeesCUDA 仓库 `docs/models/history-local-review/`，未随本仓库迁移）） | 92 条记录：40 原程序一致，7 修正 oracle 一致，16 双方停在包络边界，20 原程序违反固定裂缝不变量，4 原三裂缝关闭缺陷，4 阈值平局，1 存储应变/重主值边界；“全部已分类”不等于 92 条完全相同 |
| PROGRAM_DESCRIPTION 回归统计 | 2D 内核 234848、适配器/截面 6755、RC 膜 3961660；3D 内核 37358、适配器 1980、Tcl 954；Gf 内核 62736、层壳 69094 项记录。检查数不是独立试验数 |

第四轮优选案例显示：2H16-Ref、2H20-Ref 完成，但峰值分别高估约 29.9%、56.2%；VS、Void、Truss、板和多面墙仍提前停止；三维 BS-OA1 完成，BS-A1/BS-C1 仅覆盖部分试验位移。低 NMAE 可能只描述已算的一段。

部分差结果在原 ADINA 的配对有限元模型上同样出现，为模型或建模局限提供支持；这并不能排除所有迁移错误，也不能把不收敛位移直接当作真实破坏位移。墙改用刚性加载梁后全部通过重力，证明加载表示影响很大；该改变也改变了有限元边界模型，应保留原加载方式对照。

## 6. 使用与继续复查的重点

采用原算法比较时显式记录 `committed + secant + maxStrainIncrement 0 + crushedUnloading original`。启用 Gf、iterative/hybrid、elastic 卸载或子步后另起对比组。分层壳板需要检查面内、横向剪切和钢筋贡献，不能仅用弯曲曲线相似来宣称冲切能力正确。

下一轮若修改程序，优先用已有 oracle 和完整结构提交后力平衡检查定位差异；针对不收敛点提取单点历史与活动分支，再区分本构、积分器、截面包装及约束载荷因素。详见 后续复查流程（合并前记录：`COMPARISON_AND_REVIEW_CN.md`；参见 历史核查（OpenSeesCUDA 仓库 `docs/models/history-local-review/`，未随本仓库迁移））。
