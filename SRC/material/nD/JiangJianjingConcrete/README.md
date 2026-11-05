# 江见鲸混凝土模型（JiangJianjingConcrete）——教学用非线性弹性全量模型

来源：清华大学研究生课程《钢筋混凝土有限元》第七讲“混凝土非线性弹性本构关系”（陆新征）
及其 MSC.Marc 用户子程序 `JiangJJ2007.for`（HYPELA2）。本目录把该教学程序移植为 OpenSees 的
三维 `nDMaterial`，可直接用于实体单元，经 `PlaneStress` / `PlateFiber` 包装用于平面应力单元和分层壳。

**这是教学模型**：没有开裂、受拉软化和卸载规则，达到破坏面后应力冻结，实际构件算例中通常很早就无法继续
（见第 6 节）。请用于理解“非线性指标 + 等效单轴关系 + 割线刚度”的全量模型框架，不要用于工程分析。

---

## 1. 理论（课件第 25–39 页）

1. **破坏准则**（江见鲸四参数准则，ft/fc = 0.1）
   a·J2/fc² + (b + c·cosθ)·√J2/fc + d·I1/fc − 1 = 0，a = 1.2856，b = 1.4268，c = 10.2551，d = 3.2128；
   I1 为第一应力不变量，J2、J3 为偏应力不变量，cos3θ = 4J3/r³，r = √(4J2/3)。
2. **非线性指标（√J2 法）**：保持 I1、θ 不变，求破坏面上的 √J2f（上式对 √J2f 的一元二次方程正根），
   β = √J2 / √J2f。
3. **等效单轴关系**（中国规范二次抛物线 σ = fc[2ε/ε0 − (ε/ε0)²]，E0 = 2fc/ε0）：
   割线模量 Es = E0/2·(1 + √(1−β))，切线模量 Et = E0·√(1−β)；割线泊松比 νs = ν0（程序中关闭了课件第 31 页的
   νs 随 β 增大公式）。
4. **全量本构**：σ = D(Es, νs)·ε（各向同性割线矩阵，课件第 37 页）。
5. **程序附加规则**（`Get_DS`）：
   - β 由**增量开始时的应力**计算（显式），β > 1 时取 1；
   - β ≥ 0.99999 时应力保持为增量开始时的值（“达到极限应力后应力不变”）；
   - 压碎：εm = min(平均应变, 0)，若 min(εii − εm) < εcu（默认 −0.0033），应力清零、刚度取 1e−6·D(E0)，不可恢复；
   - 返回给求解器的迭代矩阵为 D(Et, ν) + 1e−6·D(E0, ν)。

## 2. 命令

```tcl
nDMaterial JiangJianjingConcrete $tag $fc $E0 $nu <$epsCrush> <-criterion $a $b $c $d> <-tangent original|secant> <-beta committed|current> <-rho $rho>
```

| 参数 | 含义 | 默认 / 说明 |
|---|---|---|
| `fc` | 单轴抗压强度，**输入正值**（与原程序 `C%fc=30.` 一致） | 必需 |
| `E0`、`nu` | 初始弹性模量、泊松比 | 必需 |
| `epsCrush` | 压碎应变（负值） | −0.0033 |
| `-criterion a b c d` | 破坏准则常数 | 1.2856 1.4268 10.2551 3.2128（ft/fc = 0.1） |
| `-tangent` | `original`：D(Et)+1e−6 D(E0)，同原程序；`secant`：D(Es)，本算法应力更新的**精确导数** | original |
| `-beta` | `committed`：β 取自增量开始（上一收敛步）的应力，同原程序；`current`（2026-10-07 新增，非原程序）：β 取自正在计算的应力，即课件全量割线关系 σ = D(Es(β(σ)))·ε 在每次调用中精确求解（沿 D(E0)ε 方向的一元方程），破坏前与路径和步长无关；试算越过破坏面时取该方向上破坏面上的应力，下一增量起与原程序同样保持；压碎不变。与 `-tangent secant` 配合时迭代矩阵为该更新的精确（非对称）导数 | committed |
| `-rho` | 密度 | 0 |

单位自洽（例如 N、mm、MPa 或 N、m、Pa）。应变 Voigt 顺序 [xx, yy, zz, xy, yz, zx]，工程剪应变。

用于二维与分层壳：
```tcl
nDMaterial JiangJianjingConcrete 1 30.0 30000.0 0.18 -tangent secant
nDMaterial PlaneStress 2 1        ;# quad ... PlaneStress 2
nDMaterial PlateFiber  3 1        ;# section LayeredShell ... 3 thickness ...
```

输出：`eleResponse $e material $ip state`（或 `nonlinearIndex`）返回
crushed、β、βmax、Es、Et、indexSolvable；以及标准的 `stress`、`strain`、`tangent`。存盘/重启支持。

### 求解建议

- 原程序的切线 Et 不是本算法的导数：β > 8/9（Es/Et > 2）后完整 Newton 发散。用以下任一设置：
  `algorithm ModifiedNewton -initial`（初始刚度，收敛因子 1 − Es/E0 ≤ 0.5，始终收敛），或 `-tangent secant`
  配合 `algorithm Newton`。`PlaneStress`/`PlateFiber` 包装器内部也用材料切线迭代，**与包装器配合时必须用
  `-tangent secant`**（单元平面应力试验中 original 只能到 21 MPa，secant 到 30.06 MPa）。
- β 取自上一收敛步，结果与步长有关；应变控制下峰值可超过 fc 一个增量（试验中 30.22 MPa）。

## 3. 文件

| 文件 | 内容 |
|---|---|
| `JiangJJ2007.for` | 原 MSC.Marc 程序（未改动） |
| `JiangJianjing.h/.cpp` | 与 OpenSees 无关的可重入材料点内核（`nonlinearIndex`、`crushes`、`elasticStiffness`、`update`） |
| `JiangJianjingConcrete.h/.cpp` | OpenSees 适配器（命令解析、试算/提交/回滚、存盘、响应） |
| 注册 | `classTags.h`（ND_TAG_JiangJianjingConcrete 7024）、`FEM_ObjectBrokerAllClasses.cpp`、Tcl `TclModelBuilderNDMaterialCommand.cpp`、Python `OpenSeesNDMaterialCommands.cpp`、新运行时 `nDMaterial.cpp`、`CMakeLists.txt` |
| 测试 | `tests/materials/jiang_jianjing/`：`run_tests.py`、`oracle_driver.f90`、`kernel_vs_fortran.cpp`、`kernel_checks.cpp`、`tcl_checks.tcl`；例子：`EXAMPLES/JiangJianjingConcrete/` |

## 4. 与原 Fortran 的对应和差异

| 原程序（Get_DS） | C++ | 说明 |
|---|---|---|
| `Get_D` | `elasticStiffness` | 运算顺序相同 |
| 压碎判别（增量开始应变） | `crushes(committed.strain)` | 相同（因此压碎在越过阈值的下一增量才生效） |
| 应力不变量、θ、β | `nonlinearIndex` | 相同；**D1：J2 截断为 ≥ 0**（见下） |
| 割线/切线模量、应力、保持规则 | `update` | 相同 |
| 迭代矩阵 | `update(..., tangent)` | 默认相同；可选 `-tangent secant` |
| `My_ConA/My_ConB`、`New_Step`/`Cut_Step`（Marc 中保存/恢复增量开始状态） | 适配器的 committed/trial、`commitState`/`revertToLastCommit` | 每次试算从收敛态计算，与 HYPELA2 的调用方式相同 |

发现的问题（迁移时已处理并记录）：

1. **原程序缺陷 D1（静水/近静水应力）**：偏应力的 J2 可因舍入为 −1e−30 级负数，`sqrt` 得 NaN，判别式检查失败，
   程序写出 “Wrong” 并令 β = 1，应力从此冻结。三轴等压路径在第 2 步即“破坏”（应力停在 −2.875 MPa 不再增加）。
   C++ 把 J2 截断为 0（精确值）。
2. **单精度常数**：`1.2856` 等在 Fortran 中是默认 REAL*4 常数，以约 1e−8 的相对误差参与双精度运算。C++ 默认使用
   精确的十进制常数；对比测试中按 Fortran 的实际取值逐步比较。
3. **matmul 求和顺序**：libgfortran 的 `matmul` 与顺序求和相差 1 ulp（Marc 所用的 Fortran 运行库又可能不同）。
4. **零增量再计算**：β 取自增量开始应力，若在已收敛应变处再算一次，应力会改变（用新的割线模量乘同一应变）。
   OpenSees 会在不变位移下重算单元（提交后、数据库恢复、查询），适配器在试算应变恰等于收敛应变时直接返回收敛态。
5. **课件第 39 页**：yz 分量结果印为 2.07，但按其给出的应变增量（2×10⁻⁵）应为 7.20×10³ × 2.98×10⁻⁴ = 2.15；
   其余分量与课件一致（β = 0.5093，√J2f = 8.998，Es = 17.00 GPa，σ = −6.53, −6.38, −12.69, 2.07, 1.07）。

## 5. 验证结果（`python tests/materials/jiang_jianjing/run_tests.py --build <新目录>`）

| 测试 | 结果 |
|---|---|
| 与原 Fortran（未改动，按 HYPELA2 方式驱动）逐步对比 | 45 条应变路径（单轴/双轴/三轴等压/围压后加载/受拉/纯剪/一般剪切/卸载再加载/6 条随机路径 × 3 组参数），4140 步中 3494 步逐位相同；破坏前（β < 0.98）最大应力差 4.9e−8·fc；压碎标志全部一致 |
| 差异来源 | 破坏前差异只出现在轴对称路径（cos3θ = ±1 处 dθ/d(cos3θ) 无穷大，1 ulp 的输入差导致 √eps 级差异）；破坏后 √(1−β) 与 0.99999 阈值使舍入放大，最大 1.8e−2·fc；6 条静水/围压路径触发原程序缺陷 D1（未修正的原程序差 3·fc） |
| 内核检查 | 90 项：课件例题、单轴抗压强度 = fc、单轴抗拉强度 = 0.1fc、β 旋转不变、首步线弹性、压碎、可重入、`-tangent secant` 等于有限差分导数（2.4e−11·E0）而原切线差 0.215·E0、参数校验 |
| OpenSees Tcl 检查 | 31 项：命令解析与参数校验；单个 stdBrick 单轴；应力控制立方体（见下表）；PlaneStress 包装 quad；PlateFiber + LayeredShell + ShellNLDKGQ；数据库存盘/恢复/续算逐位一致 |

应力控制立方体（fc = 30 MPa，原切线 + ModifiedNewton -initial，荷载步 0.05 MPa，取最后收敛值），与课件第 77 页 Marc 结果对比：

| 工况 | OpenSees | 课件（Marc） |
|---|---|---|
| 单轴 | 29.50 | 30 |
| 双轴 σ2 = σ3 | 35.65 | 36 |
| 围压 2 MPa | 43.55 | 44 |
| 围压 6 MPa | 67.00 | 66.6 |

原切线 + 完整 Newton：单轴只到 20.65 MPa（β > 8/9 后发散）；`-tangent secant` + Newton：29.70 MPa。

## 6. 构件算例（教学演示：效果不好是预期的）

（以下为原开发仓库中的构件算例结果，所需试验数据与驱动脚本未包含在本仓库中。）

加载和严格驱动（力容差 1e−4×试验峰值，失败二分，不跳步），只替换混凝土为
`JiangJianjingConcrete fc E0 nu -tangent secant`。覆盖率 = 计算到的挠度 / 试验最大挠度。

| 算例 | 单元 | 停止位置（覆盖率） | 已算最大力 / 试验峰值 |
|---|---|---|---|
| 2H16-Ref / 2H20-Ref 梁 | 平面应力 quad | 0.16 / 0.20 | 0.72 / 0.87 |
| B07 / B08 / B13 / B14（Bresler–Scordelis） | 平面应力 quad | 0.20 / 0.09 / 0.10 / 0.10 | 0.64 / 0.42 / 0.34 / 0.47 |
| VS-C3 / VS-OA1 | 平面应力 quad | 0.06 / 0.06 | 0.38 / 0.45 |
| 深梁 W1 | 平面应力 quad | 0.25 | 1.69（初始刚度为试验 3.6 倍） |
| W28 / W31 墙 | ShellNLDKGQ + LayeredShell | 0.09 mm / 0.20 mm 处停止 | — |
| BS-B1 / C1 / A1 / OA1 | 三维 stdBrick | 2.2 / 2.9 / 1.8 / 2.0 mm（13–26%） | 151 / 126 / 164 / 187 kN（试验 443 / 311 / 467 / 334） |

所有构件都在试验挠度的 6–26% 处停止。B07 停止时 592 个积分点中 57 个达到 β ≥ 0.99999 并冻结应力，其中 54 个在受拉区：
模型没有开裂和受拉软化，受拉区一达到 ft（= 0.1fc）就整点冻结（各方向刚度只剩 1e−6·E0），无法把拉力转给钢筋。
这正是课件所说全量非线性弹性模型的局限（“简单、难以模拟加卸载”），不是迁移错误。

## 7. 复现

```
# 只依赖 C++17 编译器的内核检查（含课件例题）
g++ -std=c++17 -O2 -ISRC/material/nD/JiangJianjingConcrete tests/materials/jiang_jianjing/kernel_checks.cpp \
    SRC/material/nD/JiangJianjingConcrete/JiangJianjing.cpp -o kernel_checks && ./kernel_checks
# 适配器检查：与 OpenSees 的 NDMaterial / PlaneStress / PlateFiber 等类链接（28 项）
#   源文件列表见 tests/materials/jiang_jianjing/adapter_checks.cpp 文件头
# 完整验证（需要 gfortran；--opensees 为编译好的 OpenSees 可执行文件，可选）
python tests/materials/jiang_jianjing/run_tests.py --build <新目录> --opensees <OpenSees>
```
