# ADINA 三维混凝土迁入与验证

2026-10-08 文档更正：当前原模式默认 `secant + committed`，
`-maxStrainIncrement 0`（不分子步）、`-crushedUnloading original`。
第三轮原 ADINA 调用链复核已更正早期 iterative 默认。存盘为 v3（198 项），兼容 v1/v2。
断裂能工厂对未显式指定的选项采用 `algorithmic + committed`。
最新证据及使用限制见 [统一 ADINA 说明](../../../docs/models/ADINA_CONCRETE_CN.md)及
[第四轮程序说明](../adina_review/PROGRAM_DESCRIPTION.md)的 R4/R3 节。
下文及早期侧向自由压碎试验记录按其各自参数、版本解释；本次没有修改程序。

本目录验证 `nDMaterial AdinaConcrete3D`。数值实现对应 ADINA 的材料模式 5，
不是已有二维 PSUMAT 的三维包装。数值核心、OpenSees 状态管理和命令解析分别位于：

- `SRC/material/nD/adinaConcrete3D/Concrete3D.{h,cpp}`：独立的 C++17 材料点算法。
- `SRC/material/nD/AdinaConcrete3D.{h,cpp}`：NDMaterial 接口、复制、回退和序列化。
- `SRC/material/nD/AdinaConcrete3DFactory.cpp`：Tcl/Python 公用参数解析。

已完成的构建、测试数量、误差及复现环境见 [验证记录](VALIDATION.md)。

## 理论和原程序对应关系

依据 `C:/Work/Adina/ADINA_ADINAT 使用手册  自动动态增量非线性分析有限元程序.pdf`：
印刷页 105–108（PDF 第 116–119 页）给出混凝土应力应变、强度包络、裂缝和卸载规则；
印刷页 129–132（PDF 第 140–143 页）给出三维材料模式 5 的输入参数。
核对的原程序为 `C:/Work/Adina/all_adina.for` 中的四个子程序。

|原子程序|迁入实现|含义|
|---|---|---|
|CMOD3D|`step`、`replay`、`integrate`|加载判断、增量应力更新、压碎历史|
|PRNCP3|`principalStress`、`strengthEnvelope`、张量旋转函数|主应力/主方向、三轴强度插值|
|CRAK3D|`updateCracks`、`rotateTransverse`|固定裂缝形成、正交裂缝、闭合和重开|
|DCRAK3|`localStiffness`、`releaseStress`|开裂前/后的局部刚度、应力释放|

拉应力、拉应变为正；压应力、压应变为负。六分量顺序为 OpenSees 的
`xx, yy, zz, xy, yz, xz`，后三个应变是工程剪应变 γ，而应力是 τ。
原 ADINA 的剪分量顺序是 `xy, xz, yz`；对照测试在输入输出边界交换这两项。

单轴压缩曲线为

```
σ = E ε / (1 + A x + B x² + C x³),  x = ε / εc
dσ/dε = E (1 − B x² − 2 C x³) / (1 + A x + B x² + C x³)²
r = εu / εc, Es = fc / εc, Eu = fu / εu
A = [E/Eu + (r−2) r² E/Es − (2r+1)(r−1)²] / [r(r−1)²]
B = 2E/Es − 3 − 2A,  C = 2 − E/Es + A
```

多轴受压通过 SP1、SP31、SP32、SP33 插值得到当前压缩强度及峰值应变；
`gamma` 同时控制受约束压缩时应变能力的增长。加载函数为
`sqrt(J2) + alpha * trace(stress)`。加载段采用三点 Gauss 积分计算增量平均模量，
小压应力区使用应力加权的各向同性模量，超过 `kappa` 阈值后使用方向模量。
裂缝固定在形成时的方向；新裂缝与已有裂缝正交。

达到压碎强度时保留预测应力并永久记录压碎状态；达到缩放后的极限压应变时
释放全部应力。不能把“达到峰值强度”实现成“立即清零”。

## 使用

```tcl
model BasicBuilder -ndm 3 -ndf 3
nDMaterial AdinaConcrete3D 1 30000 0.2 3 -30 -0.002 -6 -0.006 \
    -rho 2.4e-9 -retention 0.0001 0.5 \
    -controls 0.5 1.0 0.7 0.0 -maxStrainIncrement 0.00002
# 配合实体单元，例如：element stdBrick tag n1 n2 n3 n4 n5 n6 n7 n8 1
```

示例采用 N、mm、s 单位：模量和强度为 N/mm²，密度为质量/mm³。
Python 共用相同工厂和关键字：

```python
ops.nDMaterial('AdinaConcrete3D', 1, 30000., .2, 3., -30., -.002, -6., -.006)
```

|参数|含义和约束|默认值|
|---|---|---|
|E、nu|初始模量 E>0；−1<nu<0.5|必填|
|ft、fc、epsc、fu、epsu|ft>0；fc≤fu<0；epsu<epsc<0；fc/epsc≤E|必填|
|`-rho`|非负质量密度|0|
|`-retention normal shear`|法向 STIFAC、剪切 SHEFAC，均在 [0,1]|0.0001、0.5|
|`-controls beta gamma kappa alpha`|0<beta<1；gamma>0；0≤kappa≤1；alpha 有限|0.5、1、0.7、0|
|`-maxStrainIncrement`|每个内部子步的最大应变张量 Frobenius 范数增量，>0|2e−5|
|`-envelope` 后 24 个数|SP1 六项，SP31 六项，SP32 六项，SP33 六项|下表|

|曲线|六个默认归一化输入点|
|---|---|
|SP1|0, 0.25, 0.5, 0.75, 1, 1.2|
|SP31|1, 1.4, 1.7, 2.2, 2.5, 2.8|
|SP32|1.3, 1.5, 2, 2.3, 2.7, 3.2|
|SP33|1.25, 1.45, 1.95, 2.25, 2.65, 3.15|

这些默认包络与提供的原程序一致，应使用试验数据校准具体混凝土。
SP1 应从 0 严格递增；输入还必须使插值分母非零、各曲线非递减。
最小压应力的约束比达到或超过 SP1 最后一项时返回错误，不外推。
与旧输入卡的“0 表示采用默认值”约定不同，本接口显式给出的 0 会被保留或
按参数约束拒绝：例如 `-retention 0 0` 表示完全关闭相应残余刚度。

材料可通过实体单元的材料点响应输出：

```tcl
eleResponse 1 material 1 stress
eleResponse 1 material 1 strain
eleResponse 1 material 1 crackState
eleResponse 1 material 1 crackFrame
```

`crackState` 的九项依次为：开启裂缝数、压碎标志、最大加载函数、曾形成裂缝标志、
三个裂缝形成应变、压碎强度、压碎峰值应变。
`crackFrame` 是按行排列的 3×3 矩阵，列向量为裂缝局部坐标轴。

## 现代化改动及明确的差异

1. 数值核心不依赖 OpenSees，不含 COMMON、静态可变缓存、I/O 或全局材料状态。
   参数和两个完整历史状态均使用值语义，各积分点副本独立。
2. 主应力采用对称 Jacobi 特征分解；保留完整九个方向余弦。
   旧程序仅存四个余弦并使用 ±51、100 编码，恢复第三分量时固定取正平方根；
   新实现支持含负方向余弦的任意三维旋转。
3. 第二裂缝形成时，即使横向剪应力为零，也重新排列横向主方向。
   原 CRAK3D 的 `IF (SIGP(6).EQ.0) GO TO 112` 会跳过这个操作。
   在顺序 X/Y 拉伸的第 136 个对照步，原程序 Y 向应力仍为 3.0167567，
   迁入结果正确释放 Y 向应力为零。这项差异单独验证和记录，不计入“旧程序一致”的指标。
   裂缝闭合时同步重排坐标和形成应变，使开启平面始终排列在前；因此三条裂缝中
   任意一条闭合都能恢复该平面的压应力，而不会错误释放另一个平面的应力。
   这是对旧程序“仅减少裂缝计数，不同步交换方向”的修正，分别测试三个闭合方向。
4. 默认 `getTangent()` 返回 ADINA 的正值材料迭代矩阵（`secant`）。显式
   `-tangent algorithmic` 返回实际应力更新的数值算法切线，从同一更新起点重放。
   `-update committed` 以最后提交态为起点，`iterative` 以上一次试算态为起点。
   旧程序供求解器使用的正残余刚度不能直接视为应力映射的导数。
   差分使用相同子步数；光滑分支采用中心差分，开裂、闭合、压碎切换处采用匹配
   实际试算活动分支的单边差分，必要时缩小扰动。非光滑点不能视为具有唯一导数。
   二维、三维共用无 OpenSees 依赖的 `BranchTangent.h` 数值模块，无法解析的试算明确返回失败。
   通常需要 12 次扰动重放，切换点可能更多；大型三维模型需要考虑这项计算成本。
5. 子步数超过 10000、非有限输入、非法材料参数或包络越界均返回错误。
   失败试算不污染先前状态；有理曲线分母通过端点及导数根检查，避免漏检窄奇点。
   子步划分采用应变张量的旋转不变范数（工程剪应变平方的权重为 1/2）；
   闭合判断使用机器精度容差，避免坐标旋转的舍入误差提前触发闭合。
   抗拉阈值比较同样采用共用的应力尺度舍入判据（64 倍机器精度），并保留
   原算法首次开裂严格比较、后续指定分支包含等号比较的区别。新增精确
   抗拉阈值的随机旋转及随后关闭路径，检查应力和开裂历史的坐标一致性。
6. Channel 数据使用 v3 的 198 项数值载荷，保存参数、提交/试算状态、切线、
   求解选项与压碎卸载历史；可读 v1（173 项）和 v2（193 项）。
   接收时检查枚举、有限性和裂缝坐标正交性，校验成功后才替换对象。

数据库重启测试还修复了已有 `Domain::recvSelf` 的头数据不一致：
发送端为 17 项，接收端原为 15 项；EQ 约束计数和数据库标识还误用了 MP 的位置。
现按发送格式接收，并通过包含 equalDOF 约束的空模型数据库重建验证。
同时修复 `LoadPattern::recvSelf` 在添加节点荷载后几何标识递增，导致后续
SP 约束元数据查找错误的问题；所有元数据均按保存时传入的几何标识读取。

当前迁入范围为小应变、等温三维混凝土；速率参数不参与本构。
没有迁入温度依赖的 THERM3、有限变形 Cauchy 应力转换或重力初应力初始化。
默认保留原固定裂缝/压碎模型；可选 `-fractureEnergy Gf Lch` 仅对 Mode-I 拉伸
进行正则化，受压能量正则化未实现，软化分析仍应检查网格敏感性。
OpenSees 现有 [ASDConcrete3D 官方说明](https://opensees.github.io/OpenSeesDocumentation/user/manual/material/ndMaterials/ASDConcrete3D.html)
提供另一种带特征长度正则化的模型，应按分析目的选择。

## 可复现测试

不依赖 pytest，只需要 Python 标准库、C++17 编译器；原程序对照另需 gfortran。

```powershell
# 在仓库根目录运行，确保编译器运行库在 PATH 中。
python tests/materials/adina_3d/run_tests.py --cxx g++ --fc gfortran `
    --adina C:/Work/Adina/all_adina.for

# 已构建 OpenSees 后，增加真实接口和实体单元/数据库测试。
python tests/materials/adina_3d/run_tests.py --adina C:/Work/Adina/all_adina.for `
    --opensees-build bm --opensees bm/OpenSees.exe
```

`--adina` 按子程序边界只读提取用户原文件，在临时测试目录编译；
原程序代码未复制进迁入源文件。原始文件采用 `-fno-automatic`，
对应其依赖跨调用保留局部变换矩阵的旧式存储行为。

|测试层|检验内容|
|---|---|
|独立公式|曲线峰值/极限点、压缩切线、J2、包络节点和受约束峰值应变|
|随机张量|500 组应力；重复主值、重构、正交性、应力应变功共轭|
|材料点历程|六分量弹性、1/2/3 裂缝、闭合、剪切保留、卸载、极限压碎|
|方向与积分|20 个任意三维旋转的完整三向开裂路径；四级子步收敛|
|原 Fortran|7 条历程，2090 个增量点；一致性与已说明的横向方向修正分开判定|
|真实 NDMaterial|提交/回退、重复试算、增量接口、独立副本、Channel 重建及损坏数据拒绝|
|真实 OpenSees|stdBrick 八积分点六分量仿射应变、力控压缩的 Newton 平衡与反力、裂缝响应、数据库重启后续路径、命令错误输入|

测试结果保存在指定 `--build` 目录；默认 `tmp/adinaConcrete3D-tests`。
`oracle_results.json` 记录每条原程序历程的步数和最大误差。
测试成功表明这些公式、路径和接口契约通过验证，不等于已完成结构试验标定。
