# ADINA 二维混凝土及分层壳接口

2026-10-08 文档更正：当前原模式默认
`-tangent secant -update committed -maxStrainIncrement 0 -crushedUnloading original`。
第三轮对原 ADINA 调用链的复核已更正此前的 iterative 默认。显式开启断裂能时，
工厂默认切换为 `algorithmic + committed`，显式指定的选项优先。
存盘为 v3（118 项），兼容 v1/v2。最新测试与限制见
[当前统一说明](../../../docs/models/ADINA_CONCRETE_CN.md)及
[第四轮程序说明](../adina_review/PROGRAM_DESCRIPTION.md)的 R4/R3 节；
下文早期验证说明保留历史背景，本次没有修改程序。

本实现迁移 ADINA 模式 5 的原生小应变**平面应力**算法（A06.FOR 中的
`CDMOD / PRNCPL / CRAKID / DCRACK`，`ITYP2D >= 2`）。理论依据为
`C:/Work/Adina/ADINA_ADINAT 使用手册  自动动态增量非线性分析有限元程序.pdf`
印刷页 105–108，对应 PDF 物理页 116–119。

模型保留有理压缩曲线、多轴压缩破坏包络、加载函数、加权各向同性/正交异性切换、
固定正交裂缝、独立闭合/重新张开、压碎和极限应变释放应力。
具体核对结果和适用边界见 [VALIDATION.md](VALIDATION.md)。

## 命令与符号

```tcl
nDMaterial AdinaConcrete2D tag E nu ft fc epsc fu epsu \
    <-G transverseG> <-rho density> <-retention normal shear> \
    <-controls beta gamma kappa alpha> <-envelope 24_values> \
    <-maxStrainIncrement value>

# 可显式创建五分量材料；输入参数相同。
nDMaterial AdinaConcrete2DPlateFiber tag E nu ft fc epsc fu epsu <options>
```

- 拉应力/拉应变取正；`fc, epsc, fu, epsu` 为负。
- `E>0`，`-1<nu<0.5`；`fc<=fu<0`，`epsu<epsc<0`，`fc/epsc<=E`。
- `ft>0`；`beta` 在 (0,1)，`gamma>0`，`kappa` 在 [0,1]。
- `-retention` 默认 `0.0001 0.5`，两者均在 [0,1]。
  它们是开裂后的增量刚度保留系数。开裂法向应力释放为零；剪切保留系数
  不大于 0.001 时也释放面内裂缝剪应力。
- `-controls` 默认 `0.5 1.0 0.7 0.0`，分别对应 BETA/GAMA/KAPA/ALFA。
- `-envelope` 按 SP1(6)、SP31(6)、SP32(6)、SP33(6) 顺序输入无量纲应力比。
  默认值与 ADINA 源程序一致；参数校验拒绝非单调、分母奇异或有理曲线奇异的数据。
- `-maxStrainIncrement` 默认 `2e-5`，约束应变张量增量的 Frobenius 范数，
  即 `sqrt(de_xx²+de_yy²+dg_xy²/2)`；单次试算最多 10000 子步。
- `-G` 默认 `E/(2(1+nu))`，必须正且有限；用于壳层两个横向剪切分量。
  该选项也可在三分量材料上设置，以供分层壳创建副本时使用。
- `rho` 是体积质量密度，默认零；分层截面的面密度为各层 `rho*厚度` 之和。

所有参数均按用户选择的一致单位制输入。零值在本实现中是明确的参数值，
不会采用旧卡片输入中“零代表默认值”的约定。

## 分层壳用法

```tcl
model BasicBuilder -ndm 3 -ndf 6
nDMaterial AdinaConcrete2D 1 30000 0.2 3 -30 -0.002 -6 -0.006 -G 12500

# 直接使用二维材料，截面会请求它的 PlateFiber 副本。
section LayeredShell 10 6 1 0.12

node 1 0 0 0
node 2 1 0 0
node 3 1 1 0
node 4 0 1 0
element ShellMITC4 1 1 2 3 4 10
```

三分量接口是 `[eps_xx, eps_yy, gamma_xy]`；五分量接口是
`[eps_xx, eps_yy, gamma_xy, gamma_yz, gamma_xz]`。
所有剪应变均为**工程剪应变**，剪应力为张量分量。
当前 `LayeredShellFiberSection` 的第 7/8 个截面变形为 xz/yz 剪切，因此层材料
最后两个分量会交换；回归测试逐项验证了这个顺序。

各层使用 `eps_layer=eps_membrane-z*curvature`。材料仅处理壳局部平面内的本构，
空间壳基底旋转由壳单元处理；固定裂缝角度在各自局部平面内保存。
横向剪切是独立的弹性扩展 `tau=G*gamma`，没有从二维理论推导三维剪切破坏。
剪切修正系数如需采用，应由选定壳/截面约定或显式 G 反映，材料不额外重复施加。

钢筋混凝土可将 `PlateRebar` 层与混凝土层交替放入 `LayeredShell`。
混凝土采用 ADINA 原始固定裂缝机制，钢筋采用独立材料，不需要 MCFT 专用的钢筋状态反馈。
完整的膜、弯曲、剪切、钢筋壳循环和重启案例在 [shell_tests.tcl](shell_tests.tcl) 中。

## 模块与状态

| 模块 | 职责 |
|---|---|
| `SRC/material/nD/adinaConcrete2D/Concrete2D.h/.cpp` | 无 OpenSees 依赖的三分量积分核心、张量变换、平面应力凝聚、裂缝状态机、子步和切线 |
| `adinaConcrete3D/Concrete3D.h/.cpp` 的纯函数 | 两种维度共用的已校验参数、有理压缩标量曲线和破坏包络；二维从不调用三维积分器 |
| `AdinaConcrete2D.h/.cpp` | NDMaterial 生命周期、三/五分量接口、独立剪切状态、拷贝、响应及 Channel 序列化 |
| `AdinaConcrete2DFactory.cpp` | 参数解析与工厂，供 Tcl、Python 共用解释器和新版 runtime 调用 |

试算起点由 `-update iterative|committed` 指定；原模式默认以上一次 Newton 试算态
为起点，committed 模式从最后提交态重放。失败保持旧试算和切线不变。
每个对象保存自己的应力、应变、矩阵和历史，无材料级共享可写静态数组。
`commitState/revertToLastCommit/revertToStart` 同时管理面内历史和两个横向剪应变。
副本保持完整状态但获得独立 Channel 身份；支持已提交与尚未提交状态的重启。

默认切线为 ADINA 的材料迭代矩阵（secant）；显式 algorithmic 选项返回完整
应力更新映射的数值导数，每个扰动都从同一更新起点重放，并固定子步数。
这包括包络变化、刚度切换与裂缝状态；软化分支可以产生负或非对称切线。
光滑分支采用中心差分。开裂、闭合或压碎切换处采用与实际试算活动分支一致的单边差分，
避免把有限应力突降除以微小扰动误当作巨大刚度；必要时缩小扰动，无法解析则明确拒绝试算。
二维、三维共用无 OpenSees 依赖的 `adinaConcrete3D/BranchTangent.h` 数值模块。
结构算例采用通用矩阵求解器与足够小的增量；这项切线修复没有引入应力平滑或断裂能参数。

额外响应：

```tcl
eleResponse 1 section 1 fiber 1 stress
eleResponse 1 section 1 fiber 1 strain
eleResponse 1 section 1 fiber 1 crackState
```

`crackState` 的 10 项依次为：open1、open2、crushed、maximumLoading、hasCrackFrame、
angleRadians、crackStrain1、crackStrain2、crushingStress、crushingStrain。
两个开闭标志对应固定裂缝坐标轴，不用裂缝数量推断方向。

## 复现验证

```powershell
$env:PATH='C:/msys64/mingw64/bin;'+$env:PATH
python tests/materials/adina_2d/run_tests.py `
    --adina C:/Work/Adina/all_adina.for `
    --opensees-build bm --opensees bm/OpenSees.exe `
    --build tmp/adinaConcrete2D-verified
```

只给 `--adina` 可运行独立 C++ 核心及原程序对照；`--opensees-build` 加入真正
OpenSees 对象/分层截面的 C++ 测试；`--opensees` 加入 Tcl 壳单元、钢筋层与数据库测试。
测试脚本使用 Python 标准库、C++17 编译器和 gfortran，不依赖 pytest。
原 ADINA 程序只读；原例程和诊断修正版本均临时生成，不随迁入代码分发。
