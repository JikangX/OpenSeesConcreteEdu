# AdinaConcrete3D 验证记录

> 2026-10-08：这是早期版本的验证记录。当前默认、三维全路径分类与后续优选算例，
> 请先读 [统一 ADINA 说明](../../../docs/models/ADINA_CONCRETE_CN.md)及第四轮说明的 R4/R3 节。
> 早期 iterative 默认已经更正；本次只整理文档。

验证日期：2026-10-05。结果：构建和完整测试脚本均成功，退出码为 0。

OpenSees 基线提交为 `7123c1b150fd458e82a35cbc1fc1018dd283ffcd`，
本次修改位于当前工作目录，未创建提交。已保留用户原有工作目录改动。

## 执行环境及结果

- Windows x64；OpenSees Tcl 3.8.0；CMake Release；MinGW GCC/gfortran。
- 独立数值核心按 C++17、`-O2 -Wall -Wextra -pedantic` 编译。
- OpenSees 使用现有 `bm` 构建树，目标 `OpenSees` 成功。
- 数值核心：**22,583 项标量/状态检查通过**。
- 真实 NDMaterial 和 Channel 接口：**663 项检查通过**。
- 真实 stdBrick、Newton 平衡、空模型数据库重建、参数解析：**954 项检查通过**。
- 合计 **24,200 项检查**；此外完成 **7 条原 Fortran 历程、2,090 个增量点**的对比。

下面的误差单位与输入应力单位一致；本组 E=30000、fc=−30，按 N/mm² 可解释为 MPa。
相对误差以各历程原程序的最大绝对应力归一化。

|历程|增量数|最大绝对应力差|相对峰值误差|判定|
|---|---:|---:|---:|---|
|三维弹性混合应变|20|1.91e−14|1.23e−14|一致|
|侧向约束单轴压缩应变历程|550|6.11e−7|1.25e−8|一致|
|双轴压缩|400|1.76e−12|1.93e−14|一致|
|三轴压缩|300|4.12e−13|6.53e−15|一致|
|循环压缩|220|2.43e−9|9.02e−11|一致|
|顺序三向拉伸|300|3.0167567|1.0|按已说明的旧方向缺陷修正验证|
|单裂缝闭合|300|1.78e−15|3.01e−16|一致|

一致路径要求相对峰值误差 ≤2e−6，实际最大值为 1.25e−8。
顺序拉伸的前 135 步仍须符合一致性阈值；第 136 步开始，原程序在第二裂缝
形成时释放了错误方向的应力。迁入结果须释放正确方向应力、全程满足拉应力截断，
并在最终三方向开裂后应力为零。该历程的差异没有伪装为数值舍入误差。

另外分别验证了三条裂缝中任意一条闭合、20 种三维方向下的开裂/闭合路径，
以及含不同提交/试算裂缝坐标的完整序列化。有关两项旧方向缺陷的修正说明见
[实现说明](README.md)。

## 重启的配套修复

实体单元重启测试发现并修复了原有框架问题：

1. `Domain::recvSelf` 的头长度改为与发送端一致的 17 项；EQ 约束字段使用自身位置。
2. `LoadPattern::recvSelf` 按传入的几何标识读取元数据，避免重建期间标识变化。

测试先保存开裂的实体单元，再继续加载取得参考结果；随后清空模型、重新打开数据库，
由 FEM_ObjectBroker 重建节点、单元和各积分点材料，再加载到相同状态。
应力、九项历史响应及 equalDOF 约束位移均通过校验。

力控压缩还完成了 20 步实际 Newton 平衡求解，每步检查八个积分点的目标应力，
最终检查支座反力平衡。这项测试使用材料切线参与求解。

## 复现命令与结果文件

在 `C:/Work/OpenSees` 执行：

```powershell
$env:PATH='C:\msys64\mingw64\bin;C:\Users\Administrator\.cache\codex-runtimes\codex-primary-runtime\dependencies\python;'+$env:PATH
& C:/Users/Administrator/.cache/codex-runtimes/codex-primary-runtime/dependencies/python/python.exe `
    tests/materials/adina_3d/run_tests.py --adina C:/Work/Adina/all_adina.for `
    --opensees-build bm --opensees bm/OpenSees.exe --build tmp/adinaConcrete3D-verified
```

- 完整日志：`tmp/pdfs/adina3d/verification-final.log`。
- 原程序对比指标：`tmp/adinaConcrete3D-verified/oracle_results.json`。
- 生成的测试可执行文件、只读提取的原子程序及数据库文件均在临时测试目录。
- 最后成功构建日志：`tmp/pdfs/adina3d/build-restart2.log`。

理论手册 SHA256：
`FD8D0D2F874933A2D46AE19BC10B548B376B8A089C96B8A1AAE8C8D6E9995CE1`。

`all_adina.for` SHA256：
`9A0E2731AF05AD4EB76B8FABD890E9492896A296AA3770A560FD7FDE180FFB1B`。

当前验证覆盖小应变、等温三维模型。Python 命令注册已接入公用材料工厂，
本次运行验证使用 Tcl 可执行文件，未单独构建或运行 OpenSeesPy。
测试未包含温度耦合、有限变形、实验标定或结构层面的网格无关性验证；
原模型也未提供断裂能正则化。这些范围不能由上述通过结果推断得到。
