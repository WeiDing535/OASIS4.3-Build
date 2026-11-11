# OASIS4.3 精简开发版

NSF-D2 experimental standalone CLI: Rfree37.42% ->31.28%, still above30%. Options remain off by default. [Study report](docs/NSF_D2_IMPROVEMENT_20261026.md).

当前新增 dev5 自带引擎模式：无需另装套件，SEC17 完整隔离复验 Rfree 27.42%。详见 [独立版验证记录](docs/SELF_CONTAINED_PPLUS_20261026.md)。仅供已有授权组件的本机评估。

**两种定相方法 → 对应密度修饰 → 可选 Buccaneer / REFMAC。**

建模精修默认开启，可手动关闭；关闭后成功输出最终相位/FOM 与分别保存的定相结果和纯观测 HL，无需配置外部建模精修运行资源。

- 连续 SAD + MLHL-best 密度修饰。
- 直接法 P+ + Legacy 密度修饰。
- 统一输入、明确参数、五阶段进度、逐轮精修统计及分组结果。

双击 `OASIS4.3.command` 打开原生桌面窗口，或在 Python 3.12 专用 conda 环境中执行 `python -B run_desktop.py --project .`。

界面用 Tk（随 Python 运行时一起提供）绘制，不需要浏览器、不需要端口转发、不联网；科学计算组件全部保留在包内。窗口第一页是“任务记录”，保存跑过的每个任务并可回看日志与终点；多个任务可同时运行（各自独立输出目录、日志与进程组，上限 4 个）。“新建任务”清空表单并建议新目录，不重要的参数默认收在“显示高级参数”里，“导出结果…”复制整个输出目录，“结构 3D”页可直接旋转查看最终模型的 Cα 主链。可选的辅助程序（旧版 Qt 界面、绘图等）不属于计算核心，精简包内不再附带，请按[使用说明](docs/OASIS4.3.md)用 `pip` / `conda` 自行安装。

完整输入要求、所有选项、结果说明、CLI 恢复与验证范围见 [OASIS4.3 使用说明](docs/OASIS4.3.md)。

连续 SAD 的定相后验包含最终三重相位关系与部分结构场；自 4.3.0.dev2 起另存纯观测似然 HL 供 REFMAC 使用，DM 仍使用原定相结果。P+ 的观测 HL 数值保持不变。详见使用说明。恢复基线与备份见 [开发恢复报告](development/RECOVERY_REPORT_20261023.md)；最新修正、真实数据对照及模型局限见 [胰岛素诊断报告](docs/INSULIN_DIAGNOSIS_20261023.md)。

P+ 的可选研究模块、四组消融试验及补充先验对照见 [P+ 优化记录](docs/PPLUS_OPTIMIZATION_20261023.md)。该模块尚未接入 GUI 或默认流程；50 残基候选仍缺原子，并存在数值敏感性，不能视为稳健的默认算法升级。

专用精简源码包由 `python -B package_oasis.py --output /new/path.zip` 生成；只收录此流程所需代码与资源，排除旧产品入口、旧界面、训练脚本和训练模型。源码包不包含 Buccaneer/REFMAC 运行时，当前也不是新的完整安装包。

从源码树直接运行时，P+ 流程默认用 CCP4 Parrot 做密度修饰；找不到 `cparrot` 时不再静默改道：命令行把回退原因写到标准错误，界面在工具栏下显示红字警告，结果页把密度修饰引擎标为“内置 RESOLVE/MLHL（Parrot 不可用，已回退）”，`pipeline_report.json` 里保留 `dm_engine_requested`、`dm_engine_fallback` 两个字段。自带完整运行时内含 Parrot，不受此影响。

“免安装”指的是**完整包**：`OASIS4.3-dev17-portable` 内自带 `runtime/`（Phaser、Parrot、Buccaneer、REFMAC 与 OASIS 密度修饰），双击即用，不需要另装 CCP4 或 Phenix。源码树不含这些二进制，从源码运行时要由本机提供：macOS 上可指向已安装的 CCP4 9（`/Applications/ccp4-9`），或把完整包的 `runtime/` 放到源码树旁并用 `--runtime` 指定；Linux 目前还没有完整运行时。本仓库记录的十案例基准是在自带运行时下跑出的，因此源码树缺 Parrot 时的那次运行只能算内置引擎的参考结果，不能与基准直接对比。

完整包目前没有公开下载地址（第三方引擎许可）；怎样用本机已有的 CCP4 组装出等价的私有免安装包、以及跑通十案例协议的完整命令，见 [免安装包获取说明](docs/INSTALL_BUNDLE.md)。

为什么现在还绕不开 CCP4——四个模块（Phaser、Parrot、REFMAC、Buccaneer）各自能自研到什么程度、实测差距是多少、发布路线怎么分三档，见 [CCP4 替换评估](docs/NO_CCP4_REPLACEMENT.md)。其中定相与密度修饰已有仓库内实现（`--method pplus_internal`），精修有可调用的原型（`run_refine.py --engine internal`），建链仍只有未接入的原型。

现存 Intel Mac 0.8.3 完整包保持原样。Linux 完整运行时仍待提供。开发仓库保留历史模块用于回溯与测试；本版产品不提供其他独立功能入口。旧版说明存档见 [历史 Build 说明](docs/LEGACY_BUILD_README_20261023.md)。

原始输入、原自由集、定相结果、独立观测 HL、密度修饰结果和每轮模型分别保留；固定最后一轮结果。算法属于需实验验证的研究实现，真实数据完整流程效果尚待验证。

新增三组实验数据、六组密度修饰追加试验、数值修正及全部正负结果见 [P+ 多案例验证报告](docs/PPLUS_MULTICASE_VALIDATION_20261024.md)。研究方案尚未升级为默认流程。

后续发现旧 Parrot 对照实际读取纯观测 HL，不能将其改善归因于 P+ 后验。更正与 18 组相同先验配对结果见 [DM 先验交接审计](docs/DM_PRIOR_AUDIT_20261024.md)。

P+ 初始图与原始 SAD 双峰约束的分离、条件分支种子及三例固定预算结果见 [双通道 P+ 验证](docs/DUAL_CHANNEL_PPLUS_20261024.md)。有部分改善，尚未证明稳定优于纯观测对照，默认流程未改动。
