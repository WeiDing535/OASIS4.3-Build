# OASIS4.3 精简开发版

**两种定相方法 → 对应密度修饰 → Buccaneer → REFMAC。**

- 连续 SAD + MLHL-best 密度修饰。
- 直接法 P+ + Legacy 密度修饰。
- 统一输入、明确参数、五阶段进度、逐轮精修统计及分组结果。

双击 `OASIS4.3.command`，或在 Python 3.12 专用 conda 环境中执行 `python -B run_gui.py --project .`。

完整输入要求、所有选项、结果说明、CLI 恢复与验证范围见 [OASIS4.3 使用说明](docs/OASIS4.3.md)。

专用精简源码包由 `python -B package_oasis.py --output /new/path.zip` 生成；只收录此流程所需代码与资源，排除旧产品入口、旧界面、训练脚本和训练模型。源码包不包含 Buccaneer/REFMAC 运行时，当前也不是新的完整安装包。

现存 Intel Mac 0.8.3 完整包保持原样。Linux 完整运行时仍待提供。开发仓库保留历史模块用于回溯与测试；本版产品不提供其他独立功能入口。旧版说明存档见 [历史 Build 说明](docs/LEGACY_BUILD_README_20261023.md)。

原始输入、原自由集、原始实验 HL、密度修饰结果和每轮模型分别保留；固定最后一轮结果。算法属于需实验验证的研究实现，真实数据完整流程效果尚待验证。
