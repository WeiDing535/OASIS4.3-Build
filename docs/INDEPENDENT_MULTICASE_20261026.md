# 自带引擎版其他案例实测

先测试已交付 dev5。P9 因 ANISOU 读取不兼容失败后，追加完整各向异性张量支持的兼容修正，并以新候选包从输入重新复跑；原失败保留为单独一行。全部试验均未重新划分自由集、按中间 Rfree 选轮次或借用旧蛋白模型。所有病例均通过 sandbox-exec 屏蔽 /Applications、原科学环境、旧引擎包和研究时复制的 SDK，并禁止网络；计算只调用新包内第三方引擎。输入为已有实验 MTZ、序列和重原子子结构，不包含从零搜索重原子的测试。

| 案例 | 分辨率 Å | 最终 Rwork | 最终 Rfree | 原始第5轮 Rfree | Cα / 链 | 自由反射 | 清理 |
|---|---:|---:|---:|---:|---:|---:|---|
| p9 | 1.744 | — | — | — | — / — | 1071 | failed |
| nsf_d2 | 2.400 | 25.68% | 37.42% | 37.42% | 300 / 8 | 642 | refused |
| insulin | 1.701 | 21.70% | 24.94% | 24.94% | 51 / 2 | 422 | refused |
| p9_anisotropic_fix | 1.744 | 24.47% | 27.35% | 27.35% | 138 / 1 | 1071 | refused |

预算与 SEC17 成功分支相同：P+ KL=0.0005、最多100次、1440点积分；Wang14 + classic40 (10/20/10)；固定五轮建模，每轮十周期 REFMAC。仅唯一完整精确序列链允许清理后固定三十周期振幅精修，权重0.005；其他情况保留第五轮，未隐藏失败或拒绝清理的结果。

P9 的实验波长0.9792 Å来自冻结的官方教程说明，胰岛素采用 MTZ 的 CuKa 数据集记录1.54 Å。NSF-D2 沿用此前基准的名义0.9792 Å，尚未独立核实真实实验波长，结果应带此条件解读。P9/NSF-D2 的自由集原先因 Scalepack 未提供而按固定 HKL 哈希划分，此次完整保留；胰岛素保留公开 MTZ 原标记。

这些是既有开发数据，而非全新盲测；预处理和后续模型地图未认证为全程仅工作集，Rfree 不代表严格端到端独立盲测。不同数据分辨率不同，不能仅凭 Rfree 跨样品排名。Cα/已命名残基只是输出计数，未与真实坐标比对，不能当作正确覆盖率。

包内包括来自已有授权发行版的 Phaser、RESOLVE、CCTBX、Buccaneer、REFMAC，用户无需另装 Phenix/CCP4。不是全部算法由 OASIS 原创；仍仅供本机授权评估，未解决公开再分发条件。

逐例记录：

## p9

状态：failed；异常：Stage process exited with code 1。
清理：None；原因：None。
P+ 更新：`{}`。
原观测审计：0 组精修输入 F/SIGF/FreeR_flag 逐 HKL 一致；初始输入与完成阶段文件哈希通过。

## nsf_d2

状态：completed；异常：None。
清理：refused；原因：Cleanup requires exactly one complete, exact full-sequence chain。
P+ 更新：`{"iterations": 57, "converged": true, "mean_abs_probability_change": 0.012694915785104148, "original_observation_HL_modified": false, "posterior_HL_exported": false}`。
原观测审计：5 组精修输入 F/SIGF/FreeR_flag 逐 HKL 一致；初始输入与完成阶段文件哈希通过。

## insulin

状态：completed；异常：None。
清理：refused；原因：Single-copy cleanup requires exactly one sequence record。
P+ 更新：`{"iterations": 40, "converged": true, "mean_abs_probability_change": 0.012501696155145453, "original_observation_HL_modified": false, "posterior_HL_exported": false}`。
原观测审计：5 组精修输入 F/SIGF/FreeR_flag 逐 HKL 一致；初始输入与完成阶段文件哈希通过。

## p9_anisotropic_fix

状态：completed；异常：None。
清理：refused；原因：Cleanup requires exactly one complete, exact full-sequence chain。
P+ 更新：`{"iterations": 43, "converged": true, "mean_abs_probability_change": 0.013818261891618177, "original_observation_HL_modified": false, "posterior_HL_exported": false}`。
原观测审计：5 组精修输入 F/SIGF/FreeR_flag 逐 HKL 一致；初始输入与完成阶段文件哈希通过。

## P9 兼容性修正验证

修正保留 Phaser 输出的完整 Cartesian U 张量，随每个对称算子旋转，使用 exp(-2*pi^2*q.U.q) 衰减；没有投影为各向同性，也没有修改原始 Phaser HL。各向同性分支直接调用原函数，维持原数值路径。561 项测试通过（另1项按原配置排除）；新增8项测试包括 P1/P4/I4、非正交 P21/P31、无效张量拒绝和各向同性逐位不变，并与独立 Gemmi 结构因子实现核对。P9 的全部工作反射另见 P9_ANISOTROPIC_FACTOR_AUDIT.json。SEC17 新候选包只复跑相位/DM，原始似然、P+ 种子和 DM 地图逐位不变；这不是重新构建 SEC17 的新 Rfree 测试。

机器可读：AUDITED_RESULTS.json、RESULTS.csv；固定协议：PROTOCOL.json；各病例 COMMAND.json、RUN.log、run/pipeline_report.json。原发布包和旧案例结果未修改；源代码新增各向异性重原子结构因子模块及读取接口修正，修改前快照和独立差异文件已保存。

## 本轮交付与版本区别

原 dev5 在 P9 上于 P+ 开始前失败，不具有该例 Rfree；27.35% 来自新增各向异性支持的独立候选目录 `/Volumes/Fast_4T/GPT_WORK/OASIS4.3-Build/dist/OASIS4.3-dev5a-portable local`。候选包不需要额外安装；原 dev5 目录和 ZIP 保留。候选修正只改变各向异性重原子参数支持；胰岛素和 NSF-D2 的本轮结果来自原 dev5。

本轮结论：胰岛素 24.94%、P9（兼容修正后）27.35%，均低于30%；NSF-D2 37.42%，且300个 Cα/8条链超过预期273个残基/1条链，仍有碎片与多余模型，不能当成完整正确求解。P9 输出138个残基/1条链，序列预期139；清理按规则拒绝，未追加精修以挑低 Rfree。胰岛素输出51个残基/2条链，符合计数预期，但未与真实坐标比对。
