# 本轮结论

**尚未达到所有测试案例 Rfree 低于或接近 30% 的目标。** 本轮完成 25 组固定终点建模比较，保留了全部试验和失败适配记录。

在完全相同原 AOE HL 和建模协议下，尺度种子使 P9 从 33.53% 降到 28.24%，NSF-D2 从 39.33% 降到 33.53%；Sec17 从 60.61% 变为 60.84%，没有改善。这一共同协议对照比按案例拼接最好值更能说明新增方法的实际效果。

P9 的另一路校准 HL＋原 P+ 种子得到 26.00%，校准 HL 对照为 26.62%；同一校准 HL 配合尺度种子为 27.01%，因此该新种子并非普遍优于原种子。

NSF-D2 在原 HL 完全固定时，原 P+ 为 39.33%，尺度种子为 33.53%；原 HL 无种子对照为 33.59%。改进主要是消除了原 P+ 种子的损失，还不能说明胜过无种子基线。单独的校准 HL 路径表现更差，不应推广为默认算法。

Sec17 的所有协议仍失败，未得到可报告为成功解出的模型。更强的工作数据 LLG、较高 DM FOM、序列赋值比例或较好地图偏度都不足以证明位点和相位正确。传统 DM、单纯退出 HL 约束也没有解决此例。

## 全部第五轮终点

下面每个格子都是明确定义的相位/似然协议；不能将各案例最小值拼成一个统一算法的性能。补充组是在看到部分主实验结果后增加的探索性对照。原胰岛素 28.23% 属于历史不同协议，本轮没有复测，不纳入这些表格。

| 固定第五轮终点 | P9 | NSF-D2 | Sec17 |
|---|---:|---:|---:|
| 原 HL，无种子 | 33.94% | 33.59% | 59.68% |
| 校准 HL，无种子 | 26.62% | 38.07% | 56.48% |
| 校准 HL＋原 P+ | 26.00% | 38.87% | 59.82% |
| 校准 HL＋尺度 P+ | 27.01% | 37.47% | 60.43% |
| 独立 Phaser HL 参考 | 26.36% | 34.58% | 57.90% |
| 原 HL＋原 P+（探索） | 33.53% | 39.33% | 60.61% |
| 原 HL＋尺度 P+（探索） | 28.24% | 33.53% | 60.84% |


## Sec17 传统 DM 探索对照

| 协议 | Rwork | Rfree |
|---|---:|---:|
| classic_calibrated_hl_control_runtime_v2 | 51.70% | 59.22% |
| classic_calibrated_hl_pplus_runtime_v2 | 50.84% | 56.19% |
| classic_calibrated_hl_observed_scale_pplus_runtime_v2 | 46.05% | 59.15% |
| classic_phaser_reference_runtime_v2 | 45.35% | 57.87% |


## 模型完整性与立体化学

返回码为零表示验证工具正常运行，不代表模型没有问题。CA 数目和序列覆盖率均不是正确残基的真值计数。

| 案例/协议 | CA/已命名CA | 唯一序列覆盖 | Clashscore | Rama离群 | Cβ偏差 |
|---|---:|---:|---:|---:|---:|
| p9/legacy_hl_control | 128/128 | 91.4% | 14.18 | 0.82% | 3.0 |
| p9/calibrated_hl_control | 136/136 | 97.8% | 6.0 | 2.99% | 5.0 |
| p9/calibrated_hl_pplus | 134/134 | 96.4% | 2.83 | 0.0% | 5.0 |
| p9/calibrated_hl_observed_scale_pplus | 137/137 | 98.6% | 7.32 | 0.0% | 5.0 |
| p9/phaser_reference | 133/133 | 95.7% | 2.85 | 0.0% | 5.0 |
| p9/legacy_hl_pplus | 134/126 | 89.2% | 15.98 | 2.46% | 9.0 |
| p9/legacy_hl_observed_scale_pplus | 131/131 | 94.2% | 12.04 | 0.78% | 3.0 |
| nsf_d2/legacy_hl_control | 274/249 | 90.1% | 27.99 | 2.93% | 17.0 |
| nsf_d2/calibrated_hl_control | 309/258 | 88.3% | 39.55 | 3.31% | 18.0 |
| nsf_d2/calibrated_hl_pplus | 328/248 | 83.5% | 36.39 | 2.98% | 11.0 |
| nsf_d2/calibrated_hl_observed_scale_pplus | 311/261 | 89.4% | 37.52 | 5.26% | 12.0 |
| nsf_d2/phaser_reference | 283/233 | 85.3% | 25.46 | 1.76% | 12.0 |
| nsf_d2/legacy_hl_pplus | 285/230 | 82.8% | 41.15 | 5.56% | 12.0 |
| nsf_d2/legacy_hl_observed_scale_pplus | 304/252 | 86.4% | 26.55 | 4.13% | 13.0 |
| sec17/legacy_hl_control | 316/138 | 30.5% | 50.67 | 20.34% | 9.0 |
| sec17/calibrated_hl_control | 219/23 | 7.9% | 37.87 | 31.58% | 2.0 |
| sec17/calibrated_hl_pplus | 361/136 | 28.4% | 57.18 | 17.59% | 4.0 |
| sec17/calibrated_hl_observed_scale_pplus | 348/120 | 32.5% | 61.47 | 21.43% | 9.0 |
| sec17/phaser_reference | 415/115 | 32.5% | 63.55 | 12.63% | 5.0 |
| sec17/legacy_hl_pplus | 307/159 | 42.5% | 64.81 | 17.99% | 11.0 |
| sec17/legacy_hl_observed_scale_pplus | 351/152 | 28.1% | 57.35 | 21.77% | 11.0 |
| sec17/classic_calibrated_hl_control_runtime_v2 | 198/80 | 22.9% | 48.66 | 14.29% | 4.0 |
| sec17/classic_calibrated_hl_pplus_runtime_v2 | 271/135 | 27.7% | 60.3 | 22.22% | 1.0 |
| sec17/classic_calibrated_hl_observed_scale_pplus_runtime_v2 | 449/153 | 27.7% | 64.54 | 13.39% | 7.0 |
| sec17/classic_phaser_reference_runtime_v2 | 453/190 | 34.9% | 56.84 | 20.24% | 2.0 |


## 单独的相位约束退出诊断

固定 20 循环终点，不用它替代原第五轮，也不按两者 Rfree 较低者选择模型。

| 案例/协议 | 门控 | Rfree |
|---|---|---:|
| p9/legacy_hl_control | completed | 34.04% |
| p9/calibrated_hl_control | completed | 26.74% |
| p9/calibrated_hl_pplus | completed | 26.34% |
| p9/calibrated_hl_observed_scale_pplus | completed | 26.94% |
| p9/phaser_reference | completed | 26.53% |
| p9/legacy_hl_pplus | completed | 34.11% |
| p9/legacy_hl_observed_scale_pplus | completed | 27.44% |
| nsf_d2/legacy_hl_control | completed | 32.31% |
| nsf_d2/calibrated_hl_control | completed | 38.06% |
| nsf_d2/calibrated_hl_pplus | completed | 38.98% |
| nsf_d2/calibrated_hl_observed_scale_pplus | completed | 37.80% |
| nsf_d2/phaser_reference | completed | 33.84% |
| nsf_d2/legacy_hl_pplus | completed | 39.44% |
| nsf_d2/legacy_hl_observed_scale_pplus | completed | 34.39% |
| sec17/legacy_hl_control | completed | 60.42% |
| sec17/calibrated_hl_control | skipped_gate | — |
| sec17/calibrated_hl_pplus | skipped_gate | — |
| sec17/calibrated_hl_observed_scale_pplus | completed | 60.94% |
| sec17/phaser_reference | completed | 58.69% |
| sec17/legacy_hl_pplus | completed | 60.21% |
| sec17/legacy_hl_observed_scale_pplus | skipped_gate | — |
| sec17/classic_calibrated_hl_control_runtime_v2 | skipped_gate | — |
| sec17/classic_calibrated_hl_pplus_runtime_v2 | skipped_gate | — |
| sec17/classic_calibrated_hl_observed_scale_pplus_runtime_v2 | skipped_gate | — |
| sec17/classic_phaser_reference_runtime_v2 | completed | 60.32% |

## 已落实的开发判断与后续顺序

1. 保留 PHCOMB/FOM 交接修正与审计记录。尺度种子和新的 HL 模型继续作为独立研究选项；本轮没有将未经普遍验证的组合设成默认 P++。
2. Sec17 优先核验重原子子结构及误差模型。用工作集内部再留出的异常验证反射比较位点增删和尺度/B/占有率耦合，避免只用训练 LLG 越高越好来扩充位点。已有候选完整保留，下一阶段无需重新寻找丢失代码。
3. 在冻结的观测似然下研究 P+ 置信度和双峰噪声积分的一致校准，并加入可信部分模型反馈。当前新增尺度模块刻意保留旧误差权重和 Sim，仅仅改变几何，尚未实现这些完整方案。
4. NSF-D2 的后续重点为可验证的侧链/构象修复与工作集确定的精修协议；现有模型尚有明显 Cβ 和碰撞问题。单纯增加精修次数或取消相位约束并没有普遍改善。
5. 选定共同算法后，用未参与调试的新公开案例重新进行冻结协议验证；本轮示例属于开发数据，不能作为独立泛化证明。

这些是后续开发任务，不是已经获得的效果。没有采用更换自由集、剔除困难案例、最低轮次选择、参考目标结构拟合或单报最好结果来达到 30%。


---

# P++ 优化：完整分辨率、相位交接与异常尺度校准

实验根目录：`/Volumes/Fast_4T/GPT_WORK/OASIS4.3-pplus-target30-20261024-111000`。全部结构化数据见该目录的 `RESULTS.json` 与 `RESULTS.csv`。

本轮目标是在固定自由集下改善公开 SAD 测试案例的建模与 Rfree。研究比较包含所有预先声明的组和后续明确标记的探索性对照；每组报告第五轮建模后的最终 REFMAC 数值，不挑最低轮次。是否达到约 30%，由结果判断，不能预先保证。

## 本轮具体修正

1. **Buccaneer 相位交接**：REFMAC 的 `PHCOMB/FOM` 是一组质心相位与置信度，`FWT/PHWT` 是地图系数。此前反复把 `PHWT` 与 `FOM` 配对。现在有唯一、正确类型的 `PHCOMB` 时使用它，同时保留独立地图通道；缺少该列时显式记录原有回退契约。每轮保存输入指纹和采用的列名。七个新增测试及完整测试套件验证该修正。
2. **P+ 异常尺度种子消融**：新增 `pplus_observed_scale.py`，只用工作集非中心反射拟合异常信号的有效尺度，改变双峰几何；原有误差权重、Sim、关系图及原 AOE 模型的 HL 保持原样。这不是联合重校准的完整似然，也不是实测 f″。五个新增测试涵盖尺度恢复、无信号、冻结数组与中心单值反射等契约。
3. **外部 Parrot 真正双输入**：固定 `colin-hl`，另传 `colin-fc`。相同输入下改变种子的黑盒对照证明该通道影响输出；同一输入的重复计算及去除 HL 的对照保留在 `parrot_dual_input_probe/`。保留真实输出 F_phi，初轮也分别传地图相位与概率质心相位。该研究适配器位于实验目录，尚未接入默认产品。

## 数据与困难案例

- P9 使用完整观测分辨率约 1.744 Å，NSF-D2 约 2.400 Å，Sec17 约 3.303 Å。此前前两例按 2.5 Å 筛选。恢复完整数据会重新拟合 French–Wilson，重叠反射的振幅也可能改变，不能将前后差异只归因于 P+。
- 保留全部原有 HKL 自由标记。新出现的 HKL 使用相同 SHA256 命名空间规则扩展；P9 21557 行/1071 自由行，NSF-D2 13235/642，Sec17 8288/421。原始 Scalepack 没有提供已发表的原始自由集，这些是开发过程中建立的固定测试集合。
- Sec17 工作集 xtriage 没有检出冰环、孪晶或平移 NCS；异常可测信号主要集中在低分辨率。依据其工作数据信号截止，单独进行 5.2 Å 子结构搜索，保留全部候选。最终共识经全数据工作集 Phaser 精修完成的 LLG 为 200.219，旧搜索约 23.426；最高原生三位点候选完成后为 178.783。LLG 改善不证明位点或蛋白地图正确。
- 两手均经相同 Phaser→Parrot 工作数据协议。选手依据真实地图标准化偏度，未强制采用文件头中的对映空间群。该统计仍不是已知真值。
- P9 有一个各向异性重原子，而原生 P+ 的原子模型仅支持各向同性 ADP。原始 Phaser 输出完整保存；P+ 使用单独命名的迹等效 B 因子副本，所有相关组共用这一近似。Phaser 独立对照仍保留其各向异性相位计算。

## 对照的含义

主实验每例五组：原 AOE 观测 HL 对照；校准 HL 对照；校准 HL＋原 P+ 种子；校准 HL＋观测尺度 P+ 种子；独立 Phaser 参考。三个“校准 HL”组逐字节共享相同 HL。这里的校准 HL 使用非中心工作反射的异常方差尺度拟合，以及 `hypot(原 sigmaD,0.1*|alpha Q|)`；它仍是小异常信号近似，不等同于 Phaser 的联合 Bijvoet 似然。

在首批 NSF 结果显示校准 HL 没有普遍改善后，另加了每例两组探索性对照：原 AOE HL＋原种子、原 AOE HL＋尺度种子。它们与原 HL 对照同样逐字节共享 HL，避免把似然改变与种子几何改变混淆。

Sec17 另有四组传统 CCP4 DM 对照，统一 10 个 DM 循环，再采用相同五轮建模。这是独立引擎试验，不能与五循环 Parrot 作单因素归因。首轮由于旧式全局 harvesting 路径已存在而失败，原日志保留；后续显式 `NOHARVEST` 正常结束。DM 在四组输出中均遗漏了两个位于取整高分辨率边界的地图行；标准化适配时只令其地图/FOM/HL 权重为零，全部原始 native F、sigma 与自由标记保持不变。引擎原输出也完整保留。

所有成功建模组另做相同的探索性“退出原始相位约束”诊断：只有唯一序列覆盖率至少 30% 时，才运行固定 20 循环的仅振幅 REFMAC 精修。覆盖率使用已命名残基的局部序列比对，并对重复位置去重。它不是原始论文规则的完整实现，也不能证明残基正确；诊断结果单独列出，不替换主终点。该 30% 是模型序列覆盖率，与目标 Rfree 30% 含义不同。[Wu 等，2009](https://cryst.iphy.ac.cn/Publication/PDF_files/2009/SAD-phasing_in_OASIS4.pdf)描述了模型逐渐完整后停止输入双峰 HL 的原则。

## 结果读取原则与局限

Buccaneer 负责建模，表中 Rfree 由 REFMAC 精修输出。每组有 5 轮建模，各轮 10 次精修；首轮 Buccaneer 3 次循环，其后各 2 次。所有成功组做 MolProbity 验证，低 Rfree 不代表缺失残基、侧链或立体化学问题已经解决。

工作集用于子结构、初始相位及 DM。French–Wilson 等强度预处理使用全部测量值；后续外部引擎产生的地图没有被认证为完全排除自由集信息。这些公开示例又已经用于开发比较，因此本报告数值不是端到端独立盲测的认证结果。未使用完整目标蛋白坐标，也未以参考蛋白进行分子置换。

## 可复现与恢复

`experiment_protocol.json`、各补充协议、输入/输出 SHA256、完整外部命令、失败尝试、逐轮报告与验证输出均位于本实验目录。开始前快照为 `source_snapshot/`，新实现与测试的单独证据在 `build_handoff_fix/` 和 `observed_scale_fix/`。原始发布 ZIP、Git 历史及已有未提交修改得到保留；上游 31 个受校验文件没有变化。

完整测试结果为 **429 passed，1 个 reference_data 测试未运行**；`git diff --check` 通过。数值源文件在本轮全部建模期间保持冻结。

相关原始文档：[Phenix SAD 示例](https://phenix-online.org/documentation/tutorials/sad_phasing.html)、[CCP4 Parrot](https://www.ccp4.ac.uk/html/parrot.html)、[CCP4 DM](https://www.ccp4.ac.uk/html/dm.html)。REFMAC 的 PHCOMB 契约另有已安装官方文档证据，见 `build_handoff_fix/PHASE_CHANNEL_EVIDENCE.json`。
