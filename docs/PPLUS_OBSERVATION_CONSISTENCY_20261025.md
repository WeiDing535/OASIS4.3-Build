# P+ 观测概率一致性与分支积分修复（2026-10-25）

本轮交付了可选的观测一致性 P+ 研究实现、概率相位交接修复，以及独立发现的半圆积分偏差修正。真实数据共25次分组尝试：16个固定端点完成、7个研究接口/启动失败保留、2个按用户指令暂缓。后半程仅继续SEC17，追加闭合误差与分支形状三组；同一18组合成条件也复测。全部结果与失败尝试保留。**默认 GUI/数值 P+ 未切换；这不是所有案例 Rfree≤30%的完成声明。**

Rfree 来自 REFMAC 固定第五轮，Buccaneer负责建模。新 P+ 没有使用目标、同源蛋白坐标或外部实验相位；继承的 HA/手性已在以前的工作数据研究中选定。真实数据是反复研究过的诊断和回归案例，不能当成新的盲测。

## 同预算真实数据结果

A=旧AOE实验HL+当前旧P+种子；B=新纯观测HL+其自身积分种子；C=同一新HL+旧P+种子；D1=同一新HL+新条件P+；D2=同一新HL+只修正半圆积分的新条件P+。B/C/D1/D2的实验HL文件逐字相同。

| 案例 | A 旧实现 | B 纯观测 | C 新HL+旧种子 | D1 新条件P+ | D2 积分修正版 |
|---|---:|---:|---:|---:|---:|
| p9 | 31.22% | 27.13% | 27.30% | 26.68% | 用户暂缓 |
| nsf_d2 | 35.24% | 35.82% | 38.08% | 35.25% | 用户暂缓 |
| sec17 | 57.40% | 56.77% | 56.34% | 58.81% | 58.41% |

每组均为 Phenix classic40 DM（10初始+20收缩+10最终，scale_flip=False，无NCS）与5轮Buccaneer/每轮10周期REFMAC；第一轮3周期Buccaneer、以后2周期。使用原native全分辨率观测及固定FreeR_flag=0，不按Rfree选手性、轮次、参数或提前停止。P9/NSF-D2/SEC17分辨率分别约1.744/2.400/3.303 Å。

初始DM只读工作反射，Fo×种子FOM只加权一次；建模使用实际DM FWT，逐工作HKL乘Fnative/FDM还原至原native单位，PHWT仍为实际地图相位，PHCOMB/FOM为概率中心。初始构建Fo/系数的自由集行为NaN。每一轮精修输入F/SIGF/自由标志均逐HKL与原native逐值核验。后续外部模型地图尚未认证端到端只用工作观测，因此这里是诊断Rfree。

不能把这张表与此前5周期DM或不同幅度交接的结果作单因素比较；当前A是本轮共同接口和预算下的旧算法对照。此前SEC17独立AutoSol的28.23%/29.61%属于外部流程，未作为本轮P+成绩。

P9原四组在第一轮建模前发生研究适配器错误：native HKL已在预期ASU，但DM使用不同的对称等价HKL；原始逐HKL交集仅4行，其中2行可作工作图。修正为对内存DM调用Gemmi ensure_asu，连同P型相位和A型HL按对称转换，原native不动；交集恢复21557行、有效构建工作图20241行。全部四组恢复运行使用逐字相同的原40步DM产品，并要求准备DM输入SHA完全相同，未重新定相或改DM参数。对比Fourier地图及标量列确认仅对称表示与浮点绕相舍入改变。该问题属于研究适配器交接失败，不能算P+定相失败。[审计与官方语义链接](/Volumes/Fast_4T/GPT_WORK/OASIS4.3-pplus-calibrated-20261024-221600/ASU_HANDOFF_AUDIT.md)。

## 几何、过建和序列一致性

| 案例/分支 | Rwork/Rfree | CA/命名CA | 去重序列匹配 | 重复赋值 | Rama离群% | Clashscore |
|---|---:|---:|---:|---:|---:|---:|
| p9/legacy | 28.59%/31.22% | 127/126 | 126/139 | 0 | 0.82 | 16.06 |
| p9/observation | 24.42%/27.13% | 137/137 | 137/139 | 0 | 0.74 | 6.48 |
| p9/new_hl_old_seed | 25.80%/27.30% | 129/129 | 129/139 | 0 | 0.8 | 7.88 |
| p9/conditional_pplus | 23.87%/26.68% | 132/132 | 132/139 | 0 | 0.77 | 7.19 |
| p9/conditional_pplus_v2 | 用户暂缓 | — | — | — | — | — |
| nsf_d2/legacy | 26.26%/35.24% | 247/242 | 242/273 | 0 | 2.54 | 27.78 |
| nsf_d2/observation | 28.06%/35.82% | 275/223 | 222/273 | 1 | 1.38 | 21.79 |
| nsf_d2/new_hl_old_seed | 28.08%/38.08% | 330/228 | 228/273 | 0 | 3.12 | 34.13 |
| nsf_d2/conditional_pplus | 26.15%/35.25% | 318/247 | 240/273 | 0 | 2.9 | 26.9 |
| nsf_d2/conditional_pplus_v2 | 用户暂缓 | — | — | — | — | — |
| sec17/legacy | 50.57%/57.40% | 211/55 | 49/292 | 0 | 12.77 | 55.75 |
| sec17/observation | 52.00%/56.77% | 159/55 | 48/292 | 0 | 22.45 | 56.92 |
| sec17/new_hl_old_seed | 48.78%/56.34% | 241/83 | 69/292 | 14 | 14.49 | 44.51 |
| sec17/conditional_pplus | 47.71%/58.81% | 296/58 | 41/292 | 3 | 29.17 | 51.93 |
| sec17/conditional_pplus_v2 | 48.68%/58.41% | 288/55 | 48/292 | 7 | 22.45 | 45.53 |

所有完成端点均运行MolProbity。去重序列匹配采用每条生成链一次最佳局部比对（匹配+2、错配−3、gap−2；至少8个精确匹配且identity≥80%），跨链目标位置取并集。它衡量序列一致性，不能验证三维坐标；一条链上的多个分离岛可能被此启发式低估。CA或命名CA超过目标长度不代表完整率。低Rfree必须连同过建、重复赋值和几何离群解释。

## SEC17专项追加：闭合误差及分支形状

按用户要求P9/NSF-D2积分修正版复测已中断并留存中间产物，不当作算法失败。此后只执行SEC17。新增全局非负闭合误差σeff²=σD²+τ²，报告σD未减少、固定工作反射无新增剔除；α/τ的内部留出平均NLL从4.383405降到2.095028（旧Gaussian方差2.268943），但预测评分改善没有转化成可用结构。

| 同一闭合误差HL/固定预算 | Rwork/Rfree | CA/命名CA | 去重序列匹配 | 重复赋值 | Rama离群% | Clashscore |
|---|---:|---:|---:|---:|---:|---:|
| closure_observation | 44.37%/54.31% | 219/77 | 56/292 | 21 | 24.64 | 83.29 |
| closure_fixed_pplus | 49.11%/58.68% | 160/71 | 65/292 | 0 | 24.59 | 69.87 |
| closure_continuous_shapes | 46.41%/58.82% | 255/139 | 106/292 | 9 | 22.22 | 66.5 |

全工作集α=0.024385、τ=1.733736（F单位），纯观测FOM仅0.04156。三组纯HL字节相同，分别比较自身积分种子、固定分支P+、连续分支形状控制；Sim/关系图完全复用，权重/阻尼/轮数不变。连续控制不冒充原双点P+或新的已发表方法。

在同HL与预算下，加当前场的两种种子都没有改善SEC17；仍需重建散射/HA不确定性与场校准，而非以更高FOM为目标。研究模型10项数学/接口测试通过，6组合成观测α/τ回收均在预注册界内；此为观测模型测试，不是6个蛋白结构。

首次专项启动在冻结检查前停止、三组模式名断言在DM前停止，原记录均保留；独立execution_profile_fix仅修正允许的预注册模式，经三个实际MTZ输入契约后重跑，数值实现及已拟合参数未改变。研究代码仅放在研究目录，90份产品源码仍保持积分修正版冻结、默认算法未切换。

[SEC17专项完整报告](/Volumes/Fast_4T/GPT_WORK/OASIS4.3-Build/docs/SEC17_PPLUS_CLOSURE_20261025.md)；[专项协议/冻结](/Volumes/Fast_4T/GPT_WORK/OASIS4.3-pplus-calibrated-20261024-221600/sec17_closure/RESEARCH_FROZEN.json)；[三个固定端点](/Volumes/Fast_4T/GPT_WORK/OASIS4.3-pplus-calibrated-20261024-221600/sec17_closure/execution_profile_fix/REAL_ENDPOINT_AUDIT.json)。

## 实现与接口

- [pplus_scale.py](/Volumes/Fast_4T/GPT_WORK/OASIS4.3-Build/src/phasebuild/pplus_scale.py)：固定HA、均匀相位下规范化的Gaussian-DANO边际似然，拟合有效异常尺度α；显式保留α=0、分盆搜索及训练网格加倍精度检查。预测评分不重新拟合α，积分不达标显式失败。α不等于实测f″，不能给HA/手性排名。
- [pplus_conditional.py](/Volumes/Fast_4T/GPT_WORK/OASIS4.3-Build/src/phasebuild/pplus_conditional.py)：积分实验HL在两个半圆内的质量及一/二阶圆矩；固定条件形状，只更新P+分支质量。logit p=log(m+/m−)+Re[T·conj(z+−z−)]；观测只用一次，第一轮仅用观测，精确双点极限回到原tanh形式。它是受限的两状态变分近似，并非完整连续后验。
- [pplus_calibrated.py](/Volumes/Fast_4T/GPT_WORK/OASIS4.3-Build/src/phasebuild/pplus_calibrated.py)：冻结观测、归一化、Sim/关系图及身份指纹；观测先验和更新种子分开导出，种子不携带HL。模块名不表示已完整校准。库API为prepare_observation_pplus、observation_pplus_prior、phase_observation_pplus，尚未替换默认GUI算法。
- [hybrid_recycle.py](/Volumes/Fast_4T/GPT_WORK/OASIS4.3-Build/src/phasebuild/hybrid_recycle.py)：概率相位按PHCOMB→PHIC→明确的旧PHWT回退，验证列唯一性/类型。实际FWT/PHWT地图保持独立。
- [sad_bridge.py](/Volumes/Fast_4T/GPT_WORK/OASIS4.3-Build/src/phasebuild/sad_bridge.py) 与 [sad_pipeline.py](/Volumes/Fast_4T/GPT_WORK/OASIS4.3-Build/src/phasebuild/sad_pipeline.py)：为DM概率中心显式写PHCOMB，防止旧native PHIC抢占种子，并提前拒绝输入中保留列冲突；Fo×FOM一次的产品适配合同未变。本轮研究实际DM幅度交接另有显式单位审计。

公开OASIS双解概率更新是这项研究扩展的背景；这里没有声称逐项复现发表版或Phaser。联合Bijvoet及散射误差处理仍属于后续工作。[OASIS原始方法](https://cryst.iphy.ac.cn/Publication/PDF_files/2006/Dual-space_fragment_extension.pdf)；[Phaser SAD似然原始论文](https://pmc.ncbi.nlm.nih.gov/articles/PMC3069749/)。

## 独立发现并修正的数值偏差

D1在绝对角网格上硬切半圆，即使完整圆周矩已收敛，分支质量仍有一阶切口误差。小信号观测似然关于HA参考轴镜像对称，连续积分应有m+=m−=0.5；D1在SEC17读回HL上最大偏离约0.01745。

D2在参考相对坐标用成对半圆中点，仍保持总1440节点（每侧720）及原centric精确两态路径。它消除了硬切口，一般平滑半圆积分为二阶收敛；不能保证所有任意尖锐HL在1440节点均已达到微小误差。高精度独立adaptive-quad测试使用16384节点，另用低节点验证约四倍误差下降。

| 案例 | D1重放与保存种子最大一阶矩差 | 最大概率差 | D1半圆质量偏离0.5 | D2半圆质量偏离0.5 |
|---|---:|---:|---:|---:|
| p9 | 2.55e-07 | 3.08e-08 | 0.00501109 | 3.55e-15 |
| nsf_d2 | 2.49e-07 | 2.97e-08 | 0.00483501 | 1.42e-14 |
| sec17 | 2.59e-07 | 2.99e-08 | 0.0174486 | 2.84e-14 |

D2没有重新拟合α；重建双精度HL的float32值必须与原纯HL逐项完全相同，旧D1必须在预注册2e−6最大一阶矩/概率误差内重放，否则停止并保留失败。DM/精修仍复制同一纯HL字节。新旧版分别冻结，原12次尝试及4次P9接口恢复完成、质量审计后才改live源码，未覆盖旧结果。修复根据独立对称性与积分检验制定，非按Rfree追调。

## 合成数据与置信度

18组为6个随机原子云（P1/P21/P212121各2个种子）×3个噪声等级，不是18种独立蛋白。完整结构只生成/评分，方法仅接收观测、3Se、组成。D2复测同一批数据用于数值回归，不能作为新泛化证据。

| 版本/方法 | 条件平均一阶圆矩MSE↓ | 条件平均Brier↓ | 反射加权Brier↓ | 高FOM组预测/实际cos |
|---|---:|---:|---:|---:|
| D1/legacy_Pplus | 0.621777 | 0.178261 | 0.178627 | 0.8459/0.8900 |
| D1/new_observation | 0.762584 | 0.250005 | 0.250010 | 0.8966/0.8755 |
| D1/new_conditional_Pplus | 0.578036 | 0.177671 | 0.182016 | 0.9425/0.8039 |
| D2/legacy_Pplus | 0.621777 | 0.178261 | 0.178627 | 0.8459/0.8900 |
| D2/new_observation | 0.762584 | 0.250000 | 0.250000 | 0.8966/0.8755 |
| D2/new_conditional_Pplus | 0.578038 | 0.177667 | 0.182012 | 0.9425/0.8040 |

D1条件P+的一阶矩MSE优于旧版14/18组。


D2条件P+的一阶矩MSE优于旧版14/18组。


D1的4组一阶矩退化全在P212121；反射加权Brier略差、FOM>0.8组预测0.9425而实际相位cos均值0.8039，说明置信度明显偏高。D2应按表中全部指标解释，不能用一个较好的平均MSE或更高FOM宣布完全校准。未根据这些分数修改triplet weight、damping或轮数。

## SEC17仍需解决的误差模型与迭代问题

| 案例 | 内部留出相位边际NLL↓ | Gaussian方差模型NLL↓ | 说明 |
|---|---:|---:|---|
| p9 | 3.592847 | 3.617634 | 固定80/20工作反射切分，HA已见全部工作反射 |
| nsf_d2 | 2.611808 | 2.618160 | 固定80/20工作反射切分，HA已见全部工作反射 |
| sec17 | 4.383405 | 2.268943 | 固定80/20工作反射切分，HA已见全部工作反射 |

SEC17的1438个内部留出工作反射，新模型总NLL劣势约3040.60；最差20行贡献88.65%，最低分辨率四分位19.617–5.308 Å贡献93.53%。例如(1,5,1)的DANO=51.2896、σD=0.73695，而2α|Q|=2.22485，超出模型均值跨度约66.58σ。这里只定位不匹配，**没有认定这些观测错误或删除/截断它们**。

现模型只用DANO=2α|Q|cos(φ−reference)+Gaussian(报告σD)，未传播缺失位点、坐标/ADP/占有率、散射与观测协方差误差，也不等同联合Bijvoet幅度模型。报告σ与模型闭合误差必须分开；均匀相位的尺度目标不能独立验证HA几何。

D1第五次更新最大分支概率变化P9/NSF/SEC17约0.846/1.000/0.453，固定五次并非收敛证明。Sim仍沿用旧Wilson近似，有向截断三相关系存在相关性，局部变分logit不保证全局目标单调或概率校准。半圆积分修复没有自动解决这两类问题。

下一阶段优先：核查上述低分辨率尾部的原始观测与误差来源；在工作训练数据上引入规范化、可辨识的模型闭合误差并用固定预测似然检验；保留F+/F−各自误差/相关信息以支持联合Bijvoet升级；再审查关系重复度、方向/度数依赖和更新稳定性，以及是否应允许分支内形状随场更新、保持观测只用一次。新增带HA扰动/缺位点及相关观测误差的独立合成结构与未使用真实案例。参数选择使用数学或预测标准，最终未触碰测试集只用于一次终验。

## 测试与可复现证据

- D1完整测试：482 passed, 1 deselected in46.79s；D2完整测试：490 passed, 1 deselected in 48.61s。reference_data外部参考测试均明确排除，非默认核心测试失败。
- 初次D1完整检查的两个独立重算ULP精确比较失败已保留；只将测试改为有界浮点比较，数值源码未因此改动。同一冻结上下文字节级复用合同仍保留。
- 积分修复暂存阶段最初过紧的有限节点阈值/参考quad警告也保留。修正参考积分尺度并单列高精度与收敛测试，未为测试失败更改候选算法或真实案例预算。
- [全部25次尝试与逐轮/几何完整数据](/Volumes/Fast_4T/GPT_WORK/OASIS4.3-pplus-calibrated-20261024-221600/ALL_ENDPOINT_RESULTS.json)；[P9接口恢复4组](/Volumes/Fast_4T/GPT_WORK/OASIS4.3-pplus-calibrated-20261024-221600/asu_handoff_recovery/REAL_ENDPOINT_AUDIT.json)；[D2追加3端点](/Volumes/Fast_4T/GPT_WORK/OASIS4.3-pplus-calibrated-20261024-221600/branch_quadrature_v2/REAL_ENDPOINT_AUDIT.json)；[输入独立审计](/Volumes/Fast_4T/GPT_WORK/OASIS4.3-pplus-calibrated-20261024-221600/REAL_INPUT_AUDIT.md)；[数值独立审计](/Volumes/Fast_4T/GPT_WORK/OASIS4.3-pplus-calibrated-20261024-221600/NUMERICAL_REVIEW.md)；[SEC17尾部诊断](/Volumes/Fast_4T/GPT_WORK/OASIS4.3-pplus-calibrated-20261024-221600/PREDICTIVE_TAIL_DIAGNOSTIC.md)。
- [原18组合成审计](/Volumes/Fast_4T/GPT_WORK/OASIS4.3-pplus-calibrated-20261024-221600/synthetic/AUDIT_REPORT.md)；[D2合成结果](/Volumes/Fast_4T/GPT_WORK/OASIS4.3-pplus-calibrated-20261024-221600/branch_quadrature_v2/synthetic/SUMMARY.json)。
- [原协议](/Volumes/Fast_4T/GPT_WORK/OASIS4.3-pplus-calibrated-20261024-221600/EXPERIMENT_PROTOCOL.json)、[积分修复追加协议](/Volumes/Fast_4T/GPT_WORK/OASIS4.3-pplus-calibrated-20261024-221600/branch_quadrature_v2/EXPERIMENT_PROTOCOL.json)、两个SOURCE_FROZEN与source_frozen快照、全部脚本/原日志/MTZ/模型/MolProbity保留。研究脚本依赖现有本地CCP4 9与Phenix1.21.2运行时，证据ZIP不是可移植引擎安装包。
- [源码及输入保护审计](/Volumes/Fast_4T/GPT_WORK/OASIS4.3-pplus-calibrated-20261024-221600/PRESERVATION_AUDIT.json) 与最终MANIFEST记录本轮闭合指纹；数值源码修改前快照、原tracked差异、两个冻结版本、最终测试快照齐全。未安装依赖、重置Git、丢弃改动、提交或推送。
- 原源码ZIP保持SHA256 `87b08811472b1b603dc03a8b40906420b8af41385f0ffaa7802ef09e86e6ab93`。独立关闭的证据ZIP及CRC/SHA回读结果见 [归档验证记录](/Volumes/Fast_4T/GPT_WORK/OASIS4.3-pplus-calibrated-20261025-archive-verification.json)。
