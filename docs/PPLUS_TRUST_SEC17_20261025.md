# 观测似然兼容、限制过度自信的 P+ 分支

新增可选模块 `phasebuild.pplus_trust`。现有GUI和`--method pplus`的旧默认路径没有被自动替换；该新分支需要显式调用。

它从同一套原始SAD观测HL计算两个半圆分支的质量、条件首矩和二阶矩，再用对称性三相位关系场更新P+质量。逐反射限制更新分布相对原观测质量的Bernoulli KL距离，并使用阻尼。条件形状保持固定，原观测似然不累积；后验只作为地图种子，不冒充新的实验HL。

```python
from phasebuild.pplus_trust import phase_trusted_pplus

seed = phase_trusted_pplus(
    data, original_observation_hl, heavy_atom_phase_reference, relations,
    samples=1440, kl_budget=0.0005, cycles=100, damping=0.25,
)
```

输入HL必须与下游DM/精修使用的原始观测先验相同、HKL顺序一致。`heavy_atom_phase_reference`是异常重原子结构因子的相位，不是蛋白参考结构相位。关系图必须从同一work数据及组成生成。正KL预算是强制的；调用者还应检查实际概率/相位变化，不能把无更新结果计为P+成功。

SEC17验证路线使用**Phaser提供初始SAD观测似然，然后真正执行P+更新**。没有额外叠加Sim重原子场，避免再次加入已在该似然中的亚结构信息。这是组合前端下的新P+优化分支，不是旧AOE似然P+已被证明独立达标。

用输出`phase/fom/probability`保存P+种子；不要构造后验HL并送回作原始实验先验。RESOLVE的起始地图应是只含HKL/FWT/PHWT的独立文件，FWT=原始work Fo×种子FOM一次，PHWT=种子相位。保留概率文件中的原始HLA–D作为DM/精修约束。不要把FOM、HL混入起始地图文件，因为RESOLVE可按名字自动读取它们。

SEC17开发结果：KL预算0.0005时，4939个非中心反射的P+概率变化超过0.01，初相位变化RMS10.90°，Fo加权地图系数相对RMS改变4.84%。从无蛋白模型开始的第五轮Buccaneer/REFMAC得到27.40% Rfree；按事先固定的单拷贝清理与保守几何精修后27.42%、292残基、Ramachandran异常2.41%、clashscore12.09。完整重跑重新计算种子、DM和模型，结果及所检查的结构因子/地图列一致。

KL预算0.0025的另一P+分支最终29.33%；更强或不限制的更新在该案例反而退化。0.0005是这个开发案例验证过的设置，并非通用校准常数。算法是固定条件形状的两状态近似，不是精确连续后验，也不保证任意关系图上的全局收敛。

原始观测/Freeset未修改，没有蛋白参考坐标、旧蛋白模型或最终AutoSol相位作为求解输入。SEC17被反复用于开发，Rfree不是新的独立盲测，仍需新案例验证；2.41%主链异常也意味着还需局部重建。本次没有证明新P+的最终Rfree全面优于Phaser对照（对照同样清理后27.28%），P+候选的几何较好。

完整记录在独立工作目录 `OASIS4.3-sec17-pplus-trust-20261025-180000`，包括全部6组终点、一个完整P+复跑和输入/自由集审计。
