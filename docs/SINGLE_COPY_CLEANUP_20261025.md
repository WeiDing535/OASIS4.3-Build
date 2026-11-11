# 单拷贝完整序列链的保守清理

可选入口为 `run_single_copy_cleanup.py`。它不改变求相位或默认建模流程，不运行精修，不依据 Rfree 删除残基。

仅适用于已知非对称单元中只有目标蛋白的一份、且模型中存在唯一一条完整精确匹配输入序列的链的情况。异源复合物、多拷贝或序列不完整时不要使用。工具要求显式给出 `--asu-copies 1`；多个完整匹配或没有完整匹配都会拒绝执行。

```bash
.conda-unified/bin/python -B run_single_copy_cleanup.py \
  --model existing_model.cif --sequence sequence.fasta \
  --asu-copies 1 --output new_cleanup_directory
```

输出为新目录中的 `model.cif` 与 `cleanup_report.json`。已有输出目录会被拒绝，不覆盖输入。保留完整匹配链及全部非蛋白残基；只移除其他链中的氨基酸残基，包括 UNK。保留的原子属性、残基编号、MSE 化学名称及晶体信息通过克隆保存。报告列出每个移除残基，并记录输入、序列和输出文件哈希。

该规则适用于额外追踪出的未知/重复片段，而非普遍的结构正确性判断。完整序列匹配不证明坐标正确；清理后必须单独精修并检查几何、密度与原始自由集守恒。该工具不会把清理后的模型用于重新求初相位。

SEC17 开发案例中，原模型的一条链完整匹配292残基，其余11条链多出87残基。清理得到292残基；与独立实验脚本产生的清理模型逐字节一致。较均衡的后续手动几何权重控制得到27.28% Rfree，仍有5.17% Ramachandran异常，不能称为出版级结构。此案例反复用于开发，Rfree不是独立盲测成绩；原求相位基准采用Phaser观测概率和自身种子，不能将该结果宣称为P+核心算法改善。

完整研究记录位于工作区独立目录 `OASIS4.3-sec17-model-cleanup-20261025-174000`。使用的REFMAC固定权重含义见[官方说明](https://www.ccp4.ac.uk/html/refmac5/keywords/xray-principal.html)。
