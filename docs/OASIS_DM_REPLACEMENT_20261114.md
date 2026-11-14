# oasis-dm v2 替换旧密度修饰引擎:接入方式、验收数据与 `--cmbwt` 评估

日期:2026-11-14。范围:OASIS4.3 全部密度修饰路径。结论先行:

1. **替换已完成**。`--dm-engine oasis-dm` 是默认值,旧的外部密度修饰引擎
   已从程序中彻底删除(模块、命令行选项、打包标签、测试全部移除)。整条链
   不再需要任何外部密度修饰程序,也不再需要 CCP4 的 `hist.lib`。
2. **效果**:胰岛素 2.1 Å 上,DM 图从 CC 0.3646 / wMPE 55.4° 提升到
   CC 0.7856–0.8049 / wMPE 36.0–37.1°,下游建库与旧引擎统计上不可区分。
3. **`--cmbwt` 没有可靠的“算前评估”手段**。原因是结构性的:它只改变交给
   下游的 HL 权重,几乎不改变本轮地图,所以一切基于地图质量的指标(包括
   留出集 CC)都看不见它。实测留出集 CC 全距 0.005,且 CC 最高的 0.5 下游
   R-free 最差。默认值保持 0.4;`auto` 只在证据明确时才覆盖默认值。

---

## 1. 改了什么

### 1.1 新引擎与随仓库资源

| 文件 | 作用 |
| --- | --- |
| `src/phasebuild/vendor/oasis_dm/oasis_dm_core.py` | 上游 v1(未改动),提供 `Config`/`DensityModifier`/`hl_centroid` |
| `src/phasebuild/vendor/oasis_dm/oasis_dm_v2.py` | 上游 v2(未改动),`DensityModifierV2` |
| `src/phasebuild/vendor/oasis_dm/PROVENANCE.md` | 上游路径、sha256、`hist.lib` 许可说明、覆盖方式 |
| `src/phasebuild/oasis_dm.py` | **全部 OASIS 适配都在这里**:定位 Python、暂存输入、合并模型 HL、写 `03_dm` 列、cmbwt 标定 |
| `lib/data/hist.lib` | CCP4 直方图库副本(13373 字节),oasis-dm 不再依赖 CCP4 安装目录 |

两个上游文件的 sha256 与丁玮工作目录中的副本一致(`oasis_dm_core.py`
`6647fce3…fdd6c`,`oasis_dm_v2.py` `3588a1cf…83e8`),`PROVENANCE.md` 里可逐字复核。

### 1.2 默认值与开关

| 参数 | 旧 | 新 |
| --- | --- | --- |
| `dm_engine`(profile / P+) | 外部引擎 | `oasis-dm` |
| `dm_cycles` | 40(仅外部引擎) | 12(oasis-dm 自停) |
| `dm_cmbwt` | — | `0.4`(数值,或 `auto`) |
| `dm_holdout` | — | `0.05` |
| `dm_hist_lib` | 硬编码 `/Applications/ccp4-9/lib/data/hist.lib` | 仓库 `lib/data/hist.lib`,可覆盖 |

`--dm-engine` 二选一:`oasis-dm`(默认)、`internal`;其它取值一律报
「Unknown density-modification engine」。两条 P+ 流程(`pplus_bundled`、
`pplus_internal`)只接受 `oasis-dm`,给 `internal` 直接报错 —— 旧的内部 MLHL
密度修饰只作为 `likelihood`/`pplus` 的回退,不再冒充 P+ 的密度修饰。

回退规则:oasis-dm 需要 gemmi/numpy/scipy 的 Python,`oasis_dm.locate()` 依次
尝试 `PHASEBUILD_OASIS_DM_PYTHON`、当前解释器、`<runtime>/bin/python`。
`likelihood`/`pplus`/`pplus_internal` 找不到时回退到 `internal` 并在
`pipeline_report.json` 写 `dm_engine_requested` / `dm_engine_fallback`;
P+ 两条流程找不到时直接失败,不静默降级。

### 1.3 调用链

`sad_pipeline` → `run_trace_feedback.py` → `trace_feedback.run(dm_engine=…)` →
`oasis_dm.run(...)`。`oasis_dm.run` 的签名与旧的外部适配器对齐,所以
trace/Fcalc 反馈循环一行没改:每轮 DM 之后 trace 出模型,下一轮把模型
`Fc/φc` 并入 HL 再修饰。并入时有一个单位问题必须处理:实验 HL 四元组是
归一化量(`|HL|` 约 1–2),而 Sim 权重 `X = 2 Fo D Fc/(1-D²Fc²/Fo²)` 是观测振幅
单位;直接相加会让 12 个残基的碎片比 6 个硫位点还重 7 倍,地图当场塌掉(实测
留出集 CC 0.844 → 0.700 → 0.528)。`combine_model_hl` 因此按壳层把模型贡献
缩放到实验 `|HL|` 的量级,再乘以模型解释的散射强度比例
`<Fc²>/<Fo²>`(空模型 0,完整模型约 1)。修正后同一轮模型只贡献 0.30,实验
HL 是 1.64。`merge_phase_map` 仍然要求 DM 输出里是裸
`HLA..HLD`,oasis-dm 写的 `DM_HLADM..DM_HLDDM` 由 `oasis_dm.write_schema`
改名并按 `hypot(A,B)` 重算 FOM。

---

## 2. 验收:胰岛素 2.1 Å 的 A/B

同一份输入(`runs/ins_internal_bp3sites/`),同一套 Buccaneer+REFMAC 建库
(`ccp4_pipeline.run`,5 周期,`free_value=0`),只换 DM 引擎。

### 2.1 相位质量(evaluate.py,4509 反射,work 4291 / free 218)

| 变体 | wMPE | mean cosΔφ | σA C | 图 CC(全胞) | 图 CC(高密度) | 平均 FOM |
| --- | --- | --- | --- | --- | --- | --- |
| 初始 HL | 55.43° | 0.4791 | 0.5582 | 0.3646 | 0.3949 | 0.3214 |
| oasis-dm v2 默认 | 37.15° | 0.6987 | 0.7992 | 0.7856 | 0.7841 | 0.7738 |
| oasis-dm v2 cmbwt 0.3 | 36.80° | 0.7040 | 0.8049 | 0.7948 | 0.7938 | 0.7339 |
| oasis-dm v2 cmbwt 0.5 | **35.98°** | **0.7116** | **0.8107** | **0.8049** | **0.8016** | 0.8213 |
| 旧引擎 20 周期 | 34.29° | 0.7302 | 0.8498 | 0.8456 | 0.8534 | — |
| 旧引擎 40 周期 | 34.32° | 0.7301 | 0.8439 | 0.8435 | 0.8529 | — |
| 旧引擎 80 周期 | 34.23° | 0.7300 | 0.8435 | 0.8461 | 0.8569 | — |

旧引擎的相位仍然更好(CC 0.8456 vs 最好 0.8049,wMPE 差 1.7°),10 个壳层
里赢 8 个,最大差距在 2.6–2.8 Å(5–6.5°)。**oasis-dm 的价值是独立与自包含,
不是精度更高。**

### 2.2 下游建库

| 变体 | R_work / R_free | 残基 | CA<1 Å | CA RMSD |
| --- | --- | --- | --- | --- |
| oasis-dm 默认(0.4) | 0.2888 / 0.3295 | 52 | 50/51 | 0.21 Å |
| oasis-dm cmbwt 0.3 | **0.2634 / 0.3012** | 52 | 50/51 | 0.21 Å |
| oasis-dm cmbwt 0.5 | 0.4006 / 0.4217 | — | — | — |
| 旧引擎 20 周期 | **0.2433 / 0.3131** | 53 | 50/51 | 0.17 Å |
| 旧引擎 40 周期 | 0.3781 / 0.3896 | — | — | — |
| 生产链(3 轮旧引擎 + 2 轮反馈) | 0.2758 / 0.3181 | 50 | 48/51 | 0.22 Å |
| 生产链(3 轮 oasis-dm + 2 轮反馈,`--dm-engine oasis-dm`) | **0.2131 / 0.2662** | 52 | 51/51 | 0.12 Å |

最后一行的两组生产数字来自**同一条 `pplus_internal` 命令、同一套输入**,只差
密度修饰引擎。oasis-dm 的 R_free 比旧引擎低 5.2 个百分点,建出 52 个 CA(51 个
真值内)。这与上表的结论不矛盾:生产链多了 trace/Fcalc 反馈与固定观测先验,
单次建库的 R 不能当引擎排序依据(见本节末与 §4.3)。

引擎层面的统计(10 个 oasis-dm 变体 vs 5 个旧引擎变体):R_free
0.3599 ± 0.043 vs 0.3449 ± 0.030,差值 0.015,合并标准误约 0.019 —— **不可区分**。
单次建库的 R 因子不能用来给引擎排序:同一张图换随机种子,中位相位差 0.027°
就能让 R_free 摆动 3 个百分点(见 §4.3)。

---

## 3. `--cmbwt` 是什么

`oasis_dm_v2.py:259-260`:

```python
HL_mod = column_stack([cmbwt * X * cos(phi_mod),
                       cmbwt * X * sin(phi_mod), 0, 0])
```

`X`、`phi_mod` 来自 DM 图;这个四元组加到实验 HL 上,交给下一轮和最终输出。
`cmbwt=1.0` 是完全信任 DM 相位,`0.4` 是 v1 标定出的常数:简化的 omit ΣQ
偏乐观(留出集只有 5%),直接全信会虚高 FOM,下游 REFMAC 会过度信任。

关键含义:**它只影响“下一轮”和“输出”的权重,几乎不影响本轮 `Fc_omit`
和本轮地图**。因此单次运行内部无法用它自己的自由集指标来选它,必须做多次
短跑标定。

---

## 4. `--cmbwt` 的评估手段:实测与结论

### 4.1 试过的方案:留出集 CC 标定扫描

`oasis_dm.calibrate_cmbwt()` 对 {0.3, 0.4, 0.5, 0.6, 0.8} 各跑 3 轮,
用**标定专用种子 7717** 重新划留出集(与生产运行的 20261023 不同,避免在
同一批反射上调参),按留出集振幅 CC 打分:

| cmbwt | 0.3 | 0.4 | 0.5 | 0.6 | 0.8 |
| --- | --- | --- | --- | --- | --- |
| 留出集 CC(3 轮) | 0.7963 | 0.7987 | 0.7966 | 0.7996 | 0.8009 |

全距 **0.0046(0.6%)**,5 个候选 × 3 轮共 53 秒。跑满生产长度(9–12 轮)后
留出集 CC 也只从 0.9339(0.3)走到 0.9360(0.5)。

### 4.2 为什么这个指标不行(结构性原因)

把两个事实放在一起:

- 留出集 CC 排序:**0.5 > 0.4 > 0.3**,差距 0.002;
- 下游建库 R_free 排序:**0.3 (0.3320) < 0.4 (0.3560) < 0.5 (0.4367)**,
  差距 10 个百分点以上。

两者方向相反。原因就是 §3 说的:地图由 `Fc + DM 密度` 决定,`cmbwt` 只缩放
写出去的那一组 HL,所以任何地图质量指标(CC、free CC、R_holdout)对它的
敏感度都接近零;而真正被它影响的量 —— 交给 REFMAC 的 FOM 是否诚实 ——
恰恰是这些指标测不到的。

相位指标反而随 cmbwt 单调变好(§2.1:wMPE 37.15° → 36.80° → 35.98°),
但相位更好不等于建库更好:0.5 的相位最准、FOM 最虚高,R_free 反而最差。

### 4.3 多次重复后的真实效应量

单次 DM + 5 周期建库,同一 cmbwt 换随机种子的散布(括号内为 R_free):

| cmbwt | 各次 R_free | 中位数 |
| --- | --- | --- |
| 0.2 | 0.3992 | 0.399 |
| 0.3 | 0.3012 / 0.3285 / 0.3355 / 0.4405 | **0.332** |
| 0.4 | 0.3295 / 0.3546 / 0.3574 / 0.3821 | 0.356 |
| 0.5 | 0.4217 / 0.4367 / 0.5048 | 0.437 |
| 0.6 | 0.4179 | 0.418 |

结论:`cmbwt` 有真实效应(0.5 以上明显有害),但 0.3 与 0.4 之间的差距
(约 2 个百分点)小于同一取值的种子散布(约 ±6 个百分点)。**单次运行选不出
0.3 还是 0.4,需要至少 4 次重复才能分辨。**

### 4.4 因此 `auto` 被改成“只有证据明确才覆盖”

`--dm-cmbwt auto` 保留,但增加一道闸门:扫描全距小于 `MIN_CMBWT_SPREAD=0.01`
时**不覆盖**默认 0.4,并在日志里说明“扫描是平的,保持 0.4”。这避免了用噪声
去覆盖一个标定过的常数。日志示例:

```
cmbwt sweep is flat (0.7963..0.8009, spread 0.0046 < 0.01); keeping 0.40
```

`pipeline_report.json` 的 `dm_info` 里保留 `cmbwt`、`cmbwt_auto`、
`cmbwt_table`(每个候选的留出集 CC 与实跑轮数)、`best_free_cc`。

生产链上这道闸门的实际表现(胰岛素,3 轮 DM + 2 轮反馈,每轮都重新标定):

| 轮 | 扫描留出集 CC(0.3/0.4/0.5/0.6/0.8) | argmax | 全距 | 选中 |
| --- | --- | --- | --- | --- |
| 1 | 0.8355 / 0.8443 / 0.8480 / 0.8482 / 0.8430 | 0.6 | 0.0127 | 0.4 |
| 2 | 0.8118 / 0.8157 / 0.8175 / 0.8178 / 0.8164 | 0.6 | 0.0060 | 0.4 |
| 3 | 0.8483 / 0.8522 / 0.8544 / 0.8560 / 0.8601 | 0.8 | 0.0118 | 0.4 |

第 1、3 轮的全距其实过了 0.01 的闸门,但按“从 0.4 往上、落在最优 0.01 以内
就停”的省俭规则仍然选回 0.4。**argmax 在 0.6/0.6/0.8 之间游走**,正是 §4.2
说的盲区;而 §4.3 已经知道 0.5–0.8 在建库上是有害的。所以这次 `--dm-cmbwt auto`
跑出来的 `03_dm/map_phases.mtz` 与固定 `0.4` 的那次逐字节相同(sha256
`3913681a…`),R 也同为 0.2131 / 0.2662 —— 闸门没有帮倒忙,但也没带来收益。

### 4.5 建议

- **默认 0.4,不要动**。它保守、可复现,且在 0.3–0.4 这个平坦区间里。
- 想要更好的建库结果,用 `--dm-cmbwt 0.3`(本案例 R_free 0.332 vs 0.356),
  代价是相位略差、FOM 略低。
- **不要**用 `auto` 当作“自动调参”卖点:它的判据只在地图确实被 cmbwt 影响时
  才有分辨力,而那不是当前引擎的行为。真正的选法只有两种:固定 0.4,或者
  用端到端建库做多次重复扫描(每次约 1 分钟 DM + 1 分钟建库,至少 4 个种子)。
- 若以后要让 `auto` 真正有效,需要给 `oasis_dm` 加一个**相位敏感**的留出集
  判据,例如:留出集上 `Fc_omit` 与观测的振幅拟合只给出 ΣQ,拿不到相位误差;
  可行方向是把留出集反射从 DM 中彻底剔除后,用两套独立起始相位各跑一次,
  比较留出集上的相位一致性(半集一致性)。这属于新研究项,不在本次范围内。

---

## 5. 复现命令

同一条 `pplus_internal` 命令换 `--dm-engine` 就能复现上面的生产对照（输入是
`/Volumes/For_data/IPCAS-test/ins/` 下的胰岛素 2.1 Å SAD 数据）:

```bash
INS=/Volumes/For_data/IPCAS-test/ins
REF=$PWD/validation/unified_verified_tmp/synthetic-reference0/reference-1tqw.mtz
./.conda-unified/bin/python -u run_sad_to_model.py \
  --sad $INS/ins_wave2.0_reso2.1A_sad.mtz --native $INS/ins_wave2.0_reso2.1A_sad.mtz \
  --sequence $INS/3INS.pir --heavy /private/tmp/bp3/bp3out-1.pdb \
  --reference $REF --fpp S=1.0 \
  --output /tmp/prod_hdm --method pplus_internal --dm-engine oasis-dm \
  --dm-cycles 3 --feedback-cycles 2 --dm-solvent 0.5 \
  --input-columns --no-auto-hand --build-model --dm-cmbwt 0.4
```

`--dm-cmbwt auto` 会把 `0.4` 换成标定结果;在本例里标定三次都选回 `0.4`
(见 §4.4),所以 `03_dm/map_phases.mtz` 与固定 `0.4` 的那次**逐字节相同**
(sha256 `3913681a…`),R 也相同。


### 5.1 单个案例(默认 oasis-dm)

```bash
./.conda-unified/bin/python run_sad_to_model.py \
  --sad SAD.mtz --native NATIVE.mtz --sequence SEQ.pir --heavy SITES.pdb \
  --reference REFERENCE.mtz --fpp S=1.0 \
  --method pplus_internal --dm-engine oasis-dm \
  --dm-cycles 3 --feedback-cycles 2 --dm-solvent 0.5 \
  --input-columns --no-auto-hand --build-model --output OUT
```

### 5.2 关于旧十案例表

`docs/STATUS_20261112.md` 那张十案例表是用已删除的旧引擎跑出的历史记录,
**不能再逐点复现**;现在的 `--dm-engine` 只接受 `oasis-dm`/`internal`,换引擎
后数值属于不同引擎。要重跑十案例,按 §5.1 的 oasis-dm 命令并重新标定终点。

### 5.3 只跑 DM 并看标定表

```bash
./.conda-unified/bin/python -c "
from phasebuild.oasis_dm import calibrate_cmbwt
w, table = calibrate_cmbwt('work.mtz', 'hl.mtz', 'cal', 0.5)
print(w, table)"
```

### 5.4 测试

```bash
PYTHONPATH=src:tests ./.conda-unified/bin/python -m pytest -p no:randomly -q
# 870 passed, 1 deselected
```

---

## 6. 已知限制

- 手写引擎没有各向异性校正(上游 oasis-dm 没有这个功能),`--dm-anisotropy`
  对它无效,`dm_info['anisotropy_correction']` 恒为 `False`。
- `--sigmaq-source omit` 在留出集少于 `min_free=500` 时会自动退回 `blend`
  并把 cmbwt 压回 0.4;胰岛素只有 218–225 个留出反射,所以该开关在本案例
  上是惰性的(与默认逐字节相同)。
- `--no-auto-stop` 会忽略 `--max-cycles`/`--min-cycles`,实际按上游
  `Config.cycles=5` 跑。这是上游 CLI 的不一致,已在 `oasis_dm` 外层记录,
  但没有改上游代码。
- 十案例表是用已删除的旧引擎跑出的;换成 oasis-dm 后属于不同引擎,数值不能逐点对比。
  目前只完成了胰岛素一个案例的 A/B,十案例的 oasis-dm 全量对照尚未做。
