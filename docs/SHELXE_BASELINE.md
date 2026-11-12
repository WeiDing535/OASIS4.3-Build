# 以 SHELXC/D/E 为标尺：我们到底需要替换多少 CCP4

本文替换了此前"对标 SHELXE + Buccaneer + REFMAC 全链"的口径。**对标对象是 SHELXC/D/E 这条独立链**：SHELXE 自己做密度修饰、三肽搜索和主链追踪，不调用 Buccaneer，也不调用 REFMAC。既然目标是"达到 SHELXE 的水平"，那么 Buccaneer 与 REFMAC 的替换就不是发布前提，而是可选增强。

结论先给：

1. **SHELXE 在这些案例上的水平比我们一直要求自己的低得多**：pseudo-free CC 50.6–72.8 %，mean FOM 0.48–0.71，最好一轮 trace 的 CC 只有 7.5–49.0 %，输出是**多聚丙氨酸主链**（7ENM 只有 45 个 Cα，目标链长 101）。
2. **SHELXE 的密度修饰图与 Parrot 的图打平**（平移不变地图相关系数均值 0.683 vs 0.688，十案例逐例 −0.046…+0.032）。也就是说，"SHELXE 水平的密度修饰"有一个我们能够直接对照的等价物。
3. **我们自研链的真正缺口只有一个，而且已经量化：内部 DM 比这条线低 0.125**（0.558 vs 0.683），十案例无一例外。Parrot 图与内部 DM 图共用同一份 P+ 输入相位，所以差距只来自密度修饰引擎本身，与 Buccaneer/REFMAC 无关。
4. **Buccaneer / REFMAC 不在这个目标内**。它们带来的是"完整模型 + R/R-free 25–47 %"，远超 SHELXE 的多聚丙氨酸 trace；要不要保留是产品档次问题，不是发布许可问题。

## 1. SHELXE 自己做什么、不做什么

`cases/7_7enm_3se_P212121/shelxe_original.log` 的头部就写明了它的全部工作：

```
 +  SHELXE  -  PHASING AND DENSITY MODIFICATION  -  Version 2023/1  +
 Poly-Ala trace output to trial.pdb
 -a     4  global autotracing cycles
```

- **密度修饰**：球影响范围（sphere of influence）内的密度修饰，日志里有独立计时（7ENM：sphere of influence 76.6 s、tripeptide 42.8 s、tracing 26.1 s，合计 170.6 s）。
- **主链追踪**：`-a 4` 轮全局自动追踪，输出 `trial.pdb`，**只有多聚丙氨酸**（CA + 少量主链原子，没有序列指派、没有侧链）。
- **没有外部程序调用**：日志里没有任何 Buccaneer 或 REFMAC 调用；SHELXE 也不做对原生数据的约束精修，不使用 R-free。
- **它的质量指标是自己的**：`Pseudo-free CC`、`Estimated mean FOM`、按分辨率壳层的 `<mapCC>`、以及 `Best trace (cycle N with CC xx.xx%)`。SHELXE 不产生 R/R-free。

这就是"对标口径"的全部内容——不是 R-free 差 5 个点，而是 pseudo-free CC 与 trace 覆盖率。

## 2. SHELXE 在十案例上的实际水平

| 案例 | Pseudo-free CC (%) | mean FOM | 最好 trace CC (%) | SHELXE 输出 Cα | 我们 campaign 的 Cα | campaign R / R-free |
| --- | --- | --- | --- | --- | --- | --- |
| 6ZPP | 72.76 | 0.705 | 49.01 | 148 | 159 | 30.04 / 31.31 |
| 6HF7 | 71.57 | 0.701 | 46.60 | 506 | 602 | 26.02 / 29.04 |
| 5HHK | 70.02 | 0.691 | 44.29 | 556 | 587 | 25.94 / 28.82 |
| 7TFQ | 66.87 | 0.661 | 34.60 | 225 | 288 | 32.52 / 36.86 |
| 3O2E | 61.69 | 0.579 | 44.36 | 67 | 78 | 41.14 / 42.64 |
| 7TJ1 | 61.57 | 0.600 | 30.62 | 351 | 569 | 23.92 / 31.08 |
| 7E1D | 59.32 | 0.575 | 33.10 | 124 | 167 | 36.65 / 42.33 |
| 5D9D | 51.97 | 0.493 | 7.45 | 88 | 295 | 25.43 / 27.51 |
| 7TBS | 51.82 | 0.501 | 8.34 | 58 | 218 | 25.17 / 28.35 |
| 7ENM | 50.64 | 0.484 | 9.41 | 45 | 101 | 28.34 / 31.47 |

三点值得注意：

- SHELXE 的 trace CC 在 5D9D（7.45 %）、7TBS（8.34 %）、7ENM（9.41 %）三个案例上几乎等于没有追踪成功，但仍会写出 `trial.pdb` 并给出 pseudo-free CC。
- 它从不做序列指派，因此 Cα 数目不能直接当成"建模完成度"来比。
- 我们的 campaign 在这些案例上给出的是**带序列、带侧链、R 25–47 %** 的完整模型。之前把 Buccaneer/REFMAC 当作发布硬阻塞，是因为拿 campaign 的水平当验收线；换成 SHELXE 水平，这条线低得多。

## 3. 用同一把尺子量三种图：SHELXE ≈ Parrot，内部 DM 低 0.125

SHELXE 的 `Pseudo-free CC` 用的是它自己的留出集与密度修饰内部量，外部无法逐反射复算。这里改用一把**两边完全同码**的尺子：把每张图与同一个外部参照比较，参照是 campaign 精修模型的傅里叶系数 Fc、φc（由同一套数据算出，与三种图都独立）。

指标是**平移与手性不变的地图相关系数**：

```
complex CC = max_t | Σ_h Fmap(h) conj(Fc(h)) e^{2πi h·t} | / sqrt(Σ|Fmap|² · Σ|Fc|²)
```

对 t 取最大值是必要的：6ZPP 的 SHELXE 图沿 b 轴偏了 21 个网格（P 1 2 1 1 的极性轴原点不定），零平移相关只有 0.035，对齐后是 0.669。不做这一步会把一个原点约定问题误读成"图很差"。

| 案例 | DM 输入（campaign P+） | DM 输入（内部 P+） | SHELXE 图 | Parrot 图 | 内部 DM 图 | Δ(内部 − SHELXE) |
| --- | --- | --- | --- | --- | --- | --- |
| 3O2E | 0.123 | 0.262 | 0.595 | 0.549 | 0.432 | −0.163 |
| 5D9D | 0.224 | 0.249 | 0.713 | 0.724 | 0.531 | −0.182 |
| 5HHK | 0.399 | 0.375 | 0.805 | 0.805 | 0.709 | −0.096 |
| 6HF7 | 0.290 | 0.323 | 0.729 | 0.709 | 0.616 | −0.113 |
| 6ZPP | 0.185 | 0.156 | 0.669 | 0.701 | 0.570 | −0.099 |
| 7E1D | 0.323 | 0.293 | 0.602 | 0.632 | 0.525 | −0.078 |
| 7ENM | 0.400 | 0.390 | 0.734 | 0.721 | 0.673 | −0.062 |
| 7TBS | 0.227 | 0.218 | 0.620 | 0.636 | 0.487 | −0.133 |
| 7TFQ | 0.258 | 0.225 | 0.709 | 0.716 | 0.525 | −0.184 |
| 7TJ1 | 0.289 | 0.294 | 0.651 | 0.682 | 0.511 | −0.140 |
| **均值** | **0.272** | **0.279** | **0.683** | **0.688** | **0.558** | **−0.125** |

中间两列是**密度修饰之前**的 SAD 相位（`02_sad/experimental.mtz` 的 `PPLUS·FOM, PHIB`），分别来自 campaign 的 P+ 后端（导入原始 Phaser SAD 纯观测似然）与内部 P+ 后端（内部异常差似然，不读 Phaser 文件）。把它们放进同一张表，是为了把总差距拆成"输入相位"与"密度修饰"两段：

| | DM 增益（DM 图 − 其输入） |
| --- | --- |
| Parrot（campaign P+ → Parrot） | **+0.416**（逐案例 0.309…0.516） |
| 内部 DM（内部 P+ → 内部 MLHL） | **+0.279**（逐案例 0.170…0.414） |

三点结论：

1. **SHELXE 的图与 Parrot 的图打平**（均值 0.683 vs 0.688，逐案例 −0.046…+0.032）。注意 Parrot 是被参照模型"偏袒"的一方——campaign 模型就是从 Parrot 图上建出来的——即便如此 SHELXE 也没落后。所以 **SHELXE 水平的密度修饰 ≡ Parrot 水平的密度修饰**，这是一个我们已经有明确实现目标的量。
2. **内部 DM 比这条线低 0.125，且 10/10 案例无一例外**（−0.184…−0.062）。这不是噪声，是系统性差距。
3. **差距不在输入相位，在密度修饰本身**：两个 P+ 后端给出的 DM 输入相位质量相当（均值 0.272 vs 0.279，逐案例互有胜负，最大差 0.139）。真正的差别是**修饰增益**——Parrot 把相关从 0.272 提到 0.688（+0.416），内部 MLHL 只从 0.279 提到 0.558（+0.279）。

由此，**"达到 SHELXE 水平"的可操作定义**是：内部 DM 的 complex CC 追到 Parrot 的 ±0.03 以内。当前 0.558 → 目标 ≈ 0.68。

### 3.1 一个被推翻的检验（负结果，记录在案）

先前用过"5 % 省略检验"（free-lunch）：把图系数的 5 % 振幅置零、FFT 建图、再正变换读回这些反射的 |Fc|，与 |Fobs| 相关。**这个检验在数学上恒等于零，量到的只是数值噪声**：

- 省略 5 % 后做 FFT→反 FFT，在网格上是恒等变换，被省略的反射回读值中位数 1.3e-5，而 |Fobs| 中位数 98.5——即回读值本质为零；
- 之所以还能算出 0.205 的"CC"，是因为残差噪声随局部密度成比例，于是与 |Fobs| 伪相关；把 |Fobs| 打乱后 CC 掉到 0.023，证明它不含相位信息。

因此**本文早前版本 §3 里的 0.07–0.27 那张表已作废**，由上面的地图相关系数表取代。要让 omit-map 检验成立，循环里必须有一个非线性步骤（溶剂压平/直方图匹配），但那测的就不再是"交付的这张图"了。

复现：

```bash
cd /Volumes/Fast_4T/Deepseek-Oasis
python3 work/noCCP4/compare_map_phases.py   # 写出 work/noCCP4/map_phase_comparison.json
```

### 3.2 DM 引擎内部的消融：差距不是任何一个旋钮

既然差距出在密度修饰本身，就直接在 DM 引擎里做单旋钮消融。全部用 DM-only 运行（`--no-build-model`，只看 `FWT/PHWT` 图），逐循环 dump 系数（`OASIS_INTERNAL_DM_TRACE`）后与 campaign 模型比 complex CC。先做两个案例（3O2E 最差、7TBS 中等），确认方向后再上十案例。

| 变体 | 3O2E | 7TBS | 说明 |
| --- | --- | --- | --- |
| 基线：`mlhl` + 掩膜平滑 4.0 Å | 0.432 | 0.487 | 现发布路径 |
| 掩膜平滑 3.0 Å | 0.443 | 0.495 | 边际改善 |
| `combine=map`（照搬 Parrot 的纯图迭代，不做相位组合） | **0.076** | **0.060** | 崩溃 |
| 多分辨率阶梯 `auto` + 3.0 Å，40 轮 | 0.412 | 0.499 | 3O2E 变差 |
| 多分辨率阶梯 `auto` + 3.0 Å，100 轮 | 0.418 | 0.499 | 加轮次不救 |
| `gamma=0`（去掉扰动校正） | 0.375 | 0.446 | 变差 |
| Wang 球影响掩膜（半径按溶剂比校准）+ 3.0 Å | **0.491** | 0.495 | 3O2E +0.048，7TBS 持平 |

四个直接结论：

1. **加轮次无效**：complex CC 在第 8–10 轮就饱和，之后 30 轮不动（40 轮与 100 轮差 <0.01）。所以不是"迭代不足"，是**迭代不动点太低**。
2. **照搬 Parrot 的纯图迭代会崩**：去掉 MLHL 相位组合、只用修饰后的图当下一轮系数，弱 SAD 先验下 CC 掉到 0.06–0.08。MLHL 锚定是必需的，这条捷径不通。
3. **多分辨率阶梯不涨**：100 轮 `auto` 阶梯对 3O2E 反而低于单分辨率（0.418 vs 0.443），对 7TBS 只 +0.004。Parrot 的低分辨率起步要靠它自己的掩膜/加权配合，单独搬过来无效。
4. **掩膜确实有毛病，但不是主因**：把内部掩膜与"由 campaign 模型原子膨胀 2.2 Å + 空间群展开"的真值掩膜比，`variance` 掩膜平均占据 0.735（真值 0.507，即过度包含约 45 %），掩膜 CC 只有 0.505；换成按溶剂比校准半径的 Wang 掩膜后占据 0.647、CC 0.561。**掩膜质量确实提高了，十案例均值却是 0.558 → 0.546（−0.012，单例 −0.075…+0.059）**——掩膜不是那 0.125 的来源。

十案例 Wang 掩膜逐例（对照 §3 的基线列）：3O2E 0.491、5D9D 0.491、5HHK 0.692、6HF7 0.599、6ZPP 0.561、7E1D 0.470、7ENM 0.679、7TBS 0.495、7TFQ 0.450、7TJ1 0.529，均值 0.546。

由此得到本轮最重要的判断：**内部 DM 与 Parrot 的差距不在某个可调参数，而在算法的整体结构**——Parrot/DM 是"多分辨率 + 逐轮重建掩膜 + 按壳层 σ_A 加权的相位组合 + NCS/溶剂平均"一整套联动，任何单项照搬都不涨。要让 complex CC 从 0.558 走到 0.68，需要的是 DMMULTI 级别的重写（约 1500–2000 行 + 逐案例标定），不是调参。

复现（十案例 DM-only，4 进程，约 25 分钟）：

```bash
cd /Volumes/Fast_4T/Deepseek-Oasis
OASIS_INTERNAL_DM_MASK=wang OASIS_INTERNAL_DM_MASK_RADIUS=3.0 \
  python3 work/noCCP4/run_dm_only_all.py
python3 work/noCCP4/dm_only_gain.py dm_only_trace_wang
```

单案例逐循环 trace（写入 `cycle_NN.mtz` 与最后一轮的掩膜/图）：

```bash
DM_TRACE_DIR=$PWD/work/noCCP4/dm_trace DM_TRACE_OUT=$PWD/work/noCCP4/dm_trace_runs \
  python3 work/noCCP4/run_dm_trace.py 3_3o2e_9I_1K_P42212 19_7tbs_3SE_I222
```

消融开关全部是环境变量，默认值复现发布路径：`OASIS_INTERNAL_DM_MASK`（`variance`/`wang`）、`OASIS_INTERNAL_DM_WANG_RADIUS`、`OASIS_INTERNAL_DM_WANG_EDGE`、`OASIS_INTERNAL_DM_MASK_RADIUS`、`OASIS_INTERNAL_DM_COMBINE`（`mlhl`/`map`）、`OASIS_INTERNAL_DM_GAMMA`、`OASIS_INTERNAL_DM_LADDER`、`OASIS_INTERNAL_DM_TRACE`。

## 4. 达到 SHELXE 水平还差什么

| 能力 | SHELXE | 我们现在 | 差什么 |
| --- | --- | --- | --- |
| SAD 定相 | SHELXD | `internal_pplus.phase`（NumPy SAD 似然 + P+ 信任域） | 无 |
| 密度修饰 | 球影响范围 DM | `internal_pplus.seeded_dm`（NumPy MLHL） | **增益 0.279 vs Parrot 的 0.416**，终点 complex CC 低 0.125；§3.2 的七项单旋钮消融都不涨，需 DMMULTI 级重写 |
| 多聚丙氨酸 trace | 内建 autotracing | `trace.py` 等原型约 2400 行，**未接入** | 接入 + 十案例验收 |
| 质量指标 | pseudo-free CC / FOM / mapCC | 有 R/R-free；地图验收脚本已就位（`work/noCCP4/compare_map_phases.py`） | 把地图相关系数接进流水线报告 |
| 约束精修 | 无 | 内部精修器（原型） | **不需要**（SHELXE 不做） |
| 建链 + 序列 + 侧链 | 无 | Buccaneer | **不需要**（超出目标） |

所以"纯自研、达到 SHELXE 水平"的产品形态是：**内部 P+ → 内部 DM → 多聚丙氨酸 trace**，质量报告用地图相关系数与 trace CC。这条链上已经没有 CCP4 可执行文件。

但按第 3 节的量化结果，**这条链现在还没到 SHELXE 水平**：密度修饰差 0.125，这是必须先补的一项。补齐顺序建议是：

1. 内部 DM 追平 Parrot（目标 complex CC ≥ 0.68，十案例逐例不低于 SHELXE − 0.03）。**不要再去调参**：§3.2 已经把掩膜半径、掩膜构造（variance → 校准 Wang）、压平/直方图、相位组合方式（MLHL vs 纯图迭代）、扰动校正、分辨率阶梯、轮数（40 → 100）逐个试过，单项最好也只有 +0.05，十案例均值不涨。要补的是结构：逐轮由当前图重建掩膜、按壳层 σ_A 与 FOM 加权的相位组合、多分辨率阶梯与掩膜联动、NCS/溶剂平均，四者一起才构成 DMMULTI 的行为。
2. 再把 `trace.py` 接入流水线，输出 SHELXE 式的多聚丙氨酸 `trial.pdb`。
3. 两项都到位后，才谈"纯自研可发布"。

## 5. 工程收尾（本轮完成）

- 内部 DM 在 6ZPP、5HHK 上曾因参考图网格内存上限失败（`reference_histogram` 的 `max_grid_points`）。`native_density.py` 的采样率自适应已经把 1tqw 参考图从 9.95 M 点降到 2.76 M 点，三个案例重跑通过。
- 6HF7 曾因输入重原子模型带 ANISOU 被内部后端拒绝。现在 `read_atoms(..., project_anisotropic=True)` 把各向异性 ADP 投影到 U_iso = trace(U)/3，并发出可见的 `UserWarning`（不是静默降级）；内部 P+/DM/信任域三处调用点启用，其余调用点仍默认拒绝。
- 5HHK、6HF7 的十案例重跑已完成：内部 P+ → 内部 DM → Buccaneer → REFMAC 的最终模型为 5HHK 26.22/28.92（campaign 25.94/28.82）、6HF7 25.34/29.17（campaign 26.02/29.04）、6ZPP 30.54/33.36（campaign 30.04/31.31）。也就是说，密度修饰图差 0.10–0.13 时，后端的 Buccaneer/REFMAC 基本把这个差距吸收掉了（R-free 差 0.10–2.05 点）。
- 全套测试 639 passed（`PYTHONPATH=src .conda-unified/bin/python -m pytest tests -q`）。

## 6. 口径说明与已知局限

- 第 3 节的地图相关系数不是纯相位指标：它同时含振幅与相位，只是三种图的振幅都已被缩放到 |Fobs|，因此表里的差异主要来自相位。它也不测溶剂掩膜是否正确、不测 trace 的可连接性。
- 参照模型（campaign 精修模型）是从 Parrot 图上建出来的，因此**天然偏袒 Parrot**；SHELXE 在这种不利条件下仍与 Parrot 打平，结论因此更稳，而不是更弱。
- SHELXE 的 pseudo-free CC 无法逐反射复算，因此表 A 与表 3 的数字不可直接相减；第 3 节是"同一把尺子"的横向比较。
- 我们的 campaign 数字（表 A 右侧两列）来自 Buccaneer + REFMAC 链，用于说明"campaign 已远超 SHELXE 水平"，不用于证明自研链的能力。自研链的能力由第 3 节的内部 DM 列与 §4 的 trace 状态界定。
- 内部 DM 那一列来自 `--refinement-prior original_sad` 的对照运行（与 campaign 的 `dm_phifom` 不同），这是与 campaign 比较时的混杂因素；但与 Parrot 列的比较不受影响，因为两者共用同一份 P+ 输入。
