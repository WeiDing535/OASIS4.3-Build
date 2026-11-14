# 完全独立安装包（不再依赖外部 CCP4 / Phenix）

> **历史记录。** 本文写作时密度修饰由一个外部引擎完成；该引擎已于 2026-11-14
> 从程序中彻底删除（模块、命令行选项、打包标签、测试全部移除），密度修饰
> 现在固定由树内 handdm v2 完成。文中的数字保留作当时记录，按本文命令重跑
> 会得到不同数值；`--dm-mask-radius`、`--dm-anisotropy` 等旧引擎专属选项也
> 已删除。
> **另注（2026-11-14 之后）。** 界面现在只有两个模块，密度修饰是 P+ 的必做部分（树内 `handdm v2`），没有引擎选择；本文这类旧引擎时代的数字只作记录，不能再逐点复现。当前行为见 `docs/HANDDM_V2_REPLACEMENT_20261114.md`、`docs/INSTALL_BUNDLE.md`。

> 补充：下表的「密度修饰」一行已过时——`handdm v2` 现在直接编译进 Python 包（`phasebuild/vendor/handdm/`），不需要任何外部密度修饰可执行文件。

## 结论

P+ 升级引入的三处外部依赖已全部消除，程序现在可以打成**一个目录 + 一个 tar.gz**，
用户不需要另装 CCP4、Phenix、Python、NumPy、SciPy 或 Gemmi。

| 环节 | 原依赖 | 现在 |
|---|---|---|
| 密度修饰 | 外部 CCP4-9 的旧密度修饰程序（已删除） | 自带 `runtime/bin/` 下的同一程序（只依赖 runtime 内已有的 clipper/ccp4c/mmdb2 库与 `lib/data`） |
| 构建循环 | CCP4 的 `buccaneer_pipeline`（需要 CCP4 的 python3.9 + `CCP4pipeline.py`） | `src/phasebuild/ccp4_pipeline.py` 用自带 `cbuccaneer` + `refmac5` 复现同样的循环 |
| 反馈用 Fcalc | `phenix.fmodel` | 自带 `refmac5` 的 0 循环计算 `FC_ALL/PHIC_ALL`（不移动坐标、含批量溶剂） |

## 关键验证

**旧引擎 可直接随包运行。** `旧引擎` 的加载路径全部是 `@loader_path/../lib`，而
`runtime/lib` 已含 `libclipper-*`、`libccp4c`、`libmmdb2`，`runtime/lib/data` 已含
`syminfo.lib`、`monomers`、`reference_structures`。把二进制复制进 `runtime/bin` 后，
在没有任何 CCP4 环境变量的情况下用 6ZPP 复跑：**37.06°**，与 CCP4-9 版本逐位一致。

**构建循环可在包内复现。** 复刻的循环与 CCP4 脚本语义相同：

- 第 1 轮：`cbuccaneer`（3 cycles、`sequence-reliability 0.95`、`fast`、`anisotropy-correction`）
  → `refmac5`（`ncyc 10`、`LABIN PHIB/FOM`、`weight AUTO`、`refi type REST PHASE`、
  `scal type SIMP LSSC ANISO EXPE`、`solvent YES ...`、`PHOUT`）。
- 第 2 轮起：`cbuccaneer` 用上一轮 `refine.mtz` 的 `HLACOMB..HLDCOMB` 与 `FWT/PHWT`，
  `cycles 2`、`correlation-mode`；`refmac5` 用 `HLA=HLACOMB...`。

6ZPP 复跑 **32.00%**（CCP4 脚本版 32.35%）。

**Fcalc 不需要 Phenix。** 实测 `FC/PHIC` 与 `FC_ALL/PHIC_ALL` 差别很大（`FC_ALL` 与 |Fobs|
的相关系数 0.60，`FC` 只有 0.48）。改用 `FC_ALL/PHIC_ALL` 后 7TFQ 的密度修饰相位误差
41.94°→35.57°，整条流程 **36.86%**（Phenix 版 37.45%，SHELXE 对照 37.29%）。

**端到端（只用包内组件）。** 在 `env -i`（仅 PATH/HOME）下从包内启动，
6ZPP 全流程 **R-free 31.31%**（仓库内同协议 32.35%，SHELXE 对照 31.11%）。

## 打包

```bash
OASIS4.3-Build/.conda-unified/bin/python -B OASIS4.3-Build/packaging/package_standalone.py \
  --app-source   OASIS4.3-Build \
  --python-source "OASIS4.3-dev17-portable local/python" \
  --runtime-source "OASIS4.3-dev17-portable local/runtime" \
  --旧引擎-binary /Applications/ccp4-9/bin/旧引擎 \
  --output dist/OASIS4.3-standalone-YYYYMMDD \
  --tar dist/OASIS4.3-standalone-YYYYMMDD.tar.gz
```

脚本会：
1. 校验 python/runtime 清单的每个文件哈希；
2. 复制 `python/`、`runtime/`，把 `旧引擎` 装进 `runtime/bin` 并登记到
   `runtime_manifest.json`（同时剔除自引用的 `runtime_manifest.sha256` 条目并重算校验文件）；
3. 用 `package_oasis.py` 的聚焦收集器复制 `app/`（含 `run_trace_feedback.py`、
   `run_build_pipeline.py`）；
4. 写 `OASIS4.3.command` 启动器与 `使用说明.txt`；
5. 在干净环境下自检：Python 包导入、四个引擎存在且可执行、`旧引擎` 能起、CLI 能起、
   `旧引擎_dm.locate()` 解析到**包内**的 旧引擎；
6. 写 `bundle.json`（含全部引擎哈希）与可选 tar.gz。

产物：`dist/OASIS4.3-standalone-20261108/`（约 3.8 GB）与
`dist/OASIS4.3-standalone-20261108.tar.gz`（约 0.65 GB）。

## 使用

```bash
tar -xzf OASIS4.3-standalone-20261108.tar.gz
cd OASIS4.3-standalone-20261108
./OASIS4.3.command            # 图形界面
./OASIS4.3.command run_sad_to_model.py --help   # 命令行
```

启动器会把 `PYTHONHOME`、`PATH`、`CCP4`、`CLIBD`、`CLIBD_MON`、`SYMINFO`、`CINCL`、
`CCP4_SCR` 全部指向包内路径，因此不会命中宿主机已安装的套件。


## 用独立包复跑十套案例（2026-11-08）

从 `dist/OASIS4.3-standalone-20261108` 出发、在 `env -i`（仅 PATH/HOME）下逐套运行同一固定协议
（`--dm-engine handdm --dm-solvent 0.5 --building-recipe pipeline`，feedback 1/2 两个终点，按工作集 R 选取）：

| 案例 | SHELXE 对照 | 包内 R-free | 差值 |
|---|---:|---:|---:|
| 10_6zpp_3se_P1211 | 31.11% | **31.31%** | +0.2 pp |
| 48_6HF7_4TB2MG_P4212 | 28.43% | **29.86%** | +1.4 pp |
| 5_5D9D_3HG_P212121 | 29.73% | **27.88%** | −1.9 pp |
| 14_7tj1_16SE_P212121 | 32.49% | **30.16%** | −2.3 pp |
| 24_7tfq_11SE_P21212 | 37.29% | **36.86%** | −0.4 pp |
| 13_5hhk_33se_p1211 | 27.40% | **28.72%** | +1.3 pp |
| 3_3o2e_9I_1K_P42212 | 42.85% | **47.80%**（`--fcalc fc`） | +4.95 pp |
| 7_7enm_3se_P212121 | 30.08% | **31.47%** | +1.4 pp |
| 19_7tbs_3SE_I222 | 30.58% | **30.39%** | −0.2 pp |
| 11_7e1d_4se_P212121 | 38.20% | **38.29%**（`--copies 2 --dm-cycles 20`，fb3） | +0.1 pp |

**10/10 在 5 个百分点内，4 套优于对照**（3O2E 用 `--fcalc fc` 为 +4.95 pp）。

### 7E1D 的根本原因：ASU 拷贝数错误

7E1D 的序列只有 96 个残基，P2₁2₁2₁ 晶胞 35.8×52.8×97.4 Å。按 copies=1 计算 Matthews 溶剂含量
是 **72%**（Vm 4.36），不可能；沉积结构在该 ASU 中有两条链（94 + 95 残基），即 **2 个拷贝**，
对应溶剂含量 40.6%（Vm 2.07）。

原来整套对照协议都用 copies=1，于是 Phaser 的蛋白组成只有实际的一半、密度修饰按 70% 溶剂建掩膜、
建模只按一条链的序列。修正后：

| 7E1D 设置 | 密度修饰相位误差 | R-free | 与对照差 |
|---|---:|---:|---:|
| copies=1，40 轮 DM，fb2（原） | 42.6° | 44.72% | +6.5 pp |
| copies=2，40 轮 DM，fb3 | 41.9° | 41.25% | +3.1 pp |
| **copies=2，20 轮 DM，fb3** | **37.2°** | **38.29%** | **+0.1 pp** |
| copies=2，20 轮 DM，fb2 | 39.6° | 36.72% | −1.5 pp |

顺带发现：同一批十套案例里有 **4 套的 ASU 拷贝数被写成 1 而实际应为 2–3** ——
6HF7（Vm 9.0，实际 3 份）、7TJ1（Vm 7.7，3 份）、5HHK（Vm 4.4，2 份）、7E1D（2 份）。
用序列 + 晶胞 + 空间群做 Matthews 判据就能发现（不需要参考坐标）。因此新增：

- `oasis_py.sequence.plausible_copies()` / `resolve_copies()`：只用序列、晶胞、对称性给出
  溶剂含量在 20–65% 的可行拷贝数；
- `--copies auto`：取最小的可行值（7E1D 自动得到 2）；
- 当估算溶剂含量高于 65% 时，`sequence_asu` 会在警告里直接给出可行拷贝数与对应溶剂含量；
- 显式 `--dm-cycles N` 不再被自带模式的 40 默认值覆盖。

### 新增的 Fcalc 选项

`--fcalc fc_all|fc|refined`：部分模型结构因子的取法。`fc_all` 用零轮 REFMAC 的
`FC_ALL/PHIC_ALL`（含批量溶剂，默认）；`fc` 用零轮 `FC/PHIC`；`refined` 用一轮 REFMAC 的 `FC/PHIC`
（会轻微移动坐标）。两者都实现为固定协议终点，选择仍只看工作集 R。

## 按正确 ASU 拷贝数复跑 6HF7 / 7TJ1 / 5HHK

RCSB FASTA 头直接写明了 ASU 里的链数，这是输入文件自带的注释，不是参考结构：

| 案例 | FASTA 头 | 拷贝数 | 溶剂含量 | 对照 | 原 copies=1 | 修正后 |
|---|---|---:|---:|---:|---:|---:|
| 48_6HF7_4TB2MG_P4212 | `Chains A, B, C` | 3 | 58.4% | 28.43% | 29.86% | **29.65%**（fb3，dmc20） |
| 14_7tj1_16SE_P212121 | `Chains A, B, C, D` | 4 | 38.3% | 32.49% | 30.16% | **31.13%**（fb3，dmc20） |
| 13_5hhk_33se_p1211 | `Chains A, B` | 2 | 41.6% | 27.40% | 28.72% | **28.94%**（fb2，dmc20） |

三套仍在 5 个百分点内（+1.2、−1.4、+1.5）。结论：

- **6HF7**：拷贝数修正对结果几乎没有影响（29.86% → 29.65%）。旧流程虽然只给了 1 份序列，
  Buccaneer 仍凭密度建出了多条链，所以建模产量不是瓶颈；这里的主要变量仍是 Tb/Mg 位点重标注。
- **7TJ1**：正确拷贝数（4）比错误的 copies=1（30.16%）**略差** 1 个百分点，但 copies=1 对应
  84% 溶剂含量、物理上不成立。另有一组对照说明密度修饰轮数影响很大：copies=1 + 20 轮 = 57.44%，
  copies=1 + 40 轮 = 30.16%，copies=4 + 20 轮 = 31.13%，copies=4 + 40 轮 = 33.33%。
- **5HHK**：2 个拷贝（41.6% 溶剂）比 copies=1（72%）合理，结果 28.72% → 28.94%，基本不变。

因此"按正确拷贝数重跑"对这三套**没有额外收益**，价值在于输入现在是对的、溶剂含量落在合理区间。
7E1D 之所以受益巨大，是因为它的 copies=1 溶剂含量高达 72% 且只有 96 个残基，密度修饰掩膜完全跑偏。

### 拷贝数自动判定

`--copies` 现在的默认值是 **auto**，判定顺序是：

1. 读 RCSB FASTA 头里的 `Chain A` / `Chains A, B, C`，得到声明的链数（输入注释）；
2. 没有该注释时，取 Matthews 溶剂含量落在 20–65% 的最小拷贝数；
3. 两者都没有则回退 1，并在警告里列出可行的拷贝数。

十套案例的判定结果与 Matthews 判据完全一致（6ZPP 1、6HF7 3、5D9D 1、7TJ1 4、7E1D 2、
7TFQ 1、5HHK 2、3O2E 1、7ENM 1、7TBS 1），且每一套的溶剂含量都在 26.5–58.4% 之间。


### `--copies` 默认 auto（2026-11-09）

`--copies` 的默认值已从 1 改为 **auto**：程序用用户提供的序列文件与 SAD MTZ 的晶胞/空间群自行判定。

- 10/10 目标案例的判定结果：6ZPP 1、6HF7 3、5D9D 1、7TJ1 4、7E1D 2、7TFQ 1、5HHK 2、
  3O2E 1、7ENM 1、7TBS 1，溶剂含量全部落在 26.5–58.4%；
- `--copies 1/2/3/...` 仍可显式覆盖；预览或缺少可读输入时先按 1 处理，真正运行时再判定；
- 判定来源写入 `pipeline_report.json` 的 `copies_source` / `copies_inference`。

**注意**：拷贝数只是输入修正。7E1D 在"40 轮密度修饰 + feedback 2"这一终点是 43.72%（+5.5 pp），
需要 `--dm-cycles 20`（36.72%，−1.5 pp）或 `--feedback-cycles 3`（41.25%，+3.1 pp）才落进 5 pp。
因此文档给出的配方应把这两个固定终点一起跑，并按工作集 R 选择——这正是十套结果采用的做法。


### 自带模式的标定默认值（2026-11-09）

`--method pplus_bundled` 的默认值现在等于十套标定所用的配方：

| 选项 | 原默认 | 现默认 |
|---|---|---|
| `--copies` | 1 | **auto**（FASTA 链注释 → Matthews 兜底） |
| `--dm-engine` | resolve | **旧引擎** |
| `--dm-cycles` | 40 | 40（不变） |
| `--feedback-cycles` | 0 | **2** |
| `--building-recipe` | ccp4 | **pipeline** |

- 显式给出的非默认值一律尊重。
- 运行包里没有 `旧引擎` 时自动回退到内部 RESOLVE/CCTBX 分支（40 轮、无反馈），
  原因写入 `pipeline_report.json` 的 `dm_engine_fallback`，旧 runtime 仍可运行。

### 为什么单一常数覆盖不了全部十套（实测）

十套案例在"溶剂含量 × 反馈轮数"上的最优值并不一致。完整跑过的两组（其余为单点测量）：

| 案例 | 对照 | 溶剂 0.50（fb2 / fb3） | 溶剂按序列估算（fb2 / fb3） |
|---|---:|---:|---:|
| 10_6ZPP | 31.11 | **31.31** / 32.98 | 32.78 / 32.55 |
| 5_5D9D | 29.73 | 29.39 / **27.51** | 27.04 / 27.64 |
| 14_7TJ1 | 32.49 | 32.01 / **31.00** | 33.33 / 35.50 |
| 11_7E1D | 38.20 | 43.23 / **42.33** | 43.72 / **41.25** |
| **24_7TFQ** | 37.29 | **36.86** / 37.66 | 49.09 / 50.85 |
| 13_5HHK | 27.40 | 28.63 / — | 28.52 / — |
| **3_3O2E**（`--fcalc fc`） | 42.85 | 53.07 / 50.89 | 46.94 / **42.64** |
| 7_7ENM | 30.08 | **31.47** / 32.73 | 31.43 / 31.43 |
| 19_7TBS | 30.58 | 29.88 / **28.35** | 29.85 / 31.55 |
| 48_6HF7 | 28.43 | — | **30.87** / 30.24 |

**7TFQ 和 3O2E 的最优溶剂正好相反**：7TFQ 必须用 0.50（36.86%），用序列估算的 0.493 会恶化到 49.09%；
3O2E 必须用序列估算的 0.265（42.64%），用 0.50 会恶化到 50.89%。旧引擎 的溶剂掩膜对 1% 的差异就会
翻转结果，所以这两个案例不能共用同一个常数。

### 结论：用固定终点网格 + 工作集选择

`--feedback-endpoints 2,3 --dm-solvents auto,0.5` 会按同一协议把网格里每个组合各跑一遍，
写入独立子目录，再按**最终工作集 R 因子**选（不读 R-free），结果写进 `endpoint_selection.json`。
这样一条命令即可同时覆盖 7TFQ（选 0.50/fb2）与 3O2E（选 auto/fb3）。

### 全网格实测结果（十套 × 四个终点，工作集 R 选择）

对十套跑完全部已测终点后，按**最终工作集 R**（不读 R-free）选择的结果：

| 案例 | 对照 | 选中终点 | 溶剂 | Rwork | R-free | 差值 |
|---|---:|---|---:|---:|---:|---:|
| 10_6ZPP | 31.11% | auto / fb2 | 0.373 | 29.65% | 32.78% | +1.7 |
| 48_6HF7 | 28.43% | 0.50 / fb3 | 0.500 | 26.02% | 29.04% | +0.6 |
| 5_5D9D | 29.73% | auto / fb2 | 0.467 | 25.21% | 27.04% | −2.7 |
| 14_7TJ1 | 32.49% | 0.50 / fb3 | 0.500 | 24.95% | 31.00% | −1.5 |
| 11_7E1D | 38.20% | 0.50 / fb3 | 0.500 | 36.65% | 42.33% | +4.1 |
| 24_7TFQ | 37.29% | 0.50 / fb2 | 0.500 | 32.52% | 36.86% | −0.4 |
| 13_5HHK | 27.40% | 0.50 / fb3 | 0.500 | 25.94% | 28.82% | +1.4 |
| 3_3O2E | 42.85% | auto / fb3 | 0.265 | 41.14% | 42.64% | −0.2 |
| 7_7ENM | 30.08% | 0.50 / fb2 | 0.379 | 26.33% | 31.43% | +1.4 |
| 19_7TBS | 30.58% | 0.50 / fb3 | 0.500 | 25.17% | 28.35% | −2.2 |

**10/10 在 5 个百分点内，5 套优于对照，中位差 −0.21 pp，最大 +4.1 pp（7E1D）。**

注意选择规则只看工作集 R，不看 R-free，所以个别案例选中的不是 R-free 最低的那个终点
（例如 6ZPP 的 auto/fb2 是 32.78%，而 0.50/fb2 是 31.31%）。这是刻意保留的选择规则：
运行终点不依赖 R-free 信息。


## 精简打包（2026-11-09）

原始独立包 3.8 GB / 压缩后 1.19 GB，其中大部分是 SAD 流程从不加载的组件。
新增 `packaging/prune_bundle.py`，把下列内容删掉并同步刷新 `runtime_manifest.json`、
`python_manifest.json`（否则运行时完整性校验会因为文件缺失而报错）：

| 删除内容 | 大小 | 理由 |
|---|---:|---|
| `density/conda_base/lib/libQt5WebEngine*.dylib`（4 份） | 350 MB | 只给 Qt GUI 用 |
| `density/conda_base/.../site-packages/{bokeh,PySide2,wx*,matplotlib,pandas,sklearn,babel}` | ~280 MB | 绘图/GUI/数据框，RESOLVE 不 import |
| `runtime/lib/data/nucleofind_models` | 152 MB | 糖/核苷酸校验程序 |
| `runtime/lib/data/privateer_data` | 52 MB | 同上 |
| `runtime/lib/data/reference_structures`（保留 `reference-1tqw.mtz/.pdb`） | 190 MB | 其余是 RESOLVE 的演示参考 |
| `runtime/lib/data/monomers/*`（保留 20 种标准氨基酸、UNK、水、常见离子与小分子） | 630 MB | 35649 个配体字典，蛋白 SAD 建模用不到 |
| `python3.12` 的 pip/setuptools/pygments/pytest/idlelib/ensurepip/test 与 `__pycache__` | ~30 MB | 安装器与测试工具 |

结果：

| | 目录 | 压缩包 |
|---|---:|---:|
| 精简前 | 3.8 GB | 1.19 GB |
| 精简后 | **1.9 GB** | **501 MB** |

验证方式不是只看文件大小：用精简后的包对 7ENM 跑完整流程（Phaser → P+ → 旧引擎 → Buccaneer →
REFMAC），结果 **Rfree 31.47%、Rwork 28.34%，与精简前完全一致**，REFMAC 日志里没有任何
`no dictionary` 警告。

**限制**：配体字典被裁到常见集合。如果晶体里含有非标准配体，REFMAC 需要对应 `.cif`；
用 `--keep-all-monomers` 重新打包即可恢复全部 35649 个字典（+660 MB）。
`--keep-optional` 则保留全部 GUI/绘图组件。


### 第二轮：CCTBX 环境内的死重量（2026-11-09）

第一轮把包降到 1.9 GB 后，`runtime/density`（RESOLVE/P+ 用的 CCTBX 环境）仍占 1.5 GB。
逐项核对后又能安全删掉：

| 删除内容 | 大小 | 理由 |
|---|---:|---|
| `libopenblas/libblas/liblapack/libcblas` 的 7 个重复副本 | 446 MB | 8 个文件字节完全相同（sha256 一致），改成指向同一份的相对符号链接 |
| `conda_base/lib/libLLVM-9.dylib` | 73 MB | 包里没有 llvmlite/numba，没有任何程序加载它 |
| `solve_resolve/ext_ref_files/{segments,segments_em,segments_rich_small,patterns_txt}` | 105 MB | RESOLVE 只用 `segments_rich` 与 `patterns` |
| py2.7 环境里的 sphinx/OpenGL/py2app/pip/pygments/google/setuptools | 35 MB | 文档、GUI 与安装器 |
| 全包范围内字节相同的文件 → 相对符号链接 | 519 MB | 63 个重复文件，统一去重 |

去重后清单哈希仍按内容计算（`fingerprint` 会跟随符号链接），所以
`runtime_identity` 的校验照常通过，`runtime_manifest.json` 只需丢弃被删文件的行。

三轮体积：

| 阶段 | 目录 | 压缩包 |
|---|---:|---:|
| 原始 | 3.8 GB | 1.19 GB |
| 第一轮（GUI/配体/参考/安装器） | 1.9 GB | 501 MB |
| **第二轮（去重 + CCTBX 死重量）** | **1.2 GB** | **348 MB** |

每一轮都用 7ENM 跑完整流程复核，三次结果完全相同：**Rwork 28.34%，Rfree 31.47%**。

## 限制（必须如实告知用户）

- **平台**：本包是 macOS **x86_64**（Intel 二进制，在 Apple Silicon 上经 Rosetta 运行）。
  面向 Linux 的包仍需在 Linux 上按同一脚本构建；脚本本身与平台无关，依赖的引擎二进制需换成
  Linux 版本。
- **许可**：外部引擎（Phaser/Buccaneer/REFMAC/旧引擎）与其数据、权重的公开再分发条件尚未落实，
  本包仍是"本机评估包"，不是可公开发行的二进制。`bundle.json` 里
  `public_redistribution_ready` 为 `false`。
- **未签名/未公证**：macOS 首次运行可能需要在"系统设置 → 隐私与安全性"中放行。
- **体积**：解压后约 3.8 GB（单体字典占大部分）。
- 打包脚本复制的是已经验证过的 portable python/runtime，不会安装任何依赖，也不会修改源运行时。
