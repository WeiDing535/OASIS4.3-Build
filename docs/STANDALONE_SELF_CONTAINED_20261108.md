# 完全独立安装包（不再依赖外部 CCP4 / Phenix）

## 结论

P+ 升级引入的三处外部依赖已全部消除，程序现在可以打成**一个目录 + 一个 tar.gz**，
用户不需要另装 CCP4、Phenix、Python、NumPy、SciPy 或 Gemmi。

| 环节 | 原依赖 | 现在 |
|---|---|---|
| 密度修饰 | 外部 CCP4-9 的 `cparrot` | 自带 `runtime/bin/cparrot`（只依赖 runtime 内已有的 clipper/ccp4c/mmdb2 库与 `lib/data`） |
| 构建循环 | CCP4 的 `buccaneer_pipeline`（需要 CCP4 的 python3.9 + `CCP4pipeline.py`） | `src/phasebuild/ccp4_pipeline.py` 用自带 `cbuccaneer` + `refmac5` 复现同样的循环 |
| 反馈用 Fcalc | `phenix.fmodel` | 自带 `refmac5` 的 0 循环计算 `FC_ALL/PHIC_ALL`（不移动坐标、含批量溶剂） |

## 关键验证

**Parrot 可直接随包运行。** `cparrot` 的加载路径全部是 `@loader_path/../lib`，而
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
python3 OASIS4.3-Build/packaging/package_standalone.py \
  --app-source   OASIS4.3-Build \
  --python-source "OASIS4.3-dev17-portable local/python" \
  --runtime-source "OASIS4.3-dev17-portable local/runtime" \
  --parrot-binary /Applications/ccp4-9/bin/cparrot \
  --output dist/OASIS4.3-standalone-YYYYMMDD \
  --tar dist/OASIS4.3-standalone-YYYYMMDD.tar.gz
```

脚本会：
1. 校验 python/runtime 清单的每个文件哈希；
2. 复制 `python/`、`runtime/`，把 `cparrot` 装进 `runtime/bin` 并登记到
   `runtime_manifest.json`（同时剔除自引用的 `runtime_manifest.sha256` 条目并重算校验文件）；
3. 用 `package_oasis.py` 的聚焦收集器复制 `app/`（含 `run_trace_feedback.py`、
   `run_build_pipeline.py`）；
4. 写 `OASIS4.3.command` 启动器与 `使用说明.txt`；
5. 在干净环境下自检：Python 包导入、四个引擎存在且可执行、`cparrot` 能起、CLI 能起、
   `parrot_dm.locate()` 解析到**包内**的 cparrot；
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

## 限制（必须如实告知用户）

- **平台**：本包是 macOS **x86_64**（Intel 二进制，在 Apple Silicon 上经 Rosetta 运行）。
  面向 Linux 的包仍需在 Linux 上按同一脚本构建；脚本本身与平台无关，依赖的引擎二进制需换成
  Linux 版本。
- **许可**：外部引擎（Phaser/Buccaneer/REFMAC/Parrot）与其数据、权重的公开再分发条件尚未落实，
  本包仍是"本机评估包"，不是可公开发行的二进制。`bundle.json` 里
  `public_redistribution_ready` 为 `false`。
- **未签名/未公证**：macOS 首次运行可能需要在"系统设置 → 隐私与安全性"中放行。
- **体积**：解压后约 3.8 GB（单体字典占大部分）。
- 打包脚本复制的是已经验证过的 portable python/runtime，不会安装任何依赖，也不会修改源运行时。
