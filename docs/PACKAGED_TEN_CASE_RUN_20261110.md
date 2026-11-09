# 打包产物十案例实测（2026-11-10）

## 1. 测的是什么

直接解压发布用的压缩包，用包内启动器跑完十个困难案例，检验"压缩后的程序"本身能否复现开发树的结果。

| 项目 | 值 |
| --- | --- |
| 产物 | `dist/OASIS4.3-standalone-20261108.tar.gz` |
| SHA-256 | `485799f2f22f98b57294ea328874476d5ce2906b54b7b98bf9ddc3d9c36f44c4` |
| 大小 / 条目 | 347,770,231 字节 / 28,707 |
| 解压位置 | `work/exp/pkg10/extract/OASIS4.3-standalone-20261108/` |
| 启动方式 | `./OASIS4.3.command run_sad_to_model.py …` |
| 运行环境 | `env -i PATH=/usr/bin:/bin HOME=/tmp TMPDIR=/tmp LANG=C LC_ALL=C` |
| 代码版本 | `OASIS4.3-Build` 分支 `focus/oasis43-two-methods-20261023`，`b2fe55d` + 本文件的提交 |
| 运行记录 | `work/exp/pkg10/batch.log`、各案例 `*.out`、每案例 `endpoint_selection.json` |
| 结果汇总 | `work/PKG10_RESULTS.json` |

每个案例一个目录，内含 `01_inputs`、`02_sad`、`03_dm`、`04_bridge`、`05_build`、`06_finish` 与各阶段日志，可逐轮复查。

## 2. 自包含性证据

- 运行环境用 `env -i` 清空，`PATH` 只有 `/usr/bin:/bin`，不含任何 CCP4、Phenix、conda 路径。
- 启动器把 `CCP4`、`CLIBD`、`CLIBD_MON`、`SYMINFO`、`CINCL` 全部指向包内 `runtime/`。
- 每次运行的 `pipeline_report.json` 记录 `runtime.root = <包>/runtime`，即 Phaser / Buccaneer / REFMAC / cparrot 均来自包内。
- 包内符号链接 1,194 个，全部相对链接、无断链；`runtime/bin` 含 `phaser`、`refmac5`、`cparrot`。
- 启动器现导出 `PHASEBUILD_STANDALONE=1`，包内缺文件会直接报错，不会悄悄回退到本机 CCP4。

## 3. 结果

参考值为 SHELXE/Buccaneer/REFMAC 对照流程的独立 R-free；P+ 侧只用工作集 R 在预声明端点里选终点，从不读 R-free。

| 案例 | P+ Rwork | P+ R-free | 对照 R-free | 差值 | 端点数 |
| --- | --- | --- | --- | --- | --- |
| 6ZPP | 30.04 | **31.31** | 31.11 | +0.20 | 2 |
| 6HF7 | 26.02 | **29.04** | 28.43 | +0.61 | 2 |
| 5D9D | 25.43 | **27.51** | 29.73 | −2.22 | 2 |
| 7TJ1 | 23.92 | **31.08** | 32.49 | −1.41 | 2 |
| 7E1D | 36.65 | **42.33** | 38.20 | +4.13 | 2 |
| 7TFQ | 32.52 | **36.86** | 37.29 | −0.43 | 2 |
| 5HHK | 25.94 | **28.82** | 27.40 | +1.42 | 2 |
| 3O2E | 41.14 | **42.64** | 42.85 | −0.21 | 4 |
| 7ENM | 28.34 | **31.47** | 30.08 | +1.39 | 2 |
| 7TBS | 25.17 | **28.35** | 30.58 | −2.23 | 2 |

- 10/10 落在对照 ±5 个百分点内；5 个优于对照，中位差 +0.20 pp，最大差 +4.13 pp（7E1D）。
- 与开发树 `work/exp/batch7` 的同协议结果逐例一致（6ZPP 31.31、7E1D 42.33、7TJ1 31.08、7TBS 28.35、7TFQ 36.86、7ENM 31.47、5HHK 28.82、6HF7 29.04、5D9D 27.51），差异只在 0.0–0.1 pp，属端点选择并列时的抖动。
- 3O2E 用的是它的记录在案的专用参数：`--dm-solvents auto,0.5 --fcalc fc`（4 个端点）。若只给固定 0.5 溶剂，它停在 53.07，这是案例参数问题而非打包问题。

## 4. 耗时

00:14 同时启动十个案例（20 核机器），完成时刻：

| 案例 | 完成 | 案例 | 完成 |
| --- | --- | --- | --- |
| 3O2E | 00:25 | 5D9D | 00:47 |
| 7E1D | 00:29 | 5HHK | 01:00 |
| 7ENM | 00:31 | 6HF7 | 01:31 |
| 7TFQ | 00:37 | | |
| 6ZPP | 00:39 | | |
| 7TBS | 00:43 | | |
| 7TJ1 | 00:45 | | |

全批 77 分钟结束；最慢的 6HF7（4 个 Tb + 2 个 Mg，拷贝数 3）单案例 77 分钟，其余 11–46 分钟。

## 5. 拷贝数

十个案例全部由 `--copies auto` 自动定出，来源均为 `matthews_plausible_range`：6ZPP 1、6HF7 3、5D9D 1、7TJ1 3、7E1D 2、7TFQ 1、5HHK 2、3O2E 1、7ENM 1、7TBS 1。7TJ1 开发树按 FASTA 头标推断为 4，此处按 Matthews 得 3，R-free 31.08 与开发树的 31.00 相当，说明该案例对拷贝数不敏感。

## 6. 本次修的两处打包问题

1. 打包版 GUI 的默认溶剂分数仍是 `auto`，用界面直接开跑 7ENM 只到 38.73%，与文档协议的 31.5% 不符。已把 `dm_solvent` 默认改为标定常数 0.5。
2. 启动器没有导出 `PHASEBUILD_STANDALONE=1`，包内若缺某个二进制会静默回退到本机 `/Applications/ccp4-9`，掩盖打包不完整。已在启动器中固定该变量与 `PHASEBUILD_RUNTIME`。

两项修改后重新打包，`app/` 与上一版逐文件比较只差 `oasis_workflow.py` 一行默认值，并清除了残留的 `__pycache__`。

## 7. 边界与未决

- 包是 macOS x86_64、未签名，`public_redistribution_ready=false`，外部引擎（Buccaneer/REFMAC/Phaser/cparrot）的再分发授权仍未解决。
- 6HF7 用的是按预期占位度手工标注的 Tb/Mg 位点文件，属假设性输入修复，需实验确认。
- 3O2E 与 7TFQ 依赖案例专用参数（前者 `auto` 溶剂网格 + `--fcalc fc`，后者 `pplus_trust` + 观测先验），通用默认值在这两个案例上会明显变差。
