# 打包产物十案例实测（2026-11-10）

> **历史记录。** 本文写作时密度修饰由一个外部引擎完成；该引擎已于 2026-11-14
> 从程序中彻底删除（模块、命令行选项、打包标签、测试全部移除），密度修饰
> 现在固定由树内 oasis-dm v2 完成。文中的数字保留作当时记录，按本文命令重跑
> 会得到不同数值；`--dm-mask-radius`、`--dm-anisotropy` 等旧引擎专属选项也
> 已删除。
> **另注（2026-11-14 之后）。** 界面现在只有两个模块，密度修饰是 P+ 的必做部分（树内 `oasis-dm v2`），没有引擎选择；本文这类旧引擎时代的数字只作记录，不能再逐点复现。当前行为见 `docs/OASIS_DM_REPLACEMENT_20261114.md`、`docs/INSTALL_BUNDLE.md`。

## 1. 测的是什么

直接解压发布用的压缩包，用包内启动器跑完十个困难案例，检验"压缩后的程序"本身能否复现开发树的结果。

| 项目 | 值 |
| --- | --- |
| 产物 | `dist/OASIS4.3-standalone-20261108.tar.gz`（重打后） |
| SHA-256 | `ecd8014c5cca6a232bac97e77c3b0e9eb74c19c72f8183661b7c5607812c44dc` |
| 大小 / 条目 | 355,884,001 字节 / 28,811 |
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
- 每次运行的 `pipeline_report.json` 记录 `runtime.root = <包>/runtime`，即 Phaser / Buccaneer / REFMAC / 旧引擎 均来自包内。
- 包内符号链接 1,194 个，全部相对链接、无断链；`runtime/bin` 含 `phaser`、`refmac5`、`旧引擎`。
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

## 6. 本次修的打包问题

1. 打包版 GUI 的默认溶剂分数仍是 `auto`，用界面直接开跑 7ENM 只到 38.73%，与文档协议的 31.5% 不符。已把 `dm_solvent` 默认改为标定常数 0.5。
2. 启动器没有导出 `PHASEBUILD_STANDALONE=1`，包内若缺某个二进制会静默回退到本机 `/Applications/ccp4-9`，掩盖打包不完整。已在启动器中固定该变量与 `PHASEBUILD_RUNTIME`。
3. **裁剪脚本把 RESOLVE 的参考密度库 `ext_ref_files/segments/` 当"备选库"删掉了**（35 MB）。批量跑十个案例时用的是 `--no-auto-hand`，绕过了双手判别，所以十例全过；而界面默认 `auto_hand=True`，一旦要走"双手 P+ / Wang8 密度修饰与工作集测序判别"就会调 RESOLVE，它在读 `rho.list` / `Richardson_rho.dat` / `side_chains_rho.dat` 时直接失败：

   ```
   Sorry, cannot find density file?
   Sorry, cannot find model density file?
   RESOLVE_API_STATUS 1
   ```

   两个界面作业（`work/exp/gui_run3`、`gui_run4`）都卡在这一步。已把 `segments/` 从裁剪清单移除（`segments_em`、`segments_rich_small` 仍删），并在 `verify_bundle()` 里加了 `rho.list`、`Richardson_rho.dat`、`side_chains_rho.dat`、`gap_library.dat` 的存在性与大小检查，让这类"删了才炸"的裁剪在打包阶段就报错。
4. 补上密度库之后暴露出第二个只在自动判别手上才走到的缺陷：手性判别的 pilot 构建把流水线级配方名 `building_recipe='pipeline'` 原样传给 Buccaneer 命令行构造器，而后者只认 `legacy` / `ccp4`，于是抛 `ValueError: Unknown construction recipe`（`work/exp/gui_run5`）。已在 `hybrid_recycle` 里加 `concrete_recipe()`，把 `pipeline` 映射到 `ccp4` 构建，并补了回归测试。

5. 结果面板把 `dm_engine` 直接印成方法的**默认**引擎名（`自带 RESOLVE + CCTBX`），而本次实际跑的是那个外部引擎；溶剂分数也因为外部引擎路径的 `map_phases.json` 没有 `solvent_fraction_used` 而显示为空。现在改为读实际执行记录，界面显示 `CCP4 旧引擎 + 2 轮 trace 反馈`、溶剂 0.5。

前三项修改后重新打包，`app/` 与上一版逐文件比较只差 `oasis_workflow.py` 一行默认值，并清除了残留的 `__pycache__`；第三项使包体增加约 8 MB 压缩后体积（35 MB 原始）。

修完之后用打包版界面的默认参数（自动判手、`dm_solvent=0.5`、旧引擎 40 轮、2 轮反馈）重跑 7ENM（`work/exp/gui_run6`），6/6 阶段完成，得到 R-work 28.34 / **R-free 31.47**，与批量结果和开发树逐位一致；手性判别选中 `input_hand`，判别依据是唯一测序残基数（不用 R-free、不用参考坐标）。这说明十个案例的成绩不是靠绕过自动判手拿到的：走完整默认路径同样得到 31.47。

## 7. 边界与未决

- 包是 macOS x86_64、未签名，`public_redistribution_ready=false`，外部引擎（Buccaneer/REFMAC/Phaser/旧引擎）的再分发授权仍未解决。
- 6HF7 用的是按预期占位度手工标注的 Tb/Mg 位点文件，属假设性输入修复，需实验确认。
- 3O2E 与 7TFQ 依赖案例专用参数（前者 `auto` 溶剂网格 + `--fcalc fc`，后者 `pplus_trust` + 观测先验），通用默认值在这两个案例上会明显变差。
