# 精简版：只做模块 1（P+ 定相 + oasis-dm 密度修饰）

> **2026-11-14 起界面只有两个模块，且旧密度修饰引擎已从程序中删除。** 密度修饰
> 不再是可勾选的独立模块，而是 P+ 内置的 `oasis-dm v2`，属于模块 1、必做、没有
> 开关也没有引擎下拉框；模块 2 现在是 Buccaneer/REFMAC。本文保留 11-12 那天的
> 模块编号、旧引擎数字与检测顺序作为历史记录（文中那条旧复现命令已不能运行），凡与当前界面冲突处以本节说明为准；当前行为见
> `docs/INSTALL_BUNDLE.md`、`docs/OASIS_DM_REPLACEMENT_20261114.md`。

这一页回答两个问题：**十案例是不是模块 1 单独跑出来的**，以及**精简版（只留模块 1）能交付什么**。界面怎么排、CCP4 怎么检测，也写在同一页。

## 1. 先更正一个前提

十案例 **不是** 第一个模块单独完成的。十案例用的协议是

```
--method pplus_bundled --dm-cycles 40 \
--feedback-cycles 2/3 --trust-wang-radius 8 --dm-solvent 0.5 \
--no-auto-hand --input-columns
```

即 **模块 1（P+ 定相）+ 模块 2（旧外部引擎密度修饰，已删除）+ 模块 3（Buccaneer/REFMAC）** 三段串起来的结果。三个环节缺一不可：

| 只跑到 | 手上有什么 | 3O2E 实测 | 十案例均值 |
| --- | --- | --- | --- |
| 模块 1（P+） | 相位 PHIB、FOM、PPLUS 的 MTZ | complex CC **0.123**，mean FOM 0.482 | **0.272**（另一个 P+ 后端 0.279） |
| 模块 1+2（旧外部引擎） | 密度修饰后的图与相位 | 0.549 | **0.688** |
| 模块 1+2+3 | 模型、R/R-free、逐轮精修记录 | R/R-free 41.14 / 42.64 | 10/10 落在对照 ±5 pp 内 |

"complex CC" 是对 campaign 精修模型的平移不变复空间图相关（`work/noCCP4/compare_map_phases.py`），不是 omit-map，也不是振幅相关。

所以：**模块 1 给的是相位，不是模型**。没有模型就没有 R-free，也没有"十案例 10/10"这个说法——那个数字属于三段链。模块 1 单独的价值是把可用的初始相位做出来，让模块 2 有一个够好的起点；按图相关算，起点离成品还差 0.4 左右。

## 2. 精简版是什么

界面把流程拆成两个模块（11-12 那天是三个，见文首说明）：

| 模块 | 程序 | 依赖 | 可关闭 |
| --- | --- | --- | --- |
| ① P+ 定相 + oasis-dm 密度修饰 | 包内 Phaser 前端 + P+ 信任域更新 + 树内 `oasis-dm v2` | 只要包内运行时 | 否（必做） |
| ② Buccaneer 建模 + REFMAC 精修 | `cbuccaneer` + `refmac5` | 本机 CCP4 | 是 |

- 关掉 ②：流程在 `03_dm` 结束，输出 `02_sad/experimental.mtz`（`H K L PHIB FOM PPLUS`）、`observations.mtz` 和密度修饰后的 `03_dm/map_phases.mtz`（`FWT/PHWT/HL*`），报告里 `completed` 到 `dm` 为止。密度修饰本身关不掉。
- 精简版不做自动判手（`auto_hand` 需要试建），手性由重原子亚结构决定；界面上的"自动判手"随 ② 一起关闭。
- 11-12 那天的「关掉 ② = 停在 `02_sad`」现在只能用命令行 `--no-dm-enabled` 做到，属于诊断手段（隔离 P+ 自身的相位质量），界面不提供。

命令行写法：

```
python run_sad_to_model.py --no-build-model ...  # 界面上的「只跑模块 1」：P+ + oasis-dm
python run_sad_to_model.py --lite ...            # 诊断：连密度修饰也跳过，只到 P+ 定相
python run_sad_to_model.py --no-dm-enabled ...   # 同 --lite 的密度修饰部分，显式写出来
```

`--lite` 会同时把 `dm_enabled`、`build_model`、`auto_hand` 置为关闭，并把 `feedback_cycles` 归零——trace 反馈本来就是"用部分模型算 Fcalc 回喂密度修饰"，没有建模就没有反馈。

## 3. 界面：参数分组与模块开关

表单自上而下：

1. **任务名称**
2. **运行模块（按本机 CCP4 检测结果默认勾选）**：模块 1（P+ 定相 + oasis-dm 密度修饰）固定勾选并置灰；模块 2（Buccaneer/REFMAC）一个复选框，下面依次是「CCP4 位置（可选）」——留空用包内运行时，装在其它的地方时填安装根目录或 `ccp4.setup` / bashrc 的路径——、这一栏的说明、「重新检测 CCP4」按钮和一行检测结论。
3. **输入文件**：SAD MTZ、重原子 PDB、序列、元素、波长、份数、新输出目录。
4. **① P+ 定相参数**：P+ 后端单选、P+ 参数方案、定相迭代周期、每反射 KL 上限、每个目标的关系上限。
5. **② oasis-dm 密度修饰参数（P+ 的一部分，必做）**：一段固定说明（引擎 oasis-dm v2、不可关闭、不需要 `hist.lib`）、溶剂比例、密度修饰周期、trace 反馈轮数。没有引擎下拉框，也没有「掩膜过滤半径」（那是已删除的外部引擎的参数）。
6. **③ Buccaneer 建模 + REFMAC 精修参数**：自动判手、Buccaneer/REFMAC 轮数、每轮 REFMAC 周期、构建配方、精修相位来源。
7. **显示高级参数**（默认收起）：输入细节、P+ 掩膜与采样、运行资源与终点组合。

关闭模块 2 时，它下面的参数框整组置灰——不是隐藏，让人看得见"关了哪些"。模块 1 的密度修饰参数永远不会置灰，因为那一段没有开关。「CCP4 位置」不属于任何模块，始终可编辑，因为检测和运行都要读它。

文件、目录类字段一律排成两行——标题一行，输入框和「浏览…」一行（`_field()`）。早期版本把输入框 `side='left', expand=True` 排在按钮前面，窗口一窄按钮就被挤出可视区；现在按钮先占位、输入框吃掉剩余宽度，1180×1050 下实测每行 361 px：输入框 263 px、按钮 94 px。`tests/test_desktop_window.py::test_browse_buttons_stay_visible_and_ccp4_location_feeds_detection` 会真的建窗口量像素，防止再退化。

## 4. CCP4 检测

`src/phasebuild/module_detection.py` 只回答一个问题：**这台机器上有没有能跑的可执行文件**。找不到是正常答案，不抛异常。

按顺序在这些目录里找（命中即止）：

1. 表单里「CCP4 位置（可选）」填的值：目录 → 它自己和它的 `bin/`；`ccp4.setup` / bashrc → 解析里面的 `CCP4=`、`CCP4_ROOT=`、`CBIN=`、`CCP4_BIN=` 赋值（展开 `$HOME`），得到安装目录；
2. 包内运行时的 `bin/`（免安装包默认走这条）；
3. 环境变量 `CBIN`、`CCP4_BIN`、`OASIS_CCP4_BIN`；
4. 环境变量 `CCP4`、`CCP4_ROOT` 及其 `bin/`；
5. `/Applications/ccp4-*/bin`、`/opt/ccp4*/bin`、`/usr/local/ccp4*/bin`（双击启动的 GUI 几乎继承不到环境变量，所以显式找安装位置）；
6. `PATH` 的每一项。

判定规则（11-12 那天是前两行，现在是第二行）：

| 模块 | 判据 |
| --- | --- |
| ① oasis-dm 密度修饰 | 不需要检测：树内 `oasis-dm v2`，随程序提供，永远可用 |
| ② Buccaneer/REFMAC | `cbuccaneer` 或 `buccaneer` **且** `refmac5` 都可执行 |

结论直接写在那行提示里：

- 找到：模块 2 默认勾选，提示"已检测到 CCP4：模块 2（建模与精修）默认勾选，可随时取消。"
- 没找到：模块 2 关闭并置灰，提示"未检测到 CCP4：只做模块 1（P+ 定相 + oasis-dm 密度修饰），输出相位/FOM/MTZ。"，第二行说明装好 CCP4 或设好 `CCP4`/`CBIN` 后点「重新检测」。
- 外部密度修饰程序不参与检测：程序内已没有它，检测只看 `cbuccaneer`/`refmac5`。

重新用"沿用参数"载入一条旧记录时，记录里的开关优先；但如果记录要求一个本机现在没有的模块，复选框会被压回关闭状态——不会留一个点了必然失败的勾。

检测结果只影响**默认值**，不影响可执行性判断：用户仍可手动勾选模块 2，真正运行前的校验在 `sad_pipeline.validate()`，那里会给出具体缺哪个程序。密度修饰没有可勾的东西，也就没有"检测不到"这种状态。

### 4.1 装在别处怎么办

「CCP4 位置」这一栏同时是**运行目录**：子进程一律按 `<root>/bin/<exe>` 启动，运行环境还会清掉继承来的 `CCP4*`（`runtime_environment()`）。所以"检测到了"必须等价于"跑得起来"，两者共用同一个值：

- 填目录：`/opt/ccp4-9` 和 `/opt/ccp4-9/bin` 都认，统一归一到根目录；
- 填脚本：`ccp4.setup`、`.bashrc` 都行，解析出 `CCP4=` 后把**根目录写回这一栏**，定相、密度修饰、建模都改用它；
- 解析不出来（或填的目录里没有这些程序）：提示"`<文件名>` 里没有找到 CCP4 / CBIN 赋值"或"`<目录>` 里没找到 CCP4 程序"，并把「CCP4 位置」改回默认运行时——否则这一栏会带着一个跑不起来的目录进到运行里。提示里同时写出实际用的是哪个安装。

CCP4 自带的 `phaser` 就能当 P+ 的前端，所以一个完整的 CCP4 安装足以跑完整条链。实测（3O2E，本机 `/Applications/ccp4-9`）：

- 精简版 `--runtime /Applications/ccp4-9`：5763 反射、mean FOM 0.482，与包内运行时逐位一致；
- 完整链同一 `--runtime`：P+ → oasis-dm 密度修饰 → 2 轮 trace 反馈 → Buccaneer/REFMAC 全部从 `/Applications/ccp4-9/bin` 启动并跑通（见 4.2）。

代价是要求这个安装"什么都有"：若填的外部 CCP4 缺 `bin/phaser`，P+ 会直接报"P+ 引擎缺失：…/bin/phaser"，并提示清空该栏、改回包内运行时。

### 4.2 两种运行时的对照（3O2E）

同一台机器、同一份输入，只换 `--runtime`：

| 运行 | `--runtime` | 结果 |
|---|---|---|
| 精简版 | 包内运行时 | 5763 反射、mean FOM 0.482，`02_sad/experimental.mtz` sha256 `87a1be53…` |
| 精简版 | `/Applications/ccp4-9` | 5763 反射、mean FOM 0.482，MTZ 逐字节相同 |
| 完整链（feedback 2、溶剂 0.5） | 包内运行时（十案例记录 `endpoint_fb2_s0.50`） | R/R-free 46.97/53.07 |
| 完整链（同一协议） | `/Applications/ccp4-9` | R/R-free 46.97/53.07；`02_sad`、`03_dm` 的 MTZ 与左行逐字节相同 |

换成外部 CCP4 只改了 `phaser`、`cbuccaneer`、`refmac5` 的路径（`pipeline_report.json` 的 `config.dm_executable` 与 `runtime.root` 可见；表里这两行是旧外部引擎时代的记录，现在同样的对照只换 `phaser`/`cbuccaneer`/`refmac5`），相位、密度修饰图、模型全部一致；参考结构也直接取 `ccp4-9/lib/data/reference_structures/reference-1tqw.*`，不用另外准备。

顺带一条容易看错的：十案例表里 3O2E 的 41.14/42.64 是 feedback 3 那个终点（`endpoint_fb3`），feedback 2 是 46.97/53.07。两个数在同一台机器上都能复现，差别来自反馈轮数，不是运行时。

## 5. 实测记录（3O2E，精简版）

```
cd OASIS4.3-Build
PYTHONPATH=src .conda-unified/bin/python run_sad_to_model.py --lite \
  --sad .../01_inputs/selected_sad.mtz --native .../native.mtz \
  --sequence .../building.fasta --heavy .../02_sad/refined_heavy.pdb \
  --elements K --wavelength 1.5418 --input-columns \
  --runtime ".../OASIS4.3-dev17-portable local/runtime" \
  --output /tmp/oasis_lite_3o2e
```

结果（`/tmp/oasis_lite_3o2e/pipeline_report.json`）：

- `status=completed`，`current_name="P+ 定相完成（未启用密度修饰与建模）"`；
- `completed` 只有 `prepare`、`sad`；`dm_phases=null`、`primary_model=null`；
- 配置里 `dm_enabled=false`、`build_model=false`、`feedback_cycles=0`；
- 产物 `02_sad/experimental.mtz`：5763 个反射、分辨率 1.95 Å、列 `H K L PHIB FOM PPLUS`、mean FOM 0.482；
- 与十案例同一条路径的 `02_sad/experimental.mtz` 逐位相同（图相关同为 0.123），说明精简版没有偷偷改变模块 1 的行为。

耗时：该案例约 2 分钟（那是 `--lite`：只做 Phaser 前端 + P+ 更新，不跑密度修饰和建模。界面上的模块 1 现在含 oasis-dm 的 12 轮密度修饰，胰岛素 2.1 Å 的实测是每轮约 5 s、3 轮 16.6 s，见 `docs/OASIS_DM_REPLACEMENT_20261114.md`）。

## 6. 限制

- 精简版**没有模型、没有 R-free、没有逐轮精修记录**，不能拿来做"结构解析完成"的验收；验收指标只有相位/FOM/图相关。
- 没有本机 CCP4 时，模块 1 仍可运行（Phaser 前端在包内运行时里，密度修饰是树内 oasis-dm），但模块 2 不可用；`pplus_bundled` 本身需要包内运行时，缺运行时的裸 Python 环境只能走 `pplus_internal`（实验性）。
- 十案例的 `endpoint_*` 目录都是三段链的产物，精简版不复现它们；精简版的对照物是各案例的 `02_sad/experimental.mtz`。
- 本页所有数字都来自 macOS Intel 运行时；Apple Silicon 需要 arm64 的 CCP4。

## 7. 质量门

- `PYTHONPATH=src .conda-unified/bin/python -m pytest tests -q` → **801 passed**。
- 精简版相关测试：`tests/test_sad_pipeline.py::test_lite_plan_stops_after_phasing_and_refuses_to_build`、`tests/test_oasis_focus.py::test_lite_job_reports_two_stages_and_certifies_the_sad_phases`、`tests/test_desktop_window.py::test_density_modification_is_always_on_and_only_building_follows_ccp4`（2026-11-14 改名，替代原来的 `test_optional_modules_follow_ccp4_detection_and_gate_building`）。
