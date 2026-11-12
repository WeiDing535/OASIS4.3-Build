# 精简版：只做模块 1（P+ 定相）

这一页回答两个问题：**十案例是不是模块 1 单独跑出来的**，以及**精简版（只留模块 1）能交付什么**。界面怎么排、CCP4 怎么检测，也写在同一页。

## 1. 先更正一个前提

十案例 **不是** 第一个模块单独完成的。十案例用的协议是

```
--method pplus_bundled --dm-engine parrot --dm-cycles 40 \
--feedback-cycles 2/3 --trust-wang-radius 8 --dm-solvent 0.5 \
--no-auto-hand --input-columns
```

即 **模块 1（P+ 定相）+ 模块 2（Parrot 密度修饰）+ 模块 3（Buccaneer/REFMAC）** 三段串起来的结果。三个环节缺一不可：

| 只跑到 | 手上有什么 | 3O2E 实测 | 十案例均值 |
| --- | --- | --- | --- |
| 模块 1（P+） | 相位 PHIB、FOM、PPLUS 的 MTZ | complex CC **0.123**，mean FOM 0.482 | **0.272**（另一个 P+ 后端 0.279） |
| 模块 1+2（Parrot） | 密度修饰后的图与相位 | 0.549 | **0.688** |
| 模块 1+2+3 | 模型、R/R-free、逐轮精修记录 | R/R-free 41.14 / 42.64 | 10/10 落在对照 ±5 pp 内 |

"complex CC" 是对 campaign 精修模型的平移不变复空间图相关（`work/noCCP4/compare_map_phases.py`），不是 omit-map，也不是振幅相关。

所以：**模块 1 给的是相位，不是模型**。没有模型就没有 R-free，也没有"十案例 10/10"这个说法——那个数字属于三段链。模块 1 单独的价值是把可用的初始相位做出来，让模块 2 有一个够好的起点；按图相关算，起点离成品还差 0.4 左右。

## 2. 精简版是什么

界面把流程拆成三个可勾选的模块：

| 模块 | 程序 | 依赖 | 可关闭 |
| --- | --- | --- | --- |
| ① P+ 定相 | 包内 Phaser 前端 + P+ 信任域更新 | 只要包内运行时 | 否（必做） |
| ② Parrot 密度修饰 | `cparrot` | 本机 CCP4 | 是 |
| ③ Buccaneer 建模 + REFMAC 精修 | `cbuccaneer` + `refmac5` | 本机 CCP4 | 是 |

- 关掉 ②：流程在 `02_sad` 结束，输出 `02_sad/experimental.mtz`（`H K L PHIB FOM PPLUS`）与 `observations.mtz`，报告里 `dm_phases` 为空、`completed` 只有 `prepare`、`sad` 两段。
- ③ 依赖 ②：Buccaneer 要从密度修饰后的图开始建，直接从 P+ 原始相位建会崩。所以关掉 ② 时 ③ 自动关闭并置灰，界面给一句"模块 2 已关闭……模块 3 需要修饰后的图"。
- 精简版不做自动判手（`auto_hand` 需要试建），手性由重原子亚结构决定；界面上的"自动判手"随 ③ 一起关闭。

命令行等价写法（三种任选）：

```
python run_sad_to_model.py --lite ...            # 只做 P+ 定相
python run_sad_to_model.py --no-dm-enabled ...   # 同上，显式关模块 2
python run_sad_to_model.py --no-build-model ...  # 只到密度修饰
```

`--lite` 会同时把 `dm_enabled`、`build_model`、`auto_hand` 置为关闭，并把 `feedback_cycles` 归零——trace 反馈本来就是"用部分模型算 Fcalc 回喂密度修饰"，没有建模就没有反馈。

## 3. 界面：参数分组与模块开关

表单自上而下：

1. **任务名称**
2. **运行模块（按本机 CCP4 检测结果默认勾选）**：模块 1 固定勾选并置灰；模块 2、3 各一个复选框，旁边是「重新检测 CCP4」按钮和一行检测结论。
3. **输入文件**：SAD MTZ、重原子 PDB、序列、元素、波长、份数、新输出目录。
4. **① P+ 定相参数**：P+ 后端单选、P+ 参数方案、定相迭代周期、每反射 KL 上限、每个目标的关系上限。
5. **② Parrot 密度修饰参数**：密度修饰引擎、溶剂比例、密度修饰周期、掩膜过滤半径、trace 反馈轮数。
6. **③ Buccaneer 建模 + REFMAC 精修参数**：自动判手、Buccaneer/REFMAC 轮数、每轮 REFMAC 周期、Buccaneer/REFMAC 运行目录、构建配方、精修相位来源。
7. **显示高级参数**（默认收起）：输入细节、P+ 掩膜与采样、运行资源与终点组合。

关闭某个模块时，它下面的参数框整组置灰——不是隐藏，让人看得见"关了哪些"。"Buccaneer / REFMAC 运行目录"是例外，始终可编辑，因为检测 CCP4 时也会读它。

## 4. CCP4 检测

`src/phasebuild/module_detection.py` 只回答一个问题：**这台机器上有没有能跑的可执行文件**。找不到是正常答案，不抛异常。

按顺序在这些目录里找（命中即止）：

1. 表单里"Buccaneer / REFMAC 运行目录"指向的目录及其 `bin/`（打包运行时走这条）；
2. 环境变量 `CBIN`、`CCP4_BIN`、`OASIS_CCP4_BIN`；
3. 环境变量 `CCP4`、`CCP4_ROOT` 及其 `bin/`；
4. `/Applications/ccp4-*/bin`、`/opt/ccp4*/bin`、`/usr/local/ccp4*/bin`（双击启动的 GUI 几乎继承不到环境变量，所以显式找安装位置）；
5. `PATH` 的每一项。

判定规则：

| 模块 | 判据 |
| --- | --- |
| ② Parrot | `cparrot` 或 `parrot` 存在且可执行 |
| ③ Buccaneer/REFMAC | `cbuccaneer` 或 `buccaneer` **且** `refmac5` 都可执行 |

结论只有三种，直接写在那行提示里：

- 都找到：模块 2、3 默认勾选，提示"已检测到 CCP4：模块 2、3 默认勾选，可随时取消。"
- 只找到一部分：默认勾选可用那个，提示缺哪个。
- 都没找到：两个都关闭并置灰，提示"未检测到 CCP4（cparrot、cbuccaneer、refmac5），默认只做模块 1。"，第二行说明模块 1 输出相位/FOM/MTZ，以及装好 CCP4 或设好 `CCP4`/`CBIN` 后点「重新检测」。

重新用"沿用参数"载入一条旧记录时，记录里的开关优先；但如果记录要求一个本机现在没有的模块，复选框会被压回关闭状态——不会留一个点了必然失败的勾。

检测结果只影响**默认值**，不影响可执行性判断：用户仍可手动勾选（模块 2 可勾），真正运行前的校验在 `sad_pipeline.validate()`，那里会给出具体缺哪个程序。

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

耗时：该案例约 2 分钟（模块 1 只做 Phaser 前端 + P+ 更新，不跑 Parrot 的 40 轮和 5 轮建模）。

## 6. 限制

- 精简版**没有模型、没有 R-free、没有逐轮精修记录**，不能拿来做"结构解析完成"的验收；验收指标只有相位/FOM/图相关。
- 没有本机 CCP4 时，模块 1 仍可运行（Phaser 前端在包内运行时里），但模块 2、3 不可用；`pplus_bundled` 本身需要包内运行时，缺运行时的裸 Python 环境只能走 `pplus_internal`（实验性）。
- 十案例的 `endpoint_*` 目录都是三段链的产物，精简版不复现它们；精简版的对照物是各案例的 `02_sad/experimental.mtz`。
- 本页所有数字都来自 macOS Intel 运行时；Apple Silicon 需要 arm64 的 CCP4。

## 7. 质量门

- `PYTHONPATH=src .conda-unified/bin/python -m pytest tests -q` → **801 passed**。
- 精简版相关测试：`tests/test_sad_pipeline.py::test_lite_plan_stops_after_phasing_and_refuses_to_build`、`tests/test_oasis_focus.py::test_lite_job_reports_two_stages_and_certifies_the_sad_phases`、`tests/test_desktop_window.py::test_optional_modules_follow_ccp4_detection_and_gate_building`。
