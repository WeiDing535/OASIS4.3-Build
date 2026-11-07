# 密度修饰与 trace 反馈升级（2026-11-08）

## 结论摘要

对照流程的优势几乎全部来自 SHELXE 的自动链追踪，而不是密度修饰本身。用同一套
子结构从零重跑 SHELXE：

| 运行 | 6ZPP 加权平均相位误差 |
|---|---:|
| `-m0`（无密度修饰、无追踪） | 64.5° |
| `-m20`（20 轮密度修饰、无追踪） | 61.4° |
| `-a 4`（完整生产运行，含自动追踪） | 26.9° |

原流程的密度修饰（RESOLVE Wang14 + CCTBX classic40）把 Phaser SAD 的 61.1° 只降到
51.9°。本升级把同一先验降到 37.1°（单轮），再加一轮工作集 trace/Fcalc 反馈降到 31.4°。

## 新增内容

1. `src/phasebuild/parrot_dm.py`
   - CCP4 Parrot 适配器：各向异性标度校正、`-solvent-mask-filter-radius 3.0`、
     可配置溶剂含量、工作集限定。
   - 只接收持久化的原始观测 HL（以及反馈轮次中工作集部分模型的 Fcalc）；
     自由集反射在调用前剔除；不读取 R-free 做任何选择。
   - 输出沿用流水线的 `map_phases.mtz` 契约（F,SIGF,FWT,PHWT,PHIB,FOM,HLA..HLD）。

2. `src/phasebuild/trace_feedback.py`
   - 工作集 trace → 部分模型 Fcalc（REFMAC 批量溶剂校正）→ 与固定实验似然合并
     → 再密度修饰的循环；每轮单独保存目录与报告。

3. `run_trace_feedback.py`：流水线阶段入口，读取 `pipeline_report.json` 的配置，
   写入 `03_dm/map_phases.mtz` 与 `map_phases.json`。

4. 新增选项（默认保持旧行为）：
   - `--dm-engine resolve|parrot`（默认 `resolve`）
   - `--dm-solvent`、`--dm-mask-radius`、`--dm-anisotropy/--no-dm-anisotropy`
   - `--feedback-cycles N`、`--trace-resolution`、`--dm-executable`

## 触发条件与依赖

- `--dm-engine parrot` 需要本机可执行的 `cparrot`（CCP4）。程序会在
  `--dm-executable`、自带 runtime 的 `bin/`、`PATH` 和常见 CCP4 安装位置查找。
- 反馈轮次使用自带 runtime 的 `cbuccaneer` 与 `refmac5`；无需额外安装。
- 默认 `resolve` 引擎与 `--feedback-cycles 0` 时，行为与升级前完全一致。

## 验证（6ZPP，同一自由集，自由值 0 为自由集）

| 流程 | Rfree | 残基 |
|---|---:|---:|
| SHELXE 相位 + buccaneer_pipeline（对照） | 31.11% | 155 |
| 原 P+ 相位 + 原密度修饰 | 67.41%（历史最佳 42.33%） | 12 |
| 本升级：Parrot(各向异性)+1 轮 trace 反馈 | 33.92% | 148 |

## 局限

- 相位误差仍集中在高分辨率（1.5–1.9 Å）：31.4° 总体，但 1.5 Å 处约 46°；
  SHELXE 在各分辨率均为 27–28°。要继续缩小差距需要更强的独立链追踪器，
  而不是继续调密度修饰参数。
- 本轮只在 6ZPP 上完成端到端验证；其余九套的批量结果见 `work/` 下的实验记录。
- Parrot 是外部引擎；若要维持"自带 runtime 完全自洽"的产品形态，需要把
  `oasis_py/native_density.py` 的掩膜与各向异性处理升级到同等水平。当前内部
  实现在 6ZPP 上只有 57.3°。
