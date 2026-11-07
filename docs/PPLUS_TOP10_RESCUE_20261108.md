# P+ 十套困难案例：密度修饰与 trace 反馈升级（2026-11-08）

## 1. 问题定位

用同一套子结构、同一份数据把 SHELXE 拆成三档重跑，可以确定差距来源：

| SHELXE 运行方式 | 6ZPP 加权平均相位误差 |
|---|---:|
| `-m0`（无密度修饰、无追踪） | 64.5° |
| `-m20`（20 轮密度修饰、无追踪） | 61.4° |
| `-a 4`（生产运行，含自动追踪） | 26.9° |

也就是说：**密度修饰几乎不改变相位，全部收益来自自动链追踪与"追踪模型→相位→再密度修饰"的循环。**
原 P+ 流程的密度修饰（RESOLVE Wang14 + CCTBX classic40）只把 Phaser SAD 的 61.1° 降到 51.9°，
而且没有任何追踪反馈，因此上游相位差距一直保留到建模。

## 2. 本次改动

### 2.1 密度修饰（`src/phasebuild/parrot_dm.py`）

- 各向异性标度校正（`-anisotropy-correction`）：6ZPP 上 45.1° → 41.0°。数据确实存在各向异性，
  关掉校正后高分辨率误差从 46° 恶化到 59°。
- 溶剂掩膜滤波半径 3.0 Å（`-solvent-mask-filter-radius`）：41.0° → 40.8°，与各向异性叠加后 **37.1°**。
- 溶剂含量可显式覆盖（本批统一 0.50，序列估计值为 0.37）。
- 只喂工作集反射与原始观测 HL；自由集在调用前剔除。

### 2.2 Trace/Fcalc 反馈（`src/phasebuild/trace_feedback.py`）

每轮：Parrot → 自带 `cbuccaneer` 工作集追踪 → 部分模型 Fcalc（REFMAC 批量溶剂校正）
→ 与**固定**的原始观测 HL 合并 → 再密度修饰。每轮单独目录与报告，不累积观测似然。

6ZPP：37.1° → 31.4°（1 轮）→ 32.0°（2 轮，已收敛）。
Buccaneer 从这张图能建出的残基数由旧流程的 56 个增加到 147–162 个。

### 2.3 构建器选项（`src/phasebuild/ccp4_pipeline.py`）

`--building-recipe pipeline` 让 P+ 路线使用与对照流程**完全相同的** `buccaneer_pipeline`
（构建 + REFMAC + 其内部密度修饰循环）。这样对照的是"相位差异"而不是"构建器差异"。
自带的 `run_build_pipeline.py` 输出沿用的 `hybrid_report.json` 契约。

### 2.4 新选项（默认全部保持旧行为）

```
--dm-engine resolve|parrot
--dm-solvent --dm-mask-radius --dm-anisotropy/--no-dm-anisotropy
--feedback-cycles 0..8 --trace-resolution
--building-recipe legacy|ccp4|pipeline
```

## 3. 十套案例结果

固定协议：工作集上限 `--max-per-target 48`，溶剂 0.50，掩膜半径 3.0，各向异性校正。
每套跑 feedback=1 与 feedback=2 两个固定协议终点，按**预先固定的工作集 R 因子**（不读 R-free）
选择终点。R-free 只用于最终评估。

| 案例 | SHELXE 对照 | 原 P+ | 新 P+ R-free | 差值 | 反馈轮次 |
|---|---:|---:|---:|---:|---:|
| 10_6zpp_3se_P1211 | 31.11% | 67.41% | 31.83% | +0.7 pp | 2 |
| 48_6HF7_4TB2MG_P4212 | 28.43% | 57.23% | 30.02% | +1.6 pp | 1 |
| 5_5D9D_3HG_P212121 | 29.73% | 54.95% | 27.12% | −2.6 pp | 2 |
| 14_7tj1_16SE_P212121 | 32.49% | 54.03% | 30.42% | −2.1 pp | 1 |
| 11_7e1d_4se_P212121 | 38.20% | 58.80% | 41.07% | +2.9 pp | 2 |
| 24_7tfq_11SE_P21212 | 37.29% | 55.41% | 38.78% | +1.5 pp | 2 |
| 13_5hhk_33se_p1211 | 27.40% | 43.60% | 28.91% | +1.5 pp | 1 |
| 3_3o2e_9I_1K_P42212 | 42.85% | 56.07% | 38.98% | −3.9 pp | 2 |
| 7_7enm_3se_P212121 | 30.08% | 41.96% | 32.09% | +2.0 pp | 1 |
| 19_7tbs_3SE_I222 | 30.58% | 40.51% | 31.14% | +0.6 pp | 2 |

**10/10 在 ±5 个百分点以内，3 套优于对照**；差值中位数 +1.55 pp（原 P+ 中位数约 +11 pp）。

### 6HF7 的输入修复

6HF7 声明为 Tb/Mg 双元素，SHELXD 的 `.res` 把全部位点都标成 Mg（f″≈0.20，而 Tb≈10.60）。
用按占有率排序重标注的位点（`shelxd_sites_expected_occupancy_typed.pdb`）重跑后，
密度修饰后相位误差从 83.5° 降到 26.9°，R-free 从 57.59% 降到 30.02%。逐位点元素身份仍属假设，
必须与实验化学信息核对。

## 4. 端到端验证（6ZPP）

单条命令即可复现（见下），流水线自身输出 **R-free 32.35%**，对照 31.11%。
其中 `06_finish` 的序列清理被工作集守卫拒绝，保留固定最后一轮模型（符合既定规则）。

```
run_sad_to_model.py --sad ... --native ... --sequence ... --heavy ... \
  --method pplus_bundled --input-columns --wavelength 0.97918 --elements Se \
  --dm-engine parrot --feedback-cycles 2 --dm-solvent 0.5 \
  --trust-wang-radius 8 --max-per-target 48 --building-recipe pipeline
```

5D9D 的端到端独立验证见 `work/exp/integrated_final_5d9d`。

## 5. 局限与风险

- **Parrot 是外部引擎**（本机 CCP4-9）。自带 runtime 未含 `cparrot`；若要维持完全自洽的产品形态，
  需要把 `oasis_py/native_density.py` 的掩膜与各向异性处理升级到同等水平——当前内部实现在
  6ZPP 上只有 57.3°。`--dm-engine resolve` 保留旧路径。
- 相位误差在高分辨率仍偏高（6ZPP 1.5 Å 处约 46°，SHELXE 约 28°）。
  继续缩小差距需要更强的独立链追踪器，而不是继续调密度修饰参数。
- 元素逐位点标注（6HF7）是假设性输入修复，需实验核对。
- 每套案例的两个终点（feedback=1/2）按工作集 R 因子选择；该规则在运行前固定，不读 R-free。
- 本批只对这十套做了验证，未重跑 100 套全集。
