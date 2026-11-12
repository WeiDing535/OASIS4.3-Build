# 运行时两件套、单手续算与 insulin 2.0 Å 实测 — 2026-11-12

本轮起因：GUI 提交的作业（`runs/oasis_run`，job `7c8b11d7c0d9ee23`）报
“定相 → 密度修饰 → Buccaneer → REFMAC 失败（退出码 1）”，用户读到的
一手失败提示与实际原因不一致。

## 1. 失败的真实原因：运行时指向了裸 CCP4

`hand_selection.json` 里两手状态完全相同：

| 手 | status | error |
| --- | --- | --- |
| input_hand | failed | `[Errno 2] No such file or directory: '/Applications/ccp4-9/density/bin/oasis-density-python'` |
| opposite_hand | failed | 同上 |

即**不是“一手失败、一手有结果”，而是两手都因同一个缺文件立刻失败**；
`choose()` 当时只报“至少一手计算未完整完成”，把真正的错误藏了起来。

根因是 `cfg['runtime']` 指向 `/Applications/ccp4-9`。`pplus_bundled`
的 P+ 定相与双手试建都走包内 RESOLVE/CCTBX 树
（`<runtime>/density/bin/oasis-density-python`），这棵树只存在于便携包的
`runtime/` 里，外部 CCP4 安装没有。`oasis_workflow.defaults()` 原先只认
写死的 `dist/OASIS4.3-dev*-portable local/runtime` 列表，这些目录在开发树里
已不存在，于是静默回落到 `default_runtime()`（即 `/Applications/ccp4-9`），
GUI 把这条路径当成正常配置提交了。

修复：

- `bundled_pplus.density_root()/require_density()`：独立解析密度树
  （显式 `density_runtime` → `PHASEBUILD_DENSITY` → `<runtime>/density` →
  `<app>/runtime/density` → `<app>/../runtime/density`），解析不到就给出可执行
  的中文提示，而不是等 Phaser 跑完再报 `No such file`。
- `sad_pipeline.validate()`：`pplus_bundled` 分支在建目录、跑任何程序之前就
  调用 `require_density()`，失败信息直达 GUI。
- `oasis_workflow.defaults()`：候选运行时改为扫描 `project/runtime`、
  `project/dist/OASIS4.3*portable*/runtime`、`project.parent/OASIS4.3*portable*/runtime`，
  并优先选择同时带 `density/bin/oasis-density-python` 的那个。
- `desktop._refresh_module_hint()`：检测到的运行时不带 `density/` 时，在模块提示
  里直接写明“这是外部 CCP4 安装，P+ 需要便携包 runtime”。

同一条命令改前/改后对比（同一份输入）：

```
# 改前：Phaser 跑完后才崩
ValueError: 至少一手计算未完整完成，无法可靠判别；请查看候选日志
# 改后：validate 阶段立刻退出
ValueError: 这次运行缺“包内密度修饰引擎”bin/oasis-density-python：当前运行时
/Applications/ccp4-9 里没有 density/，它只是外部 CCP4 安装。请把“CCP4 位置”
指向便携包的 runtime 目录（例如 …/OASIS4.3-…-portable local/runtime，里面有
density/ 和 bin/），或直接启动便携包里的 OASIS4.3
```

## 2. 单手有结果时继续算（明确标注，不静默）

`hand_selection.choose()` 原先是“任一手 status 不是 completed/empty 就整单失败”。
现在分三种情况：

| 情况 | 行为 |
| --- | --- |
| 两手都完整 | 原逻辑不变：唯一测序残基差 ≥ 阈值且 ≥ margin，或弱图下 built 计数大比例胜出 |
| 一手失败、另一手唯一测序残基 ≥ 阈值 | **继续**，`policy.single_hand=True`、`hand_verified=False`，失败手的原始 error 写进 `policy.other_hand_error`；`frontend.json` 记 `hand_verified=false`，stderr 打印“本次手性未交叉验证” |
| 一手失败、另一手证据不足；或两手都失败 | 失败，并在错误信息里逐手列出真实原因 |

判据仍然是“这一手自己的证据够不够”，不是“有文件就用第一个文件”，所以单手续算
不会把无证据的图当成已判定手性。判定阈值不变（`max(3,min(15,ceil(.15*N)))`）。

## 3. 修好运行时之后：这份 insulin 数据本身在边缘

用便携包 runtime 重跑同一份输入（`runs/ins_portable`），作业不再因缺文件失败，
但停在手性判别：

| 手 | status | 唯一测序 | 建出残基 | 最长片段 |
| --- | --- | --- | --- | --- |
| input_hand | completed | 0 | 6 | 6 |
| opposite_hand | completed | 0 | 7 | 7 |

P+ 本身是给了子结构的：Phaser SAD 8 个 S 位点、FOM 0.336、LLG 150、Rfac 75.4；
把试建从 3 周期加到 10 周期仍然只建出 6 个残基，说明不是试建太短，而是这张
Wang8 + classic40 修饰图本身没有可用信息。

同一份数据、同一套 8 个位点交给 SHELXC/SHELXE（`/tmp/ins_shelx`，
`shelxc ins` → `shelxe ins ins_fa -s0.5 -m20 -h -a4`）：

| 指标 | 数值 |
| --- | --- |
| Pseudo-free CC | 57.46 % |
| estimated mean FOM | 0.549 |
| 最好一轮 trace CC | 13.57 % |
| trace 输出 | 21 个残基，3 段（目标 51） |
| 密度修饰 contrast | 0.392–0.404（成功案例通常 > 0.5） |

也就是说，这份 2.24 Å（名义 2.0 Å）S-SAD 数据上，SHELXE 同样只给出一段
13.6 % CC 的碎片主链；两边的结论一致——数据在边缘，不是哪一家的程序坏了。
这类情形下 OASIS4.3 的默认行为是明确报“证据不足”，而不是输出一个看起来
像模型的错误手性结构。
