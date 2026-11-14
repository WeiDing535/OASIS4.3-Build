# 运行日志只剩几行、进度条不消失：原因与修正（2026-11-13）

两个现象都出在界面层，与定相结果无关。用当时的作业目录 `runs/oasis_20261113_094302`（任务
`bd9221be647d3f97`，R-work 0.3166 / R-free 0.3779）复现并核对。

## 1 日志只有几行

`gui.log` 由 `jobs.Controller.run` 打开（`'w'` 模式），只写入两类内容：

- 每个“步骤”前的 `▶ <步骤名>`。SAD→模型只有一个步骤（`run_sad_to_model.py`），
  所以整轮只有一行步骤标题；
- 子进程的 stdout/stderr。而 `sad_pipeline.run_command` 把每个阶段的 stdout 重定向
  到 `<pipeline 目录>/<阶段 id>.log`，因此各阶段的叙述根本不到 `gui.log`。

该任务的 `gui.log` 877 字节 / 19 行，正是“步骤标题 + 最后那份 JSON 报告”。真正的
叙述在阶段日志里：

| 阶段日志 | 行数 | 内容 |
| --- | --- | --- |
| `02_sad/sad.log` | 0 | 当时阶段本身不打印任何东西 |
| `03_dm/dm.log` | 2 | `feedback cycle 1: 58 residues in 5 fragments`、`cycle 2: 30 …` |
| `05_build/build.log` | 1 | 最终 JSON（`run_build_pipeline.py` 只打印这一行） |
| `05_build/cycle_05/refmac.log` | 1033 | REFMAC 原始日志（HTML） |

修正分两处：

1. `desktop.stage_log_text()` 按 `pipeline_report.json` 的顺序（`completed` + `current`）
   读回各阶段日志，拼成 `── <阶段名> ──` 分块，每块保留末尾 60 行，被截断时写出
   `…（前面 N 行省略）`。`_log_view()` 把它接在 `gui.log` 之后，`_refresh()` 与
   `on_select_job()` 共用，界面与历史记录一致。
2. 补上阶段内的进度叙述，让阶段日志本身有话可说：
   - `ccp4_pipeline.run` 每轮打印 `build cycle 2/5: 45 residues, Rwork 0.3501, Rfree 0.3934`；
   - `bundled_pplus.phase` 在自动判手后打印两手各自的 `status` 与 `unique` 残基数，
     以及 `hand selection: input_hand (verified)`。

实测两轮构建（用该任务的 native/相位图/序列）：

```
build cycle 1/2: 40 residues, Rwork 0.4252, Rfree 0.4543
build cycle 2/2: 45 residues, Rwork 0.3501, Rfree 0.3934
```

## 2 完成后蓝色进度条不消失

`_refresh()` 在终态分支里把进度条设成 `value=100 if status=='success' else 0`，
随后 `poll()` 因为该任务已进 `self.finished` 不再刷新，于是满格蓝条一直停在原地，
只有 `on_select_job()`（切走再切回）会把它归零——这正是“点一下别的任务再切回来才消失”。

修正：终态一律 `value=0`；完成与否由右上角状态文字（“已完成 …”）与结果页表达。

顺带修掉一个会让人误判的细节：`jobs.get()` 里阶段总数原先硬编码为 5，而实际计划
有 6 个阶段（prepare/sad/dm/bridge/build/finish），第 6 个阶段还在跑时进度条已经是
100%。现在按 `plan['stages']` 计数，`completed/6` 才到满格。

## 3 回归测试

- `tests/test_desktop_window.py::test_the_log_pane_shows_each_stage_not_only_the_headers`
- `tests/test_desktop_window.py::test_a_long_stage_log_keeps_its_tail_and_says_so`
- `tests/test_desktop_window.py::test_a_finished_job_clears_the_progress_bar`
- `tests/test_jobs.py::test_progress_counts_the_plans_own_stages`

全量：675 passed。
