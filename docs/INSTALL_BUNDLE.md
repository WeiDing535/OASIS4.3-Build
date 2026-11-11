# 怎样获得“免安装”的 OASIS4.3 完整包

## 结论先说

**目前没有公开下载地址。** 公开仓库只有源码（约 5 MB）；完整包在 1.4–3.8 GB 之间，里面带着 Phaser、Parrot、Buccaneer、REFMAC 等第三方二进制，包内 `bundle.json` 明确写着：

```json
"public_redistribution_ready": false
```

所以它不能直接挂到公开仓库或公开 Release 上。已有 CCP4 授权的用户，可以按下面的步骤在本机组装出一模一样的包——**只要本机装了 CCP4 9，就不需要任何额外安装**。

## 方式 A：自己组装（推荐，许可上最干净）

前提：macOS Intel；一份 Python 3.12（含 numpy / scipy / gemmi）；本机已安装 CCP4 9（提供 `cparrot`、`phaser`、`cbuccaneer`、`refmac5`）。全程只做复制，不安装任何东西，也不改写来源安装。

```bash
git clone https://github.com/WeiDing535/OASIS4.3-Build.git
cd OASIS4.3-Build

# 一份 3.12 科学环境（已有可跳过，改用你自己的环境路径）
conda create -y -n oasis-bundle python=3.12 numpy scipy gemmi

# 1) 从你自己的 CCP4 生成私有 runtime（约 40 秒）：
#    复制可执行文件、动态库、完整 lib/data 与许可文本，并写出 runtime_manifest.json
#    --engines 决定带哪些引擎：默认只带 cbuccaneer/refmac5（内部 DM 路线）；
#    要跑十案例那套 Parrot 协议，把 phaser 和 cparrot 一起带上
python packaging/build_local_runtime.py \
  --source /Applications/ccp4-9 \
  --engines cbuccaneer,refmac5,phaser,cparrot \
  --output dist/runtime

# 2) 把「源码 + 这份 Python + 这份 runtime」装进一个免安装目录（约 1 分钟）
python packaging/build_installation_free.py \
  --python-source "$(conda info --base)/envs/oasis-bundle" \
  --runtime-source dist/runtime \
  --output dist/OASIS4.3-portable

# 3) 自检：脚本会用干净环境（不加载系统 PYTHONPATH / CCP4）跑一次 doctor
dist/OASIS4.3-portable/oasis doctor
```

产物结构：

```
dist/OASIS4.3-portable/
├── app/                 # 源码（package_oasis.py 选出的文件）
├── python/              # Python 3.12 + numpy/scipy/gemmi（随包）
├── runtime/             # bin/{phaser,cparrot,cbuccaneer,refmac5} + lib/data
├── oasis                # 统一入口：./oasis [gui|sad-to-model|portfolio|doctor]
├── OASIS4.3.command     # 双击开窗口
├── bundle.json          # 版本、内置引擎、许可状态
├── DOCTOR.json          # 自检结果
└── 使用说明.txt
```

双击 `OASIS4.3.command` 即开原生窗口，不需要另装 Python、CCP4 或 Phenix。

`bundle.json` 里的引擎清单和 `default_method` 是从 `runtime/runtime_manifest.json` 里实际带上的可执行文件推导出来的，所以四个引擎都在时它写 `"default_method": "pplus_bundled"`、`"bundled_engine_ids": ["phaser","cparrot","cbuccaneer","refmac5"]`；只带 Buccaneer/REFMAC 时写 `pplus_internal`。

**实测记录**（2026-11-11，macOS Intel）：

- 第 1 步：41.9 秒，1.17 GB，35 878 个文件，`bin/` 恰好是请求的四个引擎，启动检查全部通过。
- 第 2 步：1 分 05 秒，产物 1.4 GB = `python` 262 MB + `runtime` 1.2 GB + `app` 5 MB；干净环境下 `doctor` 为 `passed`。

## 方式 A 之后：跑通十案例那套 P+ 协议

组装出来的包只带四个引擎，**没有** `runtime/density/`（RESOLVE/CCTBX，Phenix 派生）。这不影响 P+ 协议，因为该协议用的是 Parrot，而 Parrot 已经取代了整个 RESOLVE/CCTBX 密度修饰；`engine_layout()` 现在只在 `dm_engine` 不是 `parrot` 时才要求 `runtime/density/`。

以 7ENM 为例（输入取自 `cases/7_7enm_3se_P212121`）：

```bash
./oasis sad-to-model \
  --method pplus_bundled --dm-engine parrot \
  --dm-cycles 40 --feedback-cycles 2 --trust-wang-radius 8 --dm-solvent 0.5 \
  --no-auto-hand --input-columns --elements Se --wavelength 0.979 \
  --sad     <case>/01_inputs/selected_sad.mtz \
  --native  <case>/01_inputs/native.mtz \
  --sequence <case>/01_inputs/building.fasta \
  --heavy   <case>/02_sad/refined_heavy.pdb \
  --reference "$PWD/runtime/lib/data/reference_structures/reference-1tqw.mtz" \
  --output /new/empty/path
```

实测（同一台机器，干净环境 `env -i`，只用包内引擎）：5 分 15 秒跑完，`primary_metric` 为 **Rwork 28.34 / Rfree 31.47**，与十案例基准的 7ENM（28.34 / 31.47）完全一致。终点仍按工作集 R 选择，全程不读 R-free。

两个容易踩的坑：

- **不要加 `--auto-hand`。** 双手判别阶段（`hand_selection`）仍然调用 RESOLVE 密度修饰，需要 `runtime/density/bin/oasis-density-python`；CCP4-only 包里没有这个文件，会报 `No such file or directory` 并终止。十案例基准本身用的就是 `--no-auto-hand`（手性由已精修的重原子子结构确定）。
- **要加 `--input-columns`。** 基准用了它，会在 `01_inputs/` 里写出 `columns.json` 和 `native.mtz`；不加的话 `native.mtz` 不进 `01_inputs`，trace/Fcalc 反馈阶段会报 `Missing stage input .../01_inputs/native.mtz`。

需要 `--dm-engine resolve`（内置 RESOLVE/MLHL 密度修饰）或 `--auto-hand` 时，才需要带 `runtime/density/` 的完整 dev17 包；那份密度树由 `packaging/build_bundled_pplus.py --density-source` 组装。

注意事项：

- 两个脚本都拒绝覆盖已存在目录，也拒绝把输出放在来源 runtime 内部。
- 目标机器必须与构建机器同架构：这里产出的是 macOS Intel x86_64；Apple Silicon 需要先有 arm64 的 CCP4 运行时。
- 这个包只在你（或你的组）已有的授权范围内使用，**不要把它当公开发行包再分发**。

## 方式 B：从维护者处拿现成构建（私有分发）

只带四个引擎的包压缩后很小，不需要分卷：

```bash
tar -czf OASIS4.3-portable-ccp4only.tar.gz OASIS4.3-portable-ccp4only-20261111
shasum -a 256 OASIS4.3-portable-ccp4only.tar.gz > OASIS4.3-portable-ccp4only.tar.gz.sha256
```

实测：1.4 GB 的目录压成 **436 MB**、45 370 个条目（55 秒），单个文件即可传，接收方：

```bash
shasum -a 256 -c OASIS4.3-portable-ccp4only.tar.gz.sha256
tar -xzf OASIS4.3-portable-ccp4only.tar.gz
cd OASIS4.3-portable-ccp4only-20261111 && ./oasis doctor
```

- 带 RESOLVE/CCTBX 的 dev17 完整包（3.8 GB）压缩后仍可能超过 2 GB 的 GitHub Release 附件上限，那才需要 `split -b 1900m … part-` 分卷，接收方 `cat part-* | tar xz`。
- 渠道用内网、移动硬盘，或**私有**仓库的 Release；不要放到公开仓库。
- 包里带着 CCP4 的二进制，按 CCP4 许可，接收方也应已持有相应授权；接收机需为 macOS Intel。

## 方式 C：只要代码

```bash
python -B package_oasis.py --output /new/path.zip
```

生成精简源码包（**不含** Buccaneer/REFMAC/Parrot 运行时）。在已装 CCP4 的机器上可直接从源码树运行；缺 Parrot 时会回退到内置引擎，并明确报出回退原因（见 [README](../README.md) 与 `docs/LINUX.md`）。

## 为什么不能直接给一个公开链接

- 第三方引擎的许可不允许把它们的二进制公开再分发；`bundle.json` 的 `public_redistribution_ready: false` 就是这个意思。
- 即使许可解决，含 RESOLVE/CCTBX 的完整包 3.8 GB，超过单个 Release 附件的 2 GB 上限，需要分卷与校验清单。
- Linux 目前还没有完整运行时；`docs/LINUX.md` 记录了在 Linux 上准备原生引擎的进展。

## 相关文件

- `packaging/build_local_runtime.py` —— 从本机 CCP4 生成私有 runtime（带 manifest 与许可文本，`--engines` 选引擎）
- `packaging/build_installation_free.py` —— 组装免安装包并跑干净环境自检
- `packaging/audit_linux_runtime.py` —— Linux 运行时审计（尚未提供完整包）
- `package_oasis.py` —— 精简源码包
