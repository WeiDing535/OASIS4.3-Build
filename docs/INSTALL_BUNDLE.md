# 怎样获得“免安装”的 OASIS4.3 完整包

## 结论先说：只装一份 CCP4 就能跑完整条链

**不需要 `density/`，也不需要我们分发任何二进制。** 本机装一份 CCP4（学术用户免费），新建任务时：

1. 「CCP4 位置」填本机安装，例如 `/Applications/ccp4-9`；
2. 方法选「连续 SAD（只需 CCP4）」；
3. 点「开始」。

跑的是：连续 SAD 定相 → 密度修饰 → Buccaneer 构建 → REFMAC 精修，全部调用 CCP4 自带的 `cparrot` / `cbuccaneer` / `refmac5`，对标 SHELXC/D/E 那条独立链。

胰岛素 2.1 Å SAD 实测（`ins_wave2.0_reso2.1A_sad.mtz` + `sg.pdb` + `3INS.pir`，λ=2.0 Å，S 边）：

| 路线 | 需要什么 | 5 轮后 Rwork / Rfree | 残基 | 用时 |
| --- | --- | --- | --- | --- |
| 连续 SAD + 内部 MLHL-best 密度修饰 | 只要 CCP4（`cbuccaneer`/`refmac5`） | **0.2796 / 0.3331** | 50 | 12 min |
| 连续 SAD + Parrot 密度修饰 | 只要 CCP4（多一个 `cparrot`） | 0.3150 / 0.3583 | 51 | 6.5 min |
| P+ 自带引擎 + `cparrot` 密度修饰（无 `density/`） | 只要 CCP4（多一个 `phaser`） | **0.3166 / 0.3779** | 48 | 17 min |
| P+ 自带引擎（完整便携包） | CCP4 + 一份 `density/` | 0.3166 / 0.3779 | 48 | 18.1 min |

三条路线都是同一份 CCP4（`/Applications/ccp4-9`，`cparrot`/`cbuccaneer`/`refmac5`）跑出来的，命令与日志都在 `runs/oasis_ccp4_only_*/` 里，可复算。

两条 CCP4 路线的差别要说清楚：

- **内部 MLHL-best（默认）**：5 轮 R 逐轮下降（0.4735→0.4061→0.3779→0.3296→0.2796），末轮 0.2796 / 0.3331，最稳。
- **Parrot + 2 轮 trace/Fcalc 反馈**：快（6.5 min 对 12 min），第 3 轮一度到 0.2670 / 0.3092，但第 4、5 轮回升到 0.3150 / 0.3583。本流程的模型选择规则是“固定末轮、不用 Rfree 挑模型”（`hybrid_report.json` 里 `selection` 字段写明），所以报告取末轮，不取第 3 轮。

所以默认推荐仍是内部 MLHL-best；Parrot 留在高级参数里，作为“想更快、或想让 CCP4 的 Parrot 参与密度修饰”的选项。P+ 的价值在十案例那一组配方，不在这一个例子——这一例里它 0.3166 / 0.3779，反而最高。

**P+ 也不再需要 `density/`。** P+ 的定相与似然精修本来就是树内 Python（`oasis_py`），`density/` 只承担两件事：双手判别用的 Wang8 试点，和旧的 RESOLVE/CCTBX 密度修饰后端。把密度修饰交给 CCP4 的 `cparrot`、手性按输入重原子文件定（不做双手判别），这条 P+ 链就只要 `phaser` + `cparrot` + `cbuccaneer` + `refmac5`。上表第四行是实测：与完整便携包（第五行）的 `experimental.mtz`、`map_phases.mtz` 逐字节相同（sha256 一致），末轮 R 也一样——也就是说这一例里去掉 `density/` 没有任何损失。GUI 在“只有 CCP4”的机器上点「开始」会自己这么走，并在红条里写明“保留 P+ 定相、改用 cparrot 做密度修饰、手性按输入重原子文件（本次不做双手判别）”。

只有想保留**双手判别**（`auto_hand`，用 CCTBX 试两只手哪只对）或旧的 RESOLVE 后端时，才需要那份 `density/`；它不能随包公开分发。没有它时，两条 CCP4 路线照样出结构。

## 以前为什么“填了 CCP4 也跑不了”

旧版把默认方法设成 `pplus_bundled`，而那条路第一步就要 `density/bin/oasis-density-python`；`/Applications/ccp4-9` 自然没有，于是「CCP4 位置」填得再对也起不来，程序还会把运行时悄悄换成便携包——用户看到的结论就是“填这个没意义”。这是设计错误，已修正：现在检测到“只有 CCP4”时，有 `phaser` 就保留 P+ 并把密度修饰换成 `cparrot`，没有 `phaser` 才改走连续 SAD；两种都在红条里写明，不偷换运行时，也不再拒绝运行。

还有一处旧行为要说明：如果便携包就放在程序旁边，`density_root()` 会把它找出来借用（这是为了让便携包用户不填路径也能跑）。现在借用时会在日志里打印“density/ 取自 …（不在所选运行时内）”，不再让它看起来只用了一台机器的运行时。

## 为什么完整包不能上传 GitHub

这是许可问题，不是技术问题，代码里绕不过去：

- `runtime/` 里是 **CCP4 的二进制**（`phaser`、`cparrot`、`cbuccaneer`、`refmac5`）加上 CCP4 的 `lib/data`（monomer library、`reference-1tqw.mtz` 等）。CCP4 对学术用户免费，但**再分发要单独授权**；`runtime/runtime_manifest.json` 与包内 `bundle.json` 都记着 `"public_redistribution_ready": false`，所以它不能进公开仓库、也不能挂公开 Release。
- `runtime/density/` 是 **Phenix 派生**的 RESOLVE/CCTBX 树（`oasis-density-python`、`conda_base`、`libtbx_env`、`ext_ref_files`）。Phenix 本身不允许再分发，学术版也要按机构申请，所以这份树同样不能公开。
- 公开仓库能放的是**源码和组装脚本**（约 5 MB）：`packaging/build_local_runtime.py`、`build_installation_free.py`、`build_bundled_pplus.py`。谁有授权，谁在本机跑一遍就能得到同样的包。

所以用户按“手上有什么”选：

| 手上有什么 | 怎么做 | 能跑到哪一步 |
| --- | --- | --- |
| 只有一份 CCP4 | 「CCP4 位置」填它，方法选「连续 SAD（只需 CCP4）」 | 完整链（SAD → 密度修饰 → Buccaneer → REFMAC），实测 0.2796 / 0.3331 |
| CCP4 + 一份 `density/` | 方法选「P+ 自带引擎」，在「包内密度引擎」填那份 `density/` | 完整链，外加 P+ 标定配方 |
| 源码 | 克隆仓库，`./oasis gui` | 同上，取决于本机有没有 CCP4 |

换句话说：**runtime 不上 GitHub，但“怎么造出 runtime”全部在 GitHub 上；而只装 CCP4 的用户已经不需要 runtime 了。**

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

### 双击 `OASIS4.3.command` 会留一个终端窗口；想要图标就生成 `.app`

Finder 把 `.command` 当终端文档：双击会开一个 Terminal 窗口，那个窗口关掉，程序也跟着退出。想要“像普通 Mac 程序那样双击图标”，在同一个目录里生成一个 `.app` 外壳：

```bash
python3 packaging/make_macos_app.py                  # 生成 ./OASIS4.3.app（源码树里）
python3 packaging/make_macos_app.py --tree dist/OASIS4.3-portable   # 给便携包也生成一个
```

`.app` 里没有第二份运行时，`Contents/MacOS/OASIS4.3` 只做三件事：找到同目录的 `OASIS4.3.command`（先看 `OASIS43_TREE` 环境变量，再看包内记录的路径，最后看 `../../../`），`cd` 过去，`exec /bin/bash OASIS4.3.command`，所以图标和脚本永远不会走成两套代码。标准错误写进 `~/Library/Logs/OASIS4.3-app.log`；找不到脚本时用 `osascript` 弹窗提示，而不是静默消失。图标是 `packaging/make_macos_app.py` 里用 Pillow 现画的（深蓝圆角底 + 衍射环 + `OASIS 4.3`），`iconutil` 编译成 `.icns`，不依赖外部美术资源。

把 `OASIS4.3.app` 拖到程序坞或“应用程序”文件夹都行；**不要**只把它单独拷走——它按相对位置找源码树，换地方后重新生成一次即可。`.app` 是生成物，已写进 `.gitignore`，不入库。

`bundle.json` 里的引擎清单和 `default_method` 是从 `runtime/runtime_manifest.json` 里实际带上的可执行文件推导出来的，所以四个引擎都在时它写 `"default_method": "pplus_bundled"`、`"bundled_engine_ids": ["phaser","cparrot","cbuccaneer","refmac5"]`；只带 Buccaneer/REFMAC 时写 `pplus_internal`。

**实测记录**（2026-11-11，macOS Intel）：

- 第 1 步：41.9 秒，1.17 GB，35 878 个文件，`bin/` 恰好是请求的四个引擎，启动检查全部通过。
- 第 2 步：1 分 05 秒，产物 1.4 GB = `python` 262 MB + `runtime` 1.2 GB + `app` 5 MB；干净环境下 `doctor` 为 `passed`。

## 方式 A 之后：跑通十案例那套 P+ 协议

组装出来的包只带四个引擎，**没有** `runtime/density/`（RESOLVE/CCTBX，Phenix 派生）。**因此它单独跑不了 `--method pplus_bundled`**：P+ 的似然精修由 `density/bin/oasis-density-python` 驱动，`bundled_pplus.phase()` 一开始就检查这份树，缺了就报 `No such file or directory` 并停在 Phaser 之前。`dm_engine=parrot` 只说明密度修饰用 Parrot 而不是 RESOLVE，它省不掉 P+ 那份树。

想用这个包跑 P+ 协议，先补一份 `density/`，两种办法：

- 本机另有一份带 `density/` 的完整包（dev17 那类）：在「CCP4 位置」填这个只带四引擎的包、在「包内密度引擎」填那份 `density/` 的路径，或用命令行 `--runtime <四引擎包>/runtime --density-runtime <完整包>/runtime/density`；
- 手头有 Phenix：用 `packaging/build_bundled_pplus.py --density-source <Phenix 安装> …` 从你自己的 Phenix 生成一份 `density/`（受 Phenix 许可约束，只在你的授权范围内使用）。

下面的 7ENM 实测是在**同时有 `density/` 的机器上**做的（`density_root()` 会从应用目录旁找到它）；只有四个引擎的包不补 `density/` 时，同一条命令会在启动阶段就停下。

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

### 窗口里的三个模块与 CCP4 检测

`./oasis gui` 的表单把流程拆成三个模块：① P+ 定相（必做）、② Parrot 密度修饰、③ Buccaneer/REFMAC。后两个是 CCP4 程序，启动时按以下顺序找可执行文件：表单里「CCP4 位置（可选）」填的目录或 setup 脚本 → 包内运行时 `bin/` → `CBIN`/`CCP4_BIN`/`OASIS_CCP4_BIN` → `CCP4`/`CCP4_ROOT` → `/Applications/ccp4-*/bin`、`/opt/ccp4*/bin`、`/usr/local/ccp4*/bin` → `PATH`。

- 找到 `cparrot`：模块 2 默认勾选；找到 `cbuccaneer` **且** `refmac5`：模块 3 默认勾选。
- 都没找到：两个模块默认关闭、参数组置灰，提示"未检测到 CCP4（cparrot、cbuccaneer、refmac5），默认只做模块 1"。装上 CCP4 或设好 `CCP4`/`CBIN` 后点「重新检测」即可。
- 装在别处、又没进 `PATH`（双击启动的 GUI 读不到 `.bashrc` 里的 `CCP4`）：在「CCP4 位置」里填安装根目录，或直接把 `ccp4.setup` / `.bashrc` 的路径粘进去；解析出的根目录会写回该栏，检测和运行都用它。这一栏同时是运行目录（子进程按 `<root>/bin/<exe>` 启动，且运行环境会清掉继承的 `CCP4*`）。
- **但只装 CCP4 跑不了 P+。** CCP4 提供的是 Phaser 前端、Parrot、Buccaneer、REFMAC；P+ 的似然精修在包内的 `density/bin/oasis-density-python` 里（RESOLVE/CCTBX 派生，见下一节），CCP4 里没有这个程序。表单因此把「CCP4 位置」和「包内密度引擎（density/）」分开：
  - 选的是便携包 runtime（自带 `density/`）：两栏都指向它，直接跑；
  - 选的是外部 CCP4：P+ 仍要用一份 `density/`。表单会先找运行时自己的 `density/`，找不到就自动改用便携包 runtime 并在提示里写明原因；你也可以在「显示高级参数 → 运行资源 → 包内密度引擎」里填一份 `density/`，这样 CCP4 继续用你自己的，只有 P+ 借用那份树；
  - 两边都没有：点「开始」会直接报错，说明缺哪个文件、以及三条解决途径，而不是让任务死在 Phaser 之前。
- 模块 3 依赖模块 2（Buccaneer 从密度修饰后的图开始建）；关掉模块 2 会自动关掉并置灰模块 3。

只勾模块 1 就是精简版：输出 `02_sad/experimental.mtz`（相位/FOM/PPLUS）与 `observations.mtz`，**没有模型、没有 R-free**。3O2E 实测约 2 分钟、5763 反射、mean FOM 0.482、图相关 0.123（三段链是 0.549）。命令行等价写法是 `--lite` 或 `--no-dm-enabled`，详见 `docs/LITE_MODULE_20261112.md`。

注意：免安装包自带的运行时里已经含 Phaser 前端与 Parrot/Buccaneer/REFMAC，所以只要包完整，检测就会通过，不需要目标机器另装 CCP4；`docs/INSTALL_BUNDLE.md` 方式 C（纯源码）才需要本机 CCP4。

## 实测：本机 CCP4 的引擎 + 包内 `density/`（2026-11-13，胰岛素 2.1 Å）

同一条胰岛素任务（`ins_wave2.0_reso2.1A_sad.mtz` + `sg.pdb` + `3INS.pir`，S 边，λ=2.0 Å，`pplus_w8_weight`，Parrot，双手自动判别）跑了三遍：

| 跑法 | 引擎来自 | P+ 来自 | 用时 | 5 轮后 Rwork / Rfree |
| --- | --- | --- | --- | --- |
| 全便携包 | 便携包 `runtime` | 便携包 `runtime/density` | 18.1 min | 31.66 / 37.79 |
| 命令行混搭 | `/Applications/ccp4-9` | dev17 便携包的 `density/` | 11.1 min | 31.66 / 37.79 |
| 窗口（填 CCP4、密度引擎留空） | 便携包 `runtime`（自动切换） | 便携包 `runtime/density` | 11.9 min | 31.66 / 37.79 |

三次的判手都是 `input_hand` / `I 21 3`，指标逐位相同，因为 `cparrot`、`cbuccaneer`、`refmac5`、`phaser` 四个文件与便携包 `runtime/bin/` 里的 sha256 完全一致——那个便携包就是用这台机器的 CCP4 组装的。换成版本不同的 CCP4，R 会变，但链路照走。

命令行混搭的完整命令（`--density-runtime` 与 `--runtime` 分开给）：

```bash
./oasis sad-to-model --method pplus_bundled --dm-engine parrot --dm-cycles 40 \
  --feedback-cycles 2 --trust-wang-radius 8 --dm-solvent 0.5 --auto-hand \
  --input-columns --elements S --wavelength 2.0 --copies 1 \
  --sad <insulin sad>.mtz --heavy sg.pdb --sequence 3INS.pir \
  --runtime /Applications/ccp4-9 \
  --density-runtime "<便携包>/runtime/density" \
  --reference /Applications/ccp4-9/lib/data/reference_structures/reference-1tqw.mtz \
  --output runs/mixed
```

窗口里的行为（2026-11-13 起）：只填 `/Applications/ccp4-9`、密度引擎留空时，点「开始」**不会**换运行时，而是把方法改成「连续 SAD」并在红条里说明「本次改用连续 SAD + 密度修饰 + Buccaneer/REFMAC，只依赖 CCP4 自带程序」；填了密度引擎路径则保留 P+ 方法。`--density-runtime` 现在只在显式指定 `density/` 时有意义，不再用来“借树”。

## 方式 B：从维护者处拿现成构建（私有分发）

只带四个引擎的包压缩后很小，不需要分卷：

```bash
tar -czf OASIS4.3-portable-ccp4only.tar.gz OASIS4.3-portable-ccp4only-20261111
shasum -a 256 OASIS4.3-portable-ccp4only.tar.gz > OASIS4.3-portable-ccp4only.tar.gz.sha256
```

实测：1.4 GB 的目录压成 **436 MB**、45 370 个条目（53 秒），单个文件即可传。2026-11-11 这一版的 tarball 是 456 700 259 字节，SHA256 `f105ff8ec243b70a72117b3e3c9cee437fc1ebff72e7edab42f5dede55d40118`。接收方：

```bash
shasum -a 256 -c OASIS4.3-portable-ccp4only.tar.gz.sha256
tar -xzf OASIS4.3-portable-ccp4only.tar.gz
cd OASIS4.3-portable-ccp4only-20261111 && ./oasis doctor
```

- 带 RESOLVE/CCTBX 的 dev17 完整包（3.8 GB）压缩后仍可能超过 2 GB 的 GitHub Release 附件上限，那才需要 `split -b 1900m … part-` 分卷，接收方 `cat part-* | tar xz`。
- 渠道用内网、移动硬盘，或**私有**仓库的 Release；不要放到公开仓库。
- 包里带着 CCP4 的二进制，按 CCP4 许可，接收方也应已持有相应授权；接收机需为 macOS Intel。
- 接收方第一次运行：从网络下载的文件会被 macOS 打上隔离属性，而包里的二进制未签名（`codesign -dv` 报 `code object is not signed at all`），双击 `OASIS4.3.command` 可能被 Gatekeeper 拦下。在终端执行一次 `xattr -dr com.apple.quarantine OASIS4.3-portable-ccp4only-20261111` 即可；从终端直接 `./oasis doctor`、`./oasis sad-to-model` 不受影响（已实测）。
- `bundle.json` 与 `DOCTOR.json` 里的路径是相对包根的（`python/bin/python3` 等），换机器、换目录都不会失效。

### 方式 B 的包能复现标定结果吗？能不装 CCP4 吗？（2026-11-11 实测）

**能复现。** 把 tarball 解到新目录，用空环境（`env -i`，不继承 CCP4、Phenix、系统 Python 的任何变量）把十套标定的选中终点逐个重跑，每套参数照抄该终点 `pipeline_report.json` 里的 `config`：

```bash
env -i PATH=/usr/bin:/bin HOME=/tmp/oasis-clean-home LANG=C \
  ./oasis sad-to-model --method pplus_bundled --dm-engine parrot --no-auto-hand \
    --input-columns --feedback-cycles 3 --dm-solvent 0.50 --dm-mask-radius 3.0 \
    --dm-anisotropy --trace-resolution 2.0 --fcalc fc_all --copies 1 --rounds 5 \
    --cycles 10 --phase-cycles 100 --dm-cycles 40 …（其余照抄 config）
```

| 案例 | 终点 | 标定 R / R-free | 包内重跑 R / R-free | 坐标 |
| --- | --- | --- | --- | --- |
| 3O2E | fb3，溶剂 auto | 41.14 / 42.64 | 41.14 / 42.64 | 逐行相同 |
| 5D9D | fb3，0.50 | 25.43 / 27.51 | 25.43 / 27.51 | 逐行相同 |
| 5HHK | fb3，0.50 | 25.94 / 28.82 | 25.94 / 28.82 | 逐行相同 |
| 6HF7 | fb3，0.50 | 26.02 / 29.04 | 26.02 / 29.04 | 逐行相同 |
| 6ZPP | fb2，0.50 | 30.04 / 31.31 | 30.04 / 31.31 | 逐行相同 |
| 7E1D | fb3，0.50 | 36.65 / 42.33 | 36.65 / 42.33 | 逐行相同 |
| 7ENM | fb2，0.50 | 28.34 / 31.47 | 28.34 / 31.47 | 逐行相同 |
| 7TBS | fb3，0.50 | 25.17 / 28.35 | 25.17 / 28.35 | 逐行相同 |
| 7TFQ | fb2，0.50（`pplus_trust`） | 32.52 / 36.86 | 32.52 / 36.86 | 逐行相同 |
| 7TJ1 | fb2，0.50 | 23.92 / 31.08 | 23.92 / 31.08 | 逐行相同 |

十套的 R 与 R-free 都与标定一致到小数点后四位，`primary_model.cif` 的坐标逐行相同。文件 sha256 仍会不同：CIF 头部 `_pdbx_database_status.recvd_initial_deposition_date` 记的是写文件日期，比对坐标前先去掉 `_` 开头的头行。

复现时最容易踩的三个坑（都是参数没照抄，不是包的问题）：

1. **3O2E 的终点用 `auto` 溶剂，不是 0.50。** 单跑时 `--dm-solvent` 只收数字，必须走端点网格 `--feedback-endpoints 3 --dm-solvents auto`；写成 0.50 会落到另一条终点（47.01 / 50.89）。
2. **7TFQ 的终点是 `pplus_trust`**，要显式 `--dm-engine parrot`（本版起允许；`--dm-engine resolve` 仍会被拒绝，因为包里没有 Phenix）。它的 `--observation-prior` 用同一套 `pplus_bundled --no-build-model` 跑出的 `02_sad/observations.mtz` 就行。
3. **每套的 `--max-per-target` 不同**（多数 48，5HHK 13，6HF7 24）；6HF7 还要 `--copies 3 --elements Tb,Mg --wavelength 1.6486`。

**不用装 CCP4。** 三条独立证据：

- `env -i ./oasis doctor` 通过，报 `requires_separate_CCP4: false`、`requires_separate_Phenix: false`，随包引擎 Phaser / Parrot / Buccaneer / REFMAC。
- 启动脚本先 `unset CCP4 CLIBD CLIBD_MON CBIN CEXAM CINCL CLIB PHENIX PHENIX_ROOT LIBTBX_BUILD` 及全部 `DYLD_*`/`LD_*`，再设 `PHASEBUILD_RUNTIME=<包>/runtime`；机器上不需要任何 CCP4 目录存在。
- `otool -L` 检查：`runtime/bin/*` 的非系统依赖全是 `@loader_path/../lib/…`（库就在包里），`python/bin/python3` 只链接 `/usr/lib/libSystem`；`./oasis gui` 也能用随包 Python 3.12 + Tk 8.6 直接起窗口。

## 方式 C：只要代码

```bash
python -B package_oasis.py --output /new/path.zip
```

生成精简源码包（**不含** Buccaneer/REFMAC/Parrot 运行时）。在已装 CCP4 的机器上可直接从源码树运行；缺 Parrot 时会回退到内置引擎，并明确报出回退原因（见 [README](../README.md) 与 `docs/LINUX.md`）。

## 包内的自研模块是实验性的（默认关闭）

包里有三个自研替换件，**都不在默认路径上**，只用于研究对比；出结构请一律用 `pplus_bundled`：

| 模块 | 开关 | 状态 |
| --- | --- | --- |
| 自研 SAD 定相 `pplus_internal` | 窗口方法选项 / `--method pplus_internal` | 可用：输入相位图相关 0.279（Phaser 0.272） |
| 自研密度修饰 `internal_pplus.seeded_dm` | 命令行/环境变量 `OASIS_INTERNAL_DM_*` | 实验性：终点图相关 0.558（Parrot 0.688） |
| 自研自动追踪 `trace.py` | 未接入流水线 | 原型：SHELXE 级地图上也只长出 3–8 残基碎片 |

- 窗口里 `内部 P+` 标着"实验性，不用于正式解析"，方法下方直接写上述数字；密度修饰下拉框只有 Parrot/RESOLVE。
- 判断依据与复现命令：`docs/NO_CCP4_REPLACEMENT.md`、`docs/SHELXE_BASELINE.md` §3.3/§3.4。
- 为什么要留在包里：它们是"能不能不依赖 CCP4"这条线的实测记录，去掉就无法复现结论；但它们**不是**可以替代 Parrot/Buccaneer 的选项。

## 还有两个非许可的障碍

- 即使许可解决，含 RESOLVE/CCTBX 的完整包 3.8 GB，超过单个 Release 附件的 2 GB 上限，需要分卷与校验清单。
- Linux 目前还没有完整运行时；`docs/LINUX.md` 记录了在 Linux 上准备原生引擎的进展。

## 相关文件

- `packaging/build_local_runtime.py` —— 从本机 CCP4 生成私有 runtime（带 manifest 与许可文本，`--engines` 选引擎）
- `packaging/build_installation_free.py` —— 组装免安装包并跑干净环境自检
- `packaging/audit_linux_runtime.py` —— Linux 运行时审计（尚未提供完整包）
- `package_oasis.py` —— 精简源码包
- `src/phasebuild/module_detection.py` —— 窗口启动时的 CCP4 模块检测（找不到是正常返回值）
