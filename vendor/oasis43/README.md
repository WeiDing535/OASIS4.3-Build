# OASIS4.3-GUI-0.1

[English](README_EN.md) | 中文

项目主页：[github.com/WeiDing535/OASIS4.3-GUI](https://github.com/WeiDing535/OASIS4.3-GUI)

OASIS4.3-GUI-0.1 是面向蛋白质 X 射线 SAD 相位计算的独立研究软件。它把 OASIS/P+ 方法的现代化 Python 实现、连续 SAD 相位似然和不依赖外部 Parrot 的密度修饰放在同一套命令行和桌面界面中。程序读取实验 MTZ、蛋白质序列与重原子坐标，自动估算 ASU 组成和溶剂含量，输出含 PHIB、FOM 与 Hendrickson–Lattman 系数的 MTZ。

本项目由中国科学院物理研究所丁玮（DingWei）维护，采用 LGPL-2.1-or-later。算法仍处于科研验证阶段；FOM 是模型置信度，不能单独证明结构已解出。

## 一键安装

需要 64 位 Python 3.12、网络连接和约 1 GB 可用空间。在源码目录运行：

```bash
python3.12 install.py
```

Windows 可使用：

```powershell
py -3.12 install.py
```

安装器会在用户目录创建隔离环境、安装 GUI 依赖，并从固定的 HTTPS 地址下载一次经过 SHA-256 校验的通用参考 MTZ。若需要离线安装，可先取得文档所列的原始 1TQW 参考文件，然后运行：

```bash
python3.12 install.py --reference-source /path/to/reference-1tqw.mtz
```

默认安装位置是 `~/.oasis-sad`。使用 `--prefix` 可指定其他位置，`--no-launch` 安装后不启动界面，`--cli-only` 只安装命令行程序。

## Conda 环境

本仓库提供可复现的 Conda 环境文件；这是源码安装方式，尚未发布到 conda-forge：

```bash
conda env create -f environment.yml
conda activate oasis43-gui
oasis-py reference-fetch
oasis-gui
```

## GUI 输入

GUI 当前重点支持蛋白质 SAD。必须提供：

- 含异常振幅和标准差的 MTZ；
- FASTA、PIR 或纯文本蛋白质序列文件；
- 整个序列文件在 ASU 中的 copies 数；
- 与实验数据同一晶胞、空间群和原点的重原子 PDB/mmCIF；
- 与元素及实验波长对应的正数 f″。

程序根据序列、copies、晶胞和空间群自动计算 Matthews 溶剂含量。每个输入项后的 `?` 提供具体说明。SIR、DMR 和更多实验选项可通过 `oasis-py --help` 查看。

## 测试

普通测试完全离线，使用合成参考数据：

```bash
python -m pip install -e '.[test,gui]' -c constraints.txt
pytest -m 'not reference_data'
```

Linux、Windows、Intel macOS 和 Apple Silicon macOS 的 GitHub Actions 已全部通过，详见[发布检查记录](docs/RELEASE_CHECKS.md)。真实参考数据集成测试需先运行 `oasis-py reference-fetch`，再显式执行 `pytest -m reference_data`。

## 引用与方法来源

软件引用信息见 [CITATION.cff](CITATION.cff)。OASIS/IPCAS 方法背景包括：Ding W. et al., *IPCAS: a direct-method-based pipeline from phasing to model building and refinement for macromolecular structure determination*, J. Appl. Cryst. 53, 253–261 (2020), DOI: [10.1107/S1600576719015115](https://doi.org/10.1107/S1600576719015115)。

Clipper 派生的数值内核、对应源码和额外许可说明见 [第三方说明](docs/THIRD_PARTY.md)。参考 MTZ 的来源、校验值和再分发处理见 [数据来源](docs/DATA_SOURCES.md)。

## 联系方式

丁玮（DingWei）：dingwei@iphy.ac.cn
