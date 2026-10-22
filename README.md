# OASIS4.3-Build

[English](README_EN.md) | 中文

从已有晶体学相位图构建蛋白质模型的研究软件，包含原生 PhaseBuild 算法、本机浏览器 GUI，以及外部 Buccaneer / REFMAC 集成。

软件内核和界面沿用 **PhaseBuild 0.8.3** 名称。本仓库首个公开源码修订为 **2026-10-22**。旧版 OASIS4.3-GUI 定相软件保留在[独立仓库](https://github.com/WeiDing535/OASIS4.3-GUI)，没有被替换。

## 完整包目标

Intel Mac 完整本机试用包已提供组装工具，包含 Python 和全部 GUI 模式所需的外部运行组件，用户无需另装 CCP4。该运行时目前只用于本机评估，尚未公开发布二进制安装包；GitHub 上仍是源码。详见 [完整包说明](docs/STANDALONE.md)。

## 本机启动

需要 Python ≥3.9、Gemmi、NumPy、SciPy。建模精修集成另外需要用户已有的 CCP4 安装及其参考资源。本仓库不附带 Python、CCP4 二进制或完整便携运行包。

macOS 已安装 CCP4 到 `/Applications/ccp4-9` 时，双击 **PhaseBuild.command**，或在仓库根目录运行：

```bash
./PhaseBuild.command
```

其他安装位置可指定：

```bash
PHASEBUILD_PYTHON=/path/to/ccp4/bin/ccp4-python \
PHASEBUILD_RUNTIME=/path/to/ccp4 ./PhaseBuild.command
```

也可直接使用已有科学 Python 环境：

```bash
python -B run_gui.py --project .
```

运行后打开终端显示的完整本机链接。关闭网页不会停止任务；通过界面停止计算，在终端按 Ctrl+C 关闭服务。输出目录必须是新目录。

## 工作方式

| 模式 | 说明 |
|---|---|
| 从相位图建模 | 原生主链搜索、序列归属、侧链和可选精修 |
| 已有模型延伸 | 从已有序列锚点继续构建 |
| 单独精修 | 调用外部 REFMAC |
| 独立 R 因子 | 固定坐标评估，工作集拟合尺度与体溶剂，无需 REFMAC |
| 无初始模型＋原 HL 精修 | 外部 Buccaneer 从已有相位构建，REFMAC 每轮保留固定实验 HL |
| Buccaneer 集成＋原 HL 精修 | 从已有种子模型继续外部构建和精修 |

集成模式需提供相位图、native MTZ、序列及原始实验 HL。原自由集数值必须与输入一致；先点击“检查输入”，再运行。详细说明见 [GUI_GUIDE.md](GUI_GUIDE.md)。

## 验证范围

发布前在现有 macOS Intel / CCP4 Python 3.9.18 环境运行了 **186 项回归测试**，全部通过，并检查 GUI 加载和本机 API。外部引擎合同测试包含模拟进程；这些结果不代表新入口已通过真实数据结构求解验证，也不代表其他操作系统已经验证。

原生算法仍为研究原型；外部集成不能作为原生算法达到相同性能的证据。此软件消费已有相位图，不包含完整上游 SAD 定相引擎。自由集标记保持不变不能消除上游地图可能已有的信息泄漏。

```bash
python -B -m pytest -q -p no:cacheprovider tests
```

测试需要已有 pytest。GitHub Actions 配置另行验证源码回归，不安装外部科学引擎。

## 源码与许可

维护者：丁玮（DingWei），中国科学院物理研究所。沿用项目的 LGPL-2.1-or-later 许可，见 [LICENSE](LICENSE)。第三方依赖、CCP4 引擎和资源的许可独立适用；本仓库不授予其再分发权。训练模板来自开发副本中的无关参考结构，见 [来源说明](docs/PROVENANCE.md)。

本仓库只包含软件、测试和小型训练模板；不上传本机任务索引、临时测试目录、目标衍射数据、目标结构或运行日志。
