"""Post-build reference evaluation and fixed-stage comparison; never a builder input."""
import json
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent / 'src'))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from phasebuild.data import fingerprint
from phasebuild.evaluate import evaluate
from phasebuild.sidechain_evaluate import evaluate_sidechains

ROOT = Path(__file__).resolve().parent
OUT = ROOT / 'runs/comparison_v0p8'
REFERENCE = Path('/Volumes/For_data/IPCAS-test/rebuilt/INS/ins_reference.pdb')


def read(path):
    return json.loads(Path(path).read_text())


def row(label, cycles, report_path):
    report_path = Path(report_path)
    report = read(report_path)
    model = report_path.parent / 'model.cif'
    backbone = evaluate(model, REFERENCE)
    sidechains = evaluate_sidechains(model, REFERENCE)
    # Apply the same current chemistry audit to every comparator, including Pro.
    import gemmi
    from phasebuild.atom_validation import atom_geometry_report
    from phasebuild.rotamers import load_library
    geometry = atom_geometry_report(gemmi.read_structure(str(model)),
                                   load_library(ROOT / 'models/rotamers_1tqw_v08.json'))
    assert report['status'] == 'success' and report['free_flags_preserved']
    assert report['reflections']['free_value'] == 0
    assert report['statistics']['reflections_used_free'] == 299
    return dict(label=label, cycles=cycles, model=fingerprint(model),
                refinement_report=fingerprint(report_path),
                r_work=report['statistics']['final']['r_work'],
                r_free=report['statistics']['final']['r_free'],
                backbone=backbone, sidechains=sidechains, geometry=geometry)


def main():
    OUT.mkdir(exist_ok=True)
    baseline = ROOT / 'runs/ins_v025_phasebuild_v07_refined/refinement_report.json'
    rows = [row('PhaseBuild 0.7 baseline', 20, baseline)]
    curves = {}
    build_audits = []
    for folder, label in [('ins_v038_v081_primary', '0.8.1 original-map restart'),
                          ('ins_v039_v081_updated_only', '0.8.1 updated-map only')]:
        stages = [rows[0]]
        for cycles, part in [(30, 'cycles_30'), (40, 'continuation/cycles_40'),
                             (50, 'continuation/cycles_50')]:
            root = ROOT / 'runs' / folder / part
            item = row(label, cycles, root / 'refine/refinement_report.json')
            rows.append(item)
            stages.append(item)
            report = read(root / 'build/build_report.json')
            assert not report['source_changed_during_run']
            assert all(fingerprint(source['path']) == source for source in report['sources'])
            build_audits.append(dict(report=fingerprint(root / 'build/build_report.json'),
                                     input_map=report['input'], extension=report['extension'],
                                     seconds=report['seconds']))
        assert not read(ROOT / 'runs' / folder / 'continuation/experiment_report.json')['source_changed_during_run']
        curves[label] = stages
    old = read(ROOT / 'runs/ins_v031_updated_map_50/experiment_report.json')
    curves['0.7 updated-map only'] = old['stages']
    rows.append(row('PhaseBuild 0.7 updated-map only', 50,
                    old['final_refinement_report']['path']))
    bc = read(ROOT / 'runs/ins_v028_buccaneer_pipeline/pipeline_report.json')
    curves['Buccaneer correlation pipeline'] = [dict(s, cycles=10*s['round']) for s in bc['stages']]
    rows.append(row('Buccaneer correlation pipeline', 50, bc['stages'][-1]['refinement_report']['path']))
    protocols = [read(r['refinement_report']['path']) for r in rows]
    assert len({r['reflections']['input']['sha256'] for r in protocols}) == 1
    assert len({tuple(r['keywords']) for r in protocols}) == 1
    result = dict(rows=rows, build_audits=build_audits,
                  map_diagnosis=read(ROOT / 'runs/diagnosis_v08/map_diagnosis.json'),
                  reference=fingerprint(REFERENCE),
                  scope='Development set, model-derived input phases. Fixed final stage; no Rfree selection. '
                        'Main arm restarts from original supplied map at cumulative cycle 20; '
                        'Buccaneer subsequently uses its own updated maps. Not identical map schedules. '
                        'Reference access is post-build evaluation only.')
    (OUT / 'comparison.json').write_text(json.dumps(result, indent=2, allow_nan=False)+'\n')
    fig, ax = plt.subplots(figsize=(8.6, 5.2))
    for label, stages in curves.items():
        x = [s.get('cycles', s.get('cumulative_cycles')) for s in stages]
        ax.plot(x, [s['r_free']*100 for s in stages], 'o-', label=label)
    ax.set(xlabel='Cumulative REFMAC cycles', ylabel='Rfree (%)', xticks=[10,20,30,40,50],
           title='INS: fixed-stage comparison (model-derived phases)')
    ax.grid(alpha=.25)
    ax.legend(fontsize=9)
    fig.tight_layout()
    for suffix in ['png','svg']:
        fig.savefig(str(OUT / ('rfree_cycles.'+suffix)), dpi=170)
    plt.close(fig)
    table = ['| 流程 | 累计精修轮数 | Rwork | Rfree | 构建残基 | N/CA/C 均 <1 Å | 正确侧链原子/已建 |',
             '|---|---:|---:|---:|---:|---:|---:|']
    for r in rows:
        b,s = r['backbone'],r['sidechains']
        table.append('| {} | {} | {:.2f}% | {:.2f}% | {} | {} | {}/{} |'.format(
            r['label'],r['cycles'],r['r_work']*100,r['r_free']*100,b['built_residues'],
            b['n_ca_c_correct_residues'],s['sidechain_atoms_within_cutoff'],s['sidechain_atoms_built']))
    report = '''# PhaseBuild 0.8.1：末端延伸与相位图对照

新流程在 INS 开发数据上将固定 50 轮的 Rfree 降至 **27.05%**，Buccaneer 对照为 **26.91%**，差 0.14 个百分点。新模型 49/49 残基严格主链正确，覆盖参考 51 残基的 96.1%；不能据此宣称在独立数据上达到或超过 Buccaneer。

## 结果

TABLE

严格正确指同一晶体对称像下 N、CA、C 分别距参考 <1 Å；侧链采用相同阈值并允许化学等价原子名互换。全部残基按一对一 CA 匹配，不能把匹配后 RMSD 当作全模型准确性。主流程最终 49 个残基类型均与参考一致，覆盖 A1–A21、B2–B29，仅缺 B1、B30。两条链分别与同一参考链、同一对称像对应，匹配残基顺序递增。当前对照为完整 correlation 模式 Buccaneer 构建流程；其他历史模式见 comparison_v0p7，未只拿较弱的默认模式比较。

![固定阶段 Rfree 曲线](rfree_cycles.png)

## 差距在哪里，解决了什么

1. 旧模型严格正确部分集中于 A1–A21 与 B4–B19，缺少 B 链两端。旧局部循环不能生成新的末端主链；增加精修轮数仍停在 38 残基。现在加入序列锚定的双向延伸、多路径扭转角搜索、末端一个残基的修复、连续优化和密度/几何/对称碰撞筛选。内部锚点在延伸阶段固定，REFMAC 随后优化全模型。Buccaneer 官方流程也包含 find/grow/join/link/sequence/correct 等多个步骤，迭代精修不等于迭代补全；见[官方说明](https://www.ccp4.ac.uk/html/cbuccaneer.html)。本程序没有移植或反编译其代码。
2. 不完整模型自身更新图会削弱尚未建区域的信息。本例参考 B20–B30 位置仅用于事后诊断采样，原始提供图与旧模型更新图的平均主链原子密度为 3.10 与 1.97，profile 为 3.49 与 2.89。主流程在 20 轮起点重新使用一次原始提供图，随后两次仅用自身更新图，最终达到 49 残基和 27.05%。全程自身更新图分支为 42 残基、35.65%。这是强烈的流程依赖证据；延伸算法与输入相位信息共同贡献了改善。
3. Gly 正 φ 的合成验证暴露末端 O 方向采样不足，现已补充。Pro 的 monomer 字典使用另一种键型列名，旧解析器漏读环连接；现在兼容两种格式并在键图缺失时报错。使用重新生成的同源 1TQW 构象库，不使用 INS 参考构象。修复字典不意味着所有 Pro 侧链已成功补全。

## 协议、质量与可复现性

所有表内精修使用同一 native F/SIGF、原始 FreeR_flag=0、299 个自由集反射以及同一 REFMAC 关键词，最终阶段固定，不按最低 Rfree 挑轮次。6148 条可用观测中 REFMAC 统计使用 6146 条（5847 work +299 free）；另两条未计入的原因尚未完全定位，但各组统计一致。保留 HKL 的自由集标记核验未改变。

主流程末次精修的键长 RMSD 为 0.0112 Å、键角 RMSD 为 1.794°；当前简化几何审计没有发现 <1.8 Å 的严重非键重原子接触。这不等同于完整 MolProbity/Ramachandran 验证。49 残基中仍有 19 个已建侧链原子未达到参考 1 Å 阈值；末端与侧链仍需改进。全部比较对象用同一修复后的化学库重新审计，详见 comparison.json。

51 项测试通过，结果见 ../verification_v0p8_tests.xml。每个构建保存实际输入地图哈希、源码哈希、序列锚定和新增事件；汇总核对当前源码与运行记录一致。运行目录 ins_v038_v081_primary 与 ins_v039_v081_updated_only 的 pipeline_report.json 均为 success。0.8.0 的探索运行 ins_v032–ins_v037 已保留并标注，不能与最终 0.8.1 混成一条轨迹。可在项目根目录运行 `ccp4-python summarize_extension_experiment.py` 重建本报告和图。

输入相位来自现成 REFMAC 模型，可能已有自由集/模型信息；本数据也用于开发。构建器不读取目标参考或 Buccaneer 坐标，参考仅由事后评估器访问。主流程与 Buccaneer 的地图使用时间表不同，累计 50 轮只是精修次数相同，不是相同总计算预算。0.14 个百分点的 Rfree 差别不足以作优劣统计判断。

## 下一步

优先在有独立实验相位的多蛋白测试集上检验“保留原始相位约束”的迭代策略，再做缺口连接、全图残基发现和侧链优化。相位约束应来自 OASIS 的实验相位/质量信息并明确自由集隔离，而不是反复依赖某个参考模型的精修地图。本轮已验证的是补全能力与地图信息保留的重要性，尚未验证 SAD 到完整模型的泛化性能。
'''.replace('TABLE','\n'.join(table))
    (OUT / 'REPORT.md').write_text(report)
    print('\n'.join(table))


if __name__ == '__main__':
    main()
