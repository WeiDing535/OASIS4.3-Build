"""Package the ten SAD cases into two self-contained bundles:

  bundles/refmac_set/     -- data for developing a REFMAC-like restrained refiner
  bundles/buccaneer_set/  -- data for developing a Buccaneer-like model builder

Every case directory is standalone: native data, sequence, models, and a
meta.json.  A manifest.json carries sizes and sha256 for integrity checks.

Frame conventions (see work/frame_fit.py):
  * the campaign DM map and the deposited model share the same hand, but differ
    by an origin shift, a space-group operator and a lattice translation;
  * the Buccaneer bundle keeps truth/built in the map frame (they are directly
    comparable);
  * the REFMAC bundle moves the campaign build back into the deposited frame so
    that RMSD against truth_model.pdb is meaningful.

Run:  ./.conda-unified/bin/python work/make_bundles.py
"""
import hashlib, json, shutil
from pathlib import Path
import gemmi
import numpy as np

ROOT = Path(__file__).resolve().parent.parent
INV = {c['pdb']: c for c in json.loads((ROOT / 'work/tencase_inventory.json').read_text())}
FRAMES = {f['pdb']: f for f in json.loads((ROOT / 'work/tencase_frames.json').read_text())}
OUT = ROOT / 'bundles'

REFMAC = ['5HHK', '6ZPP', '7ENM', '7TBS', '7TJ1']
BUCCANEER = ['3O2E', '5D9D', '6HF7', '7E1D', '7TFQ']


def strip_model(path):
    st = gemmi.read_structure(str(path))
    st.remove_alternative_conformations()
    st.remove_hydrogens()
    try:
        st.remove_ligands_and_waters()
    except RuntimeError:
        pass
    st.remove_empty_chains()
    return st


def write_pdb(st, path):
    st.write_pdb(str(path))
    return path


def perturb(st, sigma, seed):
    st2 = st.clone()
    rng = np.random.default_rng(seed)
    for ch in st2[0]:
        for r in ch:
            for at in r:
                at.pos = gemmi.Position(*(np.array([at.pos.x, at.pos.y, at.pos.z])
                                          + rng.normal(0.0, sigma, 3)))
    return st2


def map_to_deposited(st, frame):
    """Campaign build (map frame) -> deposited frame, per chain."""
    cell = st.cell
    ops = list(gemmi.SpaceGroup(st.spacegroup_hm).operations())
    shift = frame['shift_fractional']
    by_chain = {row['built_chain']: row for row in frame.get('chain_map') or []}
    out = st.clone()
    for ch in out[0]:
        row = by_chain.get(ch.name)
        if row and row.get('symmetry_op') is not None:
            op = ops[row['symmetry_op']].inverse()
            lat = row.get('lattice') or [0, 0, 0]
        else:
            op = ops[frame['symmetry_op']].inverse()
            lat = [0, 0, 0]
        for r in ch:
            for at in r:
                f = cell.fractionalize(at.pos)
                xyz = [f.x - shift[0] - lat[0],
                       f.y - shift[1] - lat[1],
                       f.z - shift[2] - lat[2]]
                xyz = op.apply_to_xyz(xyz)
                at.pos = cell.orthogonalize(gemmi.Fractional(*xyz))
    return out


def sha256(path):
    h = hashlib.sha256()
    with open(path, 'rb') as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b''):
            h.update(chunk)
    return h.hexdigest()


FREE_WORK, FREE_TEST = 0, 1


def clean_mtz(src, dest, needed, free_column=None):
    """Rewrite an MTZ for use outside OASIS.

    Two defects of the campaign files are repaired here:

      * the PDB-derived ``-sf.mtz`` carries NaN in F/SIGF (2-5% of rows) and DM
        passes some of them through, so a plain reader propagates NaN into R;
      * OASIS writes ``FreeR_flag`` with 0 marking the test set and every other
        value the work set, which is the opposite of the REFMAC/Phenix default.
        Stock ``refmac5`` therefore treats 95% of the data as free unless it is
        told ``FREE 0``.  The same test set is re-encoded as 0 = work, 1 = free.
    """
    m = gemmi.read_mtz_file(str(src))
    labels = [c.label for c in m.columns]
    data = np.array(m, copy=False)
    ok = np.ones(len(data), bool)
    for lab in needed:
        if lab in labels:
            ok &= np.isfinite(data[:, labels.index(lab)])
    sub = data[ok]
    out = gemmi.Mtz()
    out.cell = m.cell
    out.spacegroup = m.spacegroup
    out.add_dataset('HKL')
    columns = []
    for i, lab in enumerate(labels):
        col = sub[:, i]
        if lab == free_column:
            col = np.where(col == 0, FREE_TEST, FREE_WORK)
            out.add_column(lab, 'I')
        else:
            out.add_column(lab, m.columns[i].type)
        columns.append(col)
    out.set_data(np.stack(columns, axis=1).astype(np.float32))
    out.write_to_file(str(dest))
    info = dict(rows_in=int(len(data)), rows_out=int(len(sub)),
                dropped_nonfinite=int((~ok).sum()))
    if free_column:
        flags = columns[labels.index(free_column)]
        info.update(free_test_rows=int((flags == FREE_TEST).sum()),
                    free_fraction=round(float((flags == FREE_TEST).mean()), 4),
                    flag_convention='0=work, 1=free (REFMAC/Phenix default)')
    return info


def case_meta(c, frame):
    ts = c.get('truth_stats') or {}
    cfg = c.get('config') or {}
    met = c.get('metric') or {}
    return dict(
        pdb=c['pdb'], resolution=c['data']['dmin'], spacegroup=c['data']['spacegroup'],
        cell=c['data']['cell'], n_reflections=c['data']['nref'],
        asu_copies=cfg.get('copies'), chains=ts.get('chains'), residues=ts.get('residues'),
        ca=ts.get('ca'), mass_da=ts.get('mass_da'), z=ts.get('z'),
        solvent_fraction=ts.get('solvent'), heavy_elements=cfg.get('elements'),
        wavelength=cfg.get('wavelength'),
        campaign=dict(r_work=met.get('r_work'), r_free=met.get('r_free'),
                      method=cfg.get('method'), dm_solvent=cfg.get('dm_solvent'),
                      dm_cycles=cfg.get('dm_cycles'), phase_cycles=cfg.get('phase_cycles'),
                      feedback_cycles=cfg.get('feedback_cycles'),
                      building_recipe=cfg.get('building_recipe')),
        map_frame=dict(hand=frame['hand'], density_cc=frame['cc'],
                       shift_fractional=frame['shift_fractional'],
                       symmetry_op=frame['symmetry_op'],
                       symmetry_triplet=frame['symmetry_triplet'],
                       chain_map=frame.get('chain_map')),
        build_quality=frame.get('build_quality_map_frame'),
    )


def refmac_case(c, dest):
    dest.mkdir(parents=True, exist_ok=True)
    frame = FRAMES[c['pdb']]
    cleaning = {'native.mtz': clean_mtz(c['native'], dest / 'native.mtz',
                                        ['F', 'SIGF'], free_column='FreeR_flag')}
    shutil.copy2(c['fasta'], dest / 'sequence.fasta')
    shutil.copy2(c['truth'], dest / 'truth_model.pdb')
    st = strip_model(c['truth'])
    write_pdb(st, dest / 'start_model.pdb')
    write_pdb(perturb(st, 0.4, seed=sum(ord(ch) for ch in c['pdb'])),
              dest / 'start_model_perturbed.pdb')
    built = map_to_deposited(strip_model(c['model']), frame)
    write_pdb(built, dest / 'start_model_sad_built.pdb')
    meta = case_meta(c, frame)
    meta['data_cleaning'] = cleaning
    meta['starts'] = {
        'start_model.pdb': 'deposited model, waters/ligands removed (perfect start)',
        'start_model_perturbed.pdb': 'same atoms, 0.4 A Gaussian coordinate noise',
        'start_model_sad_built.pdb': 'campaign SAD build, moved back to the deposited frame',
    }
    (dest / 'meta.json').write_text(json.dumps(meta, indent=1, ensure_ascii=False))


def buccaneer_case(c, dest):
    dest.mkdir(parents=True, exist_ok=True)
    frame = FRAMES[c['pdb']]
    cleaning = {
        'map_phases.mtz': clean_mtz(c['dm'], dest / 'map_phases.mtz',
                                    ['FWT', 'PHWT', 'PHIB', 'FOM']),
        'experimental.mtz': clean_mtz(c['sad'], dest / 'experimental.mtz',
                                      ['F(+)', 'SIGF(+)', 'F(-)', 'SIGF(-)']),
        'native.mtz': clean_mtz(c['native'], dest / 'native.mtz',
                                ['F', 'SIGF'], free_column='FreeR_flag'),
    }
    shutil.copy2(c['fasta'], dest / 'sequence.fasta')
    truth_map = Path(frame['truth_in_map_frame'])
    write_pdb(strip_model(truth_map), dest / 'truth_model.pdb')
    write_pdb(strip_model(c['model']), dest / 'built_reference.pdb')
    if c.get('heavy'):
        shutil.copy2(c['heavy'], dest / 'heavy_sites.pdb')
    meta = case_meta(c, frame)
    meta['data_cleaning'] = cleaning
    meta['files_note'] = {
        'map_phases.mtz': 'density-modified map coefficients: FWT/PHWT (map), PHIB/FOM (phase+FOM)',
        'experimental.mtz': 'SAD amplitudes F(+)/F(-) with sigmas',
        'native.mtz': 'merged amplitudes + FreeR_flag (0=work, 1=free), for R/R-free scoring',
        'truth_model.pdb': 'deposited model in the frame of map_phases.mtz',
        'built_reference.pdb': 'campaign build in the same frame, with its R/R-free in meta.json',
    }
    (dest / 'meta.json').write_text(json.dumps(meta, indent=1, ensure_ascii=False))


REFMAC_README = """# REFMAC 类程序开发数据集（5 例）

用途：开发/对标一个 REFMAC 级别的**约束精修**程序（最大似然靶函数、几何约束、
溶剂/整体 B、free-R 监控）。数据来自 OASIS 十案例 SAD 战役，已裁成自包含单元，
可直接拷到另一台机器。

## 每个案例目录

| 文件 | 内容 |
| --- | --- |
| `native.mtz` | 合并振幅 `F SIGF` + `FreeR_flag`（**0=work，1=free**，约 5% free） |
| `sequence.fasta` | 序列 |
| `truth_model.pdb` | 沉积结构（含配体/水），做 RMSD、几何、水挑选的参考 |
| `start_model.pdb` | 沉积结构去掉水/配体/altloc，**完美起点**，测纯精修收敛 |
| `start_model_perturbed.pdb` | 同上原子，坐标加 0.4 Å 高斯噪声，测收敛半径 |
| `start_model_sad_built.pdb` | 战役 SAD 建模结果，已变换回沉积坐标系，**真实起点** |
| `meta.json` | 分辨率、空间群、晶胞、ASU 拷贝数、溶剂分数、战役 R/R-free、坐标系变换 |

## 为什么是这 5 例

| 案例 | 分辨率 | 空间群 | 拷贝 | 残基 | 溶剂 | 战役 R/R-free | 起点质量（CA，Å） | 作用 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 5HHK | 1.40 | P1 21 1 | 2 | 1407 | 0.36 | 25.9 / 28.8 | 0.23 | 最高分辨率，几何/ADP 精修与验证；双拷贝可测 NCS |
| 6ZPP | 1.50 | P1 21 1 | 1 | 337 | 0.42 | 30.0 / 31.3 | 0.21 | 高分辨率单链，快速迭代 |
| 7ENM | 1.70 | P21 21 21 | 1 | 251 | 0.30 | 28.3 / 31.5 | 0.07 | 小体系、低溶剂，最容易收敛的基线 |
| 7TBS | 1.96 | I222 | 1 | 336 | 0.51 | 25.2 / 28.4 | 0.08 | 高溶剂、I 心格子，测溶剂模型与整体 B |
| 7TJ1 | 2.10 | P21 21 21 | 3 | 595 | 0.50 | 23.9 / 31.1 | 0.14 | 最低分辨率 + 3 拷贝，测低分辨率几何权重与 NCS |

覆盖 1.40–2.10 Å、4 种空间群、1–3 拷贝、251–1407 残基、溶剂 0.30–0.51，
并且同时有“完美起点 / 扰动起点 / 真实 SAD 起点”三种输入。
（起点质量 = 战役 SAD 建模结果与沉积结构的 CA 最近邻距离中位数。）

## 建议的验收口径

1. `start_model.pdb` 精修后 R-free 不应劣于战役值（表内 R-free）。
2. `start_model_perturbed.pdb`（0.4 Å）应回到与 `truth_model.pdb` 的
   CA RMSD < 0.3 Å，且 R-free 与上一条相差 < 0.5 个百分点。
3. `start_model_sad_built.pdb` 精修后应达到战役 R/R-free（±0.5 点），
   并报告与 `truth_model.pdb` 的 CA RMSD / 主链几何。
4. 键长/键角 RMSD 应与 REFMAC 同量级（键长 ~0.01 Å、键角 ~1.5°）。

## 数据清洗（与战役文件不同之处）

- 源 `-sf.mtz` 里有一部分反射的 `F`/`SIGF` 是 NaN，DM 也会带出少量 NaN；
  本包已**删除 F 或 SIGF 非有限的整行**，逐例删除条数记在 `meta.json` 的
  `data_cleaning`。因此 R 的分母比战役略小（多数 <1%，5HHK/7ENM 约 3–8%），
  战役 R/R-free 仍是 REFMAC 在原始行集上算的，两者相差在 0.01 量级。
- 战役文件的 `FreeR_flag` 是 **0 表示 test set**，需要 REFMAC 加 `FREE 0`
  才正确；本包已重新编码为 REFMAC/Phenix 默认的 **0=work、1=free**，
  test set 的反射集合完全不变（`meta.json.data_cleaning.free_test_rows`），
  因此 stock `refmac5` / Phenix 不加任何关键字即可。

## 坐标系说明

`start_model_sad_built.pdb` 已用“密度互相关求原点平移 + 空间群算符/格点平移”
变换回沉积坐标系，因此与 `truth_model.pdb` 可直接算 RMSD。变换参数记在
`meta.json` 的 `map_frame` 里。
"""

BUCCANEER_README = """# Buccaneer 类程序开发数据集（5 例）

用途：开发/对标一个 Buccaneer 级别的**自动建模**程序（从实验密度图出发，
链追踪、序列对接、主链/侧链构建、迭代重建）。数据自包含。

## 每个案例目录

| 文件 | 内容 |
| --- | --- |
| `map_phases.mtz` | 密度修饰后的图系数：`FWT/PHWT`（图）、`PHIB/FOM`（相位+品质） |
| `experimental.mtz` | SAD 振幅 `F(+) SIGF(+) F(-) SIGF(-)`，可自行重算图，验证“不依赖外部程序” |
| `native.mtz` | 合并振幅 + `FreeR_flag`（**0=work，1=free**），建模后用 R/R-free 打分 |
| `sequence.fasta` | 序列 |
| `truth_model.pdb` | 沉积结构，**已变换到 `map_phases.mtz` 的坐标系**，用于覆盖率/CA RMSD |
| `built_reference.pdb` | 战役建模结果（同一坐标系），作为回归基线 |
| `heavy_sites.pdb` | 战役精修后的重原子位点（部分案例无） |
| `meta.json` | 分辨率、空间群、晶胞、拷贝数、溶剂分数、战役 R/R-free、坐标系参数、建模型质量 |

## 为什么是这 5 例

| 案例 | 分辨率 | 空间群 | 拷贝 | 残基 | 溶剂 | 战役 R/R-free | 建出 CA/真值 | 作用 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 5D9D | 1.70 | P21 21 21 | 1 | 469 | 0.46 | 25.4 / 27.5 | 307/306，98% | 顺利基线（Hg 单波长） |
| 7TFQ | 1.75 | P21 21 2 | 1 | 471 | 0.48 | 32.5 / 36.9 | 288/284，94% | 图中等，R 偏差大，测迭代重建 |
| 3O2E | 1.95 | P41 21 2 | 1 | 146 | 0.34 | 41.1 / 42.6 | 78/86，94% | 差图端：密度 CC 最低（0.53） |
| 6HF7 | 1.96 | P4 21 2 | 3 | 1122 | 0.55 | 26.0 / 29.0 | 602/573，91% | 高溶剂 + 3 拷贝 + Tb/Mg，测多拷贝 |
| 7E1D | 2.00 | P21 21 21 | 2 | 316 | 0.40 | 36.7 / 42.3 | 167/189，86% | 覆盖最不全，测不完整图下的建模 |

“建出 CA/真值”是战役模型与沉积结构的 CA 最近邻 <1 Å 的比例（见 `meta.json`
的 `build_quality`）。三个能建全的（5D9D、6HF7、7TFQ）+ 一个差图（3O2E）+
一个覆盖不全（7E1D），可做成功率与“从坏图里救回多少”的对照。

## 建议的验收口径

1. **覆盖率**：建出的 CA 数与 `truth_model.pdb` 的 CA 数之比；
2. **正确率**：建出 CA 与真值 CA 的最近邻距离 < 1 Å / < 2 Å 的比例；
3. **序列对接**：按残基编号对齐后，序列一致残基比例；
4. **R/R-free**：建模结果（或简单精修后）对 `native.mtz` 的 R/R-free，
   与 `built_reference.pdb` 的战役值比较；
5. **手性/原点**：`truth_model.pdb` 已按 `map_phases.mtz` 的坐标系放好，
   不要用沉积坐标直接比。

## 数据清洗（与战役文件不同之处）

- `map_phases.mtz`、`experimental.mtz`、`native.mtz` 里 `F`/`SIGF`/`FWT`/`PHWT`
  非有限的整行已删除（逐例条数见 `meta.json` 的 `data_cleaning`）。不删的话，
  朴素的读入会把 NaN 带进密度图和 R，源 `-sf.mtz` 本身就有 2–5% 这样的行。
- `native.mtz` 的 `FreeR_flag` 已重新编码为 REFMAC/Phenix 默认的
  **0=work、1=free**（战役文件是 0=test），test set 集合不变。

## 坐标系说明

SAD/DM 图与沉积结构**同一手性**（同手性密度 CC 0.53–0.73，反手 ≈0.03），
但原点差一个平移，并且 ASU 内的拷贝可能对应不同的空间群算符/相邻晶胞。
`truth_model.pdb` 已按“FFT 密度互相关求平移 + 空间群算符 + 格点平移”逐链
对齐到 `built_reference.pdb` 所在的位置，`meta.json.map_frame` 里记了
`shift_fractional`、`symmetry_triplet` 和 `chain_map`（逐链的算符与格点平移）。
"""


def main():
    for name, cases, writer, readme in (
            ('refmac_set', REFMAC, refmac_case, REFMAC_README),
            ('buccaneer_set', BUCCANEER, buccaneer_case, BUCCANEER_README)):
        root = OUT / name
        if root.exists():
            shutil.rmtree(root)
        (root / 'cases').mkdir(parents=True)
        for pdb in cases:
            writer(INV[pdb], root / 'cases' / pdb)
            print('packed', name, pdb, flush=True)
        (root / 'README.md').write_text(readme)
        manifest = []
        for p in sorted(root.rglob('*')):
            if p.is_file():
                manifest.append(dict(path=str(p.relative_to(root)), bytes=p.stat().st_size,
                                     sha256=sha256(p)))
        (root / 'manifest.json').write_text(json.dumps(
            dict(bundle=name, cases=cases, files=manifest), indent=1, ensure_ascii=False))
        print(name, len(manifest), 'files,',
              round(sum(f['bytes'] for f in manifest) / 1e6, 1), 'MB', flush=True)


if __name__ == '__main__':
    main()
