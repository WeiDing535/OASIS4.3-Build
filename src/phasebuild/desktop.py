"""Native desktop window for OASIS4.3.

The window is a thin front end over :mod:`phasebuild.jobs` and
:mod:`phasebuild.oasis_workflow`: it collects the same configuration the command
line accepts, starts any number of jobs (each in its own worker thread) and renders the
task list, the pipeline log and the refinement table while they run.

Tkinter ships with the bundled Python, so the desktop window needs no extra
installation.  ``Controller`` is toolkit-free and can be driven from tests.
"""
import argparse
import json
import os
import shutil
import signal
import subprocess
import sys
import time
from pathlib import Path

from . import charts
from .jobs import Jobs

STATUS_NAMES={'queued':'等待运行','running':'正在运行','stopping':'正在停止','cancelled':'已停止',
              'failed':'失败','success':'已完成','imported':'已导入（仅查看）','interrupted':'界面曾中断'}
IMPORTED_NOTE='导入的既有结果：只读目录，不会重新运行；可查看曲线/表格、导出结果、保存模型。'
METHODS=[('pplus_bundled','P+ 自带引擎（推荐）'),('pplus_trust','Phaser SAD + 外部 P+'),
         ('pplus_internal','内部 P+（实验性，不用于正式解析）'),('likelihood','连续 SAD'),
         ('pplus','旧版直接法 P+')]
METHOD_NAMES=dict(METHODS)
# Shown under the method radios; the internal chain is the one that must not be
# mistaken for a validated route (see docs/NO_CCP4_REPLACEMENT.md).
METHOD_HINTS={
 'pplus_bundled':'标定配方：P+ 定相 + Parrot 密度修饰 + Buccaneer/REFMAC，十案例 R-free 落在对照 ±5 点内。',
 'pplus_trust':'外部对照：Phaser 做 SAD，P+ 只接密度修饰之后的相位。',
 'pplus_internal':'实验性：SAD 定相已到 Phaser 水平，但密度修饰比 Parrot 低约 0.13（图相关 0.56 vs 0.69），'
                  '且没有自动追踪；结果只用于研究对比，不要用来定结构。',
 'likelihood':'连续 SAD：不经过 P+，直接用似然法定相，作为对照使用。',
 'pplus':'旧版直接法 P+：保留用于复现历史结果。'}
PROTOCOLS=[('','不使用参数方案'),('pplus_w8_weight','P+ · Wang8＋精修权重搜索（标定配方）'),
           ('pplus_w8','P+ · Wang8'),('pplus_w14','P+ · Wang14'),('stable_w14','稳定化 P+ · Wang14')]
RECIPES=[('pipeline','自带 Buccaneer/REFMAC'),('ccp4','外部 CCP4 流程'),('legacy','旧版逐轮流程')]
ENGINES=[('parrot','CCP4 Parrot（推荐）'),('resolve','RESOLVE 密度修饰')]
PRIORS=[('dm_phifom','密度修饰相位 + 模型振幅'),('original_sad','原始 SAD 相位')]

# key, title, kind, hint, options, only-methods
FIELD_SPECS=[
 ('sad','SAD 异常衍射 MTZ','file','需要 Friedel 配对振幅及各自误差；不要求已有相位。',None,None),
 ('heavy','重原子亚结构 PDB / mmCIF','file','对应本数据的 SHELXD 子结构；自带 P+ 可判手。',None,None),
 ('sequence','蛋白序列 FASTA / PIR','file','标准氨基酸；多条记录分别作为不同链。',None,None),
 ('native','原始观测 MTZ（可选）','file','留空则使用 SAD 文件并保留或新建自由集。',None,None),
 ('elements','重原子元素（必填）','text','只写元素名，例如 Se 或 Se,S；f″ 自动计算。',None,None),
 ('wavelength','实验 SAD 波长（Å，必填）','float','实际采集波长；程序不估算。',None,None),
 ('observation_prior','Phaser SAD 纯观测似然 MTZ','file','P+ 优化必填；禁止下游地图或模型相位。',None,('pplus_trust',)),
 ('composition_sequence','显式元素组成序列（可选）','file','SeMet 请写为 (MSE)；留空使用建模序列。',None,('pplus_trust','pplus_internal','pplus_bundled')),
 ('output','新的输出目录','output','目录必须尚不存在；所有结果保存在其中。',None,None),
 ('copies','整个序列文件在 ASU 中的份数','int','正整数 1–100，或 auto 按 FASTA 头 / Matthews 自动判定（默认）；不是单独一条链的拷贝数。',None,None),
 ('solvent_fraction','溶剂比例（留空自动估算）','optional','0–1 之间；留空按序列与晶胞估算。',None,None),
 ('phase_cycles','定相迭代周期','int','1–100。',None,None),
 ('dm_cycles','密度修饰周期','int','1–100；标定配方为 40。',None,None),
 ('dm_mask_radius','掩膜过滤半径（Å）','float','Parrot 溶剂掩膜的过滤半径；默认 3.0。',None,None),
 ('feedback_cycles','trace 反馈轮数','int','每轮用部分模型 Fcalc 反馈给密度修饰；0 表示不反馈，标定配方为 2–3。',None,None),
 ('rounds','Buccaneer / REFMAC 轮数','int','1–20；默认 5 轮。',None,None),
 ('cycles','每轮 REFMAC 精修周期','int','1–200；默认每轮 10 周期。',None,None),
 ('trust_wang_radius','单流程 Wang 掩膜半径（Å）','float','仅单流程使用；多流程固定 14 / 8 Å。',None,('pplus_trust','pplus_bundled')),
 ('trust_kl','P+ 优化：每反射 KL 上限','float','默认 0.0005；必须大于 0 且不超过 0.05。',None,('pplus_trust','pplus_internal','pplus_bundled')),
 ('max_per_target','每个目标的关系上限','int','默认 48；0 表示不截断。',None,None),
 ('seed','随机种子','int','0–2147483647；用于可重复的随机步骤。',None,None),
 ('labels','SAD 列映射（留空自动识别）','text','平均值：f=FP,sigf=SIGFP,dano=DANO,sigdano=SIGDANO。',None,None),
 ('free_value','原自由集数值','int','默认 0；须与原始 FreeR_flag 定义一致。',None,None),
 ('reference','密度修饰参考 MTZ','file','随运行资源提供的通用无关参考。',None,None),
 ('runtime','CCP4 位置（可选）','directory',
  '留空用包内运行时。CCP4 装在别处时填它的根目录（含 bin/），或直接粘 ccp4.setup / bashrc 的路径，点“重新检测”会解析成根目录；定相与后续模块都改用这个安装。',None,None),
 ('feedback_endpoints','反馈终点组合（留空=单终点）','optional','逗号分隔的反馈轮数，例如 2,3：每个值各跑一遍同一协议，再按工作集 R 选择终点（不读 R-free）。',None,None),
 ('dm_solvents','反馈终点溶剂组合（可选）','optional','与上一项配对，逗号分隔，auto 表示按序列估算，例如 auto,0.5；留空则沿用溶剂比例。',None,None),
]

# What an experiment needs is on the form; the calibrated constants, the run
# resources and the rare overrides wait behind “显示高级参数”.
ESSENTIAL_GROUPS=[('输入文件',('sad','heavy','sequence','elements','wavelength','copies','output')),
                  ('P+ 定相',('phase_cycles','trust_kl','max_per_target')),
                  ('Parrot 密度修饰',('solvent_fraction','dm_cycles','dm_mask_radius','feedback_cycles')),
                  ('Buccaneer / REFMAC',('rounds','cycles'))]
# Fields that live in the 「运行模块」 box rather than under a parameter group.
MODULE_FIELDS=('runtime',)
ADVANCED_GROUPS=[('输入细节',('native','observation_prior','composition_sequence','labels','free_value')),
                 ('P+ 掩膜与采样',('trust_wang_radius','seed')),
                 ('运行资源与终点组合',('reference','feedback_endpoints','dm_solvents'))]


class Controller:
    """Toolkit-free glue between the form, the job runner and the workflow."""

    def __init__(self, project, workflow=None):
        from . import oasis_workflow
        self.project=Path(project).resolve()
        self.workflow=workflow or oasis_workflow
        self.jobs=Jobs(self.project,workflow=self.workflow)
        self._metrics={}

    def defaults(self):
        return self.workflow.defaults(self.project)

    def plan(self, config):
        return self.workflow.make_plan(config,self.project)

    def describe(self, config):
        return self.workflow.describe_plan(config)

    def fallback_notice(self, config):
        """Red-line text for a run that fell back from Parrot, else ''."""
        from .sad_pipeline import fallback_notice
        return fallback_notice(config or {})

    def start(self, config, name=''):
        return self.jobs.start(config,name)

    def state(self, identifier):
        return self.jobs.get(identifier)

    def results(self, identifier):
        return self.workflow.collect_results(self.state(identifier)['output'])

    def stop(self, identifier):
        self.jobs.cancel(identifier)

    def history(self):
        """Task list for the window: job records plus their final R factors."""
        rows=self.jobs.history()
        for row in rows:
            row.update(self._metric(row))
        return rows

    def _metric(self,row):
        """Read the final statistics of one job, cached per (dir,status,update)."""
        key=(row['output'],row['status'],row.get('updated'))
        if key in self._metrics:
            return self._metrics[key]
        value={}
        if row['status'] in ('success','imported','cancelled','failed'):
            try:
                result=self.workflow.collect_results(row['output'])
                primary=result.get('primary_metric') or {}
                value=dict(method_name=result.get('method_name',''),
                           r_work=primary.get('r_work'),r_free=primary.get('r_free'),
                           stage=primary.get('stage',''),
                           endpoint=(result.get('endpoint_root') or '').split('/')[-1])
            except Exception:
                value={}
        self._metrics[key]=value
        return value

    def running(self):
        return self.jobs.live()

    def import_results(self, path):
        return self.jobs.import_results(path)

    def delete(self, identifier):
        """Drop one history record; its output directory stays on disk."""
        job=self.jobs.delete(identifier)
        self._metrics={key:value for key,value in self._metrics.items() if key[0]!=job['output']}
        return job

    def shutdown(self):
        """Stop every running job and make sure no engine child survives exit."""
        with self.jobs.lock:
            live=self.jobs.live()
            processes=dict(self.jobs.processes)
            workers=dict(self.jobs.workers)
        for identifier in live:
            self.jobs.cancel(identifier)
        for worker in workers.values():
            worker.join(timeout=2.5)
        with self.jobs.lock:
            processes={**processes,**self.jobs.processes}
        for process in processes.values():
            if process and process.poll() is None:
                try:
                    os.killpg(process.pid,9)
                except (ProcessLookupError,PermissionError):
                    pass
        for worker in workers.values():
            worker.join(timeout=5)


def _stamp():
    return time.strftime('%Y%m%d_%H%M%S')


class App:
    """Tk window; created only when a display is available."""

    def __init__(self, controller, config=None, autorun=False, exit_on_finish=False):
        import tkinter as tk
        from tkinter import ttk
        self.tk,self.ttk=tk,ttk
        self.controller=controller
        self.exit_on_finish=exit_on_finish
        self.defaults=controller.defaults()
        self.fields={}
        self.job=None
        self.result=None
        self.error_text=''
        self.finished=set()
        self.log_text=''
        self.history_signature=None
        self.root=tk.Tk()
        self.follow=tk.BooleanVar(master=self.root,value=True)
        self.root.title('OASIS4.3 · 单波长异常散射定相')
        self.root.geometry('1200x820')
        self.root.minsize(980,680)
        self.root.protocol('WM_DELETE_WINDOW',self.on_close)
        self.root.bind('<Control-n>',lambda event:self.on_new())
        self.root.bind('<Command-n>',lambda event:self.on_new())
        self._styles()
        self._layout()
        self._fill(self.defaults)
        if config:
            self._fill(config,strict=False)
        self._method_changed()
        self.root.after(300,self.poll)
        if autorun:
            self.root.after(600,self.on_start)

    # ---------------------------------------------------------------- layout
    def _styles(self):
        style=self.ttk.Style(self.root)
        if 'aqua' in style.theme_names():
            style.theme_use('aqua')
        style.configure('Hint.TLabel',foreground='#6b6b6b')
        style.configure('Status.TLabel',foreground='#1f5f9e')
        style.configure('Bad.TLabel',foreground='#a3242b')

    def _layout(self):
        tk,ttk=self.tk,self.ttk
        bar=ttk.Frame(self.root,padding=(10,8));bar.pack(side='top',fill='x')
        ttk.Label(bar,text='OASIS4.3',font=('Helvetica',16,'bold')).pack(side='left')
        ttk.Label(bar,text='单波长异常散射 → 密度修饰（可选）→ Buccaneer / REFMAC（可选）',
                  style='Hint.TLabel').pack(side='left',padx=(10,0),pady=(4,0))
        self.status=ttk.Label(bar,text='就绪',style='Status.TLabel')
        self.status.pack(side='right')

        self.warning=ttk.Label(self.root,text='',style='Bad.TLabel',wraplength=1120,justify='left')
        self.panes=ttk.Panedwindow(self.root,orient='horizontal')
        panes=self.panes
        panes.pack(fill='both',expand=True,padx=10,pady=(0,6))

        left=ttk.Frame(panes)
        canvas=tk.Canvas(left,highlightthickness=0,width=520)
        scroll=ttk.Scrollbar(left,orient='vertical',command=canvas.yview)
        canvas.configure(yscrollcommand=scroll.set)
        scroll.pack(side='right',fill='y');canvas.pack(side='left',fill='both',expand=True)
        form=ttk.Frame(canvas,padding=(0,0,8,8))
        window=canvas.create_window((0,0),window=form,anchor='nw')
        form.bind('<Configure>',lambda e:canvas.configure(scrollregion=canvas.bbox('all')))
        canvas.bind('<Configure>',lambda e:canvas.itemconfigure(window,width=e.width))
        canvas.bind_all('<MouseWheel>',lambda e:canvas.yview_scroll(int(-e.delta/40),'units'))
        self.form=form
        panes.add(left,weight=3)

        right=ttk.Frame(panes)
        self.tabs=ttk.Notebook(right)
        self.tabs.pack(fill='both',expand=True)
        self._history_tab();self._log_tab();self._results_tab();self._structure_tab();self._plan_tab()
        panes.add(right,weight=4)

        foot=ttk.Frame(self.root,padding=(10,0,10,8));foot.pack(side='bottom',fill='x')
        self.progress=ttk.Progressbar(foot,mode='determinate',maximum=100)
        self.progress.pack(side='left',fill='x',expand=True,padx=(0,10))
        self.new_button=ttk.Button(foot,text='新建任务',command=self.on_new)
        self.new_button.pack(side='left')
        self.start_button=ttk.Button(foot,text='开始运行',command=self.on_start)
        self.start_button.pack(side='left',padx=(6,0))
        self.stop_button=ttk.Button(foot,text='停止',command=self.on_stop,state='disabled')
        self.stop_button.pack(side='left',padx=6)
        ttk.Button(foot,text='打开结果目录',command=self.on_open).pack(side='left')
        self.export_button=ttk.Button(foot,text='导出结果…',command=self.on_export)
        self.export_button.pack(side='left',padx=6)
        ttk.Button(foot,text='导入已有结果…',command=self.on_import).pack(side='left')

        self._build_form()

    def _history_tab(self):
        tk,ttk=self.tk,self.ttk
        frame=ttk.Frame(self.tabs)
        tools=ttk.Frame(frame);tools.pack(fill='x',pady=(4,2))
        ttk.Label(tools,text='选中一行即在右侧各页查看；双击打开结果目录。',
                  style='Hint.TLabel').pack(side='left')
        ttk.Button(tools,text='刷新',command=lambda:self._refresh_history(force=True)).pack(side='right')
        self.delete_button=ttk.Button(tools,text='删除记录',command=self.on_delete_history)
        self.delete_button.pack(side='right',padx=6)
        self.relaunch_button=ttk.Button(tools,text='沿用参数',command=self.on_relaunch)
        self.relaunch_button.pack(side='right')
        columns=('started','name','status','method','rwork','output')
        self.history=ttk.Treeview(frame,columns=columns,show='headings',height=10)
        for key,title,width in [('started','开始时间',118),('name','任务名称',150),('status','状态',96),
                                ('method','方法',140),('rwork','最终 R-work',96),('output','输出目录',250)]:
            self.history.heading(key,text=title)
            self.history.column(key,width=width,anchor='center' if key in ('started','status','rwork') else 'w')
        y=ttk.Scrollbar(frame,orient='vertical',command=self.history.yview)
        x=ttk.Scrollbar(frame,orient='horizontal',command=self.history.xview)
        self.history.configure(yscrollcommand=y.set,xscrollcommand=x.set)
        y.pack(side='right',fill='y');x.pack(side='bottom',fill='x')
        self.history.pack(fill='both',expand=True)
        self.history.bind('<<TreeviewSelect>>',self.on_select_job)
        self.history.bind('<Double-1>',lambda event:self.on_open())
        self.tabs.add(frame,text='任务记录')

    def _log_tab(self):
        tk,ttk=self.tk,self.ttk
        frame=ttk.Frame(self.tabs)
        tools=ttk.Frame(frame);tools.pack(fill='x',pady=(4,2))
        ttk.Checkbutton(tools,text='自动滚动',variable=self.follow).pack(side='left')
        ttk.Button(tools,text='清空显示',command=lambda:(self.log.delete('1.0','end'))).pack(side='right')
        self.log=tk.Text(frame,wrap='none',height=20,font=('Menlo',11),background='#101418',
                         foreground='#d8e0e8',insertbackground='#d8e0e8')
        y=ttk.Scrollbar(frame,orient='vertical',command=self.log.yview)
        x=ttk.Scrollbar(frame,orient='horizontal',command=self.log.xview)
        self.log.configure(yscrollcommand=y.set,xscrollcommand=x.set)
        y.pack(side='right',fill='y');x.pack(side='bottom',fill='x')
        self.log.pack(fill='both',expand=True)
        self.tabs.add(frame,text='运行日志')

    def _results_tab(self):
        ttk=self.ttk
        frame=ttk.Frame(self.tabs)
        tools=ttk.Frame(frame);tools.pack(fill='x',pady=(4,2))
        ttk.Label(tools,text='R-work / R-free 为百分数；最后一行是固定后精修。',
                  style='Hint.TLabel').pack(side='left')
        self.save_model_button=ttk.Button(tools,text='保存最终模型…',command=self.on_save_model,
                                          state='disabled')
        self.save_model_button.pack(side='right')
        ttk.Button(tools,text='导出结果…',command=self.on_export).pack(side='right',padx=6)
        self.summary=self.ttk.Label(frame,text='尚无结果',justify='left',anchor='w',padding=(4,6))
        self.summary.pack(fill='x')
        self.chart=MetricChart(self,frame)
        self.chart.frame.pack(fill='x',pady=(2,6))
        columns=('stage','engine','rwork','rfree','cycles')
        self.tree=ttk.Treeview(frame,columns=columns,show='headings',height=10)
        for key,title,width in [('stage','阶段',140),('engine','引擎',220),('rwork','R-work',90),
                                ('rfree','R-free',90),('cycles','精修周期',80)]:
            self.tree.heading(key,text=title)
            self.tree.column(key,width=width,anchor='center' if key!='engine' else 'w')
        ty=ttk.Scrollbar(frame,orient='vertical',command=self.tree.yview)
        tx=ttk.Scrollbar(frame,orient='horizontal',command=self.tree.xview)
        self.tree.configure(yscrollcommand=ty.set,xscrollcommand=tx.set)
        ty.pack(side='right',fill='y');tx.pack(side='bottom',fill='x')
        self.tree.pack(fill='both',expand=True)
        self.tabs.add(frame,text='精修结果')

    def _structure_tab(self):
        frame=self.ttk.Frame(self.tabs)
        self.structure=StructureView(self,frame)
        self.tabs.add(frame,text='结构 3D')

    def _plan_tab(self):
        tk,ttk=self.tk,self.ttk
        frame=ttk.Frame(self.tabs)
        self.plan_text=tk.Text(frame,wrap='word',height=18,font=('Menlo',11))
        y=ttk.Scrollbar(frame,orient='vertical',command=self.plan_text.yview)
        self.plan_text.configure(yscrollcommand=y.set)
        y.pack(side='right',fill='y');self.plan_text.pack(fill='both',expand=True)
        self.plan_text.insert('end','填写输入后点击“检查配置”，这里显示将要执行的阶段与生效参数。')
        self.plan_text.configure(state='disabled')
        self.tabs.add(frame,text='运行计划')

    # ------------------------------------------------------------------ form
    def _section(self,title,parent=None):
        box=self.ttk.LabelFrame(parent or self.form,text=title,padding=(10,6))
        box.pack(fill='x',pady=(0,8))
        return box

    def _build_form(self):
        ttk=self.ttk
        specs={spec[0]:spec for spec in FIELD_SPECS}
        task=self._section('任务名称')
        row=ttk.Frame(task);row.pack(fill='x',pady=2)
        ttk.Label(row,text='任务名称',width=26,anchor='w').pack(side='left')
        self.task_name=self.tk.StringVar(value='')
        ttk.Entry(row,textvariable=self.task_name).pack(side='left',fill='x',expand=True)
        ttk.Label(task,text='可选：给这次任务起个名字（如“7E1D 重跑”），任务记录里显示它，'
                            '相同参数、不同数据也能一眼分清。',
                  style='Hint.TLabel',wraplength=340,justify='left').pack(fill='x',padx=(2,0))
        modules=self._section('运行模块（按本机 CCP4 检测结果默认勾选）')
        row=ttk.Frame(modules);row.pack(fill='x',pady=2)
        self.module_phase=self.tk.BooleanVar(value=True)
        ttk.Checkbutton(row,text='模块 1 · P+ 定相（必做）',variable=self.module_phase,
                        state='disabled').pack(side='left')
        self.module_parrot=self.tk.BooleanVar(value=True)
        self.parrot_check=ttk.Checkbutton(row,text='模块 2 · Parrot 密度修饰',
                                          variable=self.module_parrot,command=self._modules_changed)
        self.parrot_check.pack(side='left',padx=(18,0))
        row=ttk.Frame(modules);row.pack(fill='x',pady=(2,0))
        self.module_build=self.tk.BooleanVar(value=True)
        self.build_check=ttk.Checkbutton(row,text='模块 3 · Buccaneer 建模 + REFMAC 精修',
                                         variable=self.module_build,command=self._modules_changed)
        self.build_check.pack(side='left')
        self._field(modules,specs['runtime'])
        # The button acts on the location above it, so it sits with that field
        # rather than trailing the module checkboxes.
        row=ttk.Frame(modules);row.pack(fill='x',pady=(4,0))
        ttk.Button(row,text='重新检测 CCP4',command=self._detect_modules).pack(side='left')
        self.module_hint=self.tk.StringVar(value='正在检测本机 CCP4…')
        ttk.Label(modules,textvariable=self.module_hint,style='Hint.TLabel',wraplength=340,
                  justify='left').pack(fill='x',padx=(2,0),pady=(4,0))
        self._detect_modules(initial=True)

        inputs=self._section('输入文件')
        for key in ESSENTIAL_GROUPS[0][1]:
            self._field(inputs,specs[key])

        phase=self._section('① P+ 定相参数')
        self.method=self.tk.StringVar(value=self.defaults.get('method','pplus_bundled'))
        for key,label in METHODS:
            ttk.Radiobutton(phase,text=label,value=key,variable=self.method,
                            command=self._method_changed).pack(anchor='w')
        self.method_hint=self.tk.StringVar(value=METHOD_HINTS.get(self.method.get(),''))
        ttk.Label(phase,textvariable=self.method_hint,style='Hint.TLabel',wraplength=340,
                  justify='left').pack(fill='x',padx=(2,0),pady=(4,0))
        row=ttk.Frame(phase);row.pack(fill='x',pady=(4,0))
        ttk.Label(row,text='P+ 参数方案').pack(side='left')
        self.protocol=self.tk.StringVar(value='')
        ttk.Combobox(row,textvariable=self.protocol,values=[label for _,label in PROTOCOLS],
                     state='readonly',width=34).pack(side='left',padx=6)
        for key in ESSENTIAL_GROUPS[1][1]:
            self._field(phase,specs[key])

        dm=self._section('② Parrot 密度修饰参数')
        row=ttk.Frame(dm);row.pack(fill='x',pady=2)
        ttk.Label(row,text='密度修饰引擎',width=18,anchor='w').pack(side='left')
        self.dm_engine=self.tk.StringVar(value='parrot')
        self.dm_engine_box=ttk.Combobox(row,textvariable=self.dm_engine,values=[label for _,label in ENGINES],
                                        state='readonly',width=22)
        self.dm_engine_box.pack(side='left',fill='x',expand=True,padx=(0,2))
        ttk.Label(dm,text='Parrot 是十案例标定的配方。内部密度修饰只能由命令行/环境变量开启，'
                          '属实验性（图相关 0.56 vs Parrot 0.69），窗口里不提供。',
                  style='Hint.TLabel',wraplength=340,justify='left').pack(fill='x',padx=(2,0))
        for key in ESSENTIAL_GROUPS[2][1]:
            self._field(dm,specs[key])

        build=self._section('③ Buccaneer 建模 + REFMAC 精修参数')
        row=ttk.Frame(build);row.pack(fill='x',pady=2)
        self.auto_hand=self.tk.BooleanVar(value=True)
        self.auto_hand_box=ttk.Checkbutton(row,text='自动判手（双手试建，不用 R-free）',variable=self.auto_hand)
        self.auto_hand_box.pack(side='left')
        row=ttk.Frame(build);row.pack(fill='x',pady=(0,2))
        self.unverified_hand=self.tk.BooleanVar(value=False)
        self.unverified_hand_box=ttk.Checkbutton(row,text='证据不足也按输入手继续（结果标注未验证）',variable=self.unverified_hand)
        self.unverified_hand_box.pack(side='left')
        ttk.Label(build,text='默认关闭：两手都算完但测序证据不足时，程序按失败处理并保留两手输出。'
                             '勾选后优先沿用输入手继续建模，并在 hand_selection.json 与流水线报告里标注“手性未验证”。',
                  style='Hint.TLabel',wraplength=340,justify='left').pack(fill='x',padx=(2,0))
        for key in ESSENTIAL_GROUPS[3][1]:
            self._field(build,specs[key])
        for key,title,values,var in [('building_recipe','构建配方',RECIPES,'recipe'),
                                     ('refinement_prior','精修相位来源',PRIORS,'prior')]:
            row=ttk.Frame(build);row.pack(fill='x',pady=2)
            ttk.Label(row,text=title,width=18,anchor='w').pack(side='left')
            variable=self.tk.StringVar(value=values[0][1])
            box=ttk.Combobox(row,textvariable=variable,values=[label for _,label in values],
                             state='readonly',width=22)
            box.pack(side='left',fill='x',expand=True,padx=(0,2))
            setattr(self,var,variable);setattr(self,var+'_values',values);setattr(self,var+'_box',box)
        # One name for the module-3 switch: saved jobs and tests use build_model.
        self.build_model=self.module_build

        # The advanced parameters stay out of the way until they are asked for:
        # they are calibrated constants or run resources, not experiment inputs.
        self.advanced_box=ttk.Frame(self.form)
        for title,keys in ADVANCED_GROUPS:
            group=self._section(title,parent=self.advanced_box)
            for key in keys:
                self._field(group,specs[key])

        actions=ttk.Frame(self.form);actions.pack(fill='x',pady=(0,6))
        self.advanced=self.tk.BooleanVar(value=False)
        ttk.Checkbutton(actions,text='显示高级参数（一般保持默认）',variable=self.advanced,
                        command=self._toggle_advanced).pack(side='left')
        ttk.Button(actions,text='检查配置',command=self.on_check).pack(side='right')
        self.actions=actions

    def _toggle_advanced(self):
        if self.advanced.get():
            self.advanced_box.pack(fill='x',pady=(0,8),before=self.actions)
        else:
            self.advanced_box.pack_forget()

    def _field(self,parent,spec):
        key,title,kind,hint,_,only=spec
        ttk=self.ttk
        var=self.tk.StringVar(value='')
        if kind in ('file','directory','output'):
            # Paths get their own full-width line: the label on top, then the
            # entry and its 浏览… button.  Sharing one line with the label left
            # the entry a few characters wide and squeezed the button out.
            caption=ttk.Frame(parent);caption.pack(fill='x',pady=(4,0))
            ttk.Label(caption,text=title,anchor='w').pack(side='left')
            row=ttk.Frame(parent);row.pack(fill='x',pady=(1,0))
            ttk.Button(row,text='浏览…',
                       command=lambda:self._browse(key,var,kind)).pack(side='right')
            ttk.Entry(row,textvariable=var).pack(side='left',fill='x',expand=True,padx=(0,4))
        else:
            row=ttk.Frame(parent);row.pack(fill='x',pady=2)
            ttk.Label(row,text=title,width=26,anchor='w').pack(side='left')
            if kind=='text':
                ttk.Entry(row,textvariable=var).pack(side='left',fill='x',expand=True)
            else:
                ttk.Entry(row,textvariable=var,width=16).pack(side='left')
        ttk.Label(parent,text=hint,style='Hint.TLabel',wraplength=340,
                  justify='left').pack(fill='x',padx=(2,0))
        self.fields[key]=dict(var=var,kind=kind,title=title,only=only,row=row)

    def _browse(self,key,var,kind):
        from tkinter import filedialog
        current=var.get().strip() or str(self.controller.project)
        if kind=='directory':
            path=filedialog.askdirectory(initialdir=current,title='选择目录')
        elif kind=='output':
            path=filedialog.askdirectory(initialdir=current,title='选择新输出目录的上级目录')
            if path:
                path=str(Path(path)/('oasis_'+_stamp()))
        else:
            path=filedialog.askopenfilename(initialdir=str(Path(current).parent),
                                            title='选择文件')
        if path:
            var.set(path)

    # ------------------------------------------------------------- modules
    def _detect_modules(self,initial=False):
        """Default the optional modules from what this machine actually has."""
        from . import module_detection
        field=getattr(self,'fields',{}).get('runtime')
        typed=field['var'].get().strip() if field is not None else ''
        fallback=self.defaults.get('runtime') or None
        report=module_detection.detect(typed or fallback,location=typed or None)
        resolved=(report.get('location') or {}).get('resolved')
        self._location_warning=''
        if typed and resolved and field is not None:
            # A setup script (or a bin/ dir) is not a runtime: keep the root it
            # named so the run uses the installation the detection just proved.
            field['var'].set(resolved)
        elif typed and not resolved:
            # Nothing at the typed path.  Leaving it in the form would hand the
            # run a directory with no engines, so fall back and say so.
            if field is not None and fallback:
                field['var'].set(fallback)
            self._location_warning='“%s”里没找到 CCP4 程序，已改用 %s。'%(
                typed,fallback or '默认运行时')
            report=module_detection.detect(fallback)
        self._module_report=report
        modules=report['modules']
        self.module_parrot.set(modules['parrot']['available'])
        self.module_build.set(modules['build']['available'])
        self._modules_changed()

    def _modules_changed(self,*_):
        """Building needs a density-modified map, so module 3 follows module 2."""
        modules=(getattr(self,'_module_report',None) or {}).get('modules') or {}
        parrot_ok=bool((modules.get('parrot') or {}).get('available'))
        build_ok=bool((modules.get('build') or {}).get('available'))
        if not parrot_ok:
            self.module_parrot.set(False)      # a reused record cannot ask for what is missing
        if not self.module_parrot.get():
            self.module_build.set(False)
        self.parrot_check.configure(state='normal' if parrot_ok else 'disabled')
        self.build_check.configure(state='normal' if build_ok and self.module_parrot.get() else 'disabled')
        if not self.module_build.get():
            hand=getattr(self,'auto_hand',None)
            if hand is not None:
                hand.set(False)
        self._sync_module_state()
        self._refresh_module_hint()

    def _sync_module_state(self):
        on={'parrot':self.module_parrot.get(),'build':self.module_build.get()}
        for index,module in ((2,'parrot'),(3,'build')):
            state='normal' if on[module] else 'disabled'
            for key in ESSENTIAL_GROUPS[index][1]:
                field=getattr(self,'fields',{}).get(key)
                if not field:
                    continue
                for child in field['row'].winfo_children():
                    if child.winfo_class() not in ('TEntry','TButton'):
                        continue    # the caption label of a path row stays readable
                    try:
                        child.configure(state=state)
                    except Exception:
                        pass
        for name,module,enabled in [('dm_engine_box','parrot','readonly'),('auto_hand_box','build','normal'),('unverified_hand_box','build','normal'),
                                    ('recipe_box','build','readonly'),('prior_box','build','readonly')]:
            widget=getattr(self,name,None)
            if widget is None:
                continue
            try:
                widget.configure(state=enabled if on[module] else 'disabled')
            except Exception:
                pass
        self._method_changed()

    def _refresh_module_hint(self):
        report=getattr(self,'_module_report',None) or {}
        modules=report.get('modules') or {}
        lines=[]
        warning=getattr(self,'_location_warning','')
        if warning:lines.append(warning)
        located=(report.get('location') or {}).get('resolved') or (report.get('location') or {}).get('found_at') or ''
        if report.get('ccp4_ready'):
            lines.append('已检测到 CCP4：模块 2、3 默认勾选，可随时取消。')
            if located:lines.append('CCP4：'+located)
        elif report.get('available'):
            missing=[item['title'] for item in modules.values() if not item['available']]
            lines.append('只检测到部分 CCP4：默认勾选可用模块；未找到 '+'、'.join(missing)+'。')
            if located:lines.append('CCP4：'+located)
        else:
            lines.append('未检测到 CCP4（cparrot、cbuccaneer、refmac5），默认只做模块 1。')
            lines.append('模块 1 输出相位、FOM 与 MTZ；已装 CCP4 但不在默认位置时，在上方“CCP4 位置”里填根目录或 ccp4.setup 路径，再点“重新检测”。')
        if located and not (Path(located)/'density/bin/oasis-density-python').is_file():
            lines.append('注意：%s 是外部 CCP4 安装，没有包内密度修饰引擎（density/）；P+ 定相与双手判别需要便携包的 runtime，请把“CCP4 位置”指向 …/OASIS4.3-…-portable local/runtime。'%located)
        if not self.module_parrot.get():
            lines.append('模块 2 已关闭：跳过密度修饰，直接输出 P+ 相位；模块 3 需要修饰后的图，已一并关闭。')
        elif not self.module_build.get():
            lines.append('模块 3 已关闭：停在密度修饰后的相位与地图。')
        self.module_hint.set('\n'.join(lines))

    def _method_changed(self,*_):
        variable=getattr(self,'method',None)
        if variable is None:
            return
        method=variable.get()
        hint=getattr(self,'method_hint',None)
        if hint is not None:
            hint.set(METHOD_HINTS.get(method,''))
        for key,field in self.fields.items():
            only=field['only']
            if not only:
                continue
            widget=field['row'].winfo_children()[1]
            state='normal' if method in only else 'disabled'
            try:
                widget.configure(state=state)
            except Exception:
                pass

    # --------------------------------------------------------------- values
    def _fill(self,values,strict=True):
        method=values.get('method')
        if method in dict(METHODS):
            self.method.set(method)
        recipe=dict(self.recipe_values).get(values.get('building_recipe'),None)
        if recipe:
            self.recipe.set(recipe)
        prior=dict(self.prior_values).get(values.get('refinement_prior'),None)
        if prior:
            self.prior.set(prior)
        engine=dict(ENGINES).get(values.get('dm_engine'))
        if engine:
            self.dm_engine.set(engine)
        # Module 2/3 start from what this machine has; a reused record keeps its
        # own choice, so only a real config (strict=False) may override them.
        if not strict:
            if 'dm_enabled' in values:
                self.module_parrot.set(bool(values['dm_enabled']))
            if 'build_model' in values:
                self.build_model.set(bool(values['build_model']))
            if 'auto_hand' in values:
                self.auto_hand.set(bool(values['auto_hand']))
            if 'allow_unverified_hand' in values:
                self.unverified_hand.set(bool(values['allow_unverified_hand']))
        protocol=values.get('protocol','')
        label=dict(PROTOCOLS).get(protocol)
        if label:
            self.protocol.set(label)
        for key,field in self.fields.items():
            if key not in values:
                continue
            value=values[key]
            if isinstance(value,bool):
                value=''
            field['var'].set('' if value is None else str(value))
        if hasattr(self,'module_hint'):
            self._modules_changed()

    def collect(self):
        values=dict(self.defaults)
        for key,field in self.fields.items():
            text=field['var'].get().strip()
            kind=field['kind']
            if kind in ('int','float'):
                if not text:
                    raise ValueError('请填写：'+field['title'])
                if key=='copies' and text.lower()=='auto':
                    # ``auto`` is the documented default: FASTA header first,
                    # then the smallest Matthews-plausible number of copies.
                    values[key]='auto'
                    continue
                try:
                    values[key]=int(text) if kind=='int' else float(text)
                except ValueError:
                    raise ValueError('“%s”需要数字，当前为：%s'%(field['title'],text))
            elif kind=='optional':
                values[key]=text
            else:
                values[key]=text
        if not str(values.get('solvent_fraction') or '').strip():
            # Blank means "estimate from sequence and cell": None, not ''.
            values['solvent_fraction']=None
        values['method']=self.method.get()
        values['dm_engine']=dict((label,key) for key,label in ENGINES).get(self.dm_engine.get(),'parrot')
        values['building_recipe']=dict((label,key) for key,label in RECIPES).get(self.recipe.get(),'pipeline')
        values['refinement_prior']=dict((label,key) for key,label in PRIORS).get(self.prior.get(),'dm_phifom')
        values['dm_enabled']=bool(self.module_parrot.get())
        values['build_model']=bool(self.build_model.get()) and values['dm_enabled']
        values['auto_hand']=bool(self.auto_hand.get()) and values['build_model']
        values['allow_unverified_hand']=bool(self.unverified_hand.get()) and values['auto_hand']
        protocol=dict((label,key) for key,label in PROTOCOLS).get(self.protocol.get(),'')
        if protocol:
            values['protocol']=protocol
        else:
            values.pop('protocol',None)
        return values

    # --------------------------------------------------------------- actions
    def on_check(self):
        try:
            config=self.collect()
            plan=self.controller.describe(config)
        except Exception as error:
            self._error(str(error));return
        lines=[(stage['name'],stage.get('status') or '待运行') for stage in plan.get('stages',[])]
        text='\n'.join('· %s  [%s]'%(name,status) for name,status in lines)
        effective=plan.get('effective_options') or {}
        body='\n'.join('%s：%s'%(key,value) for key,value in effective.items())
        self._set_plan(text+('\n\n生效参数\n'+body if body else ''))
        self._warn('\n'.join(plan.get('sample_warnings') or []))
        self._say('配置检查通过',ok=True)

    def _density_runtime_fix(self,config):
        """Point a bundled-route run at a runtime that carries density/, else ''.

        A record made before the density tree was resolved separately, or a
        runtime detected as a bare CCP4, names an installation with bin/ but no
        density/.  Starting it fails before Phasing, so switch the form to the
        package runtime and say so instead of letting the job die.
        """
        from . import oasis_workflow
        if config.get('method')!='pplus_bundled':return ''
        runtime=str(config.get('runtime') or '').strip()
        if not runtime:return ''
        # Look at this runtime's own tree, not the engine's wider search: the
        # form should show the installation that actually carries density/.
        try:
            if (Path(runtime).expanduser()/'density/bin/oasis-density-python').is_file():return ''
        except (OSError,ValueError):
            return ''
        suggested=oasis_workflow.bundled_runtime(self.controller.project)
        if suggested is None or str(suggested)==runtime:return ''
        return str(suggested)

    def _use_runtime(self,path,old=''):
        """Show ``path`` in the form and explain why it replaced ``old``."""
        field=self.fields.get('runtime')
        if field is not None:field['var'].set(path)
        self._detect_modules()
        self._location_warning=('“CCP4 位置”已从 %s 改为 %s：前者没有 density/，跑不了包内 P+ 流程。'
                                %(old or '（空）',path))
        self._refresh_module_hint()

    def _apply_density_fix(self,config):
        """Apply ``_density_runtime_fix`` to the form and return the new config."""
        fixed=self._density_runtime_fix(config)
        if not fixed:return config
        self._use_runtime(fixed,str(config.get('runtime') or ''))
        return self.collect()

    def on_start(self):
        """Start one more run; earlier runs keep going in parallel."""
        try:
            config=self._apply_density_fix(self.collect())
            name=self.task_name.get().strip()
            job=self.controller.start(config,name)
        except Exception as error:
            self._error(str(error));return
        self.task_name.set('')
        self.job=job
        self._warn(self.controller.fallback_notice(job.get('config') or {}))
        self.log_text=''
        self.log.delete('1.0','end')
        self._clear_results()
        self.summary.configure(text='运行中：'+str(job['output']))
        self.stop_button.configure(state='normal')
        self.progress.configure(value=0)
        self.tabs.select(1)
        self._refresh_history(force=True)
        self._select_history(job['id'])
        self._say('已启动任务 %s%s（可继续新建任务并行运行）'
                  %(job.get('name') or job['id'],'' if not job.get('name') else '（%s）'%job['id']))
        self._show_plan_for(job)
        self._advance_output(job)

    def _advance_output(self,job):
        """Offer a fresh output directory so the next run can start at once."""
        field=self.fields.get('output')
        if not field:return
        try: current=Path(field['var'].get().strip()).expanduser().resolve()
        except (OSError,ValueError): return
        consumed=Path(job['output']).resolve()
        if current!=consumed and not current.exists(): return
        candidate=current.parent/('oasis_'+_stamp())
        suffix=2
        while candidate.exists():
            candidate=current.parent/('oasis_%s_%d'%(_stamp(),suffix));suffix+=1
        field['var'].set(str(candidate))

    def on_stop(self):
        if not self.job:return
        if self.job.get('status') not in ('queued','running','stopping'):
            self._say('该任务已经结束');return
        self.controller.stop(self.job['id'])
        self._say('正在停止 '+self.job['id']+' …')

    def on_open(self):
        if not self.job:
            self._error('还没有结果目录。');return
        path=self.job['output']
        opener='open' if sys.platform=='darwin' else 'xdg-open'
        try:
            subprocess.Popen([opener,path])
        except OSError as error:
            self._error(str(error))

    def on_new(self):
        """Clear the form for another experiment without touching the history."""
        self.job=None
        self.finished=set()
        self.log_text=''
        self.log.delete('1.0','end')
        self._clear_results()
        self._fill(self.controller.defaults())
        self.task_name.set('')
        self.fields['output']['var'].set(str(self._fresh_output()))
        self._detect_modules()
        self._method_changed()
        self.progress.configure(value=0)
        self.stop_button.configure(state='disabled')
        for item in self.history.selection():
            self.history.selection_remove(item)
        self.summary.configure(text='新任务：填写输入后点击“开始运行”')
        self._say('已新建任务：输出目录 '+self.fields['output']['var'].get())

    def _fresh_output(self):
        """A new, unused run directory next to the workflow's default one."""
        base=Path(self.controller.defaults()['output']).parent
        candidate=base/('oasis_'+_stamp())
        suffix=2
        while candidate.exists():
            candidate=base/('oasis_%s_%d'%(_stamp(),suffix));suffix+=1
        return candidate

    def on_export(self):
        """Copy the selected run directory to a place the user chooses."""
        if not self.job:
            self._error('还没有可导出的结果。');return
        source=Path(self.job['output'])
        if not source.is_dir():
            self._error('结果目录不存在：'+str(source));return
        from tkinter import filedialog
        target=filedialog.askdirectory(initialdir=str(Path.home()),title='选择导出位置')
        if not target:return
        destination=Path(target)/('oasis_export_'+_stamp())
        try:
            shutil.copytree(source,destination)
        except OSError as error:
            self._error('导出失败：'+str(error));return
        files=[item for item in destination.rglob('*') if item.is_file()]
        size=sum(item.stat().st_size for item in files)
        self._say('已导出 %d 个文件（%.1f MB）到 %s'%(len(files),size/1e6,destination),ok=True)

    def _model_path(self,result=None):
        """Absolute path of a result's final model (results store it relative)."""
        result=self.result if result is None else result
        model=(result or {}).get('final_model')
        if not model:
            return None
        path=Path(model)
        if path.is_absolute():
            return path
        base=(result or {}).get('root') or (self.job or {}).get('output')
        return Path(base)/path if base else path

    def on_save_model(self):
        """Save the final model as a single file (mmCIF or PDB)."""
        model=self._model_path()
        if not model or not model.is_file():
            self._error('还没有最终模型。');return
        from tkinter import filedialog
        path=filedialog.asksaveasfilename(defaultextension='.cif',initialfile=Path(model).name,
                                          initialdir=str(Path(model).parent),title='保存最终模型')
        if not path:return
        try:
            shutil.copyfile(model,path)
        except OSError as error:
            self._error('保存失败：'+str(error));return
        self._say('已保存模型到 '+path,ok=True)

    def on_import(self):
        from tkinter import filedialog
        path=filedialog.askdirectory(initialdir=str(self.controller.project),title='选择已有结果目录')
        if not path:return
        try:
            job=self.controller.import_results(path)
        except Exception as error:
            self._error(str(error));return
        self.job=job
        self._refresh_history(force=True)
        self._select_history(job['id'])
        self._show_results()

    def on_close(self):
        live=self.controller.running()
        if live and not self.exit_on_finish:
            from tkinter import messagebox
            question='还有 %d 个任务在运行，退出将停止它们。确定退出？'%len(live)
            if not messagebox.askokcancel('退出 OASIS4.3',question):
                return
        self.controller.shutdown()
        self.root.destroy()

    # ---------------------------------------------------------------- update
    def poll(self):
        try:
            self._refresh_history()
            if self.job and self.job.get('id') not in self.finished:
                self._refresh()
        finally:
            self.root.after(900,self.poll)

    def _refresh(self):
        identifier=self.job['id']
        try:
            job=self.controller.state(identifier)
        except KeyError:
            return
        self.job=job
        status=job.get('status')
        progress=job.get('pipeline_progress') or {}
        total=progress.get('total') or 0
        done=progress.get('completed') or 0
        if total:
            self.progress.configure(value=100*min(1,done/total))
        label=progress.get('current') or job.get('current') or ''
        self._say_status('%s %s'%(STATUS_NAMES.get(status,status),label))
        log=job.get('log') or ('（导入的既有结果没有运行日志）' if status=='imported' else '')
        if log!=self.log_text:
            self._show_log(log)
        if status in ('success','failed','cancelled','interrupted'):
            self.finished.add(identifier)
            self.stop_button.configure(state='disabled')
            self.progress.configure(value=100 if status=='success' else 0)
            if status=='success':
                self._show_results(select=False)
            elif status=='failed':
                self.error_text=job.get('error') or '任务失败，请查看日志'
                self.summary.configure(text='失败：'+self.error_text)
                self.tabs.select(1)
            if self.exit_on_finish and not self.controller.running():
                self.root.after(200,self.on_close)

    def _show_log(self,text):
        self.log_text=text
        self.log.delete('1.0','end');self.log.insert('end',text)
        if self.follow.get():
            self.log.see('end')

    def _say_status(self,label=''):
        live=self.controller.running()
        text=label or '就绪'
        if len(live)>1:
            text+='　·　共 %d 个任务运行中'%len(live)
        self._say(text)

    # ------------------------------------------------------------ task list
    def _refresh_history(self,force=False):
        try:
            rows=self.controller.history()
        except Exception:
            return
        signature=tuple((row['id'],row['status'],row.get('updated')) for row in rows)
        if not force and signature==self.history_signature:
            return
        self.history_signature=signature
        for item in self.history.get_children():
            self.history.delete(item)
        for row in rows:
            rwork=row.get('r_work')
            self.history.insert('','end',iid=row['id'],values=(
                time.strftime('%m-%d %H:%M',time.localtime(row['created'])),
                row.get('name') or row.get('label') or '',
                STATUS_NAMES.get(row['status'],row['status']),
                row.get('method_name') or METHOD_NAMES.get(row.get('method'),row.get('method') or ''),
                '' if rwork is None else '%.2f'%(100*rwork),
                row['output']))
        if self.job and self.history.exists(self.job['id']):
            self._select_history(self.job['id'])

    def on_relaunch(self):
        """Load a history record's parameters so a new run needs no retyping."""
        selection=self.history.selection()
        if not selection:
            self._error('请先在任务记录里选中一行，再点“沿用参数”。');return
        identifier=selection[0]
        try:
            job=self.controller.state(identifier)
        except KeyError:
            self._error('这条记录已经不存在。');return
        config=job.get('config')
        if not config:
            self._error('这条记录没有保存运行参数（导入的既有结果），无法沿用。');return
        self._fill(config,strict=False)
        self.task_name.set(job.get('name') or '')
        self.fields['output']['var'].set(str(self._fresh_output()))
        self._method_changed()
        # An older record may name a bare CCP4; show the runtime that can run it
        # rather than repeating a failure the user already saw.
        try:
            config=self._apply_density_fix(self.collect())
        except ValueError:
            pass
        if self._needs_advanced(config):
            self.advanced.set(True);self._toggle_advanced()
        self.job=None
        self._clear_results()
        self._show_log('')
        try:
            plan=self.controller.describe(dict(config,output=self.fields['output']['var'].get()))
            self._set_plan(json.dumps(plan,ensure_ascii=False,indent=2)[:20000])
        except Exception:
            pass
        self.tabs.select(4)
        self.summary.configure(text='已沿用任务“%s”的参数：改好输入数据后点“开始运行”，输出目录已换成新目录。'
                                    %(job.get('name') or job.get('label') or identifier))
        self._say('已沿用参数，可直接开始运行')

    def _needs_advanced(self,config):
        """True when the record used non-default advanced parameters."""
        advanced={key for _,keys in ADVANCED_GROUPS for key in keys}
        for key in advanced:
            if key not in config:
                continue
            old,new=config[key],self.defaults.get(key)
            if isinstance(old,float) and isinstance(new,(int,float)):
                if abs(old-new)>1e-9: return True
            elif str(old)!=str(new):
                return True
        return False

    def on_delete_history(self):
        """Forget selected records; the result directories stay untouched."""
        selection=list(self.history.selection())
        if not selection:
            self._error('请先在任务记录里选中要删除的行。');return
        busy={'等待运行','正在运行','正在停止'}
        targets=[i for i in selection if self.history.set(i,'status') not in busy]
        live=len(selection)-len(targets)
        if not targets:
            self._error('选中的任务还在运行，请先停止再删除记录。');return
        from tkinter import messagebox
        note='将从任务记录中删除 %d 条记录；结果目录与文件不会被删除。'%len(targets)
        if live: note+='\n正在运行的 %d 条记录会保留。'%live
        if not messagebox.askokcancel('删除任务记录',note):
            return
        removed=0
        for identifier in targets:
            try:
                self.controller.delete(identifier)
                removed+=1
            except Exception as error:
                self._error(str(error));return
        if self.job and self.job.get('id') in targets:
            self.finished.discard(self.job['id'])
            self.job=None
            self._clear_results()
            self._show_log('')
        self.history_signature=None
        self._refresh_history(force=True)
        self._say('已删除 %d 条任务记录（结果目录未改动）'%removed)

    def _select_history(self,identifier):
        if not self.history.exists(identifier):return
        self.history.selection_set(identifier)
        self.history.focus(identifier)
        self.history.see(identifier)

    def on_select_job(self,event=None):
        selection=self.history.selection()
        if not selection:return
        identifier=selection[0]
        if self.job and self.job.get('id')==identifier:return
        try:
            job=self.controller.state(identifier)
        except KeyError:
            return
        self.job=job
        status=job.get('status')
        self._show_log(job.get('log') or ('（导入的既有结果没有运行日志）' if status=='imported' else ''))
        self.stop_button.configure(state='normal' if status in ('queued','running','stopping') else 'disabled')
        self.progress.configure(value=0)
        self._show_plan_for(job)
        if status in ('success','imported'):
            if status=='success':self.finished.add(identifier)
            self._show_results(select=False)
            if status=='imported':
                self.summary.configure(text=IMPORTED_NOTE+'\n'+self.summary.cget('text'))
        elif status in ('failed','cancelled','interrupted'):
            self._clear_results()
            self.summary.configure(text=STATUS_NAMES.get(status,status)+'：'+str(job.get('error') or '查看日志'))
        else:
            self._clear_results()
            self.summary.configure(text='运行中：'+str(job.get('output')))
        self._say_status('%s %s'%(STATUS_NAMES.get(status,status),identifier))

    def _show_plan_for(self,job):
        plan=job.get('plan')
        if not plan and job.get('config'):
            try: plan=self.controller.describe(job['config'])
            except Exception: plan=None
        if plan:
            self._set_plan(json.dumps(plan,ensure_ascii=False,indent=2)[:20000])

    def _clear_results(self):
        self.result=None
        for row in self.tree.get_children():
            self.tree.delete(row)
        self.summary.configure(text='尚无结果')
        self.save_model_button.configure(state='disabled')
        self.chart.show([])
        self.structure.show(None)

    def _show_results(self,select=True):
        try:
            result=self.controller.results(self.job['id'])
        except Exception as error:
            self.summary.configure(text='读取结果失败：'+str(error));return
        for row in self.tree.get_children():
            self.tree.delete(row)
        for row in result.get('metrics') or []:
            rwork=row.get('r_work');rfree=row.get('r_free')
            self.tree.insert('','end',values=(row.get('stage',''),
                row.get('engine',''),'' if rwork is None else '%.2f'%(100*rwork),
                '' if rfree is None else '%.2f'%(100*rfree),row.get('cycles','')))
        self._warn(result.get('dm_engine_fallback') or '')
        primary=result.get('primary_metric') or {}
        lines=['方法：%s  密度修饰：%s'%(result.get('method_name',''),result.get('dm_engine',''))]
        if result.get('solvent_fraction') is not None:
            lines.append('溶剂比例：%.2f（%s）'%(result['solvent_fraction'],result.get('solvent_source') or '估算'))
        if primary:
            lines.append('最终 R-work %.2f / R-free %.2f（%s）'%(100*primary['r_work'],
                         100*primary['r_free'],primary.get('stage','')))
        endpoints=result.get('endpoints') or []
        if endpoints:
            parts=[]
            for item in endpoints:
                metric=item.get('primary_metric') or {}
                rwork=metric.get('r_work')
                solvent=item.get('solvent')
                parts.append('fb%s/溶剂%s → %s'%(item.get('feedback'),'auto' if solvent is None else '%.2f'%solvent,
                             '—' if rwork is None else 'R-work %.2f'%(100*rwork)))
            selected=(result.get('endpoint_selection') or {}).get('selected_feedback')
            lines.append('终点组合（按工作集 R 保留 fb%s）：%s'%(selected,'；'.join(parts)))
        if result.get('final_model'):
            lines.append('模型：'+result['final_model'])
        elif result.get('final_phases'):
            lines.append('相位：'+result['final_phases'])
        self.summary.configure(text='\n'.join(lines))
        self.result=result
        self.chart.show(charts.round_series(result.get('metrics')))
        model=self._model_path(result)
        self.save_model_button.configure(state='normal' if model and model.is_file() else 'disabled')
        self.structure.show(str(model) if model and model.is_file() else None,result.get('preview'))
        if select:
            self.tabs.select(2)

    def _set_plan(self,text):
        self.plan_text.configure(state='normal')
        self.plan_text.delete('1.0','end')
        self.plan_text.insert('end',text)
        self.plan_text.configure(state='disabled')

    def _say(self,text,ok=False):
        self.status.configure(text=text,style='Status.TLabel')

    def _warn(self,text=''):
        """Show one persistent red line under the toolbar; empty text hides it."""
        text=(text or '').strip()
        self.warning.configure(text=text)
        if text:
            if not self.warning.winfo_ismapped():
                self.warning.pack(side='top',fill='x',padx=14,pady=(0,4),before=self.panes)
        elif self.warning.winfo_ismapped():
            self.warning.pack_forget()

    def _error(self,text):
        self.error_text=text
        self.status.configure(text=text[:120],style='Bad.TLabel')
        if self.exit_on_finish:
            # Unattended mode: nobody can dismiss a dialog, so fail loudly instead.
            print('OASIS4.3: '+text,file=sys.stderr,flush=True)
            if not self.controller.running():
                self.root.after_idle(self.root.destroy)
            return
        from tkinter import messagebox
        messagebox.showerror('OASIS4.3',text)

    def run(self):
        self.root.mainloop()


class MetricChart:
    """Draw R-work / R-free per refinement round on a canvas."""

    WORK='#2b6cb0'
    FREE='#c53030'

    def __init__(self,app,parent):
        tk,ttk=app.tk,app.ttk
        self.rows=[]
        self.frame=ttk.Frame(parent)
        self.canvas=tk.Canvas(self.frame,height=190,background='#ffffff',
                              highlightthickness=1,highlightbackground='#dde3ea')
        self.canvas.pack(fill='x')
        self.canvas.bind('<Configure>',lambda event:self._draw())

    def show(self,rows):
        """Show the rows from :func:`phasebuild.charts.round_series`."""
        self.rows=list(rows or [])
        self._draw()

    def _draw(self):
        canvas=self.canvas
        canvas.delete('all')
        width=max(canvas.winfo_width(),160);height=max(canvas.winfo_height(),120)
        rows=self.rows
        if not rows:
            canvas.create_text(width/2,height/2,text='尚无逐轮精修数据',fill='#8a97a5')
            return
        if len(rows)<2:
            canvas.create_text(width/2,height/2,fill='#8a97a5',
                               text='只有 %d 轮数据；至少两轮才能画出曲线'%len(rows))
            return
        left,right,top,bottom=52,14,14,26
        low,high=charts.y_range(rows)
        ticks=charts.y_ticks(low,high)
        if ticks:
            low=min(low,ticks[0]);high=max(high,ticks[-1])
        plot_w=max(width-left-right,10);plot_h=max(height-top-bottom,10)

        def px(index):
            return left+plot_w*index/(len(rows)-1)

        def py(value):
            return top+plot_h*(high-value)/max(high-low,1e-9)

        for tick in ticks:
            y=py(tick)
            canvas.create_line(left,y,left+plot_w,y,fill='#e6ebf1')
            canvas.create_text(left-6,y,anchor='e',text='%.0f%%'%tick,fill='#7b8794',
                               font=('Helvetica',9))
        canvas.create_line(left,top,left,top+plot_h,fill='#9aa7b4')
        canvas.create_line(left,top+plot_h,left+plot_w,top+plot_h,fill='#9aa7b4')
        every=1 if len(rows)<=12 else 2
        for index,row in enumerate(rows):
            if index%every:
                continue
            canvas.create_text(px(index),top+plot_h+12,text=str(row['round']),fill='#7b8794',
                               font=('Helvetica',9))

        legend=(('r_work',self.WORK,None,'R-work'),('r_free',self.FREE,(4,3),'R-free'))
        for key,color,dash,name in legend:
            points=[(px(index),py(row[key])) for index,row in enumerate(rows)
                    if row[key] is not None]
            if len(points)>=2:
                flat=[value for point in points for value in point]
                canvas.create_line(*flat,fill=color,width=2,dash=dash)
            for x,y in points:
                canvas.create_oval(x-2.5,y-2.5,x+2.5,y+2.5,fill=color,outline='')

        last=[row for row in rows if row['r_work'] is not None or row['r_free'] is not None][-1]
        legend_w=112;legend_h=18*len(legend)+8
        x0=left+plot_w-legend_w;y0=top+6
        canvas.create_rectangle(x0-6,y0-2,x0+legend_w,y0+legend_h,fill='#ffffff',
                                outline='#e2e8f0')
        for index,(key,color,dash,name) in enumerate(legend):
            y=y0+6+18*index
            canvas.create_line(x0,y,x0+18,y,fill=color,width=2,dash=dash)
            value=last[key]
            canvas.create_text(x0+24,y,anchor='w',fill='#3d4852',font=('Helvetica',9),
                               text='%s %s'%(name,'—' if value is None else '%.1f%%'%value))


class StructureView:
    """Draw the built model on a canvas: drag to rotate, wheel to zoom.

    The Cα trace comes from the result preview, so the common case needs no
    file parsing; “显示所有原子” reads the model with gemmi on demand.
    """

    def __init__(self, app, parent):
        tk,ttk=app.tk,app.ttk
        self.app=app
        self.path=None
        self.preview=[]
        self.chains=[]
        self.atoms=None
        self.yaw=0.6
        self.pitch=0.35
        self.zoom=1.0
        self._drag=None
        self.frame=ttk.Frame(parent)
        self.frame.pack(fill='both',expand=True)
        tools=ttk.Frame(self.frame);tools.pack(fill='x',pady=(4,2))
        self.label=ttk.Label(tools,text='尚无模型',style='Hint.TLabel')
        self.label.pack(side='left')
        ttk.Button(tools,text='重置视角',command=self.reset).pack(side='right',padx=(6,0))
        self.all_atoms=tk.BooleanVar(master=app.root,value=False)
        ttk.Checkbutton(tools,text='显示所有原子',variable=self.all_atoms,
                        command=self._draw).pack(side='right')
        self.canvas=tk.Canvas(self.frame,background='#0f1216',highlightthickness=0)
        self.canvas.pack(fill='both',expand=True)
        self.canvas.bind('<ButtonPress-1>',self._press)
        self.canvas.bind('<B1-Motion>',self._move)
        self.canvas.bind('<ButtonRelease-1>',lambda event:setattr(self,'_drag',None))
        self.canvas.bind('<MouseWheel>',self._wheel)
        self.canvas.bind('<Configure>',lambda event:self._draw())

    def show(self, path, preview=None):
        """Display ``path`` (or clear the page when it is None)."""
        self.path=str(path) if path else None
        self.preview=preview or []
        self.atoms=None
        self.chains=[]
        if self.path:
            from . import structure_view
            self.chains=structure_view.chain_traces(self.preview)
            if not self.chains:
                try:
                    self.chains=[(name,coordinates) for name,coordinates,_ in
                                 structure_view.read_atoms(self.path)]
                except Exception:
                    self.chains=[]
            residues=sum(len(points) for _,points in self.chains)
            self.label.configure(text='%s · %d 条链 · %d 个残基'%(Path(self.path).name,
                                 len(self.chains),residues))
        else:
            self.label.configure(text='尚无模型')
        self.all_atoms.set(False)
        self.reset()

    def reset(self):
        self.yaw=0.6;self.pitch=0.35;self.zoom=1.0
        self._draw()

    def _press(self, event):
        self._drag=(event.x,event.y)

    def _move(self, event):
        if not self._drag:
            return
        x,y=self._drag
        self.yaw+=(event.x-x)*0.012
        self.pitch=max(-1.5,min(1.5,self.pitch+(event.y-y)*0.012))
        self._drag=(event.x,event.y)
        self._draw()

    def _wheel(self, event):
        step=1.1 if event.delta>0 else 1/1.1
        self.zoom=max(0.2,min(8.0,self.zoom*step))
        self._draw()
        return 'break'  # keep the form's global wheel handler out of the page

    def _atom_groups(self):
        if self.atoms is None:
            from . import structure_view
            try:
                self.atoms=[(name,coordinates) for name,coordinates,_ in
                            structure_view.read_atoms(self.path)]
            except Exception as error:
                self.atoms=[]
                self.label.configure(text='读取模型失败：'+str(error))
        return self.atoms

    def _draw(self):
        import numpy as np
        from . import structure_view
        canvas=self.canvas
        canvas.delete('all')
        width=max(canvas.winfo_width(),10);height=max(canvas.winfo_height(),10)
        if not self.path:
            canvas.create_text(width/2,height/2,text='运行完成后这里显示最终模型',
                               fill='#7f8c99')
            return
        all_atoms=self.all_atoms.get()
        groups=self._atom_groups() if all_atoms else self.chains
        points=[coordinates for _,coordinates in groups if len(coordinates)]
        if not points:
            canvas.create_text(width/2,height/2,text='模型中没有可绘制的原子',fill='#7f8c99')
            return
        matrix=structure_view.rotation(self.yaw,self.pitch)
        center,scale=structure_view.fit(np.vstack(points),width,height)
        scale*=self.zoom
        total=0
        for index,(name,coordinates) in enumerate(groups):
            if not len(coordinates):
                continue
            color=structure_view.PALETTE[index%len(structure_view.PALETTE)]
            x,y,_=structure_view.project(coordinates,matrix,center,scale,width,height)
            total+=len(coordinates)
            if all_atoms:
                for px,py in zip(x,y):
                    canvas.create_oval(px-1.4,py-1.4,px+1.4,py+1.4,fill=color,outline='')
            else:
                flat=[]
                for px,py in zip(x,y):
                    flat.extend((px,py))
                if len(flat)>=4:
                    canvas.create_line(*flat,fill=color,width=2)
        canvas.create_text(10,height-10,anchor='sw',fill='#8fa0b0',
                           text='%s · %d 个%s · 拖动旋转 / 滚轮缩放'%(
                               Path(self.path).name,total,'原子' if all_atoms else 'Cα'))


def _terminate(signum,frame):
    raise KeyboardInterrupt


def resolve_project(path):
    """Return the directory holding the pipeline scripts, or None.

    A source checkout keeps them at the root; a portable bundle keeps them in
    ``<bundle>/app`` while the runtime and the launcher sit one level above.
    """
    path=Path(path).resolve()
    for candidate in (path,path/'app'):
        if (candidate/'run_sad_to_model.py').is_file():
            return candidate
    return None


def main(argv=None):
    parser=argparse.ArgumentParser(description='OASIS4.3 桌面界面')
    parser.add_argument('--project',default=str(Path(__file__).resolve().parents[2]),
                        help='包含 run_sad_to_model.py 的源码目录')
    parser.add_argument('--config',help='预填表单的 JSON 配置')
    parser.add_argument('--autorun',action='store_true',help='打开后立即开始运行')
    parser.add_argument('--exit-on-finish',action='store_true',help='任务结束后自动关闭窗口')
    parser.add_argument('--check',action='store_true',help='只构建窗口并自检，不进入事件循环')
    args=parser.parse_args(argv)
    project=resolve_project(args.project)
    if project is None:
        parser.error('请选择包含 run_sad_to_model.py 的源码目录（--project）')
    try:
        import numpy,gemmi,scipy  # noqa: F401
    except ImportError as error:
        parser.error('缺少科学计算依赖，请使用完整便携包启动：'+str(error))
    config=None
    if args.config:
        config=json.loads(Path(args.config).read_text())
    controller=Controller(project)
    app=App(controller,config=config,autorun=args.autorun,exit_on_finish=args.exit_on_finish)
    if args.check:
        app.root.update_idletasks()
        app.root.update()
        app.root.destroy()
        print('desktop window ok')
        return 0
    # A closed terminal (SIGHUP) or an external kill must still stop the engine
    # children, so turn both into the same shutdown path as closing the window.
    for name in ('SIGTERM','SIGHUP'):
        if hasattr(signal,name):
            try:signal.signal(getattr(signal,name),_terminate)
            except (ValueError,OSError):pass
    try:
        app.run()
    except KeyboardInterrupt:
        pass
    finally:
        controller.shutdown()
    return 1 if app.error_text else 0


if __name__=='__main__':
    sys.exit(main())
