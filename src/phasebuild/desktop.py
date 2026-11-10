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

from .jobs import Jobs

STATUS_NAMES={'queued':'等待运行','running':'正在运行','stopping':'正在停止','cancelled':'已停止',
              'failed':'失败','success':'已完成','imported':'已有结果','interrupted':'界面曾中断'}
METHODS=[('pplus_bundled','P+ 自带引擎（推荐）'),('pplus_trust','Phaser SAD + 外部 P+'),
         ('pplus_internal','内部 P+（实验性）'),('likelihood','连续 SAD'),('pplus','旧版直接法 P+')]
METHOD_NAMES=dict(METHODS)
PROTOCOLS=[('','不使用参数方案'),('pplus_w8_weight','P+ · Wang8＋精修权重搜索（标定配方）'),
           ('pplus_w8','P+ · Wang8'),('pplus_w14','P+ · Wang14'),('stable_w14','稳定化 P+ · Wang14')]
RECIPES=[('pipeline','自带 Buccaneer/REFMAC 流程（推荐）'),('ccp4','外部 CCP4 流程'),('legacy','旧版逐轮流程')]
ENGINES=[('parrot','CCP4 Parrot 密度修饰（推荐）'),('resolve','RESOLVE 密度修饰')]
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
 ('rounds','Buccaneer / REFMAC 轮数','int','1–20；默认 5 轮。',None,None),
 ('cycles','每轮 REFMAC 精修周期','int','1–200；默认每轮 10 周期。',None,None),
 ('trust_wang_radius','单流程 Wang 掩膜半径（Å）','float','仅单流程使用；多流程固定 14 / 8 Å。',None,('pplus_trust','pplus_bundled')),
 ('trust_kl','P+ 优化：每反射 KL 上限','float','默认 0.0005；必须大于 0 且不超过 0.05。',None,('pplus_trust','pplus_internal','pplus_bundled')),
 ('max_per_target','每个目标的关系上限','int','默认 48；0 表示不截断。',None,None),
 ('seed','随机种子','int','0–2147483647；用于可重复的随机步骤。',None,None),
 ('labels','SAD 列映射（留空自动识别）','text','平均值：f=FP,sigf=SIGFP,dano=DANO,sigdano=SIGDANO。',None,None),
 ('free_value','原自由集数值','int','默认 0；须与原始 FreeR_flag 定义一致。',None,None),
 ('reference','密度修饰参考 MTZ','file','随运行资源提供的通用无关参考。',None,None),
 ('runtime','Buccaneer / REFMAC 运行目录','directory','目录内需有 bin/cbuccaneer、bin/refmac5。',None,None),
 ('feedback_endpoints','反馈终点组合（留空=单终点）','optional','逗号分隔的反馈轮数，例如 2,3：每个值各跑一遍同一协议，再按工作集 R 选择终点（不读 R-free）。',None,None),
 ('dm_solvents','反馈终点溶剂组合（可选）','optional','与上一项配对，逗号分隔，auto 表示按序列估算，例如 auto,0.5；留空则沿用溶剂比例。',None,None),
]

# What an experiment needs is on the form; the calibrated constants, the run
# resources and the rare overrides wait behind “显示高级参数”.
ESSENTIAL_GROUPS=[('输入文件',('sad','heavy','sequence','elements','wavelength','output')),
                  ('样本与密度修饰',('copies','solvent_fraction')),
                  ('建模与精修',('rounds','cycles'))]
ADVANCED_GROUPS=[('输入细节',('native','observation_prior','composition_sequence','labels','free_value')),
                 ('循环与采样',('phase_cycles','dm_cycles','max_per_target','seed')),
                 ('P+ 优化',('trust_wang_radius','trust_kl')),
                 ('运行资源与终点组合',('reference','runtime','feedback_endpoints','dm_solvents'))]


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

    def start(self, config):
        return self.jobs.start(config)

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
        ttk.Label(bar,text='单波长异常散射 → 密度修饰 → Buccaneer → REFMAC',
                  style='Hint.TLabel').pack(side='left',padx=(10,0),pady=(4,0))
        self.status=ttk.Label(bar,text='就绪',style='Status.TLabel')
        self.status.pack(side='right')

        panes=ttk.Panedwindow(self.root,orient='horizontal')
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
        ttk.Label(tools,text='选中一行即可查看它的日志与结果；可以同时运行多个任务。',
                  style='Hint.TLabel').pack(side='left')
        ttk.Button(tools,text='刷新',command=lambda:self._refresh_history(force=True)).pack(side='right')
        columns=('started','status','method','rwork','output')
        self.history=ttk.Treeview(frame,columns=columns,show='headings',height=10)
        for key,title,width in [('started','开始时间',130),('status','状态',86),('method','方法',150),
                                ('rwork','最终 R-work',100),('output','输出目录',330)]:
            self.history.heading(key,text=title)
            self.history.column(key,width=width,anchor='center' if key in ('started','status','rwork') else 'w')
        y=ttk.Scrollbar(frame,orient='vertical',command=self.history.yview)
        self.history.configure(yscrollcommand=y.set)
        y.pack(side='right',fill='y');self.history.pack(fill='both',expand=True)
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
        columns=('stage','engine','rwork','rfree','cycles')
        self.tree=ttk.Treeview(frame,columns=columns,show='headings',height=10)
        for key,title,width in [('stage','阶段',140),('engine','引擎',220),('rwork','R-work',90),
                                ('rfree','R-free',90),('cycles','精修周期',80)]:
            self.tree.heading(key,text=title)
            self.tree.column(key,width=width,anchor='center' if key!='engine' else 'w')
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
        box=self._section('定相方法')
        self.method=self.tk.StringVar(value=self.defaults.get('method','pplus_bundled'))
        for key,label in METHODS:
            ttk.Radiobutton(box,text=label,value=key,variable=self.method,
                            command=self._method_changed).pack(anchor='w')
        row=ttk.Frame(box);row.pack(fill='x',pady=(4,0))
        ttk.Label(row,text='P+ 参数方案').pack(side='left')
        self.protocol=self.tk.StringVar(value='')
        ttk.Combobox(row,textvariable=self.protocol,values=[label for _,label in PROTOCOLS],
                     state='readonly',width=34).pack(side='left',padx=6)

        inputs=self._section('输入文件')
        for key in ESSENTIAL_GROUPS[0][1]:
            self._field(inputs,specs[key])

        sample=self._section('样本与密度修饰')
        for key in ESSENTIAL_GROUPS[1][1]:
            self._field(sample,specs[key])
        row=ttk.Frame(sample);row.pack(fill='x',pady=2)
        ttk.Label(row,text='密度修饰引擎',width=26,anchor='w').pack(side='left')
        self.dm_engine=self.tk.StringVar(value='parrot')
        ttk.Combobox(row,textvariable=self.dm_engine,values=[label for _,label in ENGINES],
                     state='readonly',width=30).pack(side='left')

        build=self._section('建模与精修')
        row=ttk.Frame(build);row.pack(fill='x',pady=2)
        self.build_model=self.tk.BooleanVar(value=True)
        ttk.Checkbutton(row,text='Buccaneer 建模 + REFMAC 精修',variable=self.build_model).pack(side='left')
        self.auto_hand=self.tk.BooleanVar(value=True)
        ttk.Checkbutton(row,text='自动判手（双手试建，不用 R-free）',variable=self.auto_hand).pack(side='left',padx=12)
        for key in ESSENTIAL_GROUPS[2][1]:
            self._field(build,specs[key])
        for key,title,values,var in [('building_recipe','构建配方',RECIPES,'recipe'),
                                     ('refinement_prior','精修相位来源',PRIORS,'prior')]:
            row=ttk.Frame(build);row.pack(fill='x',pady=2)
            ttk.Label(row,text=title,width=26,anchor='w').pack(side='left')
            variable=self.tk.StringVar(value=values[0][1])
            ttk.Combobox(row,textvariable=variable,values=[label for _,label in values],
                         state='readonly',width=30).pack(side='left')
            setattr(self,var,variable);setattr(self,var+'_values',values)

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
        row=ttk.Frame(parent);row.pack(fill='x',pady=2)
        ttk.Label(row,text=title,width=26,anchor='w').pack(side='left')
        var=self.tk.StringVar(value='')
        if kind=='file' or kind=='directory' or kind=='output':
            ttk.Entry(row,textvariable=var).pack(side='left',fill='x',expand=True)
            ttk.Button(row,text='浏览…',width=6,
                       command=lambda:self._browse(key,var,kind)).pack(side='left',padx=(4,0))
        elif kind=='text':
            ttk.Entry(row,textvariable=var).pack(side='left',fill='x',expand=True)
        else:
            ttk.Entry(row,textvariable=var,width=16).pack(side='left')
        ttk.Label(parent,text=hint,style='Hint.TLabel',wraplength=520,
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

    def _method_changed(self,*_):
        method=self.method.get()
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
        if 'build_model' in values:
            self.build_model.set(bool(values['build_model']))
        if 'auto_hand' in values:
            self.auto_hand.set(bool(values['auto_hand']))
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
        values['build_model']=bool(self.build_model.get())
        values['auto_hand']=bool(self.auto_hand.get())
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
        lines=[(stage['name'],stage['status']) for stage in plan.get('stages',[])]
        text='\n'.join('· %s  [%s]'%(name,status) for name,status in lines)
        effective=plan.get('effective_options') or {}
        body='\n'.join('%s：%s'%(key,value) for key,value in effective.items())
        self._set_plan(text+('\n\n生效参数\n'+body if body else ''))
        self._say('配置检查通过',ok=True)

    def on_start(self):
        """Start one more run; earlier runs keep going in parallel."""
        try:
            config=self.collect()
            job=self.controller.start(config)
        except Exception as error:
            self._error(str(error));return
        self.job=job
        self.log_text=''
        self.log.delete('1.0','end')
        self._clear_results()
        self.summary.configure(text='运行中：'+str(job['output']))
        self.stop_button.configure(state='normal')
        self.progress.configure(value=0)
        self.tabs.select(1)
        self._refresh_history(force=True)
        self._select_history(job['id'])
        self._say('已启动任务 '+job['id']+'（可继续新建任务并行运行）')
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
        self.fields['output']['var'].set(str(self._fresh_output()))
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
        log=job.get('log') or ''
        if log!=self.log_text:
            self._show_log(log)
        if status in ('success','failed','cancelled','interrupted'):
            self.finished.add(identifier)
            self.stop_button.configure(state='disabled')
            self.progress.configure(value=100 if status=='success' else 0)
            if status=='success':
                self._show_results()
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
                STATUS_NAMES.get(row['status'],row['status']),
                row.get('method_name') or METHOD_NAMES.get(row.get('method'),row.get('method') or ''),
                '' if rwork is None else '%.2f'%(100*rwork),
                row['output']))
        if self.job and self.history.exists(self.job['id']):
            self._select_history(self.job['id'])

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
        self._show_log(job.get('log') or '')
        status=job.get('status')
        self.stop_button.configure(state='normal' if status in ('queued','running','stopping') else 'disabled')
        self.progress.configure(value=0)
        self._show_plan_for(job)
        if status=='success':
            self.finished.add(identifier)
            self._show_results()
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
        self.structure.show(None)

    def _show_results(self):
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
        model=self._model_path(result)
        self.save_model_button.configure(state='normal' if model and model.is_file() else 'disabled')
        self.structure.show(str(model) if model and model.is_file() else None,result.get('preview'))
        self.tabs.select(2)

    def _set_plan(self,text):
        self.plan_text.configure(state='normal')
        self.plan_text.delete('1.0','end')
        self.plan_text.insert('end',text)
        self.plan_text.configure(state='disabled')

    def _say(self,text,ok=False):
        self.status.configure(text=text,style='Status.TLabel')

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
