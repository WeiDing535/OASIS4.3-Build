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

STATUS_NAMES={'queued':'queued','running':'running','stopping':'stopping','cancelled':'stopped',
              'failed':'failed','success':'completed','imported':'imported (view only)','interrupted':'interface interrupted'}
IMPORTED_NOTE='imported existing results: read-only directory, will not be rerun; you can view the curve/table, export results, and save the model.'
METHODS=[('pplus_bundled','P+ external engine (recommended)'),('pplus_trust','Phaser SAD + external P+'),
         ('pplus_internal','internal P+ (experimental)'),('likelihood','continuous SAD'),('pplus','legacy direct-method P+')]
METHOD_NAMES=dict(METHODS)
PROTOCOLS=[('','no parameter recipe'),('pplus_w8_weight','P+ · Wang8 + refinement weight search (calibrated recipe)'),
           ('pplus_w8','P+ · Wang8'),('pplus_w14','P+ · Wang14'),('stable_w14','stabilized P+ · Wang14')]
RECIPES=[('pipeline','CCP4 Buccaneer/REFMAC pipeline (recommended)'),('ccp4','external CCP4 pipeline'),('legacy','legacy per-cycle pipeline')]
ENGINES=[('parrot','CCP4 Parrot density modification (external engine)'),
         ('handdm','handdm density modification (in-house engine, no CCP4 needed)')]
REFINE_ENGINES=[('refmac','CCP4 REFMAC refinement (external engine)'),
                ('minirefmac','minirefmac refinement (in-house engine, no CCP4 needed)')]
PRIORS=[('dm_phifom','density-modification phases + model amplitudes'),('original_sad','original SAD phases')]

# key, title, kind, hint, options, only-methods
FIELD_SPECS=[
 ('sad','SAD anomalous diffraction MTZ','file','requires Friedel-matched amplitudes and their errors; existing phases are not required.',None,None),
 ('heavy','heavy-atom substructure PDB / mmCIF','file','SHELXD substructure for this dataset; external P+ can determine the hand.',None,None),
 ('sequence','protein sequence FASTA / PIR','file','standard amino acids; multiple entries are treated as separate chains.',None,None),
 ('native','original observed MTZ (optional)','file','leave blank to use the SAD file and keep or create a free set.',None,None),
 ('elements','heavy-atom elements (required)','text','enter element names only, e.g. Se or Se,S; f″ is computed automatically.',None,None),
 ('wavelength','experimental SAD wavelength (Å, required)','float','actual acquisition wavelength; the program does not estimate it.',None,None),
 ('observation_prior','Phaser SAD likelihood-only MTZ','file','required for P+ optimization; downstream map or model phases are not allowed.',None,('pplus_trust',)),
 ('composition_sequence','explicit element composition sequence (optional)','file','write SeMet as (MSE); leave blank to use the model-building sequence.',None,('pplus_trust','pplus_internal','pplus_bundled')),
 ('output','new output directory','output','the directory must not exist yet; all results are saved there.',None,None),
 ('copies','number of copies of the whole sequence file in the ASU','int','positive integer 1–100, or auto to determine from the FASTA header / Matthews (default); not the copy number of a single chain.',None,None),
 ('solvent_fraction','solvent fraction (leave blank to estimate)','optional','between 0 and 1; leave blank to estimate from sequence and cell.',None,None),
 ('phase_cycles','phasing cycles','int','1–100.',None,None),
 ('dm_cycles','density-modification cycles','int','1–100; the calibrated recipe uses 40.',None,None),
 ('rounds','Buccaneer / REFMAC rounds','int','1–20; default 5 rounds.',None,None),
 ('cycles','REFMAC refinement cycles per round','int','1–200; default 10 cycles per round.',None,None),
 ('trust_wang_radius','single-pipeline Wang mask radius (Å)','float','used only in single-pipeline mode; fixed at 14 / 8 Å in multi-pipeline mode.',None,('pplus_trust','pplus_bundled')),
 ('trust_kl','P+ optimization: KL cap per reflection','float','default 0.0005; must be greater than 0 and at most 0.05.',None,('pplus_trust','pplus_internal','pplus_bundled')),
 ('max_per_target','relation cap per target','int','default 48; 0 means no truncation.',None,None),
 ('seed','random seed','int','0–2147483647; for reproducible random steps.',None,None),
 ('labels','SAD column mapping (leave blank to detect automatically)','text','mean values: f=FP,sigf=SIGFP,dano=DANO,sigdano=SIGDANO.',None,None),
 ('free_value','original free-set value','int','default 0; must match the original FreeR_flag definition.',None,None),
 ('reference','density-modification reference MTZ','file','a general unrelated reference provided by the user-supplied CCP4 run directory.',None,None),
 ('runtime','CCP4 installation root','directory','user-supplied CCP4: the directory must contain bin/phaser, bin/cparrot, bin/cbuccaneer, bin/refmac5, and lib/data. Click "detect" to find a standard installation automatically; the command-line equivalent is --ccp4-root.',None,None),
 ('feedback_endpoints','feedback endpoint set (blank = single endpoint)','optional','comma-separated feedback rounds, e.g. 2,3: each value runs the same protocol once, then the endpoint is chosen by work-set R (R-free is not read).',None,None),
 ('dm_solvents','feedback endpoint solvent set (optional)','optional','paired with the previous item, comma-separated; auto estimates from sequence, e.g. auto,0.5; leave blank to keep the solvent fraction.',None,None),
]

# What an experiment needs is on the form; the calibrated constants, the run
# resources and the rare overrides wait behind "show advanced parameters".
ESSENTIAL_GROUPS=[('input files',('sad','heavy','sequence','elements','wavelength','output')),
                  ('sample and density modification',('copies','solvent_fraction')),
                  ('model building and refinement',('rounds','cycles'))]
ADVANCED_GROUPS=[('input details',('native','observation_prior','composition_sequence','labels','free_value')),
                 ('cycles and sampling',('phase_cycles','dm_cycles','max_per_target','seed')),
                 ('P+ optimization',('trust_wang_radius','trust_kl')),
                 ('external engine paths and endpoint sets',('runtime','reference','feedback_endpoints','dm_solvents'))]


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

    def detect_ccp4(self, value=None):
        """Return (root, message) for a best-effort host CCP4 discovery."""
        from .engine_paths import detect_ccp4
        return detect_ccp4(value)

    def plan(self, config):
        return self.workflow.make_plan(config,self.project)

    def describe(self, config):
        return self.workflow.describe_plan(config)

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
        self.root.title('OASIS4.3 · single-wavelength anomalous diffraction phasing')
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
        self._detect_ccp4()
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
        ttk.Label(bar,text='single-wavelength anomalous diffraction → density modification → Buccaneer → REFMAC',
                  style='Hint.TLabel').pack(side='left',padx=(10,0),pady=(4,0))
        self.status=ttk.Label(bar,text='ready',style='Status.TLabel')
        self.status.pack(side='right')

        self.pages=ttk.Notebook(self.root)
        self.pages.pack(fill='both',expand=True,padx=10,pady=(0,6))

        inputs=ttk.Frame(self.pages)
        canvas=tk.Canvas(inputs,highlightthickness=0)
        scroll=ttk.Scrollbar(inputs,orient='vertical',command=canvas.yview)
        canvas.configure(yscrollcommand=scroll.set)
        scroll.pack(side='right',fill='y');canvas.pack(side='left',fill='both',expand=True)
        form=ttk.Frame(canvas,padding=(0,0,8,8))
        window=canvas.create_window((0,0),window=form,anchor='nw')
        form.bind('<Configure>',lambda e:canvas.configure(scrollregion=canvas.bbox('all')))
        canvas.bind('<Configure>',lambda e:canvas.itemconfigure(window,width=e.width))
        canvas.bind_all('<MouseWheel>',lambda e:canvas.yview_scroll(int(-e.delta/40),'units')
                        if self.pages.index(self.pages.select())==0 else None)
        self.form=form
        self.pages.add(inputs,text='  configuration inputs  ')

        outputs=ttk.Frame(self.pages)
        self.tabs=ttk.Notebook(outputs)
        self.tabs.pack(fill='both',expand=True)
        self._history_tab();self._log_tab();self._results_tab();self._structure_tab();self._plan_tab()
        self.pages.add(outputs,text='  run and results  ')

        foot=ttk.Frame(self.root,padding=(10,0,10,8));foot.pack(side='bottom',fill='x')
        self.progress=ttk.Progressbar(foot,mode='determinate',maximum=100)
        self.progress.pack(side='left',fill='x',expand=True,padx=(0,10))
        self.new_button=ttk.Button(foot,text='new job',command=self.on_new)
        self.new_button.pack(side='left')
        self.start_button=ttk.Button(foot,text='start',command=self.on_start)
        self.start_button.pack(side='left',padx=(6,0))
        self.stop_button=ttk.Button(foot,text='stop',command=self.on_stop,state='disabled')
        self.stop_button.pack(side='left',padx=6)
        ttk.Button(foot,text='open result directory',command=self.on_open).pack(side='left')
        self.export_button=ttk.Button(foot,text='export results…',command=self.on_export)
        self.export_button.pack(side='left',padx=6)
        ttk.Button(foot,text='import existing results…',command=self.on_import).pack(side='left')

        self._build_form()

    def _history_tab(self):
        tk,ttk=self.tk,self.ttk
        frame=ttk.Frame(self.tabs)
        tools=ttk.Frame(frame);tools.pack(fill='x',pady=(4,2))
        ttk.Label(tools,text='select a row to view it in the pages on the right; double-click to open the result directory.',
                  style='Hint.TLabel').pack(side='left')
        ttk.Button(tools,text='refresh',command=lambda:self._refresh_history(force=True)).pack(side='right')
        self.delete_button=ttk.Button(tools,text='delete record',command=self.on_delete_history)
        self.delete_button.pack(side='right',padx=6)
        self.relaunch_button=ttk.Button(tools,text='reuse parameters',command=self.on_relaunch)
        self.relaunch_button.pack(side='right')
        columns=('started','name','status','method','rwork','output')
        self.history=ttk.Treeview(frame,columns=columns,show='headings',height=10)
        for key,title,width in [('started','start time',118),('name','job name',150),('status','status',96),
                                ('method','method',140),('rwork','final R-work',96),('output','output directory',250)]:
            self.history.heading(key,text=title)
            self.history.column(key,width=width,anchor='center' if key in ('started','status','rwork') else 'w')
        y=ttk.Scrollbar(frame,orient='vertical',command=self.history.yview)
        self.history.configure(yscrollcommand=y.set)
        y.pack(side='right',fill='y');self.history.pack(fill='both',expand=True)
        self.history.bind('<<TreeviewSelect>>',self.on_select_job)
        self.history.bind('<Double-1>',lambda event:self.on_open())
        self.tabs.add(frame,text='job records')

    def _log_tab(self):
        tk,ttk=self.tk,self.ttk
        frame=ttk.Frame(self.tabs)
        tools=ttk.Frame(frame);tools.pack(fill='x',pady=(4,2))
        ttk.Checkbutton(tools,text='auto-scroll',variable=self.follow).pack(side='left')
        ttk.Button(tools,text='clear display',command=lambda:(self.log.delete('1.0','end'))).pack(side='right')
        self.log=tk.Text(frame,wrap='none',height=20,font=('Menlo',11),background='#101418',
                         foreground='#d8e0e8',insertbackground='#d8e0e8')
        y=ttk.Scrollbar(frame,orient='vertical',command=self.log.yview)
        x=ttk.Scrollbar(frame,orient='horizontal',command=self.log.xview)
        self.log.configure(yscrollcommand=y.set,xscrollcommand=x.set)
        y.pack(side='right',fill='y');x.pack(side='bottom',fill='x')
        self.log.pack(fill='both',expand=True)
        self.tabs.add(frame,text='run log')

    def _results_tab(self):
        ttk=self.ttk
        frame=ttk.Frame(self.tabs)
        tools=ttk.Frame(frame);tools.pack(fill='x',pady=(4,2))
        ttk.Label(tools,text='R-work / R-free are percentages; the last row is the fixed final refinement.',
                  style='Hint.TLabel').pack(side='left')
        self.save_model_button=ttk.Button(tools,text='save final model…',command=self.on_save_model,
                                          state='disabled')
        self.save_model_button.pack(side='right')
        ttk.Button(tools,text='export results…',command=self.on_export).pack(side='right',padx=6)
        self.summary=self.ttk.Label(frame,text='no results yet',justify='left',anchor='w',padding=(4,6))
        self.summary.pack(fill='x')
        self.chart=MetricChart(self,frame)
        self.chart.frame.pack(fill='x',pady=(2,6))
        columns=('stage','engine','rwork','rfree','cycles')
        self.tree=ttk.Treeview(frame,columns=columns,show='headings',height=10)
        for key,title,width in [('stage','stage',140),('engine','engine',220),('rwork','R-work',90),
                                ('rfree','R-free',90),('cycles','refinement cycles',80)]:
            self.tree.heading(key,text=title)
            self.tree.column(key,width=width,anchor='center' if key!='engine' else 'w')
        self.tree.pack(fill='both',expand=True)
        self.tabs.add(frame,text='refinement results')

    def _structure_tab(self):
        frame=self.ttk.Frame(self.tabs)
        self.structure=StructureView(self,frame)
        self.tabs.add(frame,text='Structure 3D')

    def _plan_tab(self):
        tk,ttk=self.tk,self.ttk
        frame=ttk.Frame(self.tabs)
        self.plan_text=tk.Text(frame,wrap='word',height=18,font=('Menlo',11))
        y=ttk.Scrollbar(frame,orient='vertical',command=self.plan_text.yview)
        self.plan_text.configure(yscrollcommand=y.set)
        y.pack(side='right',fill='y');self.plan_text.pack(fill='both',expand=True)
        self.plan_text.insert('end','after filling in the inputs, click "check configuration"; the stages to be executed and the effective parameters are shown here.')
        self.plan_text.configure(state='disabled')
        self.tabs.add(frame,text='run plan')

    # ------------------------------------------------------------------ form
    def _section(self,title,parent=None):
        box=self.ttk.LabelFrame(parent or self.form,text=title,padding=(10,6))
        box.pack(fill='x',pady=(0,8))
        return box

    def _build_form(self):
        ttk=self.ttk
        specs={spec[0]:spec for spec in FIELD_SPECS}
        task=self._section('job name')
        row=ttk.Frame(task);row.pack(fill='x',pady=2)
        ttk.Label(row,text='job name',width=26,anchor='w').pack(side='left')
        self.task_name=self.tk.StringVar(value='')
        ttk.Entry(row,textvariable=self.task_name).pack(side='left',fill='x',expand=True)
        ttk.Label(task,text='optional: give this job a name (e.g. "7E1D rerun"); it is shown in the job records so runs '
                            'with the same parameters but different data can be told apart at a glance.',
                  style='Hint.TLabel',wraplength=520,justify='left').pack(fill='x',padx=(2,0))
        box=self._section('phasing method')
        self.method=self.tk.StringVar(value=self.defaults.get('method','pplus_bundled'))
        for key,label in METHODS:
            ttk.Radiobutton(box,text=label,value=key,variable=self.method,
                            command=self._method_changed).pack(anchor='w')
        row=ttk.Frame(box);row.pack(fill='x',pady=(4,0))
        ttk.Label(row,text='P+ parameter recipe').pack(side='left')
        self.protocol=self.tk.StringVar(value='')
        ttk.Combobox(row,textvariable=self.protocol,values=[label for _,label in PROTOCOLS],
                     state='readonly',width=34).pack(side='left',padx=6)

        inputs=self._section('input files')
        for key in ESSENTIAL_GROUPS[0][1]:
            self._field(inputs,specs[key])

        sample=self._section('sample and density modification')
        for key in ESSENTIAL_GROUPS[1][1]:
            self._field(sample,specs[key])
        row=ttk.Frame(sample);row.pack(fill='x',pady=2)
        ttk.Label(row,text='density-modification engine',width=26,anchor='w').pack(side='left')
        self.dm_engine=self.tk.StringVar(value='parrot')
        ttk.Combobox(row,textvariable=self.dm_engine,values=[label for _,label in ENGINES],
                     state='readonly',width=30).pack(side='left')

        build=self._section('model building and refinement')
        row=ttk.Frame(build);row.pack(fill='x',pady=2)
        self.build_model=self.tk.BooleanVar(value=True)
        ttk.Checkbutton(row,text='Buccaneer model building + REFMAC refinement',variable=self.build_model).pack(side='left')
        self.auto_hand=self.tk.BooleanVar(value=True)
        ttk.Checkbutton(row,text='automatic hand determination (try both hands, without R-free)',variable=self.auto_hand).pack(side='left',padx=12)
        for key in ESSENTIAL_GROUPS[2][1]:
            self._field(build,specs[key])
        row=ttk.Frame(build);row.pack(fill='x',pady=2)
        ttk.Label(row,text='refinement engine',width=26,anchor='w').pack(side='left')
        self.refine_engine=self.tk.StringVar(value='refmac')
        ttk.Combobox(row,textvariable=self.refine_engine,values=[label for _,label in REFINE_ENGINES],
                     state='readonly',width=30).pack(side='left')
        for key,title,values,var in [('building_recipe','build recipe',RECIPES,'recipe'),
                                     ('refinement_prior','refinement phase source',PRIORS,'prior')]:
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

        # One-click engine discovery: a complete standard CCP4 install can be
        # detected and used to prefill the runtime/reference paths.
        if 'runtime' in self.fields:
            runtime_row=self.fields['runtime']['row']
            ttk.Button(runtime_row,text='detect',width=6,
                       command=lambda:self._detect_ccp4(say=True)).pack(side='right',padx=(0,4))
            self.ccp4_status=ttk.Label(runtime_row.master,text='',style='Hint.TLabel',
                                       wraplength=520,justify='left')
            self.ccp4_status.pack(fill='x',padx=(2,0))

        actions=ttk.Frame(self.form);actions.pack(fill='x',pady=(0,6))
        self.advanced=self.tk.BooleanVar(value=False)
        ttk.Checkbutton(actions,text='show advanced parameters (usually keep defaults)',variable=self.advanced,
                        command=self._toggle_advanced).pack(side='left')
        ttk.Button(actions,text='check configuration',command=self.on_check).pack(side='right')
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
            # Right-packed first: the browse button stays visible even when the
            # window is narrower than the entry's requested width.
            ttk.Button(row,text='browse…',width=6,
                       command=lambda:self._browse(key,var,kind)).pack(side='right',padx=(4,0))
            ttk.Entry(row,textvariable=var).pack(side='left',fill='x',expand=True)
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
            path=filedialog.askdirectory(initialdir=current,title='select directory')
        elif kind=='output':
            path=filedialog.askdirectory(initialdir=current,title='select the parent directory of the new output directory')
            if path:
                path=str(Path(path)/('oasis_'+_stamp()))
        else:
            path=filedialog.askopenfilename(initialdir=str(Path(current).parent),
                                            title='select file')
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

    def _detect_ccp4(self,say=False):
        """Detect a complete CCP4 install and prefill runtime/reference paths."""
        root,message=self.controller.detect_ccp4()
        if root:
            runtime=self.fields['runtime']['var']
            if not runtime.get().strip():
                runtime.set(root)
            reference=self.fields['reference']['var']
            if not reference.get().strip():
                reference.set(str(Path(root)/'lib/data/reference_structures/reference-1tqw.mtz'))
        if hasattr(self,'ccp4_status'):
            self.ccp4_status.configure(text=message)
        if say:
            self._say(message,ok=bool(root))
        return bool(root)

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
        refine=dict(REFINE_ENGINES).get(values.get('refine_engine'))
        if refine:
            self.refine_engine.set(refine)
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
                    raise ValueError('please provide: '+field['title'])
                if key=='copies' and text.lower()=='auto':
                    # ``auto`` is the documented default: FASTA header first,
                    # then the smallest Matthews-plausible number of copies.
                    values[key]='auto'
                    continue
                try:
                    values[key]=int(text) if kind=='int' else float(text)
                except ValueError:
                    raise ValueError('"%s" must be a number; currently: %s'%(field['title'],text))
            elif kind=='optional':
                values[key]=text
            else:
                values[key]=text
        if not str(values.get('solvent_fraction') or '').strip():
            # Blank means "estimate from sequence and cell": None, not ''.
            values['solvent_fraction']=None
        values['method']=self.method.get()
        values['dm_engine']=dict((label,key) for key,label in ENGINES).get(self.dm_engine.get(),'parrot')
        values['refine_engine']=dict((label,key) for key,label in REFINE_ENGINES).get(self.refine_engine.get(),'refmac')
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
        lines=[(stage.get('name',''),stage.get('status','planned')) for stage in plan.get('stages',[])]
        text='\n'.join('· %s  [%s]'%(name,status) for name,status in lines)
        effective=plan.get('effective_options') or {}
        body='\n'.join('%s: %s'%(key,value) for key,value in effective.items())
        self._set_plan(text+('\n\neffective parameters\n'+body if body else ''))
        self._say('configuration check passed',ok=True)
        self.pages.select(1)
        self.tabs.select(4)

    def on_start(self):
        """Start one more run; earlier runs keep going in parallel."""
        try:
            config=self.collect()
            name=self.task_name.get().strip()
            job=self.controller.start(config,name)
        except Exception as error:
            self._error(str(error));return
        self.task_name.set('')
        self.job=job
        self.log_text=''
        self.log.delete('1.0','end')
        self._clear_results()
        self.summary.configure(text='running: '+str(job['output']))
        self.stop_button.configure(state='normal')
        self.progress.configure(value=0)
        self.pages.select(1)
        self.tabs.select(1)
        self._refresh_history(force=True)
        self._select_history(job['id'])
        self._say('started job %s%s (you can start more jobs in parallel)'
                  %(job.get('name') or job['id'],'' if not job.get('name') else ' (%s)'%job['id']))
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
            self._say('this job has finished');return
        self.controller.stop(self.job['id'])
        self._say('stopping '+self.job['id']+' …')

    def on_open(self):
        if not self.job:
            self._error('no result directory yet.');return
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
        self._method_changed()
        self.progress.configure(value=0)
        self.stop_button.configure(state='disabled')
        for item in self.history.selection():
            self.history.selection_remove(item)
        self.summary.configure(text='new job: fill in the inputs, then click "start"')
        self._say('created new job: output directory '+self.fields['output']['var'].get())
        self.pages.select(0)

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
            self._error('no results available to export.');return
        source=Path(self.job['output'])
        if not source.is_dir():
            self._error('result directory does not exist: '+str(source));return
        from tkinter import filedialog
        target=filedialog.askdirectory(initialdir=str(Path.home()),title='select export location')
        if not target:return
        destination=Path(target)/('oasis_export_'+_stamp())
        try:
            shutil.copytree(source,destination)
        except OSError as error:
            self._error('export failed: '+str(error));return
        files=[item for item in destination.rglob('*') if item.is_file()]
        size=sum(item.stat().st_size for item in files)
        self._say('exported %d files (%.1f MB) to %s'%(len(files),size/1e6,destination),ok=True)

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
            self._error('no final model yet.');return
        from tkinter import filedialog
        path=filedialog.asksaveasfilename(defaultextension='.cif',initialfile=Path(model).name,
                                          initialdir=str(Path(model).parent),title='save final model')
        if not path:return
        try:
            shutil.copyfile(model,path)
        except OSError as error:
            self._error('save failed: '+str(error));return
        self._say('model saved to '+path,ok=True)

    def on_import(self):
        from tkinter import filedialog
        path=filedialog.askdirectory(initialdir=str(self.controller.project),title='select an existing result directory')
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
            question='%d jobs are still running; quitting will stop them. Quit anyway?'%len(live)
            if not messagebox.askokcancel('quit OASIS4.3',question):
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
        log=job.get('log') or ('(imported existing results have no run log)' if status=='imported' else '')
        if log!=self.log_text:
            self._show_log(log)
        if status in ('success','failed','cancelled','interrupted'):
            self.finished.add(identifier)
            self.stop_button.configure(state='disabled')
            self.progress.configure(value=100 if status=='success' else 0)
            if status=='success':
                self._show_results(select=False)
            elif status=='failed':
                self.error_text=job.get('error') or 'job failed; check the log'
                self.summary.configure(text='failed: '+self.error_text)
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
        text=label or 'ready'
        if len(live)>1:
            text+=' · %d jobs running'%len(live)
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
            self._error('please select a row in the job records, then click "reuse parameters".');return
        identifier=selection[0]
        try:
            job=self.controller.state(identifier)
        except KeyError:
            self._error('this record no longer exists.');return
        config=job.get('config')
        if not config:
            self._error('this record has no saved run parameters (imported existing results) and cannot be reused.');return
        self._fill(config,strict=False)
        self.task_name.set(job.get('name') or '')
        self.fields['output']['var'].set(str(self._fresh_output()))
        self._method_changed()
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
        self.pages.select(0)
        self.summary.configure(text='reused parameters from job "%s": update the input data, then click "start"; the output directory has been changed to a new one.'
                                    %(job.get('name') or job.get('label') or identifier))
        self._say('parameters reused; you can start the run now')

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
            self._error('please select the row to delete in the job records.');return
        busy={'queued','running','stopping'}
        targets=[i for i in selection if self.history.set(i,'status') not in busy]
        live=len(selection)-len(targets)
        if not targets:
            self._error('the selected job is still running; please stop it before deleting the record.');return
        from tkinter import messagebox
        note='%d records will be deleted from the job records; result directories and files will not be deleted.'%len(targets)
        if live: note+='\n%d running records will be kept.'%live
        if not messagebox.askokcancel('delete job records',note):
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
        self._say('deleted %d job records (result directories unchanged)'%removed)

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
        self._show_log(job.get('log') or ('(imported existing results have no run log)' if status=='imported' else ''))
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
            self.summary.configure(text=STATUS_NAMES.get(status,status)+': '+str(job.get('error') or 'check the log'))
        else:
            self._clear_results()
            self.summary.configure(text='running: '+str(job.get('output')))
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
        self.summary.configure(text='no results yet')
        self.save_model_button.configure(state='disabled')
        self.chart.show([])
        self.structure.show(None)

    def _show_results(self,select=True):
        try:
            result=self.controller.results(self.job['id'])
        except Exception as error:
            self.summary.configure(text='failed to read results: '+str(error));return
        for row in self.tree.get_children():
            self.tree.delete(row)
        for row in result.get('metrics') or []:
            rwork=row.get('r_work');rfree=row.get('r_free')
            self.tree.insert('','end',values=(row.get('stage',''),
                row.get('engine',''),'' if rwork is None else '%.2f'%(100*rwork),
                '' if rfree is None else '%.2f'%(100*rfree),row.get('cycles','')))
        primary=result.get('primary_metric') or {}
        lines=['method: %s  density modification: %s'%(result.get('method_name',''),result.get('dm_engine',''))]
        if result.get('solvent_fraction') is not None:
            lines.append('solvent fraction: %.2f (%s)'%(result['solvent_fraction'],result.get('solvent_source') or 'estimated'))
        if primary:
            lines.append('final R-work %.2f / R-free %.2f (%s)'%(100*primary['r_work'],
                         100*primary['r_free'],primary.get('stage','')))
        endpoints=result.get('endpoints') or []
        if endpoints:
            parts=[]
            for item in endpoints:
                metric=item.get('primary_metric') or {}
                rwork=metric.get('r_work')
                solvent=item.get('solvent')
                parts.append('fb%s/solvent %s → %s'%(item.get('feedback'),'auto' if solvent is None else '%.2f'%solvent,
                             '—' if rwork is None else 'R-work %.2f'%(100*rwork)))
            selected=(result.get('endpoint_selection') or {}).get('selected_feedback')
            lines.append('endpoint set (fb%s kept by work-set R): %s'%(selected,'; '.join(parts)))
        if result.get('final_model'):
            lines.append('model: '+result['final_model'])
        elif result.get('final_phases'):
            lines.append('phases: '+result['final_phases'])
        self.summary.configure(text='\n'.join(lines))
        self.result=result
        self.chart.show(charts.round_series(result.get('metrics')))
        model=self._model_path(result)
        self.save_model_button.configure(state='normal' if model and model.is_file() else 'disabled')
        self.structure.show(str(model) if model and model.is_file() else None,result.get('preview'))
        if select:
            self.pages.select(1)
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
            canvas.create_text(width/2,height/2,text='no per-cycle refinement data yet',fill='#8a97a5')
            return
        if len(rows)<2:
            canvas.create_text(width/2,height/2,fill='#8a97a5',
                               text='only %d cycle of data; at least two cycles are needed to draw a curve'%len(rows))
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
    file parsing; "show all atoms" reads the model with gemmi on demand.
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
        self.label=ttk.Label(tools,text='no model yet',style='Hint.TLabel')
        self.label.pack(side='left')
        ttk.Button(tools,text='reset view',command=self.reset).pack(side='right',padx=(6,0))
        self.all_atoms=tk.BooleanVar(master=app.root,value=False)
        ttk.Checkbutton(tools,text='show all atoms',variable=self.all_atoms,
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
            if self.chains:
                residues=sum(len(points) for _,points in self.chains)
                chain_count=len(self.chains)
            else:
                # No preview (e.g. imported result): draw all atoms but report
                # real chain/residue counts, not the atom count as residues.
                try:
                    self.chains=[(name,coordinates) for name,coordinates,_ in
                                 structure_view.read_atoms(self.path)]
                    chain_count,residues=structure_view.count_chains_residues(self.path)
                except Exception:
                    self.chains=[]
                    chain_count=residues=0
            self.label.configure(text='%s · %d chains · %d residues'%(Path(self.path).name,
                                 chain_count,residues))
        else:
            self.label.configure(text='no model yet')
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
                self.label.configure(text='failed to read model: '+str(error))
        return self.atoms

    def _draw(self):
        import numpy as np
        from . import structure_view
        canvas=self.canvas
        canvas.delete('all')
        width=max(canvas.winfo_width(),10);height=max(canvas.winfo_height(),10)
        if not self.path:
            canvas.create_text(width/2,height/2,text='the final model is shown here after the run completes',
                               fill='#7f8c99')
            return
        all_atoms=self.all_atoms.get()
        groups=self._atom_groups() if all_atoms else self.chains
        points=[coordinates for _,coordinates in groups if len(coordinates)]
        if not points:
            canvas.create_text(width/2,height/2,text='no drawable atoms in the model',fill='#7f8c99')
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
                           text='%s · %d %s · drag to rotate / scroll to zoom'%(
                               Path(self.path).name,total,'atoms' if all_atoms else 'Cα'))


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
    parser=argparse.ArgumentParser(description='OASIS4.3 desktop interface')
    parser.add_argument('--project',default=str(Path(__file__).resolve().parents[2]),
                        help='source directory containing run_sad_to_model.py')
    parser.add_argument('--config',help='JSON configuration to prefill the form')
    parser.add_argument('--autorun',action='store_true',help='start the run immediately after opening')
    parser.add_argument('--exit-on-finish',action='store_true',help='close the window automatically when the job finishes')
    parser.add_argument('--check',action='store_true',help='build the window and self-check without entering the event loop')
    args=parser.parse_args(argv)
    project=resolve_project(args.project)
    if project is None:
        parser.error('please select a source directory containing run_sad_to_model.py (--project)')
    try:
        import numpy,gemmi,scipy  # noqa: F401
    except ImportError as error:
        parser.error('missing scientific dependencies (numpy/gemmi/scipy): %s. Please run with a Python 3.12 environment that has the dependencies, '
                     'e.g. activate the oasis43-unified conda environment before launching, or specify an interpreter with PHASEBUILD_PYTHON.'%error)
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
