"""Desktop front end: toolkit-free controller, form mapping and Tk window."""
import json
from pathlib import Path
import pytest
from phasebuild import desktop,oasis_workflow as workflow
from test_sad_pipeline import inputs


def project(tmp_path):
    root=Path(__file__).resolve().parents[1]
    (tmp_path/'run_sad_to_model.py').symlink_to(root/'run_sad_to_model.py')
    return tmp_path


def test_controller_defaults_and_plan_come_from_the_workflow(tmp_path):
    controller=desktop.Controller(project(tmp_path))
    default=controller.defaults()
    assert default['method']=='pplus_bundled' and default['dm_engine']=='parrot'
    assert 'profile' not in default
    assert 'gui_workflow' not in __import__('sys').modules


def test_unverified_hand_switch_reaches_the_job_command(tmp_path,monkeypatch,inputs):
    """The opt-in must be a real argv flag, not a silently dropped form value."""
    from phasebuild import bundled_pplus
    monkeypatch.setattr(bundled_pplus,'engine_layout',lambda runtime,dm_engine=None,need_density=False:None)
    monkeypatch.setattr(bundled_pplus,'require_density',lambda cfg:Path(cfg['runtime'])/'density')
    root=project(tmp_path)
    base=dict(inputs,auto_hand=True,method='pplus_bundled',wavelength=2.0,
              output=str(tmp_path/'run'))
    for value,flag in [(False,'--no-allow-unverified-hand'),(True,'--allow-unverified-hand')]:
        cfg,stages=workflow.make_plan(dict(base,allow_unverified_hand=value),root,inspect=False)
        assert cfg['allow_unverified_hand'] is value
        assert flag in stages[0]['command']


def test_controller_rejects_a_retired_mode_without_creating_output(tmp_path):
    controller=desktop.Controller(project(tmp_path))
    with pytest.raises(ValueError):
        controller.plan({'mode':'refine'})
    assert not (tmp_path/'runs').exists()


def test_controller_import_requires_a_real_result_tree(tmp_path):
    controller=desktop.Controller(project(tmp_path))
    empty=tmp_path/'empty';empty.mkdir()
    with pytest.raises(ValueError):
        controller.import_results(empty)


def test_every_form_field_is_a_workflow_option():
    """A stale widget would make make_plan reject the whole configuration."""
    root=Path(__file__).resolve().parents[1]
    default=workflow.defaults(root)
    covered={spec[0] for spec in desktop.FIELD_SPECS}
    assert covered <= set(default)
    assert {'sad','heavy','sequence','elements','wavelength','copies','output'} <= covered


def test_every_method_carries_a_hint_and_the_internal_route_is_flagged():
    """The experimental chain must be labelled in the form, not only in the docs."""
    assert set(desktop.METHOD_HINTS)=={key for key,_ in desktop.METHODS}
    assert '实验性' in desktop.METHOD_NAMES['pplus_internal']
    assert '实验性' in desktop.METHOD_HINTS['pplus_internal']
    assert '推荐' in desktop.METHOD_NAMES['pplus_bundled']
    assert 'Parrot' in desktop.METHOD_HINTS['pplus_bundled']


def test_protocol_and_engine_labels_round_trip():
    for key,label in desktop.PROTOCOLS:
        assert dict((label,value) for value,label in desktop.PROTOCOLS)[label]==key
    assert dict(desktop.ENGINES)['parrot']=='CCP4 Parrot（推荐）'


def fill(app,values):
    for key,value in values.items():
        if key in app.fields:
            app.fields[key]['var'].set(str(value))


@pytest.mark.skipif(not __import__('importlib').util.find_spec('tkinter'),reason='Tk not installed')
def test_window_builds_and_collects_the_form(tmp_path):
    """Exercise the real widgets without entering the event loop."""
    import tkinter
    try:
        tkinter.Tk().destroy()
    except tkinter.TclError:
        pytest.skip('no display available')
    controller=desktop.Controller(project(tmp_path))
    app=desktop.App(controller)
    try:
        app.root.update_idletasks()
        assert app.start_button.winfo_exists()
        with pytest.raises(ValueError,match='波长'):
            app.collect()
        fill(app,dict(sad='sad.mtz',heavy='heavy.pdb',sequence='seq.fasta',elements='Se',
                      wavelength='0.97918',copies='2',output=str(tmp_path/'run')))
        values=app.collect()
        assert values['method']=='pplus_bundled' and values['dm_engine']=='parrot'
        assert values['build_model'] is True and values['auto_hand'] is True
        assert values['allow_unverified_hand'] is False   # default off: marginal evidence still fails closed
        app.unverified_hand.set(True)
        assert app.collect()['allow_unverified_hand'] is True
        app.unverified_hand.set(False)
        assert values['copies']==2 and values['wavelength']==0.97918
        assert values['protocol']=='pplus_w8_weight' and values['building_recipe']=='pipeline'
        assert values['solvent_fraction'] is None  # blank means estimate, not ''
        assert app.method_hint.get()==desktop.METHOD_HINTS['pplus_bundled']
        app.method.set('pplus_internal')
        app._method_changed()
        assert '实验性' in app.method_hint.get()
    finally:
        app.root.destroy()


def test_copy_number_accepts_auto_or_an_integer(tmp_path):
    """``auto`` is the documented default, so the form must not force a number."""
    pytest.importorskip('tkinter')
    import tkinter
    try:
        tkinter.Tk().destroy()
    except tkinter.TclError:
        pytest.skip('no display available')
    app=desktop.App(desktop.Controller(project(tmp_path)))
    try:
        fill(app,dict(sad='sad.mtz',heavy='heavy.pdb',sequence='seq.fasta',elements='Se',
                      wavelength='0.97918',copies='auto',output=str(tmp_path/'run')))
        assert app.collect()['copies']=='auto'
        app.fields['copies']['var'].set('3')
        assert app.collect()['copies']==3
    finally:
        app.root.destroy()


def test_main_check_builds_the_window_without_an_event_loop(tmp_path,capsys):
    pytest.importorskip('tkinter')
    import tkinter
    try:
        tkinter.Tk().destroy()
    except tkinter.TclError:
        pytest.skip('no display available')
    assert desktop.main(['--project',str(project(tmp_path)),'--check'])==0
    assert 'desktop window ok' in capsys.readouterr().out


def test_unattended_failure_is_reported_instead_of_waiting_for_a_dialog(tmp_path,capsys):
    """--exit-on-finish runs cannot click a message box, so errors must go to stderr."""
    pytest.importorskip('tkinter')
    import tkinter
    try:
        tkinter.Tk().destroy()
    except tkinter.TclError:
        pytest.skip('no display available')
    app=desktop.App(desktop.Controller(project(tmp_path)),exit_on_finish=True)
    try:
        app._error('boom')
        assert app.error_text=='boom'
        assert 'OASIS4.3: boom' in capsys.readouterr().err
    finally:
        app.root.destroy()


def test_check_button_shows_a_plan_that_has_no_status_yet(tmp_path,monkeypatch):
    """A fresh plan has no run status; the check must still render it."""
    pytest.importorskip('tkinter')
    import tkinter
    try:
        tkinter.Tk().destroy()
    except tkinter.TclError:
        pytest.skip('no display available')
    app=desktop.App(desktop.Controller(project(tmp_path)),exit_on_finish=True)
    try:
        monkeypatch.setattr(app,'collect',lambda: {'output':str(tmp_path/'out')})
        monkeypatch.setattr(app.controller,'describe',lambda config: dict(
            stages=[dict(id='dm',name='密度修饰',command=None)],
            effective_options={'dm_engine':'内置 RESOLVE/MLHL'},
            sample_warnings=['引擎已回退']))
        app.on_check()
        app.root.update_idletasks()
        assert '密度修饰' in app.plan_text.get('1.0','end')
        assert app.warning.winfo_ismapped()
        assert app.warning.cget('text')=='引擎已回退'
    finally:
        app.root.destroy()


def test_fallback_warning_is_a_red_line_that_clears(tmp_path):
    """An engine fallback must be visible in the window, not only in the log."""
    pytest.importorskip('tkinter')
    import tkinter
    try:
        tkinter.Tk().destroy()
    except tkinter.TclError:
        pytest.skip('no display available')
    app=desktop.App(desktop.Controller(project(tmp_path)))
    try:
        app.root.update_idletasks()
        assert not app.warning.winfo_ismapped()
        notice=app.controller.fallback_notice({'dm_engine_fallback':'Parrot (cparrot) was not found'})
        assert '已回退' in notice
        app._warn(notice)
        app.root.update_idletasks()
        assert app.warning.winfo_ismapped() and app.warning.cget('text')==notice
        assert app.warning.cget('style')=='Bad.TLabel'
        app._warn('')
        app.root.update_idletasks()
        assert not app.warning.winfo_ismapped()
        assert app.controller.fallback_notice({'dm_engine':'resolve'})==''
    finally:
        app.root.destroy()


def test_project_resolution_accepts_a_source_tree_or_a_portable_bundle(tmp_path):
    """The window must find its scripts from a checkout and from a bundle root."""
    (tmp_path/'run_sad_to_model.py').write_text('')
    assert desktop.resolve_project(tmp_path)==tmp_path.resolve()
    bundle=tmp_path/'bundle'
    (bundle/'app').mkdir(parents=True)
    (bundle/'app'/'run_sad_to_model.py').write_text('')
    assert desktop.resolve_project(bundle)==(bundle/'app').resolve()
    assert desktop.resolve_project(bundle/'app')==(bundle/'app').resolve()
    (bundle/'app'/'run_sad_to_model.py').unlink()
    assert desktop.resolve_project(bundle) is None


def test_results_panel_lists_the_endpoint_grid_and_the_retained_one(tmp_path):
    """A fixed-endpoint run must show every endpoint and which one was kept."""
    pytest.importorskip('tkinter')
    import tkinter
    try:
        tkinter.Tk().destroy()
    except tkinter.TclError:
        pytest.skip('no display available')
    controller=desktop.Controller(project(tmp_path))
    controller.results=lambda identifier: dict(method_name='P+ Wang8 + 精修权重搜索',dm_engine='parrot',metrics=[],
                solvent_fraction=0.5,solvent_source='user',final_model='00_sad_to_model/endpoint_fb2_s0.50/primary_model.cif',
                primary_metric={'stage':'完整链后精修','r_work':0.2834,'r_free':0.3147},
                endpoint_selection={'selected_feedback':2},
                endpoints=[{'feedback':2,'solvent':0.5,'primary_metric':{'r_work':0.2834,'r_free':0.3147}},
                           {'feedback':3,'solvent':0.5,'primary_metric':{'r_work':0.2911,'r_free':0.3273}}])
    app=desktop.App(controller)
    try:
        app.job={'id':'x','output':str(tmp_path)}
        app._show_results()
        text=app.summary.cget('text')
        assert '保留 fb2' in text and 'fb3/溶剂0.50' in text and 'R-work 29.11' in text
        assert '最终 R-work 28.34 / R-free 31.47' in text
    finally:
        app.root.destroy()


@pytest.mark.skipif(not __import__('importlib').util.find_spec('tkinter'),reason='Tk not installed')
def test_task_list_shows_every_job_and_switches_the_log(tmp_path):
    """The window must remember earlier runs and show the selected one's log."""
    import tkinter
    try:
        tkinter.Tk().destroy()
    except tkinter.TclError:
        pytest.skip('no display available')
    controller=desktop.Controller(project(tmp_path))
    rows=[dict(id='new',created=200,updated=210,status='running',output=str(tmp_path/'b'),
               completed=1,total=6,error=None,method='pplus_bundled',label='b.mtz',
               method_name='P+ 自带引擎',r_work=None,r_free=None),
          dict(id='old',created=100,updated=110,status='success',output=str(tmp_path/'a'),
               completed=6,total=6,error=None,method='pplus_bundled',label='a.mtz',
               method_name='P+ Wang8 + 精修权重搜索',r_work=0.2834,r_free=0.3147)]
    controller.history=lambda:rows
    controller.running=lambda:['new']
    controller.state=lambda identifier:dict(rows[0] if identifier=='new' else rows[1],
                                            log='log of '+identifier,pipeline_progress={'total':6,'completed':1})
    app=desktop.App(controller)
    try:
        app._refresh_history(force=True)
        assert set(app.history.get_children())=={'new','old'}
        assert app.history.set('old','rwork')=='28.34'
        assert app.history.set('old','status')=='已完成'
        app.history.selection_set('old')
        app.on_select_job()
        assert app.job['id']=='old' and 'log of old' in app.log.get('1.0','end')
        assert app.tabs.index(app.tabs.select())==0   # the record page keeps its buttons usable
        assert app.history.set('old','method')=='P+ Wang8 + 精修权重搜索'
    finally:
        app.root.destroy()


@pytest.mark.skipif(not __import__('importlib').util.find_spec('tkinter'),reason='Tk not installed')
def test_imported_result_is_shown_read_only(tmp_path):
    """An imported run must show its results instead of looking like a live job."""
    import tkinter
    try:
        tkinter.Tk().destroy()
    except tkinter.TclError:
        pytest.skip('no display available')
    controller=desktop.Controller(project(tmp_path))
    output=tmp_path/'portfolio';output.mkdir()
    rows=[dict(id='imported',created=100,updated=110,status='imported',output=str(output),
               completed=0,total=0,error=None,method='pplus_bundled',label='',
               method_name='P+ 自带引擎',r_work=0.3156,r_free=0.3829)]
    controller.history=lambda:rows
    controller.running=lambda:[]
    controller.state=lambda identifier:dict(rows[0],log='')
    pairs=[(0.3903,0.4266),(0.3156,0.3829)]
    metrics=[dict(stage='第 %d 轮'%index,round=index,r_work=work,r_free=free)
             for index,(work,free) in enumerate(pairs,1)]
    controller.results=lambda identifier:dict(metrics=metrics,final_model=None,preview=[],
                                              method_name='P+ 自带引擎',dm_engine='parrot',
                                              primary_metric=dict(stage='第 2 轮',r_work=0.3156,r_free=0.3829))
    app=desktop.App(controller)
    try:
        app._refresh_history(force=True)
        assert app.history.set('imported','status')=='已导入（仅查看）'
        assert app.history.set('imported','rwork')=='31.56'
        app.history.selection_set('imported')
        app.on_select_job()
        summary=app.summary.cget('text')
        assert '不会重新运行' in summary and '最终 R-work 31.56' in summary
        assert '运行中' not in summary
        assert str(app.stop_button['state'])=='disabled'
        assert len(app.tree.get_children())==2
        assert '没有运行日志' in app.log.get('1.0','end')
        assert len(app.chart.rows)==2
        app._refresh()   # the 900 ms poll must not wipe the placeholder
        assert '没有运行日志' in app.log.get('1.0','end')
        assert '不会重新运行' in app.summary.cget('text')
    finally:
        app.root.destroy()


def test_two_runs_start_from_the_window_without_stopping_the_first(tmp_path):
    """Starting a second run must work while the first is still refining."""
    import tkinter
    try:
        tkinter.Tk().destroy()
    except tkinter.TclError:
        pytest.skip('no display available')
    controller=desktop.Controller(project(tmp_path))
    started=[];rows=[]
    def start(config,name=''):
        identifier='job%d'%len(started)
        started.append((identifier,config,name))
        Path(config['output']).mkdir(parents=True,exist_ok=True)
        row=dict(id=identifier,created=100+len(started),updated=100+len(started),status='running',
                 output=str(config['output']),completed=0,total=6,error=None,
                 method=config['method'],label='sad.mtz',method_name='P+ 自带引擎',r_work=None,r_free=None)
        rows.insert(0,row)
        return dict(row,config=config,plan={'stages':[]},steps=[])
    controller.start=start
    controller.history=lambda:rows
    controller.running=lambda:[row['id'] for row in rows if row['status']=='running']
    app=desktop.App(controller)
    try:
        fill(app,dict(sad='sad.mtz',heavy='heavy.pdb',sequence='seq.fasta',elements='Se',
                      wavelength='0.97918',copies='2',output=str(tmp_path/'first')))
        app.task_name.set('第一个任务')
        app.on_start()
        assert app.job['id']=='job0'
        second=app.fields['output']['var'].get()
        assert second!=str(tmp_path/'first')  # a fresh directory is offered automatically
        fill(app,dict(output=str(tmp_path/'second')))
        app.on_start()
        assert [(identifier,name) for identifier,_,name in started]==[('job0','第一个任务'),('job1','')]
        assert len(app.history.get_children())==2
        assert app.job['id']=='job1'
        assert app.controller.running()==['job1','job0']
    finally:
        app.root.destroy()


@pytest.mark.skipif(not __import__('importlib').util.find_spec('tkinter'),reason='Tk not installed')
def test_closing_asks_once_for_all_running_jobs(tmp_path,monkeypatch):
    import tkinter
    try:
        tkinter.Tk().destroy()
    except tkinter.TclError:
        pytest.skip('no display available')
    controller=desktop.Controller(project(tmp_path))
    controller.running=lambda:['a','b']
    from tkinter import messagebox
    asked=[]
    monkeypatch.setattr(messagebox,'askokcancel',lambda *a,**k:(asked.append(a),False)[1])
    app=desktop.App(controller)
    try:
        app.on_close()
        assert asked and '2 个任务' in asked[0][1]
        assert app.root.winfo_exists()
    finally:
        app.root.destroy()


def window(tmp_path):
    """Build a real window for a test, skipping when no display exists."""
    pytest.importorskip('tkinter')
    import tkinter
    try:
        tkinter.Tk().destroy()
    except tkinter.TclError:
        pytest.skip('no display available')
    app=desktop.App(desktop.Controller(project(tmp_path)))
    app.root.update_idletasks()
    return app


def test_advanced_parameters_stay_hidden_until_asked_for(tmp_path):
    """Only experiment inputs are visible; calibrated constants wait behind a toggle."""
    app=window(tmp_path)
    try:
        assert app.advanced_box.winfo_manager()==''      # never packed at start
        essential=[key for _,keys in desktop.ESSENTIAL_GROUPS for key in keys]
        advanced=[key for _,keys in desktop.ADVANCED_GROUPS for key in keys]
        assert set(essential)|set(advanced)|set(desktop.MODULE_FIELDS)==set(app.fields)
        assert app.fields['sad']['row'].winfo_ismapped()
        assert not app.fields['native']['row'].winfo_ismapped()
        app.advanced.set(True);app._toggle_advanced();app.root.update_idletasks()
        assert app.advanced_box.winfo_manager()=='pack'
        assert app.fields['native']['row'].winfo_ismapped()
        app.advanced.set(False);app._toggle_advanced();app.root.update_idletasks()
        assert app.advanced_box.winfo_manager()==''
    finally:
        app.root.destroy()


def test_new_task_clears_the_form_and_offers_a_fresh_output(tmp_path):
    """“新建任务” must start from defaults without losing the task history."""
    app=window(tmp_path)
    try:
        fill(app,dict(sad='old.mtz',heavy='old.pdb',sequence='old.fasta',elements='Se',
                      wavelength='0.9',copies='3',output=str(tmp_path/'runs/oasis_run')))
        app.job={'id':'previous','output':str(tmp_path/'runs/oasis_run')}
        app.log.insert('end','old log')
        app.on_new()
        assert app.job is None
        assert app.fields['sad']['var'].get()=='' and app.fields['heavy']['var'].get()==''
        assert app.fields['wavelength']['var'].get()==''
        assert app.log.get('1.0','end').strip()==''
        output=app.fields['output']['var'].get()
        assert output!=str(tmp_path/'runs/oasis_run')
        assert Path(output).name.startswith('oasis_') and not Path(output).exists()
        assert Path(output).parent==tmp_path/'runs'
        assert '新任务' in app.summary.cget('text')
    finally:
        app.root.destroy()


def test_task_name_and_delete_live_in_the_job_record(tmp_path):
    """A name is stored with the record; deleting forgets it but keeps the results."""
    class Workflow:
        def make_plan(self,config,project): return dict(config),[]
        def describe_plan(self,config): return {'stages':[]}
        def save_selected_result(self,output): pass
        def enrich_job(self,job): return job
        def collect_results(self,output): raise ValueError('no results here')
    jobs=desktop.Jobs(project(tmp_path),workflow=Workflow())
    first=jobs.start(dict(output=str(tmp_path/'run1'),sad=str(tmp_path/'a.mtz'),
                          method='pplus_bundled'),'7E1D 重跑')
    jobs.workers[first['id']].join(5)
    assert jobs.history()[0]['name']=='7E1D 重跑'
    stored=json.loads((jobs.storage/(first['id']+'.json')).read_text())
    assert stored['name']=='7E1D 重跑'
    jobs.delete(first['id'])
    assert jobs.history()==[]
    assert not (jobs.storage/(first['id']+'.json')).exists()
    assert (tmp_path/'run1').is_dir()          # results stay on disk
    second=jobs.start(dict(output=str(tmp_path/'run2'),sad=str(tmp_path/'b.mtz')),'')
    jobs.workers[second['id']].join(5)
    assert jobs.history()[0]['name']==''
    jobs.items[second['id']]['status']='running'
    with pytest.raises(ValueError,match='还在运行'):
        jobs.delete(second['id'])
    assert (jobs.storage/(second['id']+'.json')).exists()


@pytest.mark.skipif(not __import__('importlib').util.find_spec('tkinter'),reason='Tk not installed')
def test_window_names_relaunches_and_deletes_history_rows(tmp_path,monkeypatch):
    """Name a run, reuse its parameters, and remove records without losing data."""
    import tkinter
    try:
        tkinter.Tk().destroy()
    except tkinter.TclError:
        pytest.skip('no display available')
    import tkinter.messagebox as box
    controller=desktop.Controller(project(tmp_path))
    config=dict(controller.defaults(),sad=str(tmp_path/'old.mtz'),heavy=str(tmp_path/'old.pdb'),
                sequence=str(tmp_path/'old.fasta'),wavelength=0.97918,copies=2,
                output=str(tmp_path/'old_run'),seed=7)
    rows=[dict(id='done',created=100,updated=110,status='success',output=str(tmp_path/'old_run'),
               name='7E1D 重跑',completed=6,total=6,error=None,method='pplus_bundled',
               label='old.mtz',method_name='P+ 自带引擎',r_work=0.28,r_free=0.31),
          dict(id='imported',created=90,updated=91,status='imported',output=str(tmp_path/'imp'),
               name='',completed=0,total=0,error=None,method='pplus_bundled',label='',
               method_name='P+ 自带引擎',r_work=0.30,r_free=0.35),
          dict(id='busy',created=80,updated=81,status='running',output=str(tmp_path/'busy'),
               name='',completed=1,total=6,error=None,method='pplus_bundled',label='',
               method_name='P+ 自带引擎',r_work=None,r_free=None)]
    deleted=[]
    controller.history=lambda:list(rows)
    controller.state=lambda identifier:dict([r for r in rows if r['id']==identifier][0],log='',
                                            config=config if identifier=='done' else None)
    def delete(identifier):
        deleted.append(identifier)
        rows[:]=[r for r in rows if r['id']!=identifier]
    controller.delete=delete
    controller.results=lambda identifier:dict(metrics=[],final_model=None,preview=[],
                                              method_name='P+ 自带引擎',dm_engine='parrot',
                                              primary_metric={})
    errors=[]
    monkeypatch.setattr(box,'showerror',lambda title,text:errors.append(text))
    monkeypatch.setattr(box,'askokcancel',lambda *args,**kwargs: True)
    app=desktop.App(controller)
    try:
        app._refresh_history(force=True)
        assert app.history.set('done','name')=='7E1D 重跑'
        assert app.history.set('imported','name')==''
        app.history.selection_set('done')
        app.on_select_job()
        assert app.tabs.index(app.tabs.select())==0
        app.on_relaunch()
        assert app.fields['sad']['var'].get()==str(tmp_path/'old.mtz')
        assert app.fields['wavelength']['var'].get()=='0.97918'
        assert app.task_name.get()=='7E1D 重跑'
        assert app.advanced.get() is True            # seed was not the default
        output=app.fields['output']['var'].get()
        assert output!=str(tmp_path/'old_run') and Path(output).name.startswith('oasis_')
        assert app.tabs.index(app.tabs.select())==4  # the plan page shows the reused setup
        assert '已沿用任务' in app.summary.cget('text')
        app.history.selection_set('imported')
        app.on_select_job()
        app.on_relaunch()
        assert '无法沿用' in errors[-1]               # imported records carry no config
        app.history.selection_set('busy')
        app.on_delete_history()
        assert deleted==[] and '还在运行' in app.status.cget('text')
        app.history.selection_set('done')
        app.on_select_job()
        app.on_delete_history()
        assert deleted==['done']
        assert 'done' not in app.history.get_children()
        assert app.job is None and '尚无结果' in app.summary.cget('text')
    finally:
        app.root.destroy()


def test_export_copies_the_whole_run_directory(tmp_path,monkeypatch):
    """The download button copies the run tree to a folder the user picks."""
    app=window(tmp_path)
    try:
        run=tmp_path/'run'
        (run/'00_sad_to_model').mkdir(parents=True)
        (run/'00_sad_to_model'/'primary_model.cif').write_text('data_x\n')
        (run/'gui.log').write_text('log\n')
        app.job={'id':'x','output':str(run)}
        from tkinter import filedialog
        target=tmp_path/'chosen';target.mkdir()
        monkeypatch.setattr(filedialog,'askdirectory',lambda **kwargs:str(target))
        app.on_export()
        exported=list(target.glob('oasis_export_*'))
        assert len(exported)==1
        assert (exported[0]/'00_sad_to_model'/'primary_model.cif').read_text()=='data_x\n'
        assert (exported[0]/'gui.log').is_file()
        assert '已导出 2 个文件' in app.status.cget('text')
    finally:
        app.root.destroy()


def test_save_model_writes_the_final_model_file(tmp_path,monkeypatch):
    app=window(tmp_path)
    try:
        model=tmp_path/'primary_model.cif';model.write_text('data_x\n')
        app.result={'final_model':str(model)}
        from tkinter import filedialog
        destination=tmp_path/'saved.cif'
        monkeypatch.setattr(filedialog,'asksaveasfilename',lambda **kwargs:str(destination))
        app.on_save_model()
        assert destination.read_text()=='data_x\n'
    finally:
        app.root.destroy()


def test_structure_page_draws_the_ca_trace_and_can_show_all_atoms(tmp_path):
    """The 3D page draws the result preview and switches to all atoms on demand."""
    app=window(tmp_path)
    try:
        assert '结构 3D' in [app.tabs.tab(index,'text') for index in range(len(app.tabs.tabs()))]
        model=tmp_path/'primary_model.pdb'
        model.write_text(
            'ATOM      1  N   MET A   1       0.000   0.000   0.000  1.00 20.00           N\n'
            'ATOM      2  CA  MET A   1       1.458   0.000   0.000  1.00 20.00           C\n'
            'ATOM      3  CA  GLY A   2       2.916   0.000   0.000  1.00 20.00           C\n'
            'END\n')
        preview=[{'chain':'A','points':[[0,0,0,'1','MET'],[3.8,0,0,'2','ALA'],[7.6,0,0,'3','GLY']]},
                 {'chain':'B','points':[[0,4,0,'1','SER'],[3.8,4,0,'2','THR']]}]
        app.tabs.select(3);app.root.update()
        app.structure.show(model,preview)
        app.structure._draw()
        assert len(app.structure.chains)==2
        trace_items=len(app.structure.canvas.find_all())
        assert trace_items>=3            # two polylines plus the caption
        app.structure.all_atoms.set(True)
        app.structure._draw()
        assert len(app.structure.canvas.find_all())>trace_items
        assert sum(len(coordinates) for _,coordinates in app.structure.atoms)==3
        app.structure.all_atoms.set(False)
        app.structure.show(None)
        assert app.structure.canvas.find_all()      # a placeholder, not an empty page
    finally:
        app.root.destroy()


def test_results_page_draws_the_r_work_and_r_free_curve(tmp_path):
    """The curve needs the per-round metrics, and one round cannot make a line."""
    app=window(tmp_path)
    try:
        pairs=[(0.3903,0.4266),(0.3639,0.4031),(0.3402,0.3931),(0.3345,0.3869),(0.3156,0.3829)]
        metrics=[dict(stage='第 %d 轮'%index,round=index,r_work=work,r_free=free)
                 for index,(work,free) in enumerate(pairs,1)]
        app.controller.results=lambda identifier:dict(metrics=metrics,final_model=None,preview=[])
        app.job={'id':'x','output':str(tmp_path)}
        app._show_results()
        app.root.update()
        app.chart._draw()
        items=app.chart.canvas.find_all()
        assert len(items)>25
        texts=[app.chart.canvas.itemcget(item,'text') for item in items
               if app.chart.canvas.type(item)=='text']
        assert any('R-work 31.6%' in text for text in texts)
        assert any('R-free 38.3%' in text for text in texts)
        assert sum(1 for text in texts if text.endswith('%'))>=3   # y axis ticks
        assert {'1','5'} <= set(texts)                             # x axis rounds

        app.controller.results=lambda identifier:dict(metrics=metrics[:1],final_model=None,preview=[])
        app._show_results()
        app.chart._draw()
        texts=[app.chart.canvas.itemcget(item,'text') for item in app.chart.canvas.find_all()
               if app.chart.canvas.type(item)=='text']
        assert any('至少两轮' in text for text in texts)
    finally:
        app.root.destroy()


@pytest.mark.skipif(not __import__('importlib').util.find_spec('tkinter'),reason='Tk not installed')
def test_optional_modules_follow_ccp4_detection_and_gate_building(tmp_path,monkeypatch):
    """Parrot and Buccaneer default on when CCP4 is found, off with a notice when not."""
    import tkinter
    try:
        tkinter.Tk().destroy()
    except tkinter.TclError:
        pytest.skip('no display available')
    from phasebuild import module_detection
    def report(parrot,build):
        modules={'parrot':dict(title='Parrot 密度修饰',available=parrot,path=None,reason=''),
                 'build':dict(title='Buccaneer 建模 + REFMAC 精修',available=build,path=None,
                              reason='',refmac=build)}
        return dict(available=parrot or build,modules=modules,summary='',
                    ccp4_ready=parrot and build,searched=[])
    monkeypatch.setattr(module_detection,'detect',lambda runtime=None,environ=None,location=None: report(False,False))
    app=desktop.App(desktop.Controller(project(tmp_path)))
    try:
        app.root.update_idletasks()
        app.fields['wavelength']['var'].set('0.97918')
        app.fields['copies']['var'].set('1')
        assert app.module_parrot.get() is False and app.module_build.get() is False
        assert str(app.build_check.cget('state'))=='disabled'
        hint=app.module_hint.get()
        assert '未检测到 CCP4' in hint and '模块 1 输出相位' in hint
        values=app.collect()
        assert values['dm_enabled'] is False and values['build_model'] is False
        monkeypatch.setattr(module_detection,'detect',lambda runtime=None,environ=None,location=None: report(True,True))
        app._detect_modules()
        assert app.module_parrot.get() is True and app.module_build.get() is True
        assert '已检测到 CCP4' in app.module_hint.get()
        assert str(app.build_check.cget('state'))=='normal'
        values=app.collect()
        assert values['dm_enabled'] is True and values['build_model'] is True
        app.module_parrot.set(False);app._modules_changed()
        assert app.module_build.get() is False
        assert str(app.build_check.cget('state'))=='disabled'
        assert '模块 2 已关闭' in app.module_hint.get()
        assert str(app.fields['dm_cycles']['row'].winfo_children()[1].cget('state'))=='disabled'
        values=app.collect()
        assert values['dm_enabled'] is False and values['build_model'] is False
        app.module_parrot.set(True);app._modules_changed()
        assert str(app.fields['dm_cycles']['row'].winfo_children()[1].cget('state'))=='normal'
    finally:
        app.root.destroy()


@pytest.mark.skipif(not __import__('importlib').util.find_spec('tkinter'),reason='Tk not installed')
def test_a_typo_in_the_ccp4_location_falls_back_and_says_so(tmp_path,monkeypatch):
    """A path with no CCP4 must not be left in the form as the run directory."""
    import tkinter
    try:
        tkinter.Tk().destroy()
    except tkinter.TclError:
        pytest.skip('no display available')
    from phasebuild import module_detection
    calls=[]
    def fake(runtime=None,environ=None,location=None):
        calls.append((runtime,location))
        modules={'parrot':dict(title='Parrot 密度修饰',available=True,path='/x/cparrot',reason=''),
                 'build':dict(title='Buccaneer 建模 + REFMAC 精修',available=True,path='/x/cbuccaneer',
                              reason='',refmac='/x/refmac5')}
        return dict(available=True,modules=modules,summary='',ccp4_ready=True,searched=[],
                    location=dict(value=str(location or ''),resolved='',found_at='/opt/ccp4-9',note=''))
    monkeypatch.setattr(module_detection,'detect',fake)
    app=desktop.App(desktop.Controller(project(tmp_path)))
    try:
        bundled=str(app.defaults.get('runtime') or '')
        app.fields['runtime']['var'].set(str(tmp_path/'typo'))
        app._detect_modules()
        # the form no longer carries a directory that cannot run anything
        assert app.fields['runtime']['var'].get()==bundled
        assert '“%s”里没找到 CCP4 程序'%str(tmp_path/'typo') in app.module_hint.get()
        assert 'CCP4：/opt/ccp4-9' in app.module_hint.get()
        assert app.module_parrot.get() and app.module_build.get()
        # the retry asked about the runtime the run will actually use
        assert calls[-1]==(bundled or None,None)
    finally:
        app.root.destroy()


@pytest.mark.skipif(not __import__('importlib').util.find_spec('tkinter'),reason='Tk not installed')
def test_browse_buttons_stay_visible_and_ccp4_location_feeds_detection(tmp_path,monkeypatch):
    """An expanding entry must not push 浏览… off the row; a typed CCP4 path is used."""
    import tkinter
    try:
        tkinter.Tk().destroy()
    except tkinter.TclError:
        pytest.skip('no display available')
    from phasebuild import module_detection
    seen={}
    def fake(runtime=None,environ=None,location=None):
        seen['runtime'],seen['location']=runtime,location
        modules={'parrot':dict(title='Parrot 密度修饰',available=True,path='/x/cparrot',reason=''),
                 'build':dict(title='Buccaneer 建模 + REFMAC 精修',available=True,path='/x/cbuccaneer',
                              reason='',refmac='/x/refmac5')}
        return dict(available=True,modules=modules,summary='',ccp4_ready=True,searched=[],
                    location=dict(value=str(location or ''),resolved='/opt/ccp4-9',note=''))
    monkeypatch.setattr(module_detection,'detect',fake)
    app=desktop.App(desktop.Controller(project(tmp_path)))
    try:
        app.root.geometry('1180x1050')
        app.root.update_idletasks();app.root.update()
        for key in ('sad','heavy','sequence','output'):
            row=app.fields[key]['row']
            buttons=[child for child in row.winfo_children() if child.winfo_class()=='TButton']
            assert len(buttons)==1,key
            button=buttons[0]
            assert button.cget('text')=='浏览…'
            assert button.winfo_ismapped() and button.winfo_width()>30,key
            # the whole button, not just its left edge, must fit inside the row
            assert button.winfo_x()+button.winfo_width()<=row.winfo_width()+1,key
            entry=[child for child in row.winfo_children() if child.winfo_class()=='TEntry'][0]
            assert entry.winfo_width()>=80,key
        # a setup script resolves to the installation root, which is then used to run
        script=tmp_path/'ccp4.setup';script.write_text('export CCP4=/opt/ccp4-9\n')
        app.fields['runtime']['var'].set(str(script))
        app._detect_modules()
        assert seen['location']==str(script)
        assert app.fields['runtime']['var'].get()=='/opt/ccp4-9'
        assert 'CCP4：/opt/ccp4-9' in app.module_hint.get()
        # the same root is what the run uses, and both CCP4 modules default on
        app.fields['wavelength']['var'].set('0.97918')
        assert app.collect()['runtime']=='/opt/ccp4-9'
        assert app.module_parrot.get() and app.module_build.get()
    finally:
        app.root.destroy()
