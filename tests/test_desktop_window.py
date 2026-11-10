"""Desktop front end: toolkit-free controller, form mapping and Tk window."""
import json
from pathlib import Path
import pytest
from phasebuild import desktop,oasis_workflow as workflow


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


def test_protocol_and_engine_labels_round_trip():
    for key,label in desktop.PROTOCOLS:
        assert dict((label,value) for value,label in desktop.PROTOCOLS)[label]==key
    assert dict(desktop.ENGINES)['parrot']=='CCP4 Parrot 密度修饰（推荐）'


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
        assert values['copies']==2 and values['wavelength']==0.97918
        assert values['protocol']=='pplus_w8_weight' and values['building_recipe']=='pipeline'
        assert values['solvent_fraction'] is None  # blank means estimate, not ''
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
        assert app.history.set('old','method')=='P+ Wang8 + 精修权重搜索'
    finally:
        app.root.destroy()


@pytest.mark.skipif(not __import__('importlib').util.find_spec('tkinter'),reason='Tk not installed')
def test_two_runs_start_from_the_window_without_stopping_the_first(tmp_path):
    """Starting a second run must work while the first is still refining."""
    import tkinter
    try:
        tkinter.Tk().destroy()
    except tkinter.TclError:
        pytest.skip('no display available')
    controller=desktop.Controller(project(tmp_path))
    started=[];rows=[]
    def start(config):
        identifier='job%d'%len(started)
        started.append((identifier,config))
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
        app.on_start()
        assert app.job['id']=='job0'
        second=app.fields['output']['var'].get()
        assert second!=str(tmp_path/'first')  # a fresh directory is offered automatically
        fill(app,dict(output=str(tmp_path/'second')))
        app.on_start()
        assert [identifier for identifier,_ in started]==['job0','job1']
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
