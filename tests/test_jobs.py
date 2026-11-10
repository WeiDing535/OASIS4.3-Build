"""Task lifecycle of the toolkit-free job runner shared by every front end."""
import json
import sys
import time
from pathlib import Path
import pytest
from phasebuild import jobs as runner


def wait_for(predicate):
    limit=time.monotonic()+8
    while time.monotonic()<limit:
        if predicate():
            return
        time.sleep(.02)
    raise AssertionError('Task state did not settle')


def stub(monkeypatch,out,steps):
    monkeypatch.setattr(runner,'make_plan',lambda *a:({'output':str(out)},steps))
    monkeypatch.setattr(runner,'save_selected_result',lambda path:path)


def test_worker_failure_preserves_output_and_stops_pipeline(monkeypatch,tmp_path):
    out=tmp_path/'output'
    steps=[dict(name='first',command=[sys.executable,'-c',"print('visible log'); raise SystemExit(7)"],output=str(out/'one'),kind='build'),
           dict(name='must not run',command=[sys.executable,'-c',"raise SystemExit(0)"],output=str(out/'two'),kind='build')]
    stub(monkeypatch,out,steps)
    jobs=runner.Jobs(tmp_path);job=jobs.start({})
    wait_for(lambda:jobs.active is None)
    result=jobs.get(job['id'])
    assert result['status']=='failed' and result['completed']==0
    assert 'visible log' in result['log'] and '退出码 7' in result['error']
    assert json.loads((out/'gui_job.json').read_text())['status']=='failed'


def test_empty_trace_stops_with_actionable_error(monkeypatch,tmp_path):
    out=tmp_path/'output';trace=out/'trace'
    code='from pathlib import Path; p=Path('+repr(str(trace))+'); p.mkdir(); (p/"build_report.json").write_text(\'{"residues":0}\')'
    steps=[dict(name='trace',command=[sys.executable,'-c',code],output=str(trace),kind='trace')]
    stub(monkeypatch,out,steps)
    jobs=runner.Jobs(tmp_path);job=jobs.start({})
    wait_for(lambda:jobs.active is None)
    assert jobs.get(job['id'])['status']=='failed'
    assert '未找到满足筛选条件' in jobs.get(job['id'])['error']


def test_stop_kills_running_group_and_prevents_next_step(monkeypatch,tmp_path):
    out=tmp_path/'output';marker=tmp_path/'should_not_exist'
    steps=[dict(name='long',command=[sys.executable,'-c','import time; print("started",flush=True); time.sleep(60)'],output=str(out/'one'),kind='build'),
           dict(name='next',command=[sys.executable,'-c','from pathlib import Path; Path('+repr(str(marker))+').touch()'],output=str(out/'two'),kind='build')]
    stub(monkeypatch,out,steps)
    jobs=runner.Jobs(tmp_path);job=jobs.start({})
    wait_for(lambda:jobs.process is not None)
    process=jobs.process
    jobs.cancel(job['id'])
    wait_for(lambda:jobs.active is None)
    assert process.poll() is not None
    assert jobs.get(job['id'])['status']=='cancelled'
    assert not marker.exists()


def test_several_jobs_run_at_the_same_time(monkeypatch,tmp_path):
    def plan(config,project):
        out=Path(config['output'])
        return ({'output':str(out)},[dict(name='wait',kind='build',output=str(out),
                command=[sys.executable,'-c','import time; print("started",flush=True); time.sleep(1.2)'])])
    monkeypatch.setattr(runner,'make_plan',plan)
    monkeypatch.setattr(runner,'save_selected_result',lambda path:path)
    jobs=runner.Jobs(tmp_path)
    first=jobs.start({'output':str(tmp_path/'first')})
    second=jobs.start({'output':str(tmp_path/'second')})
    wait_for(lambda:len(jobs.processes)==2)
    assert len(jobs.live())==2
    assert all(process.poll() is None for process in jobs.processes.values())
    wait_for(lambda:not jobs.live())
    assert jobs.get(first['id'])['status']=='success'
    assert jobs.get(second['id'])['status']=='success'
    assert {row['id'] for row in jobs.history()}=={first['id'],second['id']}
    assert jobs.active is None and jobs.process is None


def test_parallel_limit_refuses_the_next_job(monkeypatch,tmp_path):
    stub(monkeypatch,tmp_path/'output',[dict(name='wait',kind='build',output=str(tmp_path/'output'),
         command=[sys.executable,'-c','import time; time.sleep(1.5)'])])
    jobs=runner.Jobs(tmp_path,max_parallel=1)
    job=jobs.start({})
    with pytest.raises(ValueError,match='同时运行的任务已达'):
        jobs.start({})
    jobs.cancel(job['id'])
    wait_for(lambda:not jobs.live())


def test_restart_marks_unfinished_job_interrupted(tmp_path):
    jobs=runner.Jobs(tmp_path)
    data=dict(id='old',status='running',output=str(tmp_path),owned=False,created=0)
    jobs.persist(data)
    restored=runner.Jobs(tmp_path)
    assert restored.items['old']['status']=='interrupted'


def test_start_refuses_to_reuse_an_existing_output(monkeypatch,tmp_path):
    out=tmp_path/'output';out.mkdir()
    (out/'keep').write_text('preserved')
    stub(monkeypatch,out,[dict(name='never',command=[sys.executable,'-c','pass'],output=str(out),kind='build')])
    jobs=runner.Jobs(tmp_path)
    with pytest.raises(FileExistsError):
        jobs.start({})
    assert (out/'keep').read_text()=='preserved'


def test_history_is_newest_first_and_import_is_read_only(monkeypatch,tmp_path):
    out=tmp_path/'output'
    stub(monkeypatch,out,[dict(name='ok',command=[sys.executable,'-c','pass'],output=str(out),kind='build')])
    jobs=runner.Jobs(tmp_path);job=jobs.start({})
    wait_for(lambda:jobs.active is None)
    assert [item['id'] for item in jobs.history()]==[job['id']]
    assert jobs.history()[0]['status']=='success'
