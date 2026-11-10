"""Cancellable, persistent subprocess jobs shared by every front end.

The runner is deliberately toolkit-free: it owns the worker thread, the child
process group and the on-disk job records, and knows nothing about HTTP or
windows.  The desktop application drives it directly; the archived workbench
used to wrap it in a loopback HTTP server.
"""
import json
import os
import secrets
import signal
import subprocess
import threading
import time
import traceback
from pathlib import Path


def save_json(path, value):
    path = Path(path)
    temp = path.with_suffix('.tmp')
    temp.write_text(json.dumps(value,indent=2,ensure_ascii=False,allow_nan=False)+'\n')
    temp.replace(path)


# Legacy entry points are loaded only by the archived workbench. The focused
# application supplies its own workflow and does not import legacy algorithms.
def defaults(*args,**kwargs):
    from .gui_workflow import defaults as call
    return call(*args,**kwargs)
def make_plan(*args,**kwargs):
    from .gui_workflow import make_plan as call
    return call(*args,**kwargs)
def collect_results(*args,**kwargs):
    from .gui_workflow import collect_results as call
    return call(*args,**kwargs)
def save_selected_result(*args,**kwargs):
    from .gui_workflow import save_selected_result as call
    return call(*args,**kwargs)


class Jobs:
    """Owns every job of one project: threads, child groups, on-disk records.

    Several jobs may run at the same time; each one owns its own worker thread,
    child process group and ``gui.log``, so the front ends can start a new run
    while an earlier one is still refining.  ``max_parallel`` keeps a single
    machine from being oversubscribed by accident.
    """

    def __init__(self, project, workflow=None, max_parallel=4):
        self.workflow=workflow
        self.project = Path(project).resolve()
        self.storage = self.project/('.gui-oasis/jobs' if workflow else '.gui/jobs')
        self.storage.mkdir(parents=True,exist_ok=True)
        self.lock = threading.RLock()
        self.items = {}
        self.processes = {}
        self.workers = {}
        self.running = set()
        self.max_parallel = max_parallel
        self.active = None
        for index in self.storage.glob('*.json'):
            try:
                item = json.loads(index.read_text())
                if item['status'] in ('running','stopping','queued'):
                    item.update(status='interrupted',error='上次界面已退出，请检查日志并以新目录重新运行。')
                    self.persist(item)
                self.items[item['id']] = item
            except (OSError,ValueError,KeyError):
                continue

    @property
    def process(self):
        """Child process of the most recently started job (compatibility)."""
        return self.processes.get(self.active)

    @property
    def worker(self):
        """Worker thread of the latest job, running or just finished."""
        with self.lock:
            if self.active:
                return self.workers.get(self.active)
            return next(reversed(self.workers.values()),None) if self.workers else None

    def live(self):
        """Identifiers of the jobs that are queued, running or stopping."""
        with self.lock:
            return [identifier for identifier in self.running
                    if self.items.get(identifier,{}).get('status') in ('queued','running','stopping')]

    def persist(self, job):
        save_json(self.storage/(job['id']+'.json'),job)
        if job.get('owned'):
            save_json(Path(job['output'])/'gui_job.json',job)

    def start(self, config):
        with self.lock:
            if len(self.live())>=self.max_parallel:
                raise ValueError('同时运行的任务已达 %d 个，请等待其中一个完成或停止它。'%self.max_parallel)
            for identifier in [key for key,thread in self.workers.items()
                               if key not in self.running and not thread.is_alive()]:
                del self.workers[identifier]  # keep the last finished one for inspection
            cfg, steps = (self.workflow.make_plan if self.workflow else make_plan)(config,self.project)
            Path(cfg['output']).mkdir(parents=True,exist_ok=False)
            now = time.time()
            job = dict(id=secrets.token_hex(8),owned=True,created=now,updated=now,
                       status='queued',output=cfg['output'],config=cfg,steps=steps,
                       completed=0,current='',error=None,
                       label=Path(str(cfg.get('sad') or '')).name,method=cfg.get('method',''))
            if self.workflow:job['plan']=self.workflow.describe_plan(cfg)
            self.items[job['id']] = job
            self.running.add(job['id'])
            self.active = job['id']
            self.persist(job)
            worker = threading.Thread(target=self.run,args=(job,),daemon=True)
            self.workers[job['id']] = worker
            worker.start()
            return dict(job)

    def run(self, job):
        log_path = Path(job['output'])/'gui.log'
        try:
            with log_path.open('w',buffering=1) as log:
                for step in job['steps']:
                    with self.lock:
                        if job['status']=='stopping':
                            break
                        job.update(status='running',current=step['name'],updated=time.time())
                        self.persist(job)
                        log.write('\n▶ '+step['name']+'\n')
                        log.flush()
                        env = dict(os.environ,PYTHONUNBUFFERED='1')
                        process = subprocess.Popen(step['command'],cwd=str(self.project),env=env,
                                                   stdout=log,stderr=subprocess.STDOUT,
                                                   start_new_session=True)
                        self.processes[job['id']] = process
                    code = process.wait()
                    with self.lock:
                        self.processes.pop(job['id'],None)
                        if job['status']=='stopping':
                            break
                        if code:
                            raise RuntimeError('{} 失败（退出码 {}），请查看日志；已有输出已保留。'.format(step['name'],code))
                        if step['kind']=='trace':
                            report = json.loads((Path(step['output'])/'build_report.json').read_text())
                            if report.get('residues',0)<3:
                                raise RuntimeError('未找到满足筛选条件的主链片段，已停止后续步骤。请检查相位图，或增加旋转采样/种子优化数后使用新目录重试。')
                        job['completed'] += 1
                        job['updated'] = time.time()
                        self.persist(job)
                if job['status']!='stopping':
                    (self.workflow.save_selected_result if self.workflow else save_selected_result)(job['output'])
                with self.lock:
                    job.update(status='cancelled' if job['status']=='stopping' else 'success',current='')
        except Exception as error:
            with self.lock:
                job.update(status='cancelled' if job['status']=='stopping' else 'failed',error=str(error))
            with log_path.open('a') as log:
                log.write('\n'+traceback.format_exc())
        finally:
            with self.lock:
                job['updated'] = time.time()
                self.persist(job)
                self.processes.pop(job['id'],None)
                self.running.discard(job['id'])
                if self.active == job['id']:
                    # Keep ``active`` pointing at another live job, if any.
                    remaining = self.live()
                    self.active = remaining[-1] if remaining else None

    def cancel(self, identifier):
        with self.lock:
            job = self.items.get(identifier)
            if job is None or job.get('status') not in ('queued','running','stopping'):
                return
            if job['status'] == 'stopping':
                return
            job['status'] = 'stopping'
            self.persist(job)
            process = self.processes.get(identifier)
            if process:
                try:
                    os.killpg(process.pid,signal.SIGTERM)
                except ProcessLookupError:
                    pass
                def force_stop():
                    # Group leaders may exit before their REFMAC/build descendants.
                    try:
                        os.killpg(process.pid,signal.SIGKILL)
                    except ProcessLookupError:
                        pass
                timer = threading.Timer(3,force_stop)
                timer.daemon = True
                timer.start()

    def get(self, identifier):
        with self.lock:
            job = dict(self.items[identifier])
        path = Path(job['output'])/'gui.log'
        if path.is_file():
            with path.open('rb') as log:
                log.seek(max(0,path.stat().st_size-60000))
                job['log'] = log.read().decode('utf-8',errors='replace')
        else:
            job['log'] = ''
        if self.workflow and job.get('config') and not job.get('plan'):
            job['plan']=self.workflow.describe_plan(job['config'])
        if job.get('config',{}).get('mode')=='sad_to_model':
            state=Path(job['output'])/'00_sad_to_model/pipeline_report.json'
            if state.is_file():
                saved=json.loads(state.read_text())
                job['pipeline_progress']={'completed':len(saved.get('completed',[])), 'total':5,
                                          'current':saved.get('current_name','')}
        return self.workflow.enrich_job(job) if self.workflow else job

    def history(self):
        """Every known job, newest first, with the fields a task list shows."""
        with self.lock:
            rows = []
            for job in sorted(self.items.values(),key=lambda j:j['created'],reverse=True):
                config = job.get('config') or {}
                rows.append(dict(id=job['id'],created=job['created'],
                                 updated=job.get('updated',job['created']),
                                 status=job['status'],output=job['output'],
                                 completed=job.get('completed',0),
                                 total=len(job.get('steps') or []),
                                 error=job.get('error'),
                                 method=job.get('method') or config.get('method',''),
                                 label=job.get('label') or Path(str(config.get('sad') or '')).name))
            return rows

    def import_results(self, directory):
        if self.workflow:self.workflow.validate_import(directory)
        results = (self.workflow.collect_results if self.workflow else collect_results)(directory)
        if not results['artifacts']:
            raise ValueError('该目录没有找到 OASIS4.3 结果文件')
        with self.lock:
            for job in self.items.values():
                if job['output']==results['root']:
                    return dict(job)
            job = dict(id=secrets.token_hex(8),created=time.time(),updated=time.time(),
                       owned=False,status='imported',output=results['root'],completed=0,
                       steps=[],current='',error=None)
            self.items[job['id']] = job
            self.persist(job)
            return dict(job)
