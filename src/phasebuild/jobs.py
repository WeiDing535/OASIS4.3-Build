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
    def __init__(self, project, workflow=None):
        self.workflow=workflow
        self.project = Path(project).resolve()
        self.storage = self.project/('.gui-oasis/jobs' if workflow else '.gui/jobs')
        self.storage.mkdir(parents=True,exist_ok=True)
        self.lock = threading.RLock()
        self.items = {}
        self.process = None
        self.worker = None
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

    def persist(self, job):
        save_json(self.storage/(job['id']+'.json'),job)
        if job.get('owned'):
            save_json(Path(job['output'])/'gui_job.json',job)

    def start(self, config):
        with self.lock:
            if self.active:
                raise ValueError('已有任务运行中，请等待完成或停止它')
            cfg, steps = (self.workflow.make_plan if self.workflow else make_plan)(config,self.project)
            Path(cfg['output']).mkdir(parents=True,exist_ok=False)
            now = time.time()
            job = dict(id=secrets.token_hex(8),owned=True,created=now,updated=now,
                       status='queued',output=cfg['output'],config=cfg,steps=steps,
                       completed=0,current='',error=None)
            if self.workflow:job['plan']=self.workflow.describe_plan(cfg)
            self.items[job['id']] = job
            self.active = job['id']
            self.persist(job)
            self.worker = threading.Thread(target=self.run,args=(job,),daemon=True)
            self.worker.start()
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
                        self.process = subprocess.Popen(step['command'],cwd=str(self.project),env=env,
                                                        stdout=log,stderr=subprocess.STDOUT,
                                                        start_new_session=True)
                        process = self.process
                    code = process.wait()
                    with self.lock:
                        self.process = None
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
                self.active = None
                self.process = None

    def cancel(self, identifier):
        with self.lock:
            job = self.items[identifier]
            if self.active != identifier:
                return
            job['status'] = 'stopping'
            self.persist(job)
            process = self.process
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
        with self.lock:
            return [dict(id=j['id'],created=j['created'],status=j['status'],output=j['output'])
                    for j in sorted(self.items.values(),key=lambda j:j['created'],reverse=True)]

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
