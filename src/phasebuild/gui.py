"""Loopback-only local GUI and cancellable, persistent subprocess jobs."""
import argparse
import hmac
import json
import mimetypes
import os
import secrets
import signal
import subprocess
import sys
import threading
import time
import traceback
import webbrowser
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from urllib.parse import parse_qs, urlparse
from .gui_workflow import defaults, make_plan, collect_results, save_selected_result
from .gui_paths import browse_path

ASSETS = Path(__file__).with_name('gui_assets')


def save_json(path, value):
    path = Path(path)
    temp = path.with_suffix('.tmp')
    temp.write_text(json.dumps(value,indent=2,ensure_ascii=False,allow_nan=False)+'\n')
    temp.replace(path)


class Jobs:
    def __init__(self, project):
        self.project = Path(project).resolve()
        self.storage = self.project/'.gui/jobs'
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
                    item.update(status='interrupted',error='上次 GUI 服务已退出，请检查日志并以新目录重新运行。')
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
            cfg, steps = make_plan(config,self.project)
            Path(cfg['output']).mkdir(parents=True,exist_ok=False)
            now = time.time()
            job = dict(id=secrets.token_hex(8),owned=True,created=now,updated=now,
                       status='queued',output=cfg['output'],config=cfg,steps=steps,
                       completed=0,current='',error=None)
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
                    save_selected_result(job['output'])
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
        return job

    def history(self):
        with self.lock:
            return [dict(id=j['id'],created=j['created'],status=j['status'],output=j['output'])
                    for j in sorted(self.items.values(),key=lambda j:j['created'],reverse=True)]

    def import_results(self, directory):
        results = collect_results(directory)
        if not results['artifacts']:
            raise ValueError('该目录没有找到 PhaseBuild 结果文件')
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


class Server(ThreadingHTTPServer):
    daemon_threads = True


class Handler(BaseHTTPRequestHandler):
    def log_message(self, message, *args):
        pass

    def respond(self, value, code=200):
        body = json.dumps(value,ensure_ascii=False,allow_nan=False).encode('utf-8')
        self.send_response(code)
        self.headers_common('application/json; charset=utf-8',len(body))
        self.wfile.write(body)

    def headers_common(self, content_type, size):
        self.send_header('Content-Type',content_type)
        self.send_header('Content-Length',str(size))
        self.send_header('Cache-Control','no-store')
        self.send_header('X-Content-Type-Options','nosniff')
        self.send_header('Referrer-Policy','no-referrer')
        self.send_header('Content-Security-Policy',"default-src 'self'; script-src 'self'; style-src 'self'; img-src 'self' data:; connect-src 'self'; frame-ancestors 'none'; base-uri 'none'; form-action 'self'")
        self.end_headers()

    def authorized(self):
        if self.headers.get('Host') != self.server.address:
            self.respond({'error':'无效的本机地址'},403)
            return False
        origin = self.headers.get('Origin')
        if origin and origin != 'http://'+self.server.address:
            self.respond({'error':'拒绝跨站请求'},403)
            return False
        supplied = self.headers.get('X-PhaseBuild-Token','')
        # Downloads opened in another tab cannot attach the custom header.
        if self.command=='GET' and urlparse(self.path).path=='/api/download':
            supplied = parse_qs(urlparse(self.path).query).get('token',[''])[0]
        if not hmac.compare_digest(supplied,self.server.token):
            self.respond({'error':'本机访问凭据失效，请使用启动时的链接重新打开'},403)
            return False
        return True

    def do_GET(self):
        parsed = urlparse(self.path)
        route = parsed.path
        if route in ('/','/app.js','/style.css'):
            if self.headers.get('Host') != self.server.address:
                self.respond({'error':'无效地址'},403)
                return
            file = ASSETS/('index.html' if route=='/' else route[1:])
            body = file.read_bytes()
            self.send_response(200)
            self.headers_common(mimetypes.guess_type(str(file))[0]+'; charset=utf-8',len(body))
            self.wfile.write(body)
            return
        if not self.authorized():
            return
        query = parse_qs(parsed.query)
        value = lambda key,default='':query.get(key,[default])[0]
        try:
            if route=='/api/defaults':
                self.respond(dict(config=defaults(self.server.jobs.project),history=self.server.jobs.history(),
                                  python=sys.executable,version='0.8.3',project=str(self.server.jobs.project)))
            elif route=='/api/browse':
                path,selected,fallback = browse_path(value('path',str(self.server.jobs.project)),value('allow_parent')=='1')
                entries = []
                for item in sorted(path.iterdir(),key=lambda p:(not p.is_dir(),p.name.lower())):
                    if item.name.startswith('.'):
                        continue
                    if item.is_dir() or item.suffix.lower() in ('.mtz','.cif','.pdb','.pir','.fasta','.fa','.seq','.json','.npz'):
                        entries.append(dict(name=item.name,path=str(item),directory=item.is_dir()))
                    if len(entries)>=2000:
                        break
                self.respond(dict(path=str(path),parent=str(path.parent),entries=entries,
                                  selected_file=selected,parent_fallback=fallback))
            elif route=='/api/job':
                self.respond(self.server.jobs.get(value('id')))
            elif route=='/api/results':
                job = self.server.jobs.get(value('id'))
                self.respond(collect_results(job['output']))
            elif route=='/api/download':
                job = self.server.jobs.get(value('id'))
                root = Path(job['output']).resolve()
                file = (root/value('file')).resolve()
                if root not in file.parents or not file.is_file() or file.suffix.lower() not in ('.cif','.pdb','.mtz','.json','.log','.txt','.png','.svg','.ccp4','.csv'):
                    raise ValueError('无效结果文件')
                self.send_response(200)
                self.send_header('Content-Disposition',"attachment; filename*=UTF-8''"+__import__('urllib.parse',fromlist=['quote']).quote(file.name))
                self.headers_common('application/octet-stream',file.stat().st_size)
                with file.open('rb') as handle:
                    while True:
                        chunk = handle.read(1024*1024)
                        if not chunk:
                            break
                        self.wfile.write(chunk)
            else:
                self.respond({'error':'接口不存在'},404)
        except (KeyError,ValueError,OSError,RuntimeError) as error:
            self.respond({'error':str(error)},400)

    def do_POST(self):
        if not self.authorized():
            return
        try:
            size = int(self.headers.get('Content-Length','0'))
            if not 0<size<=1_000_000:
                raise ValueError('请求过大或为空')
            data = json.loads(self.rfile.read(size))
            route = urlparse(self.path).path
            if route=='/api/plan':
                cfg,steps = make_plan(data,self.server.jobs.project)
                self.respond(dict(config=cfg,steps=steps))
            elif route=='/api/start':
                self.respond(self.server.jobs.start(data))
            elif route=='/api/stop':
                self.server.jobs.cancel(data['id'])
                self.respond({'ok':True})
            elif route=='/api/import':
                self.respond(self.server.jobs.import_results(data['path']))
            elif route=='/api/audit':
                from .data import audit
                self.respond(audit(data['path']))
            else:
                self.respond({'error':'接口不存在'},404)
        except Exception as error:
            self.respond({'error':str(error)},400)


def create_server(project, port=0, token=None):
    server = Server(('127.0.0.1',port),Handler)
    server.address = '127.0.0.1:'+str(server.server_address[1])
    server.token = token or secrets.token_urlsafe(32)
    server.jobs = Jobs(project)
    return server


def main():
    p = argparse.ArgumentParser(description='PhaseBuild 本机图形界面')
    p.add_argument('--port',type=int,default=0,help='默认自动选择空闲端口')
    p.add_argument('--no-browser',action='store_true')
    p.add_argument('--project',default=str(Path(__file__).resolve().parents[2]))
    args = p.parse_args()
    project = Path(args.project).resolve()
    if not (project/'run_phasebuild.py').is_file():
        p.error('请选择包含 run_phasebuild.py 的 PhaseBuild 源码目录（--project）')
    # Fail early with a useful runtime error, before presenting an unusable UI.
    try:
        import gemmi, numpy, scipy
    except ImportError as error:
        p.error('缺少科学计算依赖，请使用完整便携包中的 PhaseBuild.command 启动：'+str(error))
    server = create_server(project,args.port)
    url = 'http://'+server.address+'/#token='+server.token
    print('PhaseBuild GUI：'+url,flush=True)
    print('仅在本机运行。按 Ctrl+C 退出服务；关闭网页不会停止计算。',flush=True)
    if not args.no_browser:
        threading.Timer(.4,lambda:webbrowser.open(url)).start()
    def request_exit(signum,frame):
        raise KeyboardInterrupt
    for name in ('SIGTERM','SIGHUP'):
        if hasattr(signal,name):
            signal.signal(getattr(signal,name),request_exit)
    try:
        server.serve_forever(poll_interval=.3)
    except KeyboardInterrupt:
        pass
    finally:
        with server.jobs.lock:
            active,process,worker = server.jobs.active,server.jobs.process,server.jobs.worker
        if active:
            server.jobs.cancel(active)
        if worker:
            # Let the CLI forward SIGTERM to its separately grouped engine
            # before resorting to SIGKILL of the CLI process group.
            worker.join(timeout=2.5)
        if process and process.poll() is None:
            # Server exit must not leave descendants behind when daemon timers end.
            try:
                os.killpg(process.pid,signal.SIGKILL)
            except ProcessLookupError:
                pass
        if worker:
            worker.join(timeout=5)
        server.server_close()


if __name__=='__main__':
    main()
