"""Real-window check: task name, relaunch ("沿用参数") and delete ("删除记录")."""
import json,os,sys,tempfile
from pathlib import Path
sys.path.insert(0,os.environ['SRC'])
from phasebuild import desktop

root=Path(os.environ['SRC']).parents[1]
source=json.loads(Path(os.environ['LIVE_JOB']).read_text())
config=dict(source['config'])
output=Path(config['output'])
assert Path(config['sad']).exists() and output.is_dir()

tmp=Path(tempfile.mkdtemp(prefix='oasis-tasktools-'))
(tmp/'run_sad_to_model.py').symlink_to(root/'run_sad_to_model.py')
storage=tmp/'.gui-oasis/jobs';storage.mkdir(parents=True)
record=dict(source,id='feed1234',output=str(output),name='7ENM 参数复用',owned=True)
(storage/'feed1234.json').write_text(json.dumps(record,ensure_ascii=False))

controller=desktop.Controller(tmp)
app=desktop.App(controller)
app.root.update_idletasks();app.root.update()
print('name entry exists:',bool(app.task_name.get()==''))
print('history columns:',app.history['columns'])
app._refresh_history(force=True)
print('row name/status:',app.history.set('feed1234','name'),'/',app.history.set('feed1234','status'))
app.tabs.select(0);app.root.update_idletasks()
print('toolbar buttons on tab 0:',app.relaunch_button.winfo_ismapped(),app.delete_button.winfo_ismapped(),
      app.relaunch_button.winfo_width(),app.delete_button.winfo_width())
app.history.selection_set('feed1234');app.on_select_job()
app.on_relaunch()
print('relaunch sad kept:',app.fields['sad']['var'].get()==config['sad'])
print('relaunch wavelength kept:',app.fields['wavelength']['var'].get()==str(config['wavelength']))
print('relaunch output is new:',app.fields['output']['var'].get()!=str(output))
print('relaunch name kept:',app.task_name.get())
print('advanced auto-shown:',app.advanced.get())
print('page after relaunch:',app.tabs.index(app.tabs.select()),'(4=运行计划)')
print('summary:',app.summary.cget('text'))
print('buttons mapped:',app.relaunch_button.winfo_ismapped(),app.delete_button.winfo_ismapped())
print('delete keeps results:')
controller.delete('feed1234')
print('  record file gone:',not (storage/'feed1234.json').exists())
print('  output dir kept:',output.is_dir())
print('  history empty:',app.controller.history()==[])
app._refresh_history(force=True)
print('  window rows:',app.history.get_children())
app.root.destroy()
