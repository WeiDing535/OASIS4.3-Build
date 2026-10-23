# SPDX-License-Identifier: LGPL-2.1-or-later
# Copyright (C) 2026 DingWei (丁玮) and contributors.
"""Native Qt desktop front end for independent SAD phasing; no external Parrot."""
from pathlib import Path
import json
import sys
import time


def main():
    try:
        from PySide6 import QtCore,QtGui,QtWidgets
    except ImportError:
        print('GUI dependency missing. Install from the project directory: python -m pip install ".[gui]"',file=sys.stderr)
        return 2
    app=QtWidgets.QApplication.instance() or QtWidgets.QApplication(sys.argv)
    window=make_window(QtCore,QtGui,QtWidgets)
    window.show()
    return app.exec()


def make_window(C,G,W):
    import gemmi
    from . import __product__, __version__
    from .desktop_job import build_sad_command
    from .reference_data import reference_info
    from .sequence import sequence_asu

    class Window(W.QMainWindow):
        def __init__(self):
            super().__init__()
            self.setWindowTitle(f'{__product__} · SAD 相位计算')
            self.resize(1000,850)
            self.process=C.QProcess(self)
            self.process.setProcessChannelMode(C.QProcess.ProcessChannelMode.SeparateChannels)
            self.process.readyReadStandardOutput.connect(self.read_stdout)
            self.process.readyReadStandardError.connect(self.read_log)
            self.process.finished.connect(self.finished)
            self.process.errorOccurred.connect(self.process_error)
            self.run_log=None;self.active_out=None;self.cancelled=False
            self.help_buttons={};self.sequence_estimate=None
            host=W.QWidget();self.setCentralWidget(host);layout=W.QVBoxLayout(host)
            title=W.QLabel('OASIS4.3-GUI-0.1  /  SAD 相位计算');title.setStyleSheet('font-size:24px;font-weight:600')
            layout.addWidget(title)
            note=W.QLabel('独立 Python 引擎 · 默认通用蛋白参考 · 无需 CCP4 / Parrot')
            layout.addWidget(note)
            scroll=W.QScrollArea();scroll.setWidgetResizable(True);layout.addWidget(scroll,3)
            self.form_host=W.QWidget();scroll.setWidget(self.form_host);form=W.QFormLayout(self.form_host)
            self.input=self.file_row(form,'实验 MTZ','open','MTZ (*.mtz)',inspect=True)
            self.sequence=self.file_row(form,'蛋白质序列文件','open','Sequence (*.fasta *.fa *.faa *.fas *.pir *.seq *.txt);;All files (*)',help_key='sequence',help_text='请选择实际结晶样品的蛋白序列文件。支持 FASTA、PIR 和无标题的单序列文本；包含标签、连接肽等实际构建序列。每个条目表示一条链，允许多个条目。不要上传 Word/PDF、结构文件或带缺口的比对。未知 X/B/Z/J 需要先明确。M 表示天然甲硫氨酸，硒代甲硫氨酸可写 (MSE)，硒代半胱氨酸写 (SEC) 或 U。')
            self.copies=W.QSpinBox();self.copies.setRange(1,10000);self.copies.setValue(1)
            self.field(form,'ASU 中的 copies',self.copies,'copies','整个序列文件在一个晶体学不对称单元（ASU）中重复的份数，必须为正整数，不是晶胞中的分子数。例如文件只有链 A，copies=2 表示 ASU 有两条 A；文件包含 A 和 B，copies=2 表示 A2B2。若文件已经列出完整 ASU 的每条链，填 1。不同链份数不同时，请在文件中按正确比例重复条目，再填写共同重复份数。')
            self.solvent=W.QLineEdit();self.solvent.setReadOnly(True);self.solvent.setPlaceholderText('选择 MTZ 和序列后自动估算')
            self.field(form,'溶剂含量（自动估算）',self.solvent,'solvent','使用序列分子量、copies、MTZ 晶胞和全部空间群对称操作数计算：Vm=Vcell/(Nsym×ASU分子量)，溶剂比例=1−1.23/Vm。显示百分比并自动传给密度修饰，不需手填。它是蛋白晶体的 Matthews 估算，不是测量值；糖基化、配体或未标记修饰未计入。异常结果请核对 copies 和结晶构建序列，软件不会自动猜测或改变 copies。')
            self.sequence_info=W.QLabel('等待序列和晶胞信息');self.sequence_info.setWordWrap(True);form.addRow(self.sequence_info)
            self.heavy=self.file_row(form,'重原子坐标','open','Coordinates (*.pdb *.cif *.mmcif)')
            self.output=self.file_row(form,'输出 MTZ','save','MTZ (*.mtz)')
            self.data_info=W.QLabel('选择实验 MTZ 后自动读取列名');self.data_info.setWordWrap(True);form.addRow(self.data_info)
            self.format=W.QComboBox();self.format.addItems(['F(+), F(-) 振幅对','F, ΔF 异常差振幅'])
            self.field(form,'数据格式',self.format,'format','选择输入的异常振幅形式：独立的 F(+)、F(−) 及各自标准差，或平均振幅 F、异常差 ΔF 及各自标准差。当前不直接接受 I 强度；I(+)/I(−) 必须先适当转换为振幅。')
            self.label_boxes={}
            label_grid=W.QWidget();grid=W.QGridLayout(label_grid);grid.setContentsMargins(0,0,0,0)
            for row,(key,label) in enumerate([('fplus','F(+)'),('sigfplus','σF(+)'),('fminus','F(-)'),('sigfminus','σF(-)'),('f','F'),('sigf','σF'),('dano','ΔF'),('sigdano','σΔF')]):
                box=W.QComboBox();box.setMinimumWidth(150);box.addItem('自动识别','')
                self.label_boxes[key]=box;grid.addWidget(W.QLabel(label),row//2,(row%2)*3);grid.addWidget(box,row//2,(row%2)*3+1)
                grid.addWidget(self.help_button('column_'+key,label,'选择实验 MTZ 中对应 '+label+' 的列。自动识别只匹配常见列名；自定义列名请手动选择。振幅列必须是振幅类型，标准差必须为匹配振幅的不确定度。灰色项目属于另一种数据格式，不参与计算。'),row//2,(row%2)*3+2)
            form.addRow('MTZ 列选择',label_grid);self.format.currentIndexChanged.connect(self.update_labels)
            self.fpp=W.QLineEdit();self.fpp.setPlaceholderText("例如 Se=3.84；必须与元素和实验波长对应，可逗号分隔")
            self.field(form,"异常散射 f''",self.fpp,"fpp","输入重原子在实验波长下的异常散射虚部，单位电子。例如 Se=3.84，多元素用逗号分隔。必须覆盖重原子坐标中的全部元素且为正数；示例值不能直接移用到其他实验。吸收边附近优先采用实验测量/可靠扫描值。")
            wave_row=W.QWidget();wave_layout=W.QHBoxLayout(wave_row);wave_layout.setContentsMargins(0,0,0,0)
            self.element=W.QLineEdit('Se');self.element.setMaximumWidth(70)
            self.wavelength=W.QLineEdit();self.wavelength.setPlaceholderText('波长 Å')
            calc=W.QPushButton("由波长估算 f''");calc.clicked.connect(self.calculate_fpp)
            wave_layout.addWidget(self.element);wave_layout.addWidget(self.help_button('element','元素','输入用于理论估算的化学元素符号，例如 Se、S、Hg；不会改变重原子文件的元素或自动推断标记方式。'))
            wave_layout.addWidget(self.wavelength);wave_layout.addWidget(self.help_button('wavelength','波长','输入采集 SAD 数据时的实际 X 射线波长，单位 Å，必须大于零。点击估算按钮会用 Cromer–Liberman 理论表计算该元素的 f 双撇；吸收边附近可能偏离实测值。'));wave_layout.addWidget(calc)
            form.addRow('可选理论估算',wave_row)
            self.method=W.QComboBox();self.method.addItem('连续 SAD + MLHL-best 密度修饰','likelihood');self.method.addItem('P+ + legacy 密度修饰','pplus')
            self.field(form,'算法',self.method,'method','连续 SAD 使用当前的连续相位似然方法，启用密度修饰时采用独立 MLHL-best 引擎；P+ 路径采用此前的 legacy 密度修饰。默认选择连续 SAD。自动组成和溶剂估算不会根据已知完整结构挑选算法。')
            self.dm=W.QCheckBox('启用独立密度修饰');self.dm.setChecked(True);self.field(form,'密度修饰',self.dm,'dm','在 SAD 初始相位后使用溶剂平坦化、直方图匹配和相位概率更新。溶剂比例来自上方自动估算。取消后仅输出 SAD 相位；全过程不调用外部 Parrot。')
            self.cycles=W.QSpinBox();self.cycles.setRange(1,100);self.cycles.setValue(5);self.field(form,'SAD 初始相位轮数',self.cycles,'cycles','初始 SAD/直接法相位计算轮数，默认 5，允许 1–100。更多轮次不保证相位更准确；它与密度修饰轮数分别设置。')
            self.dm_cycles=W.QSpinBox();self.dm_cycles.setRange(1,100);self.dm_cycles.setValue(5);self.field(form,'密度修饰轮数',self.dm_cycles,'dm_cycles','独立密度修饰迭代轮数，默认 5，允许 1–100。仅在启用密度修饰时使用。轮数增加不保证求解效果提高。')
            self.dmin=W.QLineEdit();self.dmin.setPlaceholderText('留空使用输入分辨率；复现此前基准填 2.8')
            self.field(form,'高分辨率截断 / Å',self.dmin,'dmin','单位 Å。留空保留输入数据的分辨率范围；输入 2.8 表示只使用 d≥2.8 Å 的反射。不能通过填写更小的数值创造高分辨率数据。溶剂估算取自晶胞与序列，不受截断值影响。')
            self.seed=W.QSpinBox();self.seed.setRange(0,2147483647);self.seed.setValue(20261016);self.field(form,'随机种子',self.seed,'seed','控制算法中的随机扰动，默认 20261016。相同输入、设置和依赖环境使用同一种子便于复现；改变种子不保证得到更好的结果。')
            self.reference=self.file_row(form,'可选自定义参考 MTZ','open','MTZ (*.mtz)')
            self.reference.setPlaceholderText('留空使用默认 1TQW 缓存参考，不需另装 CCP4')
            limitations=W.QLabel('当前 GUI 面向蛋白 SAD。FOM 仍是模型置信度，不能据此判断已经解出结构。输出 PHIB/FOM/HL 与运行元数据；SIR/DMR 可用命令行。')
            limitations.setWordWrap(True);form.addRow(limitations)
            buttons=W.QHBoxLayout();layout.addLayout(buttons)
            self.run_button=W.QPushButton('开始计算');self.run_button.clicked.connect(self.start)
            self.cancel_button=W.QPushButton('停止');self.cancel_button.setEnabled(False);self.cancel_button.clicked.connect(self.cancel)
            self.open_button=W.QPushButton('打开结果目录');self.open_button.setEnabled(False);self.open_button.clicked.connect(self.open_output)
            for button in [self.run_button,self.cancel_button,self.open_button]:buttons.addWidget(button)
            self.progress=W.QProgressBar();self.progress.setRange(0,1);self.progress.setValue(0);layout.addWidget(self.progress)
            self.status=W.QLabel('就绪');layout.addWidget(self.status)
            self.log=W.QPlainTextEdit();self.log.setReadOnly(True);self.log.setMaximumBlockCount(10000);layout.addWidget(self.log,2)
            self.update_labels()
            self.estimate_timer=C.QTimer(self);self.estimate_timer.setSingleShot(True);self.estimate_timer.setInterval(200)
            self.estimate_timer.timeout.connect(self.update_sequence_estimate)
            for entry in [self.input,self.sequence]:entry.textChanged.connect(lambda *_:self.estimate_timer.start())
            self.copies.valueChanged.connect(lambda *_:self.estimate_timer.start())

        def help_button(self,key,label,text):
            button=W.QToolButton();button.setText('?');button.setFixedWidth(25)
            button.setAccessibleName(label+' 帮助');button.setToolTip(text)
            button.clicked.connect(lambda checked=False:W.QMessageBox.information(self,label+' — 输入说明',text))
            self.help_buttons[key]=button
            return button

        def field(self,form,label,widget,key,text):
            row=W.QWidget();box=W.QHBoxLayout(row);box.setContentsMargins(0,0,0,0)
            box.addWidget(widget,1);box.addWidget(self.help_button(key,label,text));form.addRow(label,row)

        def update_sequence_estimate(self):
            self.sequence_estimate=None;self.solvent.clear();self.sequence_info.setStyleSheet('')
            if not self.sequence.text().strip() or not self.input.text().strip():
                self.sequence_info.setText('等待序列和晶胞信息');return
            try:
                m=gemmi.read_mtz_file(str(Path(self.input.text().strip()).expanduser()))
                info=sequence_asu(self.sequence.text().strip(),self.copies.value(),m.cell,m.spacegroup)
                self.sequence_estimate=info;self.solvent.setText(f"{info['solvent_fraction']:.2%}（Matthews 估算）")
                text=f"ASU：{info['asu_chains']} 条链 · {info['asu_residues']} 残基 · {info['asu_mass_da']/1000:.3f} kDa · Vm={info['matthews_vm']:.3f} Å³/Da"
                if info['warnings']:text+='\n'+'；'.join(info['warnings'])
                self.sequence_info.setText(text)
            except Exception as error:
                self.sequence_info.setText(str(error));self.sequence_info.setStyleSheet('color:#a33124')

        def file_row(self,form,label,mode,filter,inspect=False,help_key=None,help_text=None):
            row=W.QWidget();line=W.QLineEdit();button=W.QPushButton('选择…')
            box=W.QHBoxLayout(row);box.setContentsMargins(0,0,0,0);box.addWidget(line,1);box.addWidget(button)
            def choose():
                fn=W.QFileDialog.getSaveFileName if mode=='save' else W.QFileDialog.getOpenFileName
                path,_=fn(self,label,line.text(),filter)
                if path:
                    if mode=='save' and not path.lower().endswith('.mtz'):path+='.mtz'
                    line.setText(path)
                    if inspect:self.inspect_input()
            button.clicked.connect(choose)
            if inspect:line.editingFinished.connect(self.inspect_input)
            descriptions={
                '实验 MTZ':('input','选择实验 X 射线 SAD 振幅数据 MTZ。需要有效晶胞、空间群，以及异常振幅对或异常差和标准差。选择后自动读取列名；不会使用完整目标结构的相位。'),
                '重原子坐标':('heavy','选择重原子亚结构 PDB/mmCIF。元素、占有率、晶胞、空间群和坐标原点须与实验数据一致；不要用完整蛋白结构替代重原子文件。当前计算只接受单模型和各向同性 B 因子。'),
                '输出 MTZ':('output','选择新的 .mtz 输出文件名。同时保存同名 .json 参数记录和 .run.log 日志。不能覆盖已存在的结果，也不能覆盖输入。'),
                '可选自定义参考 MTZ':('reference','通常留空，使用首次安装时获取的 1TQW 通用参考。自定义时需要有效晶胞、空间群及 FC_BEST.F_phi.F（振幅）和 FC_BEST.F_phi.phi（相位）列。请选择独立于待求目标的合适蛋白地图，不是普通无相位实验 MTZ。')}
            key,text=descriptions.get(label,(help_key,help_text))
            self.field(form,label,row,key,text)
            return line

        def update_labels(self):
            keys={'fplus','sigfplus','fminus','sigfminus'} if self.format.currentIndex()==0 else {'f','sigf','dano','sigdano'}
            for key,box in self.label_boxes.items():box.setEnabled(key in keys)

        def inspect_input(self):
            if not self.input.text().strip():return
            try:
                m=gemmi.read_mtz_file(self.input.text().strip());cols=[c.label for c in m.columns if c.type!='H']
                defaults={'fplus':['F(+)'],'sigfplus':['SIGF(+)'],'fminus':['F(-)'],'sigfminus':['SIGF(-)'],
                    'f':['FP','F'],'sigf':['SIGFP','SIGF'],'dano':['DANO','DANOM'],'sigdano':['SIGDANO','SIGDANOM']}
                for key,box in self.label_boxes.items():
                    box.clear();box.addItem('自动识别','')
                    for name in cols:box.addItem(name,name)
                    for name in defaults[key]:
                        index=box.findData(name)
                        if index>=0:box.setCurrentIndex(index);break
                self.format.setCurrentIndex(1 if any(x in cols for x in ['DANO','DANOM']) else 0)
                self.data_info.setText(f'{m.spacegroup.hm} · {m.nreflections:,} 反射 · 最高分辨率 {m.resolution_high():.2f} Å')
                if not self.output.text():self.output.setText(str(Path(self.input.text()).with_name(Path(self.input.text()).stem+'_oasis.mtz')))
                self.update_sequence_estimate()
            except Exception as error:self.data_info.setText('MTZ 读取失败：'+str(error))

        def calculate_fpp(self):
            try:
                import math
                wave=float(self.wavelength.text());el=gemmi.Element(self.element.text().strip())
                if not math.isfinite(wave) or wave<=0 or not el.atomic_number:raise ValueError('请填写有效元素和正波长')
                value=gemmi.cromer_liberman(el.atomic_number,gemmi.hc/wave)[1]
                self.fpp.setText(f'{el.name}={value:.7g}')
                self.status.setText("已用 Cromer–Liberman 理论值估算 f''；吸收边附近请核对实验值")
            except Exception as error:W.QMessageBox.warning(self,"f'' 估算失败",str(error))

        def spec(self):
            active={'fplus','sigfplus','fminus','sigfminus'} if self.format.currentIndex()==0 else {'f','sigf','dano','sigdano'}
            return dict(input=self.input.text().strip(),heavy=self.heavy.text().strip(),output=self.output.text().strip(),
                labels={key:box.currentData() for key,box in self.label_boxes.items() if key in active and box.currentData()},
                fpp=self.fpp.text().strip(),sequence=self.sequence.text().strip(),copies=self.copies.value(),
                method=self.method.currentData(),density_modification=self.dm.isChecked(),solvent=self.solvent.text().strip(),
                cycles=self.cycles.value(),dm_cycles=self.dm_cycles.value(),dmin=self.dmin.text().strip(),
                seed=self.seed.value(),reference=self.reference.text().strip())

        def start(self):
            try:
                self.update_sequence_estimate()
                command=build_sad_command(self.spec());out=Path(self.output.text()).expanduser().resolve()
                out.parent.mkdir(parents=True,exist_ok=True)
                self.run_log=out.with_suffix('.run.log').open('x',encoding='utf-8')
                header=json.dumps(dict(command=command,sequence_asu=self.sequence_estimate,reference=reference_info() if not self.reference.text().strip() else self.reference.text().strip()),ensure_ascii=False,indent=2)
                self.log.clear()
                self.log.appendPlainText('开始 SAD 相位计算\n实验数据：'+self.input.text()+'\n序列：'+self.sequence.text()+'\nASU copies：'+str(self.copies.value())+'；估算溶剂含量：'+self.solvent.text()+'\n参考地图：'+('自定义参考' if self.reference.text().strip() else '默认 1TQW 通用蛋白参考')+'\n完整参数将保存在运行日志中。')
                self.run_log.write(header+'\n');self.run_log.flush()
                self.active_out=out;self.cancelled=False;self.started=time.monotonic()
                self.run_button.setEnabled(False);self.form_host.setEnabled(False);self.cancel_button.setEnabled(True);self.open_button.setEnabled(False)
                self.progress.setRange(0,0);self.status.setText('正在计算；完成时间随反射数和迭代轮数变化…')
                self.process.start(command[0],command[1:])
            except Exception as error:W.QMessageBox.warning(self,'输入检查失败',str(error))

        def read_stdout(self):
            text=bytes(self.process.readAllStandardOutput()).decode('utf-8',errors='replace')
            if text and self.run_log:self.run_log.write(text);self.run_log.flush()

        def read_log(self):
            text=bytes(self.process.readAllStandardError()).decode('utf-8',errors='replace')
            if text:
                self.log.moveCursor(G.QTextCursor.MoveOperation.End);self.log.insertPlainText(text)
                if self.run_log:self.run_log.write(text);self.run_log.flush()

        def finished(self,code,exit_status):
            self.read_stdout();self.read_log()
            if self.run_log:self.run_log.close();self.run_log=None
            self.run_button.setEnabled(True);self.form_host.setEnabled(True);self.cancel_button.setEnabled(False)
            self.progress.setRange(0,1)
            ok=code==0 and exit_status==C.QProcess.ExitStatus.NormalExit and self.active_out and self.active_out.is_file() and self.active_out.with_suffix('.json').is_file()
            self.progress.setValue(1 if ok else 0);self.open_button.setEnabled(bool(self.active_out))
            if ok:
                self.status.setText(f'计算完成 · {time.monotonic()-self.started:.1f} 秒 · {self.active_out}')
                self.log.appendPlainText('计算完成，已保存：\n'+str(self.active_out)+'\n'+str(self.active_out.with_suffix('.json'))+'\n'+str(self.active_out.with_suffix('.run.log')))
            elif self.cancelled:self.status.setText('计算已停止；保留运行日志，输出可能不完整，请使用新名称重跑。')
            else:self.status.setText(f'计算失败（退出码 {code}）；请查看下方日志。')

        def process_error(self,error):
            self.log.appendPlainText(self.process.errorString())
            if error==C.QProcess.ProcessError.FailedToStart:self.finished(-1,C.QProcess.ExitStatus.CrashExit)

        def cancel(self):
            self.cancelled=True;self.cancel_button.setEnabled(False);self.process.terminate()
            C.QTimer.singleShot(3000,self.kill_if_running)

        def kill_if_running(self):
            if self.cancelled and self.process.state()!=C.QProcess.ProcessState.NotRunning:self.process.kill()

        def open_output(self):
            if self.active_out:G.QDesktopServices.openUrl(C.QUrl.fromLocalFile(str(self.active_out.parent)))

        def closeEvent(self,event):
            if self.process.state()!=C.QProcess.ProcessState.NotRunning:
                answer=W.QMessageBox.question(self,'退出','计算仍在运行，停止并退出？')
                if answer!=W.QMessageBox.StandardButton.Yes:event.ignore();return
                self.process.kill();self.process.waitForFinished(3000)
            event.accept()
    return Window()

if __name__=='__main__':raise SystemExit(main())
