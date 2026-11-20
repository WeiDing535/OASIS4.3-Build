# SPDX-License-Identifier: LGPL-2.1-or-later
# Copyright (C) 2026 DingWei and contributors.
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
            self.setWindowTitle(f'{__product__} · SAD phasing')
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
            title=W.QLabel('OASIS4.3-GUI-0.1  /  SAD phasing');title.setStyleSheet('font-size:24px;font-weight:600')
            layout.addWidget(title)
            note=W.QLabel('standalone Python engine · default general protein reference · no CCP4 / Parrot needed')
            layout.addWidget(note)
            scroll=W.QScrollArea();scroll.setWidgetResizable(True);layout.addWidget(scroll,3)
            self.form_host=W.QWidget();scroll.setWidget(self.form_host);form=W.QFormLayout(self.form_host)
            self.input=self.file_row(form,'experimental MTZ','open','MTZ (*.mtz)',inspect=True)
            self.sequence=self.file_row(form,'protein sequence file','open','Sequence (*.fasta *.fa *.faa *.fas *.pir *.seq *.txt);;All files (*)',help_key='sequence',help_text='please provide the protein sequence file of the actual crystallized sample. FASTA, PIR, and untitled single-sequence text are supported; include tags, linkers, and other sequences actually built. Each entry represents one chain; multiple entries are allowed. Do not upload Word/PDF files, structure files, or gapped alignments. Unknown residues X/B/Z/J must be resolved first. M is natural methionine; selenomethionine may be written as (MSE), and selenocysteine as (SEC) or U.')
            self.copies=W.QSpinBox();self.copies.setRange(1,10000);self.copies.setValue(1)
            self.field(form,'copies in ASU',self.copies,'copies','number of times the whole sequence file repeats in one crystallographic asymmetric unit (ASU); must be a positive integer, not the number of molecules in the cell. For example, if the file contains only chain A, copies=2 means two A chains in the ASU; if it contains A and B, copies=2 means A2B2. If the file already lists every chain of the full ASU, enter 1. When chains have different copy numbers, repeat entries in the file at the correct ratio and then enter the common repetition factor.')
            self.solvent=W.QLineEdit();self.solvent.setReadOnly(True);self.solvent.setPlaceholderText('estimated automatically after selecting the MTZ and sequence')
            self.field(form,'solvent content (estimated automatically)',self.solvent,'solvent','computed from the sequence mass, copies, MTZ cell, and all space-group symmetry operations: Vm=Vcell/(Nsym×ASU mass), solvent fraction=1−1.23/Vm. The percentage is displayed and passed to density modification automatically; no manual entry is needed. It is a Matthews estimate for protein crystals, not a measured value; glycosylation, ligands, and unlabeled modifications are not counted. For unusual results, check copies and the crystallized construct sequence; the software never guesses or changes copies.')
            self.sequence_info=W.QLabel('waiting for sequence and cell information');self.sequence_info.setWordWrap(True);form.addRow(self.sequence_info)
            self.heavy=self.file_row(form,'heavy-atom coordinates','open','Coordinates (*.pdb *.cif *.mmcif)')
            self.output=self.file_row(form,'output MTZ','save','MTZ (*.mtz)')
            self.data_info=W.QLabel('column names are read automatically after selecting the experimental MTZ');self.data_info.setWordWrap(True);form.addRow(self.data_info)
            self.format=W.QComboBox();self.format.addItems(['F(+), F(-) amplitude pair','F, ΔF anomalous difference amplitudes'])
            self.field(form,'data format',self.format,'format','select the form of the input anomalous amplitudes: separate F(+) and F(−) with their standard deviations, or mean amplitude F, anomalous difference ΔF, and their standard deviations. Intensities I are not accepted directly; I(+)/I(−) must first be converted to amplitudes appropriately.')
            self.label_boxes={}
            label_grid=W.QWidget();grid=W.QGridLayout(label_grid);grid.setContentsMargins(0,0,0,0)
            for row,(key,label) in enumerate([('fplus','F(+)'),('sigfplus','σF(+)'),('fminus','F(-)'),('sigfminus','σF(-)'),('f','F'),('sigf','σF'),('dano','ΔF'),('sigdano','σΔF')]):
                box=W.QComboBox();box.setMinimumWidth(150);box.addItem('automatic','')
                self.label_boxes[key]=box;grid.addWidget(W.QLabel(label),row//2,(row%2)*3);grid.addWidget(box,row//2,(row%2)*3+1)
                grid.addWidget(self.help_button('column_'+key,label,'select the column in the experimental MTZ corresponding to '+label+'. Automatic detection matches only common column names; choose custom names manually. Amplitude columns must be of amplitude type, and standard deviations must be uncertainties of the matching amplitudes. Greyed-out items belong to the other data format and are not used.'),row//2,(row%2)*3+2)
            form.addRow('MTZ column selection',label_grid);self.format.currentIndexChanged.connect(self.update_labels)
            self.fpp=W.QLineEdit();self.fpp.setPlaceholderText("e.g. Se=3.84; must match the elements and experimental wavelength; comma-separated")
            self.field(form,"anomalous scattering f''",self.fpp,"fpp","enter the imaginary anomalous-scattering term at the experimental wavelength, in electrons. For example, Se=3.84; separate multiple elements with commas. It must cover every element in the heavy-atom coordinates and be positive; example values cannot be transferred to other experiments. Near an absorption edge, prefer measured values or a reliable scan.")
            wave_row=W.QWidget();wave_layout=W.QHBoxLayout(wave_row);wave_layout.setContentsMargins(0,0,0,0)
            self.element=W.QLineEdit('Se');self.element.setMaximumWidth(70)
            self.wavelength=W.QLineEdit();self.wavelength.setPlaceholderText('wavelength (Å)')
            calc=W.QPushButton("estimate f'' from wavelength");calc.clicked.connect(self.calculate_fpp)
            wave_layout.addWidget(self.element);wave_layout.addWidget(self.help_button('element','element','enter the chemical element symbol for theoretical estimation, e.g. Se, S, Hg; this does not change the elements in the heavy-atom file or infer labeling.'))
            wave_layout.addWidget(self.wavelength);wave_layout.addWidget(self.help_button('wavelength','wavelength','enter the actual X-ray wavelength used to collect the SAD data, in Å; it must be positive. Clicking estimate computes the f-double-prime term for the element from the Cromer–Liberman tables; near an absorption edge it may deviate from the measured value.'));wave_layout.addWidget(calc)
            form.addRow('optional theoretical estimate',wave_row)
            self.method=W.QComboBox();self.method.addItem('continuous SAD + MLHL-best density modification','likelihood');self.method.addItem('P+ + legacy density modification','pplus')
            self.field(form,'algorithm',self.method,'method','continuous SAD uses the current continuous phasing likelihood method, with the standalone MLHL-best engine for density modification; the P+ path uses the previous legacy density modification. Continuous SAD is selected by default. Automatic composition and solvent estimation do not choose an algorithm based on a known complete structure.')
            self.dm=W.QCheckBox('enable standalone density modification');self.dm.setChecked(True);self.field(form,'density modification',self.dm,'dm','apply solvent flattening, histogram matching, and phase-probability updating after the initial SAD phases. The solvent fraction comes from the automatic estimate above. If unchecked, only SAD phases are output; external Parrot is never called.')
            self.cycles=W.QSpinBox();self.cycles.setRange(1,100);self.cycles.setValue(5);self.field(form,'SAD initial phasing cycles',self.cycles,'cycles','number of initial SAD/direct-method phasing cycles, default 5, allowed 1–100. More cycles do not guarantee more accurate phases; this is set separately from the density-modification cycles.')
            self.dm_cycles=W.QSpinBox();self.dm_cycles.setRange(1,100);self.dm_cycles.setValue(5);self.field(form,'density-modification cycles',self.dm_cycles,'dm_cycles','number of standalone density-modification iterations, default 5, allowed 1–100. Used only when density modification is enabled. More cycles do not guarantee a better solution.')
            self.dmin=W.QLineEdit();self.dmin.setPlaceholderText('leave blank to use the input resolution; use 2.8 to reproduce earlier benchmarks')
            self.field(form,'high-resolution cutoff / Å',self.dmin,'dmin','in Å. Leave blank to keep the input data resolution range; entering 2.8 uses only reflections with d≥2.8 Å. Higher-resolution data cannot be created by entering a smaller number. The solvent estimate is derived from the cell and sequence and is unaffected by the cutoff.')
            self.seed=W.QSpinBox();self.seed.setRange(0,2147483647);self.seed.setValue(20261016);self.field(form,'random seed',self.seed,'seed','controls random perturbations in the algorithm, default 20261016. The same seed with the same inputs, settings, and dependencies aids reproducibility; changing it does not guarantee better results.')
            self.reference=self.file_row(form,'optional custom reference MTZ','open','MTZ (*.mtz)')
            self.reference.setPlaceholderText('leave blank to use the bundled default 1TQW reference; no extra CCP4 installation needed')
            limitations=W.QLabel('this GUI targets protein SAD. FOM is a model confidence measure and cannot be used to judge whether the structure is solved. PHIB/FOM/HL and run metadata are output; SIR/DMR are available via the command line.')
            limitations.setWordWrap(True);form.addRow(limitations)
            buttons=W.QHBoxLayout();layout.addLayout(buttons)
            self.run_button=W.QPushButton('start calculation');self.run_button.clicked.connect(self.start)
            self.cancel_button=W.QPushButton('stop');self.cancel_button.setEnabled(False);self.cancel_button.clicked.connect(self.cancel)
            self.open_button=W.QPushButton('open result directory');self.open_button.setEnabled(False);self.open_button.clicked.connect(self.open_output)
            for button in [self.run_button,self.cancel_button,self.open_button]:buttons.addWidget(button)
            self.progress=W.QProgressBar();self.progress.setRange(0,1);self.progress.setValue(0);layout.addWidget(self.progress)
            self.status=W.QLabel('ready');layout.addWidget(self.status)
            self.log=W.QPlainTextEdit();self.log.setReadOnly(True);self.log.setMaximumBlockCount(10000);layout.addWidget(self.log,2)
            self.update_labels()
            self.estimate_timer=C.QTimer(self);self.estimate_timer.setSingleShot(True);self.estimate_timer.setInterval(200)
            self.estimate_timer.timeout.connect(self.update_sequence_estimate)
            for entry in [self.input,self.sequence]:entry.textChanged.connect(lambda *_:self.estimate_timer.start())
            self.copies.valueChanged.connect(lambda *_:self.estimate_timer.start())

        def help_button(self,key,label,text):
            button=W.QToolButton();button.setText('?');button.setFixedWidth(25)
            button.setAccessibleName(label+' help');button.setToolTip(text)
            button.clicked.connect(lambda checked=False:W.QMessageBox.information(self,label+' — input help',text))
            self.help_buttons[key]=button
            return button

        def field(self,form,label,widget,key,text):
            row=W.QWidget();box=W.QHBoxLayout(row);box.setContentsMargins(0,0,0,0)
            box.addWidget(widget,1);box.addWidget(self.help_button(key,label,text));form.addRow(label,row)

        def update_sequence_estimate(self):
            self.sequence_estimate=None;self.solvent.clear();self.sequence_info.setStyleSheet('')
            if not self.sequence.text().strip() or not self.input.text().strip():
                self.sequence_info.setText('waiting for sequence and cell information');return
            try:
                m=gemmi.read_mtz_file(str(Path(self.input.text().strip()).expanduser()))
                info=sequence_asu(self.sequence.text().strip(),self.copies.value(),m.cell,m.spacegroup)
                self.sequence_estimate=info;self.solvent.setText(f"{info['solvent_fraction']:.2%} (Matthews estimate)")
                text=f"ASU: {info['asu_chains']} chains · {info['asu_residues']} residues · {info['asu_mass_da']/1000:.3f} kDa · Vm={info['matthews_vm']:.3f} Å³/Da"
                if info['warnings']:text+='\n'+'; '.join(info['warnings'])
                self.sequence_info.setText(text)
            except Exception as error:
                self.sequence_info.setText(str(error));self.sequence_info.setStyleSheet('color:#a33124')

        def file_row(self,form,label,mode,filter,inspect=False,help_key=None,help_text=None):
            row=W.QWidget();line=W.QLineEdit();button=W.QPushButton('select…')
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
                'experimental MTZ':('input','select the experimental X-ray SAD amplitude MTZ. A valid cell, space group, and anomalous amplitude pair or anomalous difference with standard deviations are required. Column names are read automatically after selection; phases from a complete target structure are never used.'),
                'heavy-atom coordinates':('heavy','select the heavy-atom substructure PDB/mmCIF. Elements, occupancies, cell, space group, and origin must match the experimental data; do not substitute a complete protein structure for the heavy-atom file. Only one model and isotropic B factors are accepted.'),
                'output MTZ':('output','choose a new .mtz output filename. A matching .json parameter record and .run.log log are saved as well. Existing results and input files cannot be overwritten.'),
                'optional custom reference MTZ':('reference','usually leave blank to use the 1TQW general reference obtained at installation. A custom reference needs a valid cell, space group, and the FC_BEST.F_phi.F (amplitude) and FC_BEST.F_phi.phi (phase) columns. Choose a suitable protein map unrelated to the target, not an ordinary experimental MTZ without phases.')}
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
                    box.clear();box.addItem('automatic','')
                    for name in cols:box.addItem(name,name)
                    for name in defaults[key]:
                        index=box.findData(name)
                        if index>=0:box.setCurrentIndex(index);break
                self.format.setCurrentIndex(1 if any(x in cols for x in ['DANO','DANOM']) else 0)
                self.data_info.setText(f'{m.spacegroup.hm} · {m.nreflections:,} reflections · highest resolution {m.resolution_high():.2f} Å')
                if not self.output.text():self.output.setText(str(Path(self.input.text()).with_name(Path(self.input.text()).stem+'_oasis.mtz')))
                self.update_sequence_estimate()
            except Exception as error:self.data_info.setText('failed to read MTZ: '+str(error))

        def calculate_fpp(self):
            try:
                import math
                wave=float(self.wavelength.text());el=gemmi.Element(self.element.text().strip())
                if not math.isfinite(wave) or wave<=0 or not el.atomic_number:raise ValueError('please provide a valid element and a positive wavelength')
                value=gemmi.cromer_liberman(el.atomic_number,gemmi.hc/wave)[1]
                self.fpp.setText(f'{el.name}={value:.7g}')
                self.status.setText("estimated with Cromer–Liberman theoretical values; check experimental values near an absorption edge")
            except Exception as error:W.QMessageBox.warning(self,"f'' estimation failed",str(error))

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
                self.log.appendPlainText('starting SAD phasing\nexperimental data: '+self.input.text()+'\nsequence: '+self.sequence.text()+'\nASU copies: '+str(self.copies.value())+'; estimated solvent content: '+self.solvent.text()+'\nreference map: '+('custom reference' if self.reference.text().strip() else 'default 1TQW general protein reference')+'\nfull parameters are saved in the run log.')
                self.run_log.write(header+'\n');self.run_log.flush()
                self.active_out=out;self.cancelled=False;self.started=time.monotonic()
                self.run_button.setEnabled(False);self.form_host.setEnabled(False);self.cancel_button.setEnabled(True);self.open_button.setEnabled(False)
                self.progress.setRange(0,0);self.status.setText('calculating; the finish time depends on the number of reflections and cycles…')
                self.process.start(command[0],command[1:])
            except Exception as error:W.QMessageBox.warning(self,'input check failed',str(error))

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
                self.status.setText(f'calculation completed · {time.monotonic()-self.started:.1f} s · {self.active_out}')
                self.log.appendPlainText('calculation completed; saved:\n'+str(self.active_out)+'\n'+str(self.active_out.with_suffix('.json'))+'\n'+str(self.active_out.with_suffix('.run.log')))
            elif self.cancelled:self.status.setText('calculation stopped; the run log is kept, output may be incomplete, please rerun with a new name.')
            else:self.status.setText(f'calculation failed (exit code {code}); check the log below.')

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
                answer=W.QMessageBox.question(self,'quit','calculation is still running; stop and quit?')
                if answer!=W.QMessageBox.StandardButton.Yes:event.ignore();return
                self.process.kill();self.process.waitForFinished(3000)
            event.accept()
    return Window()

if __name__=='__main__':raise SystemExit(main())
