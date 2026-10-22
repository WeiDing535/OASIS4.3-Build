"""Cancellation must reach Buccaneer even though it starts another session.

These tests launch only a bounded heartbeat executable; no crystallographic
program or structural input is involved.
"""
import json
import os
from pathlib import Path
import signal
import subprocess
import sys
import time

import pytest


def _eventually(predicate, timeout=6):
    deadline = time.monotonic() + timeout
    while time.monotonic() < deadline:
        if predicate():
            return
        time.sleep(0.02)
    raise AssertionError("Timed out waiting for the heartbeat process")


@pytest.mark.skipif(os.name != "posix", reason="POSIX process-group cancellation")
@pytest.mark.parametrize("stop_signal,ignore_term,gui_shutdown", [
    (signal.SIGTERM, False, False), (signal.SIGINT, False, False),
    (signal.SIGTERM, True, False), (signal.SIGTERM, True, True)])
def test_cli_cancellation_stops_independent_buccaneer_session(tmp_path, stop_signal, ignore_term, gui_shutdown):
    """The real recycle finally block terminates its separately grouped child."""
    root = tmp_path / "runtime"
    (root / "bin").mkdir(parents=True)
    library = root / "lib/data/reference_structures/reference-1tqw"
    library.parent.mkdir(parents=True)
    for suffix in (".mtz", ".pdb"):
        library.with_suffix(suffix).write_text("test resource\n")
    for name in ("map.mtz", "native.mtz", "seed.cif", "sequence.fasta"):
        (tmp_path / name).write_text("test input\n")

    # SIGTERM is observable, so stopping the heartbeat cannot be mistaken for a
    # startup failure or its independent twelve-second fail-safe expiry.
    executable = root / "bin/cbuccaneer"
    executable.write_text("#!" + sys.executable + "\n" + '''
import json, os, signal, time
from pathlib import Path
folder = Path(os.environ["HEARTBEAT_ROOT"])
def stop(signum, frame):
    (folder / "child_stopped").write_text(str(signum))
    raise SystemExit(0)
signal.signal(signal.SIGTERM, signal.SIG_IGN if os.environ["IGNORE_TERM"] == "1" else stop)
(folder / "child_ready.json").write_text(json.dumps({
    "pid": os.getpid(), "pgid": os.getpgrp(), "sid": os.getsid(0)}))
deadline = time.monotonic() + 12
counter = 0
while time.monotonic() < deadline:
    counter += 1
    (folder / "heartbeat").write_text(str(counter))
    time.sleep(0.025)
raise SystemExit("heartbeat fail-safe expired")
''')
    executable.chmod(0o755)
    source_root = Path(__file__).resolve().parents[1] / "src"
    driver = tmp_path / "driver.py"
    driver.write_text('''
import os, sys, time
from pathlib import Path
from phasebuild import hybrid_recycle as h
root = Path(os.environ["HEARTBEAT_ROOT"])
h.audit = lambda _: {"columns": [
    {"label": label, "type": kind} for label, kind in
    [("F", "F"), ("SIGF", "Q"), ("FWT", "F"), ("PHWT", "P"), ("FOM", "W")]]}
h.resolve_runtime = lambda *args, **kwargs: root / "runtime"
h.runtime_info = lambda _: {"test_stub": True}
h.runtime_environment = lambda _, **kwargs: dict(os.environ)
def unexpected_refinement(*args, **kwargs):
    raise AssertionError("Cancellation test must not reach refinement")
h.refmac_refine = unexpected_refinement
if os.environ["GUI_SHUTDOWN"] == "1" and "--engine" not in sys.argv:
    # Replace only the socket-facing server. Actual gui.main signal/finally,
    # Jobs worker, CLI process group and engine child session all execute.
    from types import SimpleNamespace
    from phasebuild import gui
    (root / "run_phasebuild.py").touch()
    jobs = gui.Jobs(root)
    gui.make_plan = lambda config, project: ({"output": str(root / "gui_result")},
        [{"name": "heartbeat engine", "kind": "hybrid", "output": str(root / "result"),
          "command": [sys.executable, __file__, "--engine"]}])
    def serve_forever(**kwargs):
        jobs.start({})
        while True:
            time.sleep(0.025)
    server = SimpleNamespace(jobs=jobs, address="127.0.0.1:0", token="test",
        serve_forever=serve_forever,
        server_close=lambda: (root / "server_closed").write_text("closed"))
    gui.create_server = lambda *args: server
    sys.argv = [__file__, "--no-browser", "--project", str(root)]
    gui.main()
else:
    h.recycle(root / "map.mtz", root / "native.mtz", root / "seed.cif",
              root / "sequence.fasta", root / "result", rounds=1, progress=lambda _: None)
''')
    environment = dict(os.environ, HEARTBEAT_ROOT=str(tmp_path),
                       IGNORE_TERM="1" if ignore_term else "0",
                       GUI_SHUTDOWN="1" if gui_shutdown else "0",
                       PYTHONPATH=str(source_root))
    child_pid = None
    with (tmp_path / "driver.log").open("w") as log:
        process = subprocess.Popen([sys.executable, str(driver)], env=environment,
                                   stdout=log, stderr=subprocess.STDOUT,
                                   start_new_session=True)
        try:
            _eventually(lambda: (tmp_path / "child_ready.json").is_file()
                        and (tmp_path / "heartbeat").is_file())
            identity = json.loads((tmp_path / "child_ready.json").read_text())
            child_pid = identity["pid"]
            assert identity["pgid"] == identity["sid"] == child_pid
            assert child_pid != process.pid
            cancelled = time.monotonic()
            os.killpg(process.pid, stop_signal)
            code = process.wait(timeout=3)
            assert code == 0 if gui_shutdown else code != 0
            assert time.monotonic() - cancelled < 3
            if gui_shutdown:
                assert (tmp_path / "server_closed").is_file()
                jobs = list((tmp_path / ".gui/jobs").glob("*.json"))
                assert len(jobs) == 1
                assert json.loads(jobs[0].read_text())["status"] == "cancelled"
            if ignore_term:
                assert not (tmp_path / "child_stopped").exists()
            else:
                assert int((tmp_path / "child_stopped").read_text()) == signal.SIGTERM
            report = json.loads((tmp_path / "result/hybrid_report.json").read_text())
            assert report["status"] == "failed"
            assert not report["stages"]
            # recycle() waited for its child; a successful cancellation leaves
            # neither a live child nor just a frozen heartbeat file as evidence.
            with pytest.raises(ProcessLookupError):
                os.kill(child_pid, 0)
            heartbeat = (tmp_path / "heartbeat").read_text()
            time.sleep(0.08)
            assert (tmp_path / "heartbeat").read_text() == heartbeat
        finally:
            if process.poll() is None:
                os.killpg(process.pid, signal.SIGKILL)
                process.wait(timeout=3)
            if child_pid is not None:
                try:
                    os.killpg(child_pid, signal.SIGKILL)
                except ProcessLookupError:
                    pass
