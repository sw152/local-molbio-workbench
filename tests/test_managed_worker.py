"""Foreground lifecycle and queue fences; real alignment plus explicit lifecycle doubles."""
import json,os,signal,subprocess,sys,threading,time
from pathlib import Path
import pytest
from localmolbio.managed_worker import ManagedAlignmentWorker
from localmolbio.database import connect
from localmolbio.config import data_dir
from localmolbio import job_queue as q
from localmolbio.cli import main
from test_sanger_history import ready
from test_alignment_evidence import minimap2
from test_fastq_alignment import setup,register,enqueue


def wait_for(check,seconds=15):
    deadline=time.monotonic()+seconds
    while time.monotonic()<deadline:
        if check():return
        time.sleep(.02)
    raise AssertionError('condition did not become true')


def status(job):
    with connect() as db:return db.execute('SELECT status FROM analysis_jobs WHERE id=?',(job,)).fetchone()[0]


def test_real_worker_runs_published_alignment_and_duplicate_lock_blocks(setup):
    client,revision,ref=setup;item=register(revision,[ref[40:500]])
    job=enqueue(revision,[item]);foreign=q.enqueue(revision,'foreign','unsupported',{})
    runner=ManagedAlignmentWorker(data_dir());runner.acquire()
    other=ManagedAlignmentWorker(data_dir())
    try:
        with pytest.raises(RuntimeError,match='already_active'):other.acquire()
        runner.start();wait_for(lambda:status(job)=='succeeded')
        assert status(foreign)=='queued'
        result=client.get(f'/api/revisions/{revision}/alignments/{job}/reads').json()
        assert result['items'][0]['alignments'][0]['reference_intervals']==[[40,500]]
        wait_for(lambda:runner.status()['completed_tasks']==1)
        runner.stop();assert not runner.process.is_alive() and runner.status()['state']=='stopped'
    finally:runner.release();other.release()
    successor=ManagedAlignmentWorker(data_dir());successor.acquire();successor.release()


def test_first_failure_stops_before_next_job_and_does_not_retry(setup):
    _,revision,ref=setup;bad=register(revision,[ref[:100]+'R'+ref[101:500]])
    failed=enqueue(revision,[bad],key='bad');good=register(revision,[ref[40:500]]);later=enqueue(revision,[good],key='good')
    runner=ManagedAlignmentWorker(data_dir());runner.acquire()
    try:
        runner.start();wait_for(lambda:runner.status()['state']=='failed')
        assert status(failed)=='failed' and status(later)=='queued'
        with connect() as db:assert db.execute('SELECT count(*) FROM job_attempts WHERE job_id=?',(failed,)).fetchone()[0]==1
    finally:runner.release()


def paused_tool(tmp_path,monkeypatch):
    # Lifecycle-only test double, not a biological alignment benchmark.
    marker=tmp_path/'tool-started';release=tmp_path/'tool-release';tool=tmp_path/'pause-minimap2'
    tool.write_text(f'#!{sys.executable}\nimport os,sys,time\nfrom pathlib import Path\nif "--version" in sys.argv:\n print("2.31-r1302");sys.exit(0)\nPath({str(marker)!r}).write_text(str(os.getpid()))\nfor _ in range(2000):\n if Path({str(release)!r}).exists():break\n time.sleep(.05)\n')
    tool.chmod(0o755);monkeypatch.setenv('MOLBIO_MINIMAP2',str(tool))
    return marker,release


def test_shutdown_finishes_current_attempt_without_claiming_next(setup,tmp_path,monkeypatch):
    _,revision,ref=setup;marker,release=paused_tool(tmp_path,monkeypatch)
    item=register(revision,[ref[40:500]]);first=enqueue(revision,[item],key='first');second=enqueue(revision,[item],key='second')
    runner=ManagedAlignmentWorker(data_dir());runner.acquire();thread=None
    try:
        runner.start();wait_for(marker.exists)
        thread=threading.Thread(target=runner.stop);thread.start()
        wait_for(lambda:runner.status()['state']=='stopping');assert status(second)=='queued'
        release.touch();thread.join(10);assert not thread.is_alive()
        assert status(first)=='succeeded' and status(second)=='queued' and not runner.process.is_alive()
    finally:release.touch();runner.release()


def test_forced_shutdown_kills_aligner_group_and_leaves_unpublished_lease(setup,tmp_path,monkeypatch):
    _,revision,ref=setup;marker,release=paused_tool(tmp_path,monkeypatch)
    item=register(revision,[ref[40:500]]);job=enqueue(revision,[item])
    runner=ManagedAlignmentWorker(data_dir(),shutdown_timeout=.1);runner.acquire()
    try:
        runner.start();wait_for(marker.exists);tool_pid=int(marker.read_text())
        runner.stop();assert not runner.process.is_alive()
        assert status(job)=='running'  # Expiring lease, never a fabricated finished report.
        with connect() as db:assert json.loads(db.execute('SELECT result_summary_json FROM analysis_jobs WHERE id=?',(job,)).fetchone()[0])=={}
        assert setup[0].get(f'/api/revisions/{revision}/alignments/{job}/report').status_code==409
        process=subprocess.run(['ps','-p',str(tool_pid),'-o','stat='],capture_output=True,text=True)
        assert not process.stdout.strip() or process.stdout.strip().startswith('Z'),process.stdout
    finally:runner.release()


def test_opt_in_missing_tool_refuses_before_initializing_or_claiming(tmp_path,monkeypatch):
    monkeypatch.delenv('MOLBIO_MINIMAP2',raising=False)
    directory=tmp_path/'unused'
    assert main(['--data-dir',str(directory),'serve','--with-alignment-worker'])==2
    assert not directory.exists()


def test_api_default_has_no_worker_or_controls(setup):
    client,revision,ref=setup
    result=client.get('/api/runtime').json()['alignment_worker']
    assert result=={'enabled':False,'state':'disabled','completed_tasks':0,'adapter':'alignment'}
    assert client.post('/api/runtime').status_code==405


def test_database_startup_failure_never_starts_worker(monkeypatch):
    import asyncio
    from types import SimpleNamespace
    from localmolbio import api
    calls=[]
    runner=SimpleNamespace(start=lambda:calls.append('start'),stop=lambda:calls.append('stop'))
    def broken():raise RuntimeError('startup failure')
    monkeypatch.setattr(api,'initialise',broken)
    async def attempt():
        async with api.lifespan(SimpleNamespace(state=SimpleNamespace(alignment_runner=runner))):pass
    with pytest.raises(RuntimeError,match='startup failure'):asyncio.run(attempt())
    assert calls==[]


def test_parent_death_aborts_inflight_child_without_publishing(setup,tmp_path,monkeypatch):
    _,revision,ref=setup;marker,release=paused_tool(tmp_path,monkeypatch)
    item=register(revision,[ref[40:500]]);job=enqueue(revision,[item])
    pidfile=tmp_path/'worker-pid'
    command=('from localmolbio.managed_worker import ManagedAlignmentWorker\n'
             'from localmolbio.config import data_dir\nfrom pathlib import Path\nimport time\n'
             f'r=ManagedAlignmentWorker(data_dir());r.acquire();r.start();Path({str(pidfile)!r}).write_text(str(r.process.pid));time.sleep(30)\n')
    parent=subprocess.Popen([sys.executable,'-c',command],stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL)
    try:
        wait_for(marker.exists);wait_for(pidfile.exists)
        parent.kill();parent.wait(timeout=5)
        wait_for(lambda:status(job)=='queued')  # Retryable parent loss, not a success.
        worker_pid=int(pidfile.read_text());tool_pid=int(marker.read_text())
        def dead(pid):
            r=subprocess.run(['ps','-p',str(pid),'-o','stat='],capture_output=True,text=True)
            return not r.stdout.strip() or r.stdout.strip().startswith('Z')
        wait_for(lambda:dead(worker_pid) and dead(tool_pid))
        assert setup[0].get(f'/api/revisions/{revision}/alignments/{job}/report').status_code==409
    finally:
        if parent.poll() is None:parent.kill();parent.wait(timeout=5)
        release.touch()


def test_worker_crash_reaps_aligner_without_waiting_for_server_shutdown(setup,tmp_path,monkeypatch):
    _,revision,ref=setup;marker,release=paused_tool(tmp_path,monkeypatch)
    item=register(revision,[ref[40:500]]);job=enqueue(revision,[item])
    runner=ManagedAlignmentWorker(data_dir());runner.acquire()
    try:
        runner.start();wait_for(marker.exists);tool_pid=int(marker.read_text())
        runner.process.kill()
        wait_for(lambda:runner.status()['state']=='failed')
        def dead():
            r=subprocess.run(['ps','-p',str(tool_pid),'-o','stat='],capture_output=True,text=True)
            return not r.stdout.strip() or r.stdout.strip().startswith('Z')
        wait_for(dead)
        assert status(job)=='running'
        assert setup[0].get(f'/api/revisions/{revision}/alignments/{job}/reads').status_code==409
    finally:runner.release()


def test_process_spawn_failure_can_release_lock(tmp_path,monkeypatch):
    runner=ManagedAlignmentWorker(tmp_path);runner.acquire()
    def fail(*args,**kwargs):raise RuntimeError('spawn failed')
    monkeypatch.setattr(runner.context.Process,'start',fail)
    try:
        with pytest.raises(RuntimeError,match='spawn failed'):runner.start()
    finally:runner.release()
    other=ManagedAlignmentWorker(tmp_path);other.acquire();other.release()
