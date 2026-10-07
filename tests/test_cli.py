"""Newcomer entry point must not guess storage, claim unintended tasks or hide failures."""
import json,socket,subprocess,sys
from pathlib import Path
import pytest
from localmolbio.cli import main,inspect_environment
from localmolbio import worker
from localmolbio.database import initialise
from test_alignment_evidence import minimap2


def test_doctor_read_only_and_explicit_storage(tmp_path,monkeypatch,capsys):
    directory=tmp_path/'not-created'/'data'
    monkeypatch.delenv('MOLBIO_DATA_DIR',raising=False);monkeypatch.delenv('MOLBIO_MINIMAP2',raising=False)
    with pytest.raises(SystemExit) as exc:main(['doctor'])
    assert exc.value.code==2 and not directory.exists()
    assert main(['--data-dir',str(directory),'doctor'])==0
    result=json.loads(capsys.readouterr().out);assert result['ready'] and not result['alignment']['ready']
    assert result['data_directory']==str(directory) and not result['database_exists']
    assert not directory.exists()
    assert main(['--data-dir',str(directory),'doctor','--require-alignment'])==2
    assert not directory.exists()


def test_doctor_validated_real_tool_and_bad_path(tmp_path,minimap2):
    r=inspect_environment(tmp_path/'data',str(minimap2));assert r['alignment']['ready']
    assert r['alignment']['version']=='2.31-r1302'
    assert not inspect_environment(tmp_path,str(tmp_path/'missing'))['alignment']['ready']
    file=tmp_path/'file';file.write_text('preserve')
    assert not inspect_environment(file)['ready'] and file.read_text()=='preserve'


def test_dependency_failure_explained_without_writes(tmp_path,monkeypatch):
    from localmolbio import cli
    real=cli.importlib.import_module
    def broken(module):
        if module=='primer3':raise OSError('ABI issue')
        return real(module)
    monkeypatch.setattr(cli.importlib,'import_module',broken)
    r=inspect_environment(tmp_path/'new');assert not r['ready']
    assert any(c['name']=='primer3-py' and not c['ready'] for c in r['checks'])
    assert not (tmp_path/'new').exists()


def test_busy_port_rejected_before_database_creation(tmp_path,monkeypatch,capsys):
    monkeypatch.delenv('MOLBIO_MINIMAP2',raising=False)
    with socket.socket() as s:
        s.bind(('127.0.0.1',0));s.listen(1)
        assert main(['--data-dir',str(tmp_path/'data'),'serve','--port',str(s.getsockname()[1])])==2
    assert 'port_unavailable' in capsys.readouterr().out and not (tmp_path/'data').exists()


@pytest.mark.parametrize('statuses,limit,calls,code',[
    (['idle'],5,1,0),(['succeeded','idle'],5,2,0),(['succeeded','succeeded'],2,2,0),
    (['attempt_failed','succeeded'],5,1,1),(['lease_lost','succeeded'],5,1,1)])
def test_finite_worker_stops_at_limit_idle_or_first_failure(tmp_path,monkeypatch,statuses,limit,calls,code):
    monkeypatch.setenv('MOLBIO_DATA_DIR',str(tmp_path));initialise()
    observed=[]
    def once(worker_id,**kw):
        observed.append(kw['adapter']);return {'status':statuses[len(observed)-1]}
    monkeypatch.setattr(worker,'run_once',once)
    assert main(['--data-dir',str(tmp_path),'work','--adapter','fastq','--max-jobs',str(limit)])==code
    assert observed==[worker.fastq_validation.ADAPTER]*calls


def test_work_preflight_does_not_claim_or_create_missing_database(tmp_path,monkeypatch):
    monkeypatch.delenv('MOLBIO_MINIMAP2',raising=False)
    def unexpected(*a,**kw):raise AssertionError('must not claim')
    monkeypatch.setattr(worker,'run_once',unexpected)
    assert main(['--data-dir',str(tmp_path/'missing'),'work','--adapter','alignment'])==2
    assert main(['--data-dir',str(tmp_path/'missing'),'work','--adapter','ab1'])==2
    assert not (tmp_path/'missing').exists()


@pytest.mark.parametrize('tail',[['work'],['work','--adapter','any'],['work','--adapter','fastq','--max-jobs','0'],['work','--adapter','fastq','--max-jobs','101'],['serve','--host','0.0.0.0'],['serve','--port','65536']])
def test_explicit_adapter_and_loopback_only(tmp_path,tail):
    with pytest.raises(SystemExit) as exc:main(['--data-dir',str(tmp_path),*tail])
    assert exc.value.code==2


def test_python_module_help():
    r=subprocess.run([sys.executable,'-m','localmolbio','--help'],capture_output=True,text=True)
    assert r.returncode==0 and 'doctor' in r.stdout and '--data-dir' in r.stdout
