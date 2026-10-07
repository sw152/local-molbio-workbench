"""Real registered-file checks, CLI execution and lease-safe failure paths."""
from hashlib import sha256
import json
import os
from pathlib import Path
import subprocess
import sys

import pytest
from localmolbio.database import connect
from localmolbio import job_queue as q
from localmolbio import input_validation as checks
from localmolbio import worker
from test_sanger_history import ready
from test_sanger_acceptance import sibling_revision


def stored(job):
    with connect() as db:
        row=dict(db.execute('SELECT * FROM analysis_jobs WHERE id=?',(job,)).fetchone())
        row['attempts']=[dict(r) for r in db.execute('SELECT * FROM job_attempts WHERE job_id=? ORDER BY number',(job,))]
    return row


def file_path(read):
    with connect() as db:return Path(db.execute('SELECT storage_path FROM sequencing_reads WHERE id=?',(read,)).fetchone()[0])


def test_success_checks_real_bytes_without_creating_alignment(ready):
    client,read,revision=ready
    original=file_path(read).read_bytes()
    job=checks.enqueue_check(revision,[read],'one')
    assert checks.enqueue_check(revision,[read],'one')==job
    outcome=worker.run_once('fixture')
    assert outcome['status']=='succeeded' and not outcome['analysis_performed']
    row=stored(job);result=json.loads(row['result_summary_json'])
    assert result['scope']=='registered_input_integrity_only' and not result['whole_reference_verified']
    assert result['files']==[{'read_id':read,'sha256':sha256(original).hexdigest(),'size_bytes':len(original)}]
    assert 'storage_path' not in row['result_summary_json']
    assert client.get(f'/api/sanger-reads/{read}/analyses').json()['total']==0
    assert file_path(read).read_bytes()==original and worker.run_once('idle')=={'status':'idle'}


@pytest.mark.parametrize('change,code',[
    ('missing','input_file_missing'),('bytes','input_hash_mismatch'),
    ('reference_bytes','reference_hash_mismatch'),('reference_identity','reference_changed_since_submission'),
    ('read_hash','read_identity_changed_since_submission'),('membership','read_missing_or_wrong_revision'),
])
def test_changed_inputs_fail_without_a_success_result(ready,change,code):
    _,read,revision=ready;job=checks.enqueue_check(revision,[read],'change')
    if change=='missing':file_path(read).unlink()
    elif change=='bytes':file_path(read).write_bytes(b'changed')
    elif change=='membership':
        other=sibling_revision(revision)
        with connect() as db:db.execute('UPDATE sequencing_reads SET sequence_revision_id=? WHERE id=?',(other,read))
    else:
        with connect() as db:
            if change=='reference_bytes':db.execute("UPDATE sequence_revisions SET sequence_text='ACGT' WHERE id=?",(revision,))
            elif change=='reference_identity':db.execute("UPDATE sequence_revisions SET topology='circular' WHERE id=?",(revision,))
            else:db.execute('UPDATE sequencing_reads SET file_sha256=? WHERE id=?',('0'*64,read))
    response=worker.run_once('fixture');row=stored(job)
    assert response['status']=='attempt_failed' and response['error']==code
    assert row['status']=='failed' and row['result_summary_json']=='{}'
    assert row['attempts'][0]['error_detail']==code


@pytest.mark.parametrize('mode,code',[('outside','input_outside_managed_reads'),('symlink','symlink_input_rejected'),('directory','input_not_regular_file'),('fifo','input_not_regular_file')])
def test_registered_paths_cannot_expand_worker_file_scope(ready,tmp_path,mode,code):
    _,read,revision=ready;path=file_path(read)
    outside=tmp_path/'outside.ab1';outside.write_bytes(path.read_bytes())
    if mode=='outside':
        with connect() as db:db.execute('UPDATE sequencing_reads SET storage_path=? WHERE id=?',(str(outside),read))
    elif mode=='symlink':path.unlink();path.symlink_to(outside)
    elif mode=='directory':path.unlink();path.mkdir()
    else:path.unlink();os.mkfifo(path)
    job=checks.enqueue_check(revision,[read],'path')
    assert worker.run_once('fixture')['error']==code
    assert stored(job)['status']=='failed'


def test_mid_read_change_is_not_published(ready,monkeypatch):
    _,read,revision=ready;job=checks.enqueue_check(revision,[read],'mutation')
    claim=q.claim('fixture',checks.ADAPTER);path=file_path(read);original=path.read_bytes();calls=[0]
    def pulse():
        calls[0]+=1
        if calls[0]==3:path.write_bytes(original+b'changed')
    with pytest.raises(checks.InputValidationError,match='input_changed_during_check'):
        checks.validate_inputs(claim,pulse)
    assert stored(job)['result_summary_json']=='{}'


def test_transient_io_retries_but_preserves_attempt_history(ready,monkeypatch):
    _,read,revision=ready;job=checks.enqueue_check(revision,[read],'io',max_attempts=2)
    def io_error(*args):raise checks.InputValidationError('input_io_error',retryable=True)
    with monkeypatch.context() as m:
        m.setattr(checks,'_hash_registered_file',io_error)
        assert worker.run_once('first')['status']=='attempt_failed'
    assert stored(job)['status']=='queued'
    assert worker.run_once('second')['status']=='succeeded'
    assert [a['status'] for a in stored(job)['attempts']]==['failed','succeeded']


@pytest.mark.parametrize('change',['cancel','expire'])
def test_lost_lease_never_publishes_or_overwrites_job_status(ready,monkeypatch,change):
    _,read,revision=ready;job=checks.enqueue_check(revision,[read],'lease')
    real_validate=worker.validate_inputs
    def lose(claim,pulse):
        result=real_validate(claim,pulse)
        if change=='cancel':q.cancel(job)
        else:
            with connect() as db:db.execute('UPDATE queued_jobs SET lease_expires=0 WHERE job_id=?',(job,))
        return result
    monkeypatch.setattr(worker,'validate_inputs',lose)
    assert worker.run_once('late')['status']=='lease_lost'
    assert stored(job)['result_summary_json']=='{}'
    assert stored(job)['status']==('cancelled' if change=='cancel' else 'running')
    if change=='expire':assert q.recover_expired()==1 and stored(job)['status']=='queued'


def test_unknown_exceptions_are_redacted_and_terminal(ready,monkeypatch):
    _,read,revision=ready;job=checks.enqueue_check(revision,[read],'error')
    def bad(*args):raise RuntimeError('private filename and sequence must not be logged')
    monkeypatch.setattr(worker,'validate_inputs',bad)
    assert worker.run_once('broken')['error']=='unexpected_input_check_error'
    row=stored(job)
    assert row['status']=='failed' and 'private' not in json.dumps(row)


def test_invalid_submission_and_unsupported_adapter_are_not_executed(ready):
    _,read,revision=ready
    for ids in [[],[read,read],['missing']]:
        with pytest.raises(checks.InputValidationError):checks.enqueue_check(revision,ids,'bad')
    with pytest.raises(checks.InputValidationError):checks.enqueue_check(sibling_revision(revision),[read],'foreign')
    job=q.enqueue(revision,'unsupported','shell-command',{'command':'not executable'})
    assert worker.run_once('restricted')=={'status':'idle'} and stored(job)['status']=='queued'


def test_cli_enqueue_and_once_in_fresh_processes(ready):
    _,read,revision=ready
    def cli(*args):
        return subprocess.run([sys.executable,'-m','localmolbio.worker',*args],capture_output=True,text=True)
    queued=cli('enqueue-check','--revision',revision,'--read-id',read,'--key','cli')
    assert queued.returncode==0,queued.stderr
    job=json.loads(queued.stdout)['job_id']
    completed=cli('once','--worker-id','cli-worker')
    assert completed.returncode==0,completed.stderr
    assert json.loads(completed.stdout)['status']=='succeeded'
    assert stored(job)['attempts'][0]['worker_id']=='cli-worker'
    assert json.loads(cli('once').stdout)['status']=='idle'
    rejected=cli('enqueue-check','--revision',revision,'--read-id','missing','--key','invalid')
    assert rejected.returncode==2 and json.loads(rejected.stdout)['status']=='rejected'

    file_path(read).write_bytes(b'changed after original check')
    failed_job=cli('enqueue-check','--revision',revision,'--read-id',read,'--key','cli-failure')
    assert failed_job.returncode==0
    failure=cli('once')
    assert failure.returncode==1 and json.loads(failure.stdout)['error']=='input_hash_mismatch'


@pytest.mark.parametrize('damaged',[False,True])
def test_multiple_inputs_are_canonical_and_never_publish_partial_success(ready,damaged):
    from abif_fixture import synthetic_ab1
    client,read,revision=ready
    with connect() as db:reference=db.execute('SELECT sequence_text FROM sequence_revisions WHERE id=?',(revision,)).fetchone()[0]
    response=client.post('/api/sanger-reads',data={'sequence_revision_id':revision},files={'file':('second.ab1',synthetic_ab1(reference[250:450],[35]*200,trace=True))})
    assert response.status_code==201
    other=response.json()['id'];ids=sorted([read,other])
    job=checks.enqueue_check(revision,[read,other],'multi')
    assert checks.enqueue_check(revision,[other,read],'multi')==job
    if damaged:file_path(ids[-1]).write_bytes(b'damaged second input')
    outcome=worker.run_once('multi-worker');row=stored(job)
    if damaged:
        assert outcome['error']=='input_hash_mismatch' and row['result_summary_json']=='{}'
    else:
        assert outcome['status']=='succeeded'
        assert [r['read_id'] for r in json.loads(row['result_summary_json'])['files']]==ids
