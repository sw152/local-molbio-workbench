"""Registered FASTQ identity, original compressed bytes and fenced queue results."""
from hashlib import sha256
import gzip
import io
import json
import os
from pathlib import Path
import subprocess
import sys

import pytest
from localmolbio.database import connect
from localmolbio.fastq_inputs import register_fastq
from localmolbio import fastq_validation as checks, input_validation as ab1, worker, job_queue as q
from test_fastq_inputs import DATA
from test_input_worker import stored
from test_sanger_history import ready
from test_sanger_acceptance import sibling_revision


def register(revision, data=DATA, compressed=False):
    content=gzip.compress(data,mtime=0) if compressed else data
    return register_fastq(revision,'input.fq.gz' if compressed else 'input.fq',io.BytesIO(content),'phred33')


def path(identity):
    with connect() as db:return Path(db.execute('SELECT storage_path FROM fastq_inputs WHERE id=?',(identity,)).fetchone()[0])


def run(name='fastq-worker'):
    return worker.run_once(name,adapter=checks.ADAPTER)


@pytest.mark.parametrize('compressed',[False,True])
def test_registered_original_bytes_and_metadata_with_no_analysis(ready,compressed):
    client,read,revision=ready;item=register(revision,compressed=compressed);original=path(item['id']).read_bytes()
    job=checks.enqueue_check(revision,[item['id']],'identity')
    assert checks.enqueue_check(revision,[item['id']],'identity')==job
    assert worker.run_once('default-ab1')=={'status':'idle'}
    assert stored(job)['status']=='queued'
    assert run()['status']=='succeeded'
    result=json.loads(stored(job)['result_summary_json'])
    assert result['scope']=='registered_fastq_integrity_only'
    assert result['analysis_performed'] is result['whole_reference_verified'] is result['summary_recomputed'] is False
    entry=result['files'][0]
    assert entry['input_id']==item['id'] and entry['sha256']==sha256(original).hexdigest()
    assert entry['size_bytes']==len(original) and entry['quality_encoding']=='phred33'
    assert entry['compression']==('gzip' if compressed else 'none')
    assert len(entry['summary_sha256'])==64
    assert 'storage_path' not in json.dumps(result) and path(item['id']).read_bytes()==original
    assert client.get(f'/api/sanger-reads/{read}/analyses').json()['total']==0
    assert run()=={'status':'idle'}


@pytest.mark.parametrize('change,code',[
    ('bytes','input_hash_mismatch'),('size','input_size_mismatch'),('missing','input_file_missing'),
    ('hash','fastq_identity_changed_since_submission'),('declared_size','fastq_identity_changed_since_submission'),
    ('compression','unsupported_fastq_registration_metadata'),('encoding','unsupported_fastq_registration_metadata'),
    ('summary','fastq_identity_changed_since_submission'),('summary_invalid','invalid_fastq_registration_summary'),
    ('membership','fastq_missing_or_wrong_revision'),('reference_bytes','reference_hash_mismatch'),
    ('reference_topology','reference_changed_since_submission'),
])
def test_changes_after_submission_never_publish_success(ready,change,code):
    _,_,revision=ready;item=register(revision);identity=item['id'];job=checks.enqueue_check(revision,[identity],'changed');p=path(identity)
    if change=='bytes':p.write_bytes(p.read_bytes().replace(b'ACGN',b'TCGN'))
    elif change=='size':p.write_bytes(p.read_bytes()+b'X')
    elif change=='missing':p.unlink()
    else:
        other=sibling_revision(revision) if change=='membership' else None
        with connect() as db:
            if change=='hash':db.execute('UPDATE fastq_inputs SET file_sha256=? WHERE id=?',('0'*64,identity))
            elif change=='declared_size':db.execute('UPDATE fastq_inputs SET size_bytes=size_bytes+1 WHERE id=?',(identity,))
            elif change=='compression':db.execute("UPDATE fastq_inputs SET compression='gzip' WHERE id=?",(identity,))
            elif change=='encoding':
                db.execute('PRAGMA ignore_check_constraints=ON')
                db.execute("UPDATE fastq_inputs SET quality_encoding='phred64' WHERE id=?",(identity,))
            elif change.startswith('summary'):
                summary=dict(item['summary']);summary['bases']+=1
                db.execute('UPDATE fastq_inputs SET summary_json=? WHERE id=?',('{' if change=='summary_invalid' else json.dumps(summary),identity))
            elif change=='membership':db.execute('UPDATE fastq_inputs SET sequence_revision_id=? WHERE id=?',(other,identity))
            elif change=='reference_bytes':db.execute("UPDATE sequence_revisions SET sequence_text='ACGT' WHERE id=?",(revision,))
            else:db.execute("UPDATE sequence_revisions SET topology='circular' WHERE id=?",(revision,))
    result=run();row=stored(job)
    assert result['error']==code and row['status']=='failed' and row['result_summary_json']=='{}'
    assert row['attempts'][0]['error_detail']==code


@pytest.mark.parametrize('mode,code',[('outside','input_outside_managed_fastq'),('symlink','symlink_input_rejected'),('directory','input_not_regular_file'),('fifo','input_not_regular_file')])
def test_fastq_paths_remain_within_managed_regular_files(ready,tmp_path,mode,code):
    _,_,revision=ready;item=register(revision);identity=item['id'];p=path(identity)
    outside=tmp_path/'outside.fq';outside.write_bytes(p.read_bytes())
    if mode=='outside':
        with connect() as db:db.execute('UPDATE fastq_inputs SET storage_path=? WHERE id=?',(str(outside),identity))
    elif mode=='symlink':p.unlink();p.symlink_to(outside)
    elif mode=='directory':p.unlink();p.mkdir()
    else:p.unlink();os.mkfifo(p)
    job=checks.enqueue_check(revision,[identity],'path');assert run()['error']==code
    assert stored(job)['result_summary_json']=='{}'


def test_canonical_multifile_request_and_no_partial_results(ready):
    _,_,revision=ready;a=register(revision);b=register(revision,compressed=True)
    ids=sorted([a['id'],b['id']]);job=checks.enqueue_check(revision,ids,'multi')
    assert checks.enqueue_check(revision,list(reversed(ids)),'multi')==job
    path(ids[-1]).write_bytes(b'broken')
    assert run()['status']=='attempt_failed' and stored(job)['result_summary_json']=='{}'
    other=checks.enqueue_check(revision,[ids[0]],'single');assert run()['status']=='succeeded'
    assert len(json.loads(stored(other)['result_summary_json'])['files'])==1


def test_summary_digest_ignores_json_spacing_but_not_values(ready):
    _,_,revision=ready;item=register(revision);identity=item['id'];job=checks.enqueue_check(revision,[identity],'stable')
    with connect() as db:db.execute('UPDATE fastq_inputs SET summary_json=? WHERE id=?',(json.dumps(item['summary'],sort_keys=True,indent=4),identity))
    assert checks.enqueue_check(revision,[identity],'stable')==job
    assert run()['status']=='succeeded'
    summary=dict(item['summary']);summary['q20_bases']-=1
    with connect() as db:db.execute('UPDATE fastq_inputs SET summary_json=? WHERE id=?',(json.dumps(summary),identity))
    with pytest.raises(q.QueueConflict):checks.enqueue_check(revision,[identity],'stable')


def test_midstream_file_or_registration_changes_are_rejected(ready,monkeypatch):
    _,_,revision=ready;item=register(revision);identity=item['id'];job=checks.enqueue_check(revision,[identity],'during')
    claim=q.claim('fixture',checks.ADAPTER);calls=[0];p=path(identity);original=p.read_bytes()
    def pulse():
        calls[0]+=1
        if calls[0]==3:p.write_bytes(original+b'changed')
    with pytest.raises(checks.InputValidationError,match='input_changed_during_check'):checks.validate_inputs(claim,pulse)
    p.write_bytes(original)
    original_hash=checks._hash_managed_file
    def change_registration(*args,**kwargs):
        result=original_hash(*args,**kwargs)
        with connect() as db:db.execute("UPDATE fastq_inputs SET original_filename='renamed.fq' WHERE id=?",(identity,))
        return result
    monkeypatch.setattr(checks,'_hash_managed_file',change_registration)
    with pytest.raises(checks.InputValidationError,match='registered_inputs_changed_during_check'):checks.validate_inputs(claim,lambda:None)
    assert stored(job)['result_summary_json']=='{}'


def test_transient_io_retry_and_attempt_budget(ready,monkeypatch):
    _,_,revision=ready;item=register(revision);job=checks.enqueue_check(revision,[item['id']],'io',max_attempts=2)
    def transient(*a,**kw):raise checks.InputValidationError('input_io_error',retryable=True)
    with monkeypatch.context() as m:
        m.setattr(checks,'_hash_managed_file',transient)
        assert run('first')['status']=='attempt_failed'
    assert stored(job)['status']=='queued';assert run('second')['status']=='succeeded'
    assert [r['status'] for r in stored(job)['attempts']]==['failed','succeeded']
    exhausted=checks.enqueue_check(revision,[item['id']],'budget',max_attempts=1)
    monkeypatch.setattr(checks,'_hash_managed_file',transient)
    assert run()['status']=='attempt_failed' and stored(exhausted)['status']=='failed'


@pytest.mark.parametrize('change',['cancel','expire'])
def test_cancelled_or_expired_worker_cannot_publish(ready,monkeypatch,change):
    _,_,revision=ready;item=register(revision);job=checks.enqueue_check(revision,[item['id']],'lease')
    original=checks.validate_inputs
    def lose(claim,pulse):
        result=original(claim,pulse)
        if change=='cancel':q.cancel(job)
        else:
            with connect() as db:db.execute('UPDATE queued_jobs SET lease_expires=0 WHERE job_id=?',(job,))
        return result
    with monkeypatch.context() as m:
        m.setattr(checks,'validate_inputs',lose)
        assert run()['status']=='lease_lost'
    assert stored(job)['result_summary_json']=='{}'
    if change=='expire':
        assert q.recover_expired()==1;assert run('replacement')['status']=='succeeded'
        assert [a['status'] for a in stored(job)['attempts']]==['expired','succeeded']
    else:assert stored(job)['status']=='cancelled' and run()=={'status':'idle'}


def test_submission_isolation_validation_and_unknown_adapter(ready):
    _,read,revision=ready;item=register(revision);other=sibling_revision(revision)
    for ids in [[],[item['id']]*2,[read],['missing']]:
        with pytest.raises(checks.InputValidationError):checks.enqueue_check(revision,ids,'bad')
    with pytest.raises(checks.InputValidationError):checks.enqueue_check(other,[item['id']],'foreign')
    ab1_job=ab1.enqueue_check(revision,[read],'ab1')
    fastq_job=checks.enqueue_check(revision,[item['id']],'fastq')
    assert run()['job_id']==fastq_job and stored(ab1_job)['status']=='queued'
    assert worker.run_once('ab1')['job_id']==ab1_job
    with pytest.raises(ValueError,match='Unsupported'):worker.run_once('bad',adapter='shell')


def test_cli_independent_processes_and_redacted_unknown_failure(ready,monkeypatch):
    _,_,revision=ready;item=register(revision,compressed=True)
    def cli(*args):return subprocess.run([sys.executable,'-m','localmolbio.worker',*args],capture_output=True,text=True)
    submitted=cli('enqueue-fastq-check','--revision',revision,'--input-id',item['id'],'--key','cli')
    assert submitted.returncode==0,submitted.stderr
    job=json.loads(submitted.stdout)['job_id']
    assert json.loads(cli('once').stdout)['status']=='idle'
    completed=cli('once','--adapter','fastq','--worker-id','fresh')
    assert completed.returncode==0,completed.stderr
    assert json.loads(completed.stdout)['status']=='succeeded' and stored(job)['attempts'][0]['worker_id']=='fresh'
    broken=checks.enqueue_check(revision,[item['id']],'unknown')
    def fail(*args):raise RuntimeError('private sequence or filename')
    monkeypatch.setattr(checks,'validate_inputs',fail)
    assert run()['error']=='unexpected_input_check_error'
    assert stored(broken)['status']=='failed' and 'private' not in json.dumps(stored(broken))
