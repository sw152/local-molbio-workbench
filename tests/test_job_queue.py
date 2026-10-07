"""Queue transitions, contention, lease fencing and durable recovery."""
from concurrent.futures import ThreadPoolExecutor
import json
import subprocess
import sys
from threading import Barrier

import pytest
from localmolbio import job_queue as q
from localmolbio.database import connect, initialise
from test_sanger_history import ready


def submit(revision,key='request',**options):
    return q.enqueue(revision,key,'fixture-v1',{'files':[{'sha256':'a'*64,'id':'fixture'}]},**options)


def state(job):
    with connect() as db:
        row=dict(db.execute('SELECT * FROM analysis_jobs WHERE id=?',(job,)).fetchone())
        attempts=[dict(r) for r in db.execute('SELECT * FROM job_attempts WHERE job_id=? ORDER BY number',(job,))]
    return row,attempts


def test_same_key_is_atomic_and_changed_requests_conflict(ready):
    _,_,revision=ready
    barrier=Barrier(6)
    def submit_once(_):
        barrier.wait()
        return submit(revision)
    with ThreadPoolExecutor(max_workers=6) as pool:ids=list(pool.map(submit_once,range(6)))
    assert len(set(ids))==1
    with pytest.raises(q.QueueConflict):submit(revision,parameters={'changed':True})
    with pytest.raises(q.QueueConflict):submit(revision,max_attempts=1)
    with pytest.raises(q.QueueConflict):q.enqueue(revision,'request','different-v1',{})
    row,attempts=state(ids[0]);assert row['status']=='queued' and attempts==[]
    manifest=json.loads(row['input_manifest_json'])
    assert manifest['reference']['id']==revision and len(manifest['reference']['sequence_sha256'])==64
    with connect() as db:assert db.execute("SELECT count(*) FROM audit_events WHERE object_id=? AND action='queued'",(ids[0],)).fetchone()[0]==1


def test_only_one_worker_claims_and_nonqueue_jobs_are_untouched(ready):
    client,read,revision=ready
    direct=client.post('/api/sanger-verifications',json={'sequencing_read_id':read}).json()
    job=submit(revision);barrier=Barrier(6)
    def take(i):
        barrier.wait()
        return q.claim(f'worker-{i}','fixture-v1')
    with ThreadPoolExecutor(max_workers=6) as pool:claims=list(pool.map(take,range(6)))
    claimed=[c for c in claims if c]
    assert len(claimed)==1 and claimed[0]['job_id']==job
    assert q.claim('wrong-adapter','uninstalled-v1') is None
    q.finish(job,claimed[0]['lease_token'],{'artifact_id':'attempt-local-result'})
    assert submit(revision)==job and q.claim('next','fixture-v1') is None
    row,attempts=state(job)
    assert row['status']=='succeeded' and len(attempts)==1
    assert json.loads(row['result_summary_json'])=={'artifact_id':'attempt-local-result'}
    assert state(direct['job_id'])[0]['status']=='succeeded'
    with pytest.raises(q.QueueConflict):q.finish(job,claimed[0]['lease_token'],{'overwritten':True})
    with pytest.raises(q.QueueConflict):q.cancel(job)


def test_expiry_boundary_fences_all_old_worker_writes(ready,monkeypatch):
    _,_,revision=ready;clock=[1000.]
    monkeypatch.setattr(q.time,'time',lambda:clock[0])
    job=submit(revision);old=q.claim('old','fixture-v1',10)
    clock[0]=1010.
    for action in [lambda:q.finish(job,old['lease_token'],{}),lambda:q.fail(job,old['lease_token'],'error',True),lambda:q.renew(job,old['lease_token'],10)]:
        with pytest.raises(q.QueueConflict):action()
    new=q.claim('new','fixture-v1',10)
    assert new['attempt']==2 and new['lease_token']!=old['lease_token']
    with pytest.raises(q.QueueConflict):q.finish(job,old['lease_token'],{'stale':True})
    q.finish(job,new['lease_token'],{'fresh':True})
    row,attempts=state(job)
    assert [a['status'] for a in attempts]==['expired','succeeded']
    assert json.loads(row['result_summary_json'])=={'fresh':True}


def test_heartbeat_extends_lease_and_recovery_is_idempotent(ready,monkeypatch):
    _,_,revision=ready;clock=[1000.]
    monkeypatch.setattr(q.time,'time',lambda:clock[0])
    job=submit(revision,max_attempts=1);claim=q.claim('one','fixture-v1',10)
    clock[0]=1009.;assert q.renew(job,claim['lease_token'],10)==1019.
    clock[0]=1010.;assert q.recover_expired()==0
    clock[0]=1019.;assert q.recover_expired()==1
    assert q.recover_expired()==0 and q.claim('two','fixture-v1') is None
    row,attempts=state(job)
    assert row['status']=='failed' and row['completed_at'] and attempts[0]['status']=='expired'


@pytest.mark.parametrize('running',[False,True])
def test_cancellation_is_idempotent_and_blocks_publication(ready,running):
    _,_,revision=ready;job=submit(revision)
    claim=q.claim('one','fixture-v1') if running else None
    q.cancel(job);q.cancel(job)
    assert state(job)[0]['status']=='cancelled' and q.claim('two','fixture-v1') is None
    if claim:
        assert state(job)[1][0]['status']=='cancelled'
        with pytest.raises(q.QueueConflict):q.finish(job,claim['lease_token'],{})
    with connect() as db:assert db.execute("SELECT count(*) FROM audit_events WHERE object_id=? AND action='cancelled'",(job,)).fetchone()[0]==1


def test_retry_budget_and_attempt_errors_survive_initialise(ready):
    _,_,revision=ready;job=submit(revision,max_attempts=2)
    one=q.claim('one','fixture-v1');q.fail(job,one['lease_token'],'transient IO',True)
    initialise();initialise()
    two=q.claim('two','fixture-v1');assert two['attempt']==2
    q.fail(job,two['lease_token'],'IO again',True)
    row,attempts=state(job)
    assert row['status']=='failed' and q.claim('three','fixture-v1') is None
    assert [a['error_detail'] for a in attempts]==['transient IO','IO again']
    assert all(a['completed_at'] for a in attempts)
    other=submit(revision,'permanent')
    claim=q.claim('one','fixture-v1');q.fail(other,claim['lease_token'],'unsupported input')
    assert state(other)[0]['status']=='failed'


def test_new_worker_process_recovers_persisted_expired_attempt(ready,monkeypatch):
    _,_,revision=ready
    # Set only this process's clock in the past; the child reads the same persisted DB.
    with monkeypatch.context() as m:
        m.setattr(q.time,'time',lambda:1000.)
        job=submit(revision);old=q.claim('crashed-process','fixture-v1',1)
    result=subprocess.run([sys.executable,'-c',
        "import json; from localmolbio.database import initialise; from localmolbio.job_queue import claim; initialise(); print(json.dumps(claim('replacement-process','fixture-v1')))"],capture_output=True,text=True,check=True)
    new=json.loads(result.stdout)
    assert new['job_id']==job and new['attempt']==2
    with pytest.raises(q.QueueConflict):q.finish(job,old['lease_token'],{})
    q.finish(job,new['lease_token'],{})
    assert [a['worker_id'] for a in state(job)[1]]==['crashed-process','replacement-process']


@pytest.mark.parametrize('options',[{'max_attempts':True},{'max_attempts':0},{'max_attempts':11},{'parameters':[]},{'parameters':{'bad':float('nan')}}])
def test_invalid_submission_leaves_no_partial_job(ready,options):
    _,_,revision=ready
    with pytest.raises(ValueError):submit(revision,**options)
    with connect() as db:assert db.execute('SELECT count(*) FROM queued_jobs').fetchone()[0]==0


def test_invalid_result_does_not_finish_attempt(ready):
    _,_,revision=ready;job=submit(revision);claim=q.claim('worker','fixture-v1')
    with pytest.raises(ValueError):q.finish(job,claim['lease_token'],{'bad':float('inf')})
    with pytest.raises(ValueError):q.renew(job,claim['lease_token'],False)
    assert state(job)[0]['status']=='running'
    q.finish(job,claim['lease_token'],{})
