"""Public job views are revision-scoped and never reveal worker lease tokens."""
import pytest
from localmolbio import job_queue as q
from localmolbio.input_validation import ADAPTER
from localmolbio.worker import run_once
from localmolbio.database import connect
from test_sanger_history import ready
from test_sanger_acceptance import sibling_revision


def url(revision):return f'/api/revisions/{revision}/input-checks'
def post(client,revision,read,key='key',**extra):return client.post(url(revision),json={'read_ids':[read],'idempotency_key':key,**extra})


def test_submit_details_execution_and_credential_redaction(ready):
    client,read,revision=ready
    queued=post(client,revision,read);assert queued.status_code==201,queued.text
    job=queued.json()['id'];assert queued.json()['status']=='queued'
    assert post(client,revision,read).json()['id']==job
    claim=q.claim('worker',ADAPTER)
    detail=client.get(url(revision)+'/'+job)
    assert detail.status_code==200 and detail.json()['status']=='running'
    assert detail.json()['read_labels']=={read:'trace.ab1'}
    assert claim['lease_token'] not in detail.text
    assert 'lease_token' not in detail.text and 'storage_path' not in detail.text and 'idempotency_key' not in detail.text
    q.fail(job,claim['lease_token'],'retry',True)
    assert run_once('replacement')['status']=='succeeded'
    detail=client.get(url(revision)+'/'+job).json()
    assert detail['result']['analysis_performed'] is False
    assert [a['status'] for a in detail['attempts']['items']]==['succeeded','failed']
    assert client.post(url(revision)+'/'+job+'/cancel').status_code==409


def test_revision_and_adapter_isolation(ready):
    client,read,revision=ready;job=post(client,revision,read).json()['id'];other=sibling_revision(revision)
    assert client.get(url(other)).json()['items']==[]
    assert client.get(url(other)+'/'+job).status_code==404
    assert client.post(url(other)+'/'+job+'/cancel').status_code==404
    assert post(client,other,read,key='foreign').status_code==422
    unknown=q.enqueue(revision,'other-adapter','unregistered',{})
    assert client.get(url(revision)+'/'+unknown).status_code==404
    assert client.post(url(revision)+'/'+unknown+'/cancel').status_code==404
    assert client.get(url(revision)).json()['total']==1
    assert client.get(url('missing')).status_code==404
    assert post(client,'missing',read).status_code==404


def test_job_and_attempt_pagination(ready):
    client,read,revision=ready
    jobs=[post(client,revision,read,key=str(i),max_attempts=10).json()['id'] for i in range(7)]
    first=client.get(url(revision)+'?limit=5').json();second=client.get(url(revision)+'?limit=5&offset=5').json()
    assert first['has_more'] and not second['has_more'] and first['total']==7
    assert {r['id'] for r in first['items']+second['items']}==set(jobs)
    # Disable the other requests so each retry addresses the same job.
    for job in jobs[1:]:q.cancel(job)
    for i in range(7):
        claim=q.claim('worker',ADAPTER);q.fail(claim['job_id'],claim['lease_token'],'retry',True)
    a=client.get(url(revision)+'/'+jobs[0]+'?limit=5').json()['attempts']
    b=client.get(url(revision)+'/'+jobs[0]+'?limit=5&offset=5').json()['attempts']
    assert a['has_more'] and not b['has_more']
    assert [r['number'] for r in a['items']+b['items']]==list(range(7,0,-1))


@pytest.mark.parametrize('running',[False,True])
def test_cancellation_revokes_publishing_and_is_idempotent(ready,running):
    client,read,revision=ready;job=post(client,revision,read).json()['id']
    claim=q.claim('worker',ADAPTER) if running else None
    for _ in range(2):assert client.post(url(revision)+'/'+job+'/cancel').json()['status']=='cancelled'
    if claim:
        with pytest.raises(q.QueueConflict):q.finish(job,claim['lease_token'],{})
    assert run_once('next')=={'status':'idle'}


def test_submission_validation_and_conflicting_key(ready):
    client,read,revision=ready
    for extra in [{'adapter':'arbitrary'},{'parameters':{}},{'max_attempts':True},{'max_attempts':0},{'read_ids':[]},{'read_ids':[read,read]},{'idempotency_key':' '}]:
        response=client.post(url(revision),json={'read_ids':[read],'idempotency_key':'bad',**extra})
        assert response.status_code==422,response.text
    assert client.get(url(revision)).json()['total']==0
    job=post(client,revision,read).json()['id']
    assert post(client,revision,read,max_attempts=2).status_code==409
    for query in ['limit=0','offset=-1','limit=101']:
        assert client.get(url(revision)+'?'+query).status_code==422
        assert client.get(url(revision)+'/'+job+'?'+query).status_code==422
