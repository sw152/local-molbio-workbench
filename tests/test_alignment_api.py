"""Public alignment evidence is scoped, bounded, provenance-linked and lease-free."""
import json
import pytest
from localmolbio import job_queue as q
from localmolbio.fastq_alignment import ADAPTER
from localmolbio.database import connect
from localmolbio.config import data_dir
from test_sanger_history import ready
from test_sanger_acceptance import sibling_revision
from test_alignment_evidence import minimap2, rc
from test_fastq_alignment import setup, register, run


def url(revision):return f'/api/revisions/{revision}/alignments'
def post(client,revision,ids,**extra):
    return client.post(url(revision),json={'input_ids':ids,'data_type':'ont-high-accuracy',
                                         'idempotency_key':'key',**extra})


def test_real_execution_public_provenance_and_read_pagination(setup):
    client,revision,ref=setup
    a=register(revision,[ref[30:330],rc(ref[200:550])]);b=register(revision,['A'*300],True)
    response=post(client,revision,[a['id'],b['id']]);assert response.status_code==201,response.text
    job=response.json()['id'];endpoint=url(revision)+'/'+job
    assert post(client,revision,[b['id'],a['id']]).json()['id']==job
    assert client.get(endpoint+'/reads').status_code==409
    claim=q.claim('/private/worker-identity',ADAPTER)
    running=client.get(endpoint);assert running.json()['result'] is None
    assert claim['lease_token'] not in running.text and '/private/' not in running.text
    q.fail(job,claim['lease_token'],'retry',True)
    assert run()['status']=='succeeded'
    detail=client.get(endpoint).json()
    assert detail['result']['read_count']==3 and detail['result']['base_quality_used'] is False
    assert detail['result']['analysis_performed'] and not detail['result']['whole_reference_verified']
    assert [a['status'] for a in detail['attempts']['items']]==['succeeded','failed']
    assert detail['input_labels']=={a['id']:'reads.fq',b['id']:'reads.fq.gz'}
    response=client.get(endpoint+'/reads?limit=2');assert response.status_code==200,response.text
    first=response.json();last=client.get(endpoint+'/reads?limit=2&offset=2').json()
    assert first['has_more'] and not last['has_more'] and first['total']==3
    evidence=first['items']+last['items'];assert [r['read_index'] for r in evidence]==[0,1,2]
    keyed={(r['source']['input_id'],r['source']['record_ordinal']):r for r in evidence}
    forward=keyed[(a['id'],1)]['alignments'][0];reverse=keyed[(a['id'],2)]['alignments'][0]
    assert forward['strand']=='+' and forward['reference_intervals']==[[30,330]]
    assert reverse['strand']=='-' and reverse['reference_intervals']==[[200,550]]
    assert keyed[(b['id'],1)]['status']=='no_alignment_reported'
    assert keyed[(b['id'],1)]['requires_review']
    for r in evidence:
        assert r['uniqueness_assessed'] is False and r['source']['query_name']==f'q{r["read_index"]}'
    assert client.get(endpoint+'/reads?offset=99').json()['items']==[]
    listing=client.get(url(revision));assert listing.json()['total']==1
    payload=json.dumps([detail,first,last,listing.json()])
    for hidden in ['lease_token','artifact_directory','relative_path','snapshots','worker_id','storage_path',
                   'idempotency_key',str(data_dir()),claim['lease_token']]:assert hidden not in payload
    assert 'reads' not in detail['result'] and 'sources' not in detail['result']
    assert 'result' not in listing.json()['items'][0]
    assert client.post(endpoint+'/cancel').status_code==409


def test_revision_adapter_and_submission_validation(setup):
    client,revision,ref=setup;item=register(revision,[ref[40:400]])
    job=post(client,revision,[item['id']]).json()['id'];other=sibling_revision(revision)
    assert client.get(url(other)).json()['total']==0
    for foreign_rev,foreign_job in [(other,job),(revision,q.enqueue(revision,'other','unknown',{}))]:
        endpoint=url(foreign_rev)+'/'+foreign_job
        assert client.get(endpoint).status_code==404
        assert client.get(endpoint+'/reads').status_code==404
        assert client.post(endpoint+'/cancel').status_code==404
    assert post(client,other,[item['id']]).status_code==422
    assert post(client,'missing',[item['id']]).status_code==404
    assert client.get(url('missing')).status_code==404
    for extra in [{'adapter':'shell'},{'parameters':{}},{'binary':'/bin/sh'}, {'data_type':'automatic'},
                  {'input_ids':[]},{'input_ids':[item['id'],item['id']]},{'max_attempts':True},
                  {'max_attempts':0},{'max_attempts':11},{'idempotency_key':' '}]:
        response=post(client,revision,[item['id']],**extra);assert response.status_code==422,response.text
    assert post(client,revision,[item['id']],data_type='short-single').status_code==409
    for suffix in ['', '/'+job, '/'+job+'/reads']:
        for query in ['limit=0','limit=101','offset=-1']:
            assert client.get(url(revision)+suffix+'?'+query).status_code==422
    assert client.get(url(revision)).json()['total']==1


def test_job_and_attempt_pages_and_failed_result(setup):
    client,revision,ref=setup;item=register(revision,[ref[40:400]])
    jobs=[post(client,revision,[item['id']],idempotency_key=str(i),max_attempts=10).json()['id'] for i in range(7)]
    first=client.get(url(revision)).json();last=client.get(url(revision)+'?offset=5').json()
    assert first['has_more'] and not last['has_more']
    assert {r['id'] for r in first['items']+last['items']}==set(jobs)
    for job in jobs[1:]:q.cancel(job)
    for _ in range(7):
        claim=q.claim('worker',ADAPTER);q.fail(claim['job_id'],claim['lease_token'],'retry',True)
    endpoint=url(revision)+'/'+jobs[0]
    first=client.get(endpoint).json();last=client.get(endpoint+'?offset=5').json()
    assert first['result'] is None
    assert [a['number'] for a in first['attempts']['items']+last['attempts']['items']]==list(range(7,0,-1))
    assert first['attempts']['has_more'] and not last['attempts']['has_more']
    claim=q.claim('worker',ADAPTER);q.fail(claim['job_id'],claim['lease_token'],'input_hash_mismatch',False)
    assert client.get(endpoint).json()['status']=='failed'
    assert client.get(endpoint+'/reads').status_code==409


@pytest.mark.parametrize('running',[False,True])
def test_cancel_fences_late_result_and_never_serves_diagnostics(setup,running):
    client,revision,ref=setup;item=register(revision,[ref[40:400]])
    job=post(client,revision,[item['id']]).json()['id'];endpoint=url(revision)+'/'+job
    claim=q.claim('worker',ADAPTER) if running else None
    for _ in range(2):assert client.post(endpoint+'/cancel').json()['status']=='cancelled'
    if claim:
        with pytest.raises(q.QueueConflict):q.finish(job,claim['lease_token'],{'diagnostic':'late'})
    assert client.get(endpoint).json()['result'] is None
    assert client.get(endpoint+'/reads').status_code==409
    assert run()=={'status':'idle'}


@pytest.mark.parametrize('change',['source','schema','hashes','ordinal','length','sequence_hash'])
def test_corrupt_evidence_fails_closed(setup,change):
    client,revision,ref=setup;item=register(revision,[ref[40:400]])
    job=post(client,revision,[item['id']]).json()['id'];assert run()['status']=='succeeded'
    with connect() as db:
        result=json.loads(db.execute('SELECT result_summary_json FROM analysis_jobs WHERE id=?',(job,)).fetchone()[0])
        if change=='source':result['sources'][0]['query_name']='q99'
        elif change=='schema':result['schema_version']=999
        elif change=='hashes':result['query_sha256']=[]
        elif change=='ordinal':result['sources'][0]['record_ordinal']=0
        elif change=='length':result['sources']=[]
        else:result['sources'][0]['sequence_sha256']='0'*64
        db.execute('UPDATE analysis_jobs SET result_summary_json=? WHERE id=?',(json.dumps(result),job))
    assert client.get(url(revision)+'/'+job+'/reads').status_code==409
