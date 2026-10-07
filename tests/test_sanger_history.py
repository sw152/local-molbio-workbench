"""Append-only reports, additive legacy migration and concurrent rerun safety."""
import io
import json
import random
import zipfile
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from threading import Barrier

import pytest
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from fastapi.testclient import TestClient
from abif_fixture import synthetic_ab1
from localmolbio.api import app
from localmolbio.database import connect, initialise
from localmolbio.sanger_verification import align_sanger_read


@pytest.fixture
def ready(tmp_path, monkeypatch):
    monkeypatch.setenv('MOLBIO_DATA_DIR', str(tmp_path/'runtime'))
    rng=random.Random(617)
    reference=''.join(rng.choice('ACGT') for _ in range(600))
    record=SeqRecord(Seq(reference),id='history',name='history',annotations={'molecule_type':'DNA','topology':'linear'})
    stream=io.StringIO();SeqIO.write(record,stream,'genbank')
    archive=io.BytesIO()
    with zipfile.ZipFile(archive,'w') as z:z.writestr('history.gb',stream.getvalue())
    with TestClient(app) as client:
        assert client.post('/api/imports/benchling',files={'file':('fixture.zip',archive.getvalue())}).status_code==201
        sequence=client.get('/api/sequences').json()[0]['id']
        revision=client.get(f'/api/sequences/{sequence}/revisions').json()[0]['id']
        read=client.post('/api/sanger-reads',data={'sequence_revision_id':revision},files={'file':('trace.ab1',synthetic_ab1(reference[50:250],[38]*200,trace=True))}).json()['id']
        yield client,read,revision


def test_rerun_preserves_history_inputs_and_job_provenance(ready):
    client,read,revision=ready
    first=client.post('/api/sanger-verifications',json={'sequencing_read_id':read}).json()
    snapshot=client.get(f'/api/sanger-reads/{read}/analyses').json()['items'][0]
    second=client.post('/api/sanger-verifications',json={'sequencing_read_id':read,'previous_alignment_id':first['id'],'direction':'forward'})
    assert second.status_code==201,second.text
    assert second.json()['run_number']==2
    history=client.get(f'/api/sanger-reads/{read}/analyses?limit=1').json()
    assert history['total']==2 and history['has_more']
    assert history['items'][0]['previous_alignment_id']==first['id']
    assert client.get(f'/api/sanger-reads/{read}/analyses?offset=1').json()['items']==[snapshot]
    assert client.get(f'/api/sanger-reads/{read}/analyses?offset=2').json()['items']==[]
    listing=client.get(f'/api/sequence-revisions/{revision}/sanger-reads').json()
    assert len(listing)==1 and listing[0]['alignment_id']==second.json()['id']
    assert listing[0]['direction']=='unknown'
    assert listing[0]['evidence']['requested_direction']=='forward'
    assert client.post('/api/sanger-verifications',json={'sequencing_read_id':read,'previous_alignment_id':first['id']}).status_code==409
    assert client.post('/api/sanger-verifications',json={'sequencing_read_id':read,'direction':'sideways'}).status_code==422
    assert client.get('/api/sanger-reads/missing/analyses').status_code==404
    assert client.get(f'/api/sanger-reads/{read}/analyses?limit=0').status_code==422
    with connect() as db:
        job=db.execute('SELECT * FROM analysis_jobs WHERE id=?',(second.json()['job_id'],)).fetchone()
        assert json.loads(job['parameters_json'])['quality_source']=='hash_checked_original_ab1'
        assert json.loads(job['input_manifest_json'])['previous_alignment_id']==first['id']
        assert len(db.execute("SELECT * FROM audit_events WHERE action='analysis-created'").fetchall())==2


def test_legacy_migration_is_additive_idempotent_and_quality_recovery_is_per_run(ready):
    client,read,revision=ready
    first=client.post('/api/sanger-verifications',json={'sequencing_read_id':read}).json()
    with connect() as db:
        run=dict(db.execute('SELECT * FROM sanger_analysis_runs').fetchone())
        report=json.loads(run['report_json'])
        # Simulate the original schema's one-report row with absent detailed evidence/quality.
        values={key:report[key] for key in ('reference_start','reference_end','wraps_origin','aligned_bases','matched_bases','mismatched_bases','inserted_bases','deleted_bases','identity_fraction')}
        values.update(id=run['id'],sequencing_read_id=read,sequence_revision_id=revision,analysis_job_id=run['analysis_job_id'],variants_json=json.dumps(report['variants']),evidence_json='{}',created_at=run['created_at'])
        db.execute(f"INSERT INTO sanger_alignments ({','.join(values)}) VALUES ({','.join('?' for _ in values)})",tuple(values.values()))
        db.execute('DELETE FROM sanger_analysis_runs')
        db.execute("UPDATE sequencing_reads SET qualities_json='[]'")
        legacy=dict(db.execute('SELECT * FROM sanger_alignments').fetchone())
    initialise();initialise()
    history=client.get(f'/api/sanger-reads/{read}/analyses').json()
    assert history['total']==1 and history['items'][0]['alignment_id']==first['id']
    assert history['items'][0]['evidence']['review_flags']==['legacy_evidence_unavailable']
    second=client.post('/api/sanger-verifications',json={'sequencing_read_id':read,'previous_alignment_id':first['id']})
    assert second.status_code==201,second.text
    assert 'quality_unavailable' not in second.json()['evidence']['review_flags']
    initialise()
    with connect() as db:
        assert dict(db.execute('SELECT * FROM sanger_alignments').fetchone())==legacy
        assert db.execute('SELECT qualities_json FROM sequencing_reads').fetchone()[0]=='[]'
        assert db.execute('SELECT count(*) FROM sanger_analysis_runs').fetchone()[0]==2


@pytest.mark.parametrize('failure',['changed','missing','calls','qualities','engine'])
def test_failed_rerun_does_not_alter_history(ready,monkeypatch,failure):
    client,read,revision=ready
    first=client.post('/api/sanger-verifications',json={'sequencing_read_id':read}).json()
    before=client.get(f'/api/sanger-reads/{read}/analyses').json()
    with connect() as db:
        stored=Path(db.execute('SELECT storage_path FROM sequencing_reads').fetchone()[0])
        if failure=='calls':db.execute("UPDATE sequencing_reads SET base_sequence='A'")
        if failure=='qualities':db.execute("UPDATE sequencing_reads SET qualities_json='[1]'")
    if failure=='changed':stored.write_bytes(b'changed')
    if failure=='missing':stored.unlink()
    if failure=='engine':
        def broken(*args,**kwargs):raise RuntimeError('simulated failure')
        monkeypatch.setattr('localmolbio.api.align_sanger_read',broken)
    response=client.post('/api/sanger-verifications',json={'sequencing_read_id':read,'previous_alignment_id':first['id']})
    assert response.status_code==(500 if failure=='engine' else 409)
    assert client.get(f'/api/sanger-reads/{read}/analyses').json()==before
    assert len(client.get('/api/jobs').json())==1


def test_simultaneous_analysis_requests_save_exactly_one_run(ready,monkeypatch):
    client,read,revision=ready
    barrier=Barrier(2)
    def synchronize(*args,**kwargs):
        barrier.wait(timeout=5)
        return align_sanger_read(*args,**kwargs)
    monkeypatch.setattr('localmolbio.api.align_sanger_read',synchronize)
    with ThreadPoolExecutor(max_workers=2) as pool:
        futures=[pool.submit(client.post,'/api/sanger-verifications',json={'sequencing_read_id':read}) for _ in range(2)]
        assert sorted(f.result().status_code for f in futures)==[201,409]
    assert client.get(f'/api/sanger-reads/{read}/analyses').json()['total']==1
    assert len(client.get('/api/jobs').json())==1
