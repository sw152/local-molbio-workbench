"""Registered bytes to real bounded minimap2 jobs, with snapshot and lease isolation."""
from dataclasses import replace
import gzip
import io
import json
from pathlib import Path
import subprocess
import sys
from hashlib import sha256

import pytest
from localmolbio import fastq_alignment as pipeline, worker, job_queue as q
from localmolbio.database import connect
from localmolbio.config import data_dir
from localmolbio.fastq_inputs import register_fastq
from test_sanger_history import ready
from test_sanger_acceptance import sibling_revision
from test_input_worker import stored
from test_alignment_evidence import minimap2, rc
from test_fastq_worker import path


@pytest.fixture
def setup(ready,minimap2,monkeypatch):
    monkeypatch.setenv('MOLBIO_MINIMAP2',str(minimap2))
    client,read,revision=ready
    with connect() as db:reference=db.execute('SELECT sequence_text FROM sequence_revisions WHERE id=?',(revision,)).fetchone()[0]
    return client,revision,reference


def register(revision,reads,compressed=False):
    # Intentionally repeated titles: provenance must use file ID + ordinal.
    data=''.join(f'@duplicate\n{r}\n+\n'+('I'*len(r))+'\n' for r in reads).encode()
    content=gzip.compress(data,mtime=0) if compressed else data
    return register_fastq(revision,'reads.fq.gz' if compressed else 'reads.fq',io.BytesIO(content),'phred33')


def enqueue(revision,items,key='align',kind='ont-high-accuracy',**kwargs):
    return pipeline.enqueue_alignment(revision,[r['id'] for r in items],kind,key,**kwargs)


def run():return worker.run_once('alignment-worker',adapter=pipeline.ADAPTER)


def attempts():return list((data_dir()/'alignment-attempts').glob('*/*'))


def test_real_multifile_snapshots_ordinals_results_and_adapter_isolation(setup):
    client,revision,reference=setup
    a=register(revision,[reference[50:350].lower(),reference[100:450]])
    b=register(revision,[rc(reference[250:550])],True)
    job=enqueue(revision,[a,b]);assert enqueue(revision,[b,a])==job
    assert worker.run_once('default')=={'status':'idle'}
    outcome=run();assert outcome['status']=='succeeded' and outcome['analysis_performed']
    result=json.loads(stored(job)['result_summary_json'])
    assert result['scope']=='registered_fastq_local_alignments_only' and result['whole_reference_verified'] is False
    assert len(result['sources'])==3 and len(result['snapshots'])==2
    expected={(a['id'],1):reference[50:350],(a['id'],2):reference[100:450],(b['id'],1):rc(reference[250:550])}
    root=data_dir()/result['artifact_directory']
    assert root in attempts()
    assert json.loads((root/'attempt.json').read_text())['attempt_number']==result['attempt_number']==1
    for i,source in enumerate(result['sources']):
        assert source['query_name']==f'q{i}'
        assert source['sequence_sha256']==sha256(expected[source['input_id'],source['record_ordinal']].encode()).hexdigest()
        assert result['reads'][i]['alignments']
    for snapshot in result['snapshots']:
        copy=root/snapshot['relative_path']
        assert copy.read_bytes()==path(snapshot['input_id']).read_bytes()
        assert copy.stat().st_ino != path(snapshot['input_id']).stat().st_ino
    assert sha256((root/'sources.json').read_bytes()).hexdigest()==result['sources_sha256']
    assert json.loads((root/'attempt-result.json').read_text())==result
    # Stored jobs never contain an absolute source storage path or a lease token.
    assert str(data_dir()) not in stored(job)['result_summary_json']
    assert stored(job)['attempts'][0]['id'] not in stored(job)['result_summary_json']
    assert len(client.get(f'/api/revisions/{revision}/fastq-inputs').json()['items'])==2


@pytest.mark.parametrize('change,code',[
    ('bytes','input_hash_mismatch'),('metadata','fastq_identity_changed_since_submission'),
    ('reference','reference_hash_mismatch'),('summary','registered_fastq_summary_does_not_match_bytes'),
    ('iupac','reads_must_be_uppercase_acgtn_within_limit'),
])
def test_invalid_input_never_produces_published_alignment(setup,change,code):
    _,revision,reference=setup;sequence=reference[50:400]
    if change=='iupac':sequence=sequence[:30]+'R'+sequence[31:]
    item=register(revision,[sequence])
    if change=='summary':
        with connect() as db:
            summary=item['summary'].copy();summary['q30_bases']-=1
            db.execute('UPDATE fastq_inputs SET summary_json=? WHERE id=?',(json.dumps(summary),item['id']))
    job=enqueue(revision,[item])
    if change=='bytes':path(item['id']).write_bytes(path(item['id']).read_bytes().replace(b'I',b'J'))
    elif change=='metadata':
        with connect() as db:db.execute('UPDATE fastq_inputs SET size_bytes=size_bytes+1 WHERE id=?',(item['id'],))
    elif change=='reference':
        with connect() as db:db.execute("UPDATE sequence_revisions SET sequence_text='ACGT' WHERE id=?",(revision,))
    result=run();assert result['error']==code
    assert stored(job)['status']=='failed' and stored(job)['result_summary_json']=='{}'


@pytest.mark.parametrize('limit,code',[
    ('decoded','decoded_size_limit_exceeded'),('records','record_count_limit_exceeded'),
    ('aggregate','alignment_aggregate_input_limit_exceeded'),('original','alignment_total_original_bytes_limit_exceeded'),
])
def test_limits_reject_whole_input_instead_of_prefix(setup,monkeypatch,limit,code):
    _,revision,reference=setup;item=register(revision,[reference[50:300],reference[100:350]])
    if limit=='original':
        monkeypatch.setattr(pipeline,'SNAPSHOT_BYTES',20)
        with pytest.raises(pipeline.InputValidationError,match=code):enqueue(revision,[item])
        assert not attempts();return
    job=enqueue(revision,[item])
    if limit=='decoded':monkeypatch.setattr(pipeline,'LIMITS',replace(pipeline.LIMITS,decoded_bytes=20))
    elif limit=='records':monkeypatch.setattr(pipeline,'LIMITS',replace(pipeline.LIMITS,records=1))
    else:monkeypatch.setattr(pipeline.core,'MAX_BASES',300)
    assert run()['error']==code and stored(job)['result_summary_json']=='{}'
    assert all(not (p/'alignment').exists() for p in attempts())


@pytest.mark.parametrize('change,code',[
    ('original','input_hash_mismatch'),('snapshot','alignment_snapshot_changed_during_run'),
    ('mapping','alignment_source_mapping_changed'),('reference','reference_hash_mismatch'),
])
def test_mutation_during_actual_alignment_prevents_publication(setup,monkeypatch,change,code):
    _,revision,reference=setup;item=register(revision,[reference[50:400]])
    job=enqueue(revision,[item]);real=pipeline.core.align
    def mutate(*args,**kwargs):
        result=real(*args,**kwargs);attempt=kwargs['output_dir'].parent
        if change=='original':path(item['id']).write_bytes(path(item['id']).read_bytes().replace(b'I',b'J'))
        elif change=='snapshot':(attempt/'inputs'/'0.fastq').write_bytes(b'changed')
        elif change=='mapping':(attempt/'sources.json').write_text('[]')
        else:
            with connect() as db:db.execute("UPDATE sequence_revisions SET sequence_text='A' WHERE id=?",(revision,))
        return result
    monkeypatch.setattr(pipeline.core,'align',mutate)
    assert run()['error']==code and stored(job)['result_summary_json']=='{}'
    assert all(not (p/'attempt-result.json').exists() for p in attempts())


@pytest.mark.parametrize('mode',['cancel','expire'])
def test_late_attempt_cannot_publish_and_recovery_uses_new_directory(setup,monkeypatch,mode):
    _,revision,reference=setup;item=register(revision,[reference[50:400]])
    job=enqueue(revision,[item]);real=pipeline.execute
    def lose(claim,pulse):
        result=real(claim,pulse)
        if mode=='cancel':q.cancel(job)
        else:
            with connect() as db:db.execute('UPDATE queued_jobs SET lease_expires=0 WHERE job_id=?',(job,))
        return result
    with monkeypatch.context() as m:
        m.setattr(pipeline,'execute',lose)
        assert run()['status']=='lease_lost'
    first=attempts()[0];assert (first/'attempt-result.json').exists()
    assert stored(job)['result_summary_json']=='{}'
    if mode=='expire':
        q.recover_expired();assert run()['status']=='succeeded'
        assert len(attempts())==2
        result=json.loads(stored(job)['result_summary_json'])
        assert data_dir()/result['artifact_directory']!=first
        assert result['attempt_number']==2
        assert [a['status'] for a in stored(job)['attempts']]==['expired','succeeded']
    else:assert stored(job)['status']=='cancelled'


def test_retry_preserves_independent_failed_attempt(setup,monkeypatch):
    _,revision,reference=setup;item=register(revision,[reference[50:400]])
    job=enqueue(revision,[item]);real=pipeline.core.align
    def io_failure(*args,**kwargs):raise OSError('private path')
    with monkeypatch.context() as m:
        m.setattr(pipeline.core,'align',io_failure)
        assert run()['error']=='alignment_artifact_io_error'
    assert stored(job)['status']=='queued';first=attempts()[0]
    assert run()['status']=='succeeded' and len(attempts())==2
    assert not (first/'attempt-result.json').exists()
    assert [a['status'] for a in stored(job)['attempts']]==['failed','succeeded']
    assert 'private path' not in json.dumps(stored(job))


def test_submission_constraints_and_tool_identity(setup,monkeypatch):
    _,revision,reference=setup;item=register(revision,[reference[50:400]])
    with pytest.raises(pipeline.InputValidationError):enqueue(sibling_revision(revision),[item])
    with pytest.raises(pipeline.InputValidationError):enqueue(revision,[item],kind='automatic')
    with pytest.raises(pipeline.InputValidationError):enqueue(revision,[item,item])
    job=enqueue(revision,[item])
    with pytest.raises(q.QueueConflict):enqueue(revision,[item],kind='short-single')
    original=pipeline._tool
    def changed():
        binary,info=original();return binary,{**info,'binary_sha256':'0'*64}
    monkeypatch.setattr(pipeline,'_tool',changed)
    assert run()['error']=='aligner_changed_since_submission' and stored(job)['result_summary_json']=='{}'
    assert not attempts()


def test_fresh_cli_enqueue_and_alignment_processes(setup):
    _,revision,reference=setup;item=register(revision,[reference[50:400]],True)
    def cli(*args):return subprocess.run([sys.executable,'-m','localmolbio.worker',*args],capture_output=True,text=True)
    submission=cli('enqueue-alignment','--revision',revision,'--input-id',item['id'],'--data-type','ont-high-accuracy','--key','cli')
    assert submission.returncode==0,submission.stderr
    job=json.loads(submission.stdout)['job_id']
    completed=cli('once','--adapter','alignment')
    assert completed.returncode==0,completed.stderr
    assert json.loads(completed.stdout)['analysis_performed'] is True
    assert json.loads(stored(job)['result_summary_json'])['reads'][0]['alignments']


def test_snapshot_parser_cancellation_stops_before_aligner(setup,monkeypatch):
    _,revision,reference=setup;item=register(revision,[reference[50:400]])
    job=enqueue(revision,[item]);real=pipeline.inspect_fastq
    def cancel_scan(*args,**kwargs):
        q.cancel(job)
        return real(*args,**kwargs)
    monkeypatch.setattr(pipeline,'inspect_fastq',cancel_scan)
    assert run()['status']=='lease_lost' and stored(job)['status']=='cancelled'
    assert stored(job)['result_summary_json']=='{}'
    assert all(not (p/'alignment').exists() for p in attempts())


def test_registered_circular_read_preserves_wrapped_coordinates(setup):
    _,revision,reference=setup
    with connect() as db:db.execute("UPDATE sequence_revisions SET topology='circular' WHERE id=?",(revision,))
    item=register(revision,[reference[-200:]+reference[:200]],True)
    job=enqueue(revision,[item]);assert run()['status']=='succeeded'
    result=json.loads(stored(job)['result_summary_json'])
    hit=next(a for a in result['reads'][0]['alignments'] if a['exact_matches']==400)
    assert hit['reference_intervals']==[[len(reference)-200,len(reference)],[0,200]]
    assert not result['whole_reference_verified']


def test_unknown_errors_are_redacted_and_configuration_is_required(setup,monkeypatch):
    _,revision,reference=setup;item=register(revision,[reference[50:400]])
    with monkeypatch.context() as m:
        m.delenv('MOLBIO_MINIMAP2')
        with pytest.raises(pipeline.InputValidationError,match='configuration_missing'):enqueue(revision,[item])
    job=enqueue(revision,[item])
    def unexpected(*args):raise RuntimeError('private input and path')
    monkeypatch.setattr(pipeline,'execute',unexpected)
    assert run()['error']=='unexpected_alignment_error'
    assert stored(job)['status']=='failed' and 'private' not in json.dumps(stored(job))
