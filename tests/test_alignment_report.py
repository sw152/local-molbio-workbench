"""Portable exports include all records, preserve boundaries, and reject corrupt identity."""
import json
import pytest
from localmolbio.database import connect
from localmolbio.config import data_dir
from localmolbio import job_queue as q
from test_alignment_evidence import minimap2, rc
from test_sanger_history import ready
from test_sanger_acceptance import sibling_revision
from test_fastq_alignment import setup,register,run


def job_for(setup):
    client,revision,ref=setup
    item=register(revision,[ref[30:400],rc(ref[100:500]),ref[50:200]+ref[207:500],'A'*300])
    root=f'/api/revisions/{revision}/alignments'
    job=client.post(root,json={'input_ids':[item['id']],'data_type':'ont-high-accuracy','idempotency_key':'report'}).json()['id']
    return client,revision,ref,item,job,root+'/'+job+'/report'


def test_complete_report_real_alignment_provenance_safe_html_and_stored_snapshot(setup):
    client,revision,ref,item,job,url=job_for(setup)
    assert client.get(url).status_code==409
    assert run()['status']=='succeeded'
    label='<img src="https://invalid.test/secret" onerror="alert(1)"> & sample'
    with connect() as db:db.execute('UPDATE fastq_inputs SET original_filename=? WHERE id=?',(label,item['id']))
    r=client.get(url);assert r.status_code==200,r.text
    report=r.json();assert report['schema']=='localmolbio.alignment-review-report'
    assert report['selection']=='all_task_records' and len(report['reads'])==4
    assert [x['source']['record_ordinal'] for x in report['reads']]==[1,2,3,4]
    assert report['reads'][1]['alignments'][0]['strand']=='-'
    assert report['reads'][2]['alignments'][0]['deleted_bases']==7
    assert report['reads'][3]['alignments']==[]
    assert report['coverage']['deletion_evidence_bases']==7 and report['review_regions']['deletion']==[[200,207]]
    assert report['review_regions']['unpaired']==[[0,30],[500,600]]
    assert report['inputs'][0]['label_at_export']==label
    assert report['reference']['id']==revision and report['job']['id']==job
    assert report['source_verification_at_export']=='not_performed'
    for flag in ['base_quality_used','pairing_used','consensus_performed','whole_reference_verified','candidate_search_exhaustive']:
        assert report['evidence'][flag] is False
    html=client.get(url+'?format=html');assert html.status_code==200,html.text
    assert '<img ' not in html.text and '&lt;img ' in html.text and '<script' not in html.text
    for phrase in ['No whole-plasmid verdict','Record 4','[200, 207)','← Reverse','all task records','Original files were not rechecked']:
        assert phrase.lower() in html.text.lower(),phrase
    for response in [r,html]:
        assert 'attachment;' in response.headers['content-disposition']
        assert 'sample' not in response.headers['content-disposition']
        assert response.headers['cache-control']=='no-store'
        assert response.headers['x-content-type-options']=='nosniff'
        assert "default-src 'none'" in response.headers['content-security-policy']
        for hidden in ['artifact_directory','lease_token','snapshots','storage_path',str(data_dir()),ref,'worker_id']:
            assert hidden not in response.text,hidden
    # Export reproduces the committed snapshot even if the current file changes; no reanalysis.
    with connect() as db:
        path=db.execute('SELECT storage_path FROM fastq_inputs WHERE id=?',(item['id'],)).fetchone()[0]
        count=db.execute('SELECT count(*) FROM job_attempts WHERE job_id=?',(job,)).fetchone()[0]
    (data_dir()/path).write_bytes(b'changed after publication')
    again=client.get(url).json();assert again['reads']==report['reads'] and again['coverage']==report['coverage']
    with connect() as db:assert db.execute('SELECT count(*) FROM job_attempts WHERE job_id=?',(job,)).fetchone()[0]==count


@pytest.mark.parametrize('change',['reference','tool','input','attempt','flag','coordinate','geometry','source','schema','scope'])
def test_export_rejects_corrupt_result(setup,change):
    client,revision,ref,item,job,url=job_for(setup);assert run()['status']=='succeeded'
    with connect() as db:
        e=json.loads(db.execute('SELECT result_summary_json FROM analysis_jobs WHERE id=?',(job,)).fetchone()[0])
        if change=='reference':e['reference_sha256']='0'*64
        elif change=='tool':e['tool']['binary_sha256']='0'*64
        elif change=='input':e['input_identities'][0]['sha256']='0'*64
        elif change=='attempt':e['attempt_number']+=1
        elif change=='flag':e['whole_reference_verified']=True
        elif change=='coordinate':e['coordinate_convention']='one_based'
        elif change=='geometry':e['reads'][0]['alignments'][0]['cigar']='10M'
        elif change=='source':e['sources'][0]['record_ordinal']=99
        elif change=='scope':e['scope']='whole_plasmid_verified'
        else:e['schema_version']=99
        db.execute('UPDATE analysis_jobs SET result_summary_json=? WHERE id=?',(json.dumps(e),job))
    for fmt in ['html','json']:assert client.get(url+'?format='+fmt).status_code==409


def test_export_revision_adapter_cancel_and_format_boundaries(setup):
    client,revision,ref,item,job,url=job_for(setup)
    other=sibling_revision(revision);assert client.get(url.replace(revision,other)).status_code==404
    foreign=q.enqueue(revision,'other','unknown',{});assert client.get(url.replace(job,foreign)).status_code==404
    assert client.get(url+'?format=pdf').status_code==422
    q.cancel(job)
    for fmt in ['html','json']:assert client.get(url+'?format='+fmt).status_code==409
