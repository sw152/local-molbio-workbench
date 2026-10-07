"""Group downloads must bind the full displayed view, never silently adopt new runs."""
from copy import deepcopy
from hashlib import sha256
import json
from pathlib import Path
from localmolbio.database import connect
from localmolbio.sanger_group_report import snapshot_digest
from test_sanger_history import ready


def view(client,revision):
    return client.get(f'/api/sequence-revisions/{revision}/sanger-reads?include_summary=true').json()

def export(client,revision,digest,format='json'):
    return client.get(f'/api/sequence-revisions/{revision}/sanger-report',params={'expected_snapshot':digest,'format':format})

def test_full_view_digest_and_export_are_stable_and_source_complete(ready):
    client,read,revision=ready
    client.post('/api/sanger-verifications',json={'sequencing_read_id':read})
    current=view(client,revision);digest=current['snapshot_sha256']
    expected=sha256(json.dumps({'reads':current['reads'],'summary':current['summary']},sort_keys=True,separators=(',',':'),ensure_ascii=False).encode()).hexdigest()
    assert digest==expected
    response=export(client,revision,digest)
    assert response.status_code==200,response.text
    payload=response.json()
    assert payload['snapshot']==current and payload['snapshot_sha256']==digest
    assert payload['schema']=='localmolbio.sanger-group-report'
    assert payload['snapshot']['reads'][0]['evidence']['base_mapping']['blocks']
    assert payload['snapshot']['summary']['source_runs'][0]['read_sha256']
    assert payload['source_verification_at_export']=='not_performed'
    assert 'storage_path' not in response.text and 'parser_metadata' not in response.text
    assert response.headers['cache-control']=='no-store' and response.headers['x-content-type-options']=='nosniff'
    assert export(client,revision,digest).json()['snapshot_sha256']==digest
    # A read-card field change must be detected even when the geometric summary does not change.
    changed=deepcopy(current);changed['reads'][0]['quality_summary']['mean_phred']=0
    assert snapshot_digest(changed)!=digest


def test_stale_snapshot_rejected_after_rerun_without_altering_previous_export(ready):
    client,read,revision=ready
    first=client.post('/api/sanger-verifications',json={'sequencing_read_id':read}).json()
    previous=view(client,revision);saved=export(client,revision,previous['snapshot_sha256']).json()
    second=client.post('/api/sanger-verifications',json={'sequencing_read_id':read,'previous_alignment_id':first['id'],'trim_quality_threshold':20}).json()
    stale=export(client,revision,previous['snapshot_sha256'])
    assert stale.status_code==409 and 'changed' in stale.json()['detail']
    current=view(client,revision)
    latest=export(client,revision,current['snapshot_sha256']).json()
    assert latest['snapshot']['reads'][0]['alignment_id']==second['id']
    assert saved['snapshot']['reads'][0]['alignment_id']==first['id']
    assert saved['snapshot_sha256']!=latest['snapshot_sha256']


def test_unanalyzed_and_legacy_exclusions_are_exported_not_hidden(ready):
    client,read,revision=ready
    pending=view(client,revision)
    data=export(client,revision,pending['snapshot_sha256']).json()['snapshot']
    assert data['summary']['source_runs'][0]['exclusion_reasons']==['not_analyzed']
    assert data['reads'][0]['alignment_id'] is None
    first=client.post('/api/sanger-verifications',json={'sequencing_read_id':read}).json()
    with connect() as db:
        row=db.execute('SELECT report_json FROM sanger_analysis_runs WHERE id=?',(first['id'],)).fetchone()
        report=json.loads(row[0]);report['evidence'].pop('base_mapping');report['evidence']['parameters']['evidence_version']=3
        db.execute('UPDATE sanger_analysis_runs SET report_json=? WHERE id=?',(json.dumps(report),first['id']))
    current=view(client,revision);data=export(client,revision,current['snapshot_sha256']).json()['snapshot']
    assert data['summary']['quality_review']['source_runs'][0]['reason']=='mapping_not_recorded'
    assert 'base_mapping' not in data['reads'][0]['evidence']


def test_missing_original_file_does_not_recompute_and_html_is_escaped(ready):
    client,read,revision=ready
    client.post('/api/sanger-verifications',json={'sequencing_read_id':read})
    with connect() as db:
        Path(db.execute('SELECT storage_path FROM sequencing_reads WHERE id=?',(read,)).fetchone()[0]).unlink()
        db.execute('UPDATE sequencing_reads SET original_filename=? WHERE id=?',('<script>alert(1)</script>.ab1',read))
        db.execute('UPDATE sequence_revisions SET label=? WHERE id=?',('<img src=x onerror=alert(1)>',revision))
    current=view(client,revision)
    html=export(client,revision,current['snapshot_sha256'],'html')
    assert html.status_code==200 and '<script>' not in html.text and '<img src=x' not in html.text
    assert '&lt;script&gt;' in html.text and '&lt;img' in html.text
    assert "default-src 'none'" in html.text and 'src="http' not in html.text
    assert 'Original AB1 files are not rechecked' in html.text
    assert 'Full snapshot' in html.text
    assert 'sanger-group-' in html.headers['content-disposition'] and 'alert' not in html.headers['content-disposition']


def test_group_export_requires_valid_digest_and_existing_revision(ready):
    client,read,revision=ready
    url=f'/api/sequence-revisions/{revision}/sanger-report'
    assert client.get(url).status_code==422
    assert export(client,revision,'bad').status_code==422
    digest=view(client,revision)['snapshot_sha256']
    assert export(client,revision,digest,'pdf').status_code==422
    assert export(client,'missing',digest).status_code==404
    assert export(client,revision,'0'*64).status_code==409
