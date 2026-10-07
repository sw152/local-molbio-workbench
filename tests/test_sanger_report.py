import json
from pathlib import Path
from test_sanger_history import ready
from localmolbio.database import connect


def export(client,read,analysis,format='json'):
    return client.get(f'/api/sanger-reads/{read}/analyses/{analysis}/export?format={format}')


def test_specific_run_export_is_snapshot_bound_and_contains_no_raw_files(ready):
    client,read,revision=ready
    first=client.post('/api/sanger-verifications',json={'sequencing_read_id':read}).json()
    first_export=export(client,read,first['id']).json()
    second=client.post('/api/sanger-verifications',json={'sequencing_read_id':read,'previous_alignment_id':first['id'],'trim_quality_threshold':20}).json()
    second_export=export(client,read,second['id']).json()
    assert first_export['analysis']['run_number']==1 and second_export['analysis']['run_number']==2
    assert first_export['saved_alignment']['evidence']['trimming']['method']=='none'
    assert second_export['saved_alignment']['evidence']['trimming']['original_read_start']==10
    assert second_export['parameters']['trim_quality_threshold']==20
    with connect() as db:
        run=db.execute('SELECT report_json FROM sanger_analysis_runs WHERE id=?',(first['id'],)).fetchone()[0]
        stored=Path(db.execute('SELECT storage_path FROM sequencing_reads WHERE id=?',(read,)).fetchone()[0])
        db.execute("UPDATE sequence_revisions SET label='Changed later' WHERE id=?",(revision,))
        db.execute("UPDATE sequencing_reads SET original_filename='Changed later' WHERE id=?",(read,))
    stored.unlink()  # Export remains a historical report, not a fresh source verification.
    again=export(client,read,first['id'])
    assert again.status_code==200
    data=again.json()
    assert data['saved_alignment']==json.loads(run)==first_export['saved_alignment']
    assert data['inputs']==first_export['inputs']
    assert data['inputs']['read']['filename']['source']=='analysis_input_snapshot'
    assert data['inputs']['read']['filename']['value']=='trace.ab1'
    assert data['inputs']['reference']['label']['value']=='history'
    assert len(data['inputs']['reference']['sha256'])==64
    assert data['source_verification_at_export']=='not_performed'
    assert str(stored) not in again.text
    assert 'base_sequence' not in again.text and 'sequence_text' not in again.text and 'parser_metadata' not in again.text
    assert 'attachment;' in again.headers['content-disposition']
    assert again.headers['cache-control']=='no-store'
    assert export(client,'other-read',first['id']).status_code==404
    assert export(client,read,'missing').status_code==404
    assert export(client,read,first['id'],'pdf').status_code==422


def test_html_escapes_metadata_and_never_adds_missing_legacy_hashes(ready):
    client,read,revision=ready
    result=client.post('/api/sanger-verifications',json={'sequencing_read_id':read}).json()
    hostile='<script>alert("sample")</script><img src="https://example.invalid/leak">'
    with connect() as db:
        db.execute("UPDATE analysis_jobs SET input_manifest_json='{}' WHERE id=?",(result['job_id'],))
        db.execute('UPDATE sequencing_reads SET original_filename=? WHERE id=?',(hostile,read))
        db.execute('UPDATE sequence_revisions SET label=? WHERE id=?',(hostile,revision))
        raw=json.loads(db.execute('SELECT report_json FROM sanger_analysis_runs WHERE id=?',(result['id'],)).fetchone()[0])
        raw['evidence']={}
        db.execute('UPDATE sanger_analysis_runs SET report_json=? WHERE id=?',(json.dumps(raw),result['id']))
    data=export(client,read,result['id']).json()
    assert data['inputs']['read']['sha256'] is None and data['inputs']['reference']['sha256'] is None
    assert data['inputs']['read']['filename']['source']=='linked_record_metadata'
    html=export(client,read,result['id'],'html')
    assert html.status_code==200
    assert '<script>' not in html.text and '<img ' not in html.text
    assert '&lt;script&gt;' in html.text
    assert 'Detailed evidence was not recorded' in html.text
    assert 'Coordinate graphic unavailable' in html.text
    assert 'Not recorded' in html.text
    assert 'default-src' in html.headers['content-security-policy']
    assert 'trace.ab1' not in html.headers['content-disposition']


def test_export_retains_boundary_zero_and_absent_deletion_quality(ready):
    client,read,revision=ready
    result=client.post('/api/sanger-verifications',json={'sequencing_read_id':read}).json()
    # Controlled archived-report fixture tests rendering conventions independently of alignment.
    with connect() as db:
        saved=json.loads(db.execute('SELECT report_json FROM sanger_analysis_runs WHERE id=?',(result['id'],)).fetchone()[0])
        saved['variants']=[{'kind':'insertion','position':0,'reference':'','read':'A','original_read_position':9,'phred':8},
                           {'kind':'deletion','position':599,'reference':'C','read':'','original_read_position':None,'phred':None}]
        saved['evidence']['reference_covered_intervals']=[[0,12],[580,600]]
        db.execute('UPDATE sanger_analysis_runs SET report_json=? WHERE id=?',(json.dumps(saved),result['id']))
    html=export(client,read,result['id'],'html').text
    assert 'Boundary after 0' in html
    assert '<td>10</td><td>8</td>' in html
    assert '<td>600</td><td>deletion</td><td>C → –</td><td>Not recorded</td><td>Not recorded</td>' in html
    assert '1–12; 581–600' in html
    assert export(client,read,result['id']).json()['saved_alignment']==saved
