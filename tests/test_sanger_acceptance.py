"""Cross-revision and mixed-version acceptance using real synthetic ABIF parsing."""
import json
from localmolbio.database import connect
from abif_fixture import synthetic_ab1
from test_sanger_history import ready
from test_sanger_group_report import view, export


def sibling_revision(revision):
    # Editing UI is not implemented; seed another immutable revision directly.
    with connect() as db:
        row=dict(db.execute('SELECT * FROM sequence_revisions WHERE id=?',(revision,)).fetchone())
        row.update(id='acceptance-revision-2',parent_revision_id=revision,revision_number=2,source_kind='manual-edit')
        db.execute(f"INSERT INTO sequence_revisions ({','.join(row)}) VALUES ({','.join('?' for _ in row)})",tuple(row.values()))
    return row['id']


def attach(client,revision,calls,quality,name):
    response=client.post('/api/sanger-reads',data={'sequence_revision_id':revision},files={'file':(name,synthetic_ab1(calls,quality,trace=True))})
    assert response.status_code==201,response.text
    return response.json()['id']


def test_same_sequence_revisions_keep_uploads_runs_and_exports_isolated(ready):
    client,read,revision=ready
    first=client.post('/api/sanger-verifications',json={'sequencing_read_id':read}).json()
    before=view(client,revision)
    other=sibling_revision(revision)
    empty=view(client,other)
    assert empty['reads']==[] and empty['summary']['quality_review']['callable_bases']==0
    assert export(client,other,before['snapshot_sha256']).status_code==409
    with connect() as db:
        row=db.execute('SELECT base_sequence,qualities_json FROM sequencing_reads WHERE id=?',(read,)).fetchone()
    second_read=attach(client,other,row[0],json.loads(row[1]),'same-file.ab1')
    second=client.post('/api/sanger-verifications',json={'sequencing_read_id':second_read}).json()
    after=view(client,other)
    assert [r['id'] for r in after['reads']]==[second_read]
    assert after['summary']['reference']['id']==other
    assert after['summary']['quality_review']['included_read_count']==1
    assert view(client,revision)==before
    assert before['snapshot_sha256']!=after['snapshot_sha256']
    assert export(client,other,after['snapshot_sha256']).status_code==200
    assert client.get(f'/api/sanger-reads/{read}/analyses/{second["id"]}/export').status_code==404
    assert client.get(f'/api/sanger-reads/{second_read}/analyses/{first["id"]}/export').status_code==404
    rejected=client.post('/api/sanger-verifications',json={'sequencing_read_id':second_read,'previous_alignment_id':first['id']})
    assert rejected.status_code==409
    assert view(client,other)==after


def test_mixed_current_legacy_and_pending_reads_keep_explicit_assessment(ready):
    client,read,revision=ready
    current=client.post('/api/sanger-verifications',json={'sequencing_read_id':read}).json()
    with connect() as db:
        reference=db.execute('SELECT sequence_text FROM sequence_revisions WHERE id=?',(revision,)).fetchone()[0]
    legacy=attach(client,revision,reference[150:350],[35]*200,'legacy.ab1')
    saved=client.post('/api/sanger-verifications',json={'sequencing_read_id':legacy}).json()
    pending=attach(client,revision,reference[350:550],[35]*200,'pending.ab1')
    with connect() as db:
        report=json.loads(db.execute('SELECT report_json FROM sanger_analysis_runs WHERE id=?',(saved['id'],)).fetchone()[0])
        report['evidence'].pop('base_mapping')
        report['evidence']['parameters']['evidence_version']=3
        db.execute('UPDATE sanger_analysis_runs SET report_json=? WHERE id=?',(json.dumps(report),saved['id']))
    snapshot=view(client,revision)
    q=snapshot['summary']['quality_review']
    assert q['included_read_count']==1 and q['excluded_read_count']==2
    assert q['callable_bases']==170 and not q['whole_reference_verified']
    sources={s['read_id']:s for s in q['source_runs']}
    assert sources[read]['included']
    assert sources[legacy]['reason']=='mapping_not_recorded'
    assert not sources[pending]['included']
    exported=export(client,revision,snapshot['snapshot_sha256']).json()['snapshot']
    assert exported==snapshot
    assert {r['alignment_id'] for r in exported['reads']}=={current['id'],saved['id'],None}
    html=export(client,revision,snapshot['snapshot_sha256'],'html')
    assert html.status_code==200
    assert all(name in html.text for name in ['trace.ab1','legacy.ab1','pending.ab1'])
