"""Saved base-level evidence must recover original calls, including reverse/indels."""
import random
import pytest
from Bio.Seq import Seq
from localmolbio.sanger_verification import align_sanger_read, SangerVerificationError
from test_sanger_history import ready

rng=random.Random(564)
REF=''.join(rng.choice('ACGT') for _ in range(500))

def expand(mapping):
    for block in mapping['blocks']:
        size=block['reference_end']-block['reference_start']
        assert len(block['read_bases'])==len(block['reference_bases'])==len(block['phred'])==size
        for i in range(size):
            yield (block['reference_start']+i,block['original_read_start']+i*block['original_read_step'],block['reference_bases'][i],block['read_bases'][i],block['phred'][i])

@pytest.mark.parametrize('reverse',[False,True])
@pytest.mark.parametrize('kind',['ordinary','insertion','deletion','circular'])
def test_lossless_mapping_through_trim_orientation_gaps_and_origin(reverse,kind):
    if kind=='circular':read=REF[-80:]+REF[:180]
    else:read=REF[40:160]+('AAA' if kind=='insertion' else '')+REF[163 if kind=='deletion' else 160:300]
    read=list(read);read[70]=next(b for b in 'ACGT' if b!=read[70]);read[90]='N';read=''.join(read)
    qualities=[(i%30)+10 for i in range(len(read))]
    qualities[:7]=[2]*7;qualities[7]=20;qualities[-11:]=[2]*11;qualities[-12]=20
    if reverse:read=str(Seq(read).reverse_complement());qualities=qualities[::-1]
    result=align_sanger_read(REF,read,'unknown',kind=='circular',qualities,20)
    mapping=result.evidence['base_mapping'];columns=list(expand(mapping))
    assert result.evidence['parameters']['evidence_version']==4
    assert len(columns)==mapping['aligned_base_count']==result.aligned_bases
    assert len({p for p,*_ in columns})==len(columns)
    trim=result.evidence['trimming']
    for p,original,ref,call,q in columns:
        assert trim['original_read_start']<=original<trim['original_read_end']
        assert ref==REF[p] and q==qualities[original]
        assert call==(str(Seq(read[original]).complement()) if reverse else read[original])
    for v in result.variants:
        if v['kind'] in {'substitution','ambiguous'}:
            assert (v['position'],v['original_read_position'],v['reference'],v['read'],v['phred']) in columns
        elif v['kind']=='deletion':assert v['position'] not in {c[0] for c in columns}
        else:assert v['original_read_position'] not in {c[1] for c in columns}
    expected={p for p,o,r,c,q in columns if r in 'ACGT' and c in 'ACGT' and q>=20}
    assert {p for a,b in mapping['callable_reference_intervals'] for p in range(a,b)}==expected
    assert sum(mapping['categories'].values())==len(columns)
    assert mapping['callable_base_count']==len(expected)
    assert mapping['categories']['matching']+mapping['categories']['differing']==len(expected)
    assert not mapping['whole_reference_verified']
    if kind!='ordinary':assert len(mapping['blocks'])>=2

def test_quality_categories_include_threshold_equal_and_matching_low_quality():
    read=list(REF);read[200]=next(b for b in 'ACGT' if b!=read[200]);read[220]='N';read=''.join(read)
    q=[30]*500;q[10]=19;q[11]=20;q[200]=20;q[220]=0
    result=align_sanger_read(REF,read,'forward',False,q)
    m=result.evidence['base_mapping']
    assert m['categories']=={'matching':497,'differing':1,'low_quality':1,'quality_unavailable':0,'ambiguous':1}
    assert m['callable_base_count']==498 and m['differing_reference_intervals']==[[200,201]]
    missing=align_sanger_read(REF,read,'forward',False).evidence['base_mapping']
    assert missing['categories']['quality_unavailable']==499 and missing['categories']['ambiguous']==1
    assert not missing['callable_reference_intervals']
    assert all(q is None for b in missing['blocks'] for q in b['phred'])

def test_partial_or_ambiguous_placement_is_not_promoted_to_validation():
    part=align_sanger_read(REF,'N'*100+REF[100:150]+'N'*100,'unknown',False,[35]*250)
    assert part.evidence['base_mapping']['callable_base_count']==50
    assert 'partial_read_alignment' in part.evidence['review_flags']
    repeat=align_sanger_read(REF+REF,REF[40:180],'unknown',False,[35]*140)
    assert 'ambiguous_alignment' in repeat.evidence['review_flags']
    assert repeat.evidence['base_mapping']['placement_eligibility']=='see_report_review_flags'
    assert not repeat.evidence['base_mapping']['whole_reference_verified']

def test_boolean_quality_is_rejected():
    with pytest.raises(SangerVerificationError,match='Quality'):
        align_sanger_read(REF,REF,'forward',False,[True]*len(REF))

def test_saved_mapping_survives_history_and_json_export(ready):
    client,read,revision=ready
    first=client.post('/api/sanger-verifications',json={'sequencing_read_id':read}).json()
    second=client.post('/api/sanger-verifications',json={'sequencing_read_id':read,'previous_alignment_id':first['id'],'trim_quality_threshold':20}).json()
    history=client.get(f'/api/sanger-reads/{read}/analyses').json()['items']
    assert history[1]['evidence']['base_mapping']==first['evidence']['base_mapping']
    assert history[0]['evidence']['base_mapping']==second['evidence']['base_mapping']
    assert first['evidence']['base_mapping']['categories']['low_quality']==30
    assert second['evidence']['base_mapping']['categories']['low_quality']==0
    exported=client.get(f'/api/sanger-reads/{read}/analyses/{first["id"]}/export?format=json').json()
    assert exported['saved_alignment']['evidence']['base_mapping']==first['evidence']['base_mapping']
    summary=client.get(f'/api/sequence-revisions/{revision}/sanger-reads?include_summary=true').json()['summary']
    assert summary['included_read_count']==1 and summary['covered_bases']==170

def test_v3_report_remains_unchanged_after_new_mapping_run(ready):
    import json
    from localmolbio.database import connect
    client,read,revision=ready
    first=client.post('/api/sanger-verifications',json={'sequencing_read_id':read}).json()
    # Emulate an archived v3 report; the application must not retrofit it.
    with connect() as db:
        row=db.execute('SELECT * FROM sanger_analysis_runs WHERE id=?',(first['id'],)).fetchone()
        old=json.loads(row['report_json']);old['evidence'].pop('base_mapping');old['evidence']['parameters']['evidence_version']=3
        encoded=json.dumps(old)
        db.execute('UPDATE sanger_analysis_runs SET report_json=? WHERE id=?',(encoded,first['id']))
        params=json.loads(db.execute('SELECT parameters_json FROM analysis_jobs WHERE id=?',(row['analysis_job_id'],)).fetchone()[0]);params['evidence_version']=3
        db.execute('UPDATE analysis_jobs SET parameters_json=? WHERE id=?',(json.dumps(params),row['analysis_job_id']))
    before=client.get(f'/api/sequence-revisions/{revision}/sanger-reads?include_summary=true').json()
    assert before['summary']['included_read_count']==1
    assert 'base_mapping' not in before['reads'][0]['evidence']
    result=client.post('/api/sanger-verifications',json={'sequencing_read_id':read,'previous_alignment_id':first['id']})
    assert result.status_code==201 and result.json()['evidence']['base_mapping']['schema_version']==1
    with connect() as db:
        assert db.execute('SELECT report_json FROM sanger_analysis_runs WHERE id=?',(first['id'],)).fetchone()[0]==encoded
