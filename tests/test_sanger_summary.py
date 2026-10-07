"""Coordinate union is not a quality-screened consensus or pass verdict."""
from copy import deepcopy
import pytest
from localmolbio.sanger_summary import summarize_sanger_reads
from test_sanger_history import ready

REF={'id':'revision','length_bp':100,'sha256':'reference-hash','topology':'linear'}

def read(name,intervals):
    return {'id':name,'original_filename':name+'.ab1','file_sha256':name+'-hash',
            'alignment_id':name+'-run','run_number':1,
            'manifest':{'read_sha256':name+'-hash','reference_sha256':REF['sha256']},
            'report':{'aligned_bases':sum(b-a for a,b in intervals),'evidence':{
                'scope':'local_read_alignment','parameters':{'evidence_version':3},
                'direction':'forward','review_flags':[],'reference_covered_intervals':intervals}}}

def test_empty_overlap_adjacency_and_deletion_holes():
    empty=summarize_sanger_reads(REF,[])
    assert empty['covered_bases']==0 and empty['uncovered_intervals']==[[0,100]]
    reads=[read('a',[[0,20],[21,40]]),read('b',[[30,60]]),read('c',[[60,100]])]
    result=summarize_sanger_reads(REF,reads)
    assert result['covered_intervals']==[[0,20],[21,100]]
    assert result['covered_bases']==99 and result['uncovered_intervals']==[[20,21]]
    assert result['overlap_bases']==10 and result['maximum_read_depth']==2
    assert result['depth_segments']==[{'start':0,'end':20,'depth':1},{'start':21,'end':30,'depth':1},{'start':30,'end':40,'depth':2},{'start':40,'end':100,'depth':1}]
    assert result['snapshot_sha256']==summarize_sanger_reads(REF,deepcopy(reads))['snapshot_sha256']

def test_circular_gaps_and_full_coverage():
    ref={**REF,'topology':'circular'}
    result=summarize_sanger_reads(ref,[read('a',[[10,40],[50,80]])])
    assert result['gap_regions']==[{'segments':[[80,100],[0,10]],'length_bp':30,'wraps_origin':True},{'segments':[[40,50]],'length_bp':10,'wraps_origin':False}]
    full=summarize_sanger_reads(ref,[read('a',[[0,100]])])
    assert full['covered_bases']==100 and not full['gap_regions']
    assert not full['whole_reference_verified']
    assert summarize_sanger_reads(ref,[])['gap_regions']==[{'segments':[[0,100]],'length_bp':100,'wraps_origin':False}]
    unknown=summarize_sanger_reads({**ref,'topology':'unknown'},[read('a',[[10,80]])])
    assert len(unknown['gap_regions'])==2

@pytest.mark.parametrize('intervals',[None,[],[[-1,10]],[[0,101]],[[4,4]],[[True,10]],[[0,10],[9,20]],[[20,30],[0,10]],[[0,10.0]],[[0,10,20]],[[0,9]]])
def test_invalid_intervals_are_excluded(intervals):
    item=read('a',[[0,10]])
    item['report']['evidence']['reference_covered_intervals']=intervals
    result=summarize_sanger_reads(REF,[item])
    assert result['covered_bases']==0
    assert 'invalid_coverage_intervals' in result['source_runs'][0]['exclusion_reasons']

@pytest.mark.parametrize('field,value,reason',[
    ('direction','unknown','unresolved_direction'),
    ('scope',None,'unsupported_or_missing_evidence'),
    ('parameters',{'evidence_version':1},'unsupported_or_missing_evidence'),
    ('parameters',{'evidence_version':3,'reference_topology':'circular'},'reference_topology_mismatch'),
    ('review_flags',['ambiguous_alignment'],'ambiguous_alignment'),
    ('review_flags',['candidate_search_truncated'],'candidate_search_truncated')])
def test_uncertain_placement_is_excluded(field,value,reason):
    item=read('a',[[0,10]])
    item['report']['evidence'][field]=value
    result=summarize_sanger_reads(REF,[item])
    assert result['covered_bases']==0 and reason in result['source_runs'][0]['exclusion_reasons']

@pytest.mark.parametrize('field',['reference_sha256','read_sha256'])
@pytest.mark.parametrize('value',[None,'other'])
def test_manifest_hashes_required(field,value):
    item=read('a',[[0,10]]);item['manifest'][field]=value
    result=summarize_sanger_reads(REF,[item])
    assert result['excluded_read_count']==1 and result['covered_bases']==0

def test_duplicate_sources_unanalyzed_and_quality_warnings():
    a=read('a',[[0,10]]);b=deepcopy(a);b['id']='b'
    duplicate=summarize_sanger_reads(REF,[a,b])
    assert duplicate['covered_bases']==0 and duplicate['excluded_read_count']==2
    assert all('duplicate_source_file' in r['exclusion_reasons'] for r in duplicate['source_runs'])
    a['alignment_id']=None;a['report']={}
    assert summarize_sanger_reads(REF,[a])['source_runs'][0]['exclusion_reasons']==['not_analyzed']
    b['report']['evidence']['review_flags']=['partial_read_alignment','low_quality_aligned_bases','quality_unavailable','ambiguous_bases']
    result=summarize_sanger_reads(REF,[b])
    assert result['covered_bases']==10 and len(result['source_runs'][0]['review_flags'])==4
    assert result['quality_screening']=='not_applied_to_union' and result['conflict_analysis']=='not_evaluated'
    assert not result['source_files_rechecked'] and not result['whole_reference_verified']

def test_api_latest_only_backwards_compatible_and_version_bound(ready):
    client,read_id,revision=ready
    url=f'/api/sequence-revisions/{revision}/sanger-reads'
    empty=client.get(url+'?include_summary=true').json()
    assert empty['summary']['excluded_read_count']==1 and empty['summary']['covered_bases']==0
    first=client.post('/api/sanger-verifications',json={'sequencing_read_id':read_id}).json()
    before=client.get(url+'?include_summary=true').json()['summary']
    second=client.post('/api/sanger-verifications',json={'sequencing_read_id':read_id,'previous_alignment_id':first['id'],'trim_quality_threshold':20}).json()
    response=client.get(url+'?include_summary=true')
    assert response.status_code==200,response.text
    data=response.json();summary=data['summary']
    assert data['reads']==client.get(url).json()
    assert summary['covered_intervals']==[[60,230]] and summary['covered_bases']==170
    assert summary['included_read_count']==1 and summary['overlap_bases']==0
    assert summary['source_runs'][0]['alignment_id']==second['id'] and summary['source_runs'][0]['run_number']==2
    assert summary['snapshot_sha256']!=before['snapshot_sha256']
    assert summary['reference']['id']==revision
    assert client.get('/api/sequence-revisions/missing/sanger-reads?include_summary=true').status_code==404

def test_union_matches_independent_per_position_enumeration():
    import random
    rng=random.Random(401)
    for _ in range(80):
        rows=[];depth=[0]*100
        for index in range(rng.randrange(1,15)):
            a,b=sorted(rng.sample(range(101),2))
            rows.append(read(str(index),[[a,b]]))
            for pos in range(a,b):depth[pos]+=1
        result=summarize_sanger_reads(REF,rows)
        assert result['covered_bases']==sum(d>0 for d in depth)
        assert result['overlap_bases']==sum(d>1 for d in depth)
        assert result['maximum_read_depth']==max(depth)
        reconstructed=[0]*100
        for s in result['depth_segments']:
            reconstructed[s['start']:s['end']]=[s['depth']]*(s['end']-s['start'])
        assert reconstructed==depth
        assert sum(g['length_bp'] for g in result['gap_regions'])==depth.count(0)
