"""Exact position-set oracles for geometry, ambiguity, gaps and API publication."""
import json, random
import pytest
from localmolbio.alignment_evidence import parse_paf
from localmolbio.alignment_coverage import summarize, regions, CoverageError
from localmolbio.database import connect
from localmolbio import job_queue as q
from test_alignment_evidence import line, minimap2
from test_sanger_history import ready
from test_fastq_alignment import setup,register,run
from test_sanger_acceptance import sibling_revision


def evidence(ref,reads,lines,topology='linear'):
    return dict(reference_length_bp=len(ref),topology=topology,reads=parse_paf(lines,ref,reads,topology))


def depths(summary):
    return [s['depths'] for s in summary['segments'] for _ in range(s['start'],s['end'])]


def test_insertions_deletions_n_and_mismatch_are_not_verified_matches():
    e=evidence('ACTAN',['ACGCN'],[line(5,0,5,'+',5,0,5,2,5,'2=1I1D1X1=',nn=1)])
    s=summarize(e)
    assert s['paired_bases']==4 and s['unpaired_bases']==1 and s['deletion_evidence_bases']==1
    assert regions(s,'unpaired')==regions(s,'deletion')==[[2,3]]
    assert depths(s)==[[1,0,0,0],[1,0,0,0],[0,0,1,0],[1,0,0,0],[1,0,0,0]]
    assert not s['base_quality_used'] and not s['whole_reference_verified'] and not s['uniqueness_assessed']


def test_reverse_end_and_circular_deletion_across_origin():
    e=evidence('A'*8,['T'*4],[line(4,0,4,'-',8,4,8,4,4,'4=')])
    assert regions(summarize(e),'unpaired')==[[0,4]]
    e=evidence('A'*8,['AAAA'],[line(4,0,4,'+',16,5,12,4,7,'2=3D2=')],'circular')
    s=summarize(e)
    assert regions(s,'deletion')==[[0,2],[7,8]]
    assert regions(s,'unpaired')==[[0,2],[4,5],[7,8]]
    assert s['paired_bases']==4


def test_alternatives_overlap_count_one_record_not_two_candidates():
    e=evidence('A'*10,['A'*6],[line(6,0,6,'+',10,0,6,6,6,'6='),line(6,0,6,'+',10,3,9,6,6,'6=',tp='S')])
    s=summarize(e)
    assert s['ambiguous_only_bases']==9 and s['single_reported_alignment_bases']==0
    assert s['max_record_depth']==1 and s['mean_record_depth']==pytest.approx(.9)
    assert regions(s,'ambiguous_only')==[[0,9]]
    assert not s['independent_molecules_assessed']


def test_paired_and_deletion_candidates_remain_separate_and_deduplicated():
    e=evidence('A'*10,['A'*6],[line(6,0,6,'+',10,0,8,6,8,'3=2D3='),line(6,0,6,'+',10,1,9,6,8,'2=2D4=',tp='S')])
    s=summarize(e)
    assert depths(s)[3]==[0,0,0,1] and depths(s)[4]==[0,0,0,1]
    assert max(d[1] for d in depths(s))==1
    assert regions(s,'deletion')==[[3,5]]


def test_duplicate_circular_copies_do_not_make_ambiguity_or_depth():
    e=evidence('A'*10,['A'*4],[line(4,0,4,'+',20,2,6,4,4,'4='),line(4,0,4,'+',20,12,16,4,4,'4=',tp='S')],'circular')
    s=summarize(e)
    assert s['single_reported_alignment_bases']==4 and s['max_record_depth']==1


def test_single_record_category_is_not_uniqueness_and_withheld_never_adds_positions():
    e=evidence('A'*10,['A'*15], [line(15,0,15,'+',20,0,15,15,15,'15='),line(15,0,5,'+',20,2,7,5,5,'5=',tp='S')],'circular')
    s=summarize(e)
    assert s['read_counts']['withheld']==1 and s['ambiguous_only_bases']==5
    assert s['unpaired_bases']==5 and s['single_reported_alignment_bases']==0
    e=evidence('A'*10,['A'*15],[line(15,0,15,'+',20,0,15,15,15,'15=')],'circular')
    s=summarize(e);assert s['paired_bases']==0 and s['read_counts']['no_alignment_reported']==0
    empty=summarize(evidence('A'*10,['A'],[]));assert empty['read_counts']['no_alignment_reported']==1
    assert regions(empty,'unpaired')==[[0,10]]


def test_position_set_oracle_many_overlapping_records_and_candidates():
    for seed in range(20):
        rng=random.Random(seed);lines=[];expected_single=[0]*37;expected_multi=[0]*37;reads=[]
        for i in range(12):
            size=rng.randint(2,12);reads.append('A'*size);starts=rng.sample(range(38-size),rng.randint(1,3))
            positions=set()
            for start in starts:
                lines.append(line(size,0,size,'+',37,start,start+size,size,size,f'{size}=').replace('q0\t',f'q{i}\t'))
                positions.update(range(start,start+size))
            target=expected_single if len(starts)==1 else expected_multi
            for pos in positions:target[pos]+=1
        s=summarize(evidence('A'*37,reads,lines));observed=depths(s)
        assert [d[0] for d in observed]==expected_single
        assert [d[1] for d in observed]==expected_multi
        assert sum(s[k] for k in ['single_reported_alignment_bases','ambiguous_only_bases','unpaired_bases'])==37


@pytest.mark.parametrize('change',['span','blocks','cigar','start','strand','query','deleted','topology'])
def test_invalid_geometry_rejects_entire_summary(change):
    e=evidence('A'*10,['AAAA'],[line(4,0,4,'+',10,3,7,4,4,'4=')]);h=e['reads'][0]['alignments'][0]
    if change=='span':h['reference_span_bp']=11
    elif change=='blocks':h['paired_blocks'][0]['reference_intervals']=[[3,8]]
    elif change=='cigar':h['cigar']='4M'
    elif change=='start':h['reference_intervals']=[[9,13]]
    elif change=='strand':h['strand']='?'
    elif change=='query':h['query_start']=True
    elif change=='deleted':h['deleted_bases']=1
    else:e['topology']='unknown'
    with pytest.raises(CoverageError):summarize(e)


def test_api_real_deletion_pagination_revision_and_provenance(setup):
    client,revision,ref=setup;item=register(revision,[ref[50:200]+ref[207:500]])
    root=f'/api/revisions/{revision}/alignments'
    response=client.post(root,json={'input_ids':[item['id']],'data_type':'ont-high-accuracy','idempotency_key':'cov'})
    job=response.json()['id'];url=root+'/'+job+'/coverage'
    assert client.get(url).status_code==409
    assert run()['status']=='succeeded'
    response=client.get(url+'?limit=1');assert response.status_code==200,response.text
    s=response.json()['summary'];assert s['paired_bases']==443 and s['deletion_evidence_bases']==7
    assert response.json()['regions']['total']==3
    spans=[client.get(url+f'?limit=1&offset={i}').json()['regions']['items'][0] for i in range(3)]
    assert [(r['start'],r['end']) for r in spans]==[(0,50),(200,207),(500,600)]
    assert client.get(url+'?kind=deletion').json()['regions']['items']==[{'start':200,'end':207,'length_bp':7}]
    assert client.get(url+'?kind=ambiguous_only').json()['regions']['items']==[]
    assert 'segments' not in s and not s['whole_reference_verified']
    other=sibling_revision(revision);assert client.get(url.replace(revision,other)).status_code==404
    foreign=q.enqueue(revision,'other','unknown',{})
    assert client.get(root+'/'+foreign+'/coverage').status_code==404
    for query in ['kind=unknown','limit=0','offset=-1','limit=101']:assert client.get(url+'?'+query).status_code==422
    with connect() as db:
        result=json.loads(db.execute('SELECT result_summary_json FROM analysis_jobs WHERE id=?',(job,)).fetchone()[0])
        result['sources'][0]['record_ordinal']=20
        db.execute('UPDATE analysis_jobs SET result_summary_json=? WHERE id=?',(json.dumps(result),job))
    assert client.get(url).status_code==409


def test_another_record_pairs_across_reported_deletion_without_consensus_claim():
    e=evidence('A'*8,['A'*4,'A'*8],
               [line(4,0,4,'+',8,1,7,4,6,'2=2D2='),
                line(8,0,8,'+',8,0,8,8,8,'8=').replace('q0\t','q1\t')])
    s=summarize(e)
    assert s['paired_bases']==8 and s['unpaired_bases']==0 and s['deletion_evidence_bases']==2
    assert depths(s)[3]==[1,0,1,0] and regions(s,'deletion')==[[3,5]]
    assert not s['whole_reference_verified'] and not s['independent_molecules_assessed']


def test_region_half_open_reverse_and_circular_boundaries():
    from localmolbio.alignment_coverage import select_region
    e=evidence('A'*8,['T'*4],[line(4,0,4,'-',8,4,8,4,4,'4=')])
    assert select_region(e,0,4)=={}
    assert select_region(e,7,8,'paired')[0]['paired_intervals']==[[7,8]]
    e=evidence('A'*8,['AAAA'],[line(4,0,4,'+',16,5,12,4,7,'2=3D2=')],'circular')
    assert select_region(e,0,2,'paired')=={}
    assert select_region(e,0,2,'deletion')[0]['deletion_intervals']==[[0,2]]
    assert select_region(e,7,8,'deletion')[0]['deletion_intervals']==[[7,8]]
    assert select_region(e,2,4)[0]['paired_intervals']==[[2,4]]
    assert select_region(e,4,5)=={}


def test_region_candidates_deduplicate_records_but_preserve_contrasting_evidence():
    from localmolbio.alignment_coverage import select_region
    e=evidence('A'*10,['A'*6],[line(6,0,6,'+',10,0,8,6,8,'3=2D3='),line(6,0,6,'+',10,3,9,6,6,'6=',tp='S')])
    paired=select_region(e,3,5,'paired');deleted=select_region(e,3,5,'deletion')
    assert list(paired)==list(deleted)==[0]
    assert paired[0]['matching_hit_indices']==[1] and deleted[0]['matching_hit_indices']==[0]
    assert paired[0]['paired_intervals']==paired[0]['deletion_intervals']==[[3,5]]
    assert paired[0]['placement_category']=='ambiguous_or_withheld'
    assert select_region(e,3,5)[0]['matching_hit_indices']==[0,1]
    e['reads'][0]['alignments'][0]['cigar']='bad'
    with pytest.raises(CoverageError):select_region(e,8,9,'paired')


@pytest.mark.parametrize('start,end,relation',[(0,0,'either'),(5,3,'either'),(-1,2,'either'),(0,11,'either'),(True,3,'either'),(0,3,'unknown')])
def test_region_invalid_queries(start,end,relation):
    from localmolbio.alignment_coverage import select_region
    with pytest.raises(CoverageError,match='invalid_reference_region'):
        select_region(evidence('A'*10,['A'],[]),start,end,relation)


def test_region_withheld_only_has_no_inferred_positions():
    from localmolbio.alignment_coverage import select_region
    e=evidence('A'*10,['A'*15],[line(15,0,15,'+',20,0,15,15,15,'15=')],'circular')
    assert select_region(e,0,10)=={}


def test_region_api_filters_before_paging_preserves_original_sources_and_scope(setup):
    client,revision,ref=setup
    item=register(revision,[ref[300:550],ref[50:200]+ref[207:500],ref[50:500],ref[30:480]])
    root=f'/api/revisions/{revision}/alignments'
    job=client.post(root,json={'input_ids':[item['id']],'data_type':'ont-high-accuracy','idempotency_key':'region'}).json()['id']
    url=root+'/'+job+'/reads';query='?start=200&end=207&limit=1'
    assert client.get(url+query).status_code==409
    assert run()['status']=='succeeded'
    pages=[client.get(url+query+f'&offset={i}').json() for i in range(3)]
    assert all(p['total']==3 and p['unfiltered_total']==4 for p in pages)
    assert [p['items'][0]['source']['record_ordinal'] for p in pages]==[2,3,4]
    assert [p['items'][0]['source']['query_name'] for p in pages]==['q1','q2','q3']
    assert pages[0]['items'][0]['region_evidence']['deletion_intervals']==[[200,207]]
    assert pages[1]['items'][0]['region_evidence']['paired_intervals']==[[200,207]]
    assert [p['has_more'] for p in pages]==[True,True,False]
    assert client.get(url+query+'&relation=paired').json()['total']==2
    assert client.get(url+query+'&relation=deletion').json()['total']==1
    assert client.get(url+'?start=550&end=600').json()['total']==0
    full=client.get(url).json();assert full['total']==4 and full['region'] is None
    assert 'region_evidence' not in full['items'][0] and not full['whole_reference_verified']
    other=sibling_revision(revision);assert client.get(url.replace(revision,other)+query).status_code==404
    for invalid in ['start=1','end=5','relation=paired','start=0&end=0','start=5&end=3','start=0&end=601','start=-1&end=5','start=0&end=5&relation=invalid']:
        assert client.get(url+'?'+invalid).status_code==422,invalid
    with connect() as db:
        result=json.loads(db.execute('SELECT result_summary_json FROM analysis_jobs WHERE id=?',(job,)).fetchone()[0])
        result['reads'][0]['alignments'][0]['cigar']='bad'
        db.execute('UPDATE analysis_jobs SET result_summary_json=? WHERE id=?',(json.dumps(result),job))
    assert client.get(url+query).status_code==409
