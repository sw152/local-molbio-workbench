"""No consensus from absent, malformed or ambiguous saved mappings."""
from copy import deepcopy
from dataclasses import asdict
from hashlib import sha256
import random
import pytest
from Bio.Seq import Seq
from localmolbio.sanger_verification import align_sanger_read
from localmolbio.sanger_summary import summarize_sanger_reads
from localmolbio.sanger_quality_summary import validated_columns
from test_sanger_history import ready

rng=random.Random(709)
SEQ=''.join(rng.choice('ACGT') for _ in range(420))
REF={'id':'r','length_bp':len(SEQ),'sha256':sha256(SEQ.encode()).hexdigest(),'topology':'linear'}

def read(name,start=40,end=240,mutate=None,quality=35,reverse=False,reference=SEQ,circular=False):
    calls=list((reference+reference)[start:end]);qs=[quality]*len(calls)
    if mutate is not None:calls[mutate]=next(b for b in 'ACGT' if b!=calls[mutate])
    calls=''.join(calls)
    if reverse:calls=str(Seq(calls).reverse_complement());qs=qs[::-1]
    report=asdict(align_sanger_read(reference,calls,'unknown',circular,qs))
    report['evidence']['parameters']['reference_topology']='circular' if circular else 'linear'
    return {'id':name,'alignment_id':name+'-run','run_number':1,'original_filename':name+'.ab1','file_sha256':name+'-file',
            'base_sequence':calls,'qualities':qs,'report':report,
            'manifest':{'read_sha256':name+'-file','reference_sha256':sha256(reference.encode()).hexdigest()}}

def summary(rows,ref=REF,seq=SEQ):return summarize_sanger_reads(ref,rows,seq)['quality_review']

def test_reverse_conflict_and_quality_coverage_are_distinct():
    a=read('a');b=read('b',120,320,mutate=70,reverse=True)
    q=summary([a,b])
    assert q['callable_bases']==280 and q['overlap_bases']==120
    assert q['conflict_position_count']==1 and q['conflict_intervals']==[[190,191]]
    assert q['reference_only_bases']==279
    difference=q['differences'][0]
    assert difference['position']==190 and difference['reference']==SEQ[190]
    assert difference['conflict'] and len(difference['calls'])==2
    call=next(c for c in difference['calls'] if c['read_id']=='b')
    assert call['original_read_position']==129 and call['phred']==35
    assert q['indel_conflicts']=='not_evaluated' and q['consensus']=='not_generated'
    assert not q['whole_reference_verified'] and not q['source_files_rechecked']

def test_low_quality_disagreement_not_a_q20_conflict_and_same_alternate_not_conflict():
    a=read('a');b=read('b',quality=19,mutate=70)
    q=summary([a,b]);assert q['callable_bases']==200 and q['overlap_bases']==0
    assert q['conflict_position_count']==0 and q['differences']==[]
    c=read('c',quality=20,mutate=70);d=read('d',quality=40,mutate=70)
    q=summary([c,d]);assert q['non_reference_position_count']==1 and q['conflict_position_count']==0
    assert q['differences'][0]['conflict'] is False and len(q['differences'][0]['calls'])==2
    assert q['reference_only_bases']==199

@pytest.mark.parametrize('mutation',[
    lambda b:b.update(reference_start=-1),
    lambda b:b.update(reference_end=421),
    lambda b:b.update(original_read_start=True),
    lambda b:b.update(original_read_step=-1),
    lambda b:b.update(read_bases='A'),
    lambda b:b.update(reference_bases='N'*200),
    lambda b:b.update(phred=[35]*199),
    lambda b:b.update(phred=[True]*200),
    lambda b:b.update(phred=[34]*200),
    lambda b:b.update(original_read_start=200),
])
def test_malformed_block_never_contributes(mutation):
    a=read('a');mutation(a['report']['evidence']['base_mapping']['blocks'][0])
    q=summary([a]);assert q['callable_bases']==0 and q['assessment']=='not_evaluated'
    assert q['source_runs'][0]['reason']=='invalid_base_mapping'

@pytest.mark.parametrize('case',['counts','coverage','direction','calls','reference_hash','duplicate_block','trimming','legacy','ambiguous','duplicate_file'])
def test_inconsistent_and_unavailable_sources_stay_explicit(case):
    a=read('a');rows=[a]
    if case=='counts':a['report']['evidence']['base_mapping']['callable_base_count']=999
    elif case=='coverage':a['report']['evidence']['reference_covered_intervals']=[[41,241]]
    elif case=='direction':a['report']['evidence']['direction']='reverse'
    elif case=='calls':a['base_sequence']='N'*200
    elif case=='reference_hash':a['manifest']['reference_sha256']='other'
    elif case=='duplicate_block':a['report']['evidence']['base_mapping']['blocks']*=2
    elif case=='trimming':a['report']['evidence']['trimming']['original_read_end']=3
    elif case=='legacy':a['report']['evidence'].pop('base_mapping');a['report']['evidence']['parameters']['evidence_version']=3
    elif case=='ambiguous':a['report']['evidence']['review_flags']=['ambiguous_alignment']
    else:
        b=deepcopy(a);b['id']='b';rows.append(b)
    q=summary(rows);assert q['included_read_count']==0 and q['callable_bases']==0
    assert all(s['reason'] for s in q['source_runs'])
    if case=='legacy':assert q['source_runs'][0]['reason']=='mapping_not_recorded'

def test_circular_original_order_and_holes_are_preserved():
    a=read('a',360,500,reverse=True,circular=True)
    q=summary([a],{**REF,'topology':'circular'})
    assert q['callable_intervals']==[[0,80],[360,420]] and q['callable_bases']==140
    a['report']['evidence']['base_mapping']['blocks'].reverse()
    assert summary([a],{**REF,'topology':'circular'})['included_read_count']==0

def test_reference_integrity_missing_quality_and_ambiguous_calls():
    a=read('a')
    assert summary([a],seq='N'*420)['source_runs'][0]['reason']=='reference_sequence_unavailable_or_mismatch'
    # Legacy storage may lack qualities recovered at analysis time; keep saved per-base evidence.
    a['qualities']=[]
    assert summary([a])['callable_bases']==200
    a['report']=asdict(align_sanger_read(SEQ,a['base_sequence'],'unknown',False))
    a['report']['evidence']['parameters']['reference_topology']='linear'
    assert summary([a])['callable_bases']==0
    a['base_sequence']=a['base_sequence'][:70]+'N'+a['base_sequence'][71:]
    a['qualities']=[35]*200
    a['report']=asdict(align_sanger_read(SEQ,a['base_sequence'],'unknown',False,a['qualities']))
    a['report']['evidence']['parameters']['reference_topology']='linear'
    assert summary([a])['callable_bases']==199

def test_api_quality_summary_does_not_return_full_inputs(ready):
    client,read_id,revision=ready
    client.post('/api/sanger-verifications',json={'sequencing_read_id':read_id})
    data=client.get(f'/api/sequence-revisions/{revision}/sanger-reads?include_summary=true').json()
    q=data['summary']['quality_review']
    assert q['callable_bases']==170 and q['included_read_count']==1
    assert q['callable_intervals']==[[60,230]]
    assert 'sequence_text' not in data['summary']['reference'] and 'base_sequence' not in data['reads'][0]

@pytest.mark.parametrize('kind',['insertion','deletion'])
@pytest.mark.parametrize('reverse',[False,True])
def test_indels_preserve_paired_coverage_with_only_scoped_indel_review(kind,reverse):
    a=read('a');calls=SEQ[40:160]+('AAA' if kind=='insertion' else '')+SEQ[163 if kind=='deletion' else 160:300]
    if reverse:calls=str(Seq(calls).reverse_complement())
    a['base_sequence']=calls;a['qualities']=[35]*len(calls)
    a['report']=asdict(align_sanger_read(SEQ,calls,'unknown',False,a['qualities']))
    a['report']['evidence']['parameters']['reference_topology']='linear'
    q=summary([a]);assert q['included_read_count']==1
    expected=a['report']['evidence']['reference_covered_intervals']
    assert q['callable_intervals']==expected and q['callable_bases']==a['report']['aligned_bases']
    assert q['conflict_position_count']==0 and q['indel_conflicts']=='exact_event_vs_reference_only'
    comparison=q['indel_review']['comparisons'][0]
    assert comparison['event_read_count']==1 and comparison['reference_read_count']==0
    assert not comparison['conflicting_support'] and comparison['alternative_indel_alleles']=='not_compared'

def test_mixed_legacy_mapping_is_not_silently_assessed():
    a,b=read('a'),read('b',120,320,mutate=70)
    b['report']['evidence'].pop('base_mapping');b['report']['evidence']['parameters']['evidence_version']=3
    q=summary([a,b])
    assert q['included_read_count']==1 and q['excluded_read_count']==1 and q['overlap_bases']==0
    assert q['source_runs'][1]['reason']=='mapping_not_recorded'
    assert q['callable_bases']==200
    a['qualities']={}
    assert summary([a])['source_runs'][0]['reason']=='stored_quality_invalid'
