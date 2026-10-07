"""Reference support requires a fully observed Q20 matching contiguous read path."""
from copy import deepcopy
from dataclasses import asdict
from hashlib import sha256
import pytest
from Bio.Seq import Seq
from localmolbio.sanger_summary import summarize_sanger_reads
from localmolbio.sanger_verification import align_sanger_read
from localmolbio.sanger_indel_comparison import compare_indels
from test_sanger_indel_review import fixture, REF, manual_event
from localmolbio.sanger_indel_review import review_indels


def record(read,qualities=None,reverse=False):
    qualities=qualities if qualities is not None else [35]*len(read)
    if reverse:read=str(Seq(read).reverse_complement());qualities=qualities[::-1]
    report=asdict(align_sanger_read(REF,read,'unknown',False,qualities))
    report['evidence']['parameters']['reference_topology']='linear'
    return {'base_sequence':read,'qualities':qualities,'report':report,'run_number':1}


def compare(rows):
    reference={'id':'rev','length_bp':len(REF),'sha256':sha256(REF.encode()).hexdigest(),'topology':'linear'}
    for i,r in enumerate(rows):
        r.setdefault('id',str(i));r.setdefault('alignment_id','run-'+str(i));r.setdefault('original_filename',str(i)+'.ab1')
        r.setdefault('file_sha256','file-'+str(i));r['manifest']={'read_sha256':r['file_sha256'],'reference_sha256':reference['sha256']}
    return summarize_sanger_reads(reference,rows,REF)['quality_review']['indel_review']

@pytest.mark.parametrize('kind',['insertion','deletion'])
@pytest.mark.parametrize('reverse',[False,True])
def test_event_and_contiguous_reference_have_conflicting_support(kind,reverse):
    event,_,start=fixture(kind,reverse)
    reference=record(REF[40:300],reverse=not reverse)
    result=compare([event,reference]);group=result['comparisons'][0]
    assert group['event_read_count']==1 and group['reference_read_count']==1 and group['unassessed_read_count']==0
    assert group['conflicting_support'] and result['conflicting_group_count']==1
    source=next(s for s in group['sources'] if s['support']=='reference')
    assert source['reference_span_bases']==(2 if kind=='insertion' else 5)
    assert source['minimum_span_phred']==35
    assert source['flanks']['left']['reference_position']==start-1
    assert group['scope']=='exact_event_vs_reference' and not group['whole_reference_verified']

@pytest.mark.parametrize('kind',['insertion','deletion'])
@pytest.mark.parametrize('case',['partial','q19','missing_quality','N','mismatch','legacy','absent_event_entries'])
def test_missing_or_uncertain_reference_evidence_is_unassessed(kind,case):
    event,_,start=fixture(kind)
    read=list(REF[40:300]);quality=[35]*len(read)
    # For deletions this is an internal deleted base: checking flanks alone would be wrong.
    position=start if kind=='deletion' else start-1
    if case=='partial':read=read[:start-40]
    elif case=='q19':quality[position-40]=19
    elif case=='N':read[position-40]='N'
    elif case=='mismatch':read[position-40]=next(b for b in 'ACGT' if b!=read[position-40])
    reference=record(''.join(read),quality[:len(read)])
    if case=='legacy':reference['report']['evidence'].pop('base_mapping');reference['report']['evidence']['parameters']['evidence_version']=3
    if case=='missing_quality':
        reference=record(REF[40:300]);reference['qualities']=[]
        reference['report']=asdict(align_sanger_read(REF,REF[40:300],'unknown',False))
        reference['report']['evidence']['parameters']['reference_topology']='linear'
    if case=='absent_event_entries':reference['report']['inserted_bases']=1
    group=compare([event,reference])['comparisons'][0]
    assert group['event_read_count']==1 and group['reference_read_count']==0 and group['unassessed_read_count']==1
    assert not group['conflicting_support'] and group['sources'][1]['reason']


def test_two_identical_events_count_once_per_read_and_q20_is_inclusive():
    first,_,_=fixture('insertion');second=deepcopy(first);second.update(id='second',alignment_id='second-run')
    reference=record(REF[40:300],[20]*260)
    group=compare([first,second,reference])['comparisons'][0]
    assert group['event_read_count']==2 and group['reference_read_count']==1
    assert len(group['sources'])==3


def test_other_insertion_allele_is_never_counted_as_reference():
    first,_,start=fixture('insertion');other=first['base_sequence']
    index=start-40+1;other=other[:index]+('T' if other[index]!='T' else 'C')+other[index+1:]
    second=record(other)
    groups=compare([first,second])['comparisons']
    assert len(groups)==2
    assert all(g['event_read_count']==1 and g['reference_read_count']==0 and g['unassessed_read_count']==1 for g in groups)
    assert all(not g['conflicting_support'] and g['alternative_indel_alleles']=='not_compared' for g in groups)


def test_duplicate_source_files_cannot_manufacture_reference_support():
    event,_,_=fixture('deletion');a=record(REF[40:300]);b=deepcopy(a)
    a['file_sha256']=b['file_sha256']='same-file'
    group=compare([event,a,b])['comparisons'][0]
    assert group['reference_read_count']==0 and group['unassessed_read_count']==2

@pytest.mark.parametrize('kind',['insertion','deletion'])
def test_circular_origin_reference_walk_is_explicit(kind):
    reference='ATGTCGATGC';boundary=0 if kind=='insertion' else 9;sequence='TG' if kind=='insertion' else 'CA'
    event,columns=manual_event(reference,kind,boundary,sequence,circular=True)
    events=review_indels(event,reference,'circular',columns)
    assert events[0]['eligible_for_exact_anchor_review']
    other={'id':'other','alignment_id':'other-run','run_number':1,'original_filename':'other.ab1','report':{'evidence':{'direction':'reverse'},'variants':[],'inserted_bases':0,'deleted_bases':0}}
    positions=[9,0] if kind=='insertion' else [8,9,0,1]
    mapped=[(p,20-i,reference[p],reference[p],35) for i,p in enumerate(positions)]
    sources=[{'read_id':r['id'],'included':True,'reason':None} for r in [event,other]]
    result=compare_indels(events,[event,other],{'r':columns,'other':mapped},sources,reference,'circular')[0]
    assert result['conflicting_support'] and result['sources'][1]['reference_span_bases']==len(positions)
    # A jump in original coordinates cannot be treated as spanning reference support.
    mapped[-1]=(mapped[-1][0],0,*mapped[-1][2:])
    result=compare_indels(events,[event,other],{'r':columns,'other':mapped},sources,reference,'circular')[0]
    assert result['reference_read_count']==0 and result['sources'][1]['reason']=='reference_span_not_contiguous_in_read'
