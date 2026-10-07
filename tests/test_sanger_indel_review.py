"""Exact-anchor eligibility is deliberately narrower than an indel conflict verdict."""
from dataclasses import asdict
from copy import deepcopy
import random
import pytest
from Bio.Seq import Seq
from localmolbio.sanger_verification import align_sanger_read
from localmolbio.sanger_quality_summary import validated_columns
from localmolbio.sanger_indel_review import review_indels

rng=random.Random(833)
REF=''.join(rng.choice('ACGT') for _ in range(400))


def fixture(kind,reverse=False,circular=False):
    start=150
    # Select flanks that cannot shift by one base; use a three-base insertion.
    if kind=='insertion':
        inserted=next(b for b in 'ACGT' if b!=REF[start])+'G'+next(b for b in 'ACGT' if b!=REF[start-1])
        query=REF[40:start]+inserted+REF[start:300]
    else:
        start=next(p for p in range(100,200) if REF[p-1]!=REF[p+2] and REF[p+3]!=REF[p])
        query=REF[40:start]+REF[start+3:300]
    quality=[35]*len(query)
    if reverse:query=str(Seq(query).reverse_complement());quality=quality[::-1]
    report=asdict(align_sanger_read(REF,query,'unknown',circular,quality))
    report['evidence']['parameters']['reference_topology']='circular' if circular else 'linear'
    row={'id':'r','alignment_id':'a','run_number':1,'original_filename':'synthetic.ab1','base_sequence':query,'qualities':quality,'report':report}
    columns,reason=validated_columns(row,REF)
    assert reason is None
    return row,columns,start

@pytest.mark.parametrize('kind',['insertion','deletion'])
@pytest.mark.parametrize('reverse',[False,True])
def test_groups_three_bases_and_recovers_original_flanks(kind,reverse):
    row,columns,start=fixture(kind,reverse)
    events=review_indels(row,REF,'linear',columns)
    assert len(events)==1
    event=events[0]
    assert event['kind']==kind and event['length_bp']==3 and event['boundary']==start
    assert event['eligible_for_exact_anchor_review'] and event['review_reasons']==[]
    assert event['flanks']['left']['reference_position']==start-1
    assert event['flanks']['right']['reference_position']==(start if kind=='insertion' else start+3)
    left,right=event['flanks']['left']['original_read_position'],event['flanks']['right']['original_read_position']
    assert right-left==(-1 if reverse else 1)*(4 if kind=='insertion' else 1)
    assert event['cross_read_comparison']=='not_performed'
    if kind=='deletion':assert event['inserted_phred']==[] and event['original_read_positions']==[]
    else:assert event['inserted_phred']==[35]*3

@pytest.mark.parametrize('issue',['missing_left','missing_right','low_flank','mismatched_flank','low_insert','missing_insert_quality','invalid_event','complex_gap'])
def test_uncertain_events_not_promoted(issue):
    row,columns,start=fixture('insertion')
    variants=row['report']['variants']
    if issue=='missing_left':columns=[c for c in columns if c[0]!=start-1]
    elif issue=='missing_right':columns=[c for c in columns if c[0]!=start]
    elif issue in {'low_flank','mismatched_flank'}:
        columns=[(p,o,r,('N' if issue=='mismatched_flank' else c),(19 if issue=='low_flank' else q)) if p==start-1 else (p,o,r,c,q) for p,o,r,c,q in columns]
    elif issue in {'low_insert','missing_insert_quality'}:
        v=variants[0];v['phred']=19 if issue=='low_insert' else None
        row['qualities']=[]
    elif issue=='invalid_event':variants[0]['position']=-1
    else:columns=[(p,o+1,r,c,q) if p==start else (p,o,r,c,q) for p,o,r,c,q in columns]
    events=review_indels(row,REF,'linear',columns)
    assert events and not any(e['eligible_for_exact_anchor_review'] for e in events)


def manual_event(reference,kind,boundary,sequence,reverse=False,circular=False):
    # Explicit alignment representation: independent of aligner's repeat traceback choice.
    L=len(reference);end=boundary if kind=='insertion' else boundary+len(sequence)
    left=(boundary-1)%L;right=end%L
    step=-1 if reverse else 1;left_o=20 if reverse else 4
    insertion_positions=[left_o+step*(i+1) for i in range(len(sequence))] if kind=='insertion' else []
    right_o=left_o+step*(len(sequence)+1 if kind=='insertion' else 1)
    raw=['A']*30
    raw[left_o]=str(Seq(reference[left]).complement()) if reverse else reference[left]
    raw[right_o]=str(Seq(reference[right]).complement()) if reverse else reference[right]
    variants=[]
    for i,b in enumerate(sequence):
        o=insertion_positions[i] if kind=='insertion' else None
        if o is not None:raw[o]=str(Seq(b).complement()) if reverse else b
        variants.append({'kind':kind,'position':boundary if kind=='insertion' else (boundary+i)%L,'reference':'' if kind=='insertion' else b,'read':b if kind=='insertion' else '', 'original_read_position':o,'phred':35 if o is not None else None})
    row={'id':'r','alignment_id':'a','run_number':1,'original_filename':'synthetic.ab1','base_sequence':''.join(raw),'qualities':[35]*30,'report':{'evidence':{'direction':'reverse' if reverse else 'forward'},'variants':variants}}
    columns=[(left,left_o,reference[left],reference[left],35),(right,right_o,reference[right],reference[right],35)]
    return row,columns

@pytest.mark.parametrize('kind',['insertion','deletion'])
def test_repeat_shift_withheld_and_independent_haplotype_equivalence(kind):
    reference='CGTAAAATGC';b=4;seq='AA'
    row,columns=manual_event(reference,kind,b,seq)
    event=review_indels(row,reference,'linear',columns)[0]
    assert 'repeat_shift_possible' in event['review_reasons']
    assert event['equivalent_one_base_shift']=={'left':True,'right':True}
    def hap(pos,bases):return reference[:pos]+(bases if kind=='insertion' else '')+reference[pos+(0 if kind=='insertion' else len(bases)):]
    assert hap(b,seq)==hap(b-1,seq[-1]+seq[:-1])==hap(b+1,seq[1:]+seq[0])

@pytest.mark.parametrize('reverse',[False,True])
def test_circular_origin_deletion_and_insertion_zero_boundary(reverse):
    reference='ACGTCGATGC'
    row,columns=manual_event(reference,'deletion',9,'CA',reverse,True)
    event=review_indels(row,reference,'circular',columns)[0]
    assert event['reference_segments']==[[9,10],[0,1]]
    assert event['flanks']['left']['reference_position']==8 and event['flanks']['right']['reference_position']==1
    row,columns=manual_event(reference,'insertion',0,'TG',reverse,True)
    event=review_indels(row,reference,'circular',columns)[0]
    assert event['boundary']==0 and event['flanks']['left']['reference_position']==9
    assert event['flanks']['right']['reference_position']==0


def test_linear_terminal_boundary_has_no_fabricated_other_flank():
    row,columns=manual_event('ACGTCGATGC','insertion',0,'TG')
    event=review_indels(row,'ACGTCGATGC','linear',columns)[0]
    assert event['flanks']['left'] is None and 'missing_left_flank' in event['review_reasons']
    row,columns=manual_event('ACGTCGATGC','insertion',10,'TG')
    event=review_indels(row,'ACGTCGATGC','linear',columns)[0]
    assert event['flanks']['right'] is None and 'missing_right_flank' in event['review_reasons']

@pytest.mark.parametrize('kind',['insertion','deletion'])
def test_duplicate_saved_event_or_count_mismatch_is_withheld(kind):
    row,columns,start=fixture(kind)
    row['report']['variants']*=2
    events=review_indels(row,REF,'linear',columns)
    assert events and all(not e['eligible_for_exact_anchor_review'] for e in events)
    assert all('inconsistent_event_totals_or_duplicates' in e['review_reasons'] for e in events)
    row,columns,start=fixture(kind)
    row['report']['inserted_bases' if kind=='insertion' else 'deleted_bases']=99
    assert not review_indels(row,REF,'linear',columns)[0]['eligible_for_exact_anchor_review']
