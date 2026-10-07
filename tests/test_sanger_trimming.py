import random
import pytest
from Bio.Seq import Seq
from localmolbio.sanger_verification import align_sanger_read, SangerVerificationError

rng=random.Random(234)
REFERENCE=''.join(rng.choice('ACGT') for _ in range(400))


@pytest.mark.parametrize('reverse',[False,True])
def test_trim_preserves_full_original_coordinates_and_internal_low_quality(reverse):
    read=list(REFERENCE[40:240]);read[80]=next(b for b in 'ACGT' if b!=read[80]);read=''.join(read)
    qualities=[8]*13+[38]*166+[9]*21
    qualities[13]=qualities[178]=20  # threshold-equal endpoints are retained
    qualities[80]=12  # interior low-quality difference must not disappear
    if reverse:read=str(Seq(read).reverse_complement());qualities=qualities[::-1]
    original_qualities=qualities.copy()
    before=align_sanger_read(REFERENCE,read,'unknown',False,qualities)
    after=align_sanger_read(REFERENCE,read,'unknown',False,qualities,20)
    assert after.evidence['direction']==('reverse' if reverse else 'forward')
    assert (after.reference_start,after.reference_end)==(53,219)
    assert after.variants[0]['original_read_position']==(119 if reverse else 80)
    assert after.variants[0]['original_read_position']==before.variants[0]['original_read_position']
    assert after.variants[0]['read_position']==80
    assert after.variants[0]['analysis_read_position']==67
    assert after.variants[0]['phred']==12
    assert after.evidence['low_quality_aligned_bases']==1
    assert after.evidence['read_aligned_fraction']==.83
    assert after.evidence['retained_read_aligned_fraction']==1
    assert (after.evidence['read_start'],after.evidence['read_end'])==(13,179)
    assert (after.evidence['original_read_aligned_start'],after.evidence['original_read_aligned_end'])==((21,187) if reverse else (13,179))
    trim=after.evidence['trimming']
    assert trim['original_length']==200 and trim['retained_length']==166
    assert (trim['removed_left'],trim['removed_right'])==((21,13) if reverse else (13,21))
    assert qualities==original_qualities


def test_disabled_trim_keeps_all_bases_and_low_quality_flags():
    result=align_sanger_read(REFERENCE,REFERENCE,'forward',False,[0]*400)
    assert result.reference_end==400 and result.evidence['trimming']['method']=='none'
    assert result.evidence['read_aligned_fraction']==result.evidence['retained_read_aligned_fraction']==1
    assert result.evidence['low_quality_aligned_bases']==400


def test_single_threshold_passing_call_is_not_a_usable_read():
    with pytest.raises(SangerVerificationError,match='fewer than 12'):
        align_sanger_read(REFERENCE,REFERENCE,'forward',False,[0]*200+[20]+[0]*199,20)


@pytest.mark.parametrize('threshold',[0,61,True,20.5,'20'])
def test_invalid_threshold_rejected(threshold):
    with pytest.raises(SangerVerificationError,match='integer'):
        align_sanger_read(REFERENCE,REFERENCE,'forward',False,[40]*400,threshold)


@pytest.mark.parametrize('qualities',[None,[40]*399,[0]*400])
def test_missing_inconsistent_or_all_low_quality_rejected(qualities):
    with pytest.raises(SangerVerificationError):
        align_sanger_read(REFERENCE,REFERENCE,'forward',False,qualities,20)


def test_trimmed_circular_origin_keeps_reference_coverage_and_original_fraction():
    read=REFERENCE[-60:]+REFERENCE[:100]
    result=align_sanger_read(REFERENCE,read,'forward',True,[5]*5+[35]*148+[5]*7,20)
    assert (result.reference_start,result.reference_end,result.wraps_origin)==(345,93,True)
    assert result.evidence['reference_covered_intervals']==[[0,93],[345,400]]
    assert result.evidence['read_aligned_fraction']==.925
    assert result.evidence['retained_read_aligned_fraction']==1


@pytest.mark.parametrize('kind',['insertion','deletion'])
@pytest.mark.parametrize('reverse',[False,True])
def test_trimmed_indels_preserve_anchors_and_original_read_positions(kind,reverse):
    read=REFERENCE[30:130]+('AAA' if kind=='insertion' else '')+REFERENCE[130 if kind=='insertion' else 133:250]
    quality=[3]*12+[35]*(len(read)-29)+[3]*17
    if reverse:read=str(Seq(read).reverse_complement());quality=quality[::-1]
    before=align_sanger_read(REFERENCE,read,'unknown',False,quality)
    after=align_sanger_read(REFERENCE,read,'unknown',False,quality,20)
    assert [(v['kind'],v['position'],v['original_read_position']) for v in after.variants]==[(v['kind'],v['position'],v['original_read_position']) for v in before.variants]
    assert len(after.variants)==3
    if kind=='deletion':assert all(v['phred'] is None and v['read_position'] is None for v in after.variants)


def test_partial_retained_alignment_does_not_inflate_original_read_coverage():
    read='N'*50+REFERENCE[70:130]+'N'*50
    quality=[3]*20+[35]*120+[3]*20
    result=align_sanger_read(REFERENCE,read,'unknown',False,quality,20)
    assert result.evidence['read_aligned_fraction']==.375
    assert result.evidence['retained_read_aligned_fraction']==.5
    assert 'partial_read_alignment' in result.evidence['review_flags']
    assert result.evidence['whole_reference_verified'] is False


def test_exactly_twelve_retained_bases_are_accepted_without_whole_read_claim():
    read=REFERENCE[50:150]
    result=align_sanger_read(REFERENCE,read,'forward',False,[0]*40+[20]*12+[0]*48,20)
    assert result.evidence['trimming']['retained_length']==12
    assert result.evidence['read_aligned_fraction']==.12
    assert result.evidence['retained_read_aligned_fraction']==1
    assert 'partial_read_alignment' in result.evidence['review_flags']
    assert result.evidence['whole_reference_verified'] is False


def test_internal_ambiguous_call_is_not_removed_by_end_trimming():
    read=REFERENCE[30:100]+'N'+REFERENCE[101:170]
    qualities=[3]*10+[35]*120+[3]*10
    qualities[70]=1
    result=align_sanger_read(REFERENCE,read,'forward',False,qualities,20)
    assert result.evidence['ambiguous_bases']==1
    assert result.variants[0]['kind']=='ambiguous'
    assert result.variants[0]['original_read_position']==70
    assert result.variants[0]['phred']==1
