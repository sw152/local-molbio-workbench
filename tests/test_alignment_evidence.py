"""Independent base/coordinate checks and optional real pinned-minimap2 integration."""
import json
import os
from pathlib import Path
import random
import sys

import pytest
from localmolbio.alignment_evidence import AlignmentError, align, parse_paf


def line(qlen,qs,qe,strand,tlen,ts,te,matches,block,cigar,mapq=60,tp='P',nn=0):
    return f'q0\t{qlen}\t{qs}\t{qe}\t{strand}\treference\t{tlen}\t{ts}\t{te}\t{matches}\t{block}\t{mapq}\ttp:A:{tp}\tnn:i:{nn}\tcg:Z:{cigar}\n'


def rc(s):return s.translate(str.maketrans('ACGTN','TGCAN'))[::-1]


def test_linear_endpoint_reverse_original_coordinates_and_unknown_mapq():
    ref='ACGTTGCA';read=rc(ref[4:])+'NN'
    result=parse_paf([line(6,0,4,'-',8,4,8,4,4,'4=',mapq=255)],ref,[read],'linear')[0]
    hit=result['alignments'][0]
    assert hit['reference_intervals']==[[4,8]] and not hit['crosses_origin']
    assert hit['paired_blocks']==[{'query_start':0,'query_end':4,'query_step':-1,'reference_intervals':[[4,8]]}]
    assert hit['raw_hits'][0]['mapq'] is None
    assert hit['aligned_query_fraction']==pytest.approx(4/6) and not result['uniqueness_assessed']


def test_circular_origin_and_equivalent_copy_deduplication():
    ref='ACGTTGCA'
    hits=parse_paf([line(4,0,4,'+',16,6,10,4,4,'4=')],ref,['CAAC'],'circular')[0]['alignments']
    assert hits[0]['reference_intervals']==[[6,8],[0,2]] and hits[0]['crosses_origin']
    copies=[line(4,0,4,'+',16,1,5,4,4,'4='),line(4,0,4,'+',16,9,13,4,4,'4=',mapq=0,tp='S')]
    result=parse_paf(copies,ref,[ref[1:5]],'circular')[0]
    assert len(result['alignments'])==1 and len(result['alignments'][0]['raw_hits'])==2
    assert result['status']=='alignment_reported' and not result['uniqueness_assessed']


def test_insertions_deletions_mismatch_and_ambiguous_bases_count_in_denominator():
    # AC / +G / -T / A->C / N~N: 2 exact, 1 mismatch, 1 ambiguous, 2 gap columns.
    ref='ACTAN';read='ACGCN'
    hit=parse_paf([line(5,0,5,'+',5,0,5,2,5,'2=1I1D1X1=',nn=1)],ref,[read],'linear')[0]['alignments'][0]
    assert (hit['exact_matches'],hit['mismatches'],hit['ambiguous_pairs'],hit['inserted_bases'],hit['deleted_bases'])==(2,1,1,1,1)
    assert hit['local_exact_identity']==pytest.approx(2/6)
    assert [b['reference_intervals'] for b in hit['paired_blocks']]==[[[0,2]],[[3,4]],[[4,5]]]


@pytest.mark.parametrize('bad',[
    line(4,0,5,'+',4,0,4,4,4,'4='),line(4,0,4,'?',4,0,4,4,4,'4='),
    line(4,0,4,'+',4,0,5,4,4,'4='),line(4,0,4,'+',4,0,4,4,4,'4M'),
    line(4,0,4,'+',4,0,4,4,4,'3='),line(4,0,4,'+',4,0,4,4,4,'5='),
    line(4,0,4,'+',4,0,4,4,5,'4='),line(4,0,4,'+',4,0,4,3,4,'4='),
    line(4,0,4,'+',4,0,4,4,4,'4X'),line(4,0,4,'+',4,0,4,4,4,'4=',tp='I'),
    line(4,0,4,'+',4,0,4,4,4,'4=',mapq=256),
    line(4,0,4,'+',4,0,4,4,4,'4=').replace('q0','q7'),
    line(4,0,4,'+',4,0,4,4,4,'4=').strip()+'\tcg:Z:4=\n',
    line(4,0,4,'+',4,0,4,4,4,'4=').replace('\treference\t','\tother\t'),
])
def test_inconsistent_or_approximate_paf_rejected(bad):
    with pytest.raises(AlignmentError):parse_paf([bad],'ACGT',['ACGT'],'linear')


def test_no_reported_alignment_and_multitraversal_are_not_passes():
    empty=parse_paf([],'ACGT',['ACGT'],'linear')[0]
    assert empty['status']=='no_alignment_reported' and not empty['uniqueness_assessed']
    withheld=parse_paf([line(8,0,8,'+',8,0,8,8,8,'8=')],'ACGT',['ACGTACGT'],'circular')[0]
    assert not withheld['alignments'] and withheld['withheld'][0]['reason']=='more_than_one_reference_traversal'
    assert withheld['requires_review']


@pytest.mark.parametrize('reference,reads,topology',[
    ('ACGT',['ACGT'],'unknown'),('ACGR',['ACGT'],'linear'),('ACGT',[],'linear'),
    ('ACGT',['A'*100001],'linear'),('A'*200001,['ACGT'],'linear'),
])
def test_input_bounds_are_explicit(reference,reads,topology):
    with pytest.raises(AlignmentError):parse_paf([],reference,reads,topology)


@pytest.fixture
def minimap2():
    binary=os.environ.get('MOLBIO_MINIMAP2')
    if not binary:pytest.skip('Set MOLBIO_MINIMAP2 to pinned 2.31-r1302 for actual aligner tests')
    return Path(binary).resolve(strict=True)


@pytest.fixture
def reference():return ''.join(random.Random(2718).choices('ACGT',k=4000))


@pytest.mark.parametrize('data_type',['ont-noisy','ont-high-accuracy','pacbio-hifi','short-single'])
def test_actual_tool_forward_reverse_partial_and_unreported(minimap2,tmp_path,reference,data_type):
    n=150 if data_type=='short-single' else 800
    reads=[reference[200:200+n],rc(reference[-n:]),'N'*n]
    result=align(reference,reads,topology='linear',data_type=data_type,binary=minimap2,output_dir=tmp_path/'run')
    first=result['reads'][0]['alignments'][0];second=result['reads'][1]['alignments'][0]
    assert first['reference_intervals']==[[200,200+n]] and first['exact_matches']==n
    assert second['reference_intervals']==[[len(reference)-n,len(reference)]] and second['strand']=='-'
    assert second['paired_blocks'][0]['query_start']==0 and second['paired_blocks'][0]['query_end']==n
    assert result['reads'][2]['status']=='no_alignment_reported'
    assert result['whole_reference_verified'] is result['consensus_performed'] is result['base_quality_used'] is False
    assert result['tool']['version']=='2.31-r1302' and len(result['tool']['binary_sha256'])==64
    assert (tmp_path/'run'/'alignments.paf').stat().st_size>0
    assert json.loads((tmp_path/'run'/'result.json').read_text())==result
    with pytest.raises(FileExistsError):align(reference,reads,topology='linear',data_type=data_type,binary=minimap2,output_dir=tmp_path/'run')


def test_actual_circular_crossing_copies_and_repeated_locations(minimap2,tmp_path,reference):
    read=reference[-400:]+reference[:400]
    result=align(reference,[read,rc(read),reference[100:900]],topology='circular',data_type='ont-high-accuracy',binary=minimap2,output_dir=tmp_path/'circular')
    for i in (0,1):
        hit=next(h for h in result['reads'][i]['alignments'] if h['exact_matches']==800)
        assert hit['reference_intervals']==[[3600,4000],[0,400]] and hit['crosses_origin']
        assert hit['strand']==('+' if i==0 else '-')
    interior=result['reads'][2]['alignments']
    assert len(interior)==1 and len(interior[0]['raw_hits'])==2
    repeat=reference[:1000]+reference[2000:3000]+reference[:1000]
    multi=align(repeat,[reference[100:900]],topology='linear',data_type='ont-high-accuracy',binary=minimap2,output_dir=tmp_path/'repeat')['reads'][0]
    assert multi['status']=='multiple_alignments_reported' and multi['requires_review']
    assert {tuple(h['reference_intervals'][0]) for h in multi['alignments']}=={(100,900),(2100,2900)}


def test_actual_mismatch_gaps_and_local_overhang(minimap2,tmp_path,reference):
    read=reference[500:1500]
    altered=read[:300]+('A' if read[300]!='A' else 'C')+read[301:500]+'TGCA'+read[500:700]+read[705:]
    flank=''.join(random.Random(812).choices('ACGT',k=200))
    result=align(reference,[altered,flank+read+flank,read[:450]+'N'+read[451:]],topology='linear',data_type='ont-high-accuracy',binary=minimap2,output_dir=tmp_path/'changes')
    hit=result['reads'][0]['alignments'][0]
    assert hit['mismatches']>=1 and hit['inserted_bases']>=4 and hit['deleted_bases']>=5
    assert hit['local_exact_identity']<1
    local=result['reads'][1]['alignments'][0]
    assert local['aligned_query_fraction']<.8 and not result['whole_reference_verified']
    ambiguous=result['reads'][2]['alignments'][0]
    assert ambiguous['ambiguous_pairs']>=1 and ambiguous['exact_matches']<1000


def fake_binary(tmp_path,body,version='2.31-r1302'):
    binary=tmp_path/'tool'
    binary.write_text(f'#!{sys.executable}\nimport sys,time\nif "--version" in sys.argv:\n print({version!r});sys.exit(0)\n'+body)
    binary.chmod(0o700)
    return binary


def test_process_timeout_and_cancellation_leave_no_result(tmp_path):
    binary=fake_binary(tmp_path,'time.sleep(10)\n')
    with pytest.raises(AlignmentError,match='alignment_timeout'):
        align('ACGT',['ACGT'],topology='linear',data_type='short-single',binary=binary,output_dir=tmp_path/'timeout',timeout_seconds=1)
    assert not (tmp_path/'timeout'/'result.json').exists()
    calls=[0]
    def pulse():
        calls[0]+=1
        if calls[0]>1:raise RuntimeError('lease cancelled')
    with pytest.raises(RuntimeError,match='lease cancelled'):
        align('ACGT',['ACGT'],topology='linear',data_type='short-single',binary=binary,output_dir=tmp_path/'cancel',pulse=pulse)
    assert not (tmp_path/'cancel'/'result.json').exists()


def test_process_failure_version_and_invalid_data_type(tmp_path):
    binary=fake_binary(tmp_path,'sys.exit(7)\n')
    with pytest.raises(AlignmentError,match='process_failed'):
        align('ACGT',['ACGT'],topology='linear',data_type='short-single',binary=binary,output_dir=tmp_path/'fail')
    with pytest.raises(AlignmentError,match='data_type'):
        align('ACGT',['ACGT'],topology='linear',data_type='automatic',binary=binary,output_dir=tmp_path/'invalid')
    fake_binary(tmp_path,'sys.exit(0)\n',version='unknown')
    with pytest.raises(AlignmentError,match='unvalidated_minimap2_version'):
        align('ACGT',['ACGT'],topology='linear',data_type='short-single',binary=binary,output_dir=tmp_path/'version')


def test_ambiguous_insertions_and_deletions_remain_in_our_denominator():
    hit=parse_paf([line(6,0,6,'+',6,0,6,4,4,'2=2I2D2=',nn=4)],'ACNNGT',['ACNNGT'],'linear')[0]['alignments'][0]
    assert hit['alignment_columns']==8 and hit['reported_paf_block_length']==4
    assert hit['ambiguous_gap_bases']==4 and hit['local_exact_identity']==.5
    with pytest.raises(AlignmentError,match='ambiguous_count'):
        parse_paf([line(6,0,6,'+',6,0,6,4,4,'2=2I2D2=')],'ACNNGT',['ACNNGT'],'linear')


def test_n_is_not_a_wildcard_for_eqx_validation():
    with pytest.raises(AlignmentError,match='base_disagreement'):
        parse_paf([line(2,0,2,'+',2,0,2,1,1,'2=',nn=1)],'AC',['AN'],'linear')


def test_output_guard_rejects_oversize_artifacts(tmp_path,monkeypatch):
    import localmolbio.alignment_evidence as core
    binary=fake_binary(tmp_path,"print('too much output')\n")
    monkeypatch.setattr(core,'MAX_OUTPUT',4)
    with pytest.raises(AlignmentError,match='output_limit'):
        align('ACGT',['ACGT'],topology='linear',data_type='short-single',binary=binary,output_dir=tmp_path/'oversized')
    assert not (tmp_path/'oversized'/'result.json').exists()


def test_changed_materialized_input_cannot_be_reported_as_no_hits(tmp_path):
    binary=fake_binary(tmp_path,"from pathlib import Path\nPath('reads.fa').write_text('>q0\\nTTTT\\n')\n")
    with pytest.raises(AlignmentError,match='inputs_changed'):
        align('ACGT',['ACGT'],topology='linear',data_type='short-single',binary=binary,output_dir=tmp_path/'changed')
    assert not (tmp_path/'changed'/'result.json').exists()
