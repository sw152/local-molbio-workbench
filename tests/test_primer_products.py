import random
import pytest
from Bio.Seq import Seq
from localmolbio.primer_sites import exact_reference_sites
from localmolbio.primer_products import review_pair_products


def item(plus=(),minus=(),length=10,reference_length=100,topology='linear',intended_start=0,truncated=False):
    return {'sequence':'A'*length,'binding_start':intended_start,'binding_end':intended_start+length,
            'reference_sites':{'reference_length':reference_length,'reference_topology':topology,
                               'truncated':truncated,'sites':[{'start':p,'strand':s} for s,positions in ((1,plus),(-1,minus)) for p in positions]}}


def test_extra_single_site_is_not_automatically_an_extra_product():
    left=item(plus=[10,80],intended_start=10);right=item(minus=[60],intended_start=60)
    result=review_pair_products(left,right,20,70)
    assert result['total_product_count']==1 and result['alternative_product_count']==0
    assert result['products'][0]['is_intended']
    assert result['products'][0]['length_bp']==60


def test_swapped_oligo_assignments_and_length_inclusive_boundaries():
    left=item(plus=[10],minus=[80],intended_start=10)
    right=item(plus=[30],minus=[60],intended_start=60)
    result=review_pair_products(left,right,60,60)
    assert result['total_product_count']==2 and result['alternative_product_count']==1
    assert result['products'][1]['assignments'][0]['plus_primer']=='reverse'
    assert review_pair_products(left,right,61,100)['total_product_count']==0


def test_overlap_excluded_but_adjacent_nonoverlapping_sites_count():
    result=review_pair_products(item(plus=[10]),item(minus=[15,20]),1,100)
    assert [(p['start'],p['length_bp']) for p in result['products']]==[(10,20)]


def test_linear_terminal_boundary_does_not_wrap():
    p=review_pair_products(item(plus=[70]),item(minus=[90]),1,100)['products'][0]
    assert p['end']==100 and not p['wraps_origin']


@pytest.mark.parametrize('positive,negative,span,segments',[(80,10,40,[[80,100],[0,20]]),(95,20,35,[[95,100],[0,30]]),(10,0,100,[[10,100],[0,10]])])
def test_circular_origin_and_single_full_traversal(positive,negative,span,segments):
    result=review_pair_products(item(plus=[positive],topology='circular'),item(minus=[negative],topology='circular'),1,100)
    p=result['products'][0]
    assert p['length_bp']==span and p['segments']==segments and p['wraps_origin']
    assert result['origin_checked']


def test_circular_overlap_across_origin_and_unknown_topology():
    assert review_pair_products(item(plus=[95],topology='circular'),item(minus=[98],topology='circular'),1,200)['total_product_count']==0
    result=review_pair_products(item(plus=[80],topology='unknown'),item(minus=[10],topology='unknown'),1,100)
    assert result['total_product_count']==0 and not result['origin_checked']


def test_identical_oligos_do_not_double_count_the_same_product_interval():
    a=item(plus=[10],minus=[60],intended_start=10)
    b=item(plus=[10],minus=[60],intended_start=60)
    result=review_pair_products(a,b,1,100)
    assert result['total_product_count']==1 and len(result['products'][0]['assignments'])==2
    assert result['alternative_product_count']==0


def test_truncated_sites_never_claim_total_counts_or_absence_of_alternatives():
    result=review_pair_products(item(plus=[10],truncated=True),item(minus=[60]),1,100)
    assert not result['site_search_complete']
    assert result['total_product_count'] is None and result['alternative_product_count'] is None
    assert result['observed_product_count']==1


def test_product_display_limit_preserves_complete_counts_and_intended_first():
    result=review_pair_products(item(plus=[0,10,20],intended_start=0),item(minus=[60,70,80],intended_start=80),1,100,limit=2)
    assert result['total_product_count']==9 and result['alternative_product_count']==8
    assert len(result['products'])==2 and result['product_list_truncated']
    assert result['products'][0]['is_intended']


def test_pair_geometry_matches_independent_sequence_window_oracle():
    rng=random.Random(678)
    for topology in ('linear','circular'):
        for _ in range(40):
            reference=''.join(rng.choice('ACGT') for _ in range(45))
            a=reference[rng.randrange(10):][:3];b=str(Seq(reference[30:34]).reverse_complement())
            left={'sequence':a,'binding_start':0,'binding_end':3,'reference_sites':exact_reference_sites(reference,a,topology)}
            right={'sequence':b,'binding_start':30,'binding_end':34,'reference_sites':exact_reference_sites(reference,b,topology)}
            result=review_pair_products(left,right,7,35)
            expected=set()
            for start in range(45):
                for span in range(7,36):
                    if topology=='linear' and start+span>45:continue
                    window=''.join(reference[(start+j)%45] for j in range(span))
                    if any(window.startswith(f) and window.endswith(str(Seq(r).reverse_complement())) for f,r in ((a,b),(b,a))):
                        expected.add((start,span))
            assert {(p['start'],p['length_bp']) for p in result['products']}==expected
            assert result['total_product_count']==len(expected)
