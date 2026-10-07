import pytest
from localmolbio.primer_sites import exact_reference_sites


def test_overlapping_repeats_and_exact_total_when_truncated():
    result = exact_reference_sites('AAAAAA', 'AAA', 'linear', limit=2)
    assert result['total_directional_matches'] == 4
    assert [s['start'] for s in result['sites']] == [0, 1]
    assert result['truncated']


def test_reverse_matches_and_terminal_boundary():
    r = exact_reference_sites('CCCAAGTGGACTT', 'AAGT', 'linear')
    assert [(s['start'], s['end'], s['strand']) for s in r['sites']] == [(3, 7, 1), (9, 13, -1)]
    assert all(not s['wraps_origin'] for s in r['sites'])


def test_origin_wrap_and_unknown_topology():
    r = exact_reference_sites('GTCCAA', 'AAGT', 'circular')
    assert r['sites'] == [{'start':4,'end':2,'strand':1,'wraps_origin':True,'segments':[[4,6],[0,2]]}]
    for topology in ('linear', 'unknown'):
        r = exact_reference_sites('GTCCAA', 'AAGT', topology)
        assert r['total_directional_matches'] == 0
        assert not r['origin_checked']


def test_palindrome_is_two_directions_not_two_physical_intervals():
    r = exact_reference_sites('ACGT', 'ACGT', 'circular')
    assert [(s['start'],s['end'],s['strand']) for s in r['sites']] == [(0,4,1),(0,4,-1)]


@pytest.mark.parametrize('sequence,oligo,topology,limit', [('', 'A','linear',100), ('AN','A','linear',100),('AC','ACG','linear',100),('AC','A','other',100),('AC','A','linear',0)])
def test_invalid_review_inputs(sequence,oligo,topology,limit):
    with pytest.raises(ValueError):exact_reference_sites(sequence,oligo,topology,limit)


def test_sites_equal_independent_window_enumeration():
    import random
    from Bio.Seq import Seq
    rng = random.Random(71)
    for topology in ('linear', 'circular'):
        for _ in range(50):
            reference = ''.join(rng.choice('ACGT') for _ in range(20))
            oligo = ''.join(rng.choice('ACGT') for _ in range(rng.randint(1, 5)))
            expected = set()
            for pos in range(20 if topology == 'circular' else 21-len(oligo)):
                window = ''.join(reference[(pos+j) % 20] for j in range(len(oligo)))
                for strand, needle in ((1, oligo), (-1, str(Seq(oligo).reverse_complement()))):
                    if window == needle: expected.add((pos, strand))
            result = exact_reference_sites(reference, oligo, topology)
            assert {(s['start'], s['strand']) for s in result['sites']} == expected
            assert result['total_directional_matches'] == len(expected)
