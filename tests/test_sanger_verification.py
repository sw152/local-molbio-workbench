from localmolbio.sanger_verification import align_sanger_read


def test_aligns_reverse_read_and_reports_variant() -> None:
    reference = "TTTTACGTCAGGCTAACCGTACGATCGGATCGTTAACCGGTTAA"
    forward_read = "ACGTCAGGCTAACCATACGATCGGATC"
    reverse_read = "GATCCGATCGTATGGTTAGCCTGACGT"

    result = align_sanger_read(reference, reverse_read, "reverse", circular=False)

    assert result.reference_start == 4
    assert result.mismatched_bases == 1
    assert result.identity_fraction > 0.95
    assert result.variants[0]["kind"] == "substitution"


import random
import pytest
from Bio.Seq import Seq
from localmolbio.sanger_verification import SangerVerificationError

_rng = random.Random(127)
REFERENCE = "".join(_rng.choice("ACGT") for _ in range(240))


def test_linear_terminal_coordinate_is_length_not_zero():
    result = align_sanger_read(REFERENCE, REFERENCE[-80:], "forward", False)
    assert (result.reference_start, result.reference_end) == (160, 240)
    assert result.evidence["reference_covered_intervals"] == [[160, 240]]
    assert not result.wraps_origin


def test_full_linear_read():
    result = align_sanger_read(REFERENCE, REFERENCE, "forward", False)
    assert result.reference_end == len(REFERENCE)
    assert result.identity_fraction == 1
    assert result.evidence["whole_reference_verified"] is False
    assert result.evidence["review_flags"] == ["quality_unavailable"]


def test_circular_origin_and_duplicate_second_copy():
    result = align_sanger_read(REFERENCE, REFERENCE[-40:] + REFERENCE[:40], "forward", True)
    assert (result.reference_start, result.reference_end, result.wraps_origin) == (200, 40, True)
    assert result.evidence["reference_covered_intervals"] == [[0, 40], [200, 240]]
    ordinary = align_sanger_read(REFERENCE, REFERENCE[50:130], "forward", True)
    assert "ambiguous_alignment" not in ordinary.evidence["review_flags"]


def test_unknown_direction_reverses_quality_positions():
    query = list(REFERENCE[40:160])
    query[30] = next(base for base in "ACGT" if base != query[30])
    reverse = str(Seq("".join(query)).reverse_complement())
    qualities = [40] * len(reverse)
    qualities[len(reverse) - 1 - 30] = 5
    result = align_sanger_read(REFERENCE, reverse, "unknown", False, qualities)
    assert result.evidence["direction"] == "reverse"
    assert result.variants[0]["position"] == 70
    assert result.variants[0]["phred"] == 5
    assert result.variants[0]["original_read_position"] == 89
    assert "low_quality_aligned_bases" in result.evidence["review_flags"]


@pytest.mark.parametrize("kind", ["insertion", "deletion"])
def test_indel_identity_includes_gap_and_insertion_has_boundary(kind):
    query = REFERENCE[30:130] + ("AAA" if kind == "insertion" else "") + REFERENCE[130 if kind == "insertion" else 133:220]
    result = align_sanger_read(REFERENCE, query, "forward", False)
    gaps = result.inserted_bases + result.deleted_bases
    assert gaps == 3
    assert result.identity_fraction == round(result.matched_bases / (result.aligned_bases + gaps), 6)
    assert all(v["position"] is not None for v in result.variants)
    assert result.variants[0]["kind"] == kind
    if kind == "deletion":
        assert result.evidence["reference_covered_intervals"] == [[30, 130], [133, 220]]


def test_ambiguous_base_is_not_matching_evidence():
    template = REFERENCE[:80] + "N" + REFERENCE[81:]
    result = align_sanger_read(template, template, "forward", False)
    assert result.matched_bases == len(template) - 1
    assert result.evidence["ambiguous_bases"] == 1
    assert result.identity_fraction < 1
    assert result.variants[0]["kind"] == "ambiguous"


def test_repeated_mapping_and_partial_match_are_flagged():
    motif = REFERENCE[20:60]
    repeated = align_sanger_read(motif + "A" * 15 + motif, motif, "forward", False)
    assert "ambiguous_alignment" in repeated.evidence["review_flags"]
    partial = align_sanger_read(REFERENCE, "N" * 100 + REFERENCE[50:90] + "N" * 100, "forward", False)
    assert "partial_read_alignment" in partial.evidence["review_flags"]


@pytest.mark.parametrize("read,direction,qualities", [("N" * 50, "forward", None), (REFERENCE, "sideways", None), (REFERENCE, "forward", [40])])
def test_unusable_input_rejected(read, direction, qualities):
    with pytest.raises(SangerVerificationError):
        align_sanger_read(REFERENCE, read, direction, False, qualities)


def test_interactive_size_limit():
    with pytest.raises(SangerVerificationError, match="limit"):
        align_sanger_read("A" * 100_000, "A" * 1000, "forward", False)
