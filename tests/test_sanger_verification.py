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
