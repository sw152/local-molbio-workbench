import pytest
from abif_fixture import synthetic_ab1
from localmolbio.sanger import parse_ab1, SangerReadError


def test_binary_abif_retains_base_quality(tmp_path):
    path = tmp_path / "synthetic.ab1"
    sequence = "ACGTNACGTACGTACG"
    qualities = [5, 10, 20, 30] * 4
    path.write_bytes(synthetic_ab1(sequence, qualities))
    result = parse_ab1(path)
    assert result.sequence == sequence
    assert result.qualities == qualities
    assert result.quality_summary["q20_fraction"] == 0.5
    assert "synthetic" in result.parser_metadata["SMPL1"]


def test_invalid_binary_is_rejected(tmp_path):
    path = tmp_path / "invalid.ab1"
    path.write_bytes(b"not an ABIF file")
    with pytest.raises(SangerReadError, match="Could not parse"):
        parse_ab1(path)
