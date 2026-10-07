import pytest
from abif_fixture import synthetic_ab1
from localmolbio.chromatogram import trace_window
from localmolbio.sanger import SangerReadError

SEQUENCE = "ACGT" * 8


@pytest.mark.parametrize("order", ["GATC", "TCAG", "ACGT"])
def test_channel_order_and_peak_quality_correspondence(tmp_path, order):
    path = tmp_path / "trace.ab1"
    quality = [38] * len(SEQUENCE)
    quality[7] = 9
    path.write_bytes(synthetic_ab1(SEQUENCE, quality, trace=True, order=order))
    trace = trace_window(path, start=5, count=8)
    assert trace["available"]
    assert trace["source_channel_order"] == order
    assert (trace["base_start"], trace["base_end"]) == (5,13)
    assert trace["orientation"] == "original_uploaded_read"
    assert trace["bases"][2]["phred"] == 9
    for base in trace["bases"]:
        peak_index = base["sample"] - trace["sample_start"]
        assert trace["channels"][base["base"]][peak_index] == 900
        assert all(trace["channels"][other][peak_index] == 40 for other in "ACGT" if other != base["base"])
    last = trace_window(path, start=31, count=24)
    assert len(last["bases"]) == 1 and last["base_end"] == 32


def test_missing_peak_data_is_unavailable_not_fabricated(tmp_path):
    path = tmp_path / "basecalls-only.ab1"
    path.write_bytes(synthetic_ab1(SEQUENCE, [30]*32))
    result = trace_window(path)
    assert result["available"] is False
    assert "PLOC2" in result["reason"]


@pytest.mark.parametrize("overrides,message", [
    ({"order":"AAAA"}, "dye order"),
    ({"peak_override":[12]*32}, "must increase"),
    ({"peak_override":list(range(31))+[32000]}, "within"),
    ({"peak_override":[12,24]}, "counts"),
])
def test_inconsistent_trace_metadata_rejected(tmp_path, overrides, message):
    path = tmp_path / "invalid.ab1"
    path.write_bytes(synthetic_ab1(SEQUENCE,[30]*32,trace=True,**overrides))
    with pytest.raises(SangerReadError, match=message):
        trace_window(path)


@pytest.mark.parametrize("start,count", [(-1,24),(0,81),(32,24),(0,0)])
def test_window_bounds(tmp_path,start,count):
    path=tmp_path/'trace.ab1'
    path.write_bytes(synthetic_ab1(SEQUENCE,[30]*32,trace=True))
    with pytest.raises(SangerReadError):
        trace_window(path,start,count)
