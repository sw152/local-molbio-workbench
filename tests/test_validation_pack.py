"""Independent decoding and construction checks for the public rehearsal inputs."""
import hashlib
import importlib.util
import json
from pathlib import Path
import zipfile

from Bio import SeqIO
from Bio.Seq import Seq
import pytest
from localmolbio.sanger_verification import align_sanger_read

spec = importlib.util.spec_from_file_location('validation_pack', Path(__file__).resolve().parents[1] / 'scripts' / 'build_validation_pack.py')
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)
FROZEN_MANIFEST = 'e9e44b1ffa1107850a3a2453cbcf5dc68131b2c3379adeef15859505a8afd9a6'


@pytest.fixture(scope='module')
def pack(tmp_path_factory):
    path = tmp_path_factory.mktemp('validation') / 'pack'
    module.build(path)
    return path


def read(pack, case, index=1):
    return SeqIO.read(pack / 'participant' / 'inputs' / f'{case}-r{index}.ab1', 'abi')


def reference(pack, case):
    return str(SeqIO.read(pack / 'participant' / 'inputs' / f'{case}.gb', 'genbank').seq)


def test_frozen_hashes_reproducibility_and_no_overwrite(pack, tmp_path):
    manifest = json.loads((pack / 'manifest.json').read_text())
    assert hashlib.sha256((pack / 'manifest.json').read_bytes()).hexdigest() == FROZEN_MANIFEST
    assert (pack / 'manifest.sha256').read_text().strip() == FROZEN_MANIFEST
    for item in manifest['files']:
        data = (pack / item['path']).read_bytes()
        assert len(data) == item['bytes']
        assert hashlib.sha256(data).hexdigest() == item['sha256']
    module.build(tmp_path / 'second')
    assert (tmp_path / 'second' / 'manifest.json').read_bytes() == (pack / 'manifest.json').read_bytes()
    with pytest.raises(FileExistsError):
        module.build(pack)
    assert (pack / 'manifest.sha256').read_text().strip() == FROZEN_MANIFEST


def test_decoded_input_truth_without_application_aligner(pack):
    ref = reference(pack, 'A01')
    assert str(read(pack, 'A01').seq) == ref[80:580]
    reverse = read(pack, 'A02')
    forward = str(reverse.seq.reverse_complement())
    assert [i for i, (a, b) in enumerate(zip(forward, ref[200:760])) if a != b] == [200]
    assert reverse.letter_annotations['phred_quality'][359] == 38
    low = read(pack, 'A03')
    assert str(low.seq[:200]) == ref[250:450]
    assert low.letter_annotations['phred_quality'][200] == 9
    assert low.letter_annotations['phred_quality'][:20] == [5] * 20
    a, b = read(pack, 'A04'), read(pack, 'A04', 2)
    assert [i for i in range(400) if a.seq[i] != b.seq[i]] == [220]
    repeat = reference(pack, 'A05')
    assert str(read(pack, 'A05').seq) == repeat[100:280] == repeat[700:880]
    assert str(read(pack, 'A06').seq) == ref[1020:] + ref[:220]
    assert str(read(pack, 'A07').seq) == ref[500:580]
    assert str(read(pack, 'A08').seq) == ref[300:800]
    assert reference(pack, 'A08') != ref
    with pytest.raises((ValueError, OSError)):
        read(pack, 'A09')
    insertion = str(read(pack, 'A10').seq.reverse_complement())
    assert insertion[:250] + insertion[253:] == ref[400:900]
    assert len(insertion) == 503


def test_participant_materials_and_unmeasured_records(pack):
    tasks = json.loads((pack / 'participant' / 'tasks.json').read_text())
    assert len(tasks['cases']) == 10
    assert len({c['case_id'] for c in tasks['cases']}) == 10
    for case in tasks['cases']:
        assert 'expected_review' not in case
        for relative in [case['reference'], *case['reads']]:
            assert (pack / 'participant' / relative).is_file()
    with zipfile.ZipFile(pack / 'participant' / 'reference-library.zip') as archive:
        assert len(archive.namelist()) == 10
        for name in archive.namelist():
            assert archive.read(name) == (pack / 'participant' / 'inputs' / name).read_bytes()
    session = json.loads((pack / 'session-template.json').read_text())
    assert session['observations'] == session['retention'] == []
    assert session['record_type'] is session['d0'] is session['setup_seconds'] is None
    assert all(v is None for v in json.loads((pack / 'observation-template.json').read_text()).values())
    assert all(v is None for v in json.loads((pack / 'retention-template.json').read_text()).values())
    ledger = json.loads((pack / 'evaluator' / 'construction-ledger.json').read_text())
    assert not ledger['independently_reviewed']


@pytest.mark.parametrize('case', ['A01', 'A02', 'A03', 'A04', 'A05', 'A06', 'A07', 'A08', 'A10'])
def test_sanger_scientific_boundaries_on_frozen_cases(pack, case):
    record = read(pack, case)
    result = align_sanger_read(reference(pack, case), str(record.seq), 'unknown', True, record.letter_annotations['phred_quality'])
    assert result.evidence['whole_reference_verified'] is False
    if case in ('A01', 'A06', 'A07'):
        expected = {'A01': [[80, 580]], 'A06': [[0, 220], [1020, 1200]], 'A07': [[500, 580]]}[case]
        assert result.evidence['reference_covered_intervals'] == expected
    if case == 'A02':
        assert result.evidence['direction'] == 'reverse'
        assert [(v['position'], v['original_read_position'], v['phred']) for v in result.variants] == [(400, 359, 38)]
    if case == 'A03':
        assert result.variants[0]['position'] == 450
        assert result.variants[0]['phred'] == 9
        assert 'low_quality_aligned_bases' in result.evidence['review_flags']
    if case == 'A05':
        assert 'ambiguous_alignment' in result.evidence['review_flags']
    if case == 'A08':
        assert 'partial_read_alignment' in result.evidence['review_flags']
    if case == 'A10':
        assert result.evidence['direction'] == 'reverse'
        assert result.inserted_bases == 3
        assert [v['position'] for v in result.variants] == [650, 650, 650]
        assert [v['original_read_position'] for v in result.variants] == [252, 251, 250]
