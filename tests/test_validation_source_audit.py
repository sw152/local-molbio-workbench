"""Inventory cannot promote missing traces, uncertain types or a signature into truth."""
import importlib.util
import io
import json
from pathlib import Path
import zipfile
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import pytest

spec = importlib.util.spec_from_file_location('audit_source', Path(__file__).resolve().parents[1] / 'scripts' / 'audit_validation_source.py')
audit_source = importlib.util.module_from_spec(spec)
spec.loader.exec_module(audit_source)


def archive(tmp_path, entries):
    path = tmp_path / 'source.zip'
    with zipfile.ZipFile(path, 'w') as z:
        for name, data in entries:
            z.writestr(name, data)
    return path


def genbank(molecule='DNA'):
    r = SeqRecord(Seq('ACGT' * 20), id='secret-research-name', description='private annotation', annotations={'molecule_type':molecule, 'topology':'circular'})
    out = io.StringIO()
    SeqIO.write(r, out, 'genbank')
    return out.getvalue()


def test_counts_records_not_filenames_and_never_claims_ready(tmp_path):
    path = archive(tmp_path, [('private.gb', genbank() * 2), ('protein.gp', genbank('protein'))])
    before = path.read_bytes()
    result = audit_source.audit(path)
    assert result['genbank_records'] == 3
    assert result['molecule_types'] == {'DNA':2, 'PROTEIN':1}
    assert result['abif_magic_count'] == 0 and result['readiness'] == 'not_established'
    assert path.read_bytes() == before
    assert 'private' not in json.dumps(result) and 'secret-research-name' not in json.dumps(result)
    assert 'ACGT' not in json.dumps(result)


def test_magic_detection_and_opaque_members_do_not_hide_missing_evidence(tmp_path):
    path = archive(tmp_path, [('renamed.gb', b'ABIFonly-a-signature'), ('trace.dat', b'.scfheader'), ('not-a-trace.ab1', b'not abif'), ('nested.zip', b'PK\x05\x06header'), ('bad.gb', 'not genbank')])
    result = audit_source.audit(path)
    assert result['abif_magic_count'] == result['scf_magic_count'] == 1
    assert result['opaque_members'] == 2 and result['nested_archives'] == 1
    assert result['empty_genbank_members'] == 1
    assert result['readiness'] == 'not_established'
    assert result['pairing_status'].startswith('not_assessed')


def test_duplicate_members_are_counted_separately_without_extracting(tmp_path):
    with pytest.warns(UserWarning):
        path = archive(tmp_path, [('../escape.gb', genbank()), ('../escape.gb', genbank())])
    result = audit_source.audit(path)
    assert result['duplicate_names'] == 1 and result['genbank_records'] == 2
    assert not (tmp_path.parent / 'escape.gb').exists()


@pytest.mark.parametrize('field,limit,error', [('MAX_ARCHIVE',1,'archive_size_limit'), ('MAX_MEMBER',1,'member_size_limit'), ('MAX_TOTAL',1,'decoded_size_limit'), ('MAX_MEMBERS',0,'member_count_limit')])
def test_limits_reject_before_unbounded_parse(tmp_path, monkeypatch, field, limit, error):
    path = archive(tmp_path, [('one.gb', genbank())])
    monkeypatch.setattr(audit_source, field, limit)
    with pytest.raises(ValueError, match=error): audit_source.audit(path)


def test_unknown_molecule_and_broken_text_are_not_reference_confirmation(tmp_path):
    path = archive(tmp_path, [('unknown.gb', genbank().replace('DNA', '   ')), ('broken.gp', b'\xff\xfe')])
    result = audit_source.audit(path)
    assert 'private_molecule' not in json.dumps(result)
    assert result['molecule_types'] == {'other_or_missing':1}
    assert result['parse_failures'] == 1


def test_double_stranded_dna_is_not_misclassified_as_unknown(tmp_path):
    result = audit_source.audit(archive(tmp_path, [('one.gb', genbank('ds-DNA'))]))
    assert result['molecule_types'] == {'DS-DNA':1}
