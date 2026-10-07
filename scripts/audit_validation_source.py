"""Read-only bounded ZIP inventory for trial readiness, not biological validation.

Reports counts and archive identity, never member names, sequences or annotations.
No extraction, reference matching, network request or application database writes.
"""
import argparse
from collections import Counter
import hashlib
import io
import json
from pathlib import Path
import warnings
import zipfile
from Bio import SeqIO

MAX_ARCHIVE = 128 * 1024 * 1024
MAX_TOTAL = 512 * 1024 * 1024
MAX_MEMBER = 32 * 1024 * 1024
MAX_MEMBERS = 10000


def audit(path):
    path = Path(path)
    if path.stat().st_size > MAX_ARCHIVE:
        raise ValueError('archive_size_limit')
    with path.open('rb') as stream:
        payload = stream.read(MAX_ARCHIVE + 1)
    if len(payload) > MAX_ARCHIVE:
        raise ValueError('archive_size_limit')
    report = {'schema_version': 1, 'scope': 'input_inventory_not_accuracy_validation',
              'archive_sha256': hashlib.sha256(payload).hexdigest(), 'archive_bytes': len(payload),
              'member_count': 0, 'suffix_counts': {}, 'duplicate_names': 0,
              'abif_magic_count': 0, 'scf_magic_count': 0,
              'genbank_records': 0, 'molecule_types': {}, 'topologies': {},
              'parse_warning_count': 0, 'parse_failures': 0, 'opaque_members': 0,
              'empty_genbank_members': 0, 'nested_archives': 0,
              'pairing_status': 'not_assessed_requires_explicit_reference_read_mapping',
              'readiness': 'not_established'}
    with zipfile.ZipFile(io.BytesIO(payload)) as archive:
        all_entries = archive.infolist()
        if len(all_entries) > MAX_MEMBERS:
            raise ValueError('member_count_limit')
        entries = [entry for entry in all_entries if not entry.is_dir()]
        if any(entry.flag_bits & 1 for entry in entries):
            raise ValueError('encrypted_member_not_supported')
        if any(entry.file_size > MAX_MEMBER for entry in entries):
            raise ValueError('member_size_limit')
        if sum(entry.file_size for entry in entries) > MAX_TOTAL:
            raise ValueError('decoded_size_limit')
        names = Counter(entry.filename for entry in entries)
        report['duplicate_names'] = sum(n - 1 for n in names.values())
        report['member_count'] = len(entries)
        suffixes, molecules, topologies = Counter(), Counter(), Counter()
        consumed = 0
        for entry in entries:
            suffix = Path(entry.filename).suffix.lower()
            suffixes[suffix or '(none)'] += 1
            with archive.open(entry) as handle:
                data = handle.read(MAX_MEMBER + 1)
            consumed += len(data)
            if len(data) > MAX_MEMBER or consumed > MAX_TOTAL:
                raise ValueError('decoded_size_limit')
            if data.startswith(b'ABIF'):
                report['abif_magic_count'] += 1  # Signature, not a claim of parseable ABI.
            elif data.startswith(b'.scf'):
                report['scf_magic_count'] += 1
            elif data.startswith((b'PK\x03\x04', b'PK\x05\x06')):
                report['nested_archives'] += 1
                report['opaque_members'] += 1
            elif suffix in ('.gb', '.gbk', '.genbank', '.gp', '.genpept'):
                count = 0
                with warnings.catch_warnings(record=True) as observed:
                    warnings.simplefilter('always')
                    try:
                        for record in SeqIO.parse(io.StringIO(data.decode('utf-8-sig')), 'genbank'):
                            count += 1
                            report['genbank_records'] += 1
                            molecule = str(record.annotations.get('molecule_type', '')).upper()
                            # Fixed categories prevent untrusted annotation text leaking into output.
                            molecules[molecule if molecule in ('DNA', 'DS-DNA', 'SS-DNA', 'MS-DNA', 'RNA', 'SS-RNA', 'DS-RNA', 'PROTEIN') else 'other_or_missing'] += 1
                            topology = str(record.annotations.get('topology', '')).lower()
                            topologies[topology if topology in ('circular', 'linear') else 'other_or_missing'] += 1
                    except (ValueError, UnicodeError):
                        report['parse_failures'] += 1
                report['parse_warning_count'] += len(observed)
                if count == 0:
                    report['empty_genbank_members'] += 1
            else:
                report['opaque_members'] += 1
        report['suffix_counts'] = dict(sorted(suffixes.items()))
        report['molecule_types'] = dict(sorted(molecules.items()))
        report['topologies'] = dict(sorted(topologies.items()))
    return report


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('archive', type=Path)
    parser.add_argument('--output', type=Path, required=True, help='New private JSON report; refuse overwrite')
    args = parser.parse_args()
    result = audit(args.archive)
    with args.output.open('x', encoding='utf-8') as stream:
        json.dump(result, stream, ensure_ascii=False, indent=2)
        stream.write('\n')
    print(json.dumps({key: result[key] for key in ['member_count', 'genbank_records', 'molecule_types', 'abif_magic_count', 'scf_magic_count', 'parse_failures', 'opaque_members', 'readiness']}))
