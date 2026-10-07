"""Build public synthetic V1 rehearsal inputs; never collect human measurements.

Run from a source checkout/sdist. The ABIF encoder is a test utility; its peaks
are artificial and cannot validate mixed instrument traces or base calling.
"""
import argparse
import hashlib
import io
import json
from pathlib import Path
import random
import sys
import zipfile

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'tests'))
from abif_fixture import synthetic_ab1

PACK_ID = 'v1-sanger-rehearsal-01'
SEED = 202610071


def digest(data):
    return hashlib.sha256(data).hexdigest()


def write_json(path, value):
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2) + '\n', encoding='utf-8')


def build(destination):
    destination = Path(destination)
    # Do not overwrite a frozen pack, observations, or even an existing empty directory.
    destination.mkdir(parents=True, exist_ok=False)
    inputs = destination / 'participant' / 'inputs'
    inputs.mkdir(parents=True)
    evaluator = destination / 'evaluator'
    evaluator.mkdir()
    rng = random.Random(SEED)
    dna = lambda n: ''.join(rng.choice('ACGT') for _ in range(n))
    reference = dna(1200)
    repeated = dna(1200)
    repeated = repeated[:700] + repeated[100:280] + repeated[880:]
    wrong_reference = dna(1200)
    cases = []
    truth = []
    genbanks = []

    def add(case_id, target, reads, expected, operations):
        record = SeqRecord(Seq(target), id=case_id, name=case_id,
                           description='Public synthetic validation rehearsal; not instrument data',
                           annotations={'molecule_type': 'DNA', 'topology': 'circular', 'date': '01-JAN-1980'})
        record.features = [SeqFeature(FeatureLocation(60, 330, strand=1), type='misc_feature', qualifiers={'label': ['Synthetic forward region']}),
                           SeqFeature(FeatureLocation(780, 1080, strand=-1), type='misc_feature', qualifiers={'label': ['Synthetic reverse region']})]
        stream = io.StringIO()
        SeqIO.write(record, stream, 'genbank')
        gb = inputs / f'{case_id}.gb'
        gb.write_text(stream.getvalue(), encoding='utf-8')
        genbanks.append(gb)
        files = []
        for index, (calls, qualities) in enumerate(reads, 1):
            filename = f'{case_id}-r{index}.ab1'
            (inputs / filename).write_bytes(synthetic_ab1(calls, qualities, trace=True))
            files.append('inputs/' + filename)
        cases.append({'case_id': case_id, 'reference': 'inputs/' + gb.name, 'reads': files,
                      'task': 'Review the assigned reference and every read; locate differences or insufficient evidence, trace them to original data, and save a review. State what remains unverified.'})
        truth.append({'case_id': case_id, 'reference_length': len(target), 'topology': 'circular',
                      'construction': operations, 'expected_review': expected,
                      'whole_reference_verified': False})

    def altered(calls, position):
        alternate = next(b for b in 'ACGT' if b != calls[position])
        return calls[:position] + alternate + calls[position + 1:]

    add('A01', reference, [(reference[80:580], [38] * 500)],
        {'paired_intervals': [[80, 580]], 'differences': [], 'uncovered_bases': 700},
        {'source_slice': [80, 580], 'direction': 'forward'})
    calls = altered(reference[200:760], 200)
    add('A02', reference, [(str(Seq(calls).reverse_complement()), [38] * 560)],
        {'substitution_position': 400, 'original_read_position': 359, 'phred': 38, 'direction': 'reverse'},
        {'source_slice': [200, 760], 'substitute_offset': 200, 'reverse_complement': True})
    calls = altered(reference[250:800], 200)
    quality = [5] * 20 + [38] * 510 + [8] * 20
    quality[200] = 9
    add('A03', reference, [(calls, quality)],
        {'substitution_position': 450, 'original_read_position': 200, 'phred': 9,
         'conclusion': 'Low-quality conflicting call requires review; not a confirmed mutation.'},
        {'source_slice': [250, 800], 'substitute_offset': 200, 'terminal_low_quality_bases': [20, 20]})
    calls = reference[600:1000]
    add('A04', reference, [(calls, [38] * 400), (altered(calls, 220), [38] * 400)],
        {'conflicting_position': 820, 'conclusion': 'Two reads disagree; no automatic consensus or biological mixture diagnosis.'},
        {'source_slice': [600, 1000], 'second_read_substitute_offset': 220})
    add('A05', repeated, [(repeated[100:280], [38] * 180)],
        {'equally_supported_intervals': [[100, 280], [700, 880]], 'conclusion': 'Mapping is not unique.'},
        {'copy_from': [100, 280], 'copy_to': [700, 880]})
    add('A06', reference, [(reference[1020:] + reference[:220], [38] * 400)],
        {'paired_intervals': [[0, 220], [1020, 1200]], 'crosses_origin': True, 'uncovered_bases': 800},
        {'source_slices_in_read_order': [[1020, 1200], [0, 220]]})
    add('A07', reference, [(reference[500:580], [38] * 80)],
        {'paired_intervals': [[500, 580]], 'uncovered_bases': 1120, 'conclusion': 'Short exact evidence does not verify the construct.'},
        {'source_slice': [500, 580]})
    add('A08', wrong_reference, [(reference[300:800], [38] * 500)],
        {'conclusion': 'Read comes from another reference. Reject the verification claim; chance local matches are possible.'},
        {'source_reference_case': 'A01', 'source_slice': [300, 800], 'assigned_reference_is_unrelated': True})
    add('A09', reference, [], {'conclusion': 'Invalid ABIF must be rejected, with no saved successful analysis.'}, {'invalid_file': True})
    (inputs / 'A09-r1.ab1').write_bytes(b'Public synthetic invalid ABIF fixture\n')
    cases[-1]['reads'] = ['inputs/A09-r1.ab1']
    inserted = next(b for b in 'ACGT' if b != reference[650]) + 'G' + next(b for b in 'ACGT' if b != reference[649])
    calls = reference[400:650] + inserted + reference[650:900]
    add('A10', reference, [(str(Seq(calls).reverse_complement()), [38] * 503)],
        {'insertion_boundary': 650, 'inserted_sequence_reference_orientation': inserted,
         'original_read_interval': [250, 253], 'direction': 'reverse',
         'conclusion': 'Insertion is between reference bases; do not invent a covered reference base for it.'},
        {'source_slice': [400, 900], 'insert_boundary': 650, 'inserted': inserted, 'reverse_complement': True})
    with zipfile.ZipFile(destination / 'participant' / 'reference-library.zip', 'w') as archive:
        for path in genbanks:
            archive.writestr(zipfile.ZipInfo(path.name, (1980, 1, 1, 0, 0, 0)), path.read_bytes())
    write_json(destination / 'participant' / 'tasks.json', {'pack_id': PACK_ID, 'synthetic_only': True, 'cases': cases})
    write_json(evaluator / 'construction-ledger.json', {'pack_id': PACK_ID, 'seed': SEED,
               'coordinates': 'zero-based half-open; insertion is a boundary',
               'oracle_source': 'generation operations, not either application alignment output',
               'independently_reviewed': False, 'cases': truth})
    write_json(destination / 'session-template.json', {
        'pack_id': PACK_ID, 'record_type': None, 'participant_id': None,
        'd0': None, 'tool_order': None, 'prior_experience': None,
        'software_version': None, 'parameters': None, 'manifest_sha256': None,
        'setup_seconds': None, 'observations': [], 'retention': [],
        'instructions': 'Use record_type technical_rehearsal or human_session. Null means unmeasured, never zero. Do not put agent timings in human records.'})
    write_json(destination / 'observation-template.json', {
        'case_id': None, 'tool': None, 'started_at': None, 'active_seconds': None,
        'waiting_seconds': None, 'total_seconds': None, 'interruptions': None,
        'completed': None, 'correct': None, 'help_count': None, 'takeover': None,
        'false_positive': None, 'false_negative': None, 'severe_error': None,
        'unresolved': None, 'format_supported': None, 'conversion_seconds': None,
        'evidence_path': None, 'friction': None})
    write_json(destination / 'retention-template.json', {
        'participant_id': None, 'd0': None, 'checkpoint': None, 'observed_at': None,
        'task_opportunities': None, 'self_initiated_completions': None,
        'reminder_or_help': None, 'alternative_tool': None, 'reason': None,
        'feedback_missing': None})
    files = [{'path': str(p.relative_to(destination)), 'bytes': p.stat().st_size,
              'sha256': digest(p.read_bytes())} for p in sorted(destination.rglob('*')) if p.is_file()]
    manifest = {'pack_id': PACK_ID, 'synthetic_only': True, 'status': 'technical_rehearsal_only',
                'holdout': False, 'human_sessions': 0, 'files': files}
    write_json(destination / 'manifest.json', manifest)
    (destination / 'manifest.sha256').write_text(digest((destination / 'manifest.json').read_bytes()) + '\n')
    return manifest


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True, help='New, nonexistent directory; keep observations private')
    args = parser.parse_args()
    build(args.output)
    print(json.dumps({'pack_id': PACK_ID, 'output': str(args.output), 'human_sessions': 0}))
