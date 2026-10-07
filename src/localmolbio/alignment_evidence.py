"""Bounded DNA alignment evidence core; no consensus, variants or pass verdict."""
from hashlib import sha256
import json
from pathlib import Path
import re
import subprocess
import time

PRESETS = {'ont-noisy': 'map-ont', 'ont-high-accuracy': 'lr:hq', 'pacbio-hifi': 'map-hifi', 'short-single': 'sr'}
MAX_REFERENCE = 200_000
MAX_READS = 5_000
MAX_READ = 100_000
MAX_BASES = 5_000_000
MAX_OUTPUT = 64 * 1024**2
MAX_LOG = 8 * 1024**2


class AlignmentError(ValueError):
    pass


def _validate(reference, reads, topology):
    if topology not in ('linear', 'circular'):
        raise AlignmentError('explicit_linear_or_circular_topology_required')
    if not isinstance(reference, str) or not 1 <= len(reference) <= MAX_REFERENCE or not re.fullmatch('[ACGTN]+', reference):
        raise AlignmentError('reference_must_be_uppercase_acgtn_within_limit')
    if not isinstance(reads, list) or not 1 <= len(reads) <= MAX_READS:
        raise AlignmentError('read_count_outside_core_limit')
    if any(not isinstance(r, str) or not 1 <= len(r) <= MAX_READ or not re.fullmatch('[ACGTN]+', r) for r in reads):
        raise AlignmentError('reads_must_be_uppercase_acgtn_within_limit')
    if sum(map(len, reads)) > MAX_BASES:
        raise AlignmentError('total_read_bases_outside_core_limit')


def _intervals(start, length, reference_length):
    start %= reference_length
    first = min(length, reference_length - start)
    return [[start, start + first]] + ([[0, length - first]] if first < length else [])


def parse_paf(lines, reference, reads, topology):
    """Validate base-level eqx CIGAR against supplied sequence text before reporting."""
    _validate(reference, reads, topology)
    target = reference * (2 if topology == 'circular' else 1)
    length = len(reference)
    output = [{'read_index': i, 'length_bp': len(read), 'alignments': [], 'withheld': [],
               'status': 'no_alignment_reported', 'uniqueness_assessed': False} for i, read in enumerate(reads)]
    seen = [{} for _ in reads]
    for number, line in enumerate(lines, 1):
        if number > MAX_READS * 100 or len(line) > 4 * MAX_READ:
            raise AlignmentError('paf_output_limit_exceeded')
        if not line.strip():
            continue
        fields = line.rstrip('\n').split('\t')
        try:
            if len(fields) < 12 or not re.fullmatch(r'q\d+', fields[0]):
                raise AlignmentError('invalid_paf_record')
            index = int(fields[0][1:])
            if not 0 <= index < len(reads) or fields[0] != f'q{index}':
                raise AlignmentError('unexpected_query_identity')
            qlen, qs, qe = map(int, fields[1:4])
            strand = fields[4]
            tlen, ts, te, matches, block, mapq = map(int, fields[6:12])
            if (fields[5] != 'reference' or tlen != len(target) or qlen != len(reads[index])
                    or strand not in ('+', '-') or not 0 <= qs < qe <= qlen
                    or not 0 <= ts < te <= tlen or not 0 <= matches <= block
                    or block <= 0 or not 0 <= mapq <= 255):
                raise AlignmentError('invalid_paf_coordinates_or_counts')
            tags = {}
            for tag in fields[12:]:
                key, kind, value = tag.split(':', 2)
                if key in tags:
                    raise AlignmentError('duplicate_paf_tag')
                tags[key] = (kind, value)
            if tags.get('tp') not in (('A', 'P'), ('A', 'S')) or tags.get('cg', ('',))[0] != 'Z':
                raise AlignmentError('unsupported_alignment_type_or_missing_cigar')
            cigar = tags['cg'][1]
            operations = re.findall(r'([1-9][0-9]*)([=XID])', cigar)
            if not operations or ''.join(a+b for a,b in operations) != cigar:
                raise AlignmentError('unsupported_or_invalid_cigar')
            oriented = reads[index][qs:qe]
            if strand == '-':
                oriented = oriented.translate(str.maketrans('ACGTN', 'TGCAN'))[::-1]
            qpos, tpos, exact, ambiguous, mismatch, insertions, deletions = 0, ts, 0, 0, 0, 0, 0
            ambiguous_gaps = 0
            paired_blocks = []
            for size, op in operations:
                size = int(size)
                qnext = qpos + (size if op in '=XI' else 0)
                tnext = tpos + (size if op in '=XD' else 0)
                if qnext > len(oriented) or tnext > te:
                    raise AlignmentError('cigar_exceeds_alignment_bounds')
                if op in '=X':
                    for a,b in zip(oriented[qpos:qnext], target[tpos:tnext]):
                        if (op == '=') != (a == b):
                            raise AlignmentError('cigar_base_disagreement')
                        if 'N' in (a,b):
                            ambiguous += 1
                        elif a == b:
                            exact += 1
                        else:
                            mismatch += 1
                    paired_blocks.append({'query_start': qs+qpos if strand == '+' else qe-qnext,
                                          'query_end': qs+qnext if strand == '+' else qe-qpos,
                                          'query_step': 1 if strand == '+' else -1,
                                          'target_start': tpos, 'target_end': tnext})
                elif op == 'I':
                    insertions += size
                    ambiguous_gaps += oriented[qpos:qnext].count('N')
                else:
                    deletions += size
                    ambiguous_gaps += target[tpos:tnext].count('N')
                qpos, tpos = qnext, tnext
            columns = exact + ambiguous + mismatch + insertions + deletions
            # Pinned minimap2 2.31 excludes all ambiguous bases (including gaps)
            # from PAF block length; nn records them. Keep them in our denominator.
            nn = ambiguous + ambiguous_gaps
            if qpos != qe-qs or tpos != te or columns - nn != block:
                raise AlignmentError('cigar_span_or_block_disagreement')
            if tags.get('nn', ('i', '0')) != ('i', str(nn)):
                raise AlignmentError('paf_ambiguous_count_disagreement')
            if exact != matches:
                raise AlignmentError('paf_match_count_disagreement')
            if te-ts > length:
                output[index]['withheld'].append({'reason': 'more_than_one_reference_traversal', 'query_start': qs, 'query_end': qe})
                continue
            intervals = _intervals(ts, te-ts, length)
            # Retain original query coordinates and normalized reference intervals.
            paired_blocks = [{'query_start': b['query_start'], 'query_end': b['query_end'], 'query_step': b['query_step'],
                              'reference_intervals': _intervals(b['target_start'], b['target_end']-b['target_start'], length)} for b in paired_blocks]
            key = (qs, qe, strand, ts % length, te-ts, cigar)
            raw = {'target_start': ts, 'target_end': te, 'type': tags['tp'][1], 'mapq': None if mapq == 255 else mapq}
            if key in seen[index]:
                seen[index][key]['raw_hits'].append(raw)
                continue
            hit = {'query_start': qs, 'query_end': qe, 'strand': strand,
                   'reference_intervals': intervals, 'reference_span_bp': te-ts,
                   'crosses_origin': len(intervals) == 2, 'cigar': cigar,
                   'paired_blocks': paired_blocks, 'exact_matches': exact, 'mismatches': mismatch,
                   'ambiguous_pairs': ambiguous, 'inserted_bases': insertions, 'deleted_bases': deletions,
                   'alignment_columns': columns, 'reported_paf_block_length': block,
                   'ambiguous_gap_bases': ambiguous_gaps, 'local_exact_identity': exact / columns,
                   'aligned_query_fraction': (qe-qs)/qlen, 'raw_hits': [raw]}
            seen[index][key] = hit
            output[index]['alignments'].append(hit)
        except (IndexError, TypeError, ValueError) as exc:
            if isinstance(exc, AlignmentError):
                raise
            raise AlignmentError('invalid_paf_record') from exc
    for item in output:
        count = len(item['alignments'])
        item['status'] = ('multiple_alignments_reported' if count > 1 else 'alignment_reported') if count else 'no_alignment_reported'
        if item['withheld']:
            item['requires_review'] = True
        item['requires_review'] = item.get('requires_review', False) or count != 1
    return output


def _file_hash(path):
    digest = sha256()
    with path.open('rb') as stream:
        while chunk := stream.read(1024**2):
            digest.update(chunk)
    return digest.hexdigest()


def align(reference, reads, *, topology, data_type, binary, output_dir, pulse=lambda: None, timeout_seconds=60):
    """Local bounded core. output_dir must be new; no DB job or public API is created."""
    _validate(reference, reads, topology)
    if data_type not in PRESETS:
        raise AlignmentError('explicit_supported_data_type_required')
    if type(timeout_seconds) is not int or not 1 <= timeout_seconds <= 300:
        raise AlignmentError('timeout_outside_limit')
    binary = Path(binary).resolve(strict=True)
    before = _file_hash(binary)
    try:
        version_run = subprocess.run([str(binary), '--version'], capture_output=True, timeout=5, check=True)
        version = version_run.stdout.decode().strip()
    except (OSError, subprocess.SubprocessError, UnicodeError) as exc:
        raise AlignmentError('aligner_version_unavailable') from exc
    if version != '2.31-r1302':
        raise AlignmentError('unvalidated_minimap2_version')
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=False)
    target = reference * (2 if topology == 'circular' else 1)
    (output_dir/'reference.fa').write_text('>reference\n'+target+'\n')
    (output_dir/'reads.fa').write_text(''.join(f'>q{i}\n{read}\n' for i,read in enumerate(reads)))
    input_hashes = {name: _file_hash(output_dir/name) for name in ('reference.fa', 'reads.fa')}
    args = ['-x', PRESETS[data_type], '--frag=no', '-c', '--eqx', '--secondary=yes', '-N', '20', '-p', '0.5',
            '--seed', '11', '-t', '2', '-K', '5m', 'reference.fa', 'reads.fa']
    paf, log = output_dir/'alignments.paf', output_dir/'aligner.log'
    pulse()
    with paf.open('wb') as stdout, log.open('wb') as stderr:
        process = subprocess.Popen([str(binary), *args], cwd=output_dir, stdout=stdout, stderr=stderr)
        started = time.monotonic()
        try:
            while process.poll() is None:
                pulse()
                if time.monotonic()-started > timeout_seconds:
                    raise AlignmentError('alignment_timeout')
                if paf.stat().st_size > MAX_OUTPUT or log.stat().st_size > MAX_LOG:
                    raise AlignmentError('alignment_output_limit_exceeded')
                time.sleep(0.05)
            if process.returncode != 0:
                raise AlignmentError('aligner_process_failed')
        finally:
            if process.poll() is None:
                process.kill()
            process.wait()
    pulse()
    if paf.stat().st_size > MAX_OUTPUT or log.stat().st_size > MAX_LOG:
        raise AlignmentError('alignment_output_limit_exceeded')
    if before != _file_hash(binary):
        raise AlignmentError('aligner_changed_during_run')
    if any(_file_hash(output_dir/name) != digest for name,digest in input_hashes.items()):
        raise AlignmentError('alignment_inputs_changed_during_run')
    with paf.open() as stream:
        evidence = parse_paf(stream, reference, reads, topology)
    result = {'schema': 'localmolbio.alignment-evidence-core', 'schema_version': 1,
              'scope': 'reported_local_alignments_only', 'tool': {'name': 'minimap2', 'version': version, 'binary_sha256': before},
              'arguments': args, 'data_type': data_type, 'topology': topology,
              'coordinate_convention': 'zero_based_half_open', 'reference_length_bp': len(reference),
              'limits': {'reference_bases': MAX_REFERENCE, 'reads': MAX_READS, 'read_bases': MAX_READ,
                         'total_read_bases': MAX_BASES, 'paf_bytes': MAX_OUTPUT, 'log_bytes': MAX_LOG,
                         'process_timeout_seconds': timeout_seconds},
              'reference_sha256': sha256(reference.encode()).hexdigest(),
              'query_sha256': [sha256(r.encode()).hexdigest() for r in reads],
              'artifacts_sha256': {name: _file_hash(output_dir/name) for name in ('reference.fa','reads.fa','alignments.paf','aligner.log')},
              'secondary_limit': 20, 'candidate_search_exhaustive': False,
              'base_quality_used': False, 'pairing_used': False, 'circular_mapq_recalibrated': False,
              'consensus_performed': False, 'whole_reference_verified': False, 'reads': evidence}
    if any(result['artifacts_sha256'][name] != digest for name,digest in input_hashes.items()):
        raise AlignmentError('alignment_inputs_changed_during_run')
    pulse()
    (output_dir/'result.json').write_text(json.dumps(result, indent=2))
    return result
