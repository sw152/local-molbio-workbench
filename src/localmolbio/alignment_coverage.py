"""Quality-unassessed paired-position evidence, counting each record once per track."""
from array import array
import re
from .alignment_evidence import MAX_REFERENCE, MAX_READS, _intervals


class CoverageError(ValueError):
    pass


def _integer(value, low, high):
    if type(value) is not int or not low <= value <= high:
        raise CoverageError('invalid_coverage_geometry')
    return value


def _union(intervals):
    merged = []
    for start, end in sorted(intervals):
        if merged and start <= merged[-1][1]:
            merged[-1][1] = max(end, merged[-1][1])
        else:
            merged.append([start, end])
    return merged


def _geometry(hit, read_length, length, topology):
    """Check stored blocks against CIGAR; never turn deletions into paired support."""
    span = _integer(hit['reference_span_bp'], 1, length)
    start = _integer(hit['reference_intervals'][0][0], 0, length-1)
    if topology == 'linear' and start+span > length:
        raise CoverageError('invalid_coverage_geometry')
    if hit['reference_intervals'] != _intervals(start, span, length):
        raise CoverageError('invalid_coverage_geometry')
    qs = _integer(hit['query_start'], 0, read_length-1)
    qe = _integer(hit['query_end'], qs+1, read_length)
    strand = hit['strand']
    if strand not in ('+', '-'):
        raise CoverageError('invalid_coverage_geometry')
    operations = re.findall(r'([1-9][0-9]*)([=XID])', hit['cigar'])
    if not operations or ''.join(a+b for a,b in operations) != hit['cigar']:
        raise CoverageError('invalid_coverage_geometry')
    paired, deleted, blocks = [], [], []
    qpos = tpos = insertions = deletions = 0
    for size, op in operations:
        size = int(size)
        qnext = qpos+(size if op in '=XI' else 0)
        tnext = tpos+(size if op in '=XD' else 0)
        if qnext > qe-qs or tnext > span:
            raise CoverageError('invalid_coverage_geometry')
        if op in '=X':
            intervals = _intervals(start+tpos, size, length)
            paired.extend(intervals)
            blocks.append({'query_start':qs+qpos if strand=='+' else qe-qnext,
                           'query_end':qs+qnext if strand=='+' else qe-qpos,
                           'query_step':1 if strand=='+' else -1,
                           'reference_intervals':intervals})
        elif op == 'D':
            deleted.extend(_intervals(start+tpos, size, length));deletions += size
        else:
            insertions += size
        qpos, tpos = qnext, tnext
    if (tpos != span or qpos != qe-qs or blocks != hit['paired_blocks']
            or deletions != hit['deleted_bases'] or insertions != hit['inserted_bases']):
        raise CoverageError('invalid_coverage_geometry')
    return paired, deleted


def summarize(evidence):
    """Derive exact runs; inputs must be queue-published and provenance-checked by caller."""
    try:
        length = _integer(evidence['reference_length_bp'], 1, MAX_REFERENCE)
        topology = evidence['topology']
        if topology not in ('linear','circular'):
            raise CoverageError('invalid_coverage_geometry')
        reads = evidence['reads']
        _integer(len(reads), 1, MAX_READS)
        tracks = [array('i',[0])*(length+1) for _ in range(4)]
        counts = dict(total=len(reads), single_reported_alignment=0, ambiguous_or_withheld=0,
                      no_alignment_reported=0, withheld=0)
        for read in reads:
            hits, withheld = read['alignments'], read['withheld']
            if withheld:
                counts['withheld'] += 1
            if not hits and not withheld:
                counts['no_alignment_reported'] += 1
            single = len(hits)==1 and not withheld
            if single:
                counts['single_reported_alignment'] += 1
            elif hits or withheld:
                counts['ambiguous_or_withheld'] += 1
            paired, deleted = [], []
            for hit in hits:
                p,d = _geometry(hit,read['length_bp'],length,topology)
                paired.extend(p);deleted.extend(d)
            # Multiple candidate hits and overlapping blocks of ONE record never add depth.
            for intervals, track in [(paired,0 if single else 1),(deleted,2 if single else 3)]:
                for start,end in _union(intervals):
                    tracks[track][start] += 1;tracks[track][end] -= 1
        depths = [0]*4
        segments = []
        totals = dict(single_reported_alignment_bases=0, ambiguous_only_bases=0,
                      unpaired_bases=0, deletion_evidence_bases=0)
        depth_sum = max_depth = 0
        for position in range(length):
            for i in range(4):depths[i] += tracks[i][position]
            if depths[0]:totals['single_reported_alignment_bases'] += 1
            elif depths[1]:totals['ambiguous_only_bases'] += 1
            else:totals['unpaired_bases'] += 1
            if depths[2]+depths[3]:totals['deletion_evidence_bases'] += 1
            depth = depths[0]+depths[1]
            depth_sum += depth;max_depth = max(max_depth,depth)
            if segments and segments[-1]['depths']==depths:
                segments[-1]['end'] = position+1
            else:
                segments.append(dict(start=position,end=position+1,depths=depths.copy()))
        return dict(schema='localmolbio.alignment-paired-coverage',schema_version=1,
                    scope='reported_paired_positions_quality_unassessed',reference_length_bp=length,
                    topology=topology,coordinate_convention='zero_based_half_open',read_counts=counts,
                    **totals,paired_bases=length-totals['unpaired_bases'],
                    max_record_depth=max_depth,mean_record_depth=depth_sum/length,
                    depth_fields=['single_reported_alignment','ambiguous_alignment','single_deletion','ambiguous_deletion'],
                    base_quality_used=False,uniqueness_assessed=False,independent_molecules_assessed=False,
                    whole_reference_verified=False,segments=segments)
    except (KeyError, IndexError, TypeError, ValueError) as exc:
        if isinstance(exc,CoverageError):raise
        raise CoverageError('invalid_coverage_geometry') from exc


def regions(summary, kind):
    """Linearized exact intervals; circular origin-end intervals remain separate."""
    predicates = {'unpaired':lambda d:not(d[0]+d[1]),
                  'ambiguous_only':lambda d:not d[0] and d[1]>0,
                  'deletion':lambda d:d[2]+d[3]>0}
    if kind not in predicates:
        raise CoverageError('invalid_coverage_region_kind')
    return _union([(s['start'],s['end']) for s in summary['segments'] if predicates[kind](s['depths'])])
