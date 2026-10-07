"""Conservative union of saved alignment coordinates; no consensus or QC verdict."""
from collections import Counter, defaultdict
from hashlib import sha256
import json

METHOD = 'latest_unambiguous_alignment_union_v1'


def summarize_sanger_reads(reference: dict, records: list[dict]) -> dict:
    length = reference['length_bp']
    if type(length) is not int or length < 1:
        raise ValueError('Reference must have a positive integer length')
    duplicate_hashes = Counter(r['file_sha256'] for r in records)
    events = defaultdict(int)
    sources = []
    for record in records:
        report = record.get('report') or {}
        evidence = report.get('evidence') or {}
        params = evidence.get('parameters') or {}
        manifest = record.get('manifest') or {}
        reasons = []
        flags = evidence.get('review_flags') or []
        intervals = evidence.get('reference_covered_intervals')
        if not record.get('alignment_id'):
            reasons.append('not_analyzed')
        else:
            if evidence.get('scope') != 'local_read_alignment' or params.get('evidence_version') not in {2,3}:
                reasons.append('unsupported_or_missing_evidence')
            if evidence.get('direction') not in {'forward','reverse'}:
                reasons.append('unresolved_direction')
            reasons.extend(f for f in ('ambiguous_alignment','candidate_search_truncated') if f in flags)
            if manifest.get('reference_sha256') != reference['sha256']:
                reasons.append('reference_hash_missing_or_mismatch')
            if manifest.get('read_sha256') != record['file_sha256']:
                reasons.append('read_hash_missing_or_mismatch')
            if params.get('reference_topology',reference['topology']) != reference['topology']:
                reasons.append('reference_topology_mismatch')
            valid = isinstance(intervals,list) and bool(intervals)
            previous_end = -1
            if valid:
                for interval in intervals:
                    if not (isinstance(interval,list) and len(interval)==2 and all(type(x) is int for x in interval)
                            and 0<=interval[0]<interval[1]<=length and interval[0]>=previous_end):
                        valid=False
                        break
                    previous_end=interval[1]
            if not valid or sum(b-a for a,b in intervals) != report.get('aligned_bases'):
                reasons.append('invalid_coverage_intervals')
        if duplicate_hashes[record['file_sha256']] > 1:
            reasons.append('duplicate_source_file')
        if not reasons:
            for start,end in intervals:
                events[start]+=1
                events[end]-=1
        sources.append({'read_id':record['id'],'filename':record['original_filename'],
                        'read_sha256':record['file_sha256'],'alignment_id':record.get('alignment_id'),
                        'run_number':record.get('run_number'),'included':not reasons,
                        'exclusion_reasons':sorted(set(reasons)),'review_flags':flags,
                        'direction':evidence.get('direction'),'source_reference_sha256':manifest.get('reference_sha256')})
    depth, previous = 0, 0
    segments = []
    for boundary,delta in sorted(events.items()):
        if boundary>previous and depth:
            if segments and segments[-1]['end']==previous and segments[-1]['depth']==depth:
                segments[-1]['end']=boundary
            else:
                segments.append({'start':previous,'end':boundary,'depth':depth})
        depth+=delta
        previous=boundary
    union=[]
    for segment in segments:
        if union and union[-1][1]==segment['start']:
            union[-1][1]=segment['end']
        else:
            union.append([segment['start'],segment['end']])
    gaps=[];cursor=0
    for start,end in union:
        if start>cursor:gaps.append([cursor,start])
        cursor=end
    if cursor<length:gaps.append([cursor,length])
    gap_regions=[{'segments':[interval],'length_bp':interval[1]-interval[0],'wraps_origin':False} for interval in gaps]
    if reference['topology']=='circular' and len(gaps)>1 and gaps[0][0]==0 and gaps[-1][1]==length:
        gap_regions=[{'segments':[gaps[-1],gaps[0]],'length_bp':length-gaps[-1][0]+gaps[0][1],'wraps_origin':True},*gap_regions[1:-1]]
    covered=sum(b-a for a,b in union)
    payload={'schema_version':1,'method':METHOD,'reference':reference,
             'coordinate_system':'zero-based-half-open','covered_intervals':union,'uncovered_intervals':gaps,
             'gap_regions':gap_regions,'depth_segments':segments,'covered_bases':covered,
             'uncovered_bases':length-covered,'covered_fraction':round(covered/length,6),
             'overlap_bases':sum(s['end']-s['start'] for s in segments if s['depth']>=2),
             'maximum_read_depth':max((s['depth'] for s in segments),default=0),
             'included_read_count':sum(r['included'] for r in sources),'excluded_read_count':sum(not r['included'] for r in sources),
             'source_runs':sources,'whole_reference_verified':False,
             'quality_screening':'not_applied_to_union','conflict_analysis':'not_evaluated',
             'source_files_rechecked':False}
    payload['snapshot_sha256']=sha256(json.dumps(payload,sort_keys=True,separators=(',',':')).encode()).hexdigest()
    return payload
