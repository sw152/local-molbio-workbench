"""Validated Q20 paired-base coverage and disagreements, never a consensus."""
from collections import defaultdict
from hashlib import sha256
from .sanger_mapping import build_base_mapping, intervals
from .sanger_indel_review import review_indels
from .sanger_indel_comparison import compare_indels

COMPLEMENT=str.maketrans('ACGTN','TGCAN')


def validated_columns(record, reference):
    """Reject a whole mapping if any coordinate, call or saved statistic disagrees."""
    evidence=record['report'].get('evidence',{})
    mapping=evidence.get('base_mapping')
    if not mapping:
        return None,'mapping_not_recorded'
    if evidence.get('parameters',{}).get('evidence_version')!=4 or not isinstance(mapping,dict) or mapping.get('schema_version')!=1:
        return None,'unsupported_mapping'
    read=record.get('base_sequence')
    if not isinstance(read,str) or not read or not reference:
        return None,'source_sequence_unavailable'
    qualities=record.get('qualities')
    if qualities is not None and (not isinstance(qualities,list) or (qualities and len(qualities)!=len(read))):
        return None,'stored_quality_invalid'
    try:
        blocks=mapping['blocks']
        if not isinstance(blocks,list) or not blocks:raise ValueError()
        trim=evidence['trimming']
        lo,hi=trim['original_read_start'],trim['original_read_end']
        if type(lo) is not int or type(hi) is not int or not 0<=lo<hi<=len(read) or trim['original_length']!=len(read):raise ValueError()
        direction=evidence['direction'];step=-1 if direction=='reverse' else 1
        if direction not in {'forward','reverse'}:raise ValueError()
        columns=[];seen=set();previous_original=None;previous_reference=None;travel=0
        for b in blocks:
            start,end,original,stride=b['reference_start'],b['reference_end'],b['original_read_start'],b['original_read_step']
            if any(type(n) is not int for n in (start,end,original,stride)):raise ValueError()
            if not 0<=start<end<=len(reference) or stride!=step:raise ValueError()
            ref,calls,qs=b['reference_bases'],b['read_bases'],b['phred']
            if not isinstance(ref,str) or not isinstance(calls,str) or not isinstance(qs,list):raise ValueError()
            if len(ref)!=end-start or len(calls)!=end-start or len(qs)!=end-start:raise ValueError()
            for i in range(end-start):
                p,o=start+i,original+i*step
                q=qs[i]
                if not lo<=o<hi or p in seen or ref[i]!=reference[p] or calls[i] not in 'ACGTN':raise ValueError()
                if q is not None and (type(q) is not int or q<0):raise ValueError()
                expected=read[o].translate(COMPLEMENT) if step==-1 else read[o]
                if calls[i]!=expected or (qualities and q!=qualities[o]):raise ValueError()
                if previous_original is not None and (o-previous_original)*step<=0:raise ValueError()
                if previous_reference is not None:
                    delta=p-previous_reference
                    if delta<=0:
                        if evidence.get('parameters',{}).get('reference_topology')!='circular':raise ValueError()
                        delta+=len(reference)
                    travel+=delta
                    if travel>=len(reference):raise ValueError()
                seen.add(p);previous_original=o;previous_reference=p
                columns.append((p,o,ref[i],calls[i],q))
        if intervals(seen)!=evidence['reference_covered_intervals'] or len(columns)!=record['report']['aligned_bases']:raise ValueError()
        expected=build_base_mapping(columns,direction,20)
        if any(mapping.get(k)!=v for k,v in expected.items()):raise ValueError()
        return columns,None
    except (KeyError,TypeError,ValueError,IndexError):
        return None,'invalid_base_mapping'


def summarize_quality(reference, reference_sequence, records, source_runs):
    source_by_id={r['read_id']:r for r in source_runs}
    sources=[];calls=defaultdict(list);indels=[];validated={}
    reference_valid=isinstance(reference_sequence,str) and len(reference_sequence)==reference['length_bp'] and sha256(reference_sequence.encode()).hexdigest()==reference['sha256']
    for record in records:
        source=source_by_id[record['id']]
        reason=None
        if not source['included']:reason='alignment_excluded'
        elif not reference_valid:reason='reference_sequence_unavailable_or_mismatch'
        if reason:columns=None
        else:columns,reason=validated_columns(record,reference_sequence)
        sources.append({'read_id':record['id'],'alignment_id':record.get('alignment_id'),
                        'filename':record['original_filename'],'included':reason is None,'reason':reason})
        if columns is None:continue
        validated[record['id']]=columns
        indels.extend(review_indels(record,reference_sequence,reference['topology'],columns))
        for pos,original,ref,call,q in columns:
            if ref not in 'ACGT' or call not in 'ACGT' or q is None or q<20:continue
            calls[pos].append({'read_id':record['id'],'alignment_id':record['alignment_id'],'run_number':record['run_number'],
                               'filename':record['original_filename'],'original_read_position':original,'base':call,'phred':q})
    differences=[];conflicts=[];reference_only=[]
    for pos,evidence in sorted(calls.items()):
        alleles={c['base'] for c in evidence}
        if alleles=={reference_sequence[pos]}:
            reference_only.append(pos)
        else:
            item={'position':pos,'reference':reference_sequence[pos],'conflict':len(alleles)>1,'calls':evidence}
            differences.append(item)
            if item['conflict']:conflicts.append(pos)
    comparisons=compare_indels(indels,records,validated,sources,reference_sequence,reference['topology']) if reference_valid else []
    for event in indels:
        event['cross_read_comparison']='see_indel_review_comparisons' if event['eligible_for_exact_anchor_review'] else 'not_evaluated'
    valid=sum(s['included'] for s in sources)
    overlap=sum(len(v)>=2 for v in calls.values())
    return {'schema_version':1,'method':'validated_q20_paired_base_review_v1','quality_threshold':20,
            'assessment':'paired_bases_only' if valid else 'not_evaluated',
            'included_read_count':valid,'excluded_read_count':len(sources)-valid,'source_runs':sources,
            'callable_bases':len(calls),'callable_fraction':round(len(calls)/reference['length_bp'],6),
            'callable_intervals':intervals(calls),'overlap_bases':overlap,
            'reference_only_bases':len(reference_only),'reference_only_intervals':intervals(reference_only),
            'non_reference_position_count':len(differences),'differences':differences,
            'conflict_position_count':len(conflicts),'conflict_intervals':intervals(conflicts),
            'indel_conflicts':'exact_event_vs_reference_only' if comparisons else 'not_evaluated','consensus':'not_generated',
            'indel_review':{'method':'saved_event_q20_flank_audit_v1','events':indels,
                            'eligible_event_count':sum(e['eligible_for_exact_anchor_review'] for e in indels),
                            'normalization':'not_performed','cross_read_comparison':'exact_event_vs_reference_only' if comparisons else 'not_performed',
                            'comparison_method':'q20_contiguous_span_event_vs_reference_v1','comparisons':comparisons,
                            'conflicting_group_count':sum(c['conflicting_support'] for c in comparisons)},
            'source_files_rechecked':False,'whole_reference_verified':False}
