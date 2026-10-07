"""Compare eligible exact events only with explicitly observed reference spans."""


def compare_indels(events, records, validated, source_runs, reference, topology):
    length=len(reference);circular=topology=='circular'
    source_by_id={s['read_id']:s for s in source_runs}
    event_by_read={r['id']:[] for r in records}
    groups={}
    for event in events:
        event_by_read[event['read_id']].append(event)
        if event['eligible_for_exact_anchor_review']:
            groups.setdefault((event['kind'],event['boundary'],event['sequence']),event)
    comparisons=[]
    for key,event in groups.items():
        start=event['boundary'];span=0 if event['kind']=='insertion' else event['length_bp']
        positions=list(range(start-1,start+span+1))
        if circular:positions=[p%length for p in positions]
        sources=[]
        for record in records:
            source=source_by_id[record['id']]
            result={'read_id':record['id'],'filename':record['original_filename'],
                    'alignment_id':record.get('alignment_id'),'run_number':record.get('run_number'),
                    'support':'unassessed','reason':None,'flanks':None}
            rows=validated.get(record['id'])
            candidates=event_by_read[record['id']]
            variants=record.get('report',{}).get('variants',[])
            # Even a missing/duplicate event list must not be taken as reference evidence.
            totals_ok=isinstance(variants,list) and all(isinstance(v,dict) for v in variants)
            if totals_ok:
                for kind,field in [('insertion','inserted_bases'),('deletion','deleted_bases')]:
                    report=record.get('report',{})
                    if field in report and sum(v.get('kind')==kind for v in variants)!=report[field]:totals_ok=False
            invalid=any(any(r in e['review_reasons'] for r in ('invalid_saved_event','inconsistent_event_totals_or_duplicates')) for e in candidates)
            matching=next((e for e in candidates if e['eligible_for_exact_anchor_review'] and (e['kind'],e['boundary'],e['sequence'])==key),None)
            if not source['included'] or rows is None:
                result['reason']=source['reason'] or 'mapping_not_available'
            elif not totals_ok or invalid:
                result['reason']='invalid_saved_indel_evidence'
            elif matching:
                result.update(support='event',flanks=matching['flanks'],event_number=matching['event_number'])
            else:
                mapped={c[0]:c for c in rows}
                span_rows=[mapped.get(p) for p in positions]
                if len(set(positions))!=len(positions) or any(c is None for c in span_rows):
                    result['reason']='reference_span_not_fully_observed'
                elif any(ref not in 'ACGT' or call!=ref or q is None or q<20 for p,o,ref,call,q in span_rows):
                    result['reason']='reference_span_not_q20_matching'
                else:
                    step=-1 if record['report']['evidence']['direction']=='reverse' else 1
                    if any(b[1]-a[1]!=step for a,b in zip(span_rows,span_rows[1:])):
                        result['reason']='reference_span_not_contiguous_in_read'
                    else:
                        # Reject contradictory saved gap entries touching this exact span.
                        deletion_positions=set(positions)
                        insertion_boundaries=set(positions[1:])
                        overlapping=any((v.get('kind')=='deletion' and v.get('position') in deletion_positions)
                                        or (v.get('kind')=='insertion' and v.get('position') in insertion_boundaries) for v in variants)
                        if overlapping:
                            result['reason']='overlapping_saved_indel'
                        else:
                            flanks={}
                            for side,c in [('left',span_rows[0]),('right',span_rows[-1])]:
                                p,o,ref,call,q=c
                                flanks[side]={'reference_position':p,'original_read_position':o,'reference':ref,'call':call,'phred':q}
                            result.update(support='reference',flanks=flanks,reference_span_bases=len(span_rows),
                                          minimum_span_phred=min(c[4] for c in span_rows))
            sources.append(result)
        supporting=sum(s['support']=='event' for s in sources)
        reference_support=sum(s['support']=='reference' for s in sources)
        comparisons.append({'kind':event['kind'],'boundary':start,'sequence':event['sequence'],'length_bp':event['length_bp'],
                            'reference_segments':event['reference_segments'],
                            'event_read_count':supporting,'reference_read_count':reference_support,
                            'unassessed_read_count':sum(s['support']=='unassessed' for s in sources),
                            'conflicting_support':bool(supporting and reference_support),'sources':sources,
                            'scope':'exact_event_vs_reference','alternative_indel_alleles':'not_compared',
                            'normalization':'not_performed','whole_reference_verified':False})
    return comparisons
