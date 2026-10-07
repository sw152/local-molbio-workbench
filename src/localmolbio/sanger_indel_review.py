"""Conservative event/anchor audit. No normalization, cross-read verdict or consensus."""
COMPLEMENT=str.maketrans('ACGTN','TGCAN')


def review_indels(record, reference, topology, columns):
    report=record['report'];direction=report['evidence']['direction'];step=-1 if direction=='reverse' else 1
    length=len(reference);circular=topology=='circular';mapped={c[0]:c for c in columns};paired_original={c[1] for c in columns}
    variants=report.get('variants',[]);groups=[]
    # Preserve saved traversal order and keep mixed gap kinds as separate events.
    for index,v in enumerate(variants):
        if v.get('kind') not in {'insertion','deletion'}:continue
        adjacent=False
        if groups and groups[-1][-1][0]==index-1:
            previous=groups[-1][-1][1]
            if previous['kind']==v['kind']:
                if v['kind']=='insertion':
                    adjacent=(type(v.get('original_read_position')) is int and type(previous.get('original_read_position')) is int
                              and v.get('position')==previous.get('position') and v['original_read_position']==previous['original_read_position']+step)
                elif type(v.get('position')) is int and type(previous.get('position')) is int:
                    adjacent=v['position']==((previous['position']+1)%length if circular else previous['position']+1)
        if adjacent:groups[-1].append((index,v))
        else:groups.append([(index,v)])
    events=[]
    for number,group in enumerate(groups,1):
        entries=[v for _,v in group];kind=entries[0]['kind'];start=entries[0].get('position')
        sequence=''.join(str(v.get('read' if kind=='insertion' else 'reference','')) for v in entries)
        reasons=[];original_positions=[];qualities=[];positions=[]
        valid=type(start) is int and 0<=start<(length if circular or kind=='deletion' else length+1)
        for v in entries:
            p=v.get('position');o=v.get('original_read_position');q=v.get('phred')
            if type(p) is not int:valid=False;continue
            if kind=='insertion':
                call=v.get('read');raw=record.get('base_sequence','');qs=record.get('qualities')
                if not (p==start and v.get('reference')=='' and isinstance(call,str) and len(call)==1 and call in 'ACGTN'
                        and type(o) is int and 0<=o<len(raw) and (q is None or type(q) is int and q>=0)):
                    valid=False;continue
                expected=raw[o].translate(COMPLEMENT) if step==-1 else raw[o]
                if call!=expected or (qs and qs[o]!=q) or o in paired_original:valid=False
                original_positions.append(o);qualities.append(q)
            else:
                if not (0<=p<length and v.get('read')=='' and v.get('reference')==reference[p] and o is None and q is None and p not in mapped):valid=False
                if 0<=p<length:positions.append(p)
        if len(sequence)!=len(entries) or not valid:
            reasons.append('invalid_saved_event')
        if set(sequence)-set('ACGT'):reasons.append('ambiguous_event_bases')
        end=(start+len(entries)) if type(start) is int else None
        right_index=(start if type(start) is int else None) if kind=='insertion' else end
        left_index=start-1 if type(start) is int else None
        if circular and valid:
            left_index%=length;right_index%=length
        left=mapped.get(left_index);right=mapped.get(right_index)
        flanks={}
        for side,column in [('left',left),('right',right)]:
            if column is None:reasons.append('missing_'+side+'_flank');flanks[side]=None;continue
            p,o,ref,call,q=column
            flanks[side]={'reference_position':p,'original_read_position':o,'reference':ref,'call':call,'phred':q}
            if ref not in 'ACGT' or call!=ref:reasons.append(side+'_flank_not_reference_matching')
            if q is None or q<20:reasons.append(side+'_flank_below_q20_or_missing')
        if kind=='insertion' and any(q is None or q<20 for q in qualities):reasons.append('inserted_calls_below_q20_or_missing')
        if valid and left is not None and right is not None:
            walk=[left[1],*original_positions,right[1]] if kind=='insertion' else [left[1],right[1]]
            if any(b-a!=step for a,b in zip(walk,walk[1:])):reasons.append('complex_or_inconsistent_gap')
        shift_left=shift_right=False
        if valid and sequence and not set(sequence)-set('ACGT'):
            shift_left=left_index is not None and 0<=left_index<length and reference[left_index]==sequence[-1]
            shift_right=right_index is not None and 0<=right_index<length and reference[right_index]==sequence[0]
            if shift_left or shift_right:reasons.append('repeat_shift_possible')
        segments=[]
        for p in positions:
            if segments and segments[-1][1]==p:segments[-1][1]=p+1
            else:segments.append([p,p+1])
        events.append({'read_id':record['id'],'filename':record['original_filename'],'alignment_id':record['alignment_id'],
                       'run_number':record['run_number'],'event_number':number,'kind':kind,'boundary':start,
                       'length_bp':len(entries),'sequence':sequence,'reference_segments':segments,
                       'original_read_positions':original_positions,'inserted_phred':qualities,
                       'saved_variant_indices':[i for i,_ in group],'flanks':flanks,
                       'eligible_for_exact_anchor_review':not reasons,'review_reasons':sorted(set(reasons)),
                       'equivalent_one_base_shift':{'left':bool(shift_left),'right':bool(shift_right)},
                       'normalization':'not_performed','cross_read_comparison':'not_performed'})
    # Duplicated saved gap entries or inconsistent totals must not become two plausible events.
    insert_positions=[p for e in events if e['kind']=='insertion' for p in e['original_read_positions']]
    delete_positions=[p for e in events if e['kind']=='deletion' for a,b in e['reference_segments'] for p in range(a,b)]
    inconsistent=len(set(insert_positions))!=len(insert_positions) or len(set(delete_positions))!=len(delete_positions)
    for kind,field in [('insertion','inserted_bases'),('deletion','deleted_bases')]:
        if field in report and sum(e['length_bp'] for e in events if e['kind']==kind)!=report[field]:inconsistent=True
    if inconsistent:
        for event in events:
            event['eligible_for_exact_anchor_review']=False
            event['review_reasons']=sorted(set([*event['review_reasons'],'inconsistent_event_totals_or_duplicates']))
    return events
