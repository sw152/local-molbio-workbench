"""Lossless mapping of aligned base pairs to original AB1 calls and qualities."""


def intervals(positions):
    result=[]
    for p in sorted(set(positions)):
        if result and result[-1][1]==p:
            result[-1][1]=p+1
        else:
            result.append([p,p+1])
    return result


def build_base_mapping(columns, direction, threshold=20):
    """Columns: reference position, original read position, ref base, oriented call, Q.

    Gaps have no paired bases and remain in the report's per-base variants.
    Blocks split at origin crossings and either coordinate discontinuity.
    """
    step=-1 if direction=='reverse' else 1
    blocks=[]
    counts=dict(matching=0,differing=0,low_quality=0,quality_unavailable=0,ambiguous=0)
    callable_positions=[]
    matching_positions=[]
    differing_positions=[]
    for pos,original,ref,call,quality in columns:
        if blocks and blocks[-1]['reference_end']==pos and blocks[-1]['original_read_start']+len(blocks[-1]['read_bases'])*step==original:
            block=blocks[-1]
        else:
            block={'reference_start':pos,'reference_end':pos,'original_read_start':original,
                   'original_read_step':step,'reference_bases':[],'read_bases':[],'phred':[]}
            blocks.append(block)
        block['reference_end']=pos+1
        block['reference_bases'].append(ref)
        block['read_bases'].append(call)
        block['phred'].append(quality)
        if ref not in 'ACGT' or call not in 'ACGT':
            counts['ambiguous']+=1
        elif quality is None:
            counts['quality_unavailable']+=1
        elif quality<threshold:
            counts['low_quality']+=1
        else:
            callable_positions.append(pos)
            if ref==call:
                counts['matching']+=1
                matching_positions.append(pos)
            else:
                counts['differing']+=1
                differing_positions.append(pos)
    for block in blocks:
        block['reference_bases']=''.join(block['reference_bases'])
        block['read_bases']=''.join(block['read_bases'])
    return {'schema_version':1,'method':'aligned_base_pairs_original_read_v1',
            'coordinate_system':'zero-based; reference intervals half-open',
            'read_bases_orientation':'reference-forward; reverse-complemented for reverse reads',
            'quality_threshold':threshold,'blocks':blocks,'aligned_base_count':len(columns),
            'categories':counts,'callable_base_count':len(callable_positions),
            'callable_reference_intervals':intervals(callable_positions),
            'matching_reference_intervals':intervals(matching_positions),
            'differing_reference_intervals':intervals(differing_positions),
            'placement_eligibility':'see_report_review_flags',
            'includes_insertions_or_deletions':False,'whole_reference_verified':False}
