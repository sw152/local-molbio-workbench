"""Bounded exact-site pair geometry on one reference; not a PCR success prediction."""


def review_pair_products(left: dict, right: dict, size_min: int, size_max: int, limit: int = 100) -> dict:
    """Count distinct reference start/span intervals using one of each oligo.

    Both oligo assignments are considered. Sites must face inward and not overlap.
    Circular products traverse the reference at most once. Missing site evidence
    makes counts lower bounds, never proof that there are no alternatives.
    """
    forward, reverse = left['reference_sites'], right['reference_sites']
    length, topology = forward['reference_length'], forward['reference_topology']
    if length < 1 or topology not in {'linear', 'circular', 'unknown'}:
        raise ValueError('Invalid reference metadata for product review')
    if (reverse['reference_length'], reverse['reference_topology']) != (length, topology):
        raise ValueError('Primer site evidence refers to different references')
    if not 1 <= size_min <= size_max <= 10000 or not 1 <= limit <= 1000:
        raise ValueError('Invalid product review size range or stored-product limit')
    for item in (left, right):
        sites = item['reference_sites']['sites']
        if len(sites) > 1000 or not 1 <= len(item['sequence']) <= length:
            raise ValueError('Product review exceeds stored-site or oligo-length limit')
        for site in sites:
            if not 0 <= site['start'] < length or site['strand'] not in {1, -1}:
                raise ValueError('Invalid binding site in product review')
            if topology != 'circular' and site['start'] + len(item['sequence']) > length:
                raise ValueError('Linear binding site exceeds reference boundary')
    products = {}
    for plus_name, plus, minus_name, minus in (
        ('forward', left, 'reverse', right), ('reverse', right, 'forward', left),
    ):
        for positive in plus['reference_sites']['sites']:
            if positive['strand'] != 1:
                continue
            for negative in minus['reference_sites']['sites']:
                if negative['strand'] != -1:
                    continue
                start, stop_site = positive['start'], negative['start']
                distance = (stop_site - start) % length if topology == 'circular' else stop_site - start
                span = distance + len(minus['sequence'])
                if distance < len(plus['sequence']) or not size_min <= span <= size_max:
                    continue
                if topology == 'circular' and span > length:
                    continue
                wraps = start + span > length
                end = start + span - length if wraps else start + span
                intended = (plus_name == 'forward' and start == left['binding_start']
                            and stop_site == right['binding_start']
                            and span == right['binding_end'] - left['binding_start'])
                key = (start, span)
                if key not in products:
                    products[key] = {'start':start, 'end':end, 'length_bp':span, 'wraps_origin':wraps,
                                     'segments':[[start,length],[0,end]] if wraps else [[start,end]],
                                     'is_intended':False, 'assignments':[]}
                product = products[key]
                product['is_intended'] |= intended
                assignment = {'plus_primer':plus_name, 'minus_primer':minus_name,
                              'plus_binding_start':start, 'minus_binding_start':stop_site}
                if assignment not in product['assignments']:
                    product['assignments'].append(assignment)
    ordered = sorted(products.values(), key=lambda p:(not p['is_intended'],p['length_bp'],p['start']))
    complete = not forward['truncated'] and not reverse['truncated']
    alternatives = sum(not p['is_intended'] for p in ordered)
    return {'method':'exact_opposed_site_pairs_v1', 'scope':'current_reference_only',
            'reference_length':length, 'reference_topology':topology,
            'coordinate_system':'zero-based-half-open', 'product_size_min':size_min, 'product_size_max':size_max,
            'origin_checked':topology == 'circular', 'maximum_reference_traversals':1,
            'overlapping_sites_evaluated':False, 'same_oligo_pairs_evaluated':False,
            'mismatches_evaluated':False, 'site_search_complete':complete,
            'observed_product_count':len(ordered), 'observed_alternative_count':alternatives,
            'total_product_count':len(ordered) if complete else None,
            'alternative_product_count':alternatives if complete else None,
            'intended_product_found':any(p['is_intended'] for p in ordered),
            'products':ordered[:limit], 'stored_product_limit':limit,
            'product_list_truncated':len(ordered)>limit}
