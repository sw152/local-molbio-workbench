"""Portable, complete reports of queue-published local alignment evidence."""
from datetime import datetime, timezone
from html import escape
import json
from .alignment_coverage import summarize, regions

LIMITATIONS = [
    'Local alignment evidence only. No consensus or whole-plasmid verification.',
    'Base quality, read pairing, uniqueness and independent molecules were not assessed.',
    'Paired positions include mismatches and N; deletions are separate and can overlap other paired evidence.',
    'Reported candidates are not exhaustive. Circular MAPQ is not recalibrated. Indels are not normalized.',
    'Withheld-only and unreported placements cannot be located. Empty intervals do not prove biological absence.',
    'This is a stored historical result. Original files were not rechecked at export; no realignment was performed.',
    'All task records and review intervals are included, regardless of the current UI page or filter.',
    'Original sequences, FASTQ qualities and files are not included. Labels, identifiers and hashes are included; labels reflect export-time metadata.',
]
FALSE_FLAGS = ('base_quality_used','pairing_used','candidate_search_exhaustive',
               'circular_mapq_recalibrated','consensus_performed','whole_reference_verified')
HIT_FIELDS = ('query_start','query_end','strand','reference_intervals','reference_span_bp','crosses_origin',
              'cigar','exact_matches','mismatches','ambiguous_pairs','inserted_bases','deleted_bases',
              'alignment_columns','reported_paf_block_length','ambiguous_gap_bases',
              'local_exact_identity','aligned_query_fraction')


def pick(value, fields):
    return {k:value[k] for k in fields}


def build_report(row, evidence, labels):
    """Caller checks publication and sources; independently check identity and geometry here."""
    manifest = json.loads(row['input_manifest_json'])
    if (evidence['scope'] != 'registered_fastq_local_alignments_only'
            or evidence['job_id'] != row['id'] or evidence['reference'] != manifest['reference']
            or evidence['input_identities'] != manifest['inputs']['fastq_inputs']
            or evidence['tool'] != manifest['inputs']['tool']
            or evidence['reference_sha256'] != manifest['reference']['sequence_sha256']
            or evidence['topology'] != manifest['reference']['topology']
            or evidence['data_type'] != json.loads(row['parameters_json'])['data_type']
            or evidence['coordinate_convention'] != 'zero_based_half_open'
            or evidence['analysis_performed'] is not True
            or any(evidence[k] is not False for k in FALSE_FLAGS)
            or type(evidence['attempt_number']) is not int
            or evidence['attempt_number'] != row['attempt_count']):
        raise ValueError('alignment_report_identity_or_scope_invalid')
    coverage = summarize(evidence)
    review = {kind:regions(coverage,kind) for kind in ('unpaired','ambiguous_only','deletion')}
    coverage.pop('segments')
    records = []
    for read,source in zip(evidence['reads'],evidence['sources']):
        if read['uniqueness_assessed'] is not False:
            raise ValueError('alignment_report_identity_or_scope_invalid')
        hits = []
        for hit in read['alignments']:
            hits.append({**pick(hit,HIT_FIELDS),
                         'paired_blocks':[pick(b,('query_start','query_end','query_step','reference_intervals')) for b in hit['paired_blocks']],
                         'raw_hits':[pick(h,('target_start','target_end','type','mapq')) for h in hit['raw_hits']]})
        records.append({**pick(read,('read_index','length_bp','status','uniqueness_assessed','requires_review')),
                        'source':pick(source,('query_name','input_id','record_ordinal','sequence_sha256')),
                        'alignments':hits,
                        'withheld':[pick(w,('reason','query_start','query_end')) for w in read['withheld']]})
    return dict(schema='localmolbio.alignment-review-report',schema_version=1,
                exported_at=datetime.now(timezone.utc).isoformat(),
                job={**pick(row,('id','created_at','completed_at','status')),'attempt_number':evidence['attempt_number']},
                reference={**pick(evidence['reference'],('id','sequence_sha256','topology')),'length_bp':evidence['reference_length_bp']},
                inputs=[{**pick(i,('id','sha256','size_bytes','compression','quality_encoding','summary_sha256')),
                         'label_at_export':labels.get(i['id'])} for i in evidence['input_identities']],
                tool=pick(evidence['tool'],('name','version','binary_sha256')),
                evidence={**pick(evidence,('schema','schema_version','scope','data_type','coordinate_convention',
                                          'analysis_performed','secondary_limit','sources_sha256',*FALSE_FLAGS)),
                          'limits':pick(evidence['limits'],('reference_bases','reads','read_bases','total_read_bases',
                                                          'paf_bytes','log_bytes','process_timeout_seconds'))},
                coverage=coverage,review_regions=review,reads=records,
                source_verification_at_export='not_performed',selection='all_task_records',limitations=LIMITATIONS.copy())


def text(value):
    return escape(str(value)) if value is not None else 'Not recorded'


def spans(intervals):
    return ' → '.join(f'[{text(a)}, {text(b)})' for a,b in intervals) or 'None reported'


def render_html(report):
    c=report['coverage'];length=report['reference']['length_bp']
    cards=[('Single reported hit evidence',c['single_reported_alignment_bases']),
           ('Alternative-only evidence',c['ambiguous_only_bases']),('No paired evidence',c['unpaired_bases'])]
    metrics=''.join(f'<div><span>{label}</span><strong>{count:,} bp</strong></div>' for label,count in cards)
    bar=''.join(f'<span class="band-{i}" style="width:{count/length*100}%"></span>' for i,(_,count) in enumerate(cards))
    identities=''.join(f'<article><h3>{text(i["label_at_export"] or i["id"])}</h3><p>Input ID <code>{text(i["id"])}</code></p><p>Original SHA-256 <code>{text(i["sha256"])}</code><br>Registration summary SHA-256 <code>{text(i["summary_sha256"])}</code></p><p>{i["size_bytes"]:,} bytes · {text(i["compression"])} · {text(i["quality_encoding"])}</p></article>' for i in report['inputs'])
    region_sections=''.join(f'<article><h3>{label}</h3><p>{len(report["review_regions"][kind]):,} intervals · 0-based, half-open</p><div class="intervals">'+''.join(f'<code>[{a}, {b}) · {b-a:,} bp</code>' for a,b in report['review_regions'][kind])+'</div>'+('' if report['review_regions'][kind] else '<p>None reported. This is not a plasmid pass verdict.</p>')+'</article>' for kind,label in [('unpaired','No paired evidence'),('ambiguous_only','Alternative-only evidence'),('deletion','Reported deletion positions')])
    labels={i['id']:i['label_at_export'] or i['id'] for i in report['inputs']}
    records=[]
    for read in report['reads']:
        source=read['source'];hits=[]
        for idx,h in enumerate(read['alignments'],1):
            blocks=''.join(f'<span style="left:{a/length*100}%;width:{(b-a)/length*100}%"></span>' for block in h['paired_blocks'] for a,b in block['reference_intervals'])
            hits.append(f'<section class="hit"><h4>Hit {idx} · {"← Reverse" if h["strand"]=="-" else "Forward →"}{" · crosses origin" if h["crosses_origin"] else ""}</h4><div class="track" aria-label="Paired reference positions">{blocks}</div><div class="ticks"><span>0 bp</span><span>{length:,} bp</span></div><p><b>Reference</b> {spans(h["reference_intervals"])}<br><b>Original read</b> [{h["query_start"]}, {h["query_end"]}) · coordinates {"decrease" if h["strand"]=="-" else "increase"} along reference</p><p>{h["exact_matches"]:,} exact · {h["mismatches"]:,} mismatches · {h["ambiguous_pairs"]:,} ambiguous pairs · {h["inserted_bases"]:,} inserted / {h["deleted_bases"]:,} deleted<br>{h["local_exact_identity"]*100:.2f}% local exact identity · {h["aligned_query_fraction"]*100:.2f}% read aligned</p><p>CIGAR <code>{text(h["cigar"])}</code></p><p>Reported raw hits: '+ '; '.join(f'{text(raw["type"])} / MAPQ {text(raw["mapq"])}' for raw in h['raw_hits'])+'</p></section>')
        withheld=''.join(f'<p class="warning">Withheld: {text(w["reason"].replace("_"," "))} · original read [{text(w["query_start"])}, {text(w["query_end"])})</p>' for w in read['withheld'])
        records.append(f'<article class="record"><p class="eyebrow">Record {source["record_ordinal"]} · {text(source["query_name"])}</p><h3>{text(labels[source["input_id"]])}</h3><p>{read["length_bp"]:,} bases · reported candidates: {len(read["alignments"])} · uniqueness not assessed</p><p>Input ID <code>{text(source["input_id"])}</code><br>Sequence SHA-256 <code>{text(source["sequence_sha256"])}</code></p>{withheld}'+(''.join(hits) or '<p>No usable local alignment reported.</p>')+'</article>')
    return '''<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><meta http-equiv="Content-Security-Policy" content="default-src 'none'; style-src 'unsafe-inline'; base-uri 'none'; form-action 'none'"><title>Local alignment review</title><style>
*{box-sizing:border-box}body{margin:0;background:#edf2f5;color:#213548;font:15px/1.65 system-ui,sans-serif}main{max-width:1000px;margin:40px auto;padding:36px;background:white;border:1px solid #cdd9e1;border-radius:18px}h1{font-size:34px;line-height:1.2;margin:10px 0 20px}h2{margin-top:38px;font-size:22px}h3{margin:0 0 8px;font-size:16px}h4{margin:10px 0}p{margin:8px 0}.eyebrow{font-size:11px;letter-spacing:.14em;text-transform:uppercase;color:#527386}.warning{border-left:4px solid #b77d37;background:#fff6e8;padding:14px 18px}.cards{display:grid;grid-template-columns:repeat(3,1fr);gap:14px;margin:18px 0}.cards>div{padding:15px;background:#f1f6f8;border-radius:8px}.cards span{font-size:12px}.cards strong{display:block;font-size:24px}.bar{display:flex;height:12px;border-radius:6px;overflow:hidden}.band-0{background:#619ab0}.band-1{background:#c79e62}.band-2{background:#7a8698}article{border:1px solid #d4e0e8;border-radius:10px;padding:20px;margin:14px 0}code{font:12px/1.7 ui-monospace,monospace;overflow-wrap:anywhere}p,h3{overflow-wrap:anywhere}.intervals{display:flex;gap:8px;flex-wrap:wrap}.intervals code{background:#edf3f6;padding:5px 9px;border-radius:4px}.hit{margin-top:20px;border-top:1px solid #d4e0e8;padding-top:8px}.track{height:16px;position:relative;background:#e0e8ed;margin-top:15px}.track span{position:absolute;top:2px;height:12px;background:#429991}.ticks{display:flex;justify-content:space-between;font-size:11px;color:#607782}.note{font-size:12px;color:#567080}li{margin-bottom:6px}@media(max-width:600px){main{margin:0;padding:20px 16px;border:0;border-radius:0}h1{font-size:28px}.cards{grid-template-columns:1fr}.cards>div{display:flex;justify-content:space-between;align-items:center}.cards strong{font-size:20px}article{padding:14px}}@media print{body{background:white;font-size:10pt}main{margin:0;padding:0;border:0}h2,h3,h4{break-after:avoid}.cards,.hit,.warning{break-inside:avoid}.track,.bar{print-color-adjust:exact}code{font-size:8pt}}
</style></head><body><main>'''+f'''<header><p class="eyebrow">Local MolBio Workbench · saved evidence</p><h1>Local alignment review</h1><p class="warning"><b>No whole-plasmid verdict</b><br>Quality, pairing and uniqueness unassessed. No consensus.</p><p>Complete task export · {len(report['reads']):,} records · {length:,} bp {text(report['reference']['topology'])} reference</p><p class="note">Exported {text(report['exported_at'])} · schema {text(report['schema'])} v1</p></header><h2>Paired-position overview</h2><p class="note">Base-count composition, not a spatial depth track. Paired bases include mismatches and N.</p><div class="bar" aria-hidden="true">{bar}</div><div class="cards">{metrics}</div><p>{c['deletion_evidence_bases']:,} positions with deletion evidence · records with withheld evidence: {c['read_counts']['withheld']:,} · records without reported alignments: {c['read_counts']['no_alignment_reported']:,}</p><p class="note">Deletion positions can overlap paired positions from other records. Single reported hits do not establish uniqueness.</p><h2>Review regions</h2>{region_sections}<h2>Run identity</h2><p>Job <code>{text(report['job']['id'])}</code> · published attempt {report['job']['attempt_number']}<br>Completed {text(report['job']['completed_at'])}<br>Data type: {text(report['evidence']['data_type'])}</p><p>Reference revision <code>{text(report['reference']['id'])}</code><br>Reference SHA-256 <code>{text(report['reference']['sequence_sha256'])}</code></p><p>{text(report['tool']['name'])} {text(report['tool']['version'])}<br>Binary SHA-256 <code>{text(report['tool']['binary_sha256'])}</code><br>Source mapping SHA-256 <code>{text(report['evidence']['sources_sha256'])}</code></p><h2>Input identities</h2>{identities}<h2>All original records</h2><p class="note">Coordinates are 0-based, half-open; file record numbers are 1-based. Colored tracks show paired reference positions; deletions remain gaps. All candidates are listed, with original read coordinates and direction.</p>{''.join(records)}<h2>Interpretation and limitations</h2><ul>{''.join('<li>'+text(s)+'</li>' for s in report['limitations'])}</ul><p class="note">For structured paired blocks and complete machine-readable evidence, also save the JSON report. No external assets or scripts are required to view this document.</p></main></body></html>'''
