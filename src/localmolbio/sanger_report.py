"""Portable saved-analysis reports. Export never realigns or substitutes a newer run."""
from datetime import datetime, timezone
from html import escape
import json
from .database import connect

LIMITATIONS = [
    'This is local read-alignment evidence, not a whole-plasmid verification verdict.',
    'Uncovered reference positions remain unassessed; identity describes only the aligned region.',
    'Differences are per-base and indels are not normalized. Review ambiguity and quality flags.',
    'End trimming, when enabled, retains internal low-quality calls; retained-read coverage is not original-read coverage.',
    'Export reproduces stored evidence. Original files and their hashes were not rechecked at export.',
    'Raw AB1 signals and full read/reference sequences are not included in this report.',
]


def load_analysis_export(read_id: str, analysis_id: str) -> dict | None:
    with connect() as db:
        db.execute('BEGIN')
        run = db.execute('SELECT * FROM sanger_analysis_runs WHERE id=? AND sequencing_read_id=?', (analysis_id, read_id)).fetchone()
        if run is None:
            return None
        job = db.execute('SELECT * FROM analysis_jobs WHERE id=?', (run['analysis_job_id'],)).fetchone()
        read = db.execute('SELECT original_filename,length_bp FROM sequencing_reads WHERE id=?', (read_id,)).fetchone()
        reference = db.execute('SELECT label,length(sequence_text) AS length_bp,topology FROM sequence_revisions WHERE id=?', (run['sequence_revision_id'],)).fetchone()
    manifest = json.loads(job['input_manifest_json'])
    def metadata(prefix, field, fallback):
        key = prefix + '_' + field
        return {'value':manifest[key] if key in manifest else fallback,
                'source':'analysis_input_snapshot' if key in manifest else 'linked_record_metadata'}
    return {
        'schema':'localmolbio.sanger-analysis-report', 'schema_version':1,
        'exported_at':datetime.now(timezone.utc).isoformat(),
        'analysis':{'id':run['id'],'job_id':run['analysis_job_id'],'read_id':read_id,
                    'sequence_revision_id':run['sequence_revision_id'],'run_number':run['run_number'],
                    'previous_alignment_id':run['previous_alignment_id'],'analyzed_at':run['created_at'],
                    'job_status':job['status']},
        'inputs':{
            'read':{'sha256':manifest.get('read_sha256'),
                    'filename':metadata('read','filename',read['original_filename']),
                    'length_bp':metadata('read','length_bp',read['length_bp'])},
            'reference':{'sha256':manifest.get('reference_sha256'),
                         'label':metadata('reference','label',reference['label']),
                         'length_bp':metadata('reference','length_bp',reference['length_bp']),
                         'topology':metadata('reference','topology',reference['topology'])}},
        'parameters':json.loads(job['parameters_json']),
        'saved_alignment':json.loads(run['report_json']),
        'source_verification_at_export':'not_performed',
        'limitations':LIMITATIONS.copy(),
    }


def _text(value):
    return escape(str(value)) if value is not None else 'Not recorded'


def _percent(value):
    return f'{value * 100:.1f}%' if isinstance(value,(int,float)) else 'Not recorded'


def _bar(intervals, length, label, excluded=False):
    if not isinstance(length,int) or length < 1 or not isinstance(intervals,list):
        return '<p class="note">Coordinate graphic unavailable for this saved analysis.</p>'
    if any(not isinstance(s,list) or len(s)!=2 or not all(isinstance(n,int) for n in s) or not 0<=s[0]<s[1]<=length for s in intervals):
        return '<p class="warning">Saved coordinates cannot be rendered safely; inspect the saved evidence below.</p>'
    blocks=''.join(f'<rect x="{a/length*1000}" width="{(b-a)/length*1000}" height="20" fill="#28aa9b"/>' for a,b in intervals)
    return f'<svg class="bar" viewBox="0 0 1000 20" preserveAspectRatio="none" role="img" aria-label="{_text(label)}"><rect width="1000" height="20" fill="{"#c58d46" if excluded else "#dde5ee"}"/>{blocks}</svg><div class="axis"><span>Original base 1</span><span>{length:,} bp</span></div>'


def render_analysis_html(payload: dict) -> str:
    analysis, inputs, saved = payload['analysis'], payload['inputs'], payload['saved_alignment']
    evidence=saved.get('evidence') or {}
    flags=evidence.get('review_flags') or []
    variants=saved.get('variants') or []
    trim=evidence.get('trimming')
    reference_length=inputs['reference']['length_bp']['value']
    coverage=_bar(evidence.get('reference_covered_intervals'), reference_length, 'Aligned reference positions')
    coverage=coverage.replace('Original base 1','Reference base 1')
    if trim:
        trim_html=f'<p><b>{"Trimming off" if trim["method"]=="none" else "End calls below Q"+str(trim["quality_threshold"])+" excluded"}</b> · {_text(trim["retained_length"])} / {_text(trim["original_length"])} original bases retained.</p>'
        trim_html+=_bar([[trim['original_read_start'],trim['original_read_end']]],trim['original_length'],'Retained original-read interval',excluded=True)
        trim_html+=f'<p>Original bases {trim["original_read_start"]+1}–{trim["original_read_end"]}; {trim["removed_left"]} left / {trim["removed_right"]} right excluded. Internal low-quality calls remain.</p>'
    else:
        trim_html='<p class="note">Retained-read interval was not recorded for this analysis.</p>'
    variant_rows=''.join('<tr>'+''.join(f'<td>{cell}</td>' for cell in (
        f'Boundary after {_text(v.get("position"))}' if v.get('kind')=='insertion' else _text(v['position']+1 if v.get('position') is not None else None),
        _text(v.get('kind')), _text(v.get('reference') or '–')+' → '+_text(v.get('read') or '–'),
        _text(v['original_read_position']+1 if v.get('original_read_position') is not None else None),
        _text(v.get('phred'))))+'</tr>' for v in variants)
    if not variant_rows:
        variant_rows='<tr><td colspan="5">No saved differences in the aligned region; uncovered positions remain unassessed.</td></tr>'
    flag_html='<ul>'+''.join(f'<li>{_text(flag.replace("_"," "))}</li>' for flag in flags)+'</ul>' if flags else '<p>No review flags recorded. This is not a validation verdict.</p>'
    if not evidence:
        flag_html='<p class="warning">Detailed evidence was not recorded for this legacy analysis. Coverage and quality must not be inferred from identity.</p>'
    input_rows=''
    for kind, fields in inputs.items():
        for name, value in fields.items():
            if isinstance(value,dict):
                text=f'{_text(value["value"])} <small>{_text(value["source"].replace("_"," "))}</small>'
            else:
                text=f'<code>{_text(value)}</code>'
            input_rows+=f'<tr><th>{_text(kind)} · {_text(name)}</th><td>{text}</td></tr>'
    metrics=''.join(f'<div><span>{title}</span><b>{_percent(value)}</b></div>' for title,value in (
        ('Aligned identity',saved.get('identity_fraction')),('Original read aligned',evidence.get('read_aligned_fraction')),
        ('Retained read aligned',evidence.get('retained_read_aligned_fraction')),('Reference covered',evidence.get('reference_covered_fraction'))))
    return f'''<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><meta http-equiv="Content-Security-Policy" content="default-src 'none'; style-src 'unsafe-inline'; base-uri 'none'; form-action 'none'"><title>Sanger evidence · run {analysis['run_number']}</title><style>
*{{box-sizing:border-box}}body{{margin:0;background:#edf2f7;color:#1b2b40;font:15px/1.55 system-ui,sans-serif}}main{{max-width:980px;margin:30px auto;padding:36px;background:white;border:1px solid #d5e0eb;border-radius:18px}}h1{{font-size:32px;letter-spacing:-.04em;margin:8px 0}}h2{{font-size:19px;margin:0 0 12px}}section{{margin-top:28px;padding-top:22px;border-top:1px solid #d5e0eb}}.eyebrow{{font-size:12px;letter-spacing:.12em;text-transform:uppercase;color:#526b83}}.note,small{{color:#5c7084}}small{{display:block;font-size:11px}}.warning{{padding:12px;border:1px solid #dfbd83;background:#fff8e9;border-radius:8px}}.metrics{{display:grid;grid-template-columns:repeat(4,1fr);gap:12px;margin:24px 0}}.metrics div{{padding:16px;background:#f0f5f9;border-radius:10px}}.metrics span{{display:block;color:#536c83;font-size:11px;text-transform:uppercase}}.metrics b{{display:block;font-size:25px;margin-top:8px}}.bar{{display:block;width:100%;height:18px;margin:14px 0 5px}}.axis{{display:flex;justify-content:space-between;font-size:11px;color:#587087}}.scroll{{overflow:auto}}table{{border-collapse:collapse;width:100%;font-size:13px}}th,td{{padding:10px;text-align:left;border-bottom:1px solid #dbe4ec;vertical-align:top}}th{{color:#526b83;font-weight:500}}.variants{{min-width:570px}}code,pre{{font-family:ui-monospace,monospace;font-size:11px;overflow-wrap:anywhere}}pre{{white-space:pre-wrap;background:#f1f5f9;padding:16px;border-radius:8px}}summary{{cursor:pointer;font-weight:600}}.identifiers{{overflow-wrap:anywhere;font-size:12px}}@media(max-width:620px){{main{{margin:0;padding:20px;border:0;border-radius:0}}.metrics{{grid-template-columns:repeat(2,1fr)}}h1{{font-size:28px}}th,td{{padding:8px}}}}@media print{{body{{background:white}}main{{max-width:none;border:0;padding:0;margin:0}}section{{break-inside:avoid}}.scroll{{overflow:visible}}.variants{{min-width:0}}th,td{{overflow-wrap:anywhere}}}}
</style></head><body><main><header><span class="eyebrow">MolBio Workbench · Saved analysis evidence</span><h1>Sanger analysis report</h1><p><b>Run {analysis['run_number']}</b> · {_text(evidence.get('direction'))} orientation · analyzed {_text(analysis['analyzed_at'])}</p><p class="identifiers">Analysis {_text(analysis['id'])}<br>Job {_text(analysis['job_id'])}<br>Reference revision {_text(analysis['sequence_revision_id'])}<br>Read {_text(analysis['read_id'])}<br>Previous analysis {_text(analysis['previous_alignment_id'])}</p><p class="warning">Local read evidence. This report does not certify the whole plasmid.</p></header><div class="metrics">{metrics}</div><section><h2>Input provenance</h2><table>{input_rows}</table><p class="note">Hashes are the values recorded in this analysis job; missing hashes remain “Not recorded”. Linked-record metadata is labeled separately from analysis-time snapshots. Files were not rechecked at export.</p></section><section><h2>Retained original-read interval</h2>{trim_html}</section><section><h2>Reference coverage</h2>{coverage}<p>Saved covered intervals (1-based inclusive): {_text('; '.join(f'{a+1}–{b}' for a,b in evidence.get('reference_covered_intervals',[]))) if evidence.get('reference_covered_intervals') is not None else 'Not recorded'}.</p></section><section><h2>Review flags</h2>{flag_html}</section><section><h2>Differences · {len(variants)} saved calls</h2><div class="scroll"><table class="variants"><thead><tr><th>Reference position¹</th><th>Type</th><th>Reference → read</th><th>Original read base</th><th>Phred</th></tr></thead><tbody>{variant_rows}</tbody></table></div><p class="note">¹ Base positions are 1-based. Insertions use boundaries after the indicated reference base; boundary 0 is before base 1. Read positions refer to the original uploaded orientation, including reverse alignments. Deletions have no read base or per-call Phred.</p></section><section><h2>Saved analysis parameters</h2><pre>{_text(json.dumps(payload['parameters'],indent=2,ensure_ascii=False))}</pre><details><summary>Full saved alignment evidence · raw zero-based coordinates</summary><pre>{_text(json.dumps(saved,indent=2,ensure_ascii=False))}</pre></details></section><section><h2>Interpretation limits</h2><ul>{''.join(f'<li>{_text(item)}</li>' for item in payload['limitations'])}</ul><p class="note">Exported {_text(payload['exported_at'])} · schema {payload['schema_version']} · self-contained report, no external resources.</p></section></main></body></html>'''
