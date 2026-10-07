/* Revision-scoped read evidence. Pending responses never update a different construct. */
(() => {
  const el = id => document.getElementById(id);
  const escape = value => String(value ?? '').replace(/[&<>"']/g, c => ({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[c]));
  const percent = value => Number.isFinite(value) ? `${(value * 100).toFixed(1)}%` : 'Unavailable';
  function timestamp(value) {
    const date=new Date(value);
    return Number.isNaN(date.getTime())?'Time not recorded':date.toISOString().replace('T',' ').replace(/\.\d{3}Z$/,' UTC');
  }
  const flags = {
    ambiguous_alignment: 'Multiple equally scoring alignments — placement needs review',
    candidate_search_truncated: 'Repeated sequence: search limit reached',
    partial_read_alignment: 'Less than 80% of the original read aligned',
    ambiguous_bases: 'Ambiguous base calls in the alignment',
    quality_unavailable: 'Per-base quality is unavailable',
    low_quality_aligned_bases: 'Aligned bases include calls below Q20',
    legacy_evidence_unavailable: 'Older analysis: coverage and quality evidence are unavailable',
  };
  let state = null, generation = 0, busy = false;
  function status(message, error = false) {
    el('sanger-status').textContent = message;
    el('sanger-status').className = error ? 'sanger-error' : '';
  }
  async function request(url, options) {
    const response = await fetch(url, options);
    const data = await response.json();
    if (!response.ok) throw new Error(typeof data.detail === 'string' ? data.detail : `Request failed (${response.status})`);
    return data;
  }
  function controls() {
    el('sanger-upload-button').disabled = busy || !state;
    el('sanger-file').disabled = busy || !state;
    el('sanger-direction').disabled = busy || !state;
    el('sanger-reads').querySelectorAll('[data-report-export], [data-analyze], [data-analysis-direction], [data-analysis-trim], [data-history-page], [data-view-run], [data-latest-run]').forEach(b => b.disabled = busy);
  }
  function clear() {
    generation++;
    state = null;
    busy = false;
    el('sanger-form').reset();
    el('sanger-reads').replaceChildren();
    window.sangerSummary.clear();
    status('');
    controls();
  }
  function coverage(evidence, length) {
    const intervals = evidence.reference_covered_intervals;
    if (!Array.isArray(intervals)) return '<p class="evidence-note">Coverage unavailable for this analysis.</p>';
    const blocks = intervals.map(([start, end]) => {
      const x = Math.max(0, Math.min(1, start / length)) * 1000;
      const width = Math.max(0, Math.min(1000 - x, (end - start) / length * 1000));
      return `<rect x="${x}" y="4" width="${width}" height="12" rx="2" fill="#55e5d7"><title>Reference bases ${start + 1}–${end}</title></rect>`;
    }).join('');
    return `<div class="coverage-heading"><span>Aligned reference positions</span><b>${percent(evidence.reference_covered_fraction)}</b></div><svg class="coverage-bar" viewBox="0 0 1000 20" preserveAspectRatio="none" role="img" aria-label="${percent(evidence.reference_covered_fraction)} of reference positions have aligned read bases"><rect y="4" width="1000" height="12" rx="3" fill="#29344b"/>${blocks}</svg><div class="coverage-axis"><span>1 bp</span><span>${length.toLocaleString()} bp</span></div>`;
  }
  function differences(variants) {
    if (!variants.length) return '<p class="evidence-note">No differences in the aligned region. Uncovered positions remain unassessed.</p>';
    return `<details class="variant-details"><summary>${variants.length} difference / uncertain base call${variants.length === 1 ? '' : 's'} <span>Inspect calls</span></summary><div class="variant-scroll"><table><thead><tr><th>Reference position</th><th>Type</th><th>Reference → read</th><th>Read base¹</th><th>Phred</th></tr></thead><tbody>${variants.map(v => `<tr><td>${v.kind === 'insertion' ? `Boundary after ${Number(v.position)}` : Number(v.position) + 1}</td><td>${escape(v.kind)}</td><td><code>${escape(v.reference || '–')} → ${escape(v.read || '–')}</code></td><td>${v.original_read_position == null ? '–' : `<button class="trace-jump" data-trace-position="${Number(v.original_read_position)}" title="Inspect chromatogram at this read base">${Number(v.original_read_position) + 1} ↗</button>`}</td><td>${v.phred == null ? 'Unavailable' : Number(v.phred)}</td></tr>`).join('')}</tbody></table></div><p class="evidence-note">¹ One-based position in the original uploaded read. Insertions use a boundary after the indicated reference base; boundary 0 is before base 1.</p></details>`;
  }
  function trimOptions(read) {
    const saved=read.evidence?.parameters?.trim_quality_threshold??null;
    const levels=[...new Set([10,20,30,...(saved===null?[]:[saved])])].sort((a,b)=>a-b);
    return `<option value="" ${saved===null?'selected':''}>Off · retain original read</option>${levels.map(q=>`<option value="${Number(q)}" ${q===saved?'selected':''}>Trim end calls below Q${Number(q)}</option>`).join('')}`;
  }
  function trimmingView(evidence) {
    const t=evidence.trimming;
    if(!t)return '<p class="evidence-note">Retained-read interval was not recorded for this historical analysis.</p>';
    if(t.method==='none')return `<p class="trim-off">End trimming off · all ${t.original_length} original bases supplied for alignment.</p>`;
    const x=1000*t.original_read_start/t.original_length,width=1000*t.retained_length/t.original_length;
    return `<section class="trim-evidence" aria-label="Quality trimming evidence"><header><b>End calls below Q${t.quality_threshold} excluded</b><span>${t.retained_length} / ${t.original_length} bases retained</span></header><svg viewBox="0 0 1000 20" preserveAspectRatio="none" role="img" aria-label="Original read: ${t.removed_left} bases excluded from the left, ${t.retained_length} retained, ${t.removed_right} excluded from the right"><rect width="1000" height="20" rx="3" fill="#9b652c"/><rect x="${x}" width="${width}" height="20" fill="#55e5d7"/></svg><div class="coverage-axis"><span>Original base 1</span><span>${t.original_length}</span></div><p>Retained original bases <button class="trace-jump" data-trace-position="${t.original_read_start}">${t.original_read_start+1} ↗</button>–<button class="trace-jump" data-trace-position="${t.original_read_end-1}">${t.original_read_end} ↗</button> · excluded ${t.removed_left} left / ${t.removed_right} right.</p><p>Amber: excluded ends. Cyan: supplied for alignment. Internal low-quality calls remain. This interval is not a validation verdict.</p></section>`;
  }
  function baseMappingView(evidence, length) {
    const m=evidence.base_mapping;
    if(!m || m.schema_version!==1)return '<p class="mapping-unavailable evidence-note">Base-level quality mapping was not recorded for this analysis. A new run can add it without changing this report.</p>';
    const c=m.categories;
    const blocks=(intervals,color)=>intervals.map(([a,b])=>`<rect x="${a/length*1000}" y="4" width="${(b-a)/length*1000}" height="12" fill="${color}"><title>Reference bases ${a+1}–${b}</title></rect>`).join('');
    return `<section class="base-mapping" aria-label="Base-level quality evidence"><header><b>Base-level quality evidence</b><span>Q${m.quality_threshold}+ · A/C/G/T</span></header><div class="base-mapping-counts"><div><strong>${c.matching}</strong><span>Matching calls</span></div><div><strong>${c.differing}</strong><span>Differing calls</span></div><div><strong>${m.aligned_base_count-m.callable_base_count}</strong><span>Outside quality filter</span></div></div><svg viewBox="0 0 1000 20" preserveAspectRatio="none" role="img" aria-label="${m.callable_base_count} reference positions with unambiguous calls at Q${m.quality_threshold} or above, including ${c.differing} difference${c.differing===1?'':'s'}"><rect y="4" width="1000" height="12" rx="3" fill="#29344b"/>${blocks(m.matching_reference_intervals,'#55e5d7')}${blocks(m.differing_reference_intervals,'#edb87a')}</svg><div class="coverage-axis"><span>1 bp</span><span>${length.toLocaleString()} bp</span></div><p>Cyan: matching. Amber: differing. Grey: outside this filter or not aligned. For this reported placement only; placement warnings still apply. Gaps remain in the difference table.</p><details><summary>Inspect saved coordinate mapping</summary><p>${c.low_quality} low-quality · ${c.quality_unavailable} quality unavailable · ${c.ambiguous} ambiguous paired calls. Categories are exclusive; ambiguous calls are counted first.</p><p class="mapping-scroll-hint">Scroll the table sideways to inspect direction.</p><div class="mapping-scroll" tabindex="0" role="region" aria-label="Saved coordinate mapping; scroll to inspect all columns"><table><thead><tr><th>Reference bases</th><th>Original read bases</th><th>Direction</th></tr></thead><tbody>${m.blocks.map(b=>`<tr><td>${b.reference_start+1}–${b.reference_end}</td><td>${b.original_read_start+1} → ${b.original_read_start+1+(b.reference_end-b.reference_start-1)*b.original_read_step}</td><td>${b.original_read_step===-1?'Reverse':'Forward'}</td></tr>`).join('')}</tbody></table></div><p>One-based coordinates. Blocks split at gaps and the reference origin. JSON export retains each paired base, its original read position and Phred value. This is not a consensus or whole-plasmid verdict.</p></details></section>`;
  }
  function reportView(read, report, latest) {
    const evidence=report.evidence||{},review=evidence.review_flags||[],params=evidence.parameters||{};
    return `<div class="analysis-exports"><span>Export this displayed run</span><div><button data-report-export="html" data-analysis-id="${escape(report.alignment_id)}" data-run-number="${report.run_number}">Download HTML</button><button data-report-export="json" data-analysis-id="${escape(report.alignment_id)}" data-run-number="${report.run_number}">Download JSON</button></div></div><div class="analysis-version ${latest?'':'historical'}"><b>Run ${report.run_number} · ${latest?'Latest saved analysis':'Historical analysis'}</b><time>${escape(timestamp(report.analyzed_at))}</time><span>${escape(evidence.direction||'Unrecorded')} orientation · requested ${escape(evidence.requested_direction||params.requested_direction||'unrecorded')}</span></div><div class="read-metrics"><div><span>Original mean quality</span><strong>${read.quality_summary?.mean_phred==null?'Unavailable':`Q${Number(read.quality_summary.mean_phred).toFixed(1)}`}</strong></div><div><span>Aligned identity</span><strong>${percent(report.identity_fraction)}</strong></div><div><span>Original read aligned</span><strong>${percent(evidence.read_aligned_fraction)}</strong></div><div><span>Retained read aligned</span><strong>${percent(evidence.retained_read_aligned_fraction)}</strong></div></div>${trimmingView(evidence)}${coverage(evidence,state.length)}${baseMappingView(evidence,state.length)}${review.length?`<ul class="read-flags">${review.map(flag=>`<li>${escape(flags[flag]||flag)}</li>`).join('')}</ul>`:''}${differences(report.variants||[])}<p class="analysis-method">${escape(params.algorithm||'Method not recorded')} ${escape(params.algorithm_version||'')} · evidence ${escape(params.evidence_version||'unrecorded')} · trimming ${escape(params.trim_quality_threshold?`end calls below Q${params.trim_quality_threshold}`:params.trimming||'not recorded')}<br>Job <code title="${escape(report.job_id)}">${escape(report.job_id)}</code></p>`;
  }
  function render(reads) {
    state.reads=new Map(reads.map(r=>[r.id,r]));
    el('sanger-reads').innerHTML=reads.length?reads.map(read=>{
      const analyzed=Boolean(read.alignment_id);
      return `<article class="read-card" data-read-id="${escape(read.id)}"><header class="read-heading"><div><b>${escape(read.original_filename)}</b><small>${Number(read.length_bp).toLocaleString()} original bases · upload orientation ${escape(read.direction)}</small></div><span class="badge">${analyzed?`${read.run_number} saved run${read.run_number===1?'':'s'}`:'Ready to analyze'}</span></header><div class="analysis-report">${analyzed?reportView(read,read,true):''}</div><div class="analysis-actions"><label>Direction for new analysis<select data-analysis-direction>${['unknown','forward','reverse'].map(v=>`<option value="${v}" ${v===(read.evidence?.requested_direction||read.direction)?'selected':''}>${v==='unknown'?'Detect automatically':v==='forward'?'Forward':'Reverse complement'}</option>`).join('')}</select></label><label>End quality trimming<select data-analysis-trim>${trimOptions(read)}</select></label><button class="export" data-analyze="${escape(read.id)}">${analyzed?'Run new analysis':'Analyze against this revision'}</button></div>${analyzed?`<p class="evidence-note">A new run preserves every report and the original AB1. End trimming stops at the first call meeting the threshold on each side; internal low-quality calls remain.</p><button class="analysis-history-button" data-history-page="0">Browse saved analyses</button><div class="analysis-history" aria-live="polite"></div>`:''}<button class="trace-toggle" data-trace-open>View chromatogram</button><div class="trace-view" hidden></div></article>`;
    }).join(''):'<div class="read-empty"><span>＋</span><b>No reads attached to this revision</b><p>Add an AB1 chromatogram to review its sequence evidence.</p></div>';
    controls();
  }
  el('sanger-reads').addEventListener('click',async event=>{
    const button=event.target.closest('[data-report-export]');
    if(!button||!state||busy)return;
    const token=generation,card=button.closest('.read-card'),format=button.dataset.reportExport;
    const analysisId=button.dataset.analysisId,run=button.dataset.runNumber;
    busy=true;controls();status(`Preparing run ${run} report…`);
    try{
      const response=await fetch(`/api/sanger-reads/${encodeURIComponent(card.dataset.readId)}/analyses/${encodeURIComponent(analysisId)}/export?format=${format}`);
      if(!response.ok){let detail;try{detail=(await response.json()).detail;}catch{}throw new Error(typeof detail==='string'?detail:`Export failed (${response.status})`);}
      const blob=await response.blob();
      if(token!==generation||!card.isConnected)return;
      const url=URL.createObjectURL(blob),link=document.createElement('a');
      link.href=url;link.download=`sanger-run-${run}-${analysisId}.${format}`;
      document.body.append(link);link.click();link.remove();setTimeout(()=>URL.revokeObjectURL(url),30000);
      status(`Run ${run} ${format.toUpperCase()} report download started.`);
    }catch(error){if(token===generation)status(`Could not export report: ${error.message}`,true);}
    finally{if(token===generation){busy=false;controls();}}
  });
  el('sanger-reads').addEventListener('click',async event=>{
    const card=event.target.closest('.read-card');
    if(!card||!state||busy)return;
    const read=state.reads.get(card.dataset.readId),token=generation;
    const show=event.target.closest('[data-view-run], [data-latest-run]');
    if(show){
      const report=show.hasAttribute('data-latest-run')?read:card.analysisHistory?.get(show.dataset.viewRun);
      if(report){card.querySelector('.analysis-report').innerHTML=reportView(read,report,report.alignment_id===read.alignment_id);card.querySelectorAll('[data-view-run]').forEach(b=>b.setAttribute('aria-pressed',String(b.dataset.viewRun===report.alignment_id)));}
      return;
    }
    const page=event.target.closest('[data-history-page]');
    if(!page||page.disabled||card.historyLoading)return;
    card.historyLoading=true;page.disabled=true;
    const panel=card.querySelector('.analysis-history'),offset=Number(page.dataset.historyPage);
    panel.textContent='Loading saved analyses…';
    try{
      const result=await request(`/api/sanger-reads/${encodeURIComponent(read.id)}/analyses?limit=5&offset=${offset}`);
      if(token!==generation||!card.isConnected)return;
      card.analysisHistory=new Map(result.items.map(r=>[r.alignment_id,r]));
      panel.innerHTML=`<p>${result.total} saved analyses · newest first · showing ${result.offset+1}–${result.offset+result.items.length}</p><div class="history-runs">${result.items.map(r=>`<button data-view-run="${escape(r.alignment_id)}" aria-pressed="false">Run ${r.run_number}<small>${escape(r.evidence?.direction||'Unrecorded')} · ${escape(timestamp(r.analyzed_at))}</small></button>`).join('')}</div><nav>${offset>0?`<button data-history-page="${Math.max(0,offset-5)}">← Newer runs</button>`:''}${result.has_more?`<button data-history-page="${offset+5}">Older runs →</button>`:''}<button data-latest-run>Show latest report</button></nav>`;
    }catch(error){if(token===generation&&card.isConnected){panel.innerHTML=`<p class="sanger-error">Could not load history: ${escape(error.message)}</p><button data-history-page="${offset}">Retry history</button>`;}}
    finally{if(card.isConnected){card.historyLoading=false;page.disabled=busy;}}
  });
  async function refresh(token) {
    window.sangerSummary.clear();
    const result = await request(`/api/sequence-revisions/${encodeURIComponent(state.id)}/sanger-reads?include_summary=true`);
    if (token === generation) { render(result.reads); window.sangerSummary.render(result.summary); }
  }
  async function open(id, length) {
    clear();
    if (!id) { status('This construct has no sequence revision.', true); return; }
    state = {id, length};
    const token = generation;
    controls();
    status('Loading attached reads…');
    try { await refresh(token); if (token === generation) status(''); }
    catch (error) { if (token === generation) status(`Could not load reads: ${error.message}`, true); }
  }
  el('sanger-form').addEventListener('submit', async event => {
    event.preventDefault();
    if (!state || busy || !el('sanger-file').files.length) return;
    const file = el('sanger-file').files[0], token = generation;
    if (!file.name.toLowerCase().endsWith('.ab1')) { status('Choose an AB1 chromatogram file.', true); return; }
    const body = new FormData();
    body.append('file', file);
    body.append('sequence_revision_id', state.id);
    body.append('direction', el('sanger-direction').value);
    busy = true; controls(); status(`Attaching ${file.name}…`);
    try {
      await request('/api/sanger-reads', {method:'POST', body});
      if (token !== generation) return;
      el('sanger-form').reset();
      await refresh(token);
      if (token === generation) status('Read attached. Analyze it against this revision to inspect differences.');
    } catch (error) { if (token === generation) status(`Could not attach read: ${error.message}`, true); }
    finally { if (token === generation) { busy = false; controls(); } }
  });
  el('sanger-reads').addEventListener('click', async event => {
    const button = event.target.closest('[data-analyze]');
    if (!button || !state || busy) return;
    const token = generation;
    busy = true; controls(); status('Aligning read and collecting evidence…');
    try {
      await request('/api/sanger-verifications', {method:'POST', headers:{'Content-Type':'application/json'}, body:JSON.stringify({sequencing_read_id:button.dataset.analyze,previous_alignment_id:state.reads.get(button.dataset.analyze)?.alignment_id||null,direction:button.closest('.read-card').querySelector('[data-analysis-direction]').value,trim_quality_threshold:button.closest('.read-card').querySelector('[data-analysis-trim]').value?Number(button.closest('.read-card').querySelector('[data-analysis-trim]').value):null})});
      if (token !== generation) return;
      await refresh(token);
      if (token === generation) status('Analysis saved. Review coverage and differences below.');
    } catch (error) { if (token === generation) status(`Analysis could not finish: ${error.message}`, true); }
    finally { if (token === generation) { busy = false; controls(); } }
  });
  window.sangerView = {open, clear};
  clear();
})();
