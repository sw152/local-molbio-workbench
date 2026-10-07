/* Revision-scoped read evidence. Pending responses never update a different construct. */
(() => {
  const el = id => document.getElementById(id);
  const escape = value => String(value ?? '').replace(/[&<>"']/g, c => ({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[c]));
  const percent = value => Number.isFinite(value) ? `${(value * 100).toFixed(1)}%` : 'Unavailable';
  const flags = {
    ambiguous_alignment: 'Multiple equally scoring alignments — placement needs review',
    candidate_search_truncated: 'Repeated sequence: search limit reached',
    partial_read_alignment: 'Less than 80% of this read aligned',
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
    el('sanger-reads').querySelectorAll('[data-analyze]').forEach(b => b.disabled = busy);
  }
  function clear() {
    generation++;
    state = null;
    busy = false;
    el('sanger-form').reset();
    el('sanger-reads').replaceChildren();
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
  function render(reads) {
    el('sanger-reads').innerHTML = reads.length ? reads.map(read => {
      const evidence = read.evidence || {}, analyzed = Boolean(read.alignment_id);
      const review = evidence.review_flags || [];
      return `<article class="read-card" data-read-id="${escape(read.id)}"><header class="read-heading"><div><b>${escape(read.original_filename)}</b><small>${Number(read.length_bp).toLocaleString()} bases · ${escape(analyzed ? evidence.direction || read.direction : read.direction)} orientation</small></div><span class="badge ${review.length ? 'review' : ''}">${analyzed ? review.length ? 'Review needed' : 'Aligned region' : 'Ready to analyze'}</span></header><div class="read-metrics"><div><span>Mean quality</span><strong>${read.quality_summary?.mean_phred == null ? 'Unavailable' : `Q${Number(read.quality_summary.mean_phred).toFixed(1)}`}</strong></div><div><span>Aligned identity</span><strong>${analyzed ? percent(read.identity_fraction) : '—'}</strong></div><div><span>Read aligned</span><strong>${analyzed ? percent(evidence.read_aligned_fraction) : '—'}</strong></div></div>${analyzed ? coverage(evidence, state.length) + (review.length ? `<ul class="read-flags">${review.map(flag => `<li>${escape(flags[flag] || flag)}</li>`).join('')}</ul>` : '') + differences(read.variants || []) : `<button class="export" data-analyze="${escape(read.id)}">Analyze against this revision</button>`}<button class="trace-toggle" data-trace-open>View chromatogram</button><div class="trace-view" hidden></div></article>`;
    }).join('') : '<div class="read-empty"><span>＋</span><b>No reads attached to this revision</b><p>Add an AB1 chromatogram to review its sequence evidence.</p></div>';
    controls();
  }
  async function refresh(token) {
    const reads = await request(`/api/sequence-revisions/${encodeURIComponent(state.id)}/sanger-reads`);
    if (token === generation) render(reads);
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
      await request('/api/sanger-verifications', {method:'POST', headers:{'Content-Type':'application/json'}, body:JSON.stringify({sequencing_read_id:button.dataset.analyze})});
      if (token !== generation) return;
      await refresh(token);
      if (token === generation) status('Analysis saved. Review coverage and differences below.');
    } catch (error) { if (token === generation) status(`Analysis could not finish: ${error.message}`, true); }
    finally { if (token === generation) { busy = false; controls(); } }
  });
  window.sangerView = {open, clear};
  clear();
})();
