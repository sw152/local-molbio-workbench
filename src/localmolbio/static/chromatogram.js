/* All coordinates here refer to the original uploaded read, including reverse alignments. */
(() => {
  const colors = {A:'#65df98', C:'#6ebaff', G:'#ffd16c', T:'#f58eb9', N:'#aeb9cd'};
  const escape = value => String(value).replace(/[&<>"']/g, c => ({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[c]));
  const requests = new WeakMap();
  function plot(trace) {
    const width = Math.max(600, trace.bases.length * 32), pad = 16;
    const samples = Object.values(trace.channels).flat();
    const min = samples.reduce((a,b)=>Math.min(a,b),0), max = samples.reduce((a,b)=>Math.max(a,b),1);
    const x = sample => pad + (sample-trace.sample_start) / Math.max(1, trace.sample_end-trace.sample_start-1) * (width-pad*2);
    const y = amplitude => 196 - (amplitude-min)/(max-min)*140;
    const baseline = y(0);
    const shades = trace.bases.filter(b => b.phred < 20).map(b => `<rect x="${x(b.sample)-11}" y="8" width="22" height="244" fill="#ffc36e15"/>`).join('');
    const paths = Object.entries(trace.channels).map(([base, values]) => `<path data-channel="${base}" d="${values.map((v,i) => `${i?'L':'M'}${x(trace.sample_start+i).toFixed(2)},${y(v).toFixed(2)}`).join(' ')}" fill="none" stroke="${colors[base]}" stroke-width="1.7" stroke-linejoin="round"/>`).join('');
    const calls = trace.bases.map(b => `<g data-peak-position="${b.position}"><rect class="trace-focus" x="${x(b.sample)-12}" y="4" width="24" height="272" rx="3" fill="none" stroke="#55e5d7"/><title>Original read base ${b.position+1}: ${escape(b.base)}, Phred ${b.phred}</title><text x="${x(b.sample)}" y="19" class="trace-number">${b.position+1}</text><text x="${x(b.sample)}" y="39" fill="${colors[b.base] || colors.N}" class="trace-call">${escape(b.base)}</text><rect x="${x(b.sample)-5}" y="${250-Math.min(60,b.phred)/60*34}" width="10" height="${Math.min(60,b.phred)/60*34}" rx="2" fill="${b.phred<20?'#ffc36e':'#7892b8'}"/><text x="${x(b.sample)}" y="271" class="trace-number">${b.phred}</text></g>`).join('');
    return `<svg class="trace-svg" style="width:${width}px" viewBox="0 0 ${width} 282" role="img" aria-label="Chromatogram, original read bases ${trace.base_start+1} through ${trace.base_end}; quality values below"><rect width="${width}" height="282" fill="#09101d"/>${shades}<path d="M 0 ${baseline} H ${width} M 0 207 H ${width}" stroke="#35445d" stroke-width="1"/>${paths}${calls}</svg>`;
  }
  function render(trace) {
    return `<header class="trace-heading"><strong>Chromatogram · original read</strong><div class="trace-legend">${Object.entries(colors).filter(([b])=>b!=='N').map(([base,color]) => `<span style="color:${color}">● ${base}</span>`).join('')}</div></header><p class="evidence-note">Bases ${trace.base_start+1}–${trace.base_end} of ${trace.read_length}. Original upload orientation, including for reverse alignments. Shared signal scale; shaded calls are below Q20.</p><div class="trace-scroll" tabindex="0" aria-label="Chromatogram, scroll horizontally for all bases">${plot(trace)}</div><p class="trace-caption">Top: original read position and base call · Bottom: Phred quality (bars capped at Q60). Scroll horizontally to inspect the window.</p><div class="trace-navigation"><button data-trace-start="${Math.max(0,trace.base_start-24)}" ${trace.base_start===0?'disabled':''}>← Previous</button><label>Go to read base <input class="trace-position" type="number" min="1" max="${trace.read_length}" value="${trace.base_start+1}" aria-label="Go to read base"></label><button data-trace-go>Go</button><button data-trace-start="${trace.base_end}" ${trace.base_end>=trace.read_length?'disabled':''}>Next →</button></div>`;
  }
  async function load(card, start, focus = null) {
    const slot = card.querySelector('.trace-view');
    const token = (requests.get(slot) || 0) + 1;
    requests.set(slot, token);
    slot.hidden = false;
    slot.setAttribute('aria-busy','true');
    slot.innerHTML = '<p class="evidence-note" role="status">Loading original chromatogram…</p>';
    try {
      const response = await fetch(`/api/sanger-reads/${encodeURIComponent(card.dataset.readId)}/trace?start=${start}&count=24`);
      const trace = await response.json();
      if (!response.ok) throw new Error(typeof trace.detail === 'string' ? trace.detail : `Request failed (${response.status})`);
      if (!slot.isConnected || requests.get(slot) !== token) return;
      slot.innerHTML = trace.available ? render(trace) : `<p class="evidence-note" role="status">Chromatogram unavailable: ${escape(trace.reason)}</p>`;
      if (trace.available && focus !== null) {
        const peak = slot.querySelector(`[data-peak-position="${focus}"]`);
        if (peak) {
          peak.classList.add('trace-selected');
          const scroll = slot.querySelector('.trace-scroll');
          scroll.scrollLeft += peak.getBoundingClientRect().left - scroll.getBoundingClientRect().left - scroll.clientWidth/2 + 12;
        }
      }
    } catch (error) {
      if (slot.isConnected && requests.get(slot) === token) slot.innerHTML = `<p class="trace-error" role="alert">Could not load chromatogram: ${escape(error.message)}. Use “View chromatogram” to retry.</p>`;
    } finally {
      if (slot.isConnected && requests.get(slot) === token) slot.removeAttribute('aria-busy');
    }
  }
  document.getElementById('sanger-summary').addEventListener('click', async event => {
    const button=event.target.closest('[data-quality-read]');
    if(!button)return;
    const card=Array.from(document.querySelectorAll('#sanger-reads .read-card')).find(c=>c.dataset.readId===button.dataset.qualityRead);
    const position=Number(button.dataset.qualityPosition);
    if(!card||!Number.isInteger(position)||position<0)return;
    await load(card,Math.max(0,position-12),position);
    if(card.isConnected)card.querySelector('.trace-view').scrollIntoView({block:'start'});
  });
  document.getElementById('sanger-reads').addEventListener('click', event => {
    const button = event.target.closest('[data-trace-open], [data-trace-position], [data-trace-start], [data-trace-go]');
    if (!button || button.disabled) return;
    const card = button.closest('.read-card');
    let start = 0, focus = null;
    if (button.hasAttribute('data-trace-position')) { focus = Number(button.dataset.tracePosition); start = Math.max(0, focus-12); }
    else if (button.hasAttribute('data-trace-start')) start = Number(button.dataset.traceStart);
    else if (button.hasAttribute('data-trace-go')) {
      const input = card.querySelector('.trace-position');
      if (!input.reportValidity() || !input.value) return;
      start = Number(input.value)-1;
    }
    load(card, start, focus);
  });
})();
