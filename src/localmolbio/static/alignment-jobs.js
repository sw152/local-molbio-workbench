/* Local evidence only; revision generations fence all asynchronous updates. */
(() => {
  const root=document.getElementById('alignment-jobs');
  const esc=v=>String(v??'').replace(/[&<>"']/g,c=>({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[c]));
  const n=v=>Number(v).toLocaleString();
  const types={'ont-noisy':'ONT · noisy reads','ont-high-accuracy':'ONT · high accuracy','pacbio-hifi':'PacBio HiFi','short-single':'Short reads · unpaired'};
  const labels={queued:'Queued · awaiting worker',running:'Aligning reads',succeeded:'Local evidence ready',failed:'Alignment failed',cancelled:'Cancelled'};
  let state=null,generation=0,busy=false;
  const pending=new Map();
  const base=()=>`/api/revisions/${encodeURIComponent(state.id)}`;
  async function request(url,options={}){
    const response=await fetch(url,options);let data;
    try{data=await response.json();}catch{throw new Error('Unreadable server response. Refresh before retrying.');}
    if(!response.ok)throw new Error(typeof data.detail==='string'?data.detail.replaceAll('_',' '):`Request rejected (${response.status}).`);
    return data;
  }
  function pager(kind,page,limit,previous,next){return `<nav class="ij-pagination" aria-label="${esc(kind)} pagination"><span>${page.total?n(page.offset+1)+'–'+n(page.offset+page.items.length)+' of '+n(page.total):'0 items'}</span><button data-aj-${kind}="${Math.max(0,page.offset-limit)}" ${busy||!page.offset?'disabled':''}>${previous}</button><button data-aj-${kind}="${page.offset+limit}" ${busy||!page.has_more?'disabled':''}>${next}</button></nav>`;}
  function diagram(hit,length){
    const x=v=>40+620*v/length;
    const blocks=hit.paired_blocks.flatMap(b=>b.reference_intervals);
    const bands=blocks.map(([start,end])=>`<rect x="${x(start)}" y="39" width="${Math.max(.4,x(end)-x(start))}" height="14" rx="1"/>`).join('');
    const arrows=hit.reference_intervals.map(([start,end])=>{const endX=x(hit.strand==='-'?start:end),dir=hit.strand==='-'?-1:1;return `<path d="M ${endX-dir*6} 34 L ${endX} 46 L ${endX-dir*6} 58" fill="none" stroke="currentColor" stroke-width="2"/>`;}).join('');
    return `<p class="aj-scroll-hint">Reference track · scroll sideways to inspect all positions →</p><div class="aj-diagram" tabindex="0" aria-label="Reference alignment diagram; scroll horizontally on narrow screens"><svg viewBox="0 0 700 92" role="img" aria-label="${esc(hit.strand==='-'?'Reverse':'Forward')} alignment on a ${length} base reference; colored blocks are paired bases"><line x1="40" x2="660" y1="46" y2="46" stroke="#405570"/>${[0,Math.floor(length/2),length].map((v,i)=>`<line x1="${x(v)}" x2="${x(v)}" y1="60" y2="66" stroke="#748ba6"/><text x="${x(v)}" y="81" text-anchor="${i===0?'start':i===2?'end':'middle'}">${n(v)} bp</text>`).join('')}<g fill="currentColor">${bands}</g>${arrows}</svg></div>`;
  }
  function hitHTML(h,index,length,matched=null){
    const ranges=h.reference_intervals.map(([s,e])=>`[${n(s)}, ${n(e)})`).join(' → ');
    return `<details class="aj-hit" ${matched===true||(matched===null&&index===0)?'open':''}><summary>Hit ${index+1}${matched===true?' · selected region':matched===false?' · outside selected relation':''} · ${h.strand==='-'?'← Reverse':'Forward →'}${h.crosses_origin?' · crosses origin':''} · ${n(h.exact_matches)} exact bases</summary>${diagram(h,length)}<p class="aj-coordinates"><b>Reference</b> ${ranges}<br><b>Original read</b> [${n(h.query_start)}, ${n(h.query_end)}) · ${h.strand==='-'?'read coordinates decrease along reference':'read coordinates increase along reference'}</p><div class="aj-stats"><span><b>${(h.local_exact_identity*100).toFixed(1)}%</b>Local exact identity</span><span><b>${(h.aligned_query_fraction*100).toFixed(1)}%</b>Read aligned</span><span><b>${n(h.mismatches)}</b>Mismatches</span><span><b>${n(h.inserted_bases)} / ${n(h.deleted_bases)}</b>Inserted / deleted</span><span><b>${n(h.ambiguous_pairs)}</b>Ambiguous pairs</span></div><p class="ij-note">Identity includes all ambiguous and gap columns. Colored blocks show paired bases; deletions remain gaps. Arrows indicate read direction. Coordinates are 0-based, half-open.</p><details><summary>CIGAR &amp; reported hits</summary><code class="aj-cigar">${esc(h.cigar)}</code><p>${h.raw_hits.map(r=>`${r.type==='P'?'Primary':'Secondary'} · MAPQ ${r.mapq===null?'unavailable':n(r.mapq)}`).join('<br>')}</p><p class="ij-note">Raw MAPQ is not recalibrated for circular references. Reported candidates are not exhaustive.</p></details></details>`;
  }
  function readHTML(r){
    const d=state.detail,label=d.input_labels[r.source.input_id]||r.source.input_id;
    return `<article class="aj-read"><header><div><span class="eyebrow">Record ${r.source.record_ordinal} · ${esc(r.source.query_name)}</span><h4>${esc(label)}</h4><small>${n(r.length_bp)} bases · ${r.alignments.length} reported alignment${r.alignments.length===1?'':'s'}</small></div><span class="ij-badge aj-review">${r.withheld.length?'Evidence withheld · review':r.alignments.length===0?'No alignment reported':r.alignments.length>1?'Multiple alignments · review':'Local match · uniqueness not assessed'}</span></header>${r.withheld.length?`<p class="ij-error">Withheld evidence: ${r.withheld.map(w=>esc(w.reason.replaceAll('_',' '))).join('; ')}</p>`:''}${!r.alignments.length?'<p class="ij-note">No usable local alignment is reported. Review the input and data type; this is not a whole-plasmid verdict.</p>':''}${r.region_evidence?`<div class="aj-region-evidence"><b>Evidence inside selected region</b><p>Paired positions: ${intervalText(r.region_evidence.paired_intervals)}<br>Deletion positions: ${intervalText(r.region_evidence.deletion_intervals)}</p><small>Union across reported candidates; both types can overlap. Matching hits: ${r.region_evidence.matching_hit_indices.map(i=>i+1).join(', ')}.</small></div>`:''}${r.alignments.map((h,i)=>hitHTML(h,i,d.result.reference_length_bp,r.region_evidence?r.region_evidence.matching_hit_indices.includes(i):null)).join('')}<details class="aj-source"><summary>Read provenance</summary><p>Input ID <code>${esc(r.source.input_id)}</code><br>File record (1-based) ${r.source.record_ordinal}<br>Uppercase sequence SHA-256 <code>${esc(r.source.sequence_sha256)}</code></p></details></article>`;
  }
  function coverageHTML(c){
    if(!c)return '';
    const s=c.summary,r=c.regions;
    const bands=[['Single reported hit evidence',s.single_reported_alignment_bases,'single'],['Alternative-only evidence',s.ambiguous_only_bases,'ambiguous'],['No paired evidence',s.unpaired_bases,'unpaired']];
    const kinds={unpaired:'No paired evidence',ambiguous_only:'Alternative-only evidence',deletion:'Reported deletion positions'};
    return `<section class="aj-coverage" aria-label="Paired-position overview"><h4>Paired-position overview</h4><p class="ij-note">${n(s.reference_length_bp)} reference bases · each record counted once at each position, across its candidate hits.</p><p class="ij-note">Reference base composition</p><div class="aj-coverage-bar" aria-hidden="true">${bands.map(([label,count,kind])=>`<span class="aj-cov-${kind}" style="width:${count/s.reference_length_bp*100}%"></span>`).join('')}</div><div class="aj-coverage-legend">${bands.map(([label,count,kind])=>`<div><i class="aj-cov-${kind}"></i><span>${label}</span><b>${n(count)} bp</b></div>`).join('')}</div><p class="ij-note">Paired positions include mismatches and N. Quality, uniqueness and independent molecules are not assessed. Single-hit evidence takes precedence where both categories overlap.</p><p class="aj-coverage-counts">${n(s.deletion_evidence_bases)} positions with deletion evidence · ${n(s.read_counts.withheld)} record${s.read_counts.withheld===1?'':'s'} with withheld evidence · ${n(s.read_counts.no_alignment_reported)} record${s.read_counts.no_alignment_reported===1?'':'s'} without reported alignments</p><label for="aj-region-kind">Inspect reference regions</label><select id="aj-region-kind" ${busy?'disabled':''}>${Object.entries(kinds).map(([k,v])=>`<option value="${k}" ${r.kind===k?'selected':''}>${v}</option>`).join('')}</select><p class="ij-note">${r.kind==='deletion'?'Deletion evidence can overlap paired positions from other records. These are reported CIGAR intervals; repeat-equivalent deletions are not normalized or called as consensus.':r.kind==='ambiguous_only'?'Only alternative or partly withheld placement evidence is reported here.':'No reported read base pairs with these positions. This includes uncovered deletion positions.'} Intervals are 0-based, half-open; origin ends are listed separately.</p><div class="aj-regions">${r.items.map(i=>`<div><code>[${n(i.start)}, ${n(i.end)})</code><b>${n(i.length_bp)} bp</b><button data-aj-inspect="${i.start}" data-aj-end="${i.end}" ${busy?'disabled':''}>Inspect related reads</button></div>`).join('')||'<p>No regions in this category. This is not a plasmid pass verdict.</p>'}</div>${pager('regions',r,5,'Previous regions','Next regions')}</section>`;
  }
  function intervalText(intervals){return intervals.map(([s,e])=>`[${n(s)}, ${n(e)})`).join(' · ')||'none reported';}
  function regionFilterHTML(){
    const page=state.reads,r=page?.region;
    if(!r)return '';
    return `<section class="aj-region-filter" tabindex="-1" aria-label="Read region filter"><div><span class="eyebrow">Focused read review</span><h4>Reference [${n(r.start)}, ${n(r.end)})</h4><p>${n(page.total)} matching record${page.total===1?'':'s'} / ${n(page.unfiltered_total)} total · 0-based, half-open</p></div><label for="aj-relation">Evidence overlapping this region</label><select id="aj-relation" ${busy?'disabled':''}>${Object.entries({either:'Paired or deletion evidence',paired:'Paired positions only',deletion:'Deletion positions only'}).map(([key,label])=>`<option value="${key}" ${r.relation===key?'selected':''}>${label}</option>`).join('')}</select><p class="ij-note">Every reported candidate is checked; each record appears once. Withheld-only and unreported placements cannot be located here.</p><button data-aj-all-reads ${busy?'disabled':''}>Show all reads</button><button data-aj-overview>Back to overview</button>${!page.total?'<p class="ij-note">No usable reported alignment overlaps with this relation. This does not prove absence of biological evidence.</p>':''}</section>`;
  }
  function detailHTML(d){
    const result=d.result;
    return `<header><h4>Alignment review</h4><button data-aj-close ${busy?'disabled':''}>Close alignment review</button></header><p><code>${esc(d.id)}</code> · ${esc(labels[d.status])}</p><p class="ij-note">${esc(types[d.data_type])} · ${esc(d.tool.name)} ${esc(d.tool.version)}</p><details><summary>Reference &amp; input identities</summary><p>Reference SHA-256 <code>${esc(d.reference.sequence_sha256)}</code></p>${d.input_identities.map(i=>`<p><b>${esc(d.input_labels[i.id]||i.id)}</b><br>Original SHA-256 <code>${esc(i.sha256)}</code></p>`).join('')}<p>Tool SHA-256 <code>${esc(d.tool.binary_sha256)}</code></p></details><h4 class="aj-subhead">Attempt history</h4>${d.attempts.items.length?`<ol>${d.attempts.items.map(a=>`<li><b>Attempt ${a.number} · ${esc(a.status)}</b><span>${esc(a.started_at)}</span>${a.error_detail?`<small class="ij-error">${esc(a.error_detail.replaceAll('_',' '))}</small>`:''}</li>`).join('')}</ol>`:'<p class="ij-note">No worker has claimed this task.</p>'}${pager('attempts',d.attempts,3,'Newer attempts','Older attempts')}${result?`<div class="aj-boundary"><b>Local alignments · ${n(result.read_count)} reads</b><p>Base quality and pairing were not used. No consensus or whole-plasmid verification. A single reported match does not establish uniqueness.</p></div>${coverageHTML(state.coverage)}${regionFilterHTML()}<div class="aj-read-list">${state.reads?.items.map(readHTML).join('')||''}</div>${state.reads?pager('reads',state.reads,3,'Previous reads','Next reads'):''}`:'<p class="ij-note">Read evidence appears only after the worker publishes a successful result. Refresh alignment tasks to update.</p>'}`;
  }
  function render(){
    if(!state){root.replaceChildren();return;}
    root.setAttribute('aria-busy',String(busy));const disabled=busy?'disabled':'';
    root.innerHTML=`<header class="ij-header"><div><span class="eyebrow">Sequence evidence</span><h3>FASTQ alignment</h3><p>Trace every local match to its original file and record.</p></div><span class="ij-scope">Local evidence · no plasmid verdict</span></header><div class="ij-layout"><div class="ij-select"><h4>Choose registered files <span>${state.selected.size} / 100</span></h4><div class="ij-read-list">${state.inputs.items.map(i=>`<label><input type="checkbox" data-aj-input="${esc(i.id)}" ${state.selected.has(i.id)?'checked':''} ${disabled}><span>${esc(i.original_filename)}<small>${n(i.summary.records)} reads · ${n(i.size_bytes)} bytes</small></span></label>`).join('')||'<p>Register FASTQ files above, then refresh alignment tasks.</p>'}</div>${pager('inputs',state.inputs,5,'Previous files','Next files')}<label class="aj-type-label" for="aj-type">Sequencing data type</label><select id="aj-type" ${disabled}><option value="">Choose explicitly…</option>${Object.entries(types).map(([k,v])=>`<option value="${k}" ${state.type===k?'selected':''}>${v}</option>`).join('')}</select><p>${state.selected.size} files selected across pages.</p><button data-aj-submit ${busy||!state.loaded||!state.type||!state.selected.size||state.selected.size>100?'disabled':''}>Queue alignment</button><details class="aj-limits"><summary>Input &amp; analysis limits</summary><p>Up to 32 MiB original and decoded total; 5,000 reads, 100,000 bases per read, 5 million total bases, 200,000-base reference. Oversized inputs are rejected in full. ACGTN only; other IUPAC bases are rejected. Quality is retained but not used for alignment. Files remain unpaired.</p></details></div><div class="ij-tasks"><div class="ij-task-heading"><h4>Alignment tasks</h4><button data-aj-refresh ${disabled}>Refresh alignment tasks</button></div><p class="ij-note">Submission queues a task. A separately invoked local worker executes it. Refresh to see progress.</p><p class="ij-status ${state.error?'ij-error':''}" role="status" aria-live="polite">${esc(state.message)}</p>${state.jobs.items.map(j=>`<article class="ij-job"><div><span class="ij-badge ij-${esc(j.status)}">${esc(labels[j.status])}</span><strong>${esc(types[j.data_type])} · ${j.input_count} files</strong><small>${j.attempt_count} / ${j.max_attempts} attempts · <code>${esc(j.id)}</code></small>${j.error_detail?`<p class="ij-error">${esc(j.error_detail.replaceAll('_',' '))}</p>`:''}</div><div class="ij-actions"><button data-aj-detail="${esc(j.id)}" ${disabled}>Review alignment</button>${['queued','running'].includes(j.status)?`<button data-aj-cancel="${esc(j.id)}" ${disabled}>Cancel alignment</button>`:''}</div></article>`).join('')||'<p class="ij-empty">No alignment tasks for this revision.</p>'}${pager('jobs',state.jobs,3,'Newer tasks','Older tasks')}</div></div><div class="ij-details aj-details">${state.detail?detailHTML(state.detail):''}</div>`;
  }
  async function action(work){
    if(!state||busy)return;const token=generation;busy=true;state.error=false;state.message='Loading…';render();
    try{await work(token);}catch(error){if(token===generation){state.error=true;state.message=`Request not completed: ${error.message} Displayed evidence may be out of date. Refresh to recover; retry uncertain submissions with the same selection and type.`;}}
    finally{if(token===generation){busy=false;render();}}
  }
  async function load(token,inputOffset=state.inputs.offset,jobOffset=state.jobs.offset){
    const url=base();
    const [inputs,jobs]=await Promise.all([request(`${url}/fastq-inputs?limit=5&offset=${inputOffset}`),request(`${url}/alignments?limit=3&offset=${jobOffset}`)]);
    if(token!==generation)return;state.inputs=inputs;state.jobs=jobs;state.loaded=true;
  }
  async function detail(token,id,attemptOffset=0,readOffset=0,kind=state.coverage?.regions.kind||'unpaired',regionOffset=state.coverage?.regions.offset||0,region=state.detail?.id===id?state.reads?.region:null){
    const url=`${base()}/alignments/${encodeURIComponent(id)}`;
    const d=await request(`${url}?limit=3&offset=${attemptOffset}`);
    if(token!==generation)return;
    const filter=region?`&start=${region.start}&end=${region.end}&relation=${region.relation}`:'';
    const [reads,coverage]=d.result?await Promise.all([request(`${url}/reads?limit=3&offset=${readOffset}${filter}`),request(`${url}/coverage?kind=${kind}&limit=5&offset=${regionOffset}`)]):[null,null];
    if(token!==generation)return;state.detail=d;state.reads=reads;state.coverage=coverage;state.message='Alignment evidence loaded.';
  }
  root.addEventListener('change',event=>{
    if(!state||busy)return;
    if(event.target.id==='aj-region-kind'){const kind=event.target.value;action(t=>detail(t,state.detail.id,state.detail.attempts.offset,state.reads.offset,kind,0));return;}
    if(event.target.id==='aj-relation'){const region={...state.reads.region,relation:event.target.value};action(t=>detail(t,state.detail.id,state.detail.attempts.offset,0,state.coverage.regions.kind,state.coverage.regions.offset,region));return;}
    if(event.target.id==='aj-type')state.type=event.target.value;
    const id=event.target.dataset.ajInput;
    if(id){if(event.target.checked)state.selected.add(id);else state.selected.delete(id);}
    render();
  });
  root.addEventListener('click',event=>{
    const b=event.target.closest('button');if(!b||b.disabled||busy||!state)return;
    if(b.hasAttribute('data-aj-close')){state.detail=null;state.reads=null;state.coverage=null;render();return;}
    if(b.dataset.ajDetail){action(t=>detail(t,b.dataset.ajDetail,0,0,'unpaired',0,null));return;}
    if(b.hasAttribute('data-aj-overview')){root.querySelector('.aj-coverage')?.scrollIntoView({block:'start'});return;}
    if(b.hasAttribute('data-aj-inspect')||b.hasAttribute('data-aj-all-reads')){
      const region=b.hasAttribute('data-aj-inspect')?{start:Number(b.dataset.ajInspect),end:Number(b.dataset.ajEnd),relation:'either'}:null;
      const token=generation;action(t=>detail(t,state.detail.id,state.detail.attempts.offset,0,state.coverage.regions.kind,state.coverage.regions.offset,region)).then(()=>{if(token===generation&&!state.error)root.querySelector('.aj-region-filter, .aj-read-list')?.scrollIntoView({block:'start'});});return;
    }
    if(b.hasAttribute('data-aj-attempts')){action(t=>detail(t,state.detail.id,Number(b.dataset.ajAttempts),state.reads?.offset||0));return;}
    if(b.hasAttribute('data-aj-reads')){action(t=>detail(t,state.detail.id,state.detail.attempts.offset,Number(b.dataset.ajReads)));return;}
    if(b.hasAttribute('data-aj-regions')){action(t=>detail(t,state.detail.id,state.detail.attempts.offset,state.reads.offset,state.coverage.regions.kind,Number(b.dataset.ajRegions)));return;}
    if(b.hasAttribute('data-aj-submit')){
      const ids=[...state.selected].sort(),type=state.type,signature=JSON.stringify([state.id,ids,type]);
      if(!pending.has(signature))pending.set(signature,crypto.randomUUID());
      const url=base(),body=JSON.stringify({input_ids:ids,data_type:type,idempotency_key:pending.get(signature)});
      action(async token=>{await request(`${url}/alignments`,{method:'POST',headers:{'Content-Type':'application/json'},body});pending.delete(signature);if(token!==generation)return;state.selected.clear();state.detail=null;state.reads=null;state.coverage=null;await load(token,0,0);if(token===generation)state.message='Alignment queued. Awaiting a local worker.';});return;
    }
    if(b.dataset.ajCancel){
      const url=`${base()}/alignments/${encodeURIComponent(b.dataset.ajCancel)}/cancel`;
      action(async token=>{await request(url,{method:'POST'});if(token!==generation)return;state.detail=null;state.reads=null;state.coverage=null;await load(token);if(token===generation)state.message='Alignment cancelled. Late results cannot be published.';});return;
    }
    const inputOffset=b.hasAttribute('data-aj-inputs')?Number(b.dataset.ajInputs):0;
    const jobOffset=b.hasAttribute('data-aj-jobs')?Number(b.dataset.ajJobs):0;
    const selectedDetail=state.detail?.id;
    action(async token=>{await load(token,inputOffset,jobOffset);if(token!==generation)return;if(selectedDetail)await detail(token,selectedDetail);if(token===generation)state.message='Alignment status refreshed.';});
  });
  function clear(){generation++;state=null;busy=false;root.removeAttribute('aria-busy');render();}
  function open(id){clear();const empty=()=>({items:[],total:0,offset:0,has_more:false});state={id,type:'',selected:new Set(),inputs:empty(),jobs:empty(),loaded:false,message:'',error:false,detail:null,reads:null,coverage:null};action(async token=>{await load(token);if(token===generation)state.message='Choose files and their sequencing data type.';});}
  window.alignmentJobs={open,clear};
})();
