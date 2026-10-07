/* Input integrity only. Revision generations guard every late response. */
(() => {
  const root=document.getElementById('input-jobs');
  const esc=v=>String(v??'').replace(/[&<>"']/g,c=>({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[c]));
  const date=v=>v?new Date(v).toLocaleString():'—';
  const labels={queued:'Queued · awaiting worker',running:'Checking inputs',succeeded:'Input check complete',failed:'Input check failed',cancelled:'Cancelled'};
  let state=null,generation=0,busy=false,historyVersion=0;
  const pending=new Map();
  async function request(url,options={}) {
    const response=await fetch(url,options);let body;
    try{body=await response.json();}catch{throw new Error('Unreadable server response');}
    if(!response.ok)throw new Error(typeof body.detail==='string'?body.detail:`Request failed (${response.status})`);
    return body;
  }
  const base=()=>`/api/revisions/${encodeURIComponent(state.id)}/input-checks`;
  function render() {
    if(!state){root.replaceChildren();return;}
    const disabled=busy?'disabled':'';
    root.innerHTML=`<header class="ij-header"><div><span class="eyebrow">Batch preparation</span><h3>Input integrity checks</h3><p>Confirm uploaded files and reference identity before analysis.</p></div><span class="ij-scope">No sequence analysis</span></header>
      <div class="ij-layout"><div class="ij-select"><h4>Select attached reads <span>${state.selected.size} / 100</span></h4><p>One check records the file hashes for this revision.</p><div class="ij-read-list">${state.reads.length?state.reads.map(r=>`<label><input type="checkbox" data-ij-read="${esc(r.id)}" ${state.selected.has(r.id)?'checked':''} ${disabled}><span>${esc(r.original_filename)}<small>${Number(r.length_bp).toLocaleString()} bases</small></span></label>`).join(''):'<p>Attach AB1 reads above, then refresh.</p>'}</div><button data-ij-submit ${busy||!state.selected.size||state.selected.size>100||!state.loaded?'disabled':''}>Queue input check</button></div>
      <div class="ij-tasks"><div class="ij-task-heading"><h4>Recent checks</h4><button data-ij-refresh ${disabled}>Refresh checks</button></div><p class="ij-note">Queued tasks wait for a local worker. Refresh to see changes.</p><div class="ij-status ${state.error?'ij-error':''}" role="status" aria-live="polite">${esc(state.message)}</div><div class="ij-job-list">${state.jobs.map(j=>`<article class="ij-job" data-ij-job="${esc(j.id)}"><div><span class="ij-badge ij-${esc(j.status)}">${esc(labels[j.status]||j.status)}</span><strong>${j.read_count} read${j.read_count===1?'':'s'} · ${j.attempt_count} / ${j.max_attempts} attempts</strong><small>${esc(date(j.created_at))} · <code>${esc(j.id)}</code></small>${j.error_detail?`<p class="ij-error">${esc(j.error_detail.replaceAll('_',' '))}</p>`:''}</div><div class="ij-actions"><button data-ij-detail="${esc(j.id)}" ${disabled}>Details &amp; attempts</button>${['queued','running'].includes(j.status)?`<button data-ij-cancel="${esc(j.id)}" ${disabled}>Cancel check</button>`:''}</div></article>`).join('')}${state.loaded&&!state.jobs.length?'<p class="ij-empty">No input checks submitted for this revision.</p>':''}</div><nav class="ij-pagination"><span>${state.total} check${state.total===1?'':'s'}</span><button data-ij-page="${Math.max(0,state.offset-5)}" ${busy||!state.offset?'disabled':''}>Newer checks</button><button data-ij-page="${state.offset+5}" ${busy||!state.more?'disabled':''}>Older checks</button></nav></div></div>
      <div class="ij-details">${state.detail?detailHTML(state.detail):''}</div><p class="ij-limit">A completed input check confirms file identity only. It does not validate the plasmid sequence, generate a consensus or replace Sanger review.</p>`;
  }
  function detailHTML(d) {
    const a=d.attempts;
    return `<header><h4>Check details</h4><button data-ij-close-detail>Close details</button></header><p><code>${esc(d.id)}</code> · ${esc(labels[d.status])}</p><p class="ij-note">Reference SHA-256 <code>${esc(d.inputs.reference.sequence_sha256)}</code></p>${d.result.scope?`<p class="ij-result">File identity checked · ${d.result.files?.length||0} files · No sequence analysis performed</p>`:''}<div class="ij-files">${d.inputs.inputs.reads.map(r=>`<div><b>${esc(d.read_labels?.[r.id]||'Read unavailable')}</b><small>Read <code>${esc(r.id)}</code></small><small>SHA-256 <code>${esc(r.sha256)}</code></small></div>`).join('')}</div><h4>Attempt history</h4>${a.items.length?`<ol>${a.items.map(r=>`<li><strong>Attempt ${r.number} · ${esc(r.status)}</strong><span>${esc(r.worker_id)} · ${esc(date(r.started_at))}</span>${r.error_detail?`<small class="ij-error">${esc(r.error_detail.replaceAll('_',' '))}</small>`:''}</li>`).join('')}</ol>`:'<p>No worker has claimed this check yet.</p>'}<nav class="ij-pagination"><span>${a.total} attempt${a.total===1?'':'s'} · newest first</span><button data-ij-history="${Math.max(0,a.offset-5)}" ${busy||!a.offset?'disabled':''}>Newer attempts</button><button data-ij-history="${a.offset+5}" ${busy||!a.has_more?'disabled':''}>Older attempts</button></nav>`;
  }
  function updateReads(id,reads) {
    if(state?.id!==id)return;
    state.reads=reads;state.selected=new Set([...state.selected].filter(i=>reads.some(r=>r.id===i)));render();
  }
  async function load(token,offset=state.offset) {
    const id=state.id,url=base();
    const [reads,list]=await Promise.all([request(`/api/sequence-revisions/${encodeURIComponent(id)}/sanger-reads`),request(`${url}?limit=5&offset=${offset}`)]);
    if(token!==generation)return;
    state.jobs=list.items;state.offset=list.offset;state.total=list.total;state.more=list.has_more;state.loaded=true;
    updateReads(id,reads);
  }
  async function action(work) {
    if(!state||busy)return;
    const token=generation;busy=true;state.error=false;state.message='Loading…';render();
    try{await work(token);}catch(error){if(token===generation){state.error=true;state.message=`Could not complete request: ${error.message}`;}}
    finally{if(token===generation){busy=false;render();}}
  }
  async function showDetail(job,offset=0) {
    const token=generation,version=++historyVersion;
    await action(async()=>{
      const detail=await request(`${base()}/${encodeURIComponent(job)}?limit=5&offset=${offset}`);
      if(token!==generation||version!==historyVersion)return;
      state.detail=detail;state.message='Attempt history loaded.';
    });
  }
  root.addEventListener('change',event=>{
    const id=event.target.dataset.ijRead;if(!id||busy||!state)return;
    if(event.target.checked)state.selected.add(id);else state.selected.delete(id);render();
  });
  root.addEventListener('click',event=>{
    const b=event.target.closest('button');if(!b||b.disabled||!state||busy)return;
    if(b.hasAttribute('data-ij-close-detail')){historyVersion++;state.detail=null;render();return;}
    if(b.dataset.ijDetail){showDetail(b.dataset.ijDetail);return;}
    if(b.hasAttribute('data-ij-history')){showDetail(state.detail.id,Number(b.dataset.ijHistory));return;}
    if(b.hasAttribute('data-ij-submit')){
      const ids=[...state.selected].sort(),signature=JSON.stringify([state.id,ids]);
      if(!pending.has(signature))pending.set(signature,crypto.randomUUID());
      const url=base(),key=pending.get(signature);
      action(async token=>{
        await request(url,{method:'POST',headers:{'Content-Type':'application/json'},body:JSON.stringify({read_ids:ids,idempotency_key:key})});
        pending.delete(signature);
        if(token!==generation)return;
        state.selected.clear();state.detail=null;state.message='Input check submitted. Waiting for a local worker.';
        await load(token,0);
      });return;
    }
    if(b.dataset.ijCancel){
      const url=`${base()}/${encodeURIComponent(b.dataset.ijCancel)}/cancel`;
      action(async token=>{await request(url,{method:'POST'});if(token!==generation)return;state.detail=null;state.message='Check cancelled. Late results cannot be published.';await load(token);});return;
    }
    const offset=b.hasAttribute('data-ij-page')?Number(b.dataset.ijPage):state.offset;
    action(async token=>{await load(token,offset);if(token!==generation)return;state.detail=null;state.message='Status updated.';});
  });
  function clear(){generation++;historyVersion++;state=null;busy=false;render();}
  function open(id){clear();state={id,reads:[],selected:new Set(),jobs:[],offset:0,total:0,more:false,loaded:false,detail:null,message:'',error:false};action(async token=>{await load(token,0);if(token===generation)state.message='Select reads to start an input check.';});}
  window.inputJobs={open,clear,updateReads};
})();
