(function(root){
  const esc=v=>String(v??'').replace(/[&<>"']/g,c=>({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[c]));
  function pairsFrom(candidates){
    const groups=new Map();
    for(const q of candidates){
      const m=q.metrics||{};
      if(!m.analysis_job_id||!Number.isInteger(m.pair_index))continue;
      const key=JSON.stringify([m.analysis_job_id,m.pair_index]);
      if(!groups.has(key))groups.set(key,[]);
      groups.get(key).push(q);
    }
    return [...groups].flatMap(([key,items])=>{
      const left=items.find(q=>q.direction==='forward'),right=items.find(q=>q.direction==='reverse');
      if(items.length!==2||!left||!right)return [];
      const coords=[left.binding_start,left.binding_end,right.binding_start,right.binding_end];
      if(!coords.every(Number.isInteger)||!(0<=coords[0]&&coords[0]<coords[1]&&coords[1]<=coords[2]&&coords[2]<coords[3]))return [];
      const size=coords[3]-coords[0];
      if(left.metrics.product_size!==size||right.metrics.product_size!==size)return [];
      const target=left.metrics.target;
      if(JSON.stringify(target)!==JSON.stringify(right.metrics.target))return [];
      if(target&&!(Number.isInteger(target.start)&&Number.isInteger(target.end)&&coords[1]<=target.start&&target.start<target.end&&target.end<=coords[2]))return [];
      return [{key,left,right,target,size}];
    });
  }
  function diagram(pair){
    const {left,right,target,size}=pair,start=left.binding_start,end=right.binding_end;
    const x=v=>60+640*(v-start)/size;
    const arrow=(q,reverse)=>{const a=x(q.binding_start),b=x(q.binding_end),head=Math.min(9,(b-a)/2);return reverse?`${b},103 ${a+head},103 ${a},112 ${a+head},121 ${b},121`:`${a},103 ${b-head},103 ${b},112 ${b-head},121 ${a},121`;};
    const ticks=[...new Set([start,...[1,2,3].map(i=>Math.round(start+size*i/4)),end])];
    return `<svg viewBox="0 0 760 210" role="img" aria-label="Primer pair ${esc(left.metrics.pair_index)}; forward arrow points right, reverse arrow points left; ${size} bp product"><path d="M60 40 V30 H700 V40" fill="none" stroke="#a084ff"/><text x="380" y="20" text-anchor="middle" fill="#d3c6ff">${size} bp expected product</text><line x1="60" y1="112" x2="700" y2="112" stroke="#7386a9" stroke-width="3"/>${target?`<rect x="${x(target.start)}" y="89" width="${x(target.end)-x(target.start)}" height="46" rx="3" fill="#ffc36e33" stroke="#ffc36e"/><text x="380" y="69" text-anchor="middle" fill="#ffcf8a">Target ${target.start+1}–${target.end} bp</text>`:'<text x="380" y="69" text-anchor="middle" fill="#95a3bd">Unconstrained placement</text>'}<polygon points="${arrow(left,false)}" fill="#55e5d7"/><polygon points="${arrow(right,true)}" fill="#a084ff"/><text x="60" y="92" fill="#55e5d7">F →</text><text x="700" y="92" text-anchor="end" fill="#c3b0ff">← R</text>${ticks.map(t=>`<line x1="${x(t)}" x2="${x(t)}" y1="148" y2="155" stroke="#7386a9"/><text x="${x(t)}" y="174" text-anchor="middle" fill="#bdcce3">${t}</text>`).join('')}<text x="380" y="202" text-anchor="middle" fill="#95a3bd">Reference boundaries · 0-based · product window to scale</text></svg>`;
  }
  function siteCard(q){
    const r=q.metrics.reference_sites;
    const heading=`<h4>${q.direction==='forward'?'Forward →':'← Reverse'} · ${esc(q.name)}</h4><code>${esc(q.sequence_text)}</code><p>Intended site ${q.binding_start+1}–${q.binding_end} bp · ${esc(q.selection_state)}</p>`;
    if(!r)return `<article>${heading}<p>Exact-site review not recorded for this candidate.</p></article>`;
    const n=r.total_directional_matches;
    return `<article>${heading}<strong class="${n>1?'site-warning':''}">${n} exact directional match${n===1?'':'es'}</strong><p>${r.reference_topology==='circular'?'Circular origin included':r.reference_topology==='unknown'?'Unknown topology: origin not checked':'Linear reference'}</p><details><summary>Review ${r.sites.length} stored site${r.sites.length===1?'':'s'}${r.truncated?` of ${n} (truncated)`:''}</summary><ul>${r.sites.map(s=>`<li>${s.strand===1?'→':'←'} ${s.segments.map(([a,b])=>`${a+1}–${b}`).join(' + ')} bp${s.wraps_origin?' · crosses origin':''}${s.start===q.binding_start&&s.end===q.binding_end&&s.strand===(q.direction==='forward'?1:-1)?' · intended':''}</li>`).join('')}</ul></details></article>`;
  }
  let selected=null;
  function clear(){selected=null;const panel=document.querySelector('#primer-review');panel.replaceChildren();panel.hidden=true;}
  function render(candidates){
    const panel=document.querySelector('#primer-review'),pairs=pairsFrom(candidates);
    if(!pairs.length){clear();return;}
    if(!pairs.some(p=>p.key===selected))selected=pairs[0].key;
    panel.hidden=false;
    panel.innerHTML=`<header><div><span class="eyebrow">Candidate evidence</span><h3>Primer pair & product</h3></div><label>Review pair<select id="primer-pair-choice">${pairs.map(p=>`<option value="${esc(p.key)}">${esc(p.left.name)} / ${esc(p.right.name)} · run ${esc(p.left.metrics.analysis_job_id.slice(0,8))}</option>`).join('')}</select></label></header><div class="primer-pair-detail"></div><p class="evidence-note">Exact, full-length matches on this reference only; both directions counted. Multiple sites warrant review. A single match does not establish specificity. Mismatch binding, other references and alternative PCR products are not evaluated. Site ranges are 1-based inclusive.</p>`;
    const choice=panel.querySelector('select');choice.value=selected;
    function show(){selected=choice.value;const pair=pairs.find(p=>p.key===selected);panel.querySelector('.primer-pair-detail').innerHTML=`<p class="primer-product-summary">Expected product ${pair.left.binding_start+1}–${pair.right.binding_end} bp · ${pair.size} bp long</p><p class="primer-scroll-hint">Scroll the diagram horizontally to inspect both primers.</p><div class="primer-product-scroll" tabindex="0" aria-label="Product diagram; scroll horizontally on small screens">${diagram(pair)}</div><div class="primer-site-cards">${siteCard(pair.left)}${siteCard(pair.right)}</div>`;}
    choice.onchange=show;show();
  }
  const api={pairsFrom,diagram,render,clear};
  if(typeof module!=='undefined')module.exports=api;else root.primerReview=api;
})(typeof window!=='undefined'?window:globalThis);
