/* Geometry uses zero-based boundaries; the annotation list uses one-based bases. */
(function(root) {
  const esc = value => String(value ?? '').replace(/[&<>"']/g, c => ({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[c]));
  const palette = {CDS:'#68d7c4',promoter:'#b69bff',rep_origin:'#f3c274',terminator:'#ed91b7',source:'#7689a7',primer_bind:'#82baff'};
  function color(type) { return palette[type] || '#91aed4'; }
  function overlaps(a,b,length,circular,gap=0) {
    return a.some(x=>b.some(y=>(circular?[-length,0,length]:[0]).some(shift=>x.start < y.end+shift+gap && x.end > y.start+shift-gap)));
  }
  function layout(features,length,circular) {
    const lanes=[], entries=[];
    const prepared=features.map((f,id)=>({id,feature:f,segments:(f.segments||[]).filter(s=>Number.isInteger(s.start)&&Number.isInteger(s.end)&&s.start>=0&&s.end<=length&&s.end>s.start)}));
    prepared.sort((a,b)=>b.segments.reduce((n,s)=>n+s.end-s.start,0)-a.segments.reduce((n,s)=>n+s.end-s.start,0)||a.id-b.id);
    for(const item of prepared) {
      if(!item.segments.length) continue;
      let lane=lanes.findIndex(used=>!overlaps(item.segments,used,length,circular,length*.008));
      if(lane<0) {lane=lanes.length;lanes.push([]);}
      lanes[lane].push(...item.segments);
      entries.push({...item,lane});
    }
    return {entries,lanes:lanes.length,invalid:features.filter((f,id)=>!entries.some(e=>e.id===id)||entries.find(e=>e.id===id).segments.length!==(f.segments||[]).length).length};
  }
  function ticks(length) {
    const raw=length/6, scale=10**Math.floor(Math.log10(raw));
    const step=[1,2,5,10].map(x=>x*scale).find(x=>x>=raw);
    const values=[];
    for(let p=0;p<length;p+=Math.max(1,step)) values.push(p);
    return values.filter(p=>p===0||length-p>=length*.1);
  }
  const pair=p=>p.map(n=>n.toFixed(3)).join(' ');
  function point(p,r,length,c) {const a=p/length*Math.PI*2-Math.PI/2;return[c+Math.cos(a)*r,c+Math.sin(a)*r];}
  function curve(start,end,r,length,c,clockwise) {return `A ${r} ${r} 0 ${Math.abs(end-start)>length/2?1:0} ${clockwise?1:0} ${pair(point(end,r,length,c))}`;}
  function ribbon(s,r,length,c) {
    const width=5,head=Math.min((s.end-s.start)/2,10*length/(2*Math.PI*r));
    let start=s.start,end=s.end;
    if(s.strand===1) end-=head;
    if(s.strand===-1) start+=head;
    // A complete undirected ring requires two arcs per edge.
    if(end-start===length) return `<circle cx="${c}" cy="${c}" r="${r}" fill="none" stroke="currentColor" stroke-width="10"/>`;
    let d=`M ${pair(point(start,r+width,length,c))} ${curve(start,end,r+width,length,c,true)}`;
    if(s.strand===1) d+=` L ${pair(point(s.end,r,length,c))}`;
    d+=` L ${pair(point(end,r-width,length,c))} ${curve(end,start,r-width,length,c,false)}`;
    if(s.strand===-1) d+=` L ${pair(point(s.start,r,length,c))}`;
    return `<path d="${d} Z" fill="currentColor"/>`;
  }
  function linearRibbon(s,y,length) {
    const x=p=>64+p/length*672, a=x(s.start),b=x(s.end),head=Math.min(10,(b-a)/2);
    const points=s.strand===1?[[a,y-6],[b-head,y-6],[b,y],[b-head,y+6],[a,y+6]]:s.strand===-1?[[b,y-6],[a+head,y-6],[a,y],[a+head,y+6],[b,y+6]]:[[a,y-6],[b,y-6],[b,y+6],[a,y+6]];
    return `<polygon points="${points.map(pair).join(' ')}" fill="currentColor"/>`;
  }
  function direction(s) {return s.strand===1?'+':s.strand===-1?'−':'unknown';}
  function featureGraphic(entry,body) {
    const f=entry.feature,description=`${f.label} · ${f.type} · ${entry.segments.map(s=>`${s.start+1}–${s.end} (${direction(s)})`).join(', ')}`;
    return `<g class="map-feature" tabindex="0" role="button" aria-pressed="false" data-feature="${entry.id}" data-lane="${entry.lane}" aria-label="${esc(description)}" style="color:${color(f.type)}"><title>${esc(description)}</title>${body}</g>`;
  }
  function svgMarkup(model,geometry) {
    const length=model.length_bp,circular=model.topology==='circular';
    const outer=150+Math.max(0,geometry.lanes-1)*20;
    const size=Math.max(600,(outer+90)*2),c=size/2;
    let parts=[];
    if(circular) {
      parts.push(`<circle cx="${c}" cy="${c}" r="${outer+20}" fill="#0b1525" stroke="#283955"/><circle cx="${c}" cy="${c}" r="104" fill="#101c30" stroke="#283955"/><text x="${c}" y="${c-9}" class="map-center-title">CIRCULAR</text><text x="${c}" y="${c+18}" class="map-center-detail">${length.toLocaleString()} bp</text>`);
      for(const t of ticks(length)) {
        const a=point(t,outer+25,length,c),b=point(t,outer+32,length,c),label=point(t,outer+51,length,c);
        parts.push(`<path d="M ${pair(a)} L ${pair(b)}" stroke="#7186a8"/><text x="${label[0]}" y="${label[1]}" class="map-tick">${t.toLocaleString()}</text>`);
      }
      geometry.entries.forEach(entry=>parts.push(featureGraphic(entry,entry.segments.map(s=>ribbon(s,150+entry.lane*20,length,c)).join(''))));
      return {width:size,height:size,markup:parts.join('')};
    }
    const height=Math.max(260,142+geometry.lanes*30);
    parts.push(`<text x="64" y="35" class="map-linear-title">${model.topology==='linear'?'LINEAR':'TOPOLOGY UNKNOWN · linear coordinate view'}</text><path d="M 64 80 H 736" stroke="#7186a8" stroke-width="2"/>`);
    const marks=[...ticks(length),length];
    for(const t of marks) {const x=64+t/length*672;parts.push(`<path d="M ${x} 75 V 86" stroke="#7186a8"/><text x="${x}" y="64" class="map-tick">${t.toLocaleString()}</text>`);}
    geometry.entries.forEach(entry=>parts.push(featureGraphic(entry,entry.segments.map(s=>linearRibbon(s,112+entry.lane*30,length)).join(''))));
    return {width:800,height,markup:parts.join('')};
  }
  function clear() {
    const svg=document.getElementById('plasmid-map');
    svg.replaceChildren();svg.onclick=null;svg.onkeydown=null;
    const legend=document.getElementById('feature-legend');legend.replaceChildren();legend.onclick=null;
    for(const id of ['map-zoom-in','map-zoom-out','map-fit','map-focus']) {const button=document.getElementById(id);button.disabled=true;button.onclick=null;}
    for(const id of ['map-selection','map-layout-note']) document.getElementById(id).textContent='';
  }
  function render(model) {
    for(const id of ['map-zoom-in','map-zoom-out','map-fit']) document.getElementById(id).disabled=false;
    const svg=document.getElementById('plasmid-map'),legend=document.getElementById('feature-legend');
    const geometry=layout(model.features,model.length_bp,model.topology==='circular');
    let drawing=svgMarkup(model,geometry);
    let selectedId=null,focused=false;
    svg.setAttribute('viewBox',`0 0 ${drawing.width} ${drawing.height}`);
    svg.setAttribute('aria-label',`${model.topology} sequence map, ${model.length_bp} base pairs`);
    svg.innerHTML=drawing.markup;
    svg.style.aspectRatio=`${drawing.width} / ${drawing.height}`;
    legend.innerHTML=model.features.map((f,id)=>{
      const entry=geometry.entries.find(e=>e.id===id),invalid=!entry||entry.segments.length!==(f.segments||[]).length;
      return `<li><button class="map-legend-button" data-feature="${id}" ${!entry?'disabled':''}><i class="swatch" style="background:${color(f.type)}"></i><span><b>${esc(f.label)}</b><small>${esc(f.type)} · ${(f.segments||[]).map(s=>`${Number(s.start)+1}–${Number(s.end)} (${direction(s)})`).join(', ')||'No coordinates'}${invalid?' · Unrenderable coordinates':''}</small></span></button></li>`;
    }).join('')||'<li>No features in this record.</li>';
    document.getElementById('map-layout-note').textContent=`${geometry.lanes} annotation tracks · ${model.topology==='circular'?'Clockwise':'Left to right'} coordinates · Tick labels are boundary positions (bp).${geometry.lanes>10?' Dense overview: select an annotation and use “Focus selected”.':''}${geometry.invalid?` ${geometry.invalid} feature(s) contain unavailable coordinates.`:''}`;
    document.getElementById('map-selection').textContent='Select a feature in the map or annotation list to inspect its coordinates.';
    const focusButton=document.getElementById('map-focus');
    focusButton.disabled=true;focusButton.textContent='Focus selected';focusButton.setAttribute('aria-pressed','false');
    let zoom=1;
    const applyZoom=()=>{svg.style.width=`${zoom*100}%`;document.getElementById('map-zoom-value').textContent=`${Math.round(zoom*100)}%`;};
    document.getElementById('map-zoom-in').onclick=()=>{zoom=Math.min(16,zoom+.5);applyZoom();};
    document.getElementById('map-zoom-out').onclick=()=>{zoom=Math.max(1,zoom-.5);applyZoom();};
    document.getElementById('map-fit').onclick=()=>{zoom=1;applyZoom();};
    applyZoom();
    function redrawFocus() {
      const visible=focused?{...geometry,entries:geometry.entries.filter(e=>e.id===selectedId).map(e=>({...e,lane:0})),lanes:1}:geometry;
      drawing=svgMarkup(model,visible);
      svg.setAttribute('viewBox',`0 0 ${drawing.width} ${drawing.height}`);
      svg.style.aspectRatio=`${drawing.width} / ${drawing.height}`;
      svg.innerHTML=drawing.markup;
      document.getElementById('map-layout-note').textContent=focused?`Showing 1 of ${model.features.length} annotations · Boundary positions (bp). Use “Show all” to restore every track.`:`${geometry.lanes} annotation tracks · Boundary positions (bp).${geometry.invalid?` ${geometry.invalid} feature(s) contain unavailable coordinates.`:''}${geometry.lanes>10?' Dense overview: select an annotation and use “Focus selected” for a readable view.':''}`;
      focusButton.textContent=focused?'Show all':'Focus selected';
      focusButton.setAttribute('aria-pressed',String(focused));
      if(selectedId!==null) {const selected=svg.querySelector(`[data-feature="${selectedId}"]`);selected?.classList.add('map-selected');selected?.setAttribute('aria-pressed','true');}
      zoom=1;applyZoom();
    }
    focusButton.onclick=()=>{if(selectedId!==null){focused=!focused;redrawFocus();}};
    function select(event) {
      const target=event.target.closest('[data-feature]');
      if(!target||target.disabled) return;
      const id=Number(target.dataset.feature),feature=model.features[id];
      selectedId=id;focusButton.disabled=false;
      if(focused) redrawFocus();
      for(const node of document.querySelectorAll('#plasmid-map [data-feature], #feature-legend [data-feature]')) {
        const selected=Number(node.dataset.feature)===id;
        node.classList.toggle('map-selected',selected);
        node.setAttribute('aria-pressed',String(selected));
      }
      document.getElementById('map-selection').textContent=`${feature.label} · ${feature.type} · ${feature.segments.map(s=>`${s.start+1}–${s.end} bp, strand ${direction(s)}`).join('; ')}${feature.location_text?' · Source: '+feature.location_text:''}`;
    }
    svg.onclick=legend.onclick=select;
    svg.onkeydown=e=>{if(e.key==='Enter'||e.key===' ') {e.preventDefault();select(e);}};
  }
  const api={layout,overlaps,ticks,ribbon,svgMarkup,render,clear};
  if(typeof module!=='undefined'&&module.exports) module.exports=api;
  else root.sequenceMap=api;
})(typeof window!=='undefined'?window:globalThis);
