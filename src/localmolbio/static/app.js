(() => {
  const $=s=>document.querySelector(s);
  const esc=v=>String(v ?? '').replace(/[&<>'"]/g,c=>({'&':'&amp;','<':'&lt;','>':'&gt;',"'":'&#39;','"':'&quot;'}[c]));
  let revision=null,mapVersion=0,mapController=null,primerBusy=false,primerVersion=0,dashboardVersion=0;
  const library={offset:0,limit:24,total:0,version:0,controller:null,loading:false};
  async function request(url,options) {
    const response=await fetch(url,options);
    let data;
    try {data=await response.json();} catch {throw new Error(`Unexpected server response (${response.status})`);}
    if(!response.ok) throw new Error(typeof data.detail==='string'?data.detail:Array.isArray(data.detail)?data.detail.map(e=>e.msg).join('; '):`Request failed (${response.status})`);
    return data;
  }
  async function dashboard() {
    const version=++dashboardVersion;
    const metrics=['#m-seq','#m-circ','#m-review','#m-primer'];
    try {
      const d=await request('/api/dashboard');
      if(version!==dashboardVersion)return;
      [d.sequence_count,d.circular_count,d.review_count,d.selected_primer_count].forEach((v,i)=>{$(metrics[i]).textContent=Number(v).toLocaleString();$(metrics[i]).removeAttribute('title');});
    } catch(error) {if(version!==dashboardVersion)return;metrics.forEach(id=>{$(id).textContent='—';$(id).title=`Could not refresh summary: ${error.message}`;});}
  }
  function libraryControls() {
    $('#library-previous').disabled=library.loading||library.offset===0;
    $('#library-next').disabled=library.loading||library.offset+library.limit>=library.total;
    $('#sequence-list').setAttribute('aria-busy',String(library.loading));
  }
  function invalidateLibrary() {
    library.version++;
    library.controller?.abort();
    library.loading=true;
    $('#sequence-list').replaceChildren();
    $('#library-page').textContent='';
    $('#library-message').textContent='Loading sequence records…';
    libraryControls();
  }
  async function list() {
    invalidateLibrary();
    const version=library.version,controller=new AbortController();
    library.controller=controller;
    const params=new URLSearchParams({query:$('#search').value,offset:library.offset,limit:library.limit});
    try {
      const page=await request('/api/sequence-library?'+params,{signal:controller.signal});
      if(version!==library.version)return;
      library.total=page.total;
      const records=page.items;
      $('#library-message').textContent=page.total?`Showing ${page.offset+1}–${page.offset+records.length} of ${Number(page.total).toLocaleString()} records`:'No matching sequence records.';
      $('#library-page').textContent=page.total?`Page ${Math.floor(page.offset/page.limit)+1} of ${Math.ceil(page.total/page.limit)}`:'';
      $('#sequence-list').innerHTML=records.length?records.map(s=>`<article class="seq"><header><strong title="${esc(s.display_name)}">${esc(s.display_name)}</strong><span class="badge ${s.parse_warning_count?'review':''}">${s.parse_warning_count?'Review '+s.parse_warning_count:esc(s.topology)}</span></header><div class="meta"><span>${Number(s.length_bp).toLocaleString()} bp</span><span>${s.feature_count} features</span></div><footer><span class="source" title="${esc(s.archive_member_name)}">${esc(s.archive_member_name)}</span><button class="link open" data-id="${esc(s.id)}">Open map →</button></footer></article>`).join(''):'<div class="empty">Try another name or source filename, or import a sequence library.</div>';
    } catch(error) {
      if(version!==library.version||error.name==='AbortError')return;
      library.total=0;
      $('#library-message').textContent='Sequence library could not be loaded.';
      $('#sequence-list').innerHTML=`<div class="empty library-error"><p>${esc(error.message)}</p><button class="export" data-library-retry>Retry loading</button></div>`;
    } finally {if(version===library.version){library.loading=false;libraryControls();}}
  }
  function primerControls() {
    $('#primer-form').querySelectorAll('input,button').forEach(e=>e.disabled=!revision||primerBusy);
    $('#export-primers').disabled=!revision||primerBusy;
    $('#primer-results').querySelectorAll('button').forEach(e=>e.disabled=primerBusy);
  }
  async function primers(id=revision,version=mapVersion) {
    if(!id)return;
    const requestVersion=++primerVersion;
    try {
      const candidates=await request('/api/revisions/'+encodeURIComponent(id)+'/primers');
      if(version!==mapVersion||id!==revision||requestVersion!==primerVersion)return;
      $('#primer-results').innerHTML=candidates.length?candidates.map(q=>`<div class="primer"><span><b>${esc(q.name||'Unnamed')}</b><small>${esc(q.direction)} · ${q.binding_start+1}..${q.binding_end}</small><small>${q.metrics.target?`Target ${q.metrics.target.start+1}–${q.metrics.target.end} bp`:'Unconstrained placement'} · ${q.metrics.product_size} bp product</small></span><code title="${esc(q.sequence_text)}">${esc(q.sequence_text)}</code><span>${Number(q.metrics.tm).toFixed(1)}°<small>${Number(q.metrics.gc_percent).toFixed(1)}% GC</small>${q.selection_state!=='selected'?`<button class="action" data-id="${esc(q.id)}" data-state="selected">Select</button>`:'<small style="color:var(--cyan)">Selected</small>'}${q.selection_state!=='archived'?`<button class="action" data-id="${esc(q.id)}" data-state="archived">Archive</button>`:''}</span></div>`).join(''):'<p class="evidence-note">No candidates saved for this revision.</p>';
      window.primerReview.render(candidates);
      primerControls();
    } catch(error) {
      if(version===mapVersion&&id===revision&&requestVersion===primerVersion) $('#primer-status').textContent=`Could not load primers: ${error.message}`;
    }
  }
  function resetConstruct() {
    mapVersion++;primerVersion++;mapController?.abort();revision=null;primerBusy=false;
    window.alignmentJobs.clear();window.fastqInputs.clear();window.inputJobs.clear();window.sangerView.clear();window.sequenceMap.clear();window.primerReview.clear();
    $('#primer-results').replaceChildren();$('#primer-status').textContent='';$('#primer-design-summary').textContent='';$('#primer-form [name=target_start]').value='';$('#primer-form [name=target_end]').value='';
    $('#map-title').textContent='Loading construct…';$('#map-summary').textContent='';
    primerControls();
  }
  async function openMap(id) {
    resetConstruct();
    const version=mapVersion,controller=new AbortController();mapController=controller;
    $('#map-panel').hidden=false;$('#map-message').className='alert';$('#map-message').textContent='Loading construct evidence…';
    try {
      const m=await request('/api/sequences/'+encodeURIComponent(id)+'/map',{signal:controller.signal});
      if(version!==mapVersion)return;
      revision=m.current_revision_id;
      $('#primer-form [name=target_start]').max=m.length_bp;$('#primer-form [name=target_end]').max=m.length_bp;
      $('#map-title').textContent=m.display_name;
      $('#map-summary').textContent=`${Number(m.length_bp).toLocaleString()} bp · ${m.topology} · ${m.feature_count} annotated features`;
      $('#map-message').innerHTML=window.sequenceMap.annotationReview(m);
      $('#map-message').className=m.requires_annotation_review?'alert warning':'alert';
      window.sequenceMap.render(m);window.sangerView.open(revision,m.length_bp);window.inputJobs.open(revision);window.fastqInputs.open(revision);window.alignmentJobs.open(revision);primerControls();
      await primers(revision,version);
    } catch(error) {
      if(version!==mapVersion||error.name==='AbortError')return;
      $('#map-title').textContent='Construct unavailable';$('#map-message').className='alert warning';
      $('#map-message').innerHTML=`${esc(error.message)} <button class="map-retry" data-map-retry="${esc(id)}">Retry</button>`;
    }
  }
  let searchTimer;
  $('#search').addEventListener('input',()=>{clearTimeout(searchTimer);library.offset=0;invalidateLibrary();searchTimer=setTimeout(list,250);});
  $('#page-size').addEventListener('change',()=>{clearTimeout(searchTimer);library.limit=Number($('#page-size').value);library.offset=0;list();});
  $('#library-previous').onclick=()=>{library.offset=Math.max(0,library.offset-library.limit);list();};
  $('#library-next').onclick=()=>{library.offset+=library.limit;list();};
  $('#sequence-list').onclick=e=>{const open=e.target.closest('.open');if(open)openMap(open.dataset.id);else if(e.target.closest('[data-library-retry]'))list();};
  $('#map-message').onclick=e=>{const button=e.target.closest('[data-map-retry]');if(button)openMap(button.dataset.mapRetry);};
  $('#close-map').onclick=()=>{resetConstruct();$('#map-panel').hidden=true;};
  $('#import-form').onsubmit=async e=>{
    e.preventDefault();const file=$('#archive').files[0];if(!file)return;
    const button=e.target.querySelector('button');if(button.disabled)return;
    button.disabled=true;$('#archive').disabled=true;
    const body=new FormData();body.append('file',file);$('#status').textContent=`Importing ${file.name}…`;
    try {
      const result=await request('/api/imports/benchling',{method:'POST',body});
      $('#status').textContent=`Imported ${result.imported_records} records. ${result.parse_warning_count||'No'} parser warnings retained.`;
      e.target.reset();clearTimeout(searchTimer);library.offset=0;list();dashboard();
    } catch(error) {$('#status').textContent=`Import failed: ${error.message}`;}
    finally {button.disabled=false;$('#archive').disabled=false;}
  };
  $('#primer-form').onsubmit=async e=>{
    e.preventDefault();if(!revision||primerBusy)return;
    const data=new FormData(e.target),version=mapVersion,id=revision;
    const first=data.get('target_start'),last=data.get('target_end');
    if(Boolean(first)!==Boolean(last)) {$('#primer-status').textContent='Enter both target positions or leave both empty.';return;}
    if(first && Number(first)>Number(last)) {$('#primer-status').textContent='Target first base must not exceed the last base. Origin-spanning targets are not supported.';return;}
    primerBusy=true;primerControls();$('#primer-design-summary').textContent='';$('#primer-status').textContent='Designing primers…';
    try {
      const result=await request('/api/primer-designs',{method:'POST',headers:{'Content-Type':'application/json'},body:JSON.stringify({sequence_revision_id:id,name_prefix:data.get('name_prefix')||null,product_size_min:+data.get('product_size_min'),product_size_max:+data.get('product_size_max'),num_return:+data.get('num_return'),target_start:first?Number(first)-1:null,target_end:last?Number(last):null})});
      if(version!==mapVersion)return;
      $('#primer-status').textContent=result.pairs.length?`Generated ${result.pairs.length} candidate pair(s).`:'No primer pairs met these constraints. Try changing the target or product size range.';
      $('#primer-design-summary').textContent=(result.parameters.target_start===null?'Unconstrained placement':`Flanking target ${result.parameters.target_start+1}–${result.parameters.target_end} bp`)+' · Linear reference coordinates · Specificity not evaluated';await primers(id,version);dashboard();
    } catch(error) {if(version===mapVersion)$('#primer-status').textContent=`Primer design failed: ${error.message}`;}
    finally {if(version===mapVersion){primerBusy=false;primerControls();}}
  };
  $('#primer-results').onclick=async e=>{
    const button=e.target.closest('.action');if(!button||!revision||primerBusy)return;
    const version=mapVersion,id=revision;primerBusy=true;primerControls();
    try {
      await request('/api/primers/'+encodeURIComponent(button.dataset.id),{method:'PATCH',headers:{'Content-Type':'application/json'},body:JSON.stringify({selection_state:button.dataset.state})});
      if(version!==mapVersion)return;
      $('#primer-status').textContent='Primer selection saved.';await primers(id,version);dashboard();
    } catch(error) {if(version===mapVersion)$('#primer-status').textContent=`Could not save selection: ${error.message}`;}
    finally {if(version===mapVersion){primerBusy=false;primerControls();}}
  };
  $('#export-primers').onclick=()=>{if(revision)location.href='/api/revisions/'+encodeURIComponent(revision)+'/primers.csv';};
  request('/health').then(()=>$('#health').textContent='Local service online').catch(()=>$('#health').textContent='Service unavailable');
  primerControls();list();dashboard();
})();
