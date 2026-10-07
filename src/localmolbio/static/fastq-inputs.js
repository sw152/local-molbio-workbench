(() => {
  const root=document.querySelector('#fastq-inputs');
  const $=s=>root.querySelector(s);
  const esc=v=>String(v??'').replace(/[&<>"']/g,c=>({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[c]));
  const number=v=>Number(v).toLocaleString(undefined,{maximumFractionDigits:1});
  let revision=null,generation=0,busy=false,offset=0,total=0,hasMore=false;
  const limit=3;
  function status(message,error=false){$('#fq-status').textContent=message;$('#fq-status').classList.toggle('fq-error',error);}
  function controls(){
    root.setAttribute('aria-busy',String(busy));
    root.querySelectorAll('input,button').forEach(e=>e.disabled=busy);
    $('#fq-submit').disabled=busy||!$('#fq-file').files.length||!$('#fq-encoding').checked;
    $('#fq-previous').disabled=busy||offset===0;
    $('#fq-next').disabled=busy||!hasMore;
  }
  function card(item){
    const s=item.summary,bands=[['Below Q20',s.bases-s.q20_bases,'low'],['Q20–29',s.q20_bases-s.q30_bases,'mid'],['Q30+',s.q30_bases,'high']];
    return `<article class="fq-card"><header><h4>${esc(item.original_filename)}</h4><span class="fq-badge">${item.compression==='gzip'?'GZIP':'FASTQ'} · Phred+33</span></header>
      <div class="fq-metrics"><div><strong>${number(s.records)}</strong><span>Reads</span></div><div><strong>${number(s.bases)}</strong><span>Total bases</span></div><div><strong>${number(s.mean_read_length)} <small>bp</small></strong><span>Mean read length</span></div></div>
      <div class="fq-range"><span>Read lengths <b>${number(s.min_read_length)}–${number(s.max_read_length)} bp</b></span><span>Ambiguous bases <b>${number(s.ambiguous_bases)} <small>(N: ${number(s.n_bases)})</small></b></span></div>
      <div class="fq-quality"><h5>Base quality composition</h5><div class="fq-bar" aria-hidden="true">${bands.map(([,count,kind])=>`<span class="fq-${kind}" style="width:${count/s.bases*100}%"></span>`).join('')}</div><ul>${bands.map(([label,count,kind])=>`<li><i class="fq-${kind}"></i><span>${label}</span><b>${number(count/s.bases*100)}%</b><small>${number(count)} bases</small></li>`).join('')}</ul></div>
      <details><summary>Original file &amp; validation details</summary><dl><dt>File size</dt><dd>${number(item.size_bytes)} bytes</dd><dt>SHA-256 · original bytes</dt><dd><code>${esc(item.file_sha256)}</code></dd><dt>Registered</dt><dd>${esc(item.created_at)}</dd><dt>Phred scores · min / mean / max</dt><dd>${number(s.min_phred)} / ${number(s.mean_phred)} / ${number(s.max_phred)}</dd><dt>Encoding source</dt><dd>Explicit user declaration; not auto-detected</dd><dt>Read layout</dt><dd>Single file; pairing not inferred</dd></dl></details></article>`;
  }
  const errorMessages={
    truncated_quality:'A read ends before all quality scores are present. Check the original export.',
    sequence_quality_length_mismatch:'A read has more quality scores than bases. Check the original export.',
    invalid_or_truncated_gzip:'The gzip file is damaged or incomplete. Export or transfer it again.',
    compression_does_not_match_filename:'The file content does not match its plain FASTQ or gzip extension.',
    invalid_dna_sequence:'A sequence contains unsupported characters. Only DNA IUPAC bases are accepted.',
    invalid_phred33_character:'A quality score is outside the supported Phred+33 character range.',
    invalid_fastq_header:'A record header is missing or malformed.',
    repeated_header_mismatch:'The repeated title after + does not match the read header.',
    truncated_sequence:'A sequence is incomplete or its + separator is missing.',
    empty_read:'The file contains a read with no bases.',empty_fastq:'The file contains no FASTQ records.',
    raw_size_limit_exceeded:'The original file exceeds 2 GiB.',decoded_size_limit_exceeded:'The decoded file exceeds 8 GiB.',
    record_size_limit_exceeded:'A read exceeds the 4 Mi base limit.',line_size_limit_exceeded:'A line exceeds the supported size limit.',
    record_count_limit_exceeded:'The file exceeds 10 million reads.',
    existing_fastq_storage_integrity_mismatch:'The stored copy has changed unexpectedly. It was not overwritten; inspect the stored original before retrying.',
    reference_missing:'This reference revision is no longer available. Reopen the construct.',
    invalid_filename:'The filename is missing or too long.',expected_fastq_or_fastq_gz:'Choose a FASTQ or gzip FASTQ file.',
    quality_encoding_must_be_explicit_phred33:'Explicitly confirm Phred+33 encoding before uploading.'
  };
  async function responseData(response){
    let data;try{data=await response.json();}catch{throw new Error('The server response could not be read.');}
    if(!response.ok){const error=new Error(typeof data.detail==='string'?data.detail:`Request rejected (${response.status}).`);error.code=response.status;throw error;}
    return data;
  }
  async function load(targetOffset=offset,message=''){
    if(!revision||busy)return;
    const version=generation,id=revision;
    busy=true;controls();status('Loading registered FASTQ inputs…');
    try{
      const data=await responseData(await fetch(`/api/revisions/${encodeURIComponent(id)}/fastq-inputs?limit=${limit}&offset=${targetOffset}`));
      if(version!==generation)return;
      offset=data.offset;total=data.total;hasMore=data.has_more;
      $('#fq-list').innerHTML=data.items.length?data.items.map(card).join(''):'<div class="fq-empty">No FASTQ inputs attached.<br><small>Register a file to inspect read lengths and base quality.</small></div>';
      $('#fq-page').textContent=total?`${offset+1}–${offset+data.items.length} of ${number(total)} files`:'0 files';
      status(message||'Registered input summaries. No alignment or consensus has been performed.');
    }catch(error){if(version===generation)status(`Could not refresh inputs: ${error.message} Use Refresh to retry. Displayed summaries may be out of date.`,true);}
    finally{if(version===generation){busy=false;controls();}}
  }
  async function upload(event){
    event.preventDefault();if(!revision||busy)return;
    const file=$('#fq-file').files[0];if(!file||!$('#fq-encoding').checked)return;
    if(!/\.(fastq|fq)(\.gz)?$/i.test(file.name)){status('Choose a .fastq, .fq, .fastq.gz or .fq.gz file.',true);return;}
    if(file.size>2*1024**3){status('File exceeds the 2 GiB original-file limit.',true);return;}
    const version=generation,id=revision,body=new FormData();body.append('file',file);body.append('quality_encoding','phred33');
    busy=true;controls();status('Uploading and validating the entire file… Large files may take time.');
    let message='',refresh=false;
    try{
      await responseData(await fetch(`/api/revisions/${encodeURIComponent(id)}/fastq-inputs`,{method:'POST',body}));
      if(version!==generation)return;
      $('#fq-form').reset();message='FASTQ registered. Original bytes retained; format and quality summary only.';refresh=true;
    }catch(error){
      if(version!==generation)return;
      if(error.code===409&&error.message==='fastq_already_registered_for_revision'){
        message='These original file bytes are already registered for this revision. Refreshed the latest inputs; use pagination to find the existing entry.';refresh=true;
      }else if(!error.code||error.code>=500){
        status('Upload outcome unknown. Refresh inputs before retrying; the server may already have registered the file.',true);
      }else{status(`FASTQ was not registered: ${errorMessages[error.message]||error.message}`,true);}
    }finally{if(version===generation){busy=false;controls();}}
    if(refresh&&version===generation)await load(0,message);
  }
  function clear(){generation++;revision=null;busy=false;offset=0;total=0;hasMore=false;root.replaceChildren();root.removeAttribute('aria-busy');}
  function open(id){
    clear();revision=id;
    root.innerHTML=`<header class="fq-header"><div><span class="eyebrow">Sequencing inputs</span><h3>FASTQ files</h3><p>Keep original reads with this reference revision. Inspect file-level quality before analysis.</p></div><span class="fq-badge">Input summary only</span></header>
      <div class="fq-layout"><form id="fq-form" class="fq-upload"><h4>Register reads</h4><label for="fq-file">FASTQ or compressed FASTQ</label><input id="fq-file" type="file" accept=".fastq,.fq,.fastq.gz,.fq.gz" required>
      <label class="fq-declaration"><input id="fq-encoding" type="checkbox" required><span>I confirm this file uses <b>Phred+33</b> quality encoding.</span></label><p>Check your sequencing provider’s export settings. Encoding cannot be reliably inferred from overlapping score ranges.</p><button id="fq-submit" type="submit" disabled>Upload &amp; inspect</button>
      <details class="fq-limits"><summary>Supported format &amp; limits</summary><p>DNA IUPAC bases, plain FASTQ or gzip. Wrapped records are supported. Up to 2 GiB original / 8 GiB decoded, 10 million reads, 4 Mi bases per read. Empty or malformed files are rejected as a whole.</p><p>Files remain unpaired. No platform detection, trimming, alignment, consensus or whole-plasmid verification is performed.</p></details></form>
      <div class="fq-results"><div class="fq-list-heading"><h4>Registered files</h4><button id="fq-refresh" type="button">Refresh</button></div><p id="fq-status" role="status" aria-live="polite"></p><div id="fq-list"></div><nav class="fq-pagination" aria-label="FASTQ files pagination"><span id="fq-page"></span><button id="fq-previous" type="button">← Previous</button><button id="fq-next" type="button">Next →</button></nav></div></div>`;
    $('#fq-form').onsubmit=upload;$('#fq-form').onchange=controls;
    $('#fq-refresh').onclick=()=>load();$('#fq-previous').onclick=()=>load(Math.max(0,offset-limit));$('#fq-next').onclick=()=>load(offset+limit);
    load(0);
  }
  window.fastqInputs={open,clear};
})();
