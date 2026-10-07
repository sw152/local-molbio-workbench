"""Real browser alignment submission, real minimap2 evidence and responsive review."""
from pathlib import Path
import io,json,os,random,socket,subprocess,sys,time,zipfile
import httpx
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from playwright.sync_api import sync_playwright,expect
ROOT=Path(__file__).resolve().parents[1]
artifacts=ROOT/'var'/'alignment-ui-check'/time.strftime('%Y%m%d-%H%M%S');artifacts.mkdir(parents=True)
if not os.environ.get('MOLBIO_MINIMAP2'):raise RuntimeError('Configure the pinned real minimap2 for this check')
rng=random.Random(617);ref=''.join(rng.choice('ACGT') for _ in range(900))
archive=artifacts/'library.zip'
with zipfile.ZipFile(archive,'w') as z:
    for name in ['Alignment demo A','Alignment demo B']:
        rec=SeqRecord(Seq(ref),id=name.replace(' ','_'),name=name.replace(' ','_'),annotations={'molecule_type':'DNA','topology':'circular'})
        stream=io.StringIO();SeqIO.write(rec,stream,'genbank');z.writestr(name+'.gb',stream.getvalue())
reads=[ref[100:500],str(Seq(ref[350:800]).reverse_complement()),ref[-200:]+ref[:200],'A'*350,ref[80:240]+ref[246:520],ref+ref]
data=''.join(f'@duplicate\n{r}\n+\n'+('I'*len(r))+'\n' for r in reads).encode()
fastq=artifacts/'synthetic_long_filename_for_alignment_provenance.fastq';fastq.write_bytes(data)
with socket.socket() as s:s.bind(('127.0.0.1',0));port=s.getsockname()[1]
base=f'http://127.0.0.1:{port}';log=(artifacts/'server.log').open('w')
env=dict(os.environ,MOLBIO_DATA_DIR=str(artifacts/'runtime'))
server=subprocess.Popen([sys.executable,'-m','uvicorn','localmolbio.api:app','--host','127.0.0.1','--port',str(port)],cwd=ROOT,env=env,stdout=log,stderr=log)
def worker():
    r=subprocess.run([sys.executable,'-m','localmolbio.worker','once','--adapter','alignment'],cwd=ROOT,env=env,capture_output=True,text=True)
    assert r.returncode==0,(r.stdout,r.stderr)
    assert json.loads(r.stdout)['status']=='succeeded'
try:
    for _ in range(100):
        try:
            if httpx.get(base+'/health').status_code==200:break
        except httpx.TransportError:pass
        time.sleep(.1)
    else:raise RuntimeError('Test server did not start')
    with sync_playwright() as p:
        browser=p.chromium.launch(headless=True,executable_path='/Applications/Google Chrome.app/Contents/MacOS/Google Chrome')
        page=browser.new_page(viewport={'width':1440,'height':1080},device_scale_factor=1)
        errors=[];page.on('pageerror',lambda error:errors.append(str(error)))
        page.goto(base);page.locator('#archive').set_input_files(str(archive));page.get_by_role('button',name='Import into local library').click()
        expect(page.locator('.open')).to_have_count(2);page.locator('.open').first.click()
        panel=page.locator('#alignment-jobs');status=panel.locator('.ij-status')
        expect(status).to_contain_text('Choose files')
        page.locator('#fq-file').set_input_files(str(fastq));page.locator('#fq-encoding').check();page.locator('#fq-submit').click()
        expect(page.locator('#fq-status')).to_contain_text('FASTQ registered')
        panel.locator('[data-aj-refresh]').click();expect(panel.locator('[data-aj-input]')).to_have_count(1)
        panel.locator('h3').scroll_into_view_if_needed();page.screenshot(path=str(artifacts/'alignment-form-desktop.png'))
        panel.locator('[data-aj-input]').check();expect(panel.locator('[data-aj-submit]')).to_be_disabled()
        panel.locator('#aj-type').select_option('ont-high-accuracy');expect(panel.locator('[data-aj-submit]')).to_be_enabled()
        # Lose only the response, after the real server has accepted the original request.
        page.evaluate('''() => {const real=window.fetch.bind(window);window.loseAlignment=true;window.fetch=async (url,opts)=>{
          const r=await real(url,opts);if(window.loseAlignment && opts?.method==='POST' && String(url).endsWith('/alignments')){
            window.loseAlignment=false;window.acceptedAlignment=r.status;throw new TypeError('Simulated response loss');}return r;};}''')
        panel.locator('[data-aj-submit]').click();expect(status).to_contain_text('same selection and type')
        assert page.evaluate('window.acceptedAlignment')==201
        panel.locator('[data-aj-submit]').click();expect(status).to_contain_text('Alignment queued')
        expect(panel.locator('.ij-job')).to_have_count(1)
        job=panel.locator('[data-aj-detail]').get_attribute('data-aj-detail')
        sequence=page.locator('.open').first.get_attribute('data-id')
        revision=page.request.get(base+f'/api/sequences/{sequence}/map').json()['current_revision_id']
        endpoint=base+f'/api/revisions/{revision}/alignments'
        assert page.request.get(endpoint).json()['total']==1
        worker();panel.locator('[data-aj-refresh]').click();expect(panel.locator('.ij-succeeded')).to_have_text('Local evidence ready')
        panel.locator('[data-aj-detail]').click();expect(panel.locator('.aj-read')).to_have_count(3)
        expect(panel.locator('.aj-boundary')).to_contain_text('No consensus or whole-plasmid verification')
        expect(panel.locator('.aj-read').nth(0)).to_contain_text('[100, 500)')
        expect(panel.locator('.aj-read').nth(1)).to_contain_text('Reverse')
        expect(panel.locator('.aj-read').nth(2)).to_contain_text('crosses origin')
        expect(panel.locator('.aj-read').nth(2)).to_contain_text('[700, 900) → [0, 200)')
        panel.locator('.aj-read').nth(2).scroll_into_view_if_needed();page.screenshot(path=str(artifacts/'alignment-wrap-desktop.png'))
        panel.locator('[data-aj-reads]').last.click();expect(panel.locator('.aj-read')).to_have_count(3)
        expect(panel.locator('.aj-read').first).to_contain_text('Record 4');expect(panel.locator('.aj-read').first).to_contain_text('No alignment reported')
        expect(panel.locator('.aj-read').nth(2)).to_contain_text('Evidence withheld')
        expect(panel.locator('.aj-read').nth(2)).to_contain_text('more than one reference traversal')
        deletion=panel.locator('.aj-read').nth(1)
        expect(deletion).to_contain_text('0 / 6')
        rects=deletion.locator('.aj-diagram rect').evaluate_all('(els)=>els.map(e=>({x:Number(e.getAttribute("x")),w:Number(e.getAttribute("width"))}))')
        assert len(rects)>=2 and any(b['x']>a['x']+a['w']+1 for a,b in zip(rects,rects[1:])),rects
        deletion.scroll_into_view_if_needed();page.screenshot(path=str(artifacts/'alignment-deletion-desktop.png'))
        panel.locator('[data-aj-reads]').first.click();expect(panel.locator('.aj-read')).to_have_count(3)
        # Preserve displayed evidence during refresh failure, then recover.
        page.route('**/alignments?*',lambda route:route.fulfill(status=503,json={'detail':'temporary unavailable'}))
        panel.locator('[data-aj-refresh]').click();expect(status).to_contain_text('Displayed evidence may be out of date')
        expect(panel.locator('.aj-read')).to_have_count(3)
        page.set_viewport_size({'width':390,'height':844});status.scroll_into_view_if_needed();page.screenshot(path=str(artifacts/'alignment-error-mobile.png'))
        page.unroute('**/alignments?*');panel.locator('[data-aj-refresh]').click();expect(status).to_contain_text('refreshed')
        reverse=panel.locator('.aj-read').nth(1);reverse.locator('.aj-source summary').click();reverse.scroll_into_view_if_needed()
        reverse.locator('.aj-diagram').first.hover();page.mouse.wheel(700,0)
        page.wait_for_function('document.querySelectorAll(".aj-read")[1].querySelector(".aj-diagram").scrollLeft>0')
        page.screenshot(path=str(artifacts/'alignment-reverse-mobile.png'))
        assert page.locator('#map-panel').evaluate('(e)=>e.scrollWidth<=e.clientWidth')
        assert reverse.locator('.aj-diagram').first.evaluate('(e)=>e.scrollWidth>e.clientWidth')
        panel.locator('.aj-read').nth(2).scroll_into_view_if_needed();page.screenshot(path=str(artifacts/'alignment-wrap-mobile.png'))
        # Input pages retain explicit selection; new tasks can be cancelled and paginated.
        for i in range(5):
            r=page.request.post(base+f'/api/revisions/{revision}/fastq-inputs',multipart={'quality_encoding':'phred33','file':{'name':f'extra-{i}.fq','mimeType':'application/octet-stream','buffer':(f'@extra{i}\n{ref[20:420]}\n+\n'+('I'*400)+'\n').encode()}})
            assert r.status==201,r.text()
        panel.locator('[data-aj-refresh]').click();expect(panel.locator('[data-aj-input]')).to_have_count(5)
        panel.locator('[data-aj-input]').first.check();panel.locator('[data-aj-inputs]').last.click()
        expect(panel.locator('[data-aj-input]')).to_have_count(1);panel.locator('[data-aj-input]').check()
        expect(panel.locator('.ij-select')).to_contain_text('2 files selected across pages')
        panel.locator('[data-aj-inputs]').first.click();expect(panel.locator('[data-aj-input]:checked')).to_have_count(1)
        panel.locator('[data-aj-submit]').click();expect(status).to_contain_text('Alignment queued')
        panel.locator('[data-aj-cancel]').first.click();expect(status).to_contain_text('Alignment cancelled')
        # A task with more attempts than one page, using the public submission API.
        input_id=page.request.get(base+f'/api/revisions/{revision}/fastq-inputs').json()['items'][0]['id']
        r=page.request.post(endpoint,data={'input_ids':[input_id],'data_type':'ont-high-accuracy','idempotency_key':'history','max_attempts':10})
        assert r.status==201,r.text();history=r.json()['id']
        command='from localmolbio import job_queue as q; from localmolbio.fastq_alignment import ADAPTER\nfor _ in range(4):\n c=q.claim("history-test",ADAPTER); q.fail(c["job_id"],c["lease_token"],"retry",True)\n'
        subprocess.run([sys.executable,'-c',command],cwd=ROOT,env=env,check=True)
        panel.locator('[data-aj-refresh]').click();expect(status).to_contain_text('refreshed')
        panel.locator(f'[data-aj-detail="{history}"]').click()
        expect(panel.locator('.aj-details ol li')).to_have_count(3)
        panel.locator('[data-aj-attempts]').last.click();expect(panel.locator('.aj-details ol li')).to_have_count(1)
        expect(panel.locator('.aj-details ol li')).to_contain_text('Attempt 1')
        panel.locator('[data-aj-attempts]').first.click();expect(panel.locator('.aj-details ol li').first).to_contain_text('Attempt 4')
        panel.locator('[data-aj-attempts]').last.scroll_into_view_if_needed();page.screenshot(path=str(artifacts/'alignment-attempts-mobile.png'))
        page.request.post(endpoint+'/'+history+'/cancel')
        r=page.request.post(endpoint,data={'input_ids':[input_id],'data_type':'short-single','idempotency_key':'fourth'})
        assert r.status==201
        panel.locator('[data-aj-refresh]').click();expect(status).to_contain_text('refreshed')
        panel.locator('[data-aj-jobs]').last.click();expect(panel.locator('.ij-job')).to_have_count(1)
        panel.locator('[data-aj-jobs]').first.click();expect(panel.locator('.ij-job')).to_have_count(3)
        # Old detail and submission responses must not contaminate a newly opened revision.
        page.evaluate('''() => {const real=window.fetch.bind(window);window.delayDetail=true;window.fetch=async (url,opts)=>{
          const r=await real(url,opts);if(window.delayDetail && /alignments\\/[^/]+\\?/.test(String(url))){window.delayDetail=false;
            window.detailWaiting=true;return new Promise(resolve=>setTimeout(()=>{window.oldDetailDone=true;resolve(r)},1200));}return r;};}''')
        panel.locator('[data-aj-detail]').first.click();page.wait_for_function('window.detailWaiting')
        page.locator('#close-map').click();page.locator('.open').nth(1).click();expect(status).to_contain_text('Choose files')
        page.wait_for_function('window.oldDetailDone');expect(panel.locator('.ij-job')).to_have_count(0);expect(panel.locator('.aj-read')).to_have_count(0)
        page.locator('#close-map').click();page.locator('.open').first.click();expect(status).to_contain_text('Choose files')
        panel.locator('[data-aj-input]').first.check();panel.locator('#aj-type').select_option('pacbio-hifi')
        page.evaluate('''() => {const real=window.fetch.bind(window);window.fetch=async (url,opts)=>{const r=await real(url,opts);
          if(opts?.method==='POST' && String(url).endsWith('/alignments')){window.postWaiting=r.status;
            return new Promise(resolve=>setTimeout(()=>{window.oldPostDone=true;resolve(r)},1200));}return r;};}''')
        panel.locator('[data-aj-submit]').click();page.wait_for_function('window.postWaiting===201')
        page.locator('#close-map').click();page.locator('.open').nth(1).click();expect(status).to_contain_text('Choose files')
        page.wait_for_function('window.oldPostDone');expect(panel.locator('.ij-job')).to_have_count(0)
        expect(panel.locator('#aj-type')).to_have_value('');expect(panel.locator('[data-aj-submit]')).to_be_disabled()
        panel.scroll_into_view_if_needed();page.screenshot(path=str(artifacts/'alignment-empty-mobile.png'))
        assert not errors,errors
        browser.close()
    (artifacts/'result.json').write_text(json.dumps({'passed':True,'real_minimap2':True,'synthetic_only':True,'checks':['explicit type','lost response idempotency','real worker','forward/reverse/circular origin evidence','real deletion remains a diagram gap','multi-traversal withheld is distinct from no hit','read/job/input/attempt pagination','cross-page selection','cancel','refresh failure recovery','late detail and submission revision isolation','mobile overflow'],'screenshots':[p.name for p in artifacts.glob('*.png')]},indent=2))
    print(json.dumps({'artifacts':str(artifacts),'passed':True}))
finally:
    server.terminate();server.wait(timeout=10);log.close()
