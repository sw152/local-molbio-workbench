"""Real browser alignment submission, real minimap2 evidence and responsive review."""
from pathlib import Path
import io,json,os,random,socket,subprocess,sys,time,zipfile
import httpx
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from playwright.sync_api import sync_playwright,expect
ROOT=Path(__file__).resolve().parents[1]
APP_PYTHON=os.environ.get('MOLBIO_VERIFY_PYTHON',sys.executable)
artifacts=ROOT/'var'/'alignment-ui-check'/time.strftime('%Y%m%d-%H%M%S');artifacts.mkdir(parents=True)
if not os.environ.get('MOLBIO_MINIMAP2'):raise RuntimeError('Configure the pinned real minimap2 for this check')
rng=random.Random(617);ref=''.join(rng.choice('ACGT') for _ in range(900))
archive=artifacts/'library.zip'
with zipfile.ZipFile(archive,'w') as z:
    for name in ['Alignment demo A','Alignment demo B']:
        rec=SeqRecord(Seq(ref),id=name.replace(' ','_'),name=name.replace(' ','_'),annotations={'molecule_type':'DNA','topology':'circular'})
        stream=io.StringIO();SeqIO.write(rec,stream,'genbank');z.writestr(name+'.gb',stream.getvalue())
reads=[ref[100:500],str(Seq(ref[350:800]).reverse_complement()),ref[-200:]+ref[:200],'A'*350,ref[80:240]+ref[246:520],ref+ref,ref[50:600],ref[70:570]]
# This repeat permits a one-base left-shift of the same 6-base deletion.
assert ref[80:239]+ref[245:520] == ref[80:240]+ref[246:520]
data=''.join(f'@duplicate\n{r}\n+\n'+('I'*len(r))+'\n' for r in reads).encode()
fastq=artifacts/'synthetic_long_filename_for_alignment_provenance.fastq';fastq.write_bytes(data)
with socket.socket() as s:s.bind(('127.0.0.1',0));port=s.getsockname()[1]
base=f'http://127.0.0.1:{port}';log=(artifacts/'server.log').open('w')
env=dict(os.environ,MOLBIO_DATA_DIR=str(artifacts/'runtime'))
server=subprocess.Popen([APP_PYTHON,'-m','localmolbio','serve','--port',str(port)],cwd=ROOT,env=env,stdout=log,stderr=log)
def worker():
    r=subprocess.run([APP_PYTHON,'-m','localmolbio','work','--adapter','alignment'],cwd=ROOT,env=env,capture_output=True,text=True)
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
        coverage=panel.locator('.aj-coverage')
        expect(coverage.locator('.aj-coverage-legend').locator('div').nth(0)).to_contain_text('900 bp')
        expect(coverage.locator('.aj-coverage-legend').locator('div').nth(2)).to_contain_text('0 bp')
        expect(coverage.locator('.aj-regions')).to_contain_text('No regions in this category')
        expect(coverage.locator('.aj-coverage-counts')).to_contain_text('6 positions with deletion evidence')
        coverage.locator('#aj-region-kind').select_option('deletion')
        expect(coverage.locator('.aj-regions')).to_contain_text('[239, 245)')
        coverage.scroll_into_view_if_needed();page.screenshot(path=str(artifacts/'alignment-coverage-desktop.png'))
        coverage.locator('[data-aj-inspect]').click()
        focused=panel.locator('.aj-region-filter')
        expect(focused).to_contain_text('4 matching records / 8 total')
        expect(panel.locator('.aj-read')).to_have_count(3)
        expect(panel.locator('.aj-read').nth(0)).to_contain_text('Record 1')
        expect(panel.locator('.aj-read').nth(1)).to_contain_text('Record 5')
        expect(panel.locator('.aj-read').nth(1).locator('.aj-region-evidence')).to_contain_text('Deletion positions: [239, 245)')
        panel.locator('[data-aj-reads]').last.click();expect(panel.locator('.aj-read')).to_have_count(1)
        expect(panel.locator('.aj-read')).to_contain_text('Record 8')
        expect(focused).to_contain_text('Reference [239, 245)')
        expect(panel.locator('#aj-relation')).to_have_value('either')
        # Export from a filtered second page still includes every record and all regions.
        with page.expect_download() as download:
            panel.locator('[data-aj-export="json"]').click()
        json_path=artifacts/'alignment-report.json';download.value.save_as(str(json_path))
        exported=json.loads(json_path.read_text());assert len(exported['reads'])==8
        assert exported['selection']=='all_task_records' and exported['review_regions']['deletion']==[[239,245]]
        assert exported['reads'][5]['withheld'] and not exported['evidence']['whole_reference_verified']
        with page.expect_download() as download:
            panel.locator('[data-aj-export="html"]').click()
        html_path=artifacts/'alignment-report.html';download.value.save_as(str(html_path))
        report_page=browser.new_page(viewport={'width':1440,'height':1080})
        external=[];report_page.on('request',lambda req:external.append(req.url) if req.url.startswith(('http:','https:')) else None)
        report_page.goto(html_path.as_uri());expect(report_page.locator('.record')).to_have_count(8)
        expect(report_page.get_by_role('heading',name='Local alignment review',exact=True)).to_be_visible()
        expect(report_page.locator('header')).to_contain_text('No whole-plasmid verdict')
        report_page.screenshot(path=str(artifacts/'alignment-report-desktop.png'))
        report_page.locator('.record').nth(4).scroll_into_view_if_needed()
        expect(report_page.locator('.record').nth(4)).to_contain_text('6 deleted')
        report_page.screenshot(path=str(artifacts/'alignment-report-deletion.png'))
        report_page.set_viewport_size({'width':390,'height':844});report_page.evaluate('window.scrollTo(0,0)')
        report_page.screenshot(path=str(artifacts/'alignment-report-mobile.png'))
        assert report_page.locator('body').evaluate('(e)=>e.scrollWidth<=window.innerWidth')
        report_page.emulate_media(media='print');report_page.screenshot(path=str(artifacts/'alignment-report-print.png'))
        assert not external,external
        report_page.close()
        page.route('**/report?*',lambda route:route.fulfill(status=409,json={'detail':'alignment report evidence invalid'}))
        panel.locator('[data-aj-export="html"]').click();expect(status).to_contain_text('alignment report evidence invalid')
        expect(focused).to_contain_text('4 matching records / 8 total')
        expect(panel.locator('.aj-read')).to_contain_text('Record 8')
        page.unroute('**/report?*')
        panel.locator('[data-aj-reads]').first.click();expect(panel.locator('.aj-read')).to_have_count(3)
        focused.scroll_into_view_if_needed();page.screenshot(path=str(artifacts/'alignment-region-desktop.png'))
        panel.locator('#aj-relation').select_option('deletion')
        expect(focused).to_contain_text('1 matching record / 8 total')
        expect(panel.locator('.aj-read')).to_have_count(1)
        expect(panel.locator('.aj-read')).to_contain_text('Record 5')
        panel.locator('[data-aj-refresh]').click();expect(status).to_contain_text('refreshed')
        expect(panel.locator('#aj-relation')).to_have_value('deletion')
        expect(panel.locator('.aj-read')).to_contain_text('Record 5')
        # A failed relation change retains the previously committed filter and evidence.
        page.route('**/reads?*',lambda route:route.fulfill(status=503,json={'detail':'reads temporarily unavailable'}))
        panel.locator('#aj-relation').select_option('paired');expect(status).to_contain_text('reads temporarily unavailable')
        expect(panel.locator('#aj-relation')).to_have_value('deletion')
        expect(panel.locator('.aj-read')).to_contain_text('Record 5')
        page.unroute('**/reads?*')
        panel.locator('#aj-relation').select_option('paired')
        expect(panel.locator('.aj-read')).to_have_count(3);expect(panel.locator('.aj-read').first).to_contain_text('Record 1')
        page.set_viewport_size({'width':390,'height':844});focused.scroll_into_view_if_needed()
        page.screenshot(path=str(artifacts/'alignment-region-mobile.png'))
        assert page.locator('#map-panel').evaluate('(e)=>e.scrollWidth<=e.clientWidth')
        panel.locator('[data-aj-overview]').click();expect(coverage).to_be_in_viewport(ratio=.1)
        panel.locator('[data-aj-all-reads]').click();expect(focused).to_have_count(0)
        expect(panel.locator('.aj-read')).to_have_count(3)
        page.set_viewport_size({'width':1440,'height':1080})
        coverage.locator('#aj-region-kind').select_option('ambiguous_only')
        expect(coverage.locator('.aj-regions')).to_contain_text('No regions in this category')
        coverage.locator('#aj-region-kind').select_option('unpaired')
        expect(status).to_contain_text('evidence loaded')

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
        coverage.locator('#aj-region-kind').select_option('deletion')
        expect(coverage.locator('.aj-regions')).to_contain_text('[239, 245)')
        coverage.scroll_into_view_if_needed();page.screenshot(path=str(artifacts/'alignment-coverage-mobile.png'))
        assert page.locator('#map-panel').evaluate('(e)=>e.scrollWidth<=e.clientWidth')
        # Coverage errors preserve previous committed review; next refresh can recover.
        page.route('**/coverage?*',lambda route:route.fulfill(status=503,json={'detail':'coverage temporarily unavailable'}))
        panel.locator('[data-aj-refresh]').click();expect(status).to_contain_text('coverage temporarily unavailable')
        expect(coverage.locator('.aj-regions')).to_contain_text('[239, 245)')
        page.unroute('**/coverage?*');panel.locator('[data-aj-refresh]').click();expect(status).to_contain_text('refreshed')
        reverse=panel.locator('.aj-read').nth(1);reverse.locator('.aj-source summary').click();reverse.scroll_into_view_if_needed()
        reverse.locator('.aj-diagram').first.hover();page.mouse.wheel(700,0)
        page.wait_for_function('document.querySelectorAll(".aj-read")[1].querySelector(".aj-diagram").scrollLeft>0')
        page.screenshot(path=str(artifacts/'alignment-reverse-mobile.png'))
        assert page.locator('#map-panel').evaluate('(e)=>e.scrollWidth<=e.clientWidth')
        assert reverse.locator('.aj-diagram').first.evaluate('(e)=>e.scrollWidth>e.clientWidth')
        panel.locator('.aj-read').nth(2).scroll_into_view_if_needed();page.screenshot(path=str(artifacts/'alignment-wrap-mobile.png'))
        # Do not download a stale report after the review is closed or revision changes.
        stale_downloads=[];page.on('download',lambda d:stale_downloads.append(d.suggested_filename))
        page.evaluate('''() => {const real=window.fetch.bind(window);window.delayReport=true;window.fetch=async (url,opts)=>{
          const r=await real(url,opts);if(window.delayReport && String(url).includes('/report?')){
            window.delayReport=false;window.reportWaiting=true;
            return new Promise(resolve=>setTimeout(()=>{window.oldReportDone=true;resolve(r)},1200));}return r;};}''')
        panel.locator('[data-aj-export="html"]').click();page.wait_for_function('window.reportWaiting')
        page.locator('#close-map').click();page.locator('.open').nth(1).click();expect(status).to_contain_text('Choose files')
        page.wait_for_function('window.oldReportDone');expect(panel.locator('.aj-read')).to_have_count(0)
        assert stale_downloads==[]
        page.locator('#close-map').click();page.locator('.open').first.click();expect(status).to_contain_text('Choose files')
        panel.locator('[data-aj-detail]').click();expect(panel.locator('.aj-read')).to_have_count(3)
        coverage.locator('#aj-region-kind').select_option('deletion');expect(coverage.locator('[data-aj-inspect]')).to_have_count(1)
        # A filtered response that arrives after switching revisions cannot restore old reads.
        page.evaluate('''() => {const real=window.fetch.bind(window);window.delayRegion=true;window.fetch=async (url,opts)=>{
          const r=await real(url,opts);if(window.delayRegion && String(url).includes('/reads?') && String(url).includes('&start=')){
            window.delayRegion=false;window.regionWaiting=true;
            return new Promise(resolve=>setTimeout(()=>{window.oldRegionDone=true;resolve(r)},1200));}return r;};}''')
        coverage.locator('[data-aj-inspect]').click();page.wait_for_function('window.regionWaiting')
        page.locator('#close-map').click();page.locator('.open').nth(1).click();expect(status).to_contain_text('Choose files')
        page.wait_for_function('window.oldRegionDone');expect(panel.locator('.aj-read')).to_have_count(0)
        expect(panel.locator('.aj-region-filter')).to_have_count(0)
        page.locator('#close-map').click();page.locator('.open').first.click();expect(status).to_contain_text('Choose files')
        panel.locator('#aj-type').select_option('ont-high-accuracy')
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
        subprocess.run([APP_PYTHON,'-c',command],cwd=ROOT,env=env,check=True)
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
    (artifacts/'result.json').write_text(json.dumps({'passed':True,'real_minimap2':True,'synthetic_only':True,'app_python':APP_PYTHON,'checks':['installed CLI serve/work','explicit type','lost response idempotency','real worker','forward/reverse/circular origin evidence','real deletion remains a diagram gap','multi-traversal withheld is distinct from no hit','read/job/input/attempt pagination','cross-page selection','cancel','refresh failure recovery','late detail and submission revision isolation','mobile overflow','coverage composition and overlapping deletion evidence','coverage filter and failure recovery','region read navigation, pagination and original provenance','paired versus deletion filters and refresh persistence','failed filter preserves committed evidence','mobile focused review and return to overview','complete HTML and JSON downloads from filtered page','standalone report desktop/mobile/print and no network assets','report failure preserves evidence','late report does not download after revision switch','late filtered response revision isolation'],'screenshots':[p.name for p in artifacts.glob('*.png')]},indent=2))
    print(json.dumps({'artifacts':str(artifacts),'passed':True}))
finally:
    server.terminate();server.wait(timeout=10);log.close()
