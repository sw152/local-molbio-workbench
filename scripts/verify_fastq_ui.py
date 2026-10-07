"""Real browser FASTQ uploads with synthetic bytes, failures and mobile review."""
from pathlib import Path
import gzip,io,json,os,socket,subprocess,sys,time,zipfile
import httpx
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from playwright.sync_api import sync_playwright,expect
ROOT=Path(__file__).resolve().parents[1]
artifacts=ROOT/'var'/'fastq-ui-check'/time.strftime('%Y%m%d-%H%M%S');artifacts.mkdir(parents=True)
archive=artifacts/'library.zip'
with zipfile.ZipFile(archive,'w') as z:
    for name in ['FASTQ demo A','FASTQ demo B']:
        rec=SeqRecord(Seq('ACGT'*150),id=name.replace(' ','_'),name=name.replace(' ','_'),annotations={'molecule_type':'DNA','topology':'linear'})
        stream=io.StringIO();SeqIO.write(rec,stream,'genbank');z.writestr(name+'.gb',stream.getvalue())
data=b'@first\nACGN\n+first\n!45?\n@second\nry\n+\n@+\n'
plain=artifacts/'demo.fastq';plain.write_bytes(data)
zipped=artifacts/'instrument_export_with_a_long_filename_for_mobile_review.fastq.gz';zipped.write_bytes(gzip.compress(data,mtime=0))
bad=artifacts/'truncated.fastq';bad.write_bytes(b'@bad\nACGT\n+\n!!\n')
with socket.socket() as s:s.bind(('127.0.0.1',0));port=s.getsockname()[1]
base=f'http://127.0.0.1:{port}';log=(artifacts/'server.log').open('w')
server=subprocess.Popen([sys.executable,'-m','uvicorn','localmolbio.api:app','--host','127.0.0.1','--port',str(port)],cwd=ROOT,env=dict(os.environ,MOLBIO_DATA_DIR=str(artifacts/'runtime')),stdout=log,stderr=log)
try:
    for _ in range(100):
        try:
            if httpx.get(base+'/health').status_code==200:break
        except httpx.TransportError:pass
        time.sleep(.1)
    else:raise RuntimeError('Test server did not start')
    with sync_playwright() as p:
        chrome=Path('/Applications/Google Chrome.app/Contents/MacOS/Google Chrome')
        browser=p.chromium.launch(headless=True,**({'executable_path':str(chrome)} if chrome.exists() else {}))
        page=browser.new_page(viewport={'width':1440,'height':1080},device_scale_factor=1)
        errors=[];page.on('pageerror',lambda error:errors.append(str(error)))
        page.goto(base);page.locator('#archive').set_input_files(str(archive));page.get_by_role('button',name='Import into local library').click()
        expect(page.locator('.open')).to_have_count(2);page.locator('.open').first.click()
        status=page.locator('#fq-status');expect(page.locator('.fq-empty')).to_be_visible()
        sequence_id=page.locator('.open').first.get_attribute('data-id')
        revision=page.request.get(base+f'/api/sequences/{sequence_id}/map').json()['current_revision_id']
        endpoint=base+f'/api/revisions/{revision}/fastq-inputs'
        def attach(path):
            page.locator('#fq-file').set_input_files(str(path));page.locator('#fq-encoding').check();page.locator('#fq-submit').click()
        def count():return page.request.get(endpoint).json()['total']
        page.locator('#fq-file').set_input_files(str(plain));expect(page.locator('#fq-submit')).to_be_disabled()
        page.locator('#fq-encoding').check();page.locator('#fq-submit').click()
        expect(status).to_contain_text('FASTQ registered');assert count()==1
        expect(page.locator('#fq-encoding')).not_to_be_checked();expect(page.locator('#fq-file')).to_have_value('')
        card=page.locator('.fq-card').first
        expect(card.locator('.fq-metrics strong')).to_have_text(['2','6','3 bp'])
        expect(card.locator('.fq-quality li b')).to_have_text(['50%','16.7%','33.3%'])
        assert card.locator('.fq-bar span').evaluate_all('(nodes)=>nodes.map(n=>parseFloat(n.style.width)).reduce((a,b)=>a+b,0)')==100
        attach(zipped);expect(status).to_contain_text('FASTQ registered');assert count()==2
        page.locator('.fq-card').first.locator('summary').click()
        page.locator('#fastq-inputs').scroll_into_view_if_needed();page.screenshot(path=str(artifacts/'fastq-desktop.png'))
        attach(plain);expect(status).to_contain_text('already registered');assert count()==2
        attach(bad);expect(status).to_contain_text('ends before all quality scores');assert count()==2
        page.set_viewport_size({'width':390,'height':844});status.scroll_into_view_if_needed();page.screenshot(path=str(artifacts/'fastq-error-mobile.png'))
        # Preserve selected file and declaration through failed/retried list refresh.
        def failed_list(route):
            if route.request.method=='GET':route.abort()
            else:route.continue_()
        page.route('**/fastq-inputs?*',failed_list);page.locator('#fq-refresh').click();expect(status).to_contain_text('Could not refresh')
        expect(page.locator('#fq-file')).to_have_value('C:\\fakepath\\truncated.fastq');expect(page.locator('#fq-encoding')).to_be_checked()
        page.unroute('**/fastq-inputs?*',failed_list);page.locator('#fq-refresh').click();expect(status).to_contain_text('Registered input summaries')
        def integrity_conflict(route):
            route.fulfill(status=409,content_type='application/json',body=json.dumps({'detail':'existing_fastq_storage_integrity_mismatch'}))
        page.route('**/fastq-inputs',integrity_conflict);attach(plain);expect(status).to_contain_text('stored copy has changed');assert count()==2
        page.unroute('**/fastq-inputs',integrity_conflict)
        # Server committed but browser loses the response: do not claim rejection.
        lost=artifacts/'lost-response.fq';lost.write_bytes(data.replace(b'@first',b'@third').replace(b'+first',b'+third'))
        page.evaluate("""() => {const real=window.fetch.bind(window);window.loseUpload=true;window.fetch=async (url,opts)=>{
          const response=await real(url,opts);if(window.loseUpload && opts?.method==='POST' && String(url).endsWith('/fastq-inputs')){
            window.loseUpload=false;window.lostUploadStatus=response.status;throw new TypeError('Simulated response loss');}return response;};}""")
        attach(lost);expect(status).to_contain_text('Upload outcome unknown');assert page.evaluate('window.lostUploadStatus')==201;assert count()==3
        page.locator('#fq-refresh').click();expect(page.locator('.fq-card')).to_have_count(3)
        page.locator('#fq-submit').click();expect(status).to_contain_text('already registered');assert count()==3
        # Add two inputs for actual server pagination (unique original bytes).
        for i in range(2):
            extra=artifacts/f'extra-{i}.fq';extra.write_bytes(f'@extra{i}\nA\n+\n!\n'.encode());attach(extra);expect(status).to_contain_text('FASTQ registered')
        expect(page.locator('#fq-page')).to_have_text('1–3 of 5 files');page.locator('#fq-next').click()
        expect(page.locator('#fq-page')).to_have_text('4–5 of 5 files');expect(page.locator('.fq-card')).to_have_count(2)
        page.locator('.fq-card').first.locator('summary').click();page.locator('.fq-card').first.scroll_into_view_if_needed()
        page.screenshot(path=str(artifacts/'fastq-card-mobile.png'))
        assert page.locator('#map-panel').evaluate('(e)=>e.scrollWidth<=e.clientWidth')
        page.locator('#fq-previous').click();expect(page.locator('#fq-page')).to_have_text('1–3 of 5 files')
        page.locator('#fq-next').scroll_into_view_if_needed();page.screenshot(path=str(artifacts/'fastq-pagination-mobile.png'))
        # Delay the old list response beyond closing and opening a different construct.
        page.evaluate('''() => {const real=window.fetch.bind(window);window.delayFastq=true;window.fetch=(url,opts)=>real(url,opts).then(r=>{
          if(window.delayFastq && String(url).includes('/fastq-inputs?')){window.delayFastq=false;return new Promise(resolve=>setTimeout(()=>{window.oldFastqFinished=true;resolve(r)},1000))}return r;});}''')
        page.locator('#fq-refresh').click();page.locator('#close-map').click();page.locator('.open').nth(1).click()
        expect(page.locator('.fq-empty')).to_be_visible();page.wait_for_function('window.oldFastqFinished===true')
        expect(page.locator('.fq-card')).to_have_count(0);expect(page.locator('#fq-page')).to_have_text('0 files')
        page.locator('#fastq-inputs').scroll_into_view_if_needed();page.screenshot(path=str(artifacts/'fastq-empty-mobile.png'))
        # A committed upload may finish after the user switches constructs.
        page.evaluate("""() => {const real=window.fetch.bind(window);window.fetch=async (url,opts)=>{
          const response=await real(url,opts);if(opts?.method==='POST' && String(url).endsWith('/fastq-inputs')){
            window.lateUploadStatus=response.status;return new Promise(resolve=>setTimeout(()=>{window.oldUploadFinished=true;resolve(response)},1000));}return response;};}""")
        attach(plain);page.wait_for_function('window.lateUploadStatus===201');page.locator('#close-map').click();page.locator('.open').first.click()
        expect(page.locator('#fq-page')).to_have_text('1–3 of 5 files');page.wait_for_function('window.oldUploadFinished===true')
        expect(page.locator('#fq-status')).to_contain_text('Registered input summaries');expect(page.locator('#fq-page')).to_have_text('1–3 of 5 files')
        expect(page.locator('#fq-submit')).to_be_disabled()
        assert not errors,errors
        browser.close()
    (artifacts/'result.json').write_text(json.dumps({'passed':True,'synthetic_fastq_only':True,'checks':['plain/gzip uploads','explicit encoding','exact disjoint quality bands','duplicate recovery','malformed rejection','list failure and retry retains form','lost upload response recovery','pagination','late list and upload isolation','integrity conflict is not duplicate success','mobile overflow'],'screenshots':[p.name for p in artifacts.glob('*.png')]},indent=2))
    print(json.dumps({'artifacts':str(artifacts),'passed':True}))
finally:
    server.terminate();server.wait(timeout=10);log.close()
