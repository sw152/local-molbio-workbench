"""Opt-in foreground analysis: submit in the UI without a separate worker command."""
from pathlib import Path
import io,json,os,random,socket,subprocess,sys,time,zipfile
import httpx
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from playwright.sync_api import sync_playwright,expect
ROOT=Path(__file__).resolve().parents[1]
APP_PYTHON=os.environ.get('MOLBIO_VERIFY_PYTHON',sys.executable)
artifacts=ROOT/'var'/'managed-ui-check'/time.strftime('%Y%m%d-%H%M%S');artifacts.mkdir(parents=True)
rng=random.Random(891);ref=''.join(rng.choice('ACGT') for _ in range(900))
record=SeqRecord(Seq(ref),id='managed-demo',name='managed-demo',annotations={'molecule_type':'DNA','topology':'circular'})
stream=io.StringIO();SeqIO.write(record,stream,'genbank')
archive=artifacts/'library.zip'
with zipfile.ZipFile(archive,'w') as z:z.writestr('managed-demo.gb',stream.getvalue())
reads=[ref[50:500],str(Seq(ref[200:800]).reverse_complement()),ref[-180:]+ref[:180]]
fastq=artifacts/'synthetic.fq';fastq.write_text(''.join(f'@record\n{r}\n+\n'+('I'*len(r))+'\n' for r in reads))
with socket.socket() as sock:sock.bind(('127.0.0.1',0));port=sock.getsockname()[1]
base=f'http://127.0.0.1:{port}'
env=dict(os.environ,MOLBIO_DATA_DIR=str(artifacts/'runtime'))
log=(artifacts/'server.log').open('w')
server=subprocess.Popen([APP_PYTHON,'-m','localmolbio','serve','--port',str(port),'--with-alignment-worker'],cwd=ROOT,env=env,stdout=log,stderr=log)
def wait(check,seconds=20):
    end=time.monotonic()+seconds
    while time.monotonic()<end:
        try:
            if check():return
        except httpx.TransportError:pass
        time.sleep(.1)
    raise AssertionError('condition did not become true')
try:
    wait(lambda:httpx.get(base+'/health').status_code==200)
    with sync_playwright() as p:
        browser=p.chromium.launch(headless=True,executable_path='/Applications/Google Chrome.app/Contents/MacOS/Google Chrome')
        page=browser.new_page(viewport={'width':1440,'height':1080});errors=[];page.on('pageerror',lambda e:errors.append(str(e)))
        page.goto(base);page.locator('#archive').set_input_files(str(archive));page.get_by_role('button',name='Import into local library').click()
        expect(page.locator('.open')).to_have_count(1);page.locator('.open').click()
        panel=page.locator('#alignment-jobs');status=panel.locator('.ij-status');mode=panel.locator('.aj-execution')
        expect(mode).to_contain_text('Analysis mode is on')
        page.locator('#fq-file').set_input_files(str(fastq));page.locator('#fq-encoding').check();page.locator('#fq-submit').click()
        expect(page.locator('#fq-status')).to_contain_text('FASTQ registered')
        panel.locator('[data-aj-refresh]').click();expect(panel.locator('[data-aj-input]')).to_have_count(1)
        panel.locator('[data-aj-input]').check();panel.locator('#aj-type').select_option('ont-high-accuracy');panel.locator('[data-aj-submit]').click()
        expect(status).to_contain_text('Alignment queued')
        sequence=page.locator('.open').get_attribute('data-id')
        revision=page.request.get(base+f'/api/sequences/{sequence}/map').json()['current_revision_id']
        endpoint=base+f'/api/revisions/{revision}/alignments'
        wait(lambda:httpx.get(endpoint).json()['items'][0]['status']=='succeeded')
        panel.locator('[data-aj-refresh]').click();expect(panel.locator('.ij-succeeded')).to_have_count(1)
        expect(mode).to_contain_text('Analysis mode is on')
        mode.scroll_into_view_if_needed();page.screenshot(path=str(artifacts/'managed-mode-desktop.png'))
        panel.locator('[data-aj-detail]').click();expect(panel.locator('.aj-read')).to_have_count(3)
        expect(panel.locator('.aj-read').nth(1)).to_contain_text('Reverse')
        expect(panel.locator('.aj-read').nth(2)).to_contain_text('crosses origin')
        with page.expect_download() as download:panel.locator('[data-aj-export="json"]').click()
        report=artifacts/'report.json';download.value.save_as(str(report));assert len(json.loads(report.read_text())['reads'])==3
        # Runtime status fetch failure cannot be mistaken for idle/active execution.
        page.route('**/api/runtime',lambda route:route.fulfill(status=503,json={'detail':'unavailable'}))
        panel.locator('[data-aj-refresh]').click();expect(mode).to_contain_text('Execution status unavailable')
        page.unroute('**/api/runtime');panel.locator('[data-aj-refresh]').click();expect(mode).to_contain_text('Analysis mode is on')
        # Unsupported alignment alphabet is a real terminal task failure, stopping the runner.
        bad=ref[:100]+'R'+ref[101:500]
        response=page.request.post(base+f'/api/revisions/{revision}/fastq-inputs',multipart={'quality_encoding':'phred33','file':{'name':'synthetic-unsupported.fq','mimeType':'application/octet-stream','buffer':(f'@bad\n{bad}\n+\n'+('I'*len(bad))+'\n').encode()}})
        assert response.status==201,response.text()
        item=response.json()['id']
        response=page.request.post(endpoint,data={'input_ids':[item],'data_type':'ont-high-accuracy','idempotency_key':'bad'})
        assert response.status==201
        wait(lambda:httpx.get(base+'/api/runtime').json()['alignment_worker']['state']=='failed')
        panel.locator('[data-aj-refresh]').click();expect(mode).to_contain_text('stopped after a worker failure')
        expect(panel.locator('.ij-failed')).to_have_count(1)
        good=page.request.get(base+f'/api/revisions/{revision}/fastq-inputs').json()['items']
        good_id=next(i['id'] for i in good if i['original_filename']=='synthetic.fq')
        pending=page.request.post(endpoint,data={'input_ids':[good_id],'data_type':'ont-high-accuracy','idempotency_key':'after-failure'}).json()['id']
        panel.locator('[data-aj-refresh]').click();expect(panel.locator('.ij-queued')).to_have_count(1)
        assert httpx.get(endpoint+'/'+pending).json()['status']=='queued'
        page.set_viewport_size({'width':390,'height':844});mode.scroll_into_view_if_needed()
        page.screenshot(path=str(artifacts/'managed-failure-mobile.png'))
        assert page.locator('#map-panel').evaluate('(e)=>e.scrollWidth<=e.clientWidth')
        assert not errors,errors
        browser.close()
    (artifacts/'result.json').write_text(json.dumps({'passed':True,'synthetic_only':True,'real_minimap2':True,'app_python':APP_PYTHON,'separate_worker_command':False,'checks':['UI submission to automatic execution','reverse and circular evidence','report download','runtime unavailable recovery','first failure stops execution','next task stays queued','desktop/mobile visual review']},indent=2))
    print(json.dumps({'passed':True,'artifacts':str(artifacts)}))
finally:
    server.terminate();server.wait(timeout=15);log.close()
