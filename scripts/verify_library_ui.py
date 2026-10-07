"""Pagination, out-of-order requests and recoverable failures using synthetic records."""
from pathlib import Path
import io,json,os,socket,subprocess,sys,time,zipfile
import httpx
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from playwright.sync_api import sync_playwright,expect
ROOT=Path(__file__).resolve().parents[1]
artifacts=ROOT/'var'/'library-ui-check'/time.strftime('%Y%m%d-%H%M%S');artifacts.mkdir(parents=True)
archive=artifacts/'library.zip'
with zipfile.ZipFile(archive,'w') as z:
    for name in [f'A{i:03d}' for i in range(53)]+['Slow','Fast']:
        rec=SeqRecord(Seq('ACGT'*40),id=name,name=name,annotations={'molecule_type':'DNA','topology':'linear'})
        stream=io.StringIO();SeqIO.write(rec,stream,'genbank');z.writestr(name+'.gb',stream.getvalue())
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
        page.add_init_script('''
          const realFetch=window.fetch.bind(window);let dashboardCalls=0;
          window.fetch=(url,options={})=>{
            const parsed=new URL(url,location.href);
            if(parsed.pathname==='/api/dashboard' && ++dashboardCalls===1) {
              return realFetch(url,options).then(response=>new Promise(resolve=>setTimeout(()=>{window.initialDashboardFinished=true;resolve(response);},1000)));
            }
            if(parsed.pathname==='/api/sequence-library' && parsed.searchParams.get('query')==='Slow') {
              window.slowSearchStarted=true;
              const opts={...options};delete opts.signal;
              return new Promise(resolve=>setTimeout(resolve,700)).then(()=>realFetch(url,opts)).then(response=>{window.slowSearchFinished=true;return response;});
            }
            if(window.slowPrimerRevision && parsed.pathname===`/api/revisions/${window.slowPrimerRevision}/primers`) {
              window.slowPrimerStarted=true;
              return new Promise(resolve=>setTimeout(()=>{window.slowPrimerFinished=true;resolve(new Response(JSON.stringify([{id:'old',name:'OLD_SENTINEL',direction:'forward',binding_start:1,binding_end:20,sequence_text:'ACGT',metrics:{tm:60,gc_percent:50},selection_state:'candidate'}]),{headers:{'Content-Type':'application/json'}}));},700));
            }
            return realFetch(url,options);
          };
        ''')
        errors=[];page.on('pageerror',lambda e:errors.append(str(e)))
        page.goto(base)
        page.locator('#archive').set_input_files(str(archive));page.get_by_role('button',name='Import into local library').click()
        expect(page.locator('#library-message')).to_have_text('Showing 1–24 of 55 records')
        expect(page.locator('.seq')).to_have_count(24)
        expect(page.locator('#m-seq')).to_have_text('55')
        page.wait_for_function('window.initialDashboardFinished===true')
        expect(page.locator('#m-seq')).to_have_text('55')
        page.locator('#library-next').click();expect(page.locator('#library-message')).to_have_text('Showing 25–48 of 55 records')
        expect(page.locator('.seq strong').first).to_have_text('A024')
        page.locator('#library-next').click();expect(page.locator('#library-message')).to_have_text('Showing 49–55 of 55 records')
        expect(page.locator('.seq')).to_have_count(7);expect(page.locator('#library-next')).to_be_disabled()
        page.locator('#search').fill('Slow');page.wait_for_function('window.slowSearchStarted===true')
        page.locator('#search').fill('Fast');expect(page.locator('.seq strong')).to_have_text('Fast')
        page.wait_for_function('window.slowSearchFinished===true');expect(page.locator('.seq strong')).to_have_text('Fast')
        expect(page.locator('#library-page')).to_have_text('Page 1 of 1')
        page.route('**/api/sequence-library?*',lambda route:route.abort(),times=1)
        page.locator('#search').fill('A000')
        expect(page.locator('#library-message')).to_have_text('Sequence library could not be loaded.')
        page.get_by_role('button',name='Retry loading').click();expect(page.locator('.seq strong')).to_have_text('A000')
        page.locator('#search').fill('');expect(page.locator('#library-message')).to_have_text('Showing 1–24 of 55 records')
        page.screenshot(path=str(artifacts/'library-desktop.png'),full_page=True)
        page.route('**/api/sequences/*/map',lambda route:route.abort(),times=1)
        page.locator('.seq').filter(has_text='A000').get_by_role('button',name='Open map').click()
        expect(page.locator('#map-title')).to_have_text('Construct unavailable')
        expect(page.locator('#primer-form button')).to_be_disabled()
        page.get_by_role('button',name='Retry',exact=True).click();expect(page.locator('#map-title')).to_have_text('A000')
        page.locator('#close-map').click()
        first=httpx.get(base+'/api/sequences',params={'query':'A000'}).json()[0]
        revision=httpx.get(base+f"/api/sequences/{first['id']}/revisions").json()[0]['id']
        page.evaluate('(id)=>window.slowPrimerRevision=id',revision)
        page.locator('.seq').filter(has_text='A000').get_by_role('button',name='Open map').click();page.wait_for_function('window.slowPrimerStarted===true')
        page.locator('#close-map').click()
        page.locator('.seq').filter(has_text='A001').get_by_role('button',name='Open map').click()
        expect(page.locator('#map-title')).to_have_text('A001')
        page.wait_for_function('window.slowPrimerFinished===true')
        expect(page.locator('#primer-results')).not_to_contain_text('OLD_SENTINEL')
        expect(page.locator('#primer-results')).to_contain_text('No candidates saved')
        page.locator('#close-map').click()
        page.route('**/api/imports/benchling',lambda route:route.fulfill(status=500,content_type='application/json',body='{"detail":"Import temporarily unavailable"}'),times=1)
        page.locator('#archive').set_input_files(str(archive));page.get_by_role('button',name='Import into local library').click()
        expect(page.locator('#status')).to_contain_text('Import failed: Import temporarily unavailable')
        expect(page.locator('#archive')).to_be_enabled()
        page.locator('#page-size').select_option('48');expect(page.locator('.seq')).to_have_count(48)
        expect(page.locator('#library-message')).to_have_text('Showing 1–48 of 55 records')
        page.set_viewport_size({'width':390,'height':844})
        page.locator('.library-pagination').scroll_into_view_if_needed()
        page.evaluate('new Promise(resolve=>requestAnimationFrame(()=>requestAnimationFrame(resolve)))')
        page.screenshot(path=str(artifacts/'library-mobile.png'))
        assert page.evaluate('document.documentElement.scrollWidth<=innerWidth+1')
        assert not errors,errors
        browser.close()
    print(json.dumps({'passed':True,'artifacts':str(artifacts),'checks':['55 records across 3 pages','page-size change','late search response','late primer response','network failure and retry','failed import unlocks controls','390px layout','no JS exceptions']},indent=2))
finally:
    server.terminate();server.wait(timeout=10);log.close()
