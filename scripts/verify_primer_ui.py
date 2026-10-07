"""Real Primer3 browser validation on an explicitly synthetic random template."""
from pathlib import Path
import io,json,os,random,socket,subprocess,sys,time,zipfile
import httpx
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature,FeatureLocation
from playwright.sync_api import sync_playwright,expect
ROOT=Path(__file__).resolve().parents[1]
artifacts=ROOT/'var'/'primer-ui-check'/time.strftime('%Y%m%d-%H%M%S');artifacts.mkdir(parents=True)
rng=random.Random(42);sequence=''.join(rng.choice('ACGT') for _ in range(1800))
record=SeqRecord(Seq(sequence),id='target-demo',name='target-demo',annotations={'molecule_type':'DNA','topology':'linear'})
record.features=[SeqFeature(FeatureLocation(650,850,strand=1),type='CDS',qualifiers={'label':['Target region']})]
stream=io.StringIO();SeqIO.write(record,stream,'genbank');archive=artifacts/'target.zip'
with zipfile.ZipFile(archive,'w') as z:z.writestr('target.gb',stream.getvalue())
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
        errors=[];page.on('pageerror',lambda e:errors.append(str(e)))
        page.goto(base);page.locator('#archive').set_input_files(str(archive))
        page.get_by_role('button',name='Import into local library').click()
        page.get_by_role('button',name='Open map').click()
        form=page.locator('#primer-form')
        expect(form.locator('[name=target_start]')).to_be_enabled()
        form.locator('[name=name_prefix]').fill('TargetPCR')
        form.locator('[name=product_size_min]').fill('300');form.locator('[name=product_size_max]').fill('500')
        form.locator('[name=target_start]').fill('651');form.get_by_role('button',name='Generate',exact=True).click()
        expect(page.locator('#primer-status')).to_contain_text('Enter both target positions')
        form.locator('[name=target_end]').fill('850')
        with page.expect_response('**/api/primer-designs') as response:
            form.get_by_role('button',name='Generate',exact=True).click()
        result=response.value.json()
        assert response.value.status==201 and result['pairs']
        assert result['parameters']['target_start']==650 and result['parameters']['target_end']==850
        assert all(pair['left']['binding_end']<=650 and pair['right']['binding_start']>=850 for pair in result['pairs'])
        expect(page.locator('#primer-status')).to_contain_text('Generated')
        expect(page.locator('#primer-design-summary')).to_contain_text('Flanking target 651–850 bp')
        expect(page.locator('#primer-design-summary')).to_contain_text('Specificity not evaluated')
        expect(page.locator('#primer-results .primer')).to_have_count(len(result['primers']))
        expect(page.locator('#primer-review')).to_be_visible()
        expect(page.locator('#primer-review')).to_contain_text('1 exact directional match')
        expect(page.locator('.pair-product-check')).to_contain_text('No additional product interval found in this range')
        assert all(p['product_review']['alternative_product_count']==0 for p in result['pairs'])
        choice=page.locator('#primer-pair-choice')
        expect(choice.locator('option')).to_have_count(len(result['pairs']))
        choice.select_option(index=1)
        expect(page.locator('.primer-product-scroll svg')).to_have_attribute('aria-label',f"Primer pair 2; forward arrow points right, reverse arrow points left; {result['pairs'][1]['product_size']} bp product")
        page.locator('#primer-review details').first.locator('summary').click()
        expect(page.locator('#primer-review details').first).to_contain_text('intended')
        page.locator('#primer-review').screenshot(path=str(artifacts/'pair-desktop.png'))
        page.locator('#primer-results').get_by_role('button',name='Select',exact=True).first.click()
        expect(page.locator('#primer-status')).to_contain_text('selection saved')
        with page.expect_download() as download:
            page.get_by_role('button',name='Export selected primers').click()
        download.value.save_as(str(artifacts/'selected.csv'))
        assert 'TargetPCR' in (artifacts/'selected.csv').read_text()
        page.locator('#map-panel').evaluate('(el)=>el.scrollTop=0')
        page.screenshot(path=str(artifacts/'target-desktop.png'))
        page.reload();page.get_by_role('button',name='Open map').click()
        expect(page.locator('#primer-results')).to_contain_text('Target 651–850 bp')
        form.locator('[name=target_start]').fill('1');form.locator('[name=target_end]').fill('10')
        form.get_by_role('button',name='Generate',exact=True).click()
        expect(page.locator('#primer-status')).to_contain_text('at least 18 bases')
        expect(form.get_by_role('button',name='Generate',exact=True)).to_be_enabled()
        form.locator('[name=target_start]').fill('900');form.locator('[name=target_end]').fill('850')
        form.get_by_role('button',name='Generate',exact=True).click()
        expect(page.locator('#primer-status')).to_contain_text('must not exceed')
        page.set_viewport_size({'width':390,'height':844})
        form.evaluate('(el)=>el.scrollIntoView({block:"start"})')
        page.evaluate('new Promise(r=>requestAnimationFrame(()=>requestAnimationFrame(r)))')
        page.screenshot(path=str(artifacts/'target-mobile.png'))
        assert page.locator('#map-panel').evaluate('(el)=>el.scrollWidth<=el.clientWidth+1')
        page.locator('#primer-review').screenshot(path=str(artifacts/'pair-mobile.png'))
        scroll=page.locator('.primer-product-scroll')
        scroll.focus();scroll.press('End')
        scroll.evaluate('(el)=>el.scrollLeft=el.scrollWidth')
        assert scroll.evaluate('(el)=>el.scrollLeft>0')
        page.locator('#primer-review').screenshot(path=str(artifacts/'pair-mobile-right.png'))
        assert page.locator('#primer-review').evaluate('(el)=>el.scrollWidth<=el.clientWidth+1')
        # Duplicate synthetic template: each designed pair has a second exact product.
        repeated=SeqRecord(Seq(sequence*2),id='repeat-demo',name='repeat-demo',annotations={'molecule_type':'DNA','topology':'linear'})
        repeated_stream=io.StringIO();SeqIO.write(repeated,repeated_stream,'genbank')
        repeated_zip=io.BytesIO()
        with zipfile.ZipFile(repeated_zip,'w') as z:z.writestr('repeat-demo.gb',repeated_stream.getvalue())
        imported=httpx.post(base+'/api/imports/benchling',files={'file':('repeat.zip',repeated_zip.getvalue())})
        assert imported.status_code==201,imported.text
        page.reload()
        page.locator('.seq').filter(has_text='repeat-demo').get_by_role('button',name='Open map').click()
        expect(form.locator('[name=target_start]')).to_be_enabled()
        form.locator('[name=product_size_min]').fill('300');form.locator('[name=product_size_max]').fill('500')
        form.locator('[name=target_start]').fill('651');form.locator('[name=target_end]').fill('850')
        with page.expect_response('**/api/primer-designs') as duplicate_response:
            form.get_by_role('button',name='Generate',exact=True).click()
        duplicate=duplicate_response.value.json()
        assert duplicate_response.value.status==201 and duplicate['pairs']
        assert all(p['product_review']['total_product_count']==2 and p['product_review']['alternative_product_count']==1 for p in duplicate['pairs'])
        expect(page.locator('.pair-product-check')).to_contain_text('1 alternative product interval found')
        page.locator('.pair-product-check summary').click()
        expect(page.locator('.pair-products-scroll tbody tr')).to_have_count(2)
        expect(page.locator('.pair-products-scroll tbody')).to_contain_text('Alternative')
        page.set_viewport_size({'width':1440,'height':1080})
        page.locator('.pair-product-check').screenshot(path=str(artifacts/'alternative-products-desktop.png'))
        page.set_viewport_size({'width':390,'height':844})
        page.locator('.pair-product-check').evaluate('(el)=>el.scrollIntoView({block:"start"})')
        page.evaluate('new Promise(r=>requestAnimationFrame(()=>requestAnimationFrame(r)))')
        page.screenshot(path=str(artifacts/'alternative-products-mobile.png'))
        table_scroll=page.locator('.pair-products-scroll')
        table_scroll.evaluate('(el)=>el.scrollLeft=el.scrollWidth')
        assert table_scroll.evaluate('(el)=>el.scrollLeft>0')
        assert page.locator('#map-panel').evaluate('(el)=>el.scrollWidth<=el.clientWidth+1')
        page.reload();page.locator('.seq').filter(has_text='repeat-demo').get_by_role('button',name='Open map').click()
        expect(page.locator('.pair-product-check')).to_contain_text('1 alternative product interval found')
        assert not errors,errors
        browser.close()
    print(json.dumps({'passed':True,'artifacts':str(artifacts),'checks':['real Primer3 target flanking','duplicated reference yields additional exact product','persisted product review and narrow table scrolling','exact-site evidence and pair switching','mobile diagram horizontal scroll','one-based UI conversion','selection and CSV','persisted target provenance','invalid target feedback','390px layout','no JS exceptions']},indent=2))
finally:
    server.terminate();server.wait(timeout=10);log.close()
