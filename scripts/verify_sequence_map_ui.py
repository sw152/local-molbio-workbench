"""Synthetic circular, linear and crowded-map browser regression. Outputs stay in var."""
from pathlib import Path
import io
import json
import os
import socket
import subprocess
import sys
import time
import zipfile
import httpx
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation, CompoundLocation
from playwright.sync_api import sync_playwright, expect

ROOT=Path(__file__).resolve().parents[1]
artifacts=ROOT/'var'/'map-ui-check'/time.strftime('%Y%m%d-%H%M%S')
artifacts.mkdir(parents=True)
def feature(start,end,strand,label,kind='CDS'):
    return SeqFeature(FeatureLocation(start,end,strand=strand),type=kind,qualifiers={'label':[label]})
records=[]
for name,topology in [('pCircular-demo','circular'),('linear-demo','linear'),('crowded-demo','circular'),('unknown-demo',None)]:
    record=SeqRecord(Seq('ACGT'*300),id=name,name=name,description='Synthetic map regression',annotations={'molecule_type':'DNA'})
    if topology:record.annotations['topology']=topology
    record.features=[feature(0,1200,None,'Source','source'),feature(80,420,1,'Reporter'),feature(300,650,-1,'Reverse marker'),feature(700,810,1,'Promoter','promoter'),feature(850,1100,None,'Origin','rep_origin')]
    if topology=='circular':
        record.features.append(SeqFeature(CompoundLocation([FeatureLocation(1100,1200,strand=1),FeatureLocation(0,70,strand=1)]),type='CDS',qualifiers={'label':['Across origin']}))
    if name=='crowded-demo':
        record.features += [feature(150+i,1000-i,1 if i%2 else -1,f'Overlapping feature {i+1}') for i in range(65)]
    records.append(record)
archive=artifacts/'maps.zip'
with zipfile.ZipFile(archive,'w') as z:
    for record in records:
        stream=io.StringIO();SeqIO.write(record,stream,'genbank');z.writestr(record.name+'.gb',stream.getvalue())
with socket.socket() as socket_:
    socket_.bind(('127.0.0.1',0));port=socket_.getsockname()[1]
base=f'http://127.0.0.1:{port}'
log=(artifacts/'server.log').open('w')
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
        def snapshot(name):
            page.evaluate('new Promise(resolve=>requestAnimationFrame(()=>requestAnimationFrame(resolve)))')
            page.screenshot(path=str(artifacts/name))
        page.goto(base)
        page.locator('#archive').set_input_files(str(archive))
        page.get_by_role('button',name='Import into local library').click()
        expect(page.locator('#status')).to_contain_text('Imported 4 records')
        def open_map(name):
            page.locator('.seq').filter(has_text=name).get_by_role('button',name='Open map').click()
            expect(page.locator('#map-title')).to_have_text(name)
            expect(page.locator('#map-layout-note')).to_contain_text('annotation tracks')
        open_map('pCircular-demo')
        expect(page.locator('#plasmid-map .map-feature')).to_have_count(6)
        page.locator('#feature-legend button').filter(has_text='Across origin').click()
        expect(page.locator('#map-selection')).to_contain_text('1101–1200 bp')
        expect(page.locator('#map-selection')).to_contain_text('1–70 bp')
        expect(page.locator('#plasmid-map .map-selected')).to_have_count(1)
        snapshot('circular-desktop.png')
        page.get_by_role('button',name='Zoom in',exact=True).click()
        expect(page.locator('#map-zoom-value')).to_have_text('150%')
        assert page.locator('.map-viewport').evaluate('(el)=>el.scrollWidth>el.clientWidth')
        page.get_by_role('button',name='Fit',exact=True).click()
        page.locator('#close-map').click()
        open_map('linear-demo')
        expect(page.locator('#plasmid-map circle')).to_have_count(0)
        expect(page.locator('#plasmid-map polygon')).to_have_count(5)
        target=page.locator('#plasmid-map [data-feature="2"]')
        target.focus();page.keyboard.press('Enter')
        expect(page.locator('#map-selection')).to_contain_text('strand −')
        snapshot('linear-desktop.png')
        page.set_viewport_size({'width':390,'height':844})
        page.locator('#map-panel').evaluate('(el)=>el.scrollTop=0')
        snapshot('linear-mobile.png')
        assert page.locator('#map-panel').evaluate('(el)=>el.scrollWidth<=el.clientWidth+1')
        page.locator('#close-map').click()
        open_map('crowded-demo')
        expect(page.locator('#feature-legend button')).to_have_count(71)
        expect(page.locator('#plasmid-map .map-feature')).to_have_count(71)
        page.locator('#feature-legend button').filter(has_text='Overlapping feature 65').click()
        expect(page.locator('#map-selection')).to_contain_text('Overlapping feature 65')
        page.locator('#map-panel').evaluate('(el)=>el.scrollTop=0')
        snapshot('crowded-mobile.png')
        page.set_viewport_size({'width':1440,'height':1080})
        snapshot('crowded-desktop.png')
        page.get_by_role('button',name='Focus selected',exact=True).click()
        expect(page.locator('#plasmid-map .map-feature')).to_have_count(1)
        expect(page.locator('#feature-legend button')).to_have_count(71)
        expect(page.locator('#map-layout-note')).to_contain_text('Showing 1 of 71')
        snapshot('crowded-focused-desktop.png')
        page.get_by_role('button',name='Show all',exact=True).click()
        expect(page.locator('#plasmid-map .map-feature')).to_have_count(71)
        page.locator('#close-map').click()
        open_map('unknown-demo')
        expect(page.locator('#plasmid-map')).to_contain_text('TOPOLOGY UNKNOWN')
        page.locator('#close-map').click()
        page.route('**/api/sequences/*/map',lambda route:route.fulfill(status=503,content_type='application/json',body='{"detail":"Map temporarily unavailable"}'))
        page.locator('.seq').filter(has_text='linear-demo').get_by_role('button',name='Open map').click()
        expect(page.locator('#map-message')).to_contain_text('Map temporarily unavailable')
        expect(page.locator('#plasmid-map .map-feature')).to_have_count(0)
        expect(page.locator('#map-focus')).to_be_disabled()
        expect(page.locator('#map-zoom-in')).to_be_disabled()
        assert not errors,errors
        browser.close()
    print(json.dumps({'passed':True,'artifacts':str(artifacts),'checks':['circular origin feature','linear strands','unknown topology','71 annotations','keyboard selection','zoom and fit','390px viewport','no JS errors']},indent=2))
finally:
    server.terminate();server.wait(timeout=10);log.close()
