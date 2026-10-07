"""Optional end-to-end browser check using synthetic GenBank and binary ABIF.

Run with the ui extra installed, e.g. python scripts/verify_sanger_ui.py.
Artifacts and its isolated database stay under var/ui-check (never committed).
Uses an installed Chrome or Playwright Chromium; no private browser profile.
"""
from pathlib import Path
import io
import json
import os
import random
import socket
import subprocess
import sys
import time
import zipfile

import httpx
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from playwright.sync_api import sync_playwright, expect

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "tests"))
from abif_fixture import synthetic_ab1

artifacts = ROOT / "var" / "ui-check" / time.strftime("%Y%m%d-%H%M%S")
artifacts.mkdir(parents=True)
rng = random.Random(43)
reference = "".join(rng.choice("ACGT") for _ in range(1400))
record = SeqRecord(Seq(reference), id="synthetic", name="pEvidence-demo", description="Synthetic UI verification construct", annotations={"molecule_type":"DNA", "topology":"circular"})
record.features = [SeqFeature(FeatureLocation(start, end, strand=1 if i % 2 else -1), type=kind, qualifiers={"label":[label]}) for i, (start, end, kind, label) in enumerate([(60,340,"CDS","Reporter"),(380,520,"promoter","Promoter"),(650,880,"rep_origin","Origin"),(940,1250,"CDS","Marker")])]
handle = io.StringIO()
SeqIO.write(record, handle, "genbank")
archive = artifacts / "synthetic.zip"
with zipfile.ZipFile(archive, "w") as z:
    z.writestr("synthetic.gb", handle.getvalue())
read = list(reference[150:750])
read[200] = next(base for base in "ACGT" if base != read[200])
quality = [38] * len(read)
quality[200] = 12
read_path = artifacts / "synthetic-Q12.ab1"
read_path.write_bytes(synthetic_ab1("".join(read), quality, trace=True))
trim_calls=list(reference[300:900])
trim_calls[180]=next(b for b in 'ACGT' if b!=trim_calls[180])
trim_quality=[5]*20+[38]*550+[8]*30;trim_quality[180]=12
trim_path=artifacts/'synthetic-reverse-trim.ab1'
trim_path.write_bytes(synthetic_ab1(str(Seq(''.join(trim_calls)).reverse_complement()),trim_quality[::-1],trace=True))
invalid = artifacts / "invalid.ab1"
invalid.write_bytes(b"invalid data")
with socket.socket() as s:
    s.bind(("127.0.0.1",0))
    port = s.getsockname()[1]
base = f"http://127.0.0.1:{port}"
environment = dict(os.environ, MOLBIO_DATA_DIR=str(artifacts / "runtime"))
log = (artifacts / "server.log").open("w")
server = subprocess.Popen([sys.executable,"-m","uvicorn","localmolbio.api:app","--host","127.0.0.1","--port",str(port)], cwd=ROOT, env=environment, stdout=log, stderr=log)
try:
    for _ in range(100):
        try:
            if httpx.get(base + "/health").status_code == 200:
                break
        except httpx.TransportError:
            pass
        time.sleep(.1)
    else:
        raise RuntimeError("Test server did not start")
    with sync_playwright() as p:
        chrome = Path("/Applications/Google Chrome.app/Contents/MacOS/Google Chrome")
        browser = p.chromium.launch(headless=True, **({"executable_path":str(chrome)} if chrome.exists() else {}))
        page = browser.new_page(viewport={"width":1440,"height":1080}, device_scale_factor=1)
        errors=[]
        page.on("pageerror", lambda error: errors.append(str(error)))
        page.goto(base)
        page.locator('#archive').set_input_files(str(archive))
        page.get_by_role('button',name='Import into local library').click()
        expect(page.locator('#status')).to_contain_text('Imported 1 records')
        page.get_by_role('button',name='Open map').click()
        expect(page.locator('#sanger-reads')).to_contain_text('No reads attached')
        expect(page.locator('#sanger-summary')).to_contain_text('0 / 1,400 bp aligned')
        expect(page.locator('#sanger-summary')).to_contain_text('No eligible saved reads')
        page.locator('#sanger-file').set_input_files(str(invalid))
        page.get_by_role('button',name='Attach AB1',exact=True).click()
        expect(page.locator('#sanger-status')).to_contain_text('Could not attach read')
        assert page.locator('#sanger-upload-button').is_enabled()
        page.locator('#sanger-file').set_input_files(str(read_path))
        page.get_by_role('button',name='Attach AB1',exact=True).click()
        expect(page.locator('#sanger-reads')).to_contain_text('Ready to analyze')
        page.locator('.summary-sources summary').click()
        expect(page.locator('.summary-sources')).to_contain_text('No saved analysis')
        page.get_by_role('button',name='Analyze against this revision').click()
        expect(page.locator('#sanger-status')).to_contain_text('Analysis saved')
        expect(page.locator('#sanger-reads')).to_contain_text('99.8%')
        expect(page.locator('#sanger-reads')).to_contain_text('42.9%')
        expect(page.locator('#sanger-reads')).to_contain_text('below Q20')
        page.locator('.variant-details summary').click()
        expect(page.locator('.variant-scroll tbody tr')).to_have_count(1)
        expect(page.locator('.variant-scroll tbody')).to_contain_text('351')
        expect(page.locator('.variant-scroll tbody')).to_contain_text('12')
        page.locator('.trace-jump').click()
        expect(page.locator('.trace-heading')).to_contain_text('original read')
        expect(page.locator('.trace-view')).to_contain_text('Bases 189–212 of 600')
        expect(page.locator('.trace-svg')).to_have_count(1)
        expect(page.locator('.trace-svg [data-channel]')).to_have_count(4)
        expect(page.locator('.trace-selected')).to_have_attribute('data-peak-position','200')
        page.locator('.trace-view').screenshot(path=str(artifacts/'trace-desktop.png'))
        page.get_by_role('button',name='Next →',exact=True).click()
        expect(page.locator('.trace-view')).to_contain_text('Bases 213–236 of 600')
        page.locator('.trace-position').fill('201')
        page.get_by_role('button',name='Go',exact=True).click()
        expect(page.locator('.trace-view')).to_contain_text('Bases 201–224 of 600')
        page.locator('#sanger-panel').screenshot(path=str(artifacts/'sanger-desktop.png'))
        page.locator('#map-panel').evaluate('(el) => el.scrollTop = 0')
        page.screenshot(path=str(artifacts/'construct-desktop.png'))
        # Persistence through reload and duplicate-upload feedback.
        page.reload()
        page.get_by_role('button',name='Open map').click()
        expect(page.locator('#sanger-reads')).to_contain_text('99.8%')
        page.locator('#sanger-file').set_input_files(str(read_path))
        page.get_by_role('button',name='Attach AB1',exact=True).click()
        expect(page.locator('#sanger-status')).to_contain_text('already attached')
        # A rerun creates an independent report, with the first still selectable.
        page.locator('[data-analysis-direction]').select_option('forward')
        page.get_by_role('button',name='Run new analysis',exact=True).click()
        expect(page.locator('.analysis-version')).to_contain_text('Run 2 · Latest')
        expect(page.locator('.analysis-version')).to_contain_text('requested forward')
        page.route('**/api/sanger-reads/*/analyses?*',lambda route:route.abort(),times=1)
        page.get_by_role('button',name='Browse saved analyses').click()
        expect(page.locator('.analysis-history')).to_contain_text('Could not load history')
        page.get_by_role('button',name='Retry history').click()
        expect(page.locator('.history-runs button')).to_have_count(2)
        page.locator('.history-runs button').filter(has_text='Run 1').click()
        expect(page.locator('.analysis-version')).to_contain_text('Run 1 · Historical')
        expect(page.locator('.analysis-version')).to_contain_text('requested unknown')
        with page.expect_download() as old_download:
            page.get_by_role('button',name='Download JSON',exact=True).click()
        old_download.value.save_as(str(artifacts/'historical-run-1.json'))
        old_payload=json.loads((artifacts/'historical-run-1.json').read_text())
        assert old_payload['analysis']['run_number']==1
        assert old_payload['saved_alignment']['evidence']['requested_direction']=='unknown'
        page.locator('#sanger-panel').screenshot(path=str(artifacts/'history-desktop.png'))
        page.get_by_role('button',name='Show latest report').click()
        expect(page.locator('.analysis-version')).to_contain_text('Run 2 · Latest')
        page.route('**/api/sanger-verifications',lambda route:route.fulfill(status=409,content_type='application/json',body=json.dumps({'detail':'Analysis history changed; reload before starting a new analysis'})),times=1)
        page.get_by_role('button',name='Run new analysis',exact=True).click()
        expect(page.locator('#sanger-status')).to_contain_text('history changed')
        expect(page.get_by_role('button',name='Run new analysis',exact=True)).to_be_enabled()
        expect(page.locator('.analysis-version')).to_contain_text('Run 2 · Latest')
        page.reload();page.get_by_role('button',name='Open map').click()
        expect(page.locator('.analysis-version')).to_contain_text('Run 2 · Latest')
        page.set_viewport_size({'width':390,'height':844})
        page.locator('.variant-details summary').click()
        page.locator('.trace-jump').click()
        expect(page.locator('.trace-view')).to_contain_text('Bases 189–212 of 600')
        page.locator('.trace-view').evaluate('(el) => el.scrollIntoView({block: "start"})')
        assert page.locator('.trace-selected').evaluate('(el) => {const r=el.getBoundingClientRect();const s=el.closest(".trace-scroll").getBoundingClientRect();return r.left>=s.left && r.right<=s.right;}')
        page.screenshot(path=str(artifacts/'trace-mobile.png'))
        page.locator('#sanger-panel').evaluate('(el) => el.scrollIntoView({block: "start"})')
        page.screenshot(path=str(artifacts/'sanger-mobile.png'))
        page.locator('.variant-details').evaluate('(el) => el.scrollIntoView({block: "center"})')
        page.screenshot(path=str(artifacts/'sanger-mobile-details.png'))
        assert page.locator('#sanger-panel').evaluate('(el) => el.scrollWidth <= el.clientWidth + 1')
        assert page.locator('#map-panel').evaluate('(el) => el.scrollWidth <= el.clientWidth + 1'), page.locator('#map-panel').evaluate('(el) => Array.from(el.querySelectorAll("*")).filter(e => e.getBoundingClientRect().right > innerWidth).map(e => [e.tagName,e.className,e.getBoundingClientRect().width]).slice(0,25)')
        # More runs exercise real history pagination without duplicating read cards.
        read_id=page.locator('.read-card').get_attribute('data-read-id')
        previous=httpx.get(base+f'/api/sanger-reads/{read_id}/analyses').json()['items'][0]['alignment_id']
        for _ in range(4):
            response=httpx.post(base+'/api/sanger-verifications',json={'sequencing_read_id':read_id,'previous_alignment_id':previous})
            assert response.status_code==201,response.text
            previous=response.json()['id']
        page.reload();page.get_by_role('button',name='Open map').click()
        expect(page.locator('.read-card')).to_have_count(1)
        expect(page.locator('.analysis-version')).to_contain_text('Run 6 · Latest')
        expect(page.locator('#sanger-summary')).to_contain_text('600 / 1,400 bp aligned')
        expect(page.locator('.summary-metrics')).to_contain_text('Included · 0 excluded')
        page.get_by_role('button',name='Browse saved analyses').click()
        expect(page.locator('.history-runs button')).to_have_count(5)
        page.get_by_role('button',name='Older runs').click()
        expect(page.locator('.history-runs button')).to_have_count(1)
        page.locator('.history-runs button').click()
        expect(page.locator('.analysis-version')).to_contain_text('Run 1 · Historical')
        page.locator('.analysis-version').evaluate('(el)=>el.scrollIntoView({block:"start"})')
        page.screenshot(path=str(artifacts/'history-mobile-report.png'))
        page.locator('.analysis-history').evaluate('(el)=>el.scrollIntoView({block:"center"})')
        page.screenshot(path=str(artifacts/'history-mobile-pagination.png'))
        assert page.locator('#map-panel').evaluate('(el)=>el.scrollWidth<=el.clientWidth+1')
        page.get_by_role('button',name='Newer runs').click()
        expect(page.locator('.history-runs button')).to_have_count(5)
        # Real binary ABIF: reverse-oriented read with asymmetric low-quality ends.
        page.set_viewport_size({'width':1440,'height':1080})
        page.locator('#sanger-file').set_input_files(str(trim_path))
        page.get_by_role('button',name='Attach AB1',exact=True).click()
        trim_card=page.locator('.read-card').filter(has_text=trim_path.name)
        expect(trim_card).to_have_count(1)
        expect(trim_card.locator('[data-analysis-trim]')).to_have_value('')
        trim_card.get_by_role('button',name='Analyze against this revision').click()
        expect(trim_card.locator('.analysis-version')).to_contain_text('Run 1 · Latest')
        expect(trim_card.locator('.trim-off')).to_contain_text('all 600 original bases')
        trim_card.locator('[data-analysis-trim]').select_option('20')
        with page.expect_response('**/api/sanger-verifications') as trim_response:
            trim_card.get_by_role('button',name='Run new analysis',exact=True).click()
        trim_report=trim_response.value.json()
        assert trim_response.value.status==201
        assert trim_report['evidence']['direction']=='reverse'
        assert trim_report['evidence']['read_aligned_fraction']==.916667
        assert trim_report['evidence']['retained_read_aligned_fraction']==1
        assert trim_report['variants'][0]['original_read_position']==419
        expect(trim_card.locator('.trim-evidence')).to_contain_text('550 / 600 bases retained')
        expect(trim_card.locator('.trim-evidence')).to_contain_text('excluded 30 left / 20 right')
        expect(trim_card.locator('.read-metrics')).to_contain_text('91.7%')
        expect(trim_card.locator('.read-flags')).to_contain_text('below Q20')
        trim_card.locator('.variant-details summary').click()
        expect(trim_card.locator('.variant-scroll tbody')).to_contain_text('481')
        expect(trim_card.locator('.variant-scroll tbody')).to_contain_text('420')
        trim_card.locator('.variant-scroll .trace-jump').click()
        expect(trim_card.locator('.trace-selected')).to_have_attribute('data-peak-position','419')
        trim_card.locator('.analysis-report').screenshot(path=str(artifacts/'trim-report-desktop.png'))
        trim_card.locator('.trim-evidence .trace-jump').first.click()
        expect(trim_card.locator('.trace-selected')).to_have_attribute('data-peak-position','30')
        trim_card.get_by_role('button',name='Browse saved analyses').click()
        trim_card.locator('.history-runs button').filter(has_text='Run 1').click()
        expect(trim_card.locator('.trim-off')).to_contain_text('all 600 original bases')
        expect(trim_card.locator('.trim-evidence')).to_have_count(0)
        trim_card.get_by_role('button',name='Show latest report').click()
        expect(trim_card.locator('.trim-evidence')).to_be_visible()
        page.reload();page.get_by_role('button',name='Open map').click()
        expect(trim_card.locator('[data-analysis-trim]')).to_have_value('20')
        expect(trim_card.locator('.trim-evidence')).to_contain_text('550 / 600 bases retained')
        page.set_viewport_size({'width':390,'height':844})
        trim_card.locator('.analysis-version').evaluate('(el)=>el.scrollIntoView({block:"start"})')
        page.screenshot(path=str(artifacts/'trim-report-mobile.png'))
        trim_card.locator('.variant-details summary').click()
        trim_card.locator('.variant-scroll .trace-jump').click()
        expect(trim_card.locator('.trace-selected')).to_have_attribute('data-peak-position','419')
        trim_card.locator('.trace-view').evaluate('(el)=>el.scrollIntoView({block:"start"})')
        page.screenshot(path=str(artifacts/'trim-trace-mobile.png'))
        assert trim_card.locator('.trace-selected').evaluate('(el)=>{const r=el.getBoundingClientRect();const s=el.closest(".trace-scroll").getBoundingClientRect();return r.left>=s.left&&r.right<=s.right;}')
        assert page.locator('#map-panel').evaluate('(el)=>el.scrollWidth<=el.clientWidth+1')
        # Latest-only union: six forward runs plus two reverse runs are two reads.
        summary=page.locator('#sanger-summary')
        expect(summary).to_contain_text('720 / 1,400 bp aligned')
        expect(summary).to_contain_text('51.4%')
        expect(summary.locator('.summary-metrics')).to_contain_text('430 bp')
        expect(summary).to_contain_text('Not quality-screened')
        summary.locator('.summary-gaps summary').click()
        expect(summary.locator('.summary-gaps')).to_contain_text('871–1,400 + 1–150 bp')
        expect(summary.locator('.summary-gaps')).to_contain_text('680 bp · crosses origin')
        summary.locator('.summary-sources summary').click()
        expect(summary.locator('.summary-sources')).to_contain_text('Run 6')
        expect(summary.locator('.summary-sources')).to_contain_text('Run 2')
        expect(summary.locator('.summary-sources')).to_contain_text('Includes calls below Q20')
        assert summary.evaluate('(el)=>el.scrollWidth<=el.clientWidth+1')
        summary.locator('.combined-coverage').evaluate('(el)=>el.scrollIntoView({block:"start"})')
        page.screenshot(path=str(artifacts/'combined-coverage-mobile.png'))
        summary.locator('.summary-sources').evaluate('(el)=>el.scrollIntoView({block:"start"})')
        page.screenshot(path=str(artifacts/'combined-sources-mobile.png'))
        page.set_viewport_size({'width':1440,'height':1080})
        summary.locator('.combined-coverage').screenshot(path=str(artifacts/'combined-coverage-desktop.png'))
        trim_card.get_by_role('button',name='Browse saved analyses').click()
        trim_card.locator('.history-runs button').filter(has_text='Run 1').click()
        expect(trim_card.locator('.analysis-version')).to_contain_text('Run 1 · Historical')
        expect(summary).to_contain_text('720 / 1,400 bp aligned')
        trim_card.get_by_role('button',name='Show latest report').click()
        page.set_viewport_size({'width':390,'height':844})
        # Downloads stay bound to the displayed analysis, with recoverable errors.
        page.route('**/analyses/*/export?format=html',lambda route:route.fulfill(status=503,content_type='application/json',body=json.dumps({'detail':'Simulated report export outage'})),times=1)
        trim_card.get_by_role('button',name='Download HTML',exact=True).click()
        expect(page.locator('#sanger-status')).to_contain_text('Could not export report')
        expect(trim_card.get_by_role('button',name='Download HTML',exact=True)).to_be_enabled()
        with page.expect_download() as html_download:
            trim_card.get_by_role('button',name='Download HTML',exact=True).click()
        html_path=artifacts/'trimmed-run-2.html';html_download.value.save_as(str(html_path))
        with page.expect_download() as json_download:
            trim_card.get_by_role('button',name='Download JSON',exact=True).click()
        json_download.value.save_as(str(artifacts/'trimmed-run-2.json'))
        payload=json.loads((artifacts/'trimmed-run-2.json').read_text())
        assert payload['analysis']['run_number']==2
        assert payload['parameters']['trim_quality_threshold']==20
        assert payload['saved_alignment']['variants'][0]['original_read_position']==419
        assert payload['inputs']['read']['filename']['value']==trim_path.name
        # Open the downloaded HTML offline: no dependency on a running app or network.
        offline=browser.new_context(offline=True,viewport={'width':1100,'height':1000})
        report_page=offline.new_page();report_errors=[];report_requests=[]
        report_page.on('pageerror',lambda error:report_errors.append(str(error)))
        report_page.on('request',lambda request:report_requests.append(request.url))
        report_page.goto(html_path.as_uri())
        expect(report_page.locator('h1')).to_have_text('Sanger analysis report')
        expect(report_page.locator('header')).to_contain_text('Run 2')
        expect(report_page.locator('.metrics')).to_contain_text('91.7%')
        expect(report_page.locator('.variants tbody')).to_contain_text('420')
        report_page.screenshot(path=str(artifacts/'export-report-desktop.png'),full_page=True)
        report_page.get_by_role('heading',name='Differences',exact=False).scroll_into_view_if_needed()
        report_page.screenshot(path=str(artifacts/'export-report-differences.png'))
        report_page.set_viewport_size({'width':390,'height':844})
        report_page.evaluate('window.scrollTo(0,0)')
        report_page.evaluate('new Promise(r=>requestAnimationFrame(()=>requestAnimationFrame(r)))')
        report_page.screenshot(path=str(artifacts/'export-report-mobile.png'))
        assert report_page.evaluate('document.documentElement.scrollWidth<=innerWidth+1')
        assert not report_errors,report_errors
        assert not any(url.startswith(('http:','https:')) for url in report_requests),report_requests
        offline.close()
        # Closing the construct cancels delivery of a late export response.
        pending=[];late_downloads=[]
        page.route('**/analyses/*/export?format=json',lambda route:pending.append(route),times=1)
        page.on('download',lambda download:late_downloads.append(download))
        trim_card.get_by_role('button',name='Download JSON',exact=True).click()
        for _ in range(30):
            if pending:break
            page.wait_for_timeout(20)
        assert pending
        page.locator('#close-map').click()
        with page.expect_response('**/analyses/*/export?format=json'):
            pending[0].fulfill(status=200,content_type='application/json',body=json.dumps(payload))
        page.evaluate('new Promise(r=>requestAnimationFrame(()=>requestAnimationFrame(r)))')
        assert not late_downloads
        page.get_by_role('button',name='Open map').click()
        expect(trim_card.locator('.analysis-version')).to_contain_text('Run 2 · Latest')
        # Optional externally supplied public ABI fixture; never bundled or committed.
        public_fixture = os.environ.get('MOLBIO_PUBLIC_AB1')
        if public_fixture:
            public_path = Path(public_fixture).resolve()
            page.set_viewport_size({'width':1440,'height':1080})
            page.locator('#sanger-file').set_input_files(str(public_path))
            page.get_by_role('button',name='Attach AB1',exact=True).click()
            card = page.locator('.read-card').filter(has_text=public_path.name)
            expect(card).to_have_count(1)
            card.get_by_role('button',name='View chromatogram',exact=True).click()
            expect(card.locator('.trace-svg')).to_have_count(1)
            card.locator('.trace-position').fill('101')
            card.get_by_role('button',name='Go',exact=True).click()
            expect(card.locator('.trace-view')).to_contain_text('Bases 101–124')
            card.locator('.trace-view').screenshot(path=str(artifacts/'public-trace-desktop.png'))
        # A failed list request must be shown, not reported as an empty library.
        page.locator('#close-map').click()
        page.route('**/api/sequence-revisions/*/sanger-reads*', lambda route: route.abort())
        page.get_by_role('button',name='Open map').click()
        expect(page.locator('#sanger-status')).to_contain_text('Could not load reads')
        expect(page.locator('#sanger-summary')).to_be_empty()
        assert not errors, errors
        browser.close()
    print(json.dumps({'passed':True,'artifacts':str(artifacts),'checks':['latest-only multi-read coverage, overlap and circular gap review','binary ABIF upload','invalid and duplicate upload','analysis and persisted evidence','append-only rerun and historical report switching','six-run history pagination','optional end trimming with original and retained coverage','historical and current JSON/HTML downloads','offline HTML rendering without external requests','export failure recovery and late-response cancellation','reverse trimmed variant and boundary jump to original peaks','history network failure and stale rerun recovery','Q12 variant and focused chromatogram','desktop and mobile layouts','network failure','no JS exceptions']},indent=2))
finally:
    server.terminate()
    server.wait(timeout=10)
    log.close()
