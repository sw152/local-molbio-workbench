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
high_calls=list(reference[100:500]);high_calls[90]=next(b for b in 'ACGT' if b!=high_calls[90])
high_path=artifacts/'synthetic-Q40-difference.ab1'
high_path.write_bytes(synthetic_ab1(''.join(high_calls),[40]*400,trace=True))
inserted=next(b for b in 'ACGT' if b!=reference[650])+'G'+next(b for b in 'ACGT' if b!=reference[649])
indel_calls=reference[500:650]+inserted+reference[650:800]
indel_path=artifacts/'synthetic-reverse-insertion.ab1'
indel_path.write_bytes(synthetic_ab1(str(Seq(indel_calls).reverse_complement()),[35]*len(indel_calls),trace=True))
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
        expect(page.locator('.quality-assessment')).to_contain_text('Not evaluated')
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
        expect(page.locator('.combined-coverage .summary-metrics')).to_contain_text('Included · 0 excluded')
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
        expect(trim_card.locator('.base-mapping-counts')).to_contain_text('549')
        expect(trim_card.locator('.base-mapping-counts')).to_contain_text('0Differing calls')
        trim_card.locator('.base-mapping summary').click()
        expect(trim_card.locator('.base-mapping')).to_contain_text('580 → 31')
        expect(trim_card.locator('.base-mapping')).to_contain_text('1 low-quality')
        trim_card.locator('.base-mapping').screenshot(path=str(artifacts/'base-mapping-desktop.png'))
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
        trim_card.locator('.base-mapping summary').click()
        trim_card.locator('.base-mapping').evaluate('(el)=>el.scrollIntoView({block:"start"})')
        page.screenshot(path=str(artifacts/'base-mapping-mobile.png'))
        trim_card.locator('.mapping-scroll').evaluate('(el)=>{el.scrollLeft=el.scrollWidth;}')
        assert trim_card.locator('.mapping-scroll td').last.evaluate('(el)=>{const r=el.getBoundingClientRect();const s=el.closest(".mapping-scroll").getBoundingClientRect();return r.left>=s.left&&r.right<=s.right+1;}')
        page.screenshot(path=str(artifacts/'base-mapping-mobile-direction.png'))
        # Latest-only union: six forward runs plus two reverse runs are two reads.
        summary=page.locator('#sanger-summary')
        expect(summary).to_contain_text('720 / 1,400 bp aligned')
        expect(summary).to_contain_text('51.4%')
        expect(summary.locator('.combined-coverage .summary-metrics')).to_contain_text('430 bp')
        expect(summary).to_contain_text('Not quality-screened')
        expect(page.locator('.quality-review')).to_contain_text('720 / 1,400 bp with Q20+ calls')
        expect(page.locator('.quality-review .summary-metrics')).to_contain_text('428 bp')
        expect(page.locator('.quality-conflict-count')).to_have_text('0')
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
        # High quality does not mean reference matching: exercise amber evidence with a real AB1.
        page.set_viewport_size({'width':1440,'height':1080})
        page.locator('#sanger-file').set_input_files(str(high_path))
        page.get_by_role('button',name='Attach AB1',exact=True).click()
        high_card=page.locator('.read-card').filter(has_text=high_path.name)
        high_card.get_by_role('button',name='Analyze against this revision').click()
        expect(high_card.locator('.base-mapping-counts')).to_contain_text('399Matching calls')
        expect(high_card.locator('.base-mapping-counts')).to_contain_text('1Differing calls')
        expect(high_card.locator('.base-mapping svg')).to_have_attribute('aria-label','400 reference positions with unambiguous calls at Q20 or above, including 1 difference')
        high_card.locator('.base-mapping').screenshot(path=str(artifacts/'base-mapping-difference.png'))
        quality_review=page.locator('.quality-review')
        expect(quality_review).to_contain_text('770 / 1,400 bp with Q20+ calls')
        expect(quality_review.locator('.quality-conflict-count')).to_have_text('1')
        quality_review.locator('.quality-differences summary').click()
        expect(quality_review.locator('.quality-locus')).to_have_count(1)
        expect(quality_review.locator('.quality-locus')).to_contain_text('Reference 191')
        expect(quality_review.locator('.quality-locus')).to_contain_text('Conflicting calls')
        expect(quality_review.locator('.quality-locus li')).to_have_count(2)
        quality_review.screenshot(path=str(artifacts/'quality-conflict-desktop.png'))
        quality_review.locator('.quality-locus li').filter(has_text=high_path.name).get_by_role('button',name='Original base 91').click()
        expect(high_card.locator('.trace-selected')).to_have_attribute('data-peak-position','90')
        page.set_viewport_size({'width':390,'height':844})
        quality_review.evaluate('(el)=>el.scrollIntoView({block:"start"})')
        page.screenshot(path=str(artifacts/'quality-conflict-mobile.png'))
        quality_review.locator('.quality-locus').evaluate('(el)=>el.scrollIntoView({block:"start"})')
        page.screenshot(path=str(artifacts/'quality-conflict-mobile-calls.png'))
        assert quality_review.evaluate('(el)=>el.scrollWidth<=el.clientWidth+1')
        page.set_viewport_size({'width':1440,'height':1080})

        # Simulate a v3 historical response: the UI must not invent missing base-level evidence.
        def legacy_history(route):
            response=route.fetch();data=response.json()
            for item in data['items']:
                item['evidence'].pop('base_mapping',None)
                item['evidence']['parameters']['evidence_version']=3
            route.fulfill(response=response,json=data)
        page.route('**/api/sanger-reads/*/analyses?limit=*',legacy_history,times=1)
        high_card.get_by_role('button',name='Browse saved analyses').click()
        high_card.locator('.history-runs button').click()
        expect(high_card.locator('.mapping-unavailable')).to_contain_text('was not recorded')
        expect(high_card.locator('.base-mapping')).to_have_count(0)
        high_card.get_by_role('button',name='Show latest report').click()
        expect(high_card.locator('.base-mapping')).to_have_count(1)
        # Group exports bind the displayed view and preserve all latest source runs.
        page.route('**/sanger-report?*',lambda route:route.fulfill(status=503,content_type='application/json',body=json.dumps({'detail':'Simulated group export outage'})),times=1)
        page.get_by_role('button',name='Group HTML',exact=True).click()
        expect(page.locator('#sanger-status')).to_contain_text('Could not export group')
        expect(page.get_by_role('button',name='Group HTML',exact=True)).to_be_enabled()
        with page.expect_download() as group_json:
            page.get_by_role('button',name='Group JSON',exact=True).click()
        group_json_path=artifacts/'sanger-group.json';group_json.value.save_as(str(group_json_path))
        group_payload=json.loads(group_json_path.read_text())
        assert group_payload['schema']=='localmolbio.sanger-group-report'
        group_snapshot=group_payload['snapshot']
        assert group_snapshot['summary']['quality_review']['conflict_position_count']==1
        assert group_snapshot['summary']['quality_review']['callable_bases']==770
        assert next(r for r in group_snapshot['reads'] if r['original_filename']==read_path.name)['run_number']==6
        with page.expect_download() as group_html:
            page.get_by_role('button',name='Group HTML',exact=True).click()
        group_html_path=artifacts/'sanger-group.html';group_html.value.save_as(str(group_html_path))
        group_offline=browser.new_context(offline=True,viewport={'width':1200,'height':1000})
        group_page=group_offline.new_page();group_requests=[];group_errors=[]
        group_page.on('request',lambda request:group_requests.append(request.url))
        group_page.on('pageerror',lambda error:group_errors.append(str(error)))
        group_page.goto(group_html_path.as_uri())
        expect(group_page.locator('h1')).to_have_text('Sanger group evidence')
        expect(group_page.locator('header')).to_contain_text(group_payload['snapshot_sha256'])
        expect(group_page.locator('.locus')).to_contain_text('Reference 191')
        expect(group_page.locator('.locus')).to_contain_text('Q40')
        expect(group_page.locator('.sources')).to_contain_text(read_path.name)
        group_page.screenshot(path=str(artifacts/'group-export-desktop.png'))
        group_page.locator('.locus').scroll_into_view_if_needed()
        group_page.screenshot(path=str(artifacts/'group-export-conflict.png'))
        group_page.set_viewport_size({'width':390,'height':844})
        group_page.evaluate('window.scrollTo(0,0)')
        group_page.evaluate('new Promise(r=>requestAnimationFrame(()=>requestAnimationFrame(r)))')
        group_page.screenshot(path=str(artifacts/'group-export-mobile.png'))
        group_page.locator('.sources').scroll_into_view_if_needed()
        group_page.locator('.sources .scroll').evaluate('(el)=>{el.scrollLeft=el.scrollWidth;}')
        group_page.evaluate('new Promise(r=>requestAnimationFrame(()=>requestAnimationFrame(r)))')
        group_page.screenshot(path=str(artifacts/'group-export-mobile-sources.png'))
        assert group_page.evaluate('document.documentElement.scrollWidth<=innerWidth+1')
        assert not group_errors and not any(url.startswith(('http:','https:')) for url in group_requests)
        group_offline.close()
        # A new server-side run invalidates the open snapshot, not the downloaded file.
        high_saved=next(r for r in group_snapshot['reads'] if r['original_filename']==high_path.name)
        updated=httpx.post(base+'/api/sanger-verifications',json={'sequencing_read_id':high_saved['id'],'previous_alignment_id':high_saved['alignment_id']})
        assert updated.status_code==201,updated.text
        page.get_by_role('button',name='Group JSON',exact=True).click()
        expect(page.locator('#sanger-status')).to_contain_text('Group evidence changed')
        page.locator('#close-map').click();page.get_by_role('button',name='Open map').click()
        expect(high_card.locator('.analysis-version')).to_contain_text('Run 2 · Latest')
        with page.expect_download() as refreshed_group:
            page.get_by_role('button',name='Group JSON',exact=True).click()
        refreshed_path=artifacts/'sanger-group-refreshed.json';refreshed_group.value.save_as(str(refreshed_path))
        refreshed=json.loads(refreshed_path.read_text())
        assert refreshed['snapshot_sha256']!=group_payload['snapshot_sha256']
        assert next(r for r in refreshed['snapshot']['reads'] if r['id']==high_saved['id'])['run_number']==2
        # Closing the view cancels a late group download.
        group_pending=[];group_late=[]
        page.route('**/sanger-report?*',lambda route:group_pending.append(route),times=1)
        page.on('download',lambda download:group_late.append(download))
        page.get_by_role('button',name='Group JSON',exact=True).click()
        for _ in range(30):
            if group_pending:break
            page.wait_for_timeout(20)
        assert group_pending
        page.locator('#close-map').click()
        with page.expect_response('**/sanger-report?*'):
            group_pending[0].fulfill(status=200,content_type='application/json',body=json.dumps(refreshed))
        page.evaluate('new Promise(r=>requestAnimationFrame(()=>requestAnimationFrame(r)))')
        assert not group_late
        page.get_by_role('button',name='Open map').click()
        expect(high_card.locator('.analysis-version')).to_contain_text('Run 2 · Latest')
        # Grouped indel anchors retain reverse-read coordinates without an indel verdict.
        page.locator('#sanger-file').set_input_files(str(indel_path))
        page.get_by_role('button',name='Attach AB1',exact=True).click()
        indel_card=page.locator('.read-card').filter(has_text=indel_path.name)
        indel_card.get_by_role('button',name='Analyze against this revision').click()
        expect(indel_card.locator('.analysis-version')).to_contain_text('Run 1')
        page.locator('.indel-review summary').click()
        expect(page.locator('.indel-event')).to_have_count(1)
        expect(page.locator('.indel-event')).to_contain_text('insertion · 3 bp')
        expect(page.locator('.indel-event')).to_contain_text('Flanks meet checks')
        expect(page.locator('.indel-event')).to_contain_text('Reference boundary after 650')
        expect(page.locator('.indel-review')).to_contain_text('compared only with explicitly observed reference spans')
        comparison=page.locator('.indel-comparison')
        expect(comparison).to_have_count(1)
        expect(comparison).to_contain_text('Conflicting event / reference support')
        expect(comparison.locator('.indel-support-counts')).to_contain_text('1Event-supporting reads')
        expect(comparison.locator('.indel-support-counts')).to_contain_text('2Reference-supporting reads')
        unassessed_count=1+int(bool(public_fixture))
        expect(comparison.locator('.indel-support-counts')).to_contain_text(f'{unassessed_count}Unassessed reads')
        expect(comparison.locator('li').filter(has_text=high_path.name)).to_contain_text('not fully observed')
        comparison.screenshot(path=str(artifacts/'indel-comparison-desktop.png'))
        comparison.locator('li').filter(has_text=read_path.name).get_by_role('button',name='Left · read 500').click()
        expect(page.locator('.read-card').filter(has_text=read_path.name).locator('.trace-selected')).to_have_attribute('data-peak-position','499')

        page.locator('.indel-event').screenshot(path=str(artifacts/'indel-anchors-desktop.png'))
        page.locator('.indel-flanks button').first.click()
        expect(indel_card.locator('.trace-selected')).to_have_attribute('data-peak-position','153')
        page.set_viewport_size({'width':390,'height':844})
        page.locator('.indel-event').evaluate('(el)=>el.scrollIntoView({block:"start"})')
        page.evaluate('new Promise(r=>requestAnimationFrame(()=>requestAnimationFrame(r)))')
        page.screenshot(path=str(artifacts/'indel-anchors-mobile.png'))
        assert page.locator('.indel-event').evaluate('(el)=>el.scrollWidth<=el.clientWidth+1')
        comparison.evaluate('(el)=>el.scrollIntoView({block:"start"})')
        page.screenshot(path=str(artifacts/'indel-comparison-mobile.png'))
        assert comparison.evaluate('(el)=>el.scrollWidth<=el.clientWidth+1')
        with page.expect_download() as indel_export:
            page.get_by_role('button',name='Group JSON',exact=True).click()
        indel_export_path=artifacts/'indel-group.json';indel_export.value.save_as(str(indel_export_path))
        exported_events=json.loads(indel_export_path.read_text())['snapshot']['summary']['quality_review']['indel_review']['events']
        assert len(exported_events)==1 and exported_events[0]['eligible_for_exact_anchor_review']
        assert exported_events[0]['original_read_positions']==[152,151,150]
        group_evidence=json.loads(indel_export_path.read_text())['snapshot']['summary']['quality_review']
        assert group_evidence['indel_conflicts']=='exact_event_vs_reference_only'
        assert group_evidence['indel_review']['comparisons'][0]['reference_read_count']==2
        if public_fixture:
            public_source=next(s for s in group_evidence['indel_review']['comparisons'][0]['sources'] if s['filename']==public_path.name)
            assert public_source['support']=='unassessed' and public_source['reason']
            assert group_evidence['indel_review']['comparisons'][0]['unassessed_read_count']==unassessed_count
        with page.expect_download() as indel_html:
            page.get_by_role('button',name='Group HTML',exact=True).click()
        indel_html_path=artifacts/'indel-comparison.html';indel_html.value.save_as(str(indel_html_path))
        indel_offline=browser.new_context(offline=True,viewport={'width':1100,'height':1000})
        indel_page=indel_offline.new_page();indel_page.goto(indel_html_path.as_uri())
        expect(indel_page.locator('.indel-comparison-report')).to_contain_text(f'1 event / 2 reference / {unassessed_count} unassessed')
        indel_page.locator('.indel-comparison-report').screenshot(path=str(artifacts/'indel-comparison-offline.png'))
        indel_offline.close()

        # Input checks are explicit queue tasks; the test starts one worker on this host.
        checks=page.locator('#input-jobs')
        page.set_viewport_size({'width':1440,'height':1080})
        checks.get_by_role('button',name='Refresh checks',exact=True).click()
        expect(checks.locator('.ij-status')).to_have_text('Status updated.')
        choices=checks.locator('.ij-read-list input')
        choices.nth(0).check();choices.nth(1).check()
        # Server commits but the response is lost: retry must reuse the same key.
        def lost_submission(route):
            route.fetch()
            route.abort()
        page.route('**/api/revisions/*/input-checks',lost_submission,times=1)
        checks.get_by_role('button',name='Queue input check',exact=True).click()
        expect(checks.locator('.ij-status')).to_contain_text('Could not complete request')
        checks.get_by_role('button',name='Queue input check',exact=True).click()
        expect(checks.locator('.ij-job')).to_have_count(1)
        expect(checks.locator('.ij-job')).to_contain_text('Queued · awaiting worker')
        job_id=checks.locator('.ij-job').get_attribute('data-ij-job')
        checks.screenshot(path=str(artifacts/'input-check-queued-desktop.png'))
        completed=subprocess.run([sys.executable,'-m','localmolbio.worker','once','--worker-id','browser-validation'],env=environment,cwd=ROOT,capture_output=True,text=True,check=True)
        assert json.loads(completed.stdout)['status']=='succeeded'
        checks.get_by_role('button',name='Refresh checks',exact=True).click()
        expect(checks.locator('.ij-job')).to_contain_text('Input check complete')
        checks.get_by_role('button',name='Details & attempts',exact=True).click()
        expect(checks.locator('.ij-details')).to_contain_text('Attempt 1 · succeeded')
        expect(checks.locator('.ij-details')).to_contain_text('No sequence analysis performed')
        checks.screenshot(path=str(artifacts/'input-check-complete-desktop.png'))
        page.set_viewport_size({'width':390,'height':844})
        checks.evaluate('(el)=>el.scrollIntoView({block:"start"})')
        page.screenshot(path=str(artifacts/'input-check-mobile.png'))
        assert checks.evaluate('(el)=>el.scrollWidth<=el.clientWidth+1')
        checks.locator('.ij-details').screenshot(path=str(artifacts/'input-check-history-mobile.png'))
        choices.nth(0).check()
        checks.get_by_role('button',name='Queue input check',exact=True).click()
        expect(checks.locator('.ij-job')).to_have_count(2)
        queued=checks.locator('.ij-job').filter(has_text='Queued · awaiting worker')
        queued.get_by_role('button',name='Cancel check',exact=True).click()
        expect(checks.locator('.ij-job').filter(has_text='Cancelled')).to_have_count(1)
        page.route('**/api/revisions/*/input-checks?*',lambda route:route.abort(),times=1)
        checks.get_by_role('button',name='Refresh checks',exact=True).click()
        expect(checks.locator('.ij-status')).to_contain_text('Could not complete request')
        checks.locator('.ij-status').scroll_into_view_if_needed()
        page.screenshot(path=str(artifacts/'input-check-error-mobile.png'))
        assert checks.get_by_role('button',name='Refresh checks',exact=True).is_enabled()
        checks.get_by_role('button',name='Refresh checks',exact=True).click()
        expect(checks.locator('.ij-status')).to_have_text('Status updated.')
        checks.screenshot(path=str(artifacts/'input-check-cancelled-mobile.png'))
        # A delayed detail response must not repopulate a closed construct.
        pending_details=[]
        page.route('**/api/revisions/*/input-checks/*?*',lambda route:pending_details.append(route))
        checks.locator(f'[data-ij-job="{job_id}"]').get_by_role('button',name='Details & attempts',exact=True).click()
        page.wait_for_timeout(100)
        assert len(pending_details)==1
        detail_response=pending_details[0].fetch()
        page.locator('#close-map').click()
        pending_details[0].fulfill(response=detail_response)
        page.unroute('**/api/revisions/*/input-checks/*?*')
        page.wait_for_timeout(100)
        expect(checks).to_be_empty()
        page.get_by_role('button',name='Open map').click()
        expect(checks.locator('.ij-job')).to_have_count(2)
        expect(checks.locator('.ij-details')).to_be_empty()

        # Populate real persisted attempts to exercise both history pages in the UI.
        listed=httpx.get(base+'/api/jobs').json()
        check_revision=next(j['sequence_revision_id'] for j in listed if j['id']==job_id)
        detail=httpx.get(f'{base}/api/revisions/{check_revision}/input-checks/{job_id}').json()
        history_job=httpx.post(f'{base}/api/revisions/{check_revision}/input-checks',json={'read_ids':[detail['inputs']['inputs']['reads'][0]['id']],'idempotency_key':'ui-seven-attempts','max_attempts':10}).json()['id']
        retry_code="from localmolbio import job_queue as q; from localmolbio.input_validation import ADAPTER\nfor i in range(7):\n c=q.claim('history-worker',ADAPTER); q.fail(c['job_id'],c['lease_token'],'transient_test_io',True)"
        subprocess.run([sys.executable,'-c',retry_code],env=environment,cwd=ROOT,check=True)
        checks.get_by_role('button',name='Refresh checks',exact=True).click()
        expect(checks.locator('.ij-job')).to_have_count(3)
        history_card=checks.locator(f'[data-ij-job="{history_job}"]')
        history_card.get_by_role('button',name='Details & attempts',exact=True).click()
        expect(checks.locator('.ij-details li')).to_have_count(5)
        expect(checks.locator('.ij-details li').first).to_contain_text('Attempt 7')
        checks.get_by_role('button',name='Older attempts',exact=True).click()
        expect(checks.locator('.ij-details li')).to_have_count(2)
        expect(checks.locator('.ij-details li').last).to_contain_text('Attempt 1')
        checks.get_by_role('button',name='Newer attempts',exact=True).click()
        expect(checks.locator('.ij-details li')).to_have_count(5)
        checks.locator('.ij-details').screenshot(path=str(artifacts/'input-check-retries-mobile.png'))
        checks.get_by_role('button',name='Older attempts',exact=True).scroll_into_view_if_needed()
        page.screenshot(path=str(artifacts/'input-check-history-controls-mobile.png'))
        history_card.get_by_role('button',name='Cancel check',exact=True).click()
        expect(checks.locator('.ij-job').filter(has_text='Cancelled')).to_have_count(2)

        # A failed list request must be shown, not reported as an empty library.
        page.locator('#close-map').click()
        page.route('**/api/sequence-revisions/*/sanger-reads*', lambda route: route.abort())
        page.get_by_role('button',name='Open map').click()
        expect(page.locator('#sanger-status')).to_contain_text('Could not load reads')
        expect(page.locator('#sanger-summary')).to_be_empty()
        assert not errors, errors
        browser.close()
    print(json.dumps({'passed':True,'artifacts':str(artifacts),'public_ab1_checked':bool(public_fixture),'checks':['input-check submission retry, worker completion, history, cancellation, errors and late-response isolation','exact indel vs reference support, unassessed short read, source peak links and offline group report','grouped reverse insertion, Q20 flank audit, original peak jump and snapshot export','group JSON/HTML downloads, offline rendering, stale snapshot rejection and late-response cancellation','quality-filtered multi-read overlap, conflict and original peak navigation','base-level Q20 mapping, reverse original coordinates and historical absence','high-quality difference shown separately from matching evidence','latest-only multi-read coverage, overlap and circular gap review','binary ABIF upload','invalid and duplicate upload','analysis and persisted evidence','append-only rerun and historical report switching','six-run history pagination','optional end trimming with original and retained coverage','historical and current JSON/HTML downloads','offline HTML rendering without external requests','export failure recovery and late-response cancellation','reverse trimmed variant and boundary jump to original peaks','history network failure and stale rerun recovery','Q12 variant and focused chromatogram','desktop and mobile layouts','network failure','no JS exceptions']},indent=2))
finally:
    server.terminate()
    server.wait(timeout=10)
    log.close()
