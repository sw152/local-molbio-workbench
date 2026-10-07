"""Technical rehearsal on frozen public inputs; no human usability measurements."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import socket
import subprocess
import sys
import time

import httpx
from playwright.sync_api import sync_playwright, expect

ROOT = Path(__file__).resolve().parents[1]
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--pack', type=Path, required=True)
args = parser.parse_args()
pack = args.pack.resolve()
manifest_bytes = (pack / 'manifest.json').read_bytes()
manifest_hash = hashlib.sha256(manifest_bytes).hexdigest()
assert manifest_hash == (pack / 'manifest.sha256').read_text().strip()
manifest = json.loads(manifest_bytes)
assert manifest['synthetic_only'] and manifest['status'] == 'technical_rehearsal_only'
for item in manifest['files']:
    assert hashlib.sha256((pack / item['path']).read_bytes()).hexdigest() == item['sha256']
artifacts = ROOT / 'var' / 'validation-pack-ui' / time.strftime('%Y%m%d-%H%M%S')
artifacts.mkdir(parents=True)
with socket.socket() as sock:
    sock.bind(('127.0.0.1', 0))
    port = sock.getsockname()[1]
base = f'http://127.0.0.1:{port}'
log = (artifacts / 'server.log').open('w')
app_python = os.environ.get('MOLBIO_VERIFY_PYTHON', sys.executable)
server = subprocess.Popen([app_python, '-m', 'localmolbio', 'serve', '--port', str(port)],
                          cwd=ROOT, env=dict(os.environ, MOLBIO_DATA_DIR=str(artifacts / 'runtime')), stdout=log, stderr=log)
try:
    for _ in range(100):
        try:
            if httpx.get(base + '/health').status_code == 200: break
        except httpx.TransportError: pass
        time.sleep(.1)
    else: raise RuntimeError('service did not become ready')
    with sync_playwright() as p:
        browser = p.chromium.launch(headless=True, executable_path='/Applications/Google Chrome.app/Contents/MacOS/Google Chrome')
        page = browser.new_page(viewport={'width':1440, 'height':1080})
        errors = []
        page.on('pageerror', lambda error: errors.append(str(error)))
        page.goto(base)
        page.locator('#archive').set_input_files(str(pack / 'participant' / 'reference-library.zip'))
        page.get_by_role('button', name='Import into local library').click()
        expect(page.locator('#status')).to_contain_text('Imported 10 records')
        sequences = page.request.get(base + '/api/sequences').json()
        # The list endpoint shape is read from the existing API, not inferred from titles.
        by_name = {s['display_name']: s['id'] for s in sequences}
        for case, original_position, q in [('A02', 359, 38), ('A03', 200, 9)]:
            page.reload()
            page.locator(f'.open[data-id="{by_name[case]}"]').click()
            expect(page.locator('#sanger-reads')).to_contain_text('No reads attached')
            expect(page.locator('#sanger-file')).to_be_enabled()
            page.locator('#sanger-file').set_input_files(str(pack / 'participant' / 'inputs' / f'{case}-r1.ab1'))
            page.get_by_role('button', name='Attach AB1', exact=True).click()
            expect(page.locator('#sanger-reads')).to_contain_text('Ready to analyze')
            page.get_by_role('button', name='Analyze against this revision').click()
            expect(page.locator('#sanger-status')).to_contain_text('Analysis saved')
            page.locator('.variant-details summary').click()
            expect(page.locator('.variant-scroll tbody tr')).to_have_count(1)
            page.locator('.trace-jump').click()
            expect(page.locator('.trace-selected')).to_have_attribute('data-peak-position', str(original_position))
            with page.expect_download() as download:
                page.get_by_role('button', name='Download JSON', exact=True).click()
            output = artifacts / f'{case}-report.json'
            download.value.save_as(str(output))
            report = json.loads(output.read_text())['saved_alignment']
            assert report['variants'][0]['original_read_position'] == original_position
            assert report['variants'][0]['phred'] == q
            assert report['evidence']['whole_reference_verified'] is False
            page.locator('.trace-view').screenshot(path=str(artifacts / f'{case}-trace.png'))
            page.set_viewport_size({'width':390, 'height':844})
            page.locator('.variant-details').scroll_into_view_if_needed()
            page.screenshot(path=str(artifacts / f'{case}-mobile.png'))
            assert page.locator('#map-panel').evaluate('(e) => e.scrollWidth <= e.clientWidth')
            page.set_viewport_size({'width':1440, 'height':1080})
        page.reload()
        page.locator(f'.open[data-id="{by_name["A09"]}"]').click()
        expect(page.locator('#sanger-reads')).to_contain_text('No reads attached')
        expect(page.locator('#sanger-file')).to_be_enabled()
        page.locator('#sanger-file').set_input_files(str(pack / 'participant' / 'inputs' / 'A09-r1.ab1'))
        page.get_by_role('button', name='Attach AB1', exact=True).click()
        expect(page.locator('#sanger-status')).to_contain_text('Could not attach read')
        expect(page.locator('#sanger-reads')).to_contain_text('No reads attached')
        assert not errors, errors
        browser.close()
    result = {'record_type':'technical_rehearsal', 'pack_id':manifest['pack_id'], 'manifest_sha256':manifest_hash,
              'passed':True, 'human_sessions':0, 'benchling_comparison':'not_run',
              'checks':['10 reference imports', 'reverse original peak and report', 'low-quality original peak and report', 'invalid AB1 rejected', 'desktop/mobile rendering']}
    (artifacts / 'result.json').write_text(json.dumps(result, indent=2))
    print(json.dumps({'passed':True, 'artifacts':str(artifacts)}))
finally:
    server.terminate()
    server.wait(timeout=15)
    log.close()
