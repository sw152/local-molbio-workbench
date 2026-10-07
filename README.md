# Local Molecular Biology Workbench

An open-source, local-first workbench for plasmid design, primer management, and sequencing verification.

The first vertical slice imports a Benchling GenBank ZIP safely, preserving duplicate archive members and recording enough provenance to make the import reproducible. It keeps parser warnings with each record so a questionable feature annotation is visible for review instead of being silently trusted.

## Run locally

Python 3.11 or later is recommended. The prototype also supports Python 3.9+.

```bash
python -m venv .venv
. .venv/bin/activate
pip install -e '.[dev]'
uvicorn localmolbio.api:app --host 127.0.0.1 --port 8000
```

Open `http://127.0.0.1:8000` in a browser. Runtime data is stored in `./var` by default; set `MOLBIO_DATA_DIR` to choose another location.

## Workstation deployment

The application binds to loopback only. Run it on the workstation and expose its local port privately through Tailscale Serve:

```bash
tailscale serve --bg localhost:8000
```

Do not expose the database or internal worker ports to the LAN or public internet.

For the complete Tailscale SSH, Serve, and future remote-worker setup, read [the workstation guide](docs/workstation-tailscale.zh-CN.md).

## Current scope

- Safe ZIP inspection and GenBank import
- Archive-level and member-level SHA-256 provenance
- Duplicate archive-member preservation
- Circular/linear sequence metadata and GenBank feature counts
- Per-record parser-warning retention for annotation review
- Searchable local sequence inventory
- Circular and linear maps with strand arrows, coordinate ticks, overlap-aware tracks, zoom and focused annotation review
- Primer3-backed PCR candidate generation with persisted design parameters and evidence
- Revision-scoped AB1 upload, local Sanger alignment, coverage and quality-aware difference review
- Windowed four-channel chromatograms, original-read navigation, and jumps from differences to their peak evidence

Real Benchling exports are private research data. Keep them outside the repository. The core test suite uses synthetic records. An optional browser check can also inspect a separately supplied public AB1 fixture; it is not bundled.

The product boundary, workstation deployment model, data model, and staged acceptance criteria are documented in [the Chinese architecture and roadmap](docs/architecture-and-roadmap.zh-CN.md).

The continuously executed delivery plan is in [the Chinese execution plan](docs/execution-plan.zh-CN.md).

## Review Sanger evidence

Open a sequence, attach an AB1 file, then analyze it against the current revision.
The read card shows local alignment identity, reference coverage and review flags.
“View chromatogram” displays the stored analyzed signal and Phred calls. Selecting
an original-read position in the difference table centers that base in the trace.
The trace always uses the uploaded read orientation, including for reverse alignments.
Missing or inconsistent trace metadata is reported; no replacement signal is generated.
High identity alone does not verify the whole plasmid.

For reproducible browser checks, install `.[dev,ui]` and run
`python scripts/verify_sanger_ui.py`. It uses an installed macOS Chrome or Playwright
Chromium (install the latter with `python -m playwright install chromium` if needed).
Its synthetic files, isolated database, logs and screenshots remain under `var/ui-check/`.
Set `MOLBIO_PUBLIC_AB1` to a separately acquired public fixture path to additionally
validate its upload and trace display.

## Inspect sequence maps

Select a feature in the diagram or annotation list to inspect its coordinates.
Overlapping annotations occupy separate tracks; compound features retain their
segments on one track. For a crowded map, “Focus selected” isolates one annotation
while preserving the complete list; “Show all” restores the overview. Zoom and Fit
control the diagram. Narrow linear diagrams scroll horizontally to retain readable labels.
Tick labels are zero-based boundary positions; annotation ranges are one-based base positions.
Unknown topology uses a clearly labeled linear coordinate view.

Geometry checks: `node --test tests/ui/sequence-map.test.cjs`.
Browser checks: `python scripts/verify_sequence_map_ui.py` with the UI dependencies above.
