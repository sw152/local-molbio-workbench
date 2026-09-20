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
- Circular plasmid-map preview with feature tracks and source-coordinate review status
- Primer3-backed PCR candidate generation with persisted design parameters and evidence

Real Benchling exports are private research data. Keep them outside the repository. The test suite uses synthetic records only.

The product boundary, workstation deployment model, data model, and staged acceptance criteria are documented in [the Chinese architecture and roadmap](docs/architecture-and-roadmap.zh-CN.md).

The continuously executed delivery plan is in [the Chinese execution plan](docs/execution-plan.zh-CN.md).
