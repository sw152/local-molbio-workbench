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
- Paginated sequence library with literal name/source search, total counts and recoverable loading errors
- Circular and linear maps with strand arrows, coordinate ticks, overlap-aware tracks, zoom and focused annotation review
- Primer3-backed PCR candidate generation, optional target flanking, persisted design parameters and evidence
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

“Run new analysis” saves another report without replacing previous evidence. You
can change the requested orientation for the new run; the upload metadata stays
unchanged. “Browse saved analyses” pages through history and opens old coverage
and difference reports, clearly labeled as historical. Each new report records its
method version, effective parameters, input hashes and predecessor.

End quality trimming is optional and defaults to off. When enabled, it removes
only the contiguous end calls below the chosen Phred threshold, stopping at the
first passing call on each side. Interior low-quality and ambiguous calls remain
in the retained interval and in the evidence. This is a simple endpoint rule,
not a sliding-window, Mott or mixed-peak quality algorithm. Fewer than 12 retained
bases, missing qualities or invalid thresholds are rejected without saving a report.

The report records the retained original-read interval, left/right exclusions,
method and threshold. It separately displays aligned fractions of the **original**
read and **retained** interval; a high retained fraction cannot establish full-read
or whole-plasmid verification. The trim bar always follows the original uploaded
orientation. Both retained boundary buttons and variant links open the original
chromatogram, including for reverse alignments. Original calls, qualities and files
are never trimmed in storage. Older reports without these fields remain marked as
unrecorded rather than being silently recalculated.

Every run reparses the hash-checked original AB1 and requires its calls and existing
saved qualities to agree. Missing historical per-base qualities can be recovered
for the new run. Missing or altered files, inconsistent calls and failed analyses
leave saved reports intact. Existing one-report databases migrate additively into
analysis history on startup, preserving the legacy table and identifiers.

API: `POST /api/sanger-verifications` accepts `sequencing_read_id`, optional
`direction`, optional `trim_quality_threshold` (integer 1–60; null/off by default),
and `previous_alignment_id` (the latest report ID, required for reruns).
Evidence version 3 keeps `read_position`/`read_start`/`read_end` on the full oriented
original read; `analysis_read_position` is relative to the retained, oriented read.
`original_read_position` and the trimming interval use original uploaded coordinates,
all zero-based. `read_aligned_fraction` always uses the original length as denominator;
`retained_read_aligned_fraction` uses the retained length.
Stale predecessors return 409, including when another request finishes first.
`GET /api/sanger-reads/{id}/analyses?limit=10&offset=0` returns newest-first reports
with pagination metadata. The revision read list continues to show one row per
read, with its latest report.

For reproducible browser checks, install `.[dev,ui]` and run
`python scripts/verify_sanger_ui.py`. It uses an installed macOS Chrome or Playwright
Chromium (install the latter with `python -m playwright install chromium` if needed).
Its synthetic files, isolated database, logs and screenshots remain under `var/ui-check/`.
Set `MOLBIO_PUBLIC_AB1` to a separately acquired public fixture path to additionally
validate its upload and trace display.

## Export one saved Sanger analysis

Each displayed report has **Download HTML** and **Download JSON** controls. Switch
to a historical run first to export that run; the export always uses its explicit
analysis ID. HTML is a self-contained, offline-readable report with coverage and
retained-read diagrams, differences, input provenance, parameters and limitations.
JSON preserves the complete saved alignment, including its original coordinate
conventions. Neither format contains raw AB1 signals or original source files. Evidence v4
includes aligned read/reference bases and per-base qualities; those aligned bases
may span the entire input.

New analysis jobs snapshot input names, lengths and reference topology alongside
the hashes. For older jobs, linked-record metadata is labeled as such, and missing
analysis-time hashes remain null/“Not recorded”. Export does not recheck source
files or rerun the alignment, so an archived report remains exportable if its AB1
is unavailable. It must not be interpreted as a fresh source-integrity check.

`GET /api/sanger-reads/{read_id}/analyses/{analysis_id}/export?format=json` supports
`json` (default) or `html`, with attachment headers and no-store caching. The
analysis must belong to the requested read. User-provided names and all stored
text are escaped in HTML; filenames use analysis identifiers. Reports contain no
scripts or external resources. Export schema version 1 is separate from the saved
alignment evidence version. The optional Sanger browser check validates real
downloads, historical run binding, error recovery, late-response cancellation and
offline rendering on desktop and narrow screens.

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

The library offers 24, 48 or 96 records per page. Searches reset to the first page;
late responses from an older search or construct are ignored. The paginated API is
`GET /api/sequence-library?query=...&limit=24&offset=0`, returning `items`, `total`,
`offset`, `limit` and `has_more`. The existing `/api/sequences` array response remains
available for compatibility. Search metacharacters in the paginated endpoint are literal.
Run `python scripts/verify_library_ui.py` for pagination, stale-response and error-recovery checks.

## Design primers around a target

In the PCR form, enter the first and last target bases (one-based, inclusive),
or leave both empty for unconstrained placement. Returned oligos must flank the
entire target and match the reference on their respective strands. Targets require
space for a primer on either side. Circular-origin-spanning products are not supported;
circular records currently use their stored linear coordinates for primer design.

The API accepts optional `target_start` and `target_end` as a zero-based half-open
interval. Saved candidates retain the target, design job, engine version and scope.
Selected-primer CSV exports also include target coordinates, pair index, product
length, exact reference match count and specificity status. New designs retain an
exact, full-length oligo scan on both strands of the current reference, including
origin-crossing sites for circular records. Counts are directional: palindromic
oligos can match both strands at one interval. The stored site list is capped at
100 per oligo; the total count remains complete and truncation is explicit.
Unknown topology is scanned as linear, with the origin marked unchecked.

The pair review panel shows opposing primer arrows, the requested target and
expected product at their actual reference coordinates. Its axis uses zero-based
boundaries; site labels use one-based inclusive ranges. Pairs are associated by
design job and pair index, and unlinked legacy candidates are not guessed into pairs.
Mobile diagrams scroll horizontally. Older candidates without site evidence are
marked as not reviewed.

New designs also enumerate exact, inward-facing, non-overlapping site pairs using
one of each oligo, in both F→/←R and R→/←F assignments. The search uses the requested
product size range (inclusive), and circular intervals traverse the reference at
most once. Products are deduplicated by reference start and span, with primer
assignments retained. The intended interval and additional intervals are displayed
with coordinate bars, directions and lengths. Unknown topology leaves the origin
unchecked. This product review can detect origin-crossing alternatives even though
primer design itself still uses linear coordinates.

If either stored binding-site list was truncated, product counts are explicitly
lower bounds and total/alternative counts remain null; no absence claim is made.
Up to 100 product intervals are stored, with the intended interval first when
found. When only this output list is truncated, complete counts are still retained.
Saved metrics, job parameters and selected-primer CSV retain the review scope,
size range and completeness. Historical candidates are not silently reanalyzed.

**Specificity is not established** by these checks. Mismatch binding, overlapping
sites, same-oligo pairs, other product sizes, other references and whole-genome
specificity are not evaluated. These are geometric candidates, not predictions of
amplification success or experimentally validated primers.
The strand and primer-inclusive size conventions follow the basic description in
[UCSC In-Silico PCR](https://genome.ucsc.edu/cgi-bin/hgPcr); this workbench uses its
own bounded exact-match implementation and does not run UCSC software.

Run `python scripts/verify_primer_ui.py` to check target entry, actual Primer3
execution on synthetic data, selection, CSV export and persisted target metadata.

## Review combined read coverage

The Sanger panel summarizes the latest saved analysis of each read attached to
one reference revision. It displays the union of aligned reference positions,
uncovered regions, overlap depth and every included or excluded source run.
Circular gaps at the two coordinate ends are shown as one region across the
origin. Reanalysis history never increases read depth.

This is geometric alignment coverage, **not quality-screened coverage or a
consensus**. Low-quality, ambiguous and differing calls may contribute aligned
positions; their saved warnings remain visible. The separate quality panel reviews Q20 paired-base conflicts;
even 100% geometric coverage does not mean whole-plasmid verification. Ambiguous
placements, truncated searches, missing/mismatched input hashes, unsupported
reports and invalid intervals are excluded. Duplicate source hashes exclude all
copies. Original files are not rechecked when viewing the summary.

`GET /api/sequence-revisions/{id}/sanger-reads?include_summary=true` returns
`{reads, summary}` from one database snapshot; without that option the original
array response is preserved. The summary includes source analysis IDs, hashes,
method version and a deterministic snapshot SHA-256. Changing the historical
report displayed in a read card does not change this latest-run summary.

## Inspect per-base quality mappings

New Sanger analyses use evidence v4. `evidence.base_mapping` stores compact blocks
of aligned reference/read base pairs: reference start/end, original uploaded-read
start and step (+1 or −1), paired reference/call strings, and a Phred array with
`null` for unavailable quality. Reverse calls are complemented into reference
orientation, while their coordinates and qualities still refer to the original
AB1. Blocks split at indels and circular origin crossings. Insertions/deletions
remain in the saved difference list and are not assigned invented paired bases.

The read panel separates Q20-or-higher unambiguous matching and differing calls.
Other paired calls are classified as ambiguous first, then quality unavailable,
then below Q20. Q20 itself is included. These are counts for the reported local
placement; ambiguity/truncated-search flags still apply. This filter is distinct
from optional end trimming and does not remove calls from saved evidence.

Historical v2/v3 reports remain unchanged and show that the mapping was not
recorded. JSON and standalone HTML's saved-evidence section retain the new
mapping, including aligned bases, which can encompass the entire input. The
combined-read panel retains geometric coverage alongside a separate Q20 paired-base
conflict review. Neither produces a consensus or whole-plasmid validation verdict.

## Review quality-filtered coverage and disagreements

The revision summary now adds `quality_review`, computed from the same database
snapshot as the read cards. Only mappings that pass the existing placement and
provenance filters are considered. The reference sequence hash, mapping version,
block lengths, original coordinates/direction, paired bases, saved quality values
when available, covered intervals and derived statistics must agree. A malformed
mapping is excluded in full. Older reports without a mapping remain explicitly
unassessed; they can still contribute geometric coverage.

The Q20 panel separately displays positions supported by unambiguous high-quality
calls, positions with at least two such calls, non-reference positions and
conflicting positions. Two different Q20+ alleles at the same reference position
are a conflict. Reads agreeing on an alternate allele are non-reference evidence,
not a conflict between reads. Low-quality/N calls do not establish agreement.
Each listed call retains its read/run identity, Phred and original AB1 coordinate;
its button opens that original chromatogram peak, including reverse-read positions.

This is a paired-base review, not variant normalization or consensus generation.
Indel support has a separate, limited exact-event-versus-reference review below.
Other indel alleles and repeat-shiftable events remain unassessed. Zero observed conflicts does
not assess missing, excluded or low-quality evidence. Original files are not
reopened for the summary; legacy stored quality arrays may be absent, in which
case the saved per-base analysis evidence remains the quality source.

## Export a group evidence snapshot

“Group JSON” and “Group HTML” export the current revision's complete displayed
read-group snapshot: each attached read's latest saved analysis, both coverage
layers, gaps, conflicts, excluded sources, analysis parameters and base mappings.
The standalone HTML works offline and includes expandable full evidence; JSON
is the machine-readable copy. No raw AB1 signals or source files are bundled.
Aligned bases in saved mappings can still encompass entire inputs.

The summary response now includes a top-level `snapshot_sha256` over **both**
`reads` and `summary`. It is distinct from the summary-only digest. Export via
`GET /api/sequence-revisions/{id}/sanger-report?expected_snapshot={snapshot_sha256}&format=json`
(or `format=html`). The server captures one database view and rejects a changed
snapshot with HTTP 409. Reopen the construct to inspect the new state before
exporting; previously downloaded snapshots remain unchanged. Closing or switching
constructs cancels delivery of pending downloads.

The digest is SHA-256 over UTF-8 JSON of `{reads,summary}`, with sorted keys,
compact separators and `ensure_ascii=False`; it excludes export time and the
outer digest itself. The export records that method, full digest and timestamp.
This is content integrity, not a signature or a biological validation verdict.
It exports the latest run per read rather than every historical run. Use each
read's historical report export when a different saved run is required.

## Inspect insertion and deletion anchors

`quality_review.indel_review` groups contiguous saved insertion/deletion entries
for reads whose base mappings passed validation. Each event retains its source
run and variant indices, sequence, reference boundary/segments, original inserted
positions/qualities and both mapped flanks. The UI shows a schematic and links
the flanks back to the original chromatogram, including reverse-read coordinates.

An event meets the current exact-anchor checks only when both flanks are
unambiguous reference matches at Q20+, inserted calls also meet Q20, and the
original-read coordinate walk forms a single simple gap. Missing flanks,
ambiguous/low-quality calls, inconsistent or duplicated saved entries, and
possible equivalent one-base shifts in repeats are withheld with reasons.
Deletion bases have no directly measured per-base Phred; flank quality must not
be relabeled as deletion quality. Linear ends never acquire a fabricated
opposite flank, and circular origin events preserve split reference segments.

Event/anchor inspection feeds the limited comparison described below. No event
normalization is performed. Passing the anchor checks is not a variant-confidence
score or a consensus verdict.
The complete audit is included in group JSON and the HTML's full snapshot;
original saved reports are not rewritten.

## Compare exact indels with observed reference support

Eligible events are grouped by kind, exact reference boundary and sequence.
Each attached read contributes at most once per group: **event**, **reference**
or **unassessed**. Reference support for an insertion requires its two adjacent
reference flanks; for a deletion it requires every deleted reference base plus
both flanks. Every compared call must match unambiguous reference bases at Q20+
and advance consecutively in the original read (decreasing for reverse reads).
An absent event in a report alone is never reference evidence.

The group exposes `conflicting_support` when distinct eligible reads support the
exact event and the reference span. It records source run IDs, original flanks,
minimum reference-span quality and every unassessed reason. Different insertion
sequences/other indel alleles are not counted as reference, and are not compared
to one another. Repeat-shiftable, complex or inconsistent evidence stays withheld.
The API marks this scope as `exact_event_vs_reference_only`; there is no event
normalization, consensus or whole-plasmid pass verdict. UI source links open the
original peaks, and group JSON/HTML exports retain the comparisons.
