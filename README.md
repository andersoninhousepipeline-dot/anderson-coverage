# Anderson Coverage Checker

A web tool to check **panel/capture coverage** and **sample read-depth coverage**
for NGS target (BED) panels — by gene, multiple genes, chromosome region,
transcript (NM/NR/XM/XR/ENST/ENSG), or rs ID — with one-click PDF reports.

Two layers of coverage in every report:

1. **Panel coverage** — is the gene / region / variant in the capture design, and
   how much of it (from the BED target intervals).
2. **Sample read-depth coverage** *(always included)* — real sequencing depth from
   reference BAM files over those targets: mean / median / min depth and
   **% of target bases ≥ 1× / 10× / 20× / 30× / 50× / 100×**, plus a **Combined**
   (average across all reference samples) row.

## Query types (auto-detected)

| Type         | Example                          |
|--------------|----------------------------------|
| Gene         | `BRCA1`                          |
| Multiple genes | `BRCA1, BRCA2, TP53`           |
| Chr region   | `chr17:43044295-43125483`        |
| Transcript   | `NM_007294.4`                    |
| rs ID        | `rs6265` (resolved via Ensembl REST) |

PDF filename follows the query: `BRCA1 coverage.pdf`, `listed gene coverage.pdf`,
`rs6265 coverage.pdf`, etc. (Print / Save PDF button).

## Setup

```bash
git clone https://github.com/andersoninhousepipeline-dot/anderson-coverage.git
cd anderson-coverage
pip install -r requirements.txt

# Point the app at your data (defaults shown):
export BED_DIR=/data/bed            # folder of *.bed panels
export BAM_DIR=/data/bed/reference  # folder of indexed reference *.bam

./start.sh                          # http://<host>:8100   (PORT=8101 ./start.sh to change)
```

Requires Python 3 with Flask, requests, pysam, numpy (see `requirements.txt`).
rs-ID lookups call `rest.ensembl.org` (needs internet). BAMs must be indexed (`.bai`).

## Server control

```bash
./start.sh    ./status.sh    ./restart.sh    ./stop.sh
```

`status.sh` reports RUNNING/STOPPED, uptime, memory, HTTP health and panel count.

## Panels (BED files)

Every `*.bed` under `BED_DIR` is auto-discovered and selectable in the dropdown.
Mixed vendor annotation styles are supported (Twist `Gene;NM_…`, Roche
`gene_symbol=…`, Sophia `Gene:NM:exon`, comma / tab / plain). Users can also
**upload a BED** in the UI. Region and rs-ID queries work on any BED; gene /
transcript search needs an annotated 4th column.

## Reference samples

Reference samples live as **one sub-directory per type** under `BAM_DIR`, each
holding one or more replicate BAMs (indexed `.bai` required). No BAMs or sample
identifiers are committed to the repo — only the directory layout matters:

```
BAM_DIR/
  Normal-Male/    *.bam   # +.bai
  Normal-Female/  *.bam
  Male-Inf/       *.bam
  Female-Inf/     *.bam
  AF/             *.bam
  POC/            *.bam
```

Every report also includes a **Case vs reference** panel: each case type
(Male/Female Infertility, AF, POC) is flagged when its mean depth falls below the
sex-matched Normal reference range (mean − 2·SD, with a 20% / 20× floor).
Male Infertility → Normal Male, Female Infertility → Normal Female, AF/POC →
pooled Normal. Regions not covered in the reference (e.g. chrY in females) show
**n/a** instead of a false flag.

> **Note:** BAM files and the reference data are **not** in this repository (too
> large, and private). The published GitHub Pages site is an informational landing
> page only — the tool must run on a server that has the BED panels and BAMs.

The dir → label/slug mapping is in `samples.py` (`DIR_MAP`). For each query, every
replicate's depth is computed in parallel (process pool) and the replicates of a
type are aggregated to one row: **mean of replicate means (± SD)**, mean median,
worst-case min, and mean %≥threshold. Replicate count shows as `n=…`.

**To add samples:** drop more indexed `*.bam` into the matching type folder and
`./restart.sh`. To add a new type, create a folder and add a `DIR_MAP` row.

Sample read-depth is **mandatory** (always computed for all types) so every report
shows both BED coverage and real sample coverage.

## API

- `GET /api/query?q=<term>&panel=<name>` → JSON report (BED + sample depth)
- `GET /api/bed?q=<term>&panel=<name>`   → matched intervals as BED
- `GET /api/panels` · `GET /api/samples`
- `POST /api/upload` (multipart `bed`) → register an uploaded panel

## Security / deployment notes

This is an internal lab tool. For a trusted LAN it runs as-is. For shared or
untrusted networks:

- `HOST=127.0.0.1 ./start.sh` — bind to localhost only (default is `0.0.0.0`
  for LAN access) and front it with an authenticating reverse proxy.
- `COVERAGE_UPLOAD_TOKEN=<secret>` — require header `X-Upload-Token` on the
  `/api/upload` write endpoint (unset = open).
- `MAX_UPLOAD_MB=512` — cap upload size.

rs-ID input is validated (`rs\d+`) and URL-encoded before the Ensembl call; the
rs-ID cache is size-bounded; API errors are logged server-side and return a
generic message. The dev server is fine for internal use; for heavier load run
behind gunicorn/uwsgi.

## Files

- `app.py` — Flask server + single-page UI (Anderson-branded)
- `coverage_index.py` — in-memory BED index (bisect overlap + token map)
- `samples.py` — reference sample-type registry (directory groups) + pysam depth engine
- `static/anderson.png` — brand logo (header + PDF)
- `start.sh` / `stop.sh` / `status.sh` / `restart.sh` / `server.conf`
