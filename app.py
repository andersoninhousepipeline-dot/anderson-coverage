"""
Coverage Checker — web application.

Query a capture panel (annotated BED) for coverage of a:
  • Gene name        e.g.  BRCA1
  • Chr region       e.g.  chr17:43044295-43125483   (or chr17:43044295)
  • NM / transcript  e.g.  NM_007294.4   (also ENST / ENSG / NR / XM / XR)
  • rs ID            e.g.  rs6265        (resolved to coordinates via Ensembl REST)

Run:   python3 app.py        then open http://localhost:8080
"""

import os
import re
import time
import datetime
import requests
from urllib.parse import quote
from flask import Flask, request, jsonify, render_template_string, Response

from coverage_index import Panel, merge_intervals, tokenize, _STOP_TOKENS
import samples as samplemod

BED_DIR = os.environ.get("BED_DIR", "/data/bed")
UPLOAD_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "uploads")
os.makedirs(UPLOAD_DIR, exist_ok=True)

# Friendly names for the two primary panels; everything else is auto-discovered.
PREFERRED = {
    "hg38_exome_comp_spikein_v2.0.2_targets_sorted.re_annotated.bed":
        "Twist Spikein v2.0.2 (re-annotated)",
    "hg38_exome_comp_spikein_v2.0.2_with_mito_merged.bed":
        "Twist Spikein v2.0.2 + mito (merged)",
}

PANELS = {}        # name -> Panel (lazy loaded)
PANEL_ORDER = []   # display order


def register_panel(name, path):
    if name in PANELS:
        return
    PANELS[name] = Panel(name, path)
    PANEL_ORDER.append(name)


def discover_panels():
    """Scan BED_DIR (and uploads/) for *.bed files and register them as panels."""
    PANELS.clear()
    PANEL_ORDER.clear()
    found = []
    for root, _dirs, files in os.walk(BED_DIR):
        # don't descend into the app's own folder twice / skip caches
        if os.path.basename(root) in ("__pycache__",):
            continue
        for fn in files:
            if fn.endswith(".bed"):
                found.append((root, fn))
    # preferred files first, in defined order
    pref_paths = []
    for base, label in PREFERRED.items():
        p = os.path.join(BED_DIR, base)
        if os.path.exists(p):
            register_panel(label, p)
            pref_paths.append(p)
    # the rest, alphabetically by relative path
    rest = []
    for root, fn in found:
        full = os.path.join(root, fn)
        if full in pref_paths:
            continue
        rel = os.path.relpath(full, BED_DIR)
        rest.append((rel, full))
    for rel, full in sorted(rest):
        register_panel(rel, full)
    # uploaded panels
    for fn in sorted(os.listdir(UPLOAD_DIR)):
        if fn.endswith(".bed"):
            register_panel("uploads/" + fn, os.path.join(UPLOAD_DIR, fn))


discover_panels()

ENSEMBL = "https://rest.ensembl.org"
GENOME_BUILD = "GRCh38 / hg38"

app = Flask(__name__)
# Cap upload size (panels can be ~90 MB; default 512 MB).
app.config["MAX_CONTENT_LENGTH"] = int(os.environ.get("MAX_UPLOAD_MB", "512")) * 1024 * 1024
# Optional shared secret to gate the write endpoint on untrusted networks.
# When unset (default), upload is open — suitable for a trusted LAN deployment.
UPLOAD_TOKEN = os.environ.get("COVERAGE_UPLOAD_TOKEN", "")


# ── panel access ──────────────────────────────────────────────────────────────
def get_panel(name=None):
    if not PANEL_ORDER:
        raise RuntimeError("No panel BED files found in %s" % BED_DIR)
    if not name or name not in PANELS:
        name = PANEL_ORDER[0]
    return PANELS[name].load()


# ── rs-id resolution (Ensembl REST) ───────────────────────────────────────────
_rsid_cache = {}
_RSID_CACHE_MAX = 2000

def resolve_rsid(rsid):
    rsid = rsid.strip()
    # Validate before it reaches the outbound HTTP request / cache key.
    if not re.match(r"^rs\d+$", rsid, re.I):
        return {"error": "Invalid rs ID (expected e.g. rs6265)."}
    if rsid in _rsid_cache:
        return _rsid_cache[rsid]
    # Bound the cache so an attacker cannot grow it without limit.
    if len(_rsid_cache) >= _RSID_CACHE_MAX:
        _rsid_cache.clear()
    url = f"{ENSEMBL}/variation/human/{quote(rsid, safe='')}?content-type=application/json"
    try:
        r = requests.get(url, timeout=15)
        if r.status_code != 200:
            res = {"error": f"Ensembl returned HTTP {r.status_code} for {rsid}"}
            _rsid_cache[rsid] = res
            return res
        data = r.json()
    except Exception as e:
        return {"error": f"Lookup failed for {rsid}: {e}"}
    mappings = []
    for m in data.get("mappings", []):
        loc = m.get("location", "")
        # location like "11:27654893-27654893" on GRCh38
        if m.get("assembly_name") and m["assembly_name"] != "GRCh38":
            continue
        mm = re.match(r"^([^:]+):(\d+)-(\d+)$", loc)
        if not mm:
            continue
        chrom = mm.group(1)
        if not chrom.startswith("chr"):
            chrom = "chr" + chrom
        mappings.append({
            "chrom": chrom,
            "start": int(mm.group(2)),
            "end": int(mm.group(3)),
            "allele": m.get("allele_string", ""),
        })
    res = {
        "rsid": rsid,
        "mappings": mappings,
        "consequence": data.get("most_severe_consequence", ""),
        "synonyms": data.get("synonyms", []),
    }
    _rsid_cache[rsid] = res
    return res


# ── annotation parsing for display ─────────────────────────────────────────────
def split_annot(annot):
    genes, ids, other = [], [], []
    seen_g, seen_i = set(), set()
    for tok in tokenize(annot):
        if tok.startswith(("NM_", "NR_", "XM_", "XR_", "ENST", "ENSG", "CCDS")):
            if tok not in seen_i:
                ids.append(tok); seen_i.add(tok)
        elif tok.startswith("ClinID"):
            other.append(tok)
        elif tok.isdigit() or tok.lower() in _STOP_TOKENS:
            continue
        else:
            if tok not in seen_g:
                genes.append(tok); seen_g.add(tok)
    return genes, ids, other


def rows_payload(panel, indices):
    rows = []
    for i in indices:
        chrom, start, end, annot = panel.row(i)
        genes, ids, _ = split_annot(annot)
        rows.append({
            "chrom": chrom, "start": start, "end": end,
            "length": end - start,
            "genes": genes, "ids": ids, "annot": annot,
        })
    rows.sort(key=lambda r: (r["chrom"], r["start"]))
    return rows


def summarize_rows(rows):
    if not rows:
        return {"intervals": 0, "targeted_bp": 0, "chroms": [],
                "span_start": None, "span_end": None}
    ivs = merge_intervals([(r["start"], r["end"]) for r in rows
                           if r["chrom"] == rows[0]["chrom"]])
    targeted = sum(e - s for s, e in ivs)
    chroms = sorted({r["chrom"] for r in rows})
    span_start = min(r["start"] for r in rows)
    span_end = max(r["end"] for r in rows)
    return {
        "intervals": len(rows),
        "targeted_bp": targeted,
        "chroms": chroms,
        "span_start": span_start,
        "span_end": span_end,
    }


# ── query-type detection ───────────────────────────────────────────────────────
RE_REGION = re.compile(r"^(chr[\w]+)[:\s]+([\d,]+)(?:[-\s]+([\d,]+))?$", re.I)
RE_RSID = re.compile(r"^rs\d+$", re.I)
RE_NM = re.compile(r"^(NM_|NR_|XM_|XR_|ENST|ENSG|CCDS)", re.I)

def detect_type(q):
    q = q.strip()
    if RE_RSID.match(q):
        return "rsid"
    if RE_REGION.match(q):
        return "region"
    # multiple genes/transcripts separated by comma / whitespace / newline
    tokens = [t for t in re.split(r"[,\s]+", q) if t]
    if len(tokens) > 1:
        return "genes"
    if RE_NM.match(q):
        return "transcript"
    return "gene"


# ── the query engine ───────────────────────────────────────────────────────────
def do_query(q, qtype, panel_name, sample_accs=None):
    panel = get_panel(panel_name)
    q = q.strip()
    sample_accs = sample_accs or []
    if qtype == "auto":
        qtype = detect_type(q)

    # A multi-token gene/transcript list -> multi-gene report.
    tokens = [t for t in re.split(r"[,\s]+", q) if t]
    if qtype in ("gene", "transcript", "genes") and len(tokens) > 1:
        qtype = "genes"

    result = {"query": q, "type": qtype, "panel": panel.name,
              "build": GENOME_BUILD, "ok": True, "messages": [],
              "samples_requested": sample_accs}

    if qtype == "genes":
        names = []
        seen = set()
        for t in tokens:                       # de-dupe, preserve order
            if t.upper() not in seen:
                seen.add(t.upper()); names.append(t)
        genes = []
        all_ivs = []
        for name in names:
            idxs = panel.rows_for_token(name)
            if not idxs:
                idxs = panel.free_text(name)
            rows = rows_payload(panel, idxs)
            ivs = [(r["chrom"], r["start"], r["end"]) for r in rows]
            all_ivs.extend(ivs)
            genes.append({
                "name": name,
                "found": bool(rows),
                "summary": summarize_rows(rows),
                "sample_coverage": compute_sample_coverage(ivs, sample_accs),
            })
        result.update({
            "genes": genes,
            "n_genes": len(names),
            "n_found": sum(1 for g in genes if g["found"]),
            "combined_summary": {
                "targeted_bp": sum(
                    e - s for ivs in samplemod._merge(all_ivs).values()
                    for s, e in ivs),
                "intervals": len(all_ivs),
            },
            "combined_sample_coverage": compute_sample_coverage(all_ivs, sample_accs),
        })
        return result

    if qtype == "region":
        m = RE_REGION.match(q)
        if not m:
            return {"ok": False, "error": "Could not parse region. Use chr1:65564-70008."}
        chrom = m.group(1)
        if not chrom.startswith("chr"):
            chrom = "chr" + chrom
        chrom = "chr" + chrom[3:]  # normalise case of 'chr'
        start = int(m.group(2).replace(",", ""))
        end = int(m.group(3).replace(",", "")) if m.group(3) else start + 1
        if end < start:
            start, end = end, start
        idxs = panel.overlap(chrom, start, end)
        rows = rows_payload(panel, idxs)
        # coverage of the queried window
        clipped = merge_intervals([(max(r["start"], start), min(r["end"], end))
                                   for r in rows])
        covered = sum(e - s for s, e in clipped)
        region_len = max(1, end - start)
        gaps = compute_gaps(start, end, clipped)
        # depth over captured portion of the region (fall back to whole region if no targets)
        depth_ivs = [(chrom, s, e) for s, e in clipped] or [(chrom, start, end)]
        result.update({
            "region": {"chrom": chrom, "start": start, "end": end,
                       "length": end - start},
            "rows": rows,
            "coverage": {
                "covered_bp": covered,
                "region_bp": end - start,
                "pct": round(100.0 * covered / region_len, 2),
                "n_covered_segments": len(clipped),
                "gaps": gaps,
                "fully_covered": covered >= (end - start) and (end - start) > 0,
            },
            "summary": summarize_rows(rows),
            "sample_coverage": compute_sample_coverage(depth_ivs, sample_accs),
        })
        return result

    if qtype == "rsid":
        info = resolve_rsid(q)
        if "error" in info:
            return {"ok": False, "error": info["error"]}
        hits = []
        for mp in info["mappings"]:
            idxs = panel.overlap(mp["chrom"], mp["start"] - 1, mp["end"])
            rows = rows_payload(panel, idxs)
            nearest = None
            if not rows:
                nearest = nearest_interval(panel, mp["chrom"], mp["start"])
            # read depth AT the variant base for each sample
            base_iv = [(mp["chrom"], mp["start"] - 1, mp["end"])]
            hits.append({
                "position": mp, "covered": bool(rows),
                "rows": rows, "nearest": nearest,
                "sample_coverage": compute_sample_coverage(base_iv, sample_accs),
            })
        result.update({
            "rsid": {"id": q, "consequence": info.get("consequence", ""),
                     "synonyms": info.get("synonyms", [])[:6]},
            "hits": hits,
            "covered": any(h["covered"] for h in hits),
        })
        return result

    # gene or transcript (token lookup), with free-text fallback
    idxs = panel.rows_for_token(q)
    used_fallback = False
    if not idxs:
        idxs = panel.free_text(q)
        used_fallback = bool(idxs)
        if used_fallback:
            result["messages"].append(
                "No exact gene/transcript token matched — showing substring matches.")
    rows = rows_payload(panel, idxs)
    depth_ivs = [(r["chrom"], r["start"], r["end"]) for r in rows]
    result.update({
        "rows": rows,
        "summary": summarize_rows(rows),
        "fallback": used_fallback,
        "found": bool(rows),
        "sample_coverage": compute_sample_coverage(depth_ivs, sample_accs),
    })
    return result


def compute_sample_coverage(intervals, sample_accs):
    """Run per-sample depth over the given intervals; returns list or None."""
    if not sample_accs or not intervals:
        return None
    accs = [a for a in sample_accs if a in samplemod.SAMPLES]
    if not accs:
        return None
    try:
        return samplemod.coverage_multi(accs, intervals)
    except Exception as e:
        return [{"error": f"Depth computation failed: {e}"}]


def compute_gaps(start, end, covered_segments):
    gaps = []
    cur = start
    for s, e in covered_segments:
        if s > cur:
            gaps.append({"start": cur, "end": s, "length": s - cur})
        cur = max(cur, e)
    if cur < end:
        gaps.append({"start": cur, "end": end, "length": end - cur})
    return gaps


def nearest_interval(panel, chrom, pos):
    c = panel.by_chrom.get(chrom)
    if not c:
        return None
    best = None
    best_d = None
    for s, e, i in zip(c["starts"], c["ends"], c["idx"]):
        if e <= pos:
            d = pos - e
        elif s > pos:
            d = s - pos
        else:
            d = 0
        if best_d is None or d < best_d:
            best_d = d
            best = i
        if s > pos and (best_d is not None and (s - pos) > best_d):
            # starts are sorted; once we pass pos and distance grows we can stop
            if s - pos > best_d:
                break
    if best is None:
        return None
    chrom, s, e, annot = panel.row(best)
    genes, ids, _ = split_annot(annot)
    return {"chrom": chrom, "start": s, "end": e,
            "distance": best_d, "genes": genes, "ids": ids}


# ── routes ─────────────────────────────────────────────────────────────────────
@app.route("/")
def index():
    return render_template_string(PAGE, panels=PANEL_ORDER, build=GENOME_BUILD)


@app.route("/api/upload", methods=["POST"])
def api_upload():
    # Optional auth: when COVERAGE_UPLOAD_TOKEN is set, require it on writes.
    if UPLOAD_TOKEN and request.headers.get("X-Upload-Token", "") != UPLOAD_TOKEN:
        return jsonify({"ok": False, "error": "Upload not authorized."}), 403
    f = request.files.get("bed")
    if not f or not f.filename:
        return jsonify({"ok": False, "error": "No file received."})
    name = re.sub(r"[^A-Za-z0-9._-]+", "_", os.path.basename(f.filename))
    if not name.endswith(".bed"):
        name += ".bed"
    dest = os.path.join(UPLOAD_DIR, name)
    f.save(dest)
    # quick validation: at least one parseable BED line
    ok_line = False
    with open(dest) as fh:
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) >= 3 and p[1].isdigit() and p[2].isdigit():
                ok_line = True
                break
    if not ok_line:
        os.remove(dest)
        return jsonify({"ok": False, "error": "File is not a tab-delimited BED "
                        "(need chrom<TAB>start<TAB>end)."})
    pname = "uploads/" + name
    # (re)register and load
    if pname in PANELS:
        del PANELS[pname]
        PANEL_ORDER.remove(pname)
    register_panel(pname, dest)
    PANELS[pname].load()
    return jsonify({"ok": True, "panel": pname,
                    "rows": PANELS[pname].n_rows,
                    "targeted_bp": PANELS[pname].total_bp})


@app.route("/api/samples")
def api_samples():
    return jsonify({"samples": samplemod.list_samples(),
                    "thresholds": samplemod.THRESHOLDS})


@app.route("/api/panels")
def api_panels():
    out = []
    for name in PANEL_ORDER:
        p = PANELS[name]
        out.append({"name": name, "loaded": p.loaded,
                    "rows": p.n_rows, "targeted_bp": p.total_bp})
    return jsonify(out)


@app.route("/api/query")
def api_query():
    q = request.args.get("q", "").strip()
    qtype = request.args.get("type", "auto")
    panel_name = request.args.get("panel", "")
    # Sample read-depth coverage is mandatory — always run all reference samples.
    sample_accs = list(samplemod.SAMPLE_ORDER)
    if not q:
        return jsonify({"ok": False, "error": "Empty query."})
    try:
        res = do_query(q, qtype, panel_name, sample_accs)
    except Exception:
        app.logger.exception("query failed for q=%r type=%r panel=%r", q, qtype, panel_name)
        return jsonify({"ok": False, "error": "Query failed; see server log."})
    return jsonify(res)


@app.route("/api/bed")
def api_bed():
    """Download matched intervals as a BED file."""
    q = request.args.get("q", "").strip()
    qtype = request.args.get("type", "auto")
    panel_name = request.args.get("panel", "")
    try:
        res = do_query(q, qtype, panel_name)
        rows = []
        if res.get("rows"):
            rows = res["rows"]
        elif res.get("hits"):
            for h in res["hits"]:
                rows.extend(h.get("rows", []))
        elif res.get("type") == "genes":
            panel = get_panel(panel_name)
            for tok in [t for t in re.split(r"[,\s]+", q) if t]:
                idxs = panel.rows_for_token(tok) or panel.free_text(tok)
                rows.extend(rows_payload(panel, idxs))
    except Exception:
        app.logger.exception("bed export failed for q=%r", q)
        return Response("# export failed; see server log\n", mimetype="text/plain")
    lines = ["\t".join([r["chrom"], str(r["start"]), str(r["end"]), r["annot"]])
             for r in rows]
    body = "\n".join(lines) + ("\n" if lines else "")
    fn = re.sub(r"[^A-Za-z0-9]+", "_", q) or "coverage"
    return Response(body, mimetype="text/plain",
                    headers={"Content-Disposition": f"attachment; filename={fn}.bed"})


# ── frontend (single page) ─────────────────────────────────────────────────────
PAGE = r"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1">
<title>Coverage Checker</title>
<style>
  * { box-sizing: border-box; }
  body { font-family: Arial, Helvetica, sans-serif; font-size: 13px; color: #222;
         background: #f5f6f8; margin: 0; }
  .wrap { max-width: 1080px; margin: 0 auto; padding: 0 18px 60px; }
  .report-header { background: #fff; border-bottom: 3px solid #F08020;
        box-shadow: 0 1px 4px rgba(0,0,0,.06); }
  .hdr-inner { display: flex; align-items: center; gap: 22px; padding: 14px 18px; }
  .brand-logo { height: 52px; width: auto; }
  .hdr-text { border-left: 2px solid #e2e6ec; padding-left: 22px; }
  .report-header h1 { font-size: 19px; margin: 0; font-weight: bold; color: #005FA0; }
  .report-header .meta { font-size: 11px; color: #6a7180; margin-top: 6px;
        display: flex; flex-wrap: wrap; gap: 18px; }
  .searchbar { background: #fff; border: 1px solid #d0d5de; border-radius: 8px;
        padding: 16px 18px; margin: 18px 0; box-shadow: 0 1px 3px rgba(0,0,0,.05); }
  .searchbar form { display: flex; gap: 10px; flex-wrap: wrap; align-items: center; }
  input[type=text] { flex: 1 1 320px; padding: 10px 12px; font-size: 14px;
        border: 1px solid #c4ccd8; border-radius: 6px; }
  select { padding: 9px 8px; border: 1px solid #c4ccd8; border-radius: 6px; font-size: 13px; }
  button { background: #005FA0; color: #fff; border: 0; padding: 10px 20px;
        font-size: 14px; border-radius: 6px; cursor: pointer; }
  button:hover { background: #00497d; }
  button.ghost { background: #fff; color: #005FA0; border: 1px solid #005FA0; padding: 7px 14px; font-size: 12px; }
  .ghostlbl { background: #fff; color: #005FA0; border: 1px solid #005FA0; padding: 6px 12px;
        font-size: 12px; border-radius: 6px; cursor: pointer; font-weight: bold; }
  .ghostlbl:hover { background: #eef2fa; }
  .samplebox { margin-top: 14px; border-top: 1px dashed #d0d5de; padding-top: 12px; }
  .samplebox-h { font-size: 11px; color: #005FA0; font-weight: bold; margin-bottom: 8px; }
  .samplelist { display: flex; flex-wrap: wrap; gap: 8px; }
  .schip { display: inline-flex; align-items: center; gap: 6px; border: 1px solid #c4ccd8;
        border-radius: 99px; padding: 5px 12px; font-size: 12px; cursor: pointer; background:#fff; }
  .schip input { cursor: pointer; }
  .schip.on { background: #eaf3ec; border-color: #2a7a47; color: #1a6033; font-weight: bold; }
  .depthwrap { background:#fff; border:1px solid #d0d5de; border-radius:8px; overflow:auto; margin-top:12px; }
  td.d-hi { color:#1a6033; font-weight:bold; } td.d-mid { color:#9a6a00; font-weight:bold; }
  td.d-lo { color:#8b1a1a; font-weight:bold; }
  .legend { font-size:10px; color:#777; margin-top:6px; }
  tr.combined td { background:#fff4e8; border-top:2px solid #F08020; font-weight:bold; }
  .hint { font-size: 11px; color: #6a7180; margin-top: 8px; }
  .hint code { background: #eef1f6; padding: 1px 5px; border-radius: 3px; cursor: pointer; }
  h2 { font-size: 11px; font-weight: bold; text-transform: uppercase; letter-spacing: .6px;
       color: #005FA0; border-bottom: 1.5px solid #F08020; padding-bottom: 4px; margin: 24px 0 12px; }
  .summary-bar { display: flex; flex-wrap: wrap; border: 1px solid #d0d5de;
       border-radius: 8px; background: #fff; overflow: hidden; }
  .s-cell { flex: 1 1 130px; padding: 14px 16px; border-right: 1px solid #e4e8f0; }
  .s-cell .val { font-size: 22px; font-weight: bold; }
  .s-cell .lbl { font-size: 10px; color: #666; margin-top: 2px; text-transform: uppercase; letter-spacing:.4px; }
  .val-total { color: #005FA0; } .val-green { color: #1a6033; } .val-red { color: #8b1a1a; }
  .prog-wrap { background: #fff; border: 1px solid #d0d5de; border-radius: 8px;
       padding: 14px 18px; margin-top: 14px; }
  .prog-outer { height: 10px; background: #e4e8f0; border-radius: 99px; overflow: hidden; }
  .prog-fill { height: 100%; background: #2a7a47; border-radius: 99px; transition: width .4s; }
  .table-wrap { background: #fff; border: 1px solid #d0d5de; border-radius: 8px;
       overflow: auto; margin-top: 12px; }
  table { border-collapse: collapse; width: 100%; }
  thead tr { background: #005FA0; }
  thead th { color: #fff; padding: 8px 10px; font-size: 10px; font-weight: bold;
       text-align: left; white-space: nowrap; border-right: 1px solid #1f6fae; }
  tbody tr:nth-child(even) { background: #f8f9fb; }
  tbody tr:hover { background: #eef2fa; }
  tbody td { padding: 7px 10px; font-size: 11px; vertical-align: top; border-right: 1px solid #e4e8f0; }
  .gene { font-weight: bold; color: #005FA0; }
  .mono { font-family: monospace; font-size: 10.5px; color: #444; }
  .cov-yes { font-weight: bold; color: #1a6033; }
  .cov-no  { font-weight: bold; color: #8b1a1a; }
  .pill { display: inline-block; padding: 2px 8px; border-radius: 99px; font-size: 11px; font-weight: bold; }
  .pill.yes { background: #e3f3e9; color: #1a6033; }
  .pill.no  { background: #fbe6e6; color: #8b1a1a; }
  .msg { background: #fff7e6; border: 1px solid #f0d9a8; color: #7a5a12;
       padding: 10px 14px; border-radius: 6px; font-size: 12px; margin-top: 12px; }
  .err { background: #fbe6e6; border: 1px solid #e3b5b5; color: #8b1a1a;
       padding: 12px 16px; border-radius: 6px; margin-top: 16px; }
  .spinner { display:inline-block; width:16px; height:16px; border:3px solid #cdd5e2;
       border-top-color:#005FA0; border-radius:50%; animation:spin .7s linear infinite; vertical-align:middle;}
  @keyframes spin { to { transform: rotate(360deg); } }
  .empty { color:#6a7180; font-size:13px; padding:14px 0; }
  .toolbar { display:flex; gap:8px; margin-top:10px; flex-wrap:wrap; }
  @media print {
    .searchbar, .toolbar { display:none; }
    body { background:#fff; -webkit-print-color-adjust: exact; print-color-adjust: exact; }
    .report-header { border-bottom:3px solid #F08020 !important; }
    thead th, tr.combined td, .prog-fill, .schip.on { -webkit-print-color-adjust: exact; print-color-adjust: exact; }
    h2 { page-break-after: avoid; }
    .depthwrap, .table-wrap { page-break-inside: auto; }
  }
</style>
</head>
<body>
<div class="report-header">
  <div class="wrap hdr-inner">
    <img src="/static/anderson.png" class="brand-logo" alt="Anderson Diagnostics &amp; Labs">
    <div class="hdr-text">
      <h1>Coverage Checker</h1>
      <div class="meta">
        <span>Genome build: {{ build }}</span>
        <span>Panel target &amp; sample read-depth coverage</span>
        <span id="hdr-date"></span>
      </div>
    </div>
  </div>
</div>

<div class="wrap">
  <div class="searchbar">
    <form id="qform" autocomplete="off">
      <input type="text" id="q" placeholder="Gene(s): BRCA1 or BRCA1, BRCA2, TP53 · region chr17:43044295-43125483 · NM_007294.4 · rs6265" autofocus>
      <select id="qtype">
        <option value="auto">Auto-detect</option>
        <option value="gene">Gene name</option>
        <option value="region">Chr region</option>
        <option value="transcript">NM / transcript ID</option>
        <option value="rsid">rs ID</option>
      </select>
      <select id="panel" title="Choose which BED panel to check against"></select>
      <button type="submit">Check coverage</button>
    </form>
    <div class="hint">Examples:
      <code class="ex">BRCA1</code>
      <code class="ex">BRCA1, BRCA2, TP53</code>
      <code class="ex">chr17:43044295-43125483</code>
      <code class="ex">NM_007294.4</code>
      <code class="ex">rs6265</code>
    </div>
    <div class="hint" style="margin-top:10px; display:flex; gap:8px; align-items:center; flex-wrap:wrap">
      <span><b>Panel (BED file):</b> switch the dropdown above to check the same query against a different panel.</span>
      <label class="ghostlbl">Upload BED…
        <input type="file" id="bedfile" accept=".bed,.txt" style="display:none">
      </label>
      <span id="uploadmsg"></span>
    </div>
    <div class="samplebox">
      <div class="samplebox-h">Sample read-depth coverage is included automatically with every search:</div>
      <div id="samplelist" class="samplelist"><span style="color:#888">loading…</span></div>
    </div>
  </div>
  <div id="out"></div>
</div>

<script>
const $ = s => document.querySelector(s);
document.getElementById('hdr-date').textContent = 'Generated: ' + new Date().toLocaleString();
const out = $('#out');

document.querySelectorAll('.ex').forEach(c =>
  c.addEventListener('click', () => { $('#q').value = c.textContent; $('#qform').requestSubmit(); }));

$('#qform').addEventListener('submit', e => { e.preventDefault(); run(); });

async function loadPanels(selectName){
  try {
    const r = await fetch('/api/panels');
    const list = await r.json();
    const sel = $('#panel');
    sel.innerHTML = '';
    for(const p of list){
      const o = document.createElement('option');
      o.value = p.name;
      o.textContent = p.name + (p.loaded ? '  ('+p.rows.toLocaleString()+' intervals)' : '');
      sel.appendChild(o);
    }
    if(selectName) sel.value = selectName;
  } catch(e){ /* ignore */ }
}
loadPanels();

let THRESHOLDS = [1,10,20,30,50,100];
async function loadSamples(){
  try {
    const r = await fetch('/api/samples');
    const d = await r.json();
    THRESHOLDS = d.thresholds || THRESHOLDS;
    const box = $('#samplelist');
    if(!d.samples || !d.samples.length){ box.innerHTML = '<span style="color:#888">No reference samples configured.</span>'; return; }
    box.innerHTML = d.samples.map(s => '<span class="schip on">'+esc(s.label)+
      (s.n?' <span style="opacity:.7;font-weight:normal">n='+s.n+'</span>':'')+'</span>').join('');
  } catch(e){ $('#samplelist').innerHTML = '<span style="color:#888">samples unavailable</span>'; }
}
loadSamples();

$('#bedfile').addEventListener('change', async e => {
  const f = e.target.files[0];
  if(!f) return;
  const msg = $('#uploadmsg');
  msg.innerHTML = '<span class="spinner"></span> uploading & indexing…';
  const fd = new FormData(); fd.append('bed', f);
  try {
    const r = await fetch('/api/upload', {method:'POST', body: fd});
    const d = await r.json();
    if(d.ok){
      msg.textContent = '✓ loaded ' + d.panel + ' (' + d.rows.toLocaleString() + ' intervals)';
      await loadPanels(d.panel);
    } else {
      msg.innerHTML = '<span style="color:#8b1a1a">✗ ' + esc(d.error) + '</span>';
    }
  } catch(err){
    msg.innerHTML = '<span style="color:#8b1a1a">✗ ' + esc(err.message) + '</span>';
  }
  e.target.value = '';
});

function fmt(n){ return (n==null)?'':n.toLocaleString(); }
function esc(s){ return (s||'').replace(/[&<>]/g, c => ({'&':'&amp;','<':'&lt;','>':'&gt;'}[c])); }
function reg(r){ return esc(r.chrom)+':'+fmt(r.start)+'-'+fmt(r.end); }

function rowsTable(rows){
  if(!rows || !rows.length) return '<div class="empty">No target intervals.</div>';
  let h = '<div class="table-wrap"><table><thead><tr>'+
    '<th>Chrom</th><th>Start</th><th>End</th><th>Length (bp)</th><th>Gene(s)</th><th>Transcript / IDs</th>'+
    '</tr></thead><tbody>';
  for(const r of rows){
    h += '<tr><td class="mono">'+esc(r.chrom)+'</td><td class="mono">'+fmt(r.start)+
      '</td><td class="mono">'+fmt(r.end)+'</td><td class="mono">'+fmt(r.length)+
      '</td><td class="gene">'+esc((r.genes||[]).join(', '))+
      '</td><td class="mono">'+esc((r.ids||[]).join(', '))+'</td></tr>';
  }
  return h+'</tbody></table></div>';
}

function dcell(pct){
  if(pct==null) return '<td>—</td>';
  const cls = pct>=95?'d-hi':(pct>=80?'d-mid':'d-lo');
  return '<td class="'+cls+'">'+pct+'%</td>';
}
function depthTable(cov, opts){
  opts = opts || {};
  if(!cov || !cov.length) return '';
  let h = '<h2>'+(opts.title||'Sample read-depth coverage')+'</h2><div class="depthwrap"><table><thead><tr>'+
    '<th>Sample</th><th>Target bp</th><th>Mean depth</th><th>Median</th><th>Min</th>';
  for(const t of THRESHOLDS) h += '<th>% &ge; '+t+'x</th>';
  h += '</tr></thead><tbody>';
  for(const c of cov){
    if(c.error){
      h += '<tr><td>'+esc(c.label||c.accession||'sample')+'</td><td colspan="'+(4+THRESHOLDS.length)+
           '" class="d-lo">'+esc(c.error)+'</td></tr>';
      continue;
    }
    h += '<tr><td class="gene">'+esc(c.label)+'</td>'+
      '<td class="mono">'+fmt(c.total_bases)+'</td>'+
      '<td class="mono"><b>'+c.mean+'x</b></td>'+
      '<td class="mono">'+c.median+'x</td>'+
      '<td class="mono">'+c.min+'x</td>';
    for(const t of THRESHOLDS) h += dcell(c.pct[String(t)]);
    h += '</tr>';
  }
  // Combined row: average across all valid samples
  const valid = cov.filter(c => !c.error && c.total_bases);
  if(valid.length > 1){
    const avg = a => Math.round(a.reduce((x,y)=>x+y,0)/a.length*10)/10;
    const tb = valid[0].total_bases;
    h += '<tr class="combined"><td class="gene">Overall Coverage</td>'+
      '<td class="mono">'+fmt(tb)+'</td>'+
      '<td class="mono"><b>'+avg(valid.map(c=>c.mean))+'x</b></td>'+
      '<td class="mono">'+avg(valid.map(c=>c.median))+'x</td>'+
      '<td class="mono">'+avg(valid.map(c=>c.min))+'x</td>';
    for(const t of THRESHOLDS) h += dcell(avg(valid.map(c=>c.pct[String(t)]||0)));
    h += '</tr>';
  }
  h += '</tbody></table></div><div class="legend">% of target bases at or above each depth. '+
       'Green &ge;95% · amber &ge;80% · red &lt;80%. Reads: duplicates/secondary/QC-fail excluded. '+
       '<b>Overall Coverage</b> = average across all reference samples.</div>';
  return h;
}
// Case-vs-reference: flag case sample types whose mean depth falls below the
// sex-matched Normal reference range (Normal mean - 2*SD, with a 20% floor).
const CASE_TYPES = [
  {id:'male-infertility',   ref:'normal-male'},
  {id:'female-infertility', ref:'normal-female'},
  {id:'af',                 ref:'pooled'},
  {id:'poc',                ref:'pooled'},
];
const ABS_FLOOR = 20;   // mean depth below this is always flagged
function referencePanel(cov){
  if(!cov || !cov.length) return '';
  const byId = {}; cov.forEach(c => { if(!c.error) byId[c.id] = c; });
  const nm = byId['normal-male'], nf = byId['normal-female'];
  if(!nm && !nf) return '';
  let pooled = null;
  if(nm && nf){ const m=(nm.mean+nf.mean)/2, sd=Math.max(nm.mean_sd||0,nf.mean_sd||0);
    pooled = {mean:m, sd:sd, label:'Normal (M+F)'}; }
  const refOf = key => key==='pooled' ? pooled
        : (key==='normal-male' ? (nm&&{mean:nm.mean,sd:nm.mean_sd||0,label:nm.label})
                               : (nf&&{mean:nf.mean,sd:nf.mean_sd||0,label:nf.label}));
  const rows = [];
  for(const ct of CASE_TYPES){
    const c = byId[ct.id]; const r = refOf(ct.ref);
    if(!c || !r) continue;
    const low = Math.min(r.mean - 2*r.sd, 0.8*r.mean);
    let status;                                  // 'na' | 'flag' | 'ok'
    if(r.mean < ABS_FLOOR) status = 'na';        // region not covered in reference (e.g. chrY in females)
    else status = ((c.mean < low) || (c.mean < ABS_FLOOR)) ? 'flag' : 'ok';
    rows.push({label:c.label, mean:c.mean, refLabel:r.label, refMean:Math.round(r.mean),
               low:Math.max(0,Math.round(low)), status});
  }
  if(!rows.length) return '';
  const anyFlag = rows.some(r => r.status==='flag');
  let h = '<h2>Case vs reference (Normal) range '+
    (anyFlag?'<span class="pill no">⚠ '+rows.filter(r=>r.status==='flag').length+' flagged</span>'
            :'<span class="pill yes">all within range</span>')+'</h2>'+
    '<div class="depthwrap"><table><thead><tr>'+
    '<th>Case type</th><th>Mean depth</th><th>Reference</th><th>Expected (mean − 2 SD)</th><th>Status</th>'+
    '</tr></thead><tbody>';
  for(const r of rows){
    const cls = r.status==='flag'?'d-lo':(r.status==='ok'?'d-hi':'');
    const badge = r.status==='flag' ? '<span class="pill no">⚠ BELOW RANGE</span>'
                : r.status==='ok'   ? '<span class="pill yes">within range</span>'
                :                     '<span class="empty">n/a — not covered in reference</span>';
    h += '<tr><td class="gene">'+esc(r.label)+'</td>'+
      '<td class="mono '+cls+'">'+r.mean+'x</td>'+
      '<td class="mono">'+esc(r.refLabel)+' ('+r.refMean+'x)</td>'+
      '<td class="mono">'+(r.status==='na'?'—':'&ge; '+r.low+'x')+'</td>'+
      '<td>'+badge+'</td></tr>';
  }
  h += '</tbody></table></div><div class="legend">Flagged when a case type’s mean depth '+
       'is below its sex-matched Normal reference (mean − 2·SD, with a 20% / '+ABS_FLOOR+'× floor). '+
       'AF / POC compared to pooled Normal (M+F). <b>n/a</b> = region not covered in the reference '+
       '(e.g. chrY in females).</div>';
  return h;
}

function depthAtVariant(cov){
  if(!cov || !cov.length) return '';
  let h = '<div class="depthwrap" style="margin-top:8px"><table><thead><tr>'+
    '<th>Sample</th><th>Read depth at variant base</th></tr></thead><tbody>';
  for(const c of cov){
    if(c.error){ h += '<tr><td>'+esc(c.label)+'</td><td class="d-lo">'+esc(c.error)+'</td></tr>'; continue; }
    const d = c.mean; const cls = d>=20?'d-hi':(d>=10?'d-mid':'d-lo');
    h += '<tr><td class="gene">'+esc(c.label)+'</td>'+
         '<td class="'+cls+'">'+Math.round(d)+'x</td></tr>';
  }
  return h+'</tbody></table></div>';
}

function sumBar(cells){
  let h = '<div class="summary-bar">';
  for(const c of cells)
    h += '<div class="s-cell"><div class="val '+(c.cls||'val-total')+'">'+c.val+
         '</div><div class="lbl">'+c.lbl+'</div></div>';
  return h+'</div>';
}

let REPORT_NAME = "coverage";   // base filename for the PDF (no extension)
function safeName(s){ return (s||'coverage').replace(/[^A-Za-z0-9 _-]+/g,'_').trim() || 'coverage'; }
function printReport(){
  const prev = document.title;
  document.title = REPORT_NAME;          // browsers use this as the PDF filename
  window.print();
  setTimeout(()=>{ document.title = prev; }, 500);
}
function dlButtons(){
  const q = encodeURIComponent($('#q').value.trim());
  const t = encodeURIComponent($('#qtype').value);
  const p = encodeURIComponent($('#panel').value);
  return '<div class="toolbar">'+
    '<a href="/api/bed?q='+q+'&type='+t+'&panel='+p+'"><button class="ghost" type="button">Download BED</button></a>'+
    '<button class="ghost" type="button" onclick="printReport()">Print / Save PDF</button>'+
    '</div>';
}

function render(d){
  if(!d.ok){ out.innerHTML = '<div class="err">'+esc(d.error||'Query failed.')+'</div>'; return; }
  let h = '';
  for(const m of (d.messages||[])) h += '<div class="msg">'+esc(m)+'</div>';

  if(d.type === 'genes'){
    REPORT_NAME = 'listed gene coverage';
    h += '<h2>Multiple-gene coverage — '+d.n_found+' / '+d.n_genes+' found</h2>';
    h += sumBar([
      {val: d.n_genes, lbl:'Genes queried'},
      {val: d.n_found, lbl:'Found in panel', cls: d.n_found===d.n_genes?'val-green':'val-red'},
      {val: fmt(d.combined_summary.targeted_bp), lbl:'Combined targeted bp', cls:'val-green'},
      {val: d.combined_summary.intervals, lbl:'Target intervals'},
    ]);
    h += dlButtons();
    h += depthTable(d.combined_sample_coverage, {title:'Combined sample read-depth coverage (all listed genes)'});
    h += referencePanel(d.combined_sample_coverage);
    for(const g of d.genes){
      h += '<h2>'+esc(g.name)+(g.found?'':' — <span class="cov-no">not found</span>')+'</h2>';
      if(!g.found){ h += '<div class="empty">No target intervals for '+esc(g.name)+' in this panel.</div>'; continue; }
      const s = g.summary;
      h += '<div class="empty" style="padding:4px 0">'+
           '<b>'+s.intervals+'</b> intervals · <b>'+fmt(s.targeted_bp)+'</b> targeted bp · '+
           esc(s.chroms.join(', '))+'</div>';
      h += depthTable(g.sample_coverage, {title:'Read-depth coverage — '+g.name});
    }
    out.innerHTML = h; return;
  }
  else if(d.type === 'region'){
    const cov = d.coverage, rg = d.region;
    REPORT_NAME = safeName(rg.chrom+'_'+rg.start+'-'+rg.end)+' coverage';
    h += '<h2>Region coverage — '+esc(rg.chrom)+':'+fmt(rg.start)+'-'+fmt(rg.end)+'</h2>';
    h += sumBar([
      {val: fmt(rg.length), lbl:'Region size (bp)'},
      {val: fmt(cov.covered_bp), lbl:'Covered by panel (bp)', cls:'val-green'},
      {val: cov.pct+'%', lbl:'Percent covered', cls: cov.pct>=99.9?'val-green':(cov.pct<=0.01?'val-red':'val-total')},
      {val: d.summary.intervals, lbl:'Target intervals'},
      {val: cov.gaps.length, lbl:'Uncovered gaps', cls: cov.gaps.length?'val-red':'val-green'},
    ]);
    h += '<div class="prog-wrap"><div class="prog-outer"><div class="prog-fill" style="width:'+
         Math.min(100,cov.pct)+'%"></div></div></div>';
    h += dlButtons();
    h += depthTable(d.sample_coverage, {title:'Sample read-depth coverage over captured region'});
    h += referencePanel(d.sample_coverage);
    h += '<h2>Target intervals overlapping region</h2>' + rowsTable(d.rows);
    if(cov.gaps && cov.gaps.length){
      h += '<h2>Uncovered gaps</h2><div class="table-wrap"><table><thead><tr>'+
           '<th>Chrom</th><th>Gap start</th><th>Gap end</th><th>Length (bp)</th></tr></thead><tbody>';
      for(const g of cov.gaps)
        h += '<tr><td class="mono">'+esc(rg.chrom)+'</td><td class="mono">'+fmt(g.start)+
             '</td><td class="mono">'+fmt(g.end)+'</td><td class="mono cov-no">'+fmt(g.length)+'</td></tr>';
      h += '</tbody></table></div>';
    }
  }
  else if(d.type === 'rsid'){
    REPORT_NAME = safeName(d.rsid.id)+' coverage';
    h += '<h2>Variant coverage — '+esc(d.rsid.id)+'</h2>';
    h += sumBar([
      {val: d.covered?'YES':'NO', lbl:'In panel target', cls: d.covered?'val-green':'val-red'},
      {val: esc(d.rsid.consequence||'—'), lbl:'Most severe consequence'},
      {val: d.hits.length, lbl:'Genomic mapping(s)'},
    ]);
    h += dlButtons();
    for(const hit of d.hits){
      const p = hit.position;
      h += '<h2>Position '+esc(p.chrom)+':'+fmt(p.start)+'  ('+esc(p.allele||'')+')  '+
           '<span class="pill '+(hit.covered?'yes':'no')+'">'+(hit.covered?'COVERED':'NOT COVERED')+'</span></h2>';
      h += depthAtVariant(hit.sample_coverage);
      if(hit.covered){ h += rowsTable(hit.rows); }
      else if(hit.nearest){
        const n = hit.nearest;
        h += '<div class="empty">Nearest target interval: <span class="mono">'+esc(n.chrom)+':'+
             fmt(n.start)+'-'+fmt(n.end)+'</span> ('+fmt(n.distance)+' bp away) — '+
             '<span class="gene">'+esc((n.genes||[]).join(', '))+'</span></div>';
      } else { h += '<div class="empty">No nearby target interval.</div>'; }
    }
  }
  else { // gene / transcript
    const s = d.summary;
    REPORT_NAME = safeName(d.query)+' coverage';
    h += '<h2>'+(d.type==='transcript'?'Transcript':'Gene')+' coverage — '+esc(d.query)+'</h2>';
    if(!d.found){
      h += '<div class="err">No target intervals found for <b>'+esc(d.query)+
           '</b> in this panel.</div>';
      out.innerHTML = h; return;
    }
    h += sumBar([
      {val: s.intervals, lbl:'Target intervals'},
      {val: fmt(s.targeted_bp), lbl:'Targeted bases (bp)', cls:'val-green'},
      {val: esc(s.chroms.join(', ')), lbl:'Chromosome'},
      {val: s.span_start!=null? fmt(s.span_end - s.span_start):'—', lbl:'Genomic span (bp)'},
    ]);
    h += dlButtons();
    h += depthTable(d.sample_coverage, {title:'Sample read-depth coverage over targets'});
    h += referencePanel(d.sample_coverage);
    h += '<h2>Target intervals</h2>' + rowsTable(d.rows);
  }
  out.innerHTML = h;
}

async function run(){
  const q = $('#q').value.trim();
  if(!q) return;
  out.innerHTML = '<div style="padding:20px"><span class="spinner"></span> '+
    '<span style="margin-left:8px;color:#005FA0">Checking BED + sample read-depth coverage…</span></div>';
  const params = new URLSearchParams({q, type:$('#qtype').value, panel:$('#panel').value});
  try {
    const r = await fetch('/api/query?'+params.toString());
    const d = await r.json();
    render(d);
  } catch(e){
    out.innerHTML = '<div class="err">Request failed: '+esc(e.message)+'</div>';
  }
}
</script>
</body>
</html>"""


if __name__ == "__main__":
    # warm the default panel so the first query is fast
    try:
        t = time.time()
        get_panel()
        print(f"Default panel loaded in {time.time()-t:.1f}s "
              f"({PANELS[PANEL_ORDER[0]].n_rows:,} intervals)")
    except Exception as e:
        print("WARN:", e)
    # HOST defaults to 0.0.0.0 (LAN access, by design). Set HOST=127.0.0.1 to
    # restrict to localhost when running behind a reverse proxy / on shared hosts.
    app.run(host=os.environ.get("HOST", "0.0.0.0"),
            port=int(os.environ.get("PORT", 8080)), debug=False)
