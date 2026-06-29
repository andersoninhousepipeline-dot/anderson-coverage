#!/usr/bin/env python3
"""
build_db.py — assemble the precomputed coverage SQLite DB from mosdepth outputs.

Lossless by design:
  • Coverage numbers are exact — sum_depth (uint64) and per-threshold base counts
    (uint32) are summed across replicates with no rounding; mosdepth's threshold
    counts are exact integers, so every "% ≥ Nx" is exact.
  • The FULL BED annotation (nothing dropped, ClinID included) is preserved,
    stored zlib-compressed (~9 MB) so gene / transcript / ClinID lookups and BED
    export are identical to the BAM-mode tool.

Not stored (need per-base histograms; available only in BAM/server mode):
  • per-region median, per-replicate SD.

Schema
  meta(key, value)
  sample_type(slug, label, n, ord)
  intervals(id INTEGER PRIMARY KEY, chrom, start, end)
  annot_pack(data BLOB)                       # single row: zlib("\n".join(annots))
  cov(interval_id INTEGER PRIMARY KEY, data BLOB)
      blob = per sample type in `ord`: struct '<Q6I'
             = sum_depth(uint64) + base-counts ≥1,10,20,30,50,100x (uint32)
"""

import os
import gzip
import json
import zlib
import struct
import sqlite3
import datetime
import glob

HERE = os.path.dirname(os.path.abspath(__file__))
MOSDIR = os.path.join(HERE, "mosdepth")
PANEL_BED = "/data/bed/hg38_exome_comp_spikein_v2.0.2_targets_sorted.re_annotated.bed"
PANEL_NAME = "Twist Spikein v2.0.2 (re-annotated)"
GENOME_BUILD = "GRCh38 / hg38"
THRESHOLDS = [1, 10, 20, 30, 50, 100]
OUT_DB = os.path.join(os.path.dirname(HERE), "coverage.db")

TYPES = [
    ("normal-male", "Normal Male"),
    ("normal-female", "Normal Female"),
    ("male-infertility", "Male Infertility"),
    ("female-infertility", "Female Infertility"),
    ("af", "AF"),
    ("poc", "POC"),
]

_PACK = struct.Struct("<Q6I")   # per type: sum_depth(uint64) + 6 counts(uint32)


def load_panel():
    intervals = []
    with open(PANEL_BED) as fh:
        for line in fh:
            if not line or line[0] == "#":
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 3:
                continue
            intervals.append((p[0], int(p[1]), int(p[2]),
                              p[3] if len(p) > 3 else ""))   # FULL annot, untouched
    return intervals


def replicate_files(slug):
    out = []
    for r in sorted(glob.glob(os.path.join(MOSDIR, f"{slug}__*.regions.bed.gz"))):
        thr = r.replace(".regions.bed.gz", ".thresholds.bed.gz")
        if os.path.exists(thr):
            out.append((r, thr))
    return out


# NOTE: mosdepth emits regions in BAM/karyotypic contig order, which differs from
# the BED's (lexical) order — so we align every row to the BED interval id by its
# (chrom,start,end) key rather than by line number. Coordinates are unique.

def read_means(path, key2id, n):
    means = [0.0] * n
    seen = 0
    with gzip.open(path, "rt") as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            iid = key2id.get((f[0], int(f[1]), int(f[2])))
            if iid is not None:
                means[iid] = float(f[4]); seen += 1
    assert seen == n, f"{path}: matched {seen} != {n}"
    return means


def read_thresholds(path, key2id, n):
    counts = [None] * n
    seen = 0
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            iid = key2id.get((f[0], int(f[1]), int(f[2])))
            if iid is not None:
                counts[iid] = tuple(int(x) for x in f[4:10]); seen += 1
    assert seen == n, f"{path}: matched {seen} != {n}"
    return counts


def main():
    print("Loading panel:", PANEL_BED)
    intervals = load_panel()
    n = len(intervals)
    bp = [e - s for (_c, s, e, _a) in intervals]
    key2id = {(c, s, e): i for i, (c, s, e, _a) in enumerate(intervals)}
    print(f"  {n:,} intervals")

    type_sd = {slug: [0.0] * n for slug, _ in TYPES}
    type_tc = {slug: [[0] * 6 for _ in range(n)] for slug, _ in TYPES}
    type_n = {}

    for slug, label in TYPES:
        reps = replicate_files(slug)
        type_n[slug] = len(reps)
        print(f"{label}: {len(reps)} replicate(s)")
        for regf, thrf in reps:
            means = read_means(regf, key2id, n)
            counts = read_thresholds(thrf, key2id, n)
            sd = type_sd[slug]; tc = type_tc[slug]
            for i in range(n):
                sd[i] += means[i] * bp[i]
                ci = counts[i]; ti = tc[i]
                for k in range(6):
                    ti[k] += ci[k]

    print("Writing DB:", OUT_DB)
    if os.path.exists(OUT_DB):
        os.remove(OUT_DB)
    con = sqlite3.connect(OUT_DB); cur = con.cursor()
    cur.execute("CREATE TABLE meta(key TEXT PRIMARY KEY, value TEXT)")
    cur.execute("CREATE TABLE sample_type(slug TEXT PRIMARY KEY, label TEXT, n INTEGER, ord INTEGER)")
    cur.execute("CREATE TABLE intervals(id INTEGER PRIMARY KEY, chrom TEXT, start INTEGER, end INTEGER)")
    cur.execute("CREATE TABLE annot_pack(data BLOB)")
    cur.execute("CREATE TABLE cov(interval_id INTEGER PRIMARY KEY, data BLOB)")

    cur.executemany("INSERT INTO meta VALUES(?,?)", [
        ("thresholds", json.dumps(THRESHOLDS)),
        ("panel_name", PANEL_NAME),
        ("genome_build", GENOME_BUILD),
        ("build_date", datetime.date.today().isoformat()),
        ("n_intervals", str(n)),
    ])
    cur.executemany("INSERT INTO sample_type VALUES(?,?,?,?)",
                    [(slug, label, type_n[slug], i) for i, (slug, label) in enumerate(TYPES)])
    cur.executemany("INSERT INTO intervals VALUES(?,?,?,?)",
                    [(i, intervals[i][0], intervals[i][1], intervals[i][2]) for i in range(n)])

    annot_blob = zlib.compress("\n".join(a for (_c, _s, _e, a) in intervals).encode(), 9)
    cur.execute("INSERT INTO annot_pack VALUES(?)", (annot_blob,))
    print(f"  annot pack: {len(annot_blob)/1e6:.1f} MB")

    def blobs():
        for i in range(n):
            yield (i, b"".join(_PACK.pack(int(round(type_sd[s][i])), *type_tc[s][i])
                               for s, _ in TYPES))
    cur.executemany("INSERT INTO cov VALUES(?,?)", blobs())

    con.commit(); con.execute("VACUUM"); con.close()
    print(f"Done. {OUT_DB}  ({os.path.getsize(OUT_DB)/1e6:.1f} MB)")


if __name__ == "__main__":
    main()
