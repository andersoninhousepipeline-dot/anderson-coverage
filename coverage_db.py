"""
coverage_db.py — read precomputed sample coverage from the SQLite DB.

When build/coverage.db (default ./coverage.db) exists, the app uses this instead
of the BAM files: the reference panel intervals and per-type coverage all come
from the DB, so the tool is fully self-contained (works straight from the repo,
no BAMs required).

Per query we aggregate, for a set of interval ids, each sample type's:
  mean depth  = sum_depth / (Σ bp × n_replicates)              (exact)
  % ≥ t       = 100 × bases≥t / (Σ bp × n_replicates)          (exact)
  min (tier)  = highest threshold fully covered across all bases (approx)
Median and per-replicate SD are not stored (per-base histograms would bloat the
DB); median is reported as None and SD as 0 in DB mode.
"""

import os
import json
import zlib
import struct
import sqlite3

from coverage_index import Panel

_TYPE_STRUCT = struct.Struct("<Q6I")   # per type: sum_depth(uint64) + 6 counts(uint32)
_TYPE_BYTES = _TYPE_STRUCT.size        # 32


class CoverageDB:
    def __init__(self, path):
        self.path = path
        self.con = sqlite3.connect(f"file:{path}?mode=ro", uri=True,
                                   check_same_thread=False)
        cur = self.con.cursor()
        meta = dict(cur.execute("SELECT key, value FROM meta").fetchall())
        self.thresholds = json.loads(meta.get("thresholds", "[1,10,20,30,50,100]"))
        self.panel_name = meta.get("panel_name", "Reference panel")
        self.genome_build = meta.get("genome_build", "GRCh38 / hg38")
        self.build_date = meta.get("build_date", "")
        self.types = [dict(slug=s, label=l, n=n, ord=o) for (s, l, n, o) in
                      cur.execute("SELECT slug,label,n,ord FROM sample_type ORDER BY ord")]
        self._ord = {t["slug"]: t["ord"] for t in self.types}
        self._meta_n = {t["slug"]: t["n"] for t in self.types}
        # Build the in-memory panel index from the DB intervals (id == row order).
        # Full annotation is restored losslessly from the compressed pack.
        (blob,) = cur.execute("SELECT data FROM annot_pack").fetchone()
        annots = zlib.decompress(blob).decode().split("\n")
        coords = cur.execute("SELECT chrom,start,end FROM intervals ORDER BY id").fetchall()
        rows = ((c, s, e, annots[i] if i < len(annots) else "")
                for i, (c, s, e) in enumerate(coords))
        self.panel = Panel(self.panel_name, path).load_rows(rows)

    # ── sample list (no identifiers) ─────────────────────────────────────────
    def list_samples(self):
        return [{"id": t["slug"], "label": t["label"], "n": t["n"]} for t in self.types]

    @property
    def sample_order(self):
        return [t["slug"] for t in self.types]

    # ── coverage over a set of interval ids ──────────────────────────────────
    def coverage(self, slugs, idxs):
        slugs = [s for s in slugs if s in self._ord]
        labels = {t["slug"]: t["label"] for t in self.types}
        bp_total = sum(self.panel.r_end[i] - self.panel.r_start[i] for i in idxs)

        # accumulate sum_depth and threshold counts per requested slug
        acc = {s: [0, [0] * 6] for s in slugs}     # slug -> [sum_depth, [t-counts]]
        if idxs and bp_total > 0:
            cur = self.con.cursor()
            ids = list(idxs)
            for off in range(0, len(ids), 900):
                chunk = ids[off:off + 900]
                q = "SELECT data FROM cov WHERE interval_id IN (%s)" % \
                    ",".join("?" * len(chunk))
                for (blob,) in cur.execute(q, chunk):
                    for s in slugs:
                        sd, *tc = _TYPE_STRUCT.unpack_from(blob, self._ord[s] * _TYPE_BYTES)
                        a = acc[s]
                        a[0] += sd
                        for k in range(6):
                            a[1][k] += tc[k]

        out = []
        for s in slugs:
            nrep = self._meta_n[s] or 1
            denom = bp_total * nrep
            base = {"id": s, "label": labels[s], "n": self._meta_n[s],
                    "total_bases": bp_total, "median": None, "mean_sd": 0.0}
            if denom == 0:
                base.update({"mean": 0, "min": 0,
                             "pct": {str(t): None for t in self.thresholds}})
                out.append(base); continue
            sd, tc = acc[s]
            pct = {str(self.thresholds[k]): round(100.0 * tc[k] / denom, 1)
                   for k in range(6)}
            # approx min tier: highest threshold whose bases cover ALL positions
            min_tier = 0
            for k in range(6):
                if tc[k] >= denom:
                    min_tier = self.thresholds[k]
            base.update({"mean": round(sd / denom, 1), "min": min_tier, "pct": pct})
            out.append(base)
        return out


def open_db(path):
    return CoverageDB(path) if path and os.path.exists(path) else None
