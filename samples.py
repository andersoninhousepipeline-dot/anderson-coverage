"""
samples.py — per-sample sequencing-depth coverage from BAM files.

Computes real coverage (read depth) over a set of target intervals for one or
more samples, and summarises it as mean depth + % of target bases at or above
depth thresholds (1x / 10x / 20x / 30x / 50x / 100x) — the standard panel QC view.

Depth is measured with pysam.count_coverage (sums A/C/G/T per base). Reads flagged
unmapped / secondary / QC-fail / duplicate are excluded (read_callback='all').
BAMs here are already deduplicated.
"""

import os
import numpy as np
import pysam

BAM_DIR = os.environ.get("BAM_DIR", "/data/bed/reference")

# Public reference samples (safe to commit). Order = display order.
#   (label, public slug)   — NO raw sample identifiers here.
SAMPLE_DEFS = [
    ("Normal Male",        "normal-male"),
    ("Normal Female",      "normal-female"),
    ("Male Infertility",   "male-infertility"),
    ("Female Infertility", "female-infertility"),
    ("AF",                 "af"),
    ("POC",                "poc"),
]

# Map public slug -> real BAM filename. Kept in a git-ignored local file so the
# underlying sample identifiers are never committed to the repository.
try:
    from samples_local import SAMPLE_FILES
except Exception:
    SAMPLE_FILES = {}

# Depth thresholds reported as "% of target bases >= X".
THRESHOLDS = [1, 10, 20, 30, 50, 100]

# Safety cap: don't attempt depth over more than this many target bases in one query.
MAX_DEPTH_BASES = 8_000_000


class Sample:
    def __init__(self, label, sid, filename):
        self.label = label
        self.id = sid                 # public slug (never a raw sample number)
        self.path = os.path.join(BAM_DIR, filename) if filename else ""
        self.filename = filename
        self.available = bool(filename) and os.path.exists(self.path)
        self._bam = None

    @property
    def name(self):
        return self.label             # label only — no identifiers exposed

    def bam(self):
        if self._bam is None:
            self._bam = pysam.AlignmentFile(self.path, "rb")
        return self._bam


SAMPLES = {}        # slug -> Sample
SAMPLE_ORDER = []   # display order
for _label, _sid in SAMPLE_DEFS:
    s = Sample(_label, _sid, SAMPLE_FILES.get(_sid))
    if s.available:
        SAMPLES[_sid] = s
        SAMPLE_ORDER.append(_sid)


def list_samples():
    return [{"id": a, "label": SAMPLES[a].label} for a in SAMPLE_ORDER]


def _merge(intervals):
    """Merge (chrom,start,end) list -> dict chrom -> merged [(s,e)]."""
    by_chrom = {}
    for c, s, e in intervals:
        if e <= s:
            continue
        by_chrom.setdefault(c, []).append((s, e))
    out = {}
    for c, ivs in by_chrom.items():
        ivs.sort()
        merged = [list(ivs[0])]
        for s, e in ivs[1:]:
            if s <= merged[-1][1]:
                merged[-1][1] = max(merged[-1][1], e)
            else:
                merged.append([s, e])
        out[c] = [(a, b) for a, b in merged]
    return out


def coverage_for_intervals(sid, intervals, thresholds=THRESHOLDS):
    """
    intervals: list of (chrom, start, end), BED 0-based half-open.
    Returns dict with mean/min/max depth, total bases, and pct>=t for each t.
    """
    sample = SAMPLES.get(sid)
    if sample is None:
        return {"error": f"Unknown sample {sid}"}

    merged = _merge(intervals)
    total_bases = sum(e - s for ivs in merged.values() for s, e in ivs)
    if total_bases == 0:
        return {"id": sid, "label": sample.label,
                "total_bases": 0, "mean": 0, "min": 0, "max": 0,
                "pct": {str(t): None for t in thresholds}}
    if total_bases > MAX_DEPTH_BASES:
        return {"id": sid, "label": sample.label,
                "total_bases": total_bases,
                "error": f"Region too large for live depth "
                         f"({total_bases:,} bp > {MAX_DEPTH_BASES:,} cap)."}

    bam = sample.bam()
    sum_depth = 0
    bmin = None
    bmax = 0
    at = {t: 0 for t in thresholds}
    # histogram (capped) for median
    HCAP = 2000
    hist = np.zeros(HCAP + 1, dtype=np.int64)

    for chrom, ivs in merged.items():
        for s, e in ivs:
            try:
                cov = bam.count_coverage(chrom, s, e, quality_threshold=0)
            except (ValueError, KeyError):
                # contig not in this BAM
                continue
            depth = np.asarray(cov, dtype=np.int64).sum(axis=0)
            if depth.size == 0:
                continue
            sum_depth += int(depth.sum())
            dmn = int(depth.min()); dmx = int(depth.max())
            bmin = dmn if bmin is None else min(bmin, dmn)
            bmax = max(bmax, dmx)
            for t in thresholds:
                at[t] += int(np.count_nonzero(depth >= t))
            clipped = np.clip(depth, 0, HCAP)
            hist += np.bincount(clipped, minlength=HCAP + 1)

    mean = sum_depth / total_bases
    # median from histogram
    csum = np.cumsum(hist)
    half = total_bases / 2.0
    median = int(np.searchsorted(csum, half))

    return {
        "id": sid,
        "label": sample.label,
        "total_bases": total_bases,
        "mean": round(mean, 1),
        "median": median,
        "min": bmin if bmin is not None else 0,
        "max": bmax,
        "pct": {str(t): round(100.0 * at[t] / total_bases, 1) for t in thresholds},
        "bases_at": {str(t): at[t] for t in thresholds},
    }


def coverage_multi(sids, intervals, thresholds=THRESHOLDS):
    return [coverage_for_intervals(a, intervals, thresholds) for a in sids]
