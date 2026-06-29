"""
samples.py — per-sample-type sequencing-depth coverage from BAM files.

Reference samples are organised as one sub-directory per type under BAM_DIR,
each holding one or more replicate BAMs:

    BAM_DIR/
      Normal/      *.bam   (+ .bai)
      Male-Inf/    *.bam
      Female-Inf/  *.bam
      AF/          *.bam
      POC/         *.bam

For a query, each replicate's depth is measured with pysam.count_coverage over
the target intervals (duplicate / secondary / QC-fail reads excluded), then the
replicates of a type are aggregated to a single row: mean of replicate means,
mean of replicate medians, worst-case min, and mean of replicate %≥threshold.

No raw sample identifiers are exposed — only the type label and replicate count.
Sample → file details stay in the (git-ignored) directory layout, not in code.
"""

import os
import numpy as np
import pysam
from concurrent.futures import ProcessPoolExecutor

BAM_DIR = os.environ.get("BAM_DIR", "/data/bed/reference")

# Directory name -> (display label, public slug). Order = display order.
DIR_MAP = [
    ("Normal-Male",   "Normal Male",        "normal-male"),
    ("Normal-Female", "Normal Female",      "normal-female"),
    ("Male-Inf",      "Male Infertility",   "male-infertility"),
    ("Female-Inf",    "Female Infertility", "female-infertility"),
    ("AF",            "AF",                 "af"),
    ("POC",           "POC",                "poc"),
]

# Depth thresholds reported as "% of target bases >= X".
THRESHOLDS = [1, 10, 20, 30, 50, 100]

# Safety cap: don't attempt depth over more than this many target bases in one query.
MAX_DEPTH_BASES = 8_000_000

# Parallel workers for replicate depth (pysam releases the GIL inside htslib).
_MAX_WORKERS = int(os.environ.get("DEPTH_WORKERS", "8"))


class SampleGroup:
    def __init__(self, label, sid, dirname):
        self.label = label
        self.id = sid
        self.dir = os.path.join(BAM_DIR, dirname)
        self.paths = []
        if os.path.isdir(self.dir):
            for fn in sorted(os.listdir(self.dir)):
                if fn.endswith(".bam"):
                    self.paths.append(os.path.join(self.dir, fn))
        self.available = len(self.paths) > 0
        self.n = len(self.paths)


SAMPLES = {}        # slug -> SampleGroup
SAMPLE_ORDER = []   # display order
for _dirname, _label, _sid in DIR_MAP:
    g = SampleGroup(_label, _sid, _dirname)
    if g.available:
        SAMPLES[_sid] = g
        SAMPLE_ORDER.append(_sid)


def list_samples():
    return [{"id": a, "label": SAMPLES[a].label, "n": SAMPLES[a].n}
            for a in SAMPLE_ORDER]


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


# Process-local BAM handle cache (each worker process opens its own handles).
_bam_handles = {}

def _bam(path):
    h = _bam_handles.get(path)
    if h is None:
        h = pysam.AlignmentFile(path, "rb")
        _bam_handles[path] = h
    return h


def _depth_task(args):
    """Top-level (picklable) worker: per-base depth stats for one BAM."""
    path, merged, total_bases, thresholds = args
    bam = _bam(path)
    sum_depth = 0
    bmin = None
    bmax = 0
    at = {t: 0 for t in thresholds}
    HCAP = 2000
    hist = np.zeros(HCAP + 1, dtype=np.int64)
    for chrom, ivs in merged.items():
        for s, e in ivs:
            try:
                cov = bam.count_coverage(chrom, s, e, quality_threshold=0)
            except (ValueError, KeyError):
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
            hist += np.bincount(np.clip(depth, 0, HCAP), minlength=HCAP + 1)
    csum = np.cumsum(hist)
    median = int(np.searchsorted(csum, total_bases / 2.0))
    return {
        "mean": sum_depth / total_bases,
        "median": median,
        "min": bmin if bmin is not None else 0,
        "max": bmax,
        "pct": {t: 100.0 * at[t] / total_bases for t in thresholds},
    }


# Persistent worker pool (fork-based on Linux; reused across requests).
_POOL = None

def _get_pool():
    global _POOL
    if _POOL is None:
        _POOL = ProcessPoolExecutor(max_workers=_MAX_WORKERS)
    return _POOL


def _aggregate(sid, group, reps, total_bases, thresholds):
    means = [r["mean"] for r in reps]
    return {
        "id": sid, "label": group.label, "n": group.n,
        "total_bases": total_bases,
        "mean": round(float(np.mean(means)), 1),
        "mean_sd": round(float(np.std(means)), 1) if len(means) > 1 else 0.0,
        "median": int(round(float(np.mean([r["median"] for r in reps])))),
        "min": min(r["min"] for r in reps),          # worst-case across replicates
        "max": max(r["max"] for r in reps),
        "pct": {str(t): round(float(np.mean([r["pct"][t] for r in reps])), 1)
                for t in thresholds},
    }


def coverage_multi(sids, intervals, thresholds=THRESHOLDS):
    """Aggregate depth per sample type. All replicate BAMs run concurrently
    in a process pool (pysam.count_coverage holds the GIL, so threads won't do)."""
    sids = [s for s in sids if s in SAMPLES]
    merged = _merge(intervals)
    total_bases = sum(e - s for ivs in merged.values() for s, e in ivs)

    # Degenerate cases shared by all groups (intervals identical per query).
    if total_bases == 0 or total_bases > MAX_DEPTH_BASES:
        out = []
        for s in sids:
            g = SAMPLES[s]
            base = {"id": s, "label": g.label, "n": g.n, "total_bases": total_bases}
            if total_bases == 0:
                base.update({"mean": 0, "median": 0, "min": 0, "max": 0,
                             "pct": {str(t): None for t in thresholds}})
            else:
                base["error"] = (f"Region too large for live depth "
                                 f"({total_bases:,} bp > {MAX_DEPTH_BASES:,} cap).")
            out.append(base)
        return out

    # Dispatch every replicate BAM across all requested groups at once.
    tasks, owner = [], []
    for s in sids:
        for path in SAMPLES[s].paths:
            tasks.append((path, merged, total_bases, thresholds))
            owner.append(s)
    try:
        results = list(_get_pool().map(_depth_task, tasks))
    except Exception:
        results = [_depth_task(t) for t in tasks]   # fallback: serial in-process

    by_group = {s: [] for s in sids}
    for s, r in zip(owner, results):
        by_group[s].append(r)
    return [_aggregate(s, SAMPLES[s], by_group[s], total_bases, thresholds)
            for s in sids]
