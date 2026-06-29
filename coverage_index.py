"""
coverage_index.py — Panel-coverage index for annotated BED target files.

A "panel" is a capture-target BED of the form:
    chrom \t start \t end \t Gene;ENST..;NM_..;...   (0-based half-open, BED standard)

The index supports:
  • interval overlap queries (chr region, rs-id position)   -> bisect on sorted starts
  • token lookups (gene symbol, NM/NR/XM/XR id, ENST/ENSG)  -> token -> row map
  • free-text substring fallback                            -> linear scan

Coverage here = panel/target coverage: which portion of a gene / region / variant
is included in the capture design. These BEDs carry no per-base sequencing depth.
"""

import os
import re
import bisect
from collections import defaultdict

# Tokens we deliberately keep OUT of the token map to avoid memory blow-up.
# (ClinID-* tokens are near-unique per row; CCDS adds little lookup value.)
_SKIP_TOKEN_PREFIXES = ("ClinID",)

# Token prefixes that identify a transcript/ID (everything else is treated as a gene symbol).
_ID_PREFIXES = ("NM_", "NR_", "XM_", "XR_", "ENST", "ENSG", "CCDS")

# BED 4th-column annotations come in many vendor styles, e.g.
#   Twist:     OR4F5;ENST00000641515.2;NM_001005484.2
#   Roche:     gene_symbol=OR4F5;hgnc_id=14825;ensembl_gene_id=ENSG00000186092
#   Sophia:    AGRN:NM_001305275.2:ex1;AGRN:NM_198576.4:ex1
#   Nonacus:   Target_1_NPHP4_1,Target_1_NPHP4_2
#   aviti:     OR4F5    1   (tab-separated extra columns)
# Splitting on all of these delimiters yields searchable gene/transcript tokens.
_TOKEN_SPLIT = re.compile(r"[;,:=|\t ]+")

# Keys / non-informative tokens we don't index as searchable terms.
_STOP_TOKENS = {
    "gene_symbol", "hgnc_id", "ensembl_gene_id", "ncbi_gene_id",
    "regionidxkey", "no_annotation", "target", "covered", "intron",
}


def tokenize(annot):
    """Split a BED annotation into candidate gene / transcript tokens."""
    return [t for t in _TOKEN_SPLIT.split(annot) if t]


class Panel:
    def __init__(self, name, path):
        self.name = name
        self.path = path
        self.loaded = False
        # parallel row arrays
        self.r_chrom = []
        self.r_start = []
        self.r_end = []
        self.r_annot = []
        # per-chromosome sorted index: chrom -> dict(starts, ends, idx, max_len)
        self.by_chrom = {}
        # token (upper) -> list of row indices
        self.token_map = defaultdict(list)
        self.n_rows = 0
        self.total_bp = 0

    # ── loading ────────────────────────────────────────────────────────────
    def load(self):
        if self.loaded:
            return self
        chrom_rows = defaultdict(list)  # chrom -> [(start, end, idx)]
        with open(self.path) as fh:
            i = 0
            for line in fh:
                if not line or line[0] == "#":
                    continue
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 3:
                    continue
                chrom = parts[0]
                try:
                    start = int(parts[1])
                    end = int(parts[2])
                except ValueError:
                    continue
                annot = parts[3] if len(parts) > 3 else ""
                self.r_chrom.append(chrom)
                self.r_start.append(start)
                self.r_end.append(end)
                self.r_annot.append(annot)
                self.total_bp += max(0, end - start)
                chrom_rows[chrom].append((start, end, i))
                # token map
                if annot:
                    for tok in tokenize(annot):
                        if (tok.startswith(_SKIP_TOKEN_PREFIXES)
                                or tok.isdigit()
                                or tok.lower() in _STOP_TOKENS):
                            continue
                        self.token_map[tok.upper()].append(i)
                i += 1
            self.n_rows = i

        # build per-chrom sorted arrays
        for chrom, rows in chrom_rows.items():
            rows.sort(key=lambda r: r[0])
            starts = [r[0] for r in rows]
            ends = [r[1] for r in rows]
            idx = [r[2] for r in rows]
            max_len = max((e - s for s, e in zip(starts, ends)), default=0)
            self.by_chrom[chrom] = {
                "starts": starts, "ends": ends, "idx": idx, "max_len": max_len,
            }
        self.loaded = True
        return self

    # ── core queries ─────────────────────────────────────────────────────────
    def overlap(self, chrom, qstart, qend):
        """Return list of row indices whose interval overlaps [qstart, qend)."""
        c = self.by_chrom.get(chrom)
        if not c:
            return []
        starts, ends, idx, max_len = c["starts"], c["ends"], c["idx"], c["max_len"]
        lo = bisect.bisect_left(starts, qstart - max_len)
        if lo < 0:
            lo = 0
        out = []
        n = len(starts)
        i = lo
        while i < n and starts[i] < qend:
            if ends[i] > qstart:
                out.append(idx[i])
            i += 1
        return out

    def rows_for_token(self, token):
        return self.token_map.get(token.upper(), [])

    def free_text(self, needle, limit=2000):
        """Substring scan over annotations (case-insensitive). For ClinID / partial."""
        needle = needle.upper()
        out = []
        for i, annot in enumerate(self.r_annot):
            if needle in annot.upper():
                out.append(i)
                if len(out) >= limit:
                    break
        return out

    def row(self, i):
        return (self.r_chrom[i], self.r_start[i], self.r_end[i], self.r_annot[i])


# ── helpers shared by the web layer ──────────────────────────────────────────
def merge_intervals(intervals):
    """intervals: list of (start, end) -> merged, sorted, non-overlapping."""
    if not intervals:
        return []
    s = sorted(intervals)
    merged = [list(s[0])]
    for st, en in s[1:]:
        if st <= merged[-1][1]:
            merged[-1][1] = max(merged[-1][1], en)
        else:
            merged.append([st, en])
    return [(a, b) for a, b in merged]


def classify_token(tok):
    if tok.startswith(_ID_PREFIXES):
        return "id"
    return "gene"
