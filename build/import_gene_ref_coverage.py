#!/usr/bin/env python3
"""
import_gene_ref_coverage.py — load an external per-gene "% of coding region
covered" reference table into coverage.db, as a separate, clearly-labeled
table alongside our own mosdepth-derived panel/depth data.

This does NOT replace or recompute anything from our own BAMs/BED — it is a
cross-reference from a different lab's panel (SurfSeq), used to flag genes
where coverage is known to be incomplete (e.g. multi-copy paralogs, segmental
duplications) even when our own panel shows bait design present.

Input: TSV with `gene\tpct_coding_covered` (see build/ref/*.tsv for provenance
comments — source URL, upstream commit, fetch date).

Schema added:
  gene_ref_coverage(gene TEXT PRIMARY KEY, pct_coding_covered REAL)
  meta gains: gene_ref_source, gene_ref_source_date, gene_ref_n_genes
"""

import os
import sys
import sqlite3

HERE = os.path.dirname(os.path.abspath(__file__))
DEFAULT_DB = os.path.join(os.path.dirname(HERE), "coverage.db")
DEFAULT_TSV = os.path.join(HERE, "ref", "surfseq_gene_coverage_2025-07.tsv")

SOURCE = "SurfSeq Gene Coverage (andersonbioinfo.github.io/SurfSeq_Gene_coverage)"
SOURCE_DATE = "2025-07-26"


def load_tsv(path):
    rows = []
    with open(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            gene, pct = line.rstrip("\n").split("\t")
            if gene == "gene":   # header
                continue
            rows.append((gene.strip().upper(), float(pct)))
    return rows


def main():
    db_path = sys.argv[1] if len(sys.argv) > 1 else DEFAULT_DB
    tsv_path = sys.argv[2] if len(sys.argv) > 2 else DEFAULT_TSV

    print("Loading reference TSV:", tsv_path)
    rows = load_tsv(tsv_path)
    print(f"  {len(rows):,} genes")

    print("Writing into:", db_path)
    con = sqlite3.connect(db_path)
    cur = con.cursor()
    cur.execute("DROP TABLE IF EXISTS gene_ref_coverage")
    cur.execute("CREATE TABLE gene_ref_coverage(gene TEXT PRIMARY KEY, pct_coding_covered REAL)")
    cur.executemany("INSERT INTO gene_ref_coverage VALUES(?,?)", rows)
    cur.executemany("INSERT OR REPLACE INTO meta VALUES(?,?)", [
        ("gene_ref_source", SOURCE),
        ("gene_ref_source_date", SOURCE_DATE),
        ("gene_ref_n_genes", str(len(rows))),
    ])
    con.commit()
    con.close()
    print("Done.")


if __name__ == "__main__":
    main()
