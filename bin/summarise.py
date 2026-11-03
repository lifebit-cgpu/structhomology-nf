#!/usr/bin/env python3
"""Combine per-protein outputs into run-level tables.

Writes, in --outdir:
  combined_hits.tsv        top-N structure and sequence hits of every protein
  per_protein_summary.tsv  one row per protein in the samplesheet, including
                           proteins stopped by the pLDDT gate
  plddt_qc_summary.tsv     the pLDDT QC rows of every protein
  run_summary.json         database details and run-level counts

Standard library only.
"""
import argparse
import csv
import glob
import json
import os

BEST_KEYS = [
    ("top_structure_hit", "top_struct"),
    ("top_sequence_hit", "top_seq"),
    ("best_structure_informative_hit", "best_struct_inf"),
    ("best_structure_only_informative_hit", "best_struct_only_inf"),
    ("best_sequence_informative_hit", "best_seq_inf"),
]
BRIEF_FIELDS = ["target", "description", "organism", "pfam", "ec", "category",
                "fs_evalue", "fs_prob", "fs_tm_max", "fs_fident", "seq_evalue", "seq_fident"]
COUNT_KEYS = [
    "n_structure_hits", "n_sequence_hits", "n_structure_significant", "n_sequence_significant",
    "n_shared_significant", "n_structure_only_significant", "n_sequence_only_significant",
    "n_structure_informative_pdb", "n_sequence_informative_pdb",
    "n_structure_only_informative_pdb", "n_sequence_only_informative_pdb",
    "n_structure_significant_poorly_annotated",
]


def read_tsv(path):
    with open(path, newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def write_tsv(path, cols, rows):
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", extrasaction="ignore", lineterminator="\n")
        w.writeheader()
        for r in rows:
            w.writerow({c: ("" if r.get(c) is None else r.get(c)) for c in cols})


def num(x, default=0.0):
    try:
        return float(x)
    except (TypeError, ValueError):
        return default


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--top-hits-dir", required=True)
    ap.add_argument("--summaries-dir", required=True)
    ap.add_argument("--qc-dir", required=True)
    ap.add_argument("--db-info", required=True)
    ap.add_argument("--outdir", default=".")
    a = ap.parse_args()

    try:
        with open(a.db_info) as fh:
            db = json.load(fh)
    except (OSError, ValueError):
        db = {}
    db_label = db.get("name", "") + (f" ({db['version']})" if db.get("version") else "")

    qc_rows = []
    for f in sorted(glob.glob(os.path.join(a.qc_dir, "*.tsv"))):
        qc_rows += read_tsv(f)
    qc_by = {r["accession"]: r for r in qc_rows}

    summaries = {}
    for f in sorted(glob.glob(os.path.join(a.summaries_dir, "*.json"))):
        with open(f) as fh:
            s = json.load(fh)
        if s.get("accession"):
            summaries[s["accession"]] = s

    # Combined hit table
    hits, hit_cols = [], None
    for f in sorted(glob.glob(os.path.join(a.top_hits_dir, "*.tsv"))):
        rows = read_tsv(f)
        if rows and hit_cols is None:
            hit_cols = list(rows[0].keys())
        hits += rows
    hit_cols = hit_cols or ["query", "target"]
    for h in hits:
        q = qc_by.get(h.get("query"), {})
        h["query_gene"] = q.get("gene", "")
        h["query_description"] = q.get("description", "")
        h["query_mean_plddt"] = q.get("mean_plddt", "")
        h["database"] = db_label
    lead = ["query", "query_group", "query_gene", "query_description", "query_mean_plddt"]
    combined_cols = lead + [c for c in hit_cols if c not in lead] + ["database"]

    def hit_order(h):
        fr = num(h.get("fs_rank"), 1e9)
        sr = num(h.get("seq_rank"), 1e9)
        return (h.get("query_group") == "control", h.get("query", ""), min(fr, sr), h.get("target", ""))

    hits.sort(key=hit_order)
    write_tsv(os.path.join(a.outdir, "combined_hits.tsv"), combined_cols, hits)

    # Per-protein summary (every protein in the samplesheet)
    per = []
    for q in qc_rows:
        acc = q["accession"]
        s = summaries.get(acc, {})
        row = {
            "accession": acc, "group": q.get("group", ""), "gene": q.get("gene", ""),
            "description": q.get("description", ""), "afdb_entry": q.get("afdb_entry", ""),
            "model_version": q.get("model_version", ""), "n_residues": q.get("n_residues", ""),
            "mean_plddt": q.get("mean_plddt", ""), "frac_ge70": q.get("frac_ge70", ""),
            "frac_very_low_le50": q.get("frac_very_low_le50", ""), "qc_flag": q.get("qc_flag", ""),
            "searched": "true" if s else "false",
            "structure_added_class": s.get("structure_added_class", "not_searched" if not s else ""),
            "structure_gain": s.get("structure_gain", ""),
            "same_protein_in_db": s.get("same_protein_in_db", ""),
        }
        for k in COUNT_KEYS:
            row[k] = (s.get("counts") or {}).get(k, "")
        for key, prefix in BEST_KEYS:
            b = s.get(key) or {}
            for fld in BRIEF_FIELDS:
                row[f"{prefix}_{fld}"] = b.get(fld, "")
        per.append(row)

    per.sort(key=lambda r: (r["group"] == "control", -num(r["structure_gain"], -1e9), r["accession"]))
    per_cols = list(per[0].keys()) if per else ["accession"]
    write_tsv(os.path.join(a.outdir, "per_protein_summary.tsv"), per_cols, per)

    qc_cols = list(qc_rows[0].keys()) if qc_rows else ["accession"]
    write_tsv(os.path.join(a.outdir, "plddt_qc_summary.tsv"), qc_cols, qc_rows)

    classes = {}
    for r in per:
        classes[r["structure_added_class"]] = classes.get(r["structure_added_class"], 0) + 1
    run = {
        "database": db,
        "n_proteins": len(qc_rows),
        "n_passed_plddt_gate": sum(1 for r in qc_rows if r.get("qc_flag") == "PASS"),
        "n_searched": len(summaries),
        "n_combined_hit_rows": len(hits),
        "structure_added_class_counts": classes,
    }
    with open(os.path.join(a.outdir, "run_summary.json"), "w") as fh:
        json.dump(run, fh, indent=2)
    print(f"[summarise] proteins={run['n_proteins']} searched={run['n_searched']} hit_rows={len(hits)} classes={classes}")


if __name__ == "__main__":
    main()
