#!/usr/bin/env python3
"""Positive-control acceptance test.

A control PASSES when at least one of its top 10 Foldseek hits is
  (1) a known relative: its PDB ID is in expected_pdb_ids, OR its annotation
      matches family_regex (needed because clustered databases such as PDB100
      keep one representative chain per group, which may be a different
      entry of the same family), AND
  (2) significant: E-value <= --sig-evalue OR max(qTM, tTM) >= --min-tmscore.
Controls not present in the run are reported as NOT_RUN.

Secondary (reported, not used to pass/fail): which expected PDB IDs each arm
recovers among its significant hits, and which only structure recovers.

With --fail-on-failure the script exits 1 when any control FAILS.
Standard library only.
"""
import argparse
import csv
import glob
import os
import re
import sys

TEXT_FIELDS = ["hit_description", "hit_title", "hit_uniprot_name", "hit_pfam", "hit_interpro",
               "hit_cath", "hit_scop", "hit_ecod", "target_header"]


def read_tsv(path):
    with open(path, newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def f(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--combined", required=True)
    ap.add_argument("--all-hits-dir", required=True)
    ap.add_argument("--controls", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--sig-evalue", type=float, default=1e-3)
    ap.add_argument("--min-tmscore", type=float, default=0.5)
    ap.add_argument("--top", type=int, default=10)
    ap.add_argument("--fail-on-failure", action="store_true")
    a = ap.parse_args()

    combined = read_tsv(a.combined) if os.path.getsize(a.combined) > 0 else []
    all_hits = {}
    for p in glob.glob(os.path.join(a.all_hits_dir, "*.all_hits.tsv")):
        all_hits[os.path.basename(p).split(".all_hits.tsv")[0]] = read_tsv(p)
    with open(a.controls, newline="", encoding="utf-8") as fh:
        controls = list(csv.DictReader(fh))

    cols = ["accession", "gene", "fold", "status", "matched_target", "matched_pdb_id", "matched_by",
            "matched_fs_rank", "matched_fs_evalue", "matched_tm_max", "matched_description",
            "n_top_hits_checked", "expected_found_by_structure", "expected_found_by_sequence",
            "expected_found_by_structure_only", "note"]
    out = []
    for c in controls:
        acc = c["accession"].strip()
        expected = {x.strip().upper() for x in re.split(r"[;,\s]+", c.get("expected_pdb_ids", "")) if x.strip()}
        fam = re.compile(c.get("family_regex") or r"(?!x)x", re.IGNORECASE)
        row = {k: "" for k in cols}
        row.update({"accession": acc, "gene": c.get("gene", ""), "fold": c.get("fold", "")})

        if acc not in all_hits:
            row["status"] = "NOT_RUN"
            row["note"] = "control not in this run (not in samplesheet, failed the pLDDT gate, or not searched)"
            out.append(row)
            continue

        top = [h for h in combined if h.get("query") == acc and f(h.get("fs_rank")) is not None
               and f(h["fs_rank"]) <= a.top]
        top.sort(key=lambda h: f(h["fs_rank"]))
        row["n_top_hits_checked"] = len(top)
        match = None
        for h in top:
            pid = (h.get("pdb_id") or h.get("target", "")[:4]).upper()
            text = " ".join(h.get(k, "") for k in TEXT_FIELDS)
            by = "pdb_id" if pid in expected else "family_name" if fam.search(text) else ""
            ev = f(h.get("fs_evalue"))
            tms = [x for x in (f(h.get("fs_qtmscore")), f(h.get("fs_ttmscore"))) if x is not None]
            tm = max(tms) if tms else None
            sig = (ev is not None and ev <= a.sig_evalue) or (tm is not None and tm >= a.min_tmscore)
            if by and sig:
                match = (h, pid, by, ev, tm)
                break
        if match:
            h, pid, by, ev, tm = match
            row.update({
                "status": "PASS", "matched_target": h.get("target", ""), "matched_pdb_id": pid,
                "matched_by": by, "matched_fs_rank": h.get("fs_rank", ""), "matched_fs_evalue": ev,
                "matched_tm_max": tm,
                "matched_description": h.get("hit_description") or h.get("hit_title") or h.get("target_header", ""),
            })
        else:
            row["status"] = "FAIL"
            row["note"] = f"no significant known relative among the top {a.top} Foldseek hits"

        rows = all_hits[acc]
        s_found = sorted({(h.get("pdb_id") or "").upper() for h in rows
                          if h.get("structure_significant") == "true"} & expected)
        q_found = sorted({(h.get("pdb_id") or "").upper() for h in rows
                          if h.get("sequence_significant") == "true"} & expected)
        row["expected_found_by_structure"] = ";".join(s_found)
        row["expected_found_by_sequence"] = ";".join(q_found)
        row["expected_found_by_structure_only"] = ";".join(sorted(set(s_found) - set(q_found)))
        out.append(row)

    with open(a.out, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(out)

    for r in out:
        print(f"[check_controls] {r['accession']} {r['gene']}: {r['status']} "
              f"{r['matched_target']} {r['matched_by']} E={r['matched_fs_evalue']} TM={r['matched_tm_max']}")
    failed = [r["accession"] for r in out if r["status"] == "FAIL"]
    if failed and a.fail_on_failure:
        sys.stderr.write("[check_controls] ERROR: positive controls failed: " + ", ".join(failed) +
                         ". The pipeline did not recover known structural relatives; results are not trustworthy.\n")
        sys.exit(1)


if __name__ == "__main__":
    main()
