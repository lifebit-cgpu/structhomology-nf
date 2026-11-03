#!/usr/bin/env python3
"""Per-residue pLDDT summary and search gate for one AlphaFold DB model.

pLDDT is read from the AlphaFold DB confidence JSON ('confidenceScore'). If
that file is missing or empty, it is read from the B-factor column of the
model's C-alpha atoms, where AlphaFold stores pLDDT.

Bands follow AlphaFold DB: very high > 90, confident 70-90, low 50-70,
very low <= 50. A model PASSES when mean pLDDT >= --min-mean-plddt AND the
fraction of residues with pLDDT >= 70 is >= --min-frac-confident.

Flags: PASS, FAIL_MEAN, FAIL_FRAC, NO_MODEL.

Standard library only.
"""
import argparse
import json
import os
import statistics


def from_confidence(path):
    try:
        with open(path) as fh:
            data = json.load(fh)
    except (OSError, ValueError):
        return None
    if isinstance(data, list) and data:
        data = data[0]
    if isinstance(data, dict) and data.get("confidenceScore"):
        return [float(x) for x in data["confidenceScore"]]
    return None


def from_pdb(path):
    vals = []
    with open(path, errors="replace") as fh:
        for line in fh:
            if line.startswith("ATOM") and line[12:16].strip() == "CA":
                try:
                    vals.append(float(line[60:66]))
                except ValueError:
                    pass
    return vals


def from_cif(path):
    cols, vals = [], []
    with open(path, errors="replace") as fh:
        for line in fh:
            s = line.strip()
            if s.startswith("_atom_site."):
                cols.append(s.split(".", 1)[1].split()[0])
                continue
            if not cols:
                continue
            if s.startswith("#") or s.startswith("loop_") or s.startswith("_"):
                if vals:
                    break
                continue
            parts = s.split()
            if len(parts) != len(cols):
                continue
            row = dict(zip(cols, parts))
            if row.get("group_PDB") == "ATOM" and row.get("label_atom_id") == "CA":
                try:
                    vals.append(float(row["B_iso_or_equiv"]))
                except (KeyError, ValueError):
                    pass
    return vals


def longest_run(vals, thr=70.0):
    best = cur = 0
    for v in vals:
        cur = cur + 1 if v >= thr else 0
        best = max(best, cur)
    return best


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--accession", required=True)
    ap.add_argument("--group", default="query")
    ap.add_argument("--model", required=True)
    ap.add_argument("--confidence", required=True)
    ap.add_argument("--record", required=True)
    ap.add_argument("--min-mean-plddt", type=float, default=70.0)
    ap.add_argument("--min-frac-confident", type=float, default=0.5)
    ap.add_argument("--out-prefix", required=True)
    a = ap.parse_args()

    try:
        with open(a.record) as fh:
            rec = json.load(fh)
    except (OSError, ValueError):
        rec = {}

    conf = from_confidence(a.confidence)
    model_vals = []
    if os.path.exists(a.model) and os.path.getsize(a.model) > 0:
        model_vals = from_cif(a.model) if a.model.endswith(".cif") else from_pdb(a.model)

    if conf:
        plddt, source = conf, "confidence_json"
    else:
        plddt, source = model_vals, "model_bfactor"

    maxdiff = ""
    if conf and model_vals and len(conf) == len(model_vals):
        maxdiff = round(max(abs(x - y) for x, y in zip(conf, model_vals)), 3)

    n = len(plddt)
    row = {
        "accession": a.accession,
        "group": a.group,
        "gene": rec.get("gene") or "",
        "description": rec.get("uniprot_description") or "",
        "afdb_entry": rec.get("entry_id") or "",
        "model_version": rec.get("model_version") if rec.get("model_version") is not None else "",
        "model_created_date": rec.get("model_created_date") or "",
        "afdb_status": rec.get("status") or "",
        "n_residues": n,
        "mean_plddt": "",
        "median_plddt": "",
        "frac_very_high_gt90": "",
        "frac_confident_70_90": "",
        "frac_low_50_70": "",
        "frac_very_low_le50": "",
        "frac_ge70": "",
        "longest_run_ge70": "",
        "api_mean_plddt": rec.get("api_global_metric_plddt") if rec.get("api_global_metric_plddt") is not None else "",
        "plddt_source": source if n else "",
        "plddt_json_vs_model_maxdiff": maxdiff,
        "min_mean_plddt": a.min_mean_plddt,
        "min_frac_confident": a.min_frac_confident,
        "qc_flag": "",
        "qc_reason": "",
    }

    if n == 0:
        row["qc_flag"] = "NO_MODEL"
        row["qc_reason"] = f"no pLDDT values (AlphaFold DB status: {rec.get('status', 'unknown')})"
    else:
        mean = statistics.fmean(plddt)
        frac = lambda f: round(sum(1 for v in plddt if f(v)) / n, 4)  # noqa: E731
        row.update({
            "mean_plddt": round(mean, 2),
            "median_plddt": round(statistics.median(plddt), 2),
            "frac_very_high_gt90": frac(lambda v: v > 90),
            "frac_confident_70_90": frac(lambda v: 70 < v <= 90),
            "frac_low_50_70": frac(lambda v: 50 < v <= 70),
            "frac_very_low_le50": frac(lambda v: v <= 50),
            "frac_ge70": frac(lambda v: v >= 70),
            "longest_run_ge70": longest_run(plddt),
        })
        if mean < a.min_mean_plddt:
            row["qc_flag"] = "FAIL_MEAN"
            row["qc_reason"] = f"mean pLDDT {mean:.1f} < {a.min_mean_plddt:g}"
        elif row["frac_ge70"] < a.min_frac_confident:
            row["qc_flag"] = "FAIL_FRAC"
            row["qc_reason"] = f"fraction of residues with pLDDT >= 70 is {row['frac_ge70']:.2f} < {a.min_frac_confident:g}"
        else:
            row["qc_flag"] = "PASS"
            row["qc_reason"] = "meets both pLDDT thresholds"

    pre = a.out_prefix
    with open(f"{pre}.plddt_qc.tsv", "w") as fh:
        fh.write("\t".join(row) + "\n")
        fh.write("\t".join(str(v) for v in row.values()) + "\n")
    with open(f"{pre}.plddt_qc.json", "w") as fh:
        json.dump(row, fh, indent=2)
    with open(f"{pre}.qc_flag.txt", "w") as fh:
        fh.write(row["qc_flag"])
    print(f"[plddt_qc] {a.accession}: {row['qc_flag']} mean={row['mean_plddt']} frac_ge70={row['frac_ge70']} source={row['plddt_source']}")


if __name__ == "__main__":
    main()
