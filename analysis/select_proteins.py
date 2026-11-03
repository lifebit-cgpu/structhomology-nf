#!/usr/bin/env python3
"""select_proteins.py - select reviewed human 'Uncharacterized' proteins with a confident AlphaFold DB model.
Re-implements the selection logic used for the structhomology-nf analysis run of 2026-11-03.
  python select_proteins.py --count-only          # how many entries the UniProt query returns today
  python select_proteins.py --outdir selection    # full selection (~480 AFDB API + confidence downloads, serial)
Requires python>=3.9, requests, pandas."""
import argparse, re, time, sys, os, requests, pandas as pd
QUERY = "reviewed:true AND organism_id:9606 AND protein_name:uncharacterized"
UNIPROT = "https://rest.uniprot.org/uniprotkb/search"
FIELDS = "accession,gene_names,protein_name,length,annotation_score,protein_existence"
AFDB_API = "https://alphafold.ebi.ac.uk/api/prediction/{acc}"
MIN_MEAN, MIN_FRAC = 70.0, 0.5          # pLDDT gate, same defaults as the pipeline
BANDS = [(100,200),(200,300),(300,450),(450,650),(650,1001)]; PER_BAND, N_TOTAL = 4, 20
FAMILY_WORDS = re.compile(r"domain|family|kinase|homolog|-like|containing", re.I)

def get(url, params=None, tries=8):
    for i in range(tries):
        try:
            r = requests.get(url, params=params, timeout=60)
            if r.status_code == 200 and r.content: return r
            if r.status_code == 404: return None
        except requests.RequestException: pass
        time.sleep(min(60, 2 ** i))            # AFDB throttles (HTTP 429 / empty bodies): back off
    raise RuntimeError(f"failed: {url}")

def uniprot_all():
    rows, url, params = [], UNIPROT, {"query": QUERY, "fields": FIELDS, "format": "tsv", "size": 500}
    while url:
        r = get(url, params); lines = r.text.strip().split("\n"); hdr = lines[0].split("\t")
        rows += [dict(zip(hdr, l.split("\t"))) for l in lines[1:]]
        m = re.search(r'<([^>]+)>;\s*rel="next"', r.headers.get("Link", "")); url, params = (m.group(1), None) if m else (None, None)
    return pd.DataFrame(rows)

def afdb(acc):
    recs = get(AFDB_API.format(acc=acc))
    if recs is None: return None
    recs = [x for x in recs.json() if x.get("uniprotAccession") == acc]   # canonical AF-<acc>-F1, not isoform models
    if not recs: return None
    rec = recs[0]; conf = get(rec["plddtDocUrl"]).json(); conf = conf[0] if isinstance(conf, list) else conf
    s = pd.Series(conf["confidenceScore"], dtype=float)
    return dict(afdb_version=rec.get("latestVersion"), entry_id=rec.get("entryId"), api_mean_plddt=rec.get("globalMetricValue"),
                model_start=rec.get("uniprotStart"), model_end=rec.get("uniprotEnd"), cif_url=rec.get("cifUrl"),
                mean_plddt=round(s.mean(), 2), frac_very_high=(s > 90).mean(), frac_confident=((s >= 70) & (s <= 90)).mean(),
                frac_low=((s >= 50) & (s < 70)).mean(), frac_very_low=(s < 50).mean(), frac_ge70=(s >= 70).mean())

def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--outdir", default="selection"); ap.add_argument("--count-only", action="store_true")
    a = ap.parse_args()
    if a.count_only:
        r = get(UNIPROT, {"query": QUERY, "format": "tsv", "fields": "accession", "size": 1})
        print("UniProt query:", QUERY); print("x-total-results:", r.headers.get("x-total-results")); return
    os.makedirs(a.outdir, exist_ok=True)
    up = uniprot_all(); print("UniProt entries:", len(up))
    up = up[up["Protein names"].str.contains(r"^(Putative )?[Uu]ncharacterized")].copy(); print("named Uncharacterized:", len(up))
    out = []
    for i, acc in enumerate(up.Entry):
        out.append({"accession": acc, **(afdb(acc) or {})}); time.sleep(0.5)
        if i % 50 == 0: print(i, acc, file=sys.stderr)
    df = up.rename(columns={"Entry": "accession", "Gene Names": "gene", "Protein names": "protein_name", "Length": "length",
                            "Annotation": "annotation_score", "Protein existence": "protein_existence"}).merge(pd.DataFrame(out), on="accession")
    df["length"] = df.length.astype(int); df["annotation_score"] = df.annotation_score.astype(float)
    df.to_csv(f"{a.outdir}/candidates_with_plddt.csv", index=False)
    df["reason"] = ""
    df.loc[df.mean_plddt.isna(), "reason"] = "no AFDB model"
    df.loc[(df.reason == "") & ((df.mean_plddt < MIN_MEAN) | (df.frac_ge70 < MIN_FRAC)), "reason"] = "fails pLDDT gate"
    df.loc[(df.reason == "") & ((df.length < 100) | (df.length > 1000)), "reason"] = "length outside 100-1000"
    df.loc[(df.reason == "") & ~((df.model_start == 1) & (df.model_end == df.length)), "reason"] = "model does not cover full protein"
    df.loc[(df.reason == "") & df.protein_name.str.contains(FAMILY_WORDS), "reason"] = "family-implying name"
    el = df[df.reason == ""].copy(); print("eligible:", len(el))
    el["putative"] = el.protein_name.str.startswith("Putative")
    el = el.sort_values(["putative", "annotation_score", "mean_plddt"], ascending=[True, True, False])
    pick = []
    for lo, hi in BANDS: pick += list(el[(el.length >= lo) & (el.length < hi)].accession[:PER_BAND])
    pick += [x for x in el.accession if x not in pick][:max(0, N_TOTAL - len(pick))]
    df.loc[(df.reason == "") & ~df.accession.isin(pick), "reason"] = "eligible, not picked (length spread)"
    s = df[df.accession.isin(pick)].copy(); s["uniprot_query"] = QUERY; s["access_date"] = time.strftime("%Y-%m-%d")
    s.to_csv(f"{a.outdir}/uncharacterized_selection.csv", index=False); df[df.reason != ""].to_csv(f"{a.outdir}/uncharacterized_rejected.csv", index=False)
    pd.DataFrame({"accession": s.accession, "group": "query"}).to_csv(f"{a.outdir}/samplesheet_queries.csv", index=False)
    print("selected:", len(s)); print(",".join(s.accession))

if __name__ == "__main__": main()
