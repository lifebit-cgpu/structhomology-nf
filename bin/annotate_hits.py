#!/usr/bin/env python3
"""Merge Foldseek (structure) and MMseqs2 (sequence) hits for one query that
searched the SAME target set, classify every target, and annotate.

Significance
  structure arm: E <= --sig-evalue, OR max(qTM, tTM) >= --min-tmscore over an
                 alignment of >= --min-aln-len residues
  sequence arm : E <= --sig-evalue
Category per target: shared (both significant), structure_only,
sequence_only, weak (neither significant).

Annotation (unless --skip-remote): RCSB PDB GraphQL for PDB targets (entry
title, entity description, organism, EC, Pfam, InterPro, CATH, SCOP, ECOD, GO,
UniProt mapping) and UniProt REST for the mapped UniProt entries (protein
name, annotation score, function text). AlphaFold DB targets (AF-<acc>-F1)
are annotated from UniProt directly.

A hit is 'poorly_annotated' when its name says uncharacterized / unknown
function / hypothetical / DUF, or it has no name at all. A hit is
'same_protein' when it maps to the query accession. 'informative' =
significant AND not poorly annotated AND not the query itself: only these can
suggest a function.

Inputs are BLAST-tab files with a header line (--format-mode 4).
Standard library only.
"""
import argparse
import json
import re
import sys
import time
import urllib.error
import urllib.parse
import urllib.request

UA = "structhomology-nf/1.0 (https://github.com/lifebit-cgpu/structhomology-nf)"

POORLY_ANNOTATED = re.compile(
    r"uncharacteri[sz]ed|unknown function|hypothetical protein|\bDUF\d+|unnamed protein"
    r"|putative protein\b|protein of unknown|de novo designed|designed protein",
    re.IGNORECASE,
)

GQL_FULL = """query($ids: [String!]!) {
  entries(entry_ids: $ids) {
    rcsb_id
    struct { title }
    exptl { method }
    rcsb_entry_info { resolution_combined }
    polymer_entities {
      rcsb_polymer_entity { pdbx_description rcsb_ec_lineage { id } }
      rcsb_polymer_entity_container_identifiers { auth_asym_ids uniprot_ids }
      rcsb_entity_source_organism { scientific_name }
      rcsb_polymer_entity_annotation { type name }
    }
  }
}"""

GQL_MIN = """query($ids: [String!]!) {
  entries(entry_ids: $ids) {
    rcsb_id
    struct { title }
    polymer_entities {
      rcsb_polymer_entity { pdbx_description }
      rcsb_polymer_entity_container_identifiers { auth_asym_ids uniprot_ids }
    }
  }
}"""

OUT_COLS = [
    "query", "query_group", "target", "pdb_id", "chain", "hit_kind",
    "in_structure", "in_sequence", "structure_significant", "sequence_significant", "category",
    "fs_rank", "fs_evalue", "fs_bits", "fs_prob", "fs_fident", "fs_alnlen",
    "fs_alntmscore", "fs_qtmscore", "fs_ttmscore", "fs_lddt",
    "fs_qstart", "fs_qend", "fs_tstart", "fs_tend", "fs_qcov", "fs_tcov",
    "seq_rank", "seq_evalue", "seq_bits", "seq_fident", "seq_alnlen",
    "seq_qstart", "seq_qend", "seq_qcov", "seq_tcov",
    "qlen", "tlen", "target_header",
    "hit_title", "hit_description", "hit_organism", "hit_uniprot", "hit_uniprot_name", "hit_gene",
    "hit_annotation_score", "hit_ec", "hit_pfam", "hit_interpro", "hit_cath", "hit_scop", "hit_ecod",
    "hit_go", "hit_function", "experimental_method", "resolution",
    "annotation_source", "same_protein", "poorly_annotated", "informative",
]


def http(url, data=None, headers=None, tries=5, backoff=3.0):
    hdr = {"User-Agent": UA}
    hdr.update(headers or {})
    last = None
    for attempt in range(tries):
        try:
            req = urllib.request.Request(url, data=data, headers=hdr)
            with urllib.request.urlopen(req, timeout=120) as resp:
                body = resp.read()
            if body:
                return body
            last = "empty body"
        except urllib.error.HTTPError as exc:
            if exc.code in (400, 404):
                sys.stderr.write(f"[annotate_hits] HTTP {exc.code} for {url}\n")
                return None
            last = f"HTTP {exc.code}"
        except Exception as exc:
            last = repr(exc)
        time.sleep(backoff * (2 ** attempt))
    sys.stderr.write(f"[annotate_hits] WARNING: giving up on {url}: {last}\n")
    return None


def as_float(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


def read_m8(path):
    """Read a BLAST-tab file with a header line; keep the best alignment
    (lowest E-value, then highest bit score) per target."""
    best = {}
    with open(path, encoding="utf-8", errors="replace") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        if header == [""] or "target" not in header:
            return best
        nf = len(header)
        hidx = header.index("theader") if "theader" in header else None
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) > nf and hidx is not None:
                extra = len(parts) - nf
                parts = parts[:hidx] + [" ".join(parts[hidx:hidx + extra + 1])] + parts[hidx + extra + 1:]
            if len(parts) != nf:
                continue
            row = dict(zip(header, parts))
            ev = as_float(row.get("evalue"))
            row["_ev"] = ev if ev is not None else float("inf")
            row["_bits"] = as_float(row.get("bits")) or 0.0
            cur = best.get(row["target"])
            if cur is None or (row["_ev"], -row["_bits"]) < (cur["_ev"], -cur["_bits"]):
                best[row["target"]] = row
    return best


def ranks(hits):
    order = sorted(hits, key=lambda t: (hits[t]["_ev"], -hits[t]["_bits"], t))
    return {t: i + 1 for i, t in enumerate(order)}


def parse_target(t):
    m = re.match(r"^AF-([A-Z0-9]+)-F\d+", t)
    if m:
        return {"hit_kind": "afdb", "uniprot": m.group(1), "pdb_id": "", "chain": ""}
    m = re.match(r"^([0-9][A-Za-z0-9]{3})(?:[-_.].*)?$", t)
    if m:
        chain = t.rsplit("_", 1)[1] if "_" in t else ""
        return {"hit_kind": "pdb", "uniprot": "", "pdb_id": m.group(1).upper(), "chain": chain.split("-")[0]}
    return {"hit_kind": "other", "uniprot": "", "pdb_id": "", "chain": ""}


def header_text(target, header):
    h = (header or "").strip()
    if h.startswith(target):
        h = h[len(target):].strip()
    return h


def rcsb_annotate(pdb_ids, url):
    out = {}
    ids = sorted(set(pdb_ids))
    for i in range(0, len(ids), 50):
        chunk = ids[i:i + 50]
        for query in (GQL_FULL, GQL_MIN):
            payload = json.dumps({"query": query, "variables": {"ids": chunk}}).encode()
            body = http(url, data=payload, headers={"Content-Type": "application/json"})
            if not body:
                continue
            try:
                js = json.loads(body)
            except ValueError:
                continue
            entries = (js.get("data") or {}).get("entries")
            if entries is None:
                sys.stderr.write(f"[annotate_hits] RCSB GraphQL errors: {str(js.get('errors'))[:300]}\n")
                continue
            for e in entries:
                if e and e.get("rcsb_id"):
                    out[e["rcsb_id"].upper()] = e
            break
    return out


def pick_entity(entry, chain):
    ents = entry.get("polymer_entities") or []
    for e in ents:
        ids = (e.get("rcsb_polymer_entity_container_identifiers") or {}).get("auth_asym_ids") or []
        if chain and chain in ids:
            return e
    return ents[0] if ents else {}


def entity_fields(entry, ent):
    pe = ent.get("rcsb_polymer_entity") or {}
    ci = ent.get("rcsb_polymer_entity_container_identifiers") or {}
    org = ent.get("rcsb_entity_source_organism") or []
    ann = ent.get("rcsb_polymer_entity_annotation") or []

    def names(types, n=4):
        seen = []
        for x in ann:
            nm = x.get("name")
            if x.get("type") in types and nm and nm not in seen:
                seen.append(nm)
        return "; ".join(seen[:n])

    ec_ids = {x.get("id") for x in (pe.get("rcsb_ec_lineage") or []) if x.get("id")}
    ec_leaf = sorted(i for i in ec_ids if not any(j != i and j.startswith(i + ".") for j in ec_ids))
    res = (entry.get("rcsb_entry_info") or {}).get("resolution_combined") or []
    return {
        "hit_title": (entry.get("struct") or {}).get("title") or "",
        "hit_description": pe.get("pdbx_description") or "",
        "hit_organism": "; ".join(sorted({o.get("scientific_name") for o in org if o.get("scientific_name")})),
        "hit_uniprot": ";".join(ci.get("uniprot_ids") or []),
        "hit_ec": "; ".join(ec_leaf),
        "hit_pfam": names({"Pfam"}),
        "hit_interpro": names({"InterPro"}),
        "hit_cath": names({"CATH"}),
        "hit_scop": names({"SCOP", "SCOP2", "SCOP2B"}),
        "hit_ecod": names({"ECOD"}),
        "hit_go": names({"GO"}, 6),
        "experimental_method": "; ".join(x.get("method", "") for x in (entry.get("exptl") or [])),
        "resolution": res[0] if res else "",
    }


def uniprot_annotate(accs, base):
    out = {}
    accs = sorted({a for a in accs if a})
    fields = "accession,protein_name,gene_primary,organism_name,annotation_score,ec,cc_function"

    def query(chunk):
        url = f"{base.rstrip('/')}/uniprotkb/accessions?" + urllib.parse.urlencode(
            {"accessions": ",".join(chunk), "fields": fields, "format": "tsv"})
        return http(url)

    for i in range(0, len(accs), 100):
        chunk = accs[i:i + 100]
        bodies = []
        body = query(chunk)
        if body:
            bodies.append(body)
        elif len(chunk) > 1:  # one bad accession fails the batch; retry one by one
            bodies += [b for b in (query([x]) for x in chunk) if b]
        for b in bodies:
            for line in b.decode("utf-8", "replace").splitlines()[1:]:
                p = line.split("\t")
                p += [""] * (7 - len(p))
                out[p[0]] = {
                    "hit_uniprot_name": p[1],
                    "hit_gene": p[2],
                    "hit_uniprot_organism": p[3],
                    "hit_annotation_score": p[4],
                    "hit_uniprot_ec": p[5],
                    "hit_function": re.sub(r"^FUNCTION:\s*", "", p[6])[:400],
                }
    return out


def fmt(v):
    if isinstance(v, bool):
        return "true" if v else "false"
    if v is None:
        return ""
    return str(v).replace("\t", " ").replace("\n", " ")


def brief(r):
    if r is None:
        return None
    tm = [x for x in (as_float(r.get("fs_qtmscore")), as_float(r.get("fs_ttmscore"))) if x is not None]
    return {
        "target": r["target"],
        "pdb_id": r["pdb_id"],
        "description": r["hit_description"] or r["hit_uniprot_name"] or r["hit_title"] or r["target_header"],
        "organism": r["hit_organism"],
        "pfam": r["hit_pfam"],
        "ec": r["hit_ec"],
        "category": r["category"],
        "fs_evalue": as_float(r["fs_evalue"]),
        "fs_prob": as_float(r["fs_prob"]),
        "fs_tm_max": max(tm) if tm else None,
        "fs_fident": as_float(r["fs_fident"]),
        "seq_evalue": as_float(r["seq_evalue"]),
        "seq_fident": as_float(r["seq_fident"]),
        "informative": r["informative"],
        "same_protein": r["same_protein"],
    }


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--accession", required=True)
    ap.add_argument("--group", default="query")
    ap.add_argument("--foldseek", required=True)
    ap.add_argument("--mmseqs", required=True)
    ap.add_argument("--qc", required=True)
    ap.add_argument("--out-prefix", required=True)
    ap.add_argument("--top-n", type=int, default=10)
    ap.add_argument("--sig-evalue", type=float, default=1e-3)
    ap.add_argument("--min-tmscore", type=float, default=0.5)
    ap.add_argument("--min-aln-len", type=int, default=50)
    ap.add_argument("--max-annotate", type=int, default=500,
                    help="annotate at most this many targets (best-ranked first)")
    ap.add_argument("--skip-remote", action="store_true")
    ap.add_argument("--rcsb-url", default="https://data.rcsb.org/graphql")
    ap.add_argument("--uniprot-url", default="https://rest.uniprot.org")
    a = ap.parse_args()
    acc = a.accession

    try:
        with open(a.qc) as fh:
            qc = json.load(fh)
    except (OSError, ValueError):
        qc = {}

    fs, ms = read_m8(a.foldseek), read_m8(a.mmseqs)
    fs_rank, ms_rank = ranks(fs), ranks(ms)

    rows = []
    for t in set(fs) | set(ms):
        f, m = fs.get(t, {}), ms.get(t, {})
        info = parse_target(t)
        qtm, ttm = as_float(f.get("qtmscore")), as_float(f.get("ttmscore"))
        tms = [x for x in (qtm, ttm) if x is not None]
        tm_max = max(tms) if tms else None
        fev, mev = as_float(f.get("evalue")), as_float(m.get("evalue"))
        alnlen = as_float(f.get("alnlen")) or 0
        s_sig = bool(f) and ((fev is not None and fev <= a.sig_evalue) or
                             (tm_max is not None and tm_max >= a.min_tmscore and alnlen >= a.min_aln_len))
        q_sig = bool(m) and mev is not None and mev <= a.sig_evalue
        category = "shared" if (s_sig and q_sig) else "structure_only" if s_sig else \
            "sequence_only" if q_sig else "weak"
        hdr = f.get("theader") or m.get("theader") or ""
        row = {c: "" for c in OUT_COLS}
        row.update({
            "query": acc, "query_group": a.group, "target": t,
            "pdb_id": info["pdb_id"], "chain": info["chain"], "hit_kind": info["hit_kind"],
            "in_structure": bool(f), "in_sequence": bool(m),
            "structure_significant": s_sig, "sequence_significant": q_sig, "category": category,
            "fs_rank": fs_rank.get(t, ""), "seq_rank": ms_rank.get(t, ""),
            "qlen": f.get("qlen") or m.get("qlen") or "", "tlen": f.get("tlen") or m.get("tlen") or "",
            "target_header": header_text(t, hdr),
            "hit_uniprot": info["uniprot"],
            "annotation_source": "database_header",
            "_rank": min(fs_rank.get(t, 10 ** 9), ms_rank.get(t, 10 ** 9)),
        })
        for k in ("evalue", "bits", "prob", "fident", "alnlen", "alntmscore", "qtmscore", "ttmscore",
                  "lddt", "qstart", "qend", "tstart", "tend", "qcov", "tcov"):
            row[f"fs_{k}"] = f.get(k, "")
        for k in ("evalue", "bits", "fident", "alnlen", "qstart", "qend", "qcov", "tcov"):
            row[f"seq_{k}"] = m.get(k, "")
        rows.append(row)

    rows.sort(key=lambda r: (r["_rank"], r["target"]))

    # Annotate: top-N of each arm plus every significant hit, best-ranked first, capped
    want = [r for r in rows if r["structure_significant"] or r["sequence_significant"]
            or (r["fs_rank"] != "" and r["fs_rank"] <= a.top_n)
            or (r["seq_rank"] != "" and r["seq_rank"] <= a.top_n)][:a.max_annotate]
    if not a.skip_remote and want:
        pdb = rcsb_annotate([r["pdb_id"] for r in want if r["pdb_id"]], a.rcsb_url)
        for r in want:
            e = pdb.get(r["pdb_id"])
            if e:
                r.update(entity_fields(e, pick_entity(e, r["chain"])))
                r["annotation_source"] = "rcsb"
        up_ids = {r["hit_uniprot"].split(";")[0] for r in want if r["hit_uniprot"]}
        up = uniprot_annotate(up_ids, a.uniprot_url)
        for r in want:
            u = up.get(r["hit_uniprot"].split(";")[0]) if r["hit_uniprot"] else None
            if u:
                r["hit_uniprot_name"] = u["hit_uniprot_name"]
                r["hit_gene"] = u["hit_gene"]
                r["hit_annotation_score"] = u["hit_annotation_score"]
                r["hit_function"] = u["hit_function"]
                r["hit_ec"] = r["hit_ec"] or u["hit_uniprot_ec"]
                r["hit_organism"] = r["hit_organism"] or u["hit_uniprot_organism"]
                r["annotation_source"] = "rcsb+uniprot" if r["annotation_source"] == "rcsb" else "uniprot"

    for r in rows:
        desc = r["hit_description"] or r["hit_uniprot_name"] or r["hit_title"] or r["target_header"]
        score = as_float(r["hit_annotation_score"])
        poor = (not desc.strip()) or bool(POORLY_ANNOTATED.search(desc)) or (
            bool(r["hit_uniprot_name"]) and bool(POORLY_ANNOTATED.search(r["hit_uniprot_name"])))
        if score is not None and score <= 1 and not r["hit_ec"] and not r["hit_function"]:
            poor = True
        same = acc in [x.strip() for x in r["hit_uniprot"].split(";") if x.strip()]
        r["poorly_annotated"] = poor
        r["same_protein"] = same
        r["informative"] = (r["structure_significant"] or r["sequence_significant"]) and not poor and not same

    pre = a.out_prefix

    def write(path, subset):
        with open(path, "w") as fh:
            fh.write("\t".join(OUT_COLS) + "\n")
            for r in subset:
                fh.write("\t".join(fmt(r[c]) for c in OUT_COLS) + "\n")

    write(f"{pre}.all_hits.tsv", rows)
    top = [r for r in rows if (r["fs_rank"] != "" and r["fs_rank"] <= a.top_n)
           or (r["seq_rank"] != "" and r["seq_rank"] <= a.top_n)]
    write(f"{pre}.top_hits.tsv", top)

    def uniq(sel):
        return len({r["pdb_id"] or r["target"] for r in sel})

    s_sig = [r for r in rows if r["structure_significant"]]
    q_sig = [r for r in rows if r["sequence_significant"]]
    s_inf = [r for r in s_sig if r["informative"]]
    q_inf = [r for r in q_sig if r["informative"]]
    so_inf = [r for r in s_inf if r["category"] == "structure_only"]
    qo_inf = [r for r in q_inf if r["category"] == "sequence_only"]
    counts = {
        "n_structure_hits": len(fs), "n_sequence_hits": len(ms),
        "n_structure_significant": len(s_sig), "n_sequence_significant": len(q_sig),
        "n_shared_significant": sum(1 for r in rows if r["category"] == "shared"),
        "n_structure_only_significant": sum(1 for r in rows if r["category"] == "structure_only"),
        "n_sequence_only_significant": sum(1 for r in rows if r["category"] == "sequence_only"),
        "n_structure_informative_pdb": uniq(s_inf),
        "n_sequence_informative_pdb": uniq(q_inf),
        "n_structure_only_informative_pdb": uniq(so_inf),
        "n_sequence_only_informative_pdb": uniq(qo_inf),
        "n_structure_significant_poorly_annotated": sum(1 for r in s_sig if r["poorly_annotated"]),
    }
    if not s_sig and not q_sig:
        cls = "no_significant_hits"
    elif not s_inf and not q_inf:
        cls = "only_self_or_poorly_annotated_hits"
    elif s_inf and not q_inf:
        cls = "structure_only"
    elif so_inf:
        cls = "structure_extends_sequence"
    elif s_inf:
        cls = "structure_agrees_with_sequence"
    else:
        cls = "sequence_only"

    def best(sel, key):
        sel = [r for r in sel if r[key] != ""]
        return brief(min(sel, key=lambda r: r[key])) if sel else None

    summary = {
        "accession": acc,
        "group": a.group,
        "gene": qc.get("gene", ""),
        "description": qc.get("description", ""),
        "model_version": qc.get("model_version", ""),
        "mean_plddt": qc.get("mean_plddt", ""),
        "frac_ge70": qc.get("frac_ge70", ""),
        "qc_flag": qc.get("qc_flag", ""),
        "thresholds": {"sig_evalue": a.sig_evalue, "min_tmscore": a.min_tmscore, "min_aln_len": a.min_aln_len},
        "counts": counts,
        "structure_added_class": cls,
        "structure_gain": counts["n_structure_only_informative_pdb"] - counts["n_sequence_only_informative_pdb"],
        "same_protein_in_db": any(r["same_protein"] for r in rows if r["structure_significant"] or r["sequence_significant"]),
        "top_structure_hit": best(rows, "fs_rank"),
        "top_sequence_hit": best(rows, "seq_rank"),
        "best_structure_informative_hit": best(s_inf, "fs_rank"),
        "best_structure_only_informative_hit": best(so_inf, "fs_rank"),
        "best_sequence_informative_hit": best(q_inf, "seq_rank"),
    }
    with open(f"{pre}.summary.json", "w") as fh:
        json.dump(summary, fh, indent=2)
    print(f"[annotate_hits] {acc}: structure {len(fs)} hits ({len(s_sig)} sig), sequence {len(ms)} hits "
          f"({len(q_sig)} sig), class={cls}")


if __name__ == "__main__":
    main()
