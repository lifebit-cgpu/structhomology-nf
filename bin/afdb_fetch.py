#!/usr/bin/env python3
"""Fetch the AlphaFold DB record, model, confidence file and sequence for one
UniProt accession.

The model version is resolved from the AlphaFold DB API ('latestVersion' of the
canonical record) and is never hard-coded. When the API returns several models
(for example isoform models such as Q8IXR9-2), the record whose
'uniprotAccession' equals the query exactly - the canonical AF-<acc>-F1 model -
is used.

Outputs (prefix = --out-prefix):
  <prefix>.model.pdb (or .cif)   the model; empty file if AFDB has none
  <prefix>.confidence.json       per-residue pLDDT ({} if unavailable)
  <prefix>.fasta                 sequence of the modelled region
  <prefix>.record.json           API record + resolved version + status

Standard library only.
"""
import argparse
import datetime as dt
import json
import sys
import time
import urllib.error
import urllib.request

UA = "structhomology-nf/1.0 (https://github.com/lifebit-cgpu/structhomology-nf)"


class FetchError(Exception):
    pass


def fetch(url, tries=6, backoff=5.0):
    """GET with retry and exponential backoff. Returns bytes, or None on 404.
    Retries on HTTP 429/5xx, network errors and empty bodies, which the
    AlphaFold file server returns under load."""
    last = None
    for attempt in range(tries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": UA})
            with urllib.request.urlopen(req, timeout=180) as resp:
                body = resp.read()
            if body:
                return body
            last = "empty response body"
        except urllib.error.HTTPError as exc:
            if exc.code == 404:
                return None
            last = f"HTTP {exc.code}"
        except Exception as exc:  # network errors, timeouts
            last = repr(exc)
        wait = backoff * (2 ** attempt)
        sys.stderr.write(f"[afdb_fetch] {url}: {last}; retry {attempt + 1}/{tries} in {wait:.0f}s\n")
        time.sleep(wait)
    raise FetchError(f"could not fetch {url} after {tries} attempts ({last})")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--accession", required=True)
    ap.add_argument("--api", default="https://alphafold.ebi.ac.uk/api/prediction")
    ap.add_argument("--out-prefix", required=True)
    a = ap.parse_args()

    acc, pre = a.accession.strip(), a.out_prefix
    record = {
        "accession": acc,
        "api_url": f"{a.api.rstrip('/')}/{acc}",
        "fetched_at_utc": dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds"),
    }

    try:
        body = fetch(record["api_url"])
    except FetchError as exc:
        sys.exit(f"[afdb_fetch] ERROR: {exc}")
    entries = json.loads(body) if body else []
    if isinstance(entries, dict):
        entries = [entries]
    record["n_models_returned"] = len(entries)
    record["models_returned"] = [e.get("entryId") for e in entries]

    canonical = [e for e in entries if e.get("uniprotAccession") == acc]
    pick = next((e for e in canonical if e.get("entryId") == f"AF-{acc}-F1"),
                canonical[0] if canonical else None)

    model_name = f"{pre}.model.pdb"
    seq = ""
    conf = {}
    if pick is None:
        record["status"] = "no_afdb_model" if not entries else "no_canonical_model"
        open(model_name, "wb").close()
    else:
        version = pick.get("latestVersion")
        record.update({
            "status": "ok",
            "entry_id": pick.get("entryId"),
            "model_version": version,
            "model_created_date": pick.get("modelCreatedDate"),
            "api_global_metric_plddt": pick.get("globalMetricValue"),
            "uniprot_start": pick.get("uniprotStart"),
            "uniprot_end": pick.get("uniprotEnd"),
            "uniprot_description": pick.get("uniprotDescription"),
            "gene": pick.get("gene"),
            "organism": pick.get("organismScientificName"),
            "sequence_version_date": pick.get("sequenceVersionDate"),
        })
        model_url = pick.get("pdbUrl") or pick.get("cifUrl")
        if not model_url:
            sys.exit(f"[afdb_fetch] ERROR: record {pick.get('entryId')} has neither pdbUrl nor cifUrl")
        if model_url.endswith(".cif"):
            model_name = f"{pre}.model.cif"
        if version is not None and f"_v{version}." not in model_url:
            record["version_warning"] = f"model URL does not carry _v{version}: {model_url}"
        record["model_url"] = model_url
        try:
            data = fetch(model_url)
        except FetchError as exc:
            sys.exit(f"[afdb_fetch] ERROR: {exc}")
        if data is None:
            sys.exit(f"[afdb_fetch] ERROR: model file not found (404): {model_url}")
        with open(model_name, "wb") as fh:
            fh.write(data)

        conf_url = pick.get("plddtDocUrl") or \
            f"https://alphafold.ebi.ac.uk/files/{pick.get('entryId')}-confidence_v{version}.json"
        record["confidence_url"] = conf_url
        try:
            cbody = fetch(conf_url, tries=4)
            conf = json.loads(cbody) if cbody else {}
        except (FetchError, ValueError) as exc:
            record["confidence_warning"] = f"confidence JSON unavailable ({exc}); pLDDT will be read from model B-factors"
            conf = {}

        seq = pick.get("uniprotSequence") or pick.get("sequence") or ""
        start, end = pick.get("uniprotStart"), pick.get("uniprotEnd")
        if seq and isinstance(start, int) and isinstance(end, int) and (start > 1 or end < len(seq)):
            seq = seq[start - 1:end]
        record["model_sequence_length"] = len(seq)

    with open(f"{pre}.confidence.json", "w") as fh:
        json.dump(conf, fh)
    with open(f"{pre}.fasta", "w") as fh:
        if seq:
            fh.write(f">{acc} {record.get('entry_id', '')} v{record.get('model_version', '')}\n")
            for i in range(0, len(seq), 60):
                fh.write(seq[i:i + 60] + "\n")
    with open(f"{pre}.record.json", "w") as fh:
        json.dump(record, fh, indent=2)

    print(f"[afdb_fetch] {acc}: status={record['status']} entry={record.get('entry_id')} "
          f"version=v{record.get('model_version')} models_returned={record['n_models_returned']}")


if __name__ == "__main__":
    main()
