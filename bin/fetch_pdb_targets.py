#!/usr/bin/env python3
"""Download mmCIF files for a list of PDB IDs from RCSB.

Input: a text file with PDB IDs separated by newlines, spaces or commas;
lines starting with '#' are comments. Writes <outdir>/<id>.cif.gz and a
manifest. Exits non-zero if any ID is malformed or cannot be downloaded, so a
test database is never silently incomplete.

Standard library only.
"""
import argparse
import os
import re
import sys
import time
import urllib.error
import urllib.request

UA = "structhomology-nf/1.0 (https://github.com/lifebit-cgpu/structhomology-nf)"


def fetch(url, tries=5, backoff=3.0):
    last = None
    for attempt in range(tries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": UA})
            with urllib.request.urlopen(req, timeout=120) as resp:
                body = resp.read()
            if body:
                return body
            last = "empty body"
        except urllib.error.HTTPError as exc:
            if exc.code == 404:
                return None
            last = f"HTTP {exc.code}"
        except Exception as exc:
            last = repr(exc)
        time.sleep(backoff * (2 ** attempt))
    sys.stderr.write(f"[fetch_pdb_targets] giving up on {url}: {last}\n")
    return None


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--ids", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--manifest", required=True)
    ap.add_argument("--base-url", default="https://files.rcsb.org/download")
    a = ap.parse_args()

    ids = []
    with open(a.ids) as fh:
        for line in fh:
            line = line.split("#", 1)[0]
            ids += [t for t in re.split(r"[\s,]+", line) if t]
    bad = [i for i in ids if not re.fullmatch(r"[0-9][A-Za-z0-9]{3}", i)]
    if bad:
        sys.exit(f"[fetch_pdb_targets] ERROR: not 4-character PDB IDs: {', '.join(bad)}")
    ids = sorted({i.upper() for i in ids})
    if not ids:
        sys.exit("[fetch_pdb_targets] ERROR: no PDB IDs in " + a.ids)

    os.makedirs(a.outdir, exist_ok=True)
    rows, failed = [], []
    for pid in ids:
        url = f"{a.base_url.rstrip('/')}/{pid}.cif.gz"
        data = fetch(url)
        if not data:
            failed.append(pid)
            rows.append((pid, url, "failed", 0))
            continue
        out = os.path.join(a.outdir, f"{pid.lower()}.cif.gz")
        with open(out, "wb") as fh:
            fh.write(data)
        rows.append((pid, url, "ok", len(data)))

    with open(a.manifest, "w") as fh:
        fh.write("pdb_id\turl\tstatus\tbytes\n")
        for r in rows:
            fh.write("\t".join(str(x) for x in r) + "\n")
    print(f"[fetch_pdb_targets] downloaded {len(ids) - len(failed)}/{len(ids)} entries")
    if failed:
        sys.exit(f"[fetch_pdb_targets] ERROR: could not download {', '.join(failed)}")


if __name__ == "__main__":
    main()
