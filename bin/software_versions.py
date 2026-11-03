#!/usr/bin/env python3
"""Collate versions.yml files (concatenated) into one YAML and one TSV.

Each input block looks like
    "WORKFLOW:PROCESS":
        tool: version
The output lists each process once (by its short name) plus a 'Workflow'
entry with Nextflow and the pipeline version.
Standard library only (no PyYAML needed).
"""
import argparse
import platform


def parse(path):
    out, proc = {}, None
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            if not line.strip():
                continue
            if not line.startswith((" ", "\t")):
                proc = line.strip().rstrip(":").strip().strip('"').strip("'").split(":")[-1]
                continue
            if proc and ":" in line:
                k, v = line.strip().split(":", 1)
                out.setdefault(proc, {})[k.strip()] = v.strip().strip('"').strip("'")
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--input", required=True)
    ap.add_argument("--nextflow", required=True)
    ap.add_argument("--pipeline", required=True)
    ap.add_argument("--process", default="SOFTWARE_VERSIONS")
    ap.add_argument("--out-prefix", default="software_versions")
    a = ap.parse_args()

    versions = parse(a.input)
    versions[a.process.split(":")[-1]] = {"python": platform.python_version()}
    versions["Workflow"] = {"Nextflow": a.nextflow, "structhomology-nf": a.pipeline}

    with open(f"{a.out_prefix}.yml", "w") as fh:
        for proc in sorted(versions):
            fh.write(f"{proc}:\n")
            for tool, ver in sorted(versions[proc].items()):
                fh.write(f"    {tool}: '{ver}'\n")
    with open(f"{a.out_prefix}.tsv", "w") as fh:
        fh.write("process\ttool\tversion\n")
        for proc in sorted(versions):
            for tool, ver in sorted(versions[proc].items()):
                fh.write(f"{proc}\t{tool}\t{ver}\n")
    tools = sorted({(t, v) for p in versions.values() for t, v in p.items()})
    print("[software_versions] " + ", ".join(f"{t}={v}" for t, v in tools))


if __name__ == "__main__":
    main()
