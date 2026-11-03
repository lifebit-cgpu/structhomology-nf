#!/usr/bin/env python3
"""Print the parameter table of nextflow_schema.json as Markdown.

Usage: python3 bin/schema_to_markdown.py nextflow_schema.json > params.md
The README parameter table is generated with this script.
"""
import json
import sys


def main():
    path = sys.argv[1] if len(sys.argv) > 1 else "nextflow_schema.json"
    with open(path) as fh:
        schema = json.load(fh)
    for ref in schema.get("allOf", []):
        key = ref["$ref"].split("/")[-1]
        group = schema["$defs"][key]
        required = set(group.get("required", []))
        print(f"**{group.get('title', key)}**\n")
        print("| Parameter | Type | Default | Description |")
        print("|---|---|---|---|")
        for name, p in group.get("properties", {}).items():
            typ = p.get("type", "")
            if p.get("format"):
                typ += f" ({p['format']})"
            default = p.get("default", "required" if name in required else "")
            if isinstance(default, bool):
                default = str(default).lower()
            desc = p.get("description", "").replace("|", "\\|")
            print(f"| `--{name}` | {typ} | `{default}` | {desc} |" if default != "" else
                  f"| `--{name}` | {typ} | | {desc} |")
        print()


if __name__ == "__main__":
    main()
