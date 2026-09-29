#!/usr/bin/env python
"""Enabled residue variable mods from a FragPipe workflow, one per line.

Skips disabled rows, terminal sites, and isobaric-tag masses. Exit 2 if none remain.
"""
from __future__ import annotations

import argparse
import re
import sys

TAG_MASS_PREFIXES = (
    "304.2071",
    "304.2053",
    "229.1629",
    "144.1020",
)


def read_workflow(path: str) -> dict[str, str]:
    out: dict[str, str] = {}
    with open(path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#") or "=" not in line:
                continue
            k, v = line.split("=", 1)
            out[k.strip()] = v.strip()
    return out


def is_tag_mass(mass: str) -> bool:
    return mass.startswith(TAG_MASS_PREFIXES)


def enabled_residue_mods(var_mods: str) -> list[str]:
    found: list[str] = []
    seen: set[str] = set()
    for part in var_mods.split(";"):
        bits = [b.strip() for b in part.split(",")]
        if len(bits) < 3:
            continue
        mass, sites, flag = bits[0], bits[1], bits[2].lower()
        if flag != "true" or is_tag_mass(mass):
            continue
        if not re.fullmatch(r"[A-Z]+", sites):
            continue
        token = sites.upper()
        if token in seen:
            continue
        seen.add(token)
        found.append(token)
    return found


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("workflow", help="FragPipe .workflow path")
    args = ap.parse_args()
    wf = read_workflow(args.workflow)
    mods = enabled_residue_mods(wf.get("msfragger.table.var-mods", ""))
    if not mods:
        print("ERROR: no enabled residue PTM in msfragger.table.var-mods", file=sys.stderr)
        return 2
    print("\n".join(mods))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
