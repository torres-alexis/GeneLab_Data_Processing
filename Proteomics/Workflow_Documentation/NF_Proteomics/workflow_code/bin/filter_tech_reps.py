#!/usr/bin/env python
"""First row per dup key (order preserved). LFQ: Factor Value[*]+Bioreplicate; else Sample Name+run. TMT: plex+TechRepMixture+fraction."""

import argparse
import csv
import sys


def _factor_tuple(row, fieldnames):
    return tuple(sorted(
        (k, row[k].strip())
        for k in fieldnames
        if k.startswith("Factor Value[") and (row.get(k) or "").strip()
    ))


def _dedup_key_lfq(row, fieldnames):
    ft = _factor_tuple(row, fieldnames)
    bio = str(row.get("Bioreplicate", "")).strip()
    if ft and bio:
        return ("factor_bio", ft, bio)
    return ("unique", row.get("Sample Name", ""), row.get("run", ""))


def _dedup_key_tmt(row):
    plex = (row.get("plex") or "").strip()
    tr = str(row.get("TechRepMixture", "") or "").strip() or "1"
    fr = row.get("fraction")
    frac = ("" if fr is None else str(fr)).strip()
    return ("tmt", plex, tr, frac)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--mode", choices=("lfq", "tmt"), required=True)
    ap.add_argument("--input", required=True)
    ap.add_argument("--output", required=True)
    args = ap.parse_args()

    with open(args.input, newline="", encoding="utf-8") as f:
        rd = csv.DictReader(f)
        cols = rd.fieldnames or []
        rows = list(rd)

    key_fn = _dedup_key_tmt if args.mode == "tmt" else lambda r: _dedup_key_lfq(r, cols)
    seen = set()
    out = []
    for row in rows:
        k = key_fn(row)
        if k in seen:
            continue
        seen.add(k)
        out.append(row)

    with open(args.output, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=cols, extrasaction="ignore")
        w.writeheader()
        w.writerows(out)

    print(f"filter_tech_reps: {len(rows)} -> {len(out)} ({args.mode})", file=sys.stderr)


if __name__ == "__main__":
    main()
