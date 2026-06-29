#!/usr/bin/env python
"""Keep first row per technical-replicate group (order preserved).

Has Tech Reps column optional:
  Missing column, blank, or FALSE: row is kept (unique per Sample Name + run [LFQ] or run [TMT]).
  TRUE: collapse with other TRUE rows sharing:
    LFQ: Source Name + Factor Value[*]
    TMT: plex + TechRepMixture + fraction
"""

import argparse
import csv
import shutil
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from tech_rep_utils import parse_has_tech_reps


def _factor_tuple(row, fieldnames):
    return tuple(sorted(
        (k, row[k].strip())
        for k in fieldnames
        if k.startswith("Factor Value[") and (row.get(k) or "").strip()
    ))


def _has_column(fieldnames, name: str) -> bool:
    return name in (fieldnames or [])


def _row_collapse_tech_reps(row, fieldnames) -> bool:
    """TRUE only when column present and cell is TRUE; otherwise keep row."""
    if not _has_column(fieldnames, "Has Tech Reps"):
        return False
    return parse_has_tech_reps(row.get("Has Tech Reps")) is True


def _unique_key_lfq(row):
    run = (row.get("run") or "").strip()
    return ("unique", row.get("Sample Name", ""), run)


def _unique_key_tmt(row):
    run = (row.get("run") or "").strip()
    return ("unique", run or row.get("Sample Name", ""))


def _dedup_key_lfq(row, fieldnames):
    if not _row_collapse_tech_reps(row, fieldnames):
        return _unique_key_lfq(row)

    source = (row.get("Source Name") or "").strip()
    if not source:
        sys.exit(
            f"Error: Has Tech Reps=TRUE requires Source Name (Sample Name={row.get('Sample Name', '')})"
        )
    return ("source_tech", source, _factor_tuple(row, fieldnames))


def _dedup_key_tmt(row, fieldnames):
    if not _row_collapse_tech_reps(row, fieldnames):
        return _unique_key_tmt(row)

    plex = (row.get("plex") or "").strip()
    tr = str(row.get("TechRepMixture", "") or "").strip() or "1"
    fr = row.get("fraction")
    frac = ("" if fr is None else str(fr)).strip()
    return ("tmt_tech", plex, tr, frac)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--mode", choices=("lfq", "tmt"), required=True)
    ap.add_argument("--input", required=True)
    ap.add_argument("--output", required=True)
    ap.add_argument(
        "--publish-if-changed",
        default="",
        help="Copy output to this path only when row count decreases.",
    )
    args = ap.parse_args()

    with open(args.input, newline="", encoding="utf-8") as f:
        rd = csv.DictReader(f)
        cols = rd.fieldnames or []
        rows = list(rd)

    key_fn = (
        (lambda r: _dedup_key_tmt(r, cols))
        if args.mode == "tmt"
        else (lambda r: _dedup_key_lfq(r, cols))
    )
    seen = set()
    out = []
    for row in rows:
        k = key_fn(row)
        if k in seen:
            continue
        seen.add(k)
        out.append(row)

    with open(args.output, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=cols, extrasaction="ignore", lineterminator="\n")
        w.writeheader()
        w.writerows(out)

    if args.publish_if_changed and len(out) < len(rows):
        publish_path = Path(args.publish_if_changed)
        publish_path.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy(args.output, publish_path)

    print(f"filter_tech_reps: {len(rows)} -> {len(out)} ({args.mode})", file=sys.stderr)


if __name__ == "__main__":
    main()
