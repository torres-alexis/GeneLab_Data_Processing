#!/usr/bin/env python
"""Keep first row per technical-replicate group (order preserved).

Has Tech Reps column optional:
  Missing column, blank, or FALSE: row is kept (unique per Sample Name + run
  and data_file basename [LFQ], or run [TMT]).
  TRUE: collapse with other TRUE rows sharing:
    LFQ: Source Name + Factor Value[*] + fraction (fraction optional; blank if absent)
    TMT: plex + TechRepMixture + fraction

LFQ: same Sample Name with multiple data_files and Has Tech Reps=TRUE requires
non-empty fraction on those rows (tech-rep collapse is per fraction). Without
fraction, that pattern errors. Fraction-only sheets leave Has Tech Reps unset.
"""

import argparse
import csv
import shutil
import sys
from collections import defaultdict
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


def _data_file_basename(row) -> str:
    return Path((row.get("data_file") or "").strip()).name


def _fraction_value(row) -> str:
    fr = row.get("fraction")
    return ("" if fr is None else str(fr)).strip()


def _unique_key_lfq(row):
    run = (row.get("run") or "").strip()
    if run:
        return ("unique", run)
    return ("unique", row.get("Sample Name", ""), _data_file_basename(row))


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
    return ("source_tech", source, _factor_tuple(row, fieldnames), _fraction_value(row))


def _dedup_key_tmt(row, fieldnames):
    if not _row_collapse_tech_reps(row, fieldnames):
        return _unique_key_tmt(row)

    plex = (row.get("plex") or "").strip()
    tr = str(row.get("TechRepMixture", "") or "").strip() or "1"
    return ("tmt_tech", plex, tr, _fraction_value(row))


def _guard_lfq_tech_rep_needs_fraction(rows, fieldnames) -> None:
    """Require fraction when Has Tech Reps=TRUE would otherwise merge distinct fraction files."""
    by_sample = defaultdict(list)
    for row in rows:
        if not _row_collapse_tech_reps(row, fieldnames):
            continue
        sid = (row.get("Sample Name") or "").strip()
        base = _data_file_basename(row)
        if sid and base:
            by_sample[sid].append(row)

    for sid, group in by_sample.items():
        bases = {_data_file_basename(r) for r in group}
        if len(bases) <= 1:
            continue
        missing = [r for r in group if not _fraction_value(r)]
        if missing:
            sys.exit(
                "Error: Has Tech Reps=TRUE with multiple data_files for "
                f"Sample Name={sid!r} requires a non-empty fraction on each of those rows "
                f"(data_files: {', '.join(sorted(bases))}). "
                "Set fraction per LC slice so tech-rep collapse is Source Name + factors + fraction, "
                "or leave Has Tech Reps unset/FALSE on fraction-only rows."
            )


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
    ap.add_argument(
        "--drop-log",
        default="tech_reps_dropped.tsv",
        help="TSV of kept/dropped tech-rep rows (empty groups still get a header).",
    )
    args = ap.parse_args()

    with open(args.input, newline="", encoding="utf-8") as f:
        rd = csv.DictReader(f)
        cols = rd.fieldnames or []
        rows = list(rd)

    if args.mode == "lfq":
        _guard_lfq_tech_rep_needs_fraction(rows, cols)

    key_fn = (
        (lambda r: _dedup_key_tmt(r, cols))
        if args.mode == "tmt"
        else (lambda r: _dedup_key_lfq(r, cols))
    )
    seen = set()
    out = []
    log_rows = []
    key_kept = {}
    for row in rows:
        k = key_fn(row)
        sample = (row.get("Sample Name") or row.get("run") or "").strip()
        data_file = (row.get("data_file") or "").strip()
        collapsing = _row_collapse_tech_reps(row, cols)
        if k in seen:
            log_rows.append(
                {
                    "group_key": "|".join(str(x) for x in k),
                    "action": "dropped",
                    "Sample Name": sample,
                    "data_file": data_file,
                    "kept_sample": key_kept.get(k, ""),
                    "reason": "keep-first" if collapsing else "duplicate unique key",
                }
            )
            continue
        seen.add(k)
        out.append(row)
        key_kept[k] = sample
        if collapsing:
            log_rows.append(
                {
                    "group_key": "|".join(str(x) for x in k),
                    "action": "kept",
                    "Sample Name": sample,
                    "data_file": data_file,
                    "kept_sample": sample,
                    "reason": "keep-first",
                }
            )

    with open(args.output, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=cols, extrasaction="ignore", lineterminator="\n")
        w.writeheader()
        w.writerows(out)

    if args.publish_if_changed and len(out) < len(rows):
        publish_path = Path(args.publish_if_changed)
        publish_path.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy(args.output, publish_path)

    if args.drop_log:
        log_path = Path(args.drop_log)
        log_path.parent.mkdir(parents=True, exist_ok=True)
        log_fields = [
            "group_key",
            "action",
            "Sample Name",
            "data_file",
            "kept_sample",
            "reason",
        ]
        with log_path.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=log_fields, delimiter="\t", lineterminator="\n")
            writer.writeheader()
            writer.writerows(log_rows)

    print(f"filter_tech_reps: {len(rows)} -> {len(out)} ({args.mode})", file=sys.stderr)


if __name__ == "__main__":
    main()
