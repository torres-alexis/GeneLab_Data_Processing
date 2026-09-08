#!/usr/bin/env python
"""Pick the tech-rep row with the most identified peptides (AWG 8/12).

Two-pass SOP:
  1. Run FragPipe with --tech_rep all (search all TRs).
  2. Run this on combined_peptide.tsv (or combined_protein.tsv) + the original sheet.
  3. Rerun FragPipe on the written sheet (--tech_rep first).

Per-sample peptide count = rows with Spectral Count > 0 (else Intensity > 0)
for that sample's column. Columns map to Sample Name, run, or data_file stem
(not Original Sample Name — reruns share the parent original name).
"""

from __future__ import annotations

import argparse
import csv
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from filter_tech_reps import (  # noqa: E402
    _dedup_key_lfq,
    _dedup_key_tmt,
    _row_collapse_tech_reps,
)
from decoy_contam import _open_table  # noqa: E402


def _num(val) -> float:
    try:
        return float(val)
    except (TypeError, ValueError):
        return 0.0


def per_sample_peptide_counts(quant_path: Path) -> dict[str, int]:
    text, dialect = _open_table(quant_path)
    reader = csv.DictReader(text.splitlines(), dialect=dialect)
    fields = reader.fieldnames or []
    spec_cols = [
        c
        for c in fields
        if c.endswith(" Spectral Count")
        and "Unique" not in c
        and "Total" not in c
        and "Combined" not in c
    ]
    int_cols = [
        c
        for c in fields
        if c.endswith(" Intensity") and "MaxLFQ" not in c and "Combined" not in c
    ]
    cols = spec_cols or int_cols
    if not cols:
        sys.exit(f"Error: no per-sample Spectral Count / Intensity columns in {quant_path}")
    counts = {c: 0 for c in cols}
    for row in reader:
        for col in cols:
            if _num(row.get(col)) > 0:
                counts[col] += 1
    out = {}
    for col, n in counts.items():
        sample = col
        for suffix in (" Spectral Count", " Intensity"):
            if sample.endswith(suffix):
                sample = sample[: -len(suffix)]
                break
        out[sample] = n
    return out


def _row_aliases(row: dict) -> list[str]:
    aliases = []
    for key in ("Sample Name", "run"):
        val = (row.get(key) or "").strip()
        if val:
            aliases.append(val)
    data_file = (row.get("data_file") or "").strip()
    if data_file:
        stem = Path(data_file).name
        for ext in (".mzML", ".mzml", ".raw", ".RAW"):
            if stem.endswith(ext):
                stem = stem[: -len(ext)]
                break
        aliases.append(stem)
    return aliases


def _count_for_row(row: dict, counts: dict[str, int]) -> int | None:
    for alias in _row_aliases(row):
        if alias in counts:
            return counts[alias]
    return None


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--mode", choices=("lfq", "tmt"), required=True)
    ap.add_argument("--input", required=True, help="Original runsheet / data sheet (all tech reps)")
    ap.add_argument("--quant", required=True, help="combined_peptide.tsv or combined_protein.tsv")
    ap.add_argument("--output", required=True, help="Filtered sheet keeping max-peptide tech rep")
    ap.add_argument("--log", default="tech_reps_selected.tsv")
    args = ap.parse_args()

    counts = per_sample_peptide_counts(Path(args.quant))
    with open(args.input, newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        cols = reader.fieldnames or []
        rows = list(reader)

    key_fn = (
        (lambda r: _dedup_key_tmt(r, cols))
        if args.mode == "tmt"
        else (lambda r: _dedup_key_lfq(r, cols))
    )
    groups = defaultdict(list)
    for i, row in enumerate(rows):
        groups[key_fn(row)].append((i, row))

    keep_idx = set()
    log_rows = []
    for key, members in groups.items():
        collapsing = any(_row_collapse_tech_reps(r, cols) for _, r in members)
        scored = []
        for i, row in members:
            n = _count_for_row(row, counts)
            scored.append((n if n is not None else -1, i, row))
        scored.sort(key=lambda t: (t[0], -t[1]), reverse=True)
        winner_n, winner_i, winner = scored[0]
        keep_idx.add(winner_i)
        reason = "most_peptides" if collapsing and len(members) > 1 else "only_row"
        if winner_n < 0:
            reason = "no_quant_match_keep_first"
        for n, i, row in scored:
            log_rows.append(
                {
                    "group_key": "|".join(str(x) for x in key),
                    "action": "kept" if i == winner_i else "dropped",
                    "Sample Name": (row.get("Sample Name") or row.get("run") or "").strip(),
                    "data_file": (row.get("data_file") or "").strip(),
                    "peptide_count": "" if n < 0 else str(n),
                    "reason": reason,
                }
            )

    out_rows = [row for i, row in enumerate(rows) if i in keep_idx]
    with open(args.output, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=cols, extrasaction="ignore", lineterminator="\n")
        writer.writeheader()
        writer.writerows(out_rows)

    log_path = Path(args.log)
    log_path.parent.mkdir(parents=True, exist_ok=True)
    with log_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["group_key", "action", "Sample Name", "data_file", "peptide_count", "reason"],
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(log_rows)

    print(
        f"select_tech_reps_by_peptides: {len(rows)} -> {len(out_rows)} ({args.mode})",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
