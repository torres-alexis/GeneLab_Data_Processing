#!/usr/bin/env python
"""Opt-in selection of technical runs by distinct identified peptide sequence.

Search all technical replicates first (--tech_rep all), then supply filtered
FragPipe psm.tsv files with --psm. Combined peptide tables are supported only
when their per-sample columns map unambiguously to individual sheet rows.
This helper does not change the workflow's --tech_rep first default.
"""

from __future__ import annotations

import argparse
import csv
import io
import math
import re
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from filter_tech_reps import (  # noqa: E402
    _dedup_key_lfq,
    _dedup_key_tmt,
    _guard_lfq_tech_rep_needs_fraction,
    _row_collapse_tech_reps,
)
from tech_rep_utils import parse_has_tech_reps  # noqa: E402
from decoy_contam import row_is_junk  # noqa: E402


LOG_FIELDS = [
    "group_key", "action", "Sample Name", "data_file", "peptide_count", "reason",
    "row_number", "kept_row_number", "quant_sample", "selection_metric",
    "evidence_type", "evidence_files",
]
SPECTRUM_RE = re.compile(r"^(.+)\.\d+\.\d+\.\d+$")


def _table(path: Path):
    text = path.read_text(encoding="utf-8-sig")
    if not text.strip():
        raise ValueError(f"empty table: {path}")
    dialect = csv.Sniffer().sniff(text.splitlines()[0], delimiters=",\t")
    return csv.DictReader(io.StringIO(text), dialect=dialect)


def _check_width(row, path, line):
    if None in row or any(value is None for value in row.values()):
        raise ValueError(f"{path}:{line}: malformed table row (wrong number of fields)")


def _positive(value, column: str, line: int) -> bool:
    # Blank/NA intensities represent missing observations, not detections.
    if value is not None and value.strip().lower() in ("", "na", "nan"):
        return False
    try:
        number = float(value)
    except (TypeError, ValueError):
        raise ValueError(f"invalid measurement in {column!r} at line {line}: {value!r}")
    if not math.isfinite(number) or number < 0:
        raise ValueError(f"invalid measurement in {column!r} at line {line}: {value!r}")
    return number > 0


def per_sample_peptide_counts(quant_path: Path, allow_intensity=False) -> dict[str, int]:
    reader = _table(quant_path)
    fields = reader.fieldnames or []
    if "Peptide Sequence" not in fields:
        raise ValueError("--quant requires combined_peptide.tsv with Peptide Sequence; "
                         "protein tables do not measure peptide counts")
    if len(fields) != len(set(fields)):
        raise ValueError("duplicate combined-peptide column names are ambiguous")
    cols = [c for c in fields if c.endswith(" Spectral Count")
            and not c.endswith((" Unique Spectral Count", " Total Spectral Count"))
            and c not in ("Combined Spectral Count", "Total Spectral Count", "Unique Spectral Count")]
    suffix = " Spectral Count"
    if not cols and allow_intensity:
        cols = [c for c in fields if c.endswith(" Intensity")
                and not c.endswith(" MaxLFQ Intensity")
                and c not in ("Combined Intensity", "MaxLFQ Intensity", "Total Intensity")]
        suffix = " Intensity"
    if not cols:
        raise ValueError("no per-sample Spectral Count columns; use filtered --psm files. "
                         "Intensity-only ranking requires --allow-intensity-fallback "
                         "and can include MBR detections")
    peptides = {c: set() for c in cols}
    for line, row in enumerate(reader, start=2):
        _check_width(row, quant_path, line)
        if row_is_junk(row, columns=[c for c in ("Protein", "Protein ID") if c in fields]):
            continue
        sequence = (row.get("Peptide Sequence") or "").strip()
        if not sequence:
            raise ValueError(f"missing Peptide Sequence at line {line}")
        for col in cols:
            if _positive(row.get(col), col, line):
                peptides[col].add(sequence)
    return {col[:-len(suffix)]: len(sequences) for col, sequences in peptides.items()}


def per_run_psm_counts(paths: list[Path]) -> dict[str, int]:
    """Count distinct unmodified sequences, excluding explicit decoys/contaminants.

    Spectrum is <run>.<scan>.<scan>.<charge>. Spectrum File is deliberately not
    used: Philosopher may write interact-<run>.pep.xml there, not an mzML name.
    """
    peptides = defaultdict(set)
    for path in paths:
        reader = _table(path)
        required = ("Spectrum", "Peptide", "Is Decoy", "Is Contaminant")
        if any((reader.fieldnames or []).count(col) != 1 for col in required):
            raise ValueError(f"{path}: expected one of each filtered PSM column: {', '.join(required)}")
        for line, row in enumerate(reader, start=2):
            _check_width(row, path, line)
            match = SPECTRUM_RE.fullmatch((row.get("Spectrum") or "").strip())
            sequence = (row.get("Peptide") or "").strip()
            if not match or not sequence:
                raise ValueError(f"{path}:{line}: missing peptide or invalid Spectrum identifier")
            run = Path(match.group(1)).name
            peptides[run]  # Keep observed runs with zero target peptides as zero.
            flags = [parse_has_tech_reps(row.get(col)) for col in required[2:]]
            if any(flag is None for flag in flags):
                raise ValueError(f"{path}:{line}: invalid Is Decoy / Is Contaminant flag")
            if not any(flags):
                peptides[run].add(sequence)
    return {run: len(sequences) for run, sequences in peptides.items()}


def _row_aliases(row: dict) -> list[str]:
    aliases = [(row.get(key) or "").strip() for key in ("Sample Name", "run")]
    data_file = (row.get("data_file") or "").strip()
    if data_file:
        stem = Path(data_file).name
        if Path(stem).suffix.lower() in (".mzml", ".raw"):
            stem = Path(stem).stem
        aliases.append(stem)
    return list(dict.fromkeys(alias for alias in aliases if alias))


def select_rows(rows, cols, counts, mode, metric):
    if mode not in ("lfq", "tmt"):
        raise ValueError(f"unsupported mode: {mode}")
    aliases = defaultdict(set)
    for i, row in enumerate(rows):
        value = (row.get("Has Tech Reps") or "").strip()
        if value and parse_has_tech_reps(value) is None:
            raise ValueError(f"row {i + 1}: invalid Has Tech Reps value {value!r}")
        for alias in _row_aliases(row):
            aliases[alias].add(i)
    if mode == "lfq":
        _guard_lfq_tech_rep_needs_fraction(rows, cols)
    key_fn = _dedup_key_tmt if mode == "tmt" else _dedup_key_lfq
    groups = defaultdict(list)
    for i, row in enumerate(rows):
        if not _row_collapse_tech_reps(row, cols):
            # The selection helper must never collapse unflagged biological rows.
            key = ("unflagged", i + 1)
        else:
            if mode == "tmt" and any(not (row.get(c) or "").strip()
                                     for c in ("plex", "TechRepMixture", "fraction")):
                raise ValueError(f"row {i + 1}: TMT technical runs require plex, TechRepMixture and fraction")
            key = key_fn(row, cols)
        groups[key].append(i)

    kept, audit = set(), {}
    for key, members in groups.items():
        scores, matched = {}, {}
        reason = "not_technical_replicate" if key[0] == "unflagged" else "only_row"
        winner = members[0]
        if len(members) > 1:
            for i in members:
                matches = [a for a in _row_aliases(rows[i]) if a in counts and aliases[a] == {i}]
                if len(matches) != 1:
                    raise ValueError(
                        f"row {i + 1}: expected one unambiguous run-level quant match, found {matches}. "
                        "Missing or shared experiment/bioreplicate counts cannot rank runs; "
                        "supply filtered --psm files for all candidate runs.")
                matched[i] = matches[0]
                scores[i] = counts[matches[0]]
            winner = max(members, key=lambda i: scores[i])  # Stable sheet-order tie break.
            best = scores[winner]
            ties = sum(score == best for score in scores.values())
            reason = "all_zero_keep_first" if best == 0 else "tie_keep_first" if ties > 1 else "most_peptides"
        kept.add(winner)
        for i in members:
            row = rows[i]
            audit[i] = {
                "group_key": "|".join(str(x) for x in key),
                "action": "kept" if i == winner else "dropped",
                "Sample Name": (row.get("Sample Name") or row.get("run") or "").strip(),
                "data_file": (row.get("data_file") or "").strip(),
                "peptide_count": scores.get(i, ""),
                "reason": reason,
                "row_number": i + 1,
                "kept_row_number": winner + 1,
                "quant_sample": matched.get(i, ""),
                "selection_metric": metric if scores else "not_ranked",
            }
    return [row for i, row in enumerate(rows) if i in kept], [audit[i] for i in range(len(rows))]


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--mode", choices=("lfq", "tmt"), required=True)
    ap.add_argument("--input", required=True, help="Original CSV runsheet / data sheet (all technical runs)")
    source = ap.add_mutually_exclusive_group(required=True)
    source.add_argument("--quant", help="combined_peptide.tsv with individually resolvable run columns")
    source.add_argument("--psm", nargs="+", help="Filtered FragPipe psm.tsv files from all candidate runs")
    ap.add_argument("--allow-intensity-fallback", action="store_true",
                    help="Allow --quant intensity-only detections, which may include MBR; not identified-only ranking")
    ap.add_argument("--output", required=True, help="Selected CSV sheet")
    ap.add_argument("--log", default="tech_reps_selected.tsv", help="Selection audit TSV")
    args = ap.parse_args()
    try:
        inputs = [Path(args.input)] + ([Path(args.quant)] if args.quant else [Path(p) for p in args.psm])
        outputs = [Path(args.output).resolve(), Path(args.log).resolve()]
        if outputs[0] == outputs[1] or set(outputs) & {p.resolve() for p in inputs}:
            raise ValueError("output and log must be distinct from each other and all inputs")
        if args.psm and args.allow_intensity_fallback:
            raise ValueError("--allow-intensity-fallback applies only to --quant")
        if args.psm:
            counts = per_run_psm_counts([Path(p) for p in args.psm])
            metric = "distinct_identified_peptide_sequences"
        else:
            counts = per_sample_peptide_counts(Path(args.quant), args.allow_intensity_fallback)
            fields = _table(Path(args.quant)).fieldnames or []
            metric = "distinct_identified_peptide_sequences" if any(
                c.endswith(" Spectral Count") and c[:-len(" Spectral Count")] in counts for c in fields
            ) else "distinct_intensity_detected_peptide_sequences"
        with open(args.input, newline="", encoding="utf-8-sig") as handle:
            reader = csv.DictReader(handle)
            cols = reader.fieldnames or []
            if not cols or len(cols) != len(set(cols)):
                raise ValueError("input sheet must have nonempty, unique CSV column names")
            rows = list(reader)
        if any(None in row or any(value is None for value in row.values()) for row in rows):
            raise ValueError("input sheet has a malformed row (wrong number of CSV fields)")
        out_rows, log_rows = select_rows(rows, cols, counts, args.mode, metric)
        for entry in log_rows:
            entry["evidence_type"] = "filtered_psm" if args.psm else "combined_peptide"
            entry["evidence_files"] = ";".join(str(p.resolve()) for p in inputs[1:])
        # All scientific/mapping validation completes before any output is opened.
        for path in outputs:
            path.parent.mkdir(parents=True, exist_ok=True)
        with open(args.output, "w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=cols, lineterminator="\n")
            writer.writeheader()
            writer.writerows(out_rows)
        with open(args.log, "w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=LOG_FIELDS, delimiter="\t", lineterminator="\n")
            writer.writeheader()
            writer.writerows(log_rows)
    except (OSError, ValueError, csv.Error) as exc:
        ap.exit(2, f"Error: {exc}\n")
    print(f"select_tech_reps_by_peptides: {len(rows)} -> {len(out_rows)} ({args.mode})", file=sys.stderr)


if __name__ == "__main__":
    main()
