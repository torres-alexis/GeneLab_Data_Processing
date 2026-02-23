#!/usr/bin/env python3

"""
Update FragPipe experiment_annotation.tsv with condition, replicate, and basenames from runsheet.

FragPipe generates experiment_annotation from the manifest; for LFQ with
Experiment=Sample Name, condition is parsed from first token (wrong for multi-factor designs).
This script overwrites:
- file: basename only (strip path)
- sample_name: basename for display (e.g. Std_Mars_L1); sample unchanged (must match FragPipe quant columns)
- condition: all Factor Values joined with " & ", sanitized via R make.names (FragPipe-Analyst compatible)
- replicate: per-condition (1,2,3 within each condition) for FragPipe-Analyst design matrix

Input: FragPipe experiment_annotation.tsv, runsheet CSV
Output: Updated experiment_annotation.tsv
"""

import argparse
import csv
import os
import re
import sys


def _make_names_safe(s: str) -> str:
    """Mimic R's make.names(): valid R symbol from arbitrary string."""
    if not s:
        return ""
    # Replace invalid chars with . (R make.names behavior)
    out = re.sub(r"[^a-zA-Z0-9._]", ".", s.strip())
    out = re.sub(r"\.+", ".", out).strip(".")
    if not out:
        return ""
    if out[0].isdigit():
        out = "X" + out
    return out


def _condition_from_factors(row: dict, factor_columns: list) -> str:
    """Join all Factor Values with ' & ', then sanitize for R (FragPipe-Analyst)."""
    values = [row.get(c, "").strip() for c in factor_columns if row.get(c, "").strip()]
    raw = " & ".join(values)
    return _make_names_safe(raw) if raw else ""


def main():
    parser = argparse.ArgumentParser(
        description="Update experiment_annotation with condition/replicate from runsheet"
    )
    parser.add_argument("--experiment_annotation", required=True, help="FragPipe experiment_annotation.tsv")
    parser.add_argument("--runsheet", required=True, help="Runsheet CSV")
    parser.add_argument("--output", default="experiment_annotation.tsv", help="Output path")
    args = parser.parse_args()

    # Build Sample Name -> (condition, source_name) from runsheet
    factor_columns = []
    with open(args.runsheet, "r") as f:
        reader = csv.DictReader(f)
        factor_columns = [c for c in reader.fieldnames if c.startswith("Factor Value[")]
        rows = list(reader)

    sample_to_meta = {}
    for row in rows:
        sample_name = row.get("Sample Name", "").strip()
        if not sample_name:
            continue
        condition = _condition_from_factors(row, factor_columns) or sample_name
        source_name = row.get("Source Name", "").strip()
        sample_to_meta[sample_name] = {"condition": condition, "source_name": source_name}

    # Per-condition replicate (1,2,3 within each condition) for FragPipe-Analyst design matrix
    cond_to_repl = {}
    for row in rows:
        sample_name = row.get("Sample Name", "").strip()
        if sample_name not in sample_to_meta:
            continue
        cond = sample_to_meta[sample_name]["condition"]
        if cond not in cond_to_repl:
            cond_to_repl[cond] = 0
        cond_to_repl[cond] += 1
        sample_to_meta[sample_name]["replicate"] = cond_to_repl[cond]

    # Read FragPipe experiment_annotation
    with open(args.experiment_annotation, "r") as f:
        reader = csv.DictReader(f, delimiter="\t")
        rows_anno = list(reader)
        fieldnames = list(reader.fieldnames)
    if "sample_name" not in fieldnames:
        fieldnames.append("sample_name")

    # Update file (basename), condition, replicate; sample_name = basename for display
    for row in rows_anno:
        sample = row.get("sample", "").strip()
        # sample = Experiment_Bioreplicate (e.g. Std_Mars_L1_1); Sample Name = Experiment (basename)
        sample_name = re.sub(r"_\d+$", "", sample) if sample else ""
        if "file" in row and row.get("file"):
            row["file"] = os.path.basename(row["file"])
        if sample_name in sample_to_meta:
            meta = sample_to_meta[sample_name]
            row["condition"] = meta["condition"]
            row["replicate"] = str(meta["replicate"])
        row["sample_name"] = sample_name or sample  # basename for display

    with open(args.output, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows_anno)

    print(f"Updated experiment_annotation written to {args.output}")


if __name__ == "__main__":
    main()
