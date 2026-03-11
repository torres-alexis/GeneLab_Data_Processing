#!/usr/bin/env python

"""
Convert runsheet (LFQ) or data_sheet + sample_sheet (TMT) to FragPipe manifest.tsv and experiment_annotation.tsv.

  --runsheet: LFQ input (one row per mzML). --data_sheet: TMT input. Mode: data_sheet→TMT, runsheet→LFQ.
  --sample_sheet: TMT only; required for experiment_annotation.

Outputs: manifest.tsv (Path|Experiment|Bioreplicate|Data type, no header), experiment_annotation.tsv.
LFQ: Experiment from Factor Value columns or "1"; Bioreplicate from column or sequential per condition.
TMT: Experiment=plex; Bioreplicate=TechRepMixture or "1"; plex column must match FragPipe folder (plex_Bioreplicate).
"""

import argparse
import csv
import os
import re
import sys


def _make_names_safe(s: str) -> str:
    """R-safe symbol: [^a-zA-Z0-9._] -> ., leading digit -> X prefix."""
    if not s:
        return ""
    out = re.sub(r"[^a-zA-Z0-9._]", ".", s.strip())
    out = re.sub(r"\.+", ".", out).strip(".")
    if not out:
        return ""
    if out[0].isdigit():
        out = "X" + out
    return out


def _condition_label_from_factors(row: dict, factor_columns: list) -> str:
    """Join Factor Values with ' & '."""
    values = [row.get(c, "").strip() for c in factor_columns if row.get(c, "").strip()]
    return " & ".join(values) if values else ""


def _sanitize_for_fragpipe(s: str) -> str:
    """FragPipe Experiment: [^A-Za-z0-9_] -> _."""
    if not s:
        return ""
    return re.sub(r"[^A-Za-z0-9_]", "_", s.strip())


def _condition_from_factors(row: dict, factor_columns: list) -> str:
    """FragPipe sanitize + R-safe."""
    raw = _condition_label_from_factors(row, factor_columns)
    if not raw:
        return ""
    return _make_names_safe(_sanitize_for_fragpipe(raw))


def main():
    parser = argparse.ArgumentParser(
        description="Convert runsheet CSV to FragPipe manifest TSV"
    )
    parser.add_argument("--runsheet", default="", help="LFQ: path to runsheet CSV (one row per mzML)")
    parser.add_argument("--data_sheet", default="", help="TMT: path to data sheet CSV (one row per mzML). Use --runsheet or --data_sheet.")
    parser.add_argument(
        "--assay_suffix",
        default="",
        help="Suffix to append to sample names in output filenames (e.g., _GLProteomics)",
    )
    parser.add_argument(
        "--output",
        default="manifest.tsv",
        help="Output manifest TSV file path",
    )
    parser.add_argument(
        "--sample_sheet",
        default="",
        help="TMT only: path to sample_sheet CSV (channel-centric). Required for TMT experiment_annotation.",
    )
    args = parser.parse_args()

    input_file = args.data_sheet or args.runsheet
    if not input_file:
        sys.exit("Error: Provide --runsheet (LFQ) or --data_sheet (TMT)")
    mode = "TMT" if args.data_sheet else "LFQ"

    # Read input
    with open(input_file, "r") as f:
        reader = csv.DictReader(f)
        rows = list(reader)

    if not rows:
        sys.exit("Error: Input is empty")

    factor_columns = [
        col for col in reader.fieldnames if col.startswith("Factor Value[")
    ]

    # Bioreplicate
    has_bioreplicate_col = "Bioreplicate" in (reader.fieldnames or [])
    cond_to_biorep = {}
    sample_to_biorep = {}
    sample_to_experiment = {}

    def _get_sample_id(row):
        return row.get("run", "").strip() or row.get("Sample Name", "").strip()

    for row in rows:
        sample_name = _get_sample_id(row)
        if not sample_name:
            continue
        if mode == "LFQ":
            cond_label = _condition_label_from_factors(row, factor_columns)
            if cond_label:
                sample_to_experiment[sample_name] = _sanitize_for_fragpipe(cond_label)
            else:
                sample_to_experiment[sample_name] = "1"
        if mode == "LFQ" and has_bioreplicate_col and row.get("Bioreplicate", "").strip():
            sample_to_biorep[sample_name] = str(row.get("Bioreplicate", "").strip())
        elif mode == "LFQ":
            condition = _condition_from_factors(row, factor_columns)
            if condition:
                if condition not in cond_to_biorep:
                    cond_to_biorep[condition] = 0
                cond_to_biorep[condition] += 1
                offset = sum(cond_to_biorep.get(c, 0) for c in cond_to_biorep if c != condition)
                sample_to_biorep[sample_name] = str(offset + cond_to_biorep[condition])
            else:
                sample_to_biorep[sample_name] = "1"
        elif mode == "TMT":
            plex = row.get("plex", "").strip()
            if not plex:
                sys.exit(f"Error: TMT mode requires 'plex' column in data sheet. Missing for Sample Name: {sample_name}")
            biorep = row.get("TechRepMixture", "").strip()
            sample_to_biorep[sample_name] = biorep if biorep else "1"

    with open(args.output, "w", newline="\n") as f:
        writer = csv.writer(f, delimiter="\t")

        for row in rows:
            input_file_path = row.get("data_file", "").strip()
            if not input_file_path:
                continue

            sample_name = _get_sample_id(row)
            if not sample_name:
                sys.exit(f"Error: Missing 'run' or 'Sample Name' column in input for row: {row}")

            input_file = f"{sample_name}.mzML"
            experiment = sample_to_experiment.get(sample_name, "1") if mode == "LFQ" else row.get("plex", "").strip() or ""
            data_type = "DDA"
            bioreplicate = sample_to_biorep.get(sample_name, "1")
            writer.writerow([input_file, experiment or "1", bioreplicate or "1", data_type])

    print(f"Manifest written to {args.output}")

    # Write experiment_annotation.tsv 
    base = os.path.basename(args.output)
    exp_base = base.replace("manifest", "experiment_annotation", 1) if "manifest" in base else "experiment_annotation" + os.path.splitext(base)[1]
    exp_anno_path = os.path.join(os.path.dirname(args.output), exp_base)

    if mode == "LFQ":
        exp_fieldnames = ["file", "sample", "sample_name", "condition_label", "condition", "replicate"]
        exp_rows = []
        for row in rows:
            input_file_path = row.get("data_file", "").strip()
            if not input_file_path:
                continue
            sample_id = _get_sample_id(row)
            experiment = sample_to_experiment.get(sample_id, "1")
            bioreplicate = sample_to_biorep.get(sample_id, "1")
            sample = f"{experiment}_{bioreplicate}"
            cond_label = _condition_label_from_factors(row, factor_columns) or "Experiment"
            cond_safe = _condition_from_factors(row, factor_columns) or "Experiment"
            exp_rows.append({
                "file": f"{sample_id}.mzML",
                "sample": sample,
                "sample_name": sample,
                "condition_label": cond_label,
                "condition": cond_safe,
                "replicate": bioreplicate,
            })
        with open(exp_anno_path, "w", newline="\n") as f:
            writer = csv.DictWriter(f, fieldnames=exp_fieldnames, delimiter="\t")
            writer.writeheader()
            writer.writerows(exp_rows)
        print(f"Experiment annotation written to {exp_anno_path}")

    # TMT: write experiment_annotation from sample_sheet + data_sheet (requires --sample_sheet)
    # plex column must match FragPipe folder = Experiment_Bioreplicate (e.g. TMTa_1, TMTa_2 when multiple runs)
    # Per https://fragpipe.nesvilab.org/docs/tutorial_fragpipe.html: one annotation.txt per plex folder
    if mode == "TMT" and args.sample_sheet:
        with open(args.sample_sheet, "r") as f:
            sample_reader = csv.DictReader(f)
            sample_rows = list(sample_reader)
            sample_fieldnames = sample_reader.fieldnames or []
        if not sample_rows:
            sys.exit("Error: Sample sheet is empty")
        factor_cols = [c for c in sample_fieldnames if c.startswith("Factor Value[")]

        folders_seen = set()
        folder_to_plex_base = {}  # folder -> plex base for sample_sheet lookup
        for row in rows:
            plex = row.get("plex", "").strip()
            if not plex:
                continue
            biorep = sample_to_biorep.get(_get_sample_id(row), "1")
            folder = f"{plex}_{biorep}"
            if folder not in folders_seen:
                folders_seen.add(folder)
                folder_to_plex_base[folder] = plex

        plex_base_to_channels = {}
        for row in sample_rows:
            plex = row.get("plex", "").strip()
            if not plex:
                continue
            if plex not in plex_base_to_channels:
                plex_base_to_channels[plex] = []
            plex_base_to_channels[plex].append(row)

        exp_fieldnames = ["plex", "channel", "sample", "sample_name", "condition_label", "condition", "replicate"]
        exp_rows = []
        for folder, plex_base in folder_to_plex_base.items():
            channel_rows = plex_base_to_channels.get(plex_base, [])
            for row in channel_rows:
                sample_name = row.get("Sample Name", "").strip()
                if not sample_name:
                    continue
                channel = row.get("channel", "").strip()
                replicate = row.get("Bioreplicate", "").strip() or "1"
                cond_label = _condition_label_from_factors(row, factor_cols)
                cond_safe = _condition_from_factors(row, factor_cols) if cond_label else _sanitize_for_fragpipe(sample_name)
                sample_fp = _sanitize_for_fragpipe(sample_name) or sample_name
                exp_rows.append({
                    "plex": folder,
                    "channel": channel,
                    "sample": sample_fp,
                    "sample_name": sample_fp,
                    "condition_label": cond_label or sample_name,
                    "condition": cond_safe or _make_names_safe(sample_name),
                    "replicate": replicate,
                })
        with open(exp_anno_path, "w", newline="\n") as f:
            writer = csv.DictWriter(f, fieldnames=exp_fieldnames, delimiter="\t")
            writer.writeheader()
            writer.writerows(exp_rows)
        print(f"Experiment annotation written to {exp_anno_path}")


if __name__ == "__main__":
    main()
