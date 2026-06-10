#!/usr/bin/env python

"""
Convert runsheet (LFQ) or data_sheet + sample_sheet (TMT) to downstream metadata tables:
FragPipe manifest, FragPipeAnalystR experiment_annotation, and (TMT) MSstatsTMT annotation.

  --runsheet: LFQ input (one row per mzML). TMT: pass --data_sheet first, then --sample_sheet. Mode: data_sheet→TMT, runsheet→LFQ.
  --sample_sheet: TMT only; required for experiment_annotation (paired with --data_sheet).

Outputs (cwd): manifest[.suffix].tsv, experiment_annotation[.suffix].tsv, and (TMT) MSstatsTMT_annotation[.suffix].csv.
Manifest data_type column: from runsheet/data_sheet `data_type` per row.
LFQ: Experiment from Factor Value columns or "1"; Bioreplicate from column or sequential per condition.
  LFQ experiment_annotation: sample = `{Experiment}_{Bioreplicate}` (quant match); sample_name = runsheet 'Sample Name'
  TMT: Experiment=plex; manifest Bioreplicate=TechRepMixture from data sheet (required). fraction and TechRepMixture must be set on every data sheet row. plex column must match FragPipe folder (plex_Bioreplicate).
  MSstatsTMT annotation only: sample-sheet condition Pool is written as Norm (bridge channel); experiment_annotation keeps Pool for FragPipe/FPAR.
  If one logical plex spans multiple folders (TechRepMixture / fraction batches), sample/sample_name become <folder>_<Sample Name> (batch prefix). Single folder per plex → no prefix.
"""

import argparse
import csv
import re
import sys
from collections import Counter

DATA_TYPE_CHOICES = ("DDA", "DIA")


def _row_data_type(row: dict) -> str:
    val = (row.get("data_type") or "").strip()
    if not val:
        return "DDA"
    if val not in DATA_TYPE_CHOICES:
        sys.exit(
            f"Error: invalid data_type '{val}' (expected one of {', '.join(DATA_TYPE_CHOICES)})"
        )
    return val


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


def _condition_name_from_factors(row: dict, factor_columns: list) -> str:
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
    raw = _condition_name_from_factors(row, factor_columns)
    if not raw:
        return ""
    return _make_names_safe(_sanitize_for_fragpipe(raw))


def _msstats_condition(condition: str) -> str:
    """MSstatsTMT bridge/pool channels must use Condition Norm (not in experiment_annotation)."""
    if condition.strip().lower() == "pool":
        return "Norm"
    return condition


def _require_tmt_cell(row: dict, column: str) -> str:
    """Require a non-empty TMT data sheet column value."""
    run = (row.get("run") or "").strip()
    val = (row.get(column) or "").strip()
    if not val:
        sys.exit(f"Error: data_sheet row missing {column} for run {run or row}")
    return val


def _require_sample_bioreplicate(row: dict) -> str:
    """Require a non-empty sample sheet Bioreplicate value."""
    sample_name = (row.get("Sample Name") or "").strip()
    val = (row.get("Bioreplicate") or "").strip()
    if not val:
        sys.exit(f"Error: sample_sheet row missing Bioreplicate for Sample Name {sample_name or row}")
    return val


def _write_msstats_tmt_annotation(
    data_rows: list,
    sample_rows: list,
    sample_fieldnames: list,
    output_path: str,
) -> None:
    """MSstatsTMT annotation: Run, Fraction, TechRepMixture, Mixture, Channel, BioReplicate, Condition."""
    factor_cols = [c for c in sample_fieldnames if c.startswith("Factor Value[")]

    plex_to_channels = {}
    for row in sample_rows:
        plex = (row.get("plex") or "").strip()
        channel = (row.get("channel") or "").strip()
        if not plex or not channel:
            continue
        sample_name = (row.get("Sample Name") or "").strip()
        biorep = _require_sample_bioreplicate(row)
        cond = _condition_from_factors(row, factor_cols)
        if not cond:
            cond = "Empty" if not sample_name else _make_names_safe(_sanitize_for_fragpipe(sample_name))
        plex_to_channels.setdefault(plex, []).append(
            {"channel": channel, "bioreplicate": biorep, "condition": _msstats_condition(cond)}
        )

    if not plex_to_channels:
        sys.exit("Error: sample_sheet has no plex/channel rows")

    fieldnames = [
        "Run",
        "Fraction",
        "TechRepMixture",
        "Mixture",
        "Channel",
        "BioReplicate",
        "Condition",
    ]
    out_rows = []
    for row in data_rows:
        run = (row.get("run") or "").strip()
        if not run:
            sys.exit(f"Error: data_sheet row missing run: {row}")
        plex = (row.get("plex") or "").strip()
        if not plex:
            sys.exit(f"Error: data_sheet row missing plex for run {run}")
        fraction = _require_tmt_cell(row, "fraction")
        tech_rep = _require_tmt_cell(row, "TechRepMixture")
        channels = plex_to_channels.get(plex)
        if not channels:
            sys.exit(f"Error: no sample_sheet channels for plex {plex} (run {run})")
        for ch in channels:
            out_rows.append(
                {
                    "Run": run,
                    "Fraction": fraction,
                    "TechRepMixture": tech_rep,
                    "Mixture": plex,
                    "Channel": ch["channel"],
                    "BioReplicate": ch["bioreplicate"],
                    "Condition": ch["condition"],
                }
            )

    with open(output_path, "w", newline="\n") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(out_rows)

    print(f"MSstatsTMT annotation written to {output_path} ({len(out_rows)} rows)")


def main():
    parser = argparse.ArgumentParser(
        description="Convert runsheet CSV to FragPipe manifest TSV"
    )
    parser.add_argument("--runsheet", default="", help="LFQ: path to runsheet CSV (one row per mzML)")
    parser.add_argument("--data_sheet", default="", help="TMT: path to data sheet CSV (one row per mzML). Use --runsheet or --data_sheet.")
    parser.add_argument(
        "--sample_sheet",
        default="",
        help="TMT only: path to sample_sheet CSV (channel-centric). Required for TMT experiment_annotation.",
    )
    parser.add_argument(
        "--assay_suffix",
        default="",
        help="Optional stem suffix for outputs (e.g. _GLProteomics → manifest_GLProteomics.tsv, MSstatsTMT_annotation_GLProteomics.csv). Default: manifest.tsv, experiment_annotation.tsv.",
    )
    args = parser.parse_args()

    suffix = (args.assay_suffix or "").strip()
    manifest_stem = f"manifest{suffix}" if suffix else "manifest"
    manifest_path = f"{manifest_stem}.tsv"
    exp_anno_path = f"{manifest_stem.replace('manifest', 'experiment_annotation', 1)}.tsv"
    msstats_tmt_anno_path = f"MSstatsTMT_annotation{suffix}.csv" if suffix else "MSstatsTMT_annotation.csv"

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
        r = (row.get("run") or "").strip()
        if r:
            return r
        return (row.get("Sample Name") or "").strip()

    for row in rows:
        sample_name = _get_sample_id(row)
        if not sample_name:
            continue
        if mode == "LFQ":
            cond_name = _condition_name_from_factors(row, factor_columns)
            if cond_name:
                sample_to_experiment[sample_name] = _sanitize_for_fragpipe(cond_name)
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
            if not row.get("data_file", "").strip():
                continue
            plex = row.get("plex", "").strip()
            if not plex:
                sys.exit(f"Error: TMT mode requires 'plex' column in data sheet. Missing for Sample Name: {sample_name}")
            sample_to_biorep[sample_name] = _require_tmt_cell(row, "TechRepMixture")

    with open(manifest_path, "w", newline="\n") as f:
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
            data_type = _row_data_type(row)
            if mode == "TMT":
                bioreplicate = sample_to_biorep.get(sample_name)
                if not bioreplicate:
                    sys.exit(f"Error: missing TechRepMixture for run {sample_name}")
            else:
                bioreplicate = sample_to_biorep.get(sample_name, "1")
            writer.writerow([input_file, experiment or "1", bioreplicate or "1", data_type])

    print(f"Manifest written to {manifest_path}")

    if mode == "LFQ":
        exp_fieldnames = ["file", "sample", "sample_name", "condition", "condition_name", "replicate"]
        sample_to_display_name = {}
        sample_display_from_source = {}
        for row in rows:
            input_file_path = row.get("data_file", "").strip()
            if not input_file_path:
                continue
            sample_id = _get_sample_id(row)
            experiment = sample_to_experiment.get(sample_id, "1")
            bioreplicate = sample_to_biorep.get(sample_id, "1")
            sample = f"{experiment}_{bioreplicate}"
            source_name = (row.get("Source Name") or "").strip()
            display_name = source_name or (row.get("Sample Name") or "").strip() or sample_id
            if sample not in sample_to_display_name or (source_name and not sample_display_from_source.get(sample, False)):
                sample_to_display_name[sample] = display_name
                sample_display_from_source[sample] = bool(source_name)

        exp_rows = []
        for row in rows:
            input_file_path = row.get("data_file", "").strip()
            if not input_file_path:
                continue
            sample_id = _get_sample_id(row)
            experiment = sample_to_experiment.get(sample_id, "1")
            bioreplicate = sample_to_biorep.get(sample_id, "1")
            sample = f"{experiment}_{bioreplicate}"
            sample_name_out = sample_to_display_name.get(sample, sample_id)
            cond_name = _condition_name_from_factors(row, factor_columns) or "Experiment"
            cond_safe = _condition_from_factors(row, factor_columns) or "Experiment"
            exp_rows.append({
                "file": f"{sample_id}.mzML",
                "sample": sample,
                "sample_name": sample_name_out,
                "condition": cond_safe,
                "condition_name": cond_name,
                "replicate": bioreplicate,
            })
        with open(exp_anno_path, "w", newline="\n") as f:
            writer = csv.DictWriter(f, fieldnames=exp_fieldnames, delimiter="\t")
            writer.writeheader()
            writer.writerows(exp_rows)
        print(f"Experiment annotation written to {exp_anno_path}")

    # TMT: write experiment_annotation from data_sheet + sample_sheet (requires --sample_sheet)
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

        # One logical plex in multiple folders → duplicate sample strings across folders; prefix all channels with folder_ for that plex.
        plex_needs_folder_prefix = {pb for pb, n in Counter(folder_to_plex_base.values()).items() if n > 1}

        plex_base_to_channels = {}
        for row in sample_rows:
            plex = row.get("plex", "").strip()
            if not plex:
                continue
            if plex not in plex_base_to_channels:
                plex_base_to_channels[plex] = []
            plex_base_to_channels[plex].append(row)

        exp_fieldnames = ["plex", "channel", "sample", "sample_name", "condition", "condition_name", "replicate"]
        exp_rows = []
        for folder, plex_base in folder_to_plex_base.items():
            channel_rows = plex_base_to_channels.get(plex_base, [])
            for row in channel_rows:
                sample_name = row.get("Sample Name", "").strip()
                if not sample_name:
                    continue
                channel = row.get("channel", "").strip()
                replicate = _require_sample_bioreplicate(row)
                cond_name = _condition_name_from_factors(row, factor_cols)
                cond_safe = _condition_from_factors(row, factor_cols) if cond_name else _sanitize_for_fragpipe(sample_name)
                sample_base = _sanitize_for_fragpipe(sample_name) or sample_name
                if plex_base in plex_needs_folder_prefix:
                    sample_fp = _sanitize_for_fragpipe(f"{folder}_{sample_name}")
                else:
                    sample_fp = sample_base
                exp_rows.append({
                    "plex": folder,
                    "channel": channel,
                    "sample": sample_fp,
                    "sample_name": sample_fp,
                    "condition": cond_safe,
                    "condition_name": cond_name,
                    "replicate": replicate,
                })
        with open(exp_anno_path, "w", newline="\n") as f:
            writer = csv.DictWriter(f, fieldnames=exp_fieldnames, delimiter="\t")
            writer.writeheader()
            writer.writerows(exp_rows)
        print(f"Experiment annotation written to {exp_anno_path}")

        data_rows = [r for r in rows if (r.get("data_file") or "").strip()]
        _write_msstats_tmt_annotation(data_rows, sample_rows, sample_fieldnames, msstats_tmt_anno_path)


if __name__ == "__main__":
    main()
