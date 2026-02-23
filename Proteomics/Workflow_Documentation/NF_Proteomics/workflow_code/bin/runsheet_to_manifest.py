#!/usr/bin/env python3

"""
Convert runsheet CSV to FragPipe manifest TSV format.

Manifest format (4 columns, tab-separated, no header):
1. Path (input file, e.g. SampleName_assay_suffix.mzML)
2. Experiment (group ID; FragPipe appends _Bioreplicate to form sample names)
3. Bioreplicate (per-condition: 1,2,3 within each condition; matches experiment_annotation replicate)
4. Data type (DDA, DIA, etc.)

For LFQ: Experiment = Sample Name; bioreplicate = per-condition (inferred) or from optional "Bioreplicate" runsheet column.
Optional "Bioreplicate" column: if present, use it; else infer per-condition (1,2,3...).
"""

import argparse
import csv
import re
import sys


def _make_names_safe(s: str) -> str:
    """Mimic R's make.names(): valid R symbol from arbitrary string."""
    if not s:
        return ""
    out = re.sub(r"[^a-zA-Z0-9._]", ".", s.strip())
    out = re.sub(r"\.+", ".", out).strip(".")
    if not out:
        return ""
    if out[0].isdigit():
        out = "X" + out
    return out


def _condition_from_factors(row: dict, factor_columns: list) -> str:
    """Join all Factor Values with ' & ', then sanitize for R."""
    values = [row.get(c, "").strip() for c in factor_columns if row.get(c, "").strip()]
    raw = " & ".join(values)
    return _make_names_safe(raw) if raw else ""


def main():
    parser = argparse.ArgumentParser(
        description="Convert runsheet CSV to FragPipe manifest TSV"
    )
    parser.add_argument("--runsheet", required=True, help="Path to runsheet CSV file")
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
        "--mode",
        default="LFQ",
        choices=["LFQ", "TMT"],
        help="Runsheet mode (LFQ or TMT)",
    )
    args = parser.parse_args()

    # Read runsheet
    with open(args.runsheet, "r") as f:
        reader = csv.DictReader(f)
        rows = list(reader)

    if not rows:
        sys.exit("Error: Runsheet is empty")

    # Find Factor Value columns
    factor_columns = [
        col for col in reader.fieldnames if col.startswith("Factor Value[")
    ]

    # Bioreplicate: optional runsheet column, else infer per-condition (LFQ) or from Source Name (TMT)
    has_bioreplicate_col = "Bioreplicate" in (reader.fieldnames or [])
    cond_to_biorep = {}
    sample_to_biorep = {}
    for row in rows:
        sample_name = row.get("Sample Name", "").strip()
        if not sample_name:
            continue
        if has_bioreplicate_col and row.get("Bioreplicate", "").strip():
            sample_to_biorep[sample_name] = str(row.get("Bioreplicate", "").strip())
        elif args.mode == "LFQ":
            condition = _condition_from_factors(row, factor_columns) or sample_name
            if condition not in cond_to_biorep:
                cond_to_biorep[condition] = 0
            cond_to_biorep[condition] += 1
            sample_to_biorep[sample_name] = str(cond_to_biorep[condition])
        else:
            # TMT: Source Name order (legacy)
            source_name = row.get("Source Name", "").strip()
            if source_name not in cond_to_biorep:
                cond_to_biorep[source_name] = len(cond_to_biorep) + 1
            sample_to_biorep[sample_name] = str(cond_to_biorep[source_name])

    # Write manifest
    with open(args.output, "w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")

        for row in rows:
            input_file_path = row.get("data_file", "").strip()
            if not input_file_path:
                continue

            sample_name = row.get("Sample Name", "").strip()
            if not sample_name:
                sys.exit(f"Error: Missing 'Sample Name' column in runsheet for row: {row}")

            input_file = f"{sample_name}.mzML"

            if args.mode == "LFQ":
                experiment = sample_name
            else:
                factor_values = []
                for factor_col in factor_columns:
                    value = row.get(factor_col, "").strip()
                    if value:
                        factor_values.append(value)
                experiment = " & ".join(factor_values) if factor_values else ""

            bioreplicate = sample_to_biorep.get(sample_name, "1")
            data_type = row.get("data_type", "").strip()

            writer.writerow([input_file, experiment, bioreplicate, data_type])

    print(f"Manifest written to {args.output}")


if __name__ == "__main__":
    main()

