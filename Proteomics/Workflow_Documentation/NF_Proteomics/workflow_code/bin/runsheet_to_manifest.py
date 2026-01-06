#!/usr/bin/env python3

"""
Convert runsheet CSV to FragPipe manifest TSV format.

Manifest format (4 columns):
1. input_file (path to mzML file)
2. joined_factor_values (Factor Value columns joined with __, spaces replaced with _)
3. bioreplicate
4. data_type
"""

import argparse
import csv
import sys
import os


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

    # Write manifest
    with open(args.output, "w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")

        for row in rows:
            # Column 1: input_file (basename only)
            input_file_path = row.get("input_file", "").strip()
            if not input_file_path:
                continue

            # Get basename (filename only, no path)
            input_file = os.path.basename(input_file_path)

            # Apply assay_suffix to filename if provided
            if args.assay_suffix:
                # Add suffix before .mzML extension
                if input_file.endswith(".mzML"):
                    input_file = input_file.replace(".mzML", f"{args.assay_suffix}.mzML")
                elif input_file.endswith(".mzml"):
                    input_file = input_file.replace(".mzml", f"{args.assay_suffix}.mzML")
                else:
                    input_file = f"{input_file}{args.assay_suffix}.mzML"

            # Column 2: Join Factor Value columns with __, replace spaces with _
            factor_values = []
            for factor_col in factor_columns:
                value = row.get(factor_col, "").strip()
                if value:
                    # Replace spaces with underscores
                    value = value.replace(" ", "_")
                    factor_values.append(value)

            joined_factors = "__".join(factor_values) if factor_values else ""

            # Column 3: bioreplicate
            bioreplicate = row.get("bioreplicate", "").strip()

            # Column 4: data_type
            data_type = row.get("data_type", "").strip()

            # Write row
            writer.writerow([input_file, joined_factors, bioreplicate, data_type])

    print(f"Manifest written to {args.output}")


if __name__ == "__main__":
    main()

