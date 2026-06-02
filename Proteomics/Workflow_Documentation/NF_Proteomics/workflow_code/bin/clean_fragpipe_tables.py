#!/usr/bin/env python3
"""Rename FragPipe table sample columns for published outputs."""

import argparse
import csv
import sys
from pathlib import Path


def read_sample_display_map(experiment_annotation: Path) -> dict[str, str]:
    with experiment_annotation.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if "sample" not in (reader.fieldnames or []) or "sample_name" not in (reader.fieldnames or []):
            raise SystemExit("experiment_annotation must contain 'sample' and 'sample_name' columns")

        sample_display_map = {}
        for row in reader:
            sample = (row.get("sample") or "").strip()
            sample_name = (row.get("sample_name") or "").strip() or sample
            if sample and sample not in sample_display_map:
                sample_display_map[sample] = sample_name
        return sample_display_map


def rename_column(column: str, sample_display_map: dict[str, str]) -> str:
    for sample, display in sample_display_map.items():
        if column == sample:
            return display
        if column.startswith(f"{sample} ") or column.startswith(f"{sample}."):
            return display + column[len(sample) :]
    return column


def clean_table(input_table: Path, output_dir: Path, sample_display_map: dict[str, str]) -> Path:
    with input_table.open(newline="", encoding="utf-8") as in_handle:
        reader = csv.reader(in_handle, delimiter="\t")
        try:
            header = next(reader)
        except StopIteration:
            print(f"WARNING: {input_table} is empty; skipping cleaned table", file=sys.stderr)
            return None

        output_dir.mkdir(parents=True, exist_ok=True)
        output_table = output_dir / input_table.name
        with output_table.open("w", newline="", encoding="utf-8") as out_handle:
            writer = csv.writer(out_handle, delimiter="\t", lineterminator="\n")
            writer.writerow([rename_column(column, sample_display_map) for column in header])
            writer.writerows(reader)

    return output_table


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--experiment_annotation", required=True, type=Path)
    parser.add_argument("--output_dir", default=".", type=Path)
    parser.add_argument("--type", required=True, choices=("LFQ", "TMT"))
    parser.add_argument("tables", nargs="+", type=Path)
    args = parser.parse_args()

    sample_display_map = read_sample_display_map(args.experiment_annotation)
    for table in args.tables:
        clean_table(table, args.output_dir, sample_display_map)


if __name__ == "__main__":
    main()
