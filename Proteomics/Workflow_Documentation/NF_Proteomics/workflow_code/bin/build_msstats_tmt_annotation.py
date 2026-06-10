#!/usr/bin/env python3
"""
Build MSstatsTMT annotation (7 columns) from TMT data_sheet + sample_sheet.

Output columns: Run, Fraction, TechRepMixture, Mixture, Channel, BioReplicate, Condition
  Run             = data_sheet run (must match msstats.csv Run / SpectrumName)
  Fraction        = data_sheet fraction (default 1)
  TechRepMixture  = data_sheet TechRepMixture (default 1)
  Mixture         = data_sheet plex (TMTa, TMT1, …)
  Channel         = sample_sheet channel
  BioReplicate    = sample_sheet Bioreplicate (default 1)
  Condition       = biology group from Factor Value[…] (R-safe); unused channels → Empty
"""

import argparse
import csv
import re
import sys


def _make_names_safe(s: str) -> str:
    if not s:
        return ""
    out = re.sub(r"[^a-zA-Z0-9._]", ".", s.strip())
    out = re.sub(r"\.+", ".", out).strip(".")
    if not out:
        return ""
    if out[0].isdigit():
        out = "X" + out
    return out


def _sanitize_for_fragpipe(s: str) -> str:
    if not s:
        return ""
    return re.sub(r"[^A-Za-z0-9_]", "_", s.strip())


def _condition_from_factors(row: dict, factor_columns: list) -> str:
    values = [row.get(c, "").strip() for c in factor_columns if row.get(c, "").strip()]
    if not values:
        return ""
    raw = " & ".join(values)
    return _make_names_safe(_sanitize_for_fragpipe(raw))


def main():
    parser = argparse.ArgumentParser(description="Build MSstatsTMT annotation from TMT sheets")
    parser.add_argument("--data_sheet", required=True)
    parser.add_argument("--sample_sheet", required=True)
    parser.add_argument(
        "--output",
        default="msstats_tmt_annotation.tsv",
        help="Output TSV path (default: msstats_tmt_annotation.tsv)",
    )
    args = parser.parse_args()

    with open(args.data_sheet, newline="") as f:
        data_reader = csv.DictReader(f)
        data_rows = [r for r in data_reader if (r.get("data_file") or "").strip()]
        if not data_rows:
            sys.exit("Error: data_sheet has no rows with data_file")

    with open(args.sample_sheet, newline="") as f:
        sample_reader = csv.DictReader(f)
        sample_rows = list(sample_reader)
        sample_fields = sample_reader.fieldnames or []
        if not sample_rows:
            sys.exit("Error: sample_sheet is empty")

    factor_cols = [c for c in sample_fields if c.startswith("Factor Value[")]

    plex_to_channels = {}
    for row in sample_rows:
        plex = (row.get("plex") or "").strip()
        channel = (row.get("channel") or "").strip()
        if not plex or not channel:
            continue
        sample_name = (row.get("Sample Name") or "").strip()
        biorep = (row.get("Bioreplicate") or row.get("Bioreplicate") or "1").strip() or "1"
        cond = _condition_from_factors(row, factor_cols)
        if not cond:
            cond = "Empty" if not sample_name else _make_names_safe(_sanitize_for_fragpipe(sample_name))
        plex_to_channels.setdefault(plex, []).append(
            {"channel": channel, "bioreplicate": biorep, "condition": cond}
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
        fraction = (row.get("fraction") or "1").strip() or "1"
        tech_rep = (row.get("TechRepMixture") or "1").strip() or "1"
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

    with open(args.output, "w", newline="\n") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(out_rows)

    print(f"MSstatsTMT annotation written to {args.output} ({len(out_rows)} rows)")


if __name__ == "__main__":
    main()
