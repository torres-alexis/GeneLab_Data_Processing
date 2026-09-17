#!/usr/bin/env python
"""Keep RED/HALT rows from VV_log_final (drop flag_code 20 and 30)."""

import argparse
import csv
from pathlib import Path


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--assay_suffix", required=True)
    args = ap.parse_args()
    inp = Path(f"VV_log_final{args.assay_suffix}.csv")
    out = Path(f"VV_log_final_only_issues{args.assay_suffix}.csv")
    with inp.open(newline="") as f:
        reader = csv.DictReader(f)
        fieldnames = reader.fieldnames
        rows = [r for r in reader if str(r.get("flag_code", "")) not in ("20", "30")]
    with out.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    main()
