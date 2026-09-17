#!/usr/bin/env python
"""Fail if any listed path is missing or empty."""

import argparse
import csv
import sys
from pathlib import Path

HEADER = ["component", "sample_id", "check_name", "status", "flag_code", "message", "details"]
FLAG = {"GREEN": "20", "YELLOW": "30", "RED": "50", "HALT": "80"}


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--name", required=True)
    ap.add_argument("--suffix", default="")
    ap.add_argument("paths", nargs="+")
    args = ap.parse_args()

    rows = []
    bad = []
    for raw in args.paths:
        path = Path(raw)
        sid = path.name
        if not path.is_file() or path.stat().st_size == 0:
            bad.append(raw)
            rows.append(
                {
                    "component": args.name,
                    "sample_id": sid,
                    "check_name": "check_output_existence",
                    "status": "HALT",
                    "flag_code": FLAG["HALT"],
                    "message": "empty or missing",
                    "details": raw,
                }
            )
        else:
            rows.append(
                {
                    "component": args.name,
                    "sample_id": sid,
                    "check_name": "check_output_existence",
                    "status": "GREEN",
                    "flag_code": FLAG["GREEN"],
                    "message": "ok",
                    "details": raw,
                }
            )

    out = Path(f"VV_log_{args.name}{args.suffix}.csv")
    with out.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=HEADER)
        writer.writeheader()
        writer.writerows(rows)

    if bad:
        sys.exit(f"VV_STEP {args.name}: empty/missing: {', '.join(bad)}")
    print(f"VV_STEP {args.name}: ok ({len(args.paths)} file(s))")


if __name__ == "__main__":
    main()
