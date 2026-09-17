#!/usr/bin/env python
"""Concatenate per-step VV CSVs staged as VV_in.csv*."""

import argparse
from pathlib import Path


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--assay_suffix", required=True)
    args = ap.parse_args()
    logs = sorted(Path.cwd().glob("VV_in.csv*"))
    out = Path(f"VV_log_final{args.assay_suffix}.csv")
    if not logs:
        out.write_text("component,sample_id,check_name,status,flag_code,message,details\n")
        return
    parts = []
    for i, log in enumerate(logs):
        text = log.read_text()
        parts.append(text if i == 0 else "".join(text.splitlines(True)[1:]))
    out.write_text("".join(parts))


if __name__ == "__main__":
    main()
