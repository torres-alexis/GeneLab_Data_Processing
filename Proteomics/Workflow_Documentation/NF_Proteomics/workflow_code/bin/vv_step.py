#!/usr/bin/env python
"""Fail if any listed path is missing or empty."""

import argparse
import sys
from pathlib import Path


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--name", required=True)
    ap.add_argument("--suffix", default="")
    ap.add_argument("paths", nargs="+")
    args = ap.parse_args()
    bad = []
    for raw in args.paths:
        path = Path(raw)
        if not path.is_file() or path.stat().st_size == 0:
            bad.append(raw)
    if bad:
        sys.exit(f"VV_STEP {args.name}: empty/missing: {', '.join(bad)}")
    msg = f"VV_STEP {args.name}: ok ({len(args.paths)} file(s))"
    print(msg)
    Path(f"VV_log_{args.name}{args.suffix}.log").write_text(msg + "\n")


if __name__ == "__main__":
    main()
