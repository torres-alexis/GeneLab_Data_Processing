#!/usr/bin/env python
"""Fail if any listed path is missing or empty."""

import argparse
import sys
from pathlib import Path


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--name", required=True)
    ap.add_argument("paths", nargs="+")
    args = ap.parse_args()
    bad = []
    for raw in args.paths:
        path = Path(raw)
        if not path.is_file() or path.stat().st_size == 0:
            bad.append(raw)
    if bad:
        sys.exit(f"VV_STEP {args.name}: empty/missing: {', '.join(bad)}")
    print(f"VV_STEP {args.name}: ok ({len(args.paths)} file(s))")


if __name__ == "__main__":
    main()
