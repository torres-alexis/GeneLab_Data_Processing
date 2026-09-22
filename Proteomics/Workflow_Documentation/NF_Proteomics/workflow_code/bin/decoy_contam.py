#!/usr/bin/env python
"""Drop decoy / contaminant rows from FragPipe and DE tables."""

from __future__ import annotations

import argparse
import csv
import re
import sys
from pathlib import Path

ID_COLUMNS = (
    "ProteinName",
    "Protein ID",
    "Protein.ID",
    "ProteinID",
    "Protein",
    "Entry Name",
    "Entry.Name",
    "Index",
    "name",
    "ID",
)

GROUP_SPLIT_RE = re.compile(r"[;,]+")


def group_parts(value) -> list[str]:
    if value is None:
        return []
    text = str(value).strip()
    if not text or text.lower() in {"na", "nan", "none"}:
        return []
    return [p.strip() for p in GROUP_SPLIT_RE.split(text) if p.strip()]


def token_is_junk(
    token: str,
    *,
    decoy_prefix: str = "rev_",
    contam_prefix: str = "contam_",
) -> bool:
    t = token.strip()
    if not t:
        return False
    if decoy_prefix and t.startswith(decoy_prefix):
        return True
    if contam_prefix and t.startswith(contam_prefix):
        return True
    return False


def value_is_junk(value, **kwargs) -> bool:
    parts = group_parts(value)
    if not parts:
        return False
    if all(token_is_junk(p, **kwargs) for p in parts):
        return True
    return token_is_junk(parts[0], **kwargs)


def row_id_values(row: dict, columns=None) -> list[str]:
    cols = columns or [c for c in ID_COLUMNS if c in row]
    if not cols:
        cols = [c for c in row if c]
    return [str(row[c]) for c in cols if row.get(c) not in (None, "")]


def row_is_junk(row: dict, columns=None, **kwargs) -> bool:
    return any(value_is_junk(v, **kwargs) for v in row_id_values(row, columns))


def _open_table(path: Path):
    raw = path.read_bytes()
    if raw.startswith(b"\xef\xbb\xbf"):
        raw = raw[3:]
    text = raw.decode("utf-8", errors="replace")
    dialect = csv.Sniffer().sniff(text.splitlines()[0] + "\n", delimiters=",\t")
    return text, dialect


def filter_table(in_path: Path, out_path: Path, **kwargs) -> tuple[int, int]:
    text, dialect = _open_table(in_path)
    reader = csv.DictReader(text.splitlines(), dialect=dialect)
    rows = list(reader)
    fieldnames = reader.fieldnames or []
    kept = [r for r in rows if not row_is_junk(r, **kwargs)]
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fieldnames,
            extrasaction="ignore",
            dialect=dialect,
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(kept)
    return len(rows), len(kept)


def scan_table(path: Path, **kwargs) -> list[dict]:
    text, dialect = _open_table(path)
    reader = csv.DictReader(text.splitlines(), dialect=dialect)
    hits = []
    for i, row in enumerate(reader, start=2):
        if row_is_junk(row, **kwargs):
            hits.append({"line": i, "ids": row_id_values(row)[:3]})
    return hits


def _self_test() -> None:
    assert token_is_junk("rev_sp|P12345|X_MOUSE")
    assert not token_is_junk("REV_sp|P12345|X_MOUSE")
    assert not token_is_junk("decoy_sp|P12345|X_MOUSE")
    assert token_is_junk("contam_sp|P00761|TRYP_PIG")
    assert not token_is_junk("Cont_sp|P00761|TRYP_PIG")
    assert not token_is_junk("CONTAM_sp|P00761|TRYP_PIG")
    assert not token_is_junk("CRAP_sp|P00761|TRYP_PIG")
    assert not token_is_junk("crap_sp|P00761|TRYP_PIG")
    assert not token_is_junk("sp|P00761|TRYP_PIG")
    assert not token_is_junk("P02769")
    assert not token_is_junk("ALBU_BOVIN")
    assert not token_is_junk("P02768")
    assert not value_is_junk("sp|A0JNU3|LPP60_MOUSE;contam_sp|P00761|TRYP_PIG")
    assert value_is_junk("contam_sp|P00761|TRYP_PIG;sp|A0JNU3|LPP60_MOUSE")
    print("decoy_contam.py self-test OK")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--input", help="Table to filter (CSV/TSV)")
    ap.add_argument("--output", help="Filtered table path")
    ap.add_argument("--scan", help="Scan only; exit 1 if junk rows remain")
    ap.add_argument("--decoy-prefix", default="rev_")
    ap.add_argument("--contam-prefix", default="contam_")
    ap.add_argument("--self-test", action="store_true")
    args = ap.parse_args()
    if args.self_test:
        _self_test()
        return
    kwargs = {
        "decoy_prefix": args.decoy_prefix,
        "contam_prefix": args.contam_prefix,
    }
    if args.scan:
        hits = scan_table(Path(args.scan), **kwargs)
        if hits:
            print(f"FAIL: {len(hits)} decoy/contam row(s) in {args.scan}", file=sys.stderr)
            for hit in hits[:10]:
                print(f"  line {hit['line']}: {hit['ids']}", file=sys.stderr)
            sys.exit(1)
        print(f"OK: no decoy/contam rows in {args.scan}")
        return
    if not args.input or not args.output:
        ap.error("--input and --output required unless --scan or --self-test")
    n_in, n_out = filter_table(Path(args.input), Path(args.output), **kwargs)
    print(f"decoy_contam: {n_in} -> {n_out} ({n_in - n_out} dropped)", file=sys.stderr)


if __name__ == "__main__":
    main()
