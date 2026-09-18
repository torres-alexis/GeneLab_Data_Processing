#!/usr/bin/env python3
"""Inspect decoys and unlabeled accessions; strip decoys before Philosopher --custom."""

from __future__ import annotations

import argparse
import tempfile
from pathlib import Path


def parse_fasta(path: Path):
    header = None
    seq = []
    with path.open() as fh:
        for line in fh:
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(seq)
                header = line[1:].rstrip("\n")
                seq = []
            else:
                seq.append(line.rstrip("\n"))
        if header is not None:
            yield header, "".join(seq)


def write_fasta(path: Path, records) -> int:
    n = 0
    with path.open("w") as out:
        for header, seq in records:
            out.write(f">{header}\n{seq}\n")
            n += 1
    return n


def is_decoy(header: str, decoy_prefix: str) -> bool:
    return bool(decoy_prefix) and header.startswith(decoy_prefix)


def is_contam_tagged(header: str, decoy_prefix: str, contam_prefix: str) -> bool:
    if not contam_prefix:
        return False
    if decoy_prefix and header.startswith(decoy_prefix):
        header = header[len(decoy_prefix) :]
    return header.startswith(contam_prefix)


def strip_tags(header: str, decoy_prefix: str, contam_prefix: str) -> str:
    if decoy_prefix and header.startswith(decoy_prefix):
        header = header[len(decoy_prefix) :]
    if contam_prefix and header.startswith(contam_prefix):
        header = header[len(contam_prefix) :]
    return header


def accession(header: str, decoy_prefix: str = "rev_", contam_prefix: str = "contam_") -> str | None:
    parts = strip_tags(header, decoy_prefix, contam_prefix).split("|")
    if len(parts) >= 2 and parts[1]:
        return parts[1].split()[0]
    return None


def inspect_fasta(path: Path, decoy_prefix: str, contam_prefix: str) -> dict:
    n = n_decoy = n_tagged = 0
    unlabeled = set()
    for header, _seq in parse_fasta(path):
        n += 1
        decoy = is_decoy(header, decoy_prefix)
        tagged = is_contam_tagged(header, decoy_prefix, contam_prefix)
        if decoy:
            n_decoy += 1
        if tagged:
            n_tagged += 1
        acc = accession(header, decoy_prefix, contam_prefix)
        if acc and not decoy and not tagged:
            unlabeled.add(acc)
    return {
        "n": n,
        "n_decoy": n_decoy,
        "n_tagged": n_tagged,
        "has_decoys": n_decoy > 0,
        "unlabeled_accs": unlabeled,
    }


def write_targets(src: Path, dest: Path, decoy_prefix: str) -> int:
    return write_fasta(
        dest,
        ((h, s) for h, s in parse_fasta(src) if not is_decoy(h, decoy_prefix)),
    )


def restore_study(
    src: Path, dest: Path, study_accs: set[str], decoy_prefix: str, contam_prefix: str
) -> tuple[int, int]:
    if not contam_prefix or not study_accs:
        dest.write_text(src.read_text())
        return sum(1 for _ in parse_fasta(src)), 0
    n = n_restored = 0
    records = []
    for header, seq in parse_fasta(src):
        acc = accession(header, decoy_prefix, contam_prefix)
        if acc and acc in study_accs and is_contam_tagged(header, decoy_prefix, contam_prefix):
            decoy = decoy_prefix if decoy_prefix and header.startswith(decoy_prefix) else ""
            rest = header[len(decoy) :] if decoy else header
            if rest.startswith(contam_prefix):
                header = decoy + rest[len(contam_prefix) :]
                n_restored += 1
        records.append((header, seq))
        n += 1
    write_fasta(dest, records)
    return n, n_restored


def emit_env(path: Path, need_decoy: bool) -> None:
    path.write_text(f"NEED_DECOY={'true' if need_decoy else 'false'}\n")


def write_accs(accs: set[str], path: Path) -> None:
    path.write_text("".join(f"{a}\n" for a in sorted(accs)))


def load_accs(path: Path) -> set[str]:
    return {ln.strip() for ln in path.read_text().splitlines() if ln.strip()}


def _self_test() -> None:
    assert accession("sp|P02768|ALBU_HUMAN Albumin") == "P02768"
    assert accession("contam_sp|P00761|TRYP_PIG") == "P00761"
    assert accession("rev_contam_sp|P02768|ALBU_HUMAN") == "P02768"
    assert is_decoy("rev_sp|P02768|ALBU_HUMAN", "rev_")
    assert not is_decoy("sp|P02768|ALBU_HUMAN", "rev_")
    assert not is_decoy("sp|P12345|PREV_HUMAN", "rev_")
    assert is_contam_tagged("contam_sp|P00761|TRYP_PIG", "rev_", "contam_")
    assert not is_contam_tagged("sp|P02768|ALBU_HUMAN", "rev_", "contam_")

    with tempfile.TemporaryDirectory() as td:
        td = Path(td)
        src = td / "in.fas"
        src.write_text(
            ">sp|P02768|ALBU_HUMAN Albumin OS=Homo sapiens OX=9606\nMKWVTFISLL\n"
            ">sp|P04264|K2C1_HUMAN Keratin OS=Homo sapiens OX=9606\nMSRQFSSRSF\n"
            ">rev_sp|P12345|FAKE_HUMAN decoy\nIHGFEDCBAM\n"
            ">contam_sp|P00761|TRYP_PIG Trypsin OS=Sus scrofa\nIVGGYTCAAN\n"
            ">sp|P12345|FAKE_HUMAN not cRAP\nMABCDEFGHI\n"
        )
        st = inspect_fasta(src, "rev_", "contam_")
        assert st["n"] == 5
        assert st["n_decoy"] == 1
        assert st["has_decoys"] is True
        assert "P02768" in st["unlabeled_accs"]
        assert "P00761" not in st["unlabeled_accs"]

        targets = td / "targets.fas"
        n = write_targets(src, targets, "rev_")
        text = targets.read_text()
        assert n == 4
        assert ">rev_sp|P12345|" not in text
        assert ">contam_sp|P00761|TRYP_PIG" in text
        assert ">sp|P02768|ALBU_HUMAN" in text

        phi = td / "phi.fas"
        phi.write_text(
            ">contam_sp|P02768|ALBU_HUMAN Albumin OS=Homo sapiens OX=9606\nMKWVTFISLL\n"
            ">rev_contam_sp|P02768|ALBU_HUMAN Albumin OS=Homo sapiens OX=9606\nLLSIFTVWKM\n"
            ">contam_sp|P00761|TRYP_PIG Trypsin OS=Sus scrofa\nIVGGYTCAAN\n"
            ">rev_contam_sp|P00761|TRYP_PIG Trypsin OS=Sus scrofa\nNAACTYGGVI\n"
            ">sp|P12345|FAKE_HUMAN not cRAP\nMABCDEFGHI\n"
        )
        out = td / "out.fas"
        n, restored = restore_study(phi, out, st["unlabeled_accs"], "rev_", "contam_")
        assert n == 5
        assert restored == 2
        text = out.read_text()
        assert ">sp|P02768|ALBU_HUMAN" in text
        assert ">rev_sp|P02768|ALBU_HUMAN" in text
        assert ">contam_sp|P00761|TRYP_PIG" in text
        assert ">rev_contam_sp|P00761|TRYP_PIG" in text
        assert ">contam_sp|P02768|" not in text
    print("add_decoys_contams.py self-test OK")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("command", choices=["inspect", "write-targets", "restore-study", "self-test"])
    ap.add_argument("fasta", nargs="?", type=Path)
    ap.add_argument("--decoy-prefix", default="rev_")
    ap.add_argument("--contam-prefix", default="contam_")
    ap.add_argument("--want-decoys", action="store_true")
    ap.add_argument("--emit-env", type=Path)
    ap.add_argument("--emit-study-accs", type=Path)
    ap.add_argument("--study-accs", type=Path)
    ap.add_argument("--output", type=Path)
    args = ap.parse_args()

    if args.command == "self-test":
        _self_test()
        return

    if args.fasta is None:
        ap.error("fasta required")

    decoy_prefix = args.decoy_prefix or "rev_"
    contam_prefix = args.contam_prefix or ""

    if args.command == "inspect":
        st = inspect_fasta(args.fasta, decoy_prefix, contam_prefix)
        need_decoy = args.want_decoys and not st["has_decoys"]
        print(
            f"add_decoys_contams: n={st['n']} decoys={st['n_decoy']} tagged={st['n_tagged']} "
            f"NEED_DECOY={str(need_decoy).lower()}",
            flush=True,
        )
        if args.emit_env:
            emit_env(args.emit_env, need_decoy)
        if args.emit_study_accs:
            write_accs(st["unlabeled_accs"], args.emit_study_accs)
        return

    if args.command == "write-targets":
        if not args.output:
            ap.error("--output required")
        n = write_targets(args.fasta, args.output, decoy_prefix)
        print(f"write-targets: {n} target sequences", flush=True)
        return

    if args.command == "restore-study":
        if not args.output or not args.study_accs:
            ap.error("--output and --study-accs required")
        n, restored = restore_study(
            args.fasta, args.output, load_accs(args.study_accs), decoy_prefix, contam_prefix
        )
        print(f"restore-study: {n} sequences, stripped contam_ from {restored} study accessions", flush=True)
        return


if __name__ == "__main__":
    main()
