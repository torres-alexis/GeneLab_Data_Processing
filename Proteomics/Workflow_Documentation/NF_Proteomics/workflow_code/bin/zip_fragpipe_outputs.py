#!/usr/bin/env python3
"""Stage, scrub paths, and zip FragPipe run outputs for publication.

Excluded: *.mzML (spectra under RawData/), search intermediates that are binary or
path-heavy and not needed for GeneLab deliverables (*.pepXML, *.pep.xml, *.pin, *.mzBIN).

Text configs and logs are path-scrubbed before zipping; tables, tmt-report/, and
Philosopher TSVs are included as-is.
"""

from __future__ import annotations

import argparse
import fnmatch
import re
import shutil
import sys
import zipfile
from dataclasses import dataclass
from pathlib import Path

# Glob patterns on basename (fnmatch).
EXCLUDE_PATTERNS = (
    "*.mzML",
    "*.pepXML",
    "*.pep.xml",
    "*.pin",
    "*.mzBIN",
)

EXCLUDE_DIRS: frozenset[str] = frozenset()

# Text files rewritten to strip absolute paths before zipping.
SANITIZE_BASENAMES = frozenset(
    {
        "filelist_ionquant.txt",
        "filelist_proteinprophet.txt",
        "sdrf.tsv",
        "experiment_annotation.tsv",
        "experiment_annotation_GLProteomics.tsv",
        "msbooster_params.txt",
    }
)

# Match Unix absolute paths inside log/filelist text.
_UNIX_ABS_PATH = re.compile(r"/(?:[^\s'\"\]\),]+(?:/[^\s'\"\]\),]+)*)")
_SKIP_PATH_PARTS = frozenset(
    {
        "work",
        "output",
        "fragpipe_temp",
        "tools",
        "ext",
        "RawData",
        "fragpipe_bin",
        "home",
        "git",
        "global",
    }
)


@dataclass
class ZipStats:
    kept_files: int = 0 
    kept_bytes: int = 0
    excluded_files: int = 0
    excluded_bytes: int = 0


def _scrub_absolute_path(match: re.Match[str]) -> str:
    raw = match.group(0).rstrip("],;")
    parts = [p for p in raw.split("/") if p]
    if not parts:
        return raw
    name = parts[-1]
    if name.endswith(".fas") or name.endswith(".mzML"):
        return name
    if len(parts) >= 2 and parts[-2] not in _SKIP_PATH_PARTS and not re.fullmatch(
        r"[0-9a-f]{32}", parts[-2]
    ):
        return f"{parts[-2]}/{name}"
    return name


def _scrub_line(line: str) -> str:
    line = _UNIX_ABS_PATH.sub(_scrub_absolute_path, line)
    if line.startswith("workdir="):
        suffix = "\n" if line.endswith("\n") else ""
        return f"workdir=.{suffix}"
    return line


def _matches_exclude(rel_path: Path) -> bool:
    if any(part in EXCLUDE_DIRS for part in rel_path.parts):
        return True
    name = rel_path.name
    return any(fnmatch.fnmatch(name, pat) for pat in EXCLUDE_PATTERNS)


def _should_scrub(path: Path) -> bool:
    if path.name in SANITIZE_BASENAMES:
        return True
    if path.name.startswith("log_") and path.suffix == ".txt":
        return True
    if path.name.endswith(".fp-manifest"):
        return True
    if path.suffix.lower() == ".log":
        return True
    return path.suffix in {".params", ".workflow", ".job", ".yml"}


def _scrub_experiment_annotation_line(line: str) -> str:
    if not line.strip() or line.lower().startswith("file\t"):
        return line
    parts = line.rstrip("\n").split("\t")
    if parts:
        parts[0] = Path(parts[0]).name
    return "\t".join(parts) + ("\n" if line.endswith("\n") else "")


def _scrub_file(src: Path, dest: Path) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    with src.open(encoding="utf-8", errors="replace") as handle:
        lines = list(handle)
    if src.name.startswith("experiment_annotation") and src.suffix == ".tsv":
        lines = [_scrub_experiment_annotation_line(line) for line in lines]
    else:
        lines = [_scrub_line(line) for line in lines]
    dest.write_text("".join(lines), encoding="utf-8")


def _stage_tree(input_dir: Path, staging_dir: Path) -> ZipStats:
    stats = ZipStats()
    if staging_dir.exists():
        shutil.rmtree(staging_dir)
    staging_dir.mkdir(parents=True)

    for src in sorted(input_dir.rglob("*")):
        if not src.is_file():
            continue
        rel = src.relative_to(input_dir)
        size = src.stat().st_size
        if _matches_exclude(rel):
            stats.excluded_files += 1
            stats.excluded_bytes += size
            continue
        dest = staging_dir / rel
        if _should_scrub(src):
            _scrub_file(src, dest)
        else:
            dest.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(src, dest)
        stats.kept_files += 1
        stats.kept_bytes += size

    return stats


def _make_zip(staging_dir: Path, output_zip: Path) -> None:
    output_zip.parent.mkdir(parents=True, exist_ok=True)
    if output_zip.exists():
        output_zip.unlink()
    with zipfile.ZipFile(output_zip, "w", compression=zipfile.ZIP_DEFLATED) as zf:
        for path in sorted(staging_dir.rglob("*")):
            if path.is_file():
                zf.write(path, arcname=path.relative_to(staging_dir).as_posix())


def zip_fragpipe_outputs(input_dir: Path, output_zip: Path, workflow_type: str) -> ZipStats:
    del workflow_type  # LFQ vs TMT share excludes; branch here if rules diverge.
    staging_dir = output_zip.parent / f".{output_zip.stem}_staging"
    stats = _stage_tree(input_dir, staging_dir)
    if stats.kept_files == 0:
        raise SystemExit(f"No publishable FragPipe files found under {input_dir}")
    _make_zip(staging_dir, output_zip)
    shutil.rmtree(staging_dir)
    return stats


def _format_bytes(num: int) -> str:
    if num >= 1_000_000_000:
        return f"{num / 1_000_000_000:.2f} GB"
    if num >= 1_000_000:
        return f"{num / 1_000_000:.1f} MB"
    if num >= 1_000:
        return f"{num / 1_000:.1f} KB"
    return f"{num} B"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, type=Path, help="FragPipe output/ directory")
    parser.add_argument("--output", required=True, type=Path, help="Output .zip path")
    parser.add_argument(
        "--workflow-type",
        required=True,
        choices=("LFQ", "TMT"),
        help="LFQ or TMT (reserved for workflow-specific exclude rules)",
    )
    args = parser.parse_args()

    if not args.input.is_dir():
        raise SystemExit(f"Input directory not found: {args.input}")

    stats = zip_fragpipe_outputs(args.input, args.output, args.workflow_type)
    zip_size = args.output.stat().st_size
    print(
        f"Wrote {args.output} ({_format_bytes(zip_size)}); "
        f"kept {stats.kept_files} files ({_format_bytes(stats.kept_bytes)} raw), "
        f"excluded {stats.excluded_files} files ({_format_bytes(stats.excluded_bytes)})",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
