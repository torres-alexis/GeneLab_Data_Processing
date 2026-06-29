#!/usr/bin/env python3
"""Clean filesystem paths from published proteomics outputs (post-processing).

Targets:
  - processing_info/nextflow*.txt
  - FragPipeAnalystR/**/FragPipeAnalystR_parameters*.txt
  - pmultiqc/multiqc*.html and multiqc*_data.zip
  - FragPipe/fragpipe*.zip (text configs/logs inside the archive)
"""

from __future__ import annotations

import argparse
import fnmatch
import json
import os
import re
import shutil
import sys
import tempfile
import zipfile
from pathlib import Path

# --- shared text scrubbing (fragpipe zip internals) ---

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

FRAGPIPE_SANITIZE_BASENAMES = frozenset(
    {
        "filelist_ionquant.txt",
        "filelist_proteinprophet.txt",
        "sdrf.tsv",
        "experiment_annotation.tsv",
        "experiment_annotation_GLProteomics.tsv",
        "msbooster_params.txt",
    }
)

PATH_REMOVED = "<path-removed-for-security-purposes>"


def _basename_path(path: str) -> str:
    return os.path.basename(path.rstrip("/")) or path


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


def _scrub_fragpipe_line(line: str) -> str:
    line = _UNIX_ABS_PATH.sub(_scrub_absolute_path, line)
    if line.startswith("workdir="):
        suffix = "\n" if line.endswith("\n") else ""
        return f"workdir=.{suffix}"
    return line


def _scrub_experiment_annotation_line(line: str) -> str:
    if not line.strip() or line.lower().startswith("file\t"):
        return line
    parts = line.rstrip("\n").split("\t")
    if parts:
        parts[0] = Path(parts[0]).name
    return "\t".join(parts) + ("\n" if line.endswith("\n") else "")


def _should_scrub_fragpipe_member(name: str) -> bool:
    path = Path(name)
    if path.name in FRAGPIPE_SANITIZE_BASENAMES:
        return True
    if path.name.startswith("log_") and path.suffix == ".txt":
        return True
    if path.name.endswith(".fp-manifest"):
        return True
    if path.suffix.lower() == ".log":
        return True
    return path.suffix in {".params", ".workflow", ".job", ".yml"}


def _scrub_fragpipe_text(name: str, text: str) -> str:
    lines = text.splitlines(keepends=True)
    if Path(name).name.startswith("experiment_annotation") and Path(name).suffix == ".tsv":
        return "".join(_scrub_experiment_annotation_line(line) for line in lines)
    return "".join(_scrub_fragpipe_line(line) for line in lines)


def _scrub_nextflow_log_text(text: str) -> str:
    """Port of clean_paths.sh for processing_info nextflow logs."""
    rules = [
        (r"/[A-Za-z0-9._/-]*/(OSD-[0-9][A-Za-z0-9._-]*/RawData/[^\" )]*)", r"\1"),
        (r"/[A-Za-z0-9._/-]*/(OSD-[0-9][A-Za-z0-9._-]*/results/[^\" )]*)", r"\1"),
        (r"/[A-Za-z0-9._/-]*/(GLDS-[0-9][A-Za-z0-9._-]*/RawData/[^\" )]*)", r"\1"),
        (r"/[A-Za-z0-9._/-]*/(GLDS-[0-9][A-Za-z0-9._-]*/results/[^\" )]*)", r"\1"),
        (r"/[A-Za-z0-9._/-]*/(workflow_code/[^\" )]*)", r"\1"),
        (r"/[A-Za-z0-9._/-]*/main\.nf", "main.nf"),
        (r"/[A-Za-z0-9._/-]*/work/[A-Za-z0-9][A-Za-z0-9._/-]*", "<workdir>"),
        (r'"/home/[^"]*"', f'"{PATH_REMOVED}"'),
        (r"'/home/[^']*'", f"'{PATH_REMOVED}'"),
        (r"/home/[^\" )]*", PATH_REMOVED),
        (r'"/global/[^"]*"', f'"{PATH_REMOVED}"'),
        (r"'/global/[^']*'", f"'{PATH_REMOVED}'"),
        (r"/global/[^\" )]*", PATH_REMOVED),
    ]
    for pattern, repl in rules:
        text = re.sub(pattern, repl, text)
    return text


def _scrub_parameter_line(line: str) -> str:
    if ":" not in line:
        return _scrub_fragpipe_line(line)
    prefix, _, value = line.partition(":")
    value = value.strip()
    if value.startswith(("/", "./")) or "/home/" in value or "/global/" in value:
        if "/" in value:
            value = _basename_path(value.split()[0]) + value[len(value.split()[0]) :]
    return f"{prefix}: {value}\n" if line.endswith("\n") else f"{prefix}: {value}"


# --- per-target strippers ---

def strip_processing_info(processed_dir: Path, assay_suffix: str) -> list[str]:
    touched: list[str] = []
    info_dir = processed_dir / "processing_info"
    if not info_dir.is_dir():
        return touched

    for path in sorted(info_dir.glob("nextflow*.txt")):
        original = path.read_text(encoding="utf-8", errors="replace")
        cleaned = _scrub_nextflow_log_text(original)
        if cleaned != original:
            path.write_text(cleaned, encoding="utf-8")
            touched.append(str(path))

    for path in sorted(processed_dir.rglob(f"FragPipeAnalystR_parameters_*{assay_suffix}.txt")):
        original = path.read_text(encoding="utf-8", errors="replace")
        cleaned = "".join(_scrub_parameter_line(line) for line in original.splitlines(keepends=True))
        if cleaned != original:
            path.write_text(cleaned, encoding="utf-8")
            touched.append(str(path))

    return touched


def _clean_multiqc_json(data: dict) -> dict:
    for key in ("config_output_dir", "config_script_path"):
        if key in data and isinstance(data[key], str):
            data[key] = _basename_path(data[key])
    if "config_analysis_dir_abs" in data and isinstance(data["config_analysis_dir_abs"], list):
        data["config_analysis_dir_abs"] = [_basename_path(p) for p in data["config_analysis_dir_abs"]]
    if "config_analysis_dir" in data and isinstance(data["config_analysis_dir"], list):
        data["config_analysis_dir"] = [_basename_path(p) for p in data["config_analysis_dir"]]
    if "report_multiqc_command" in data and isinstance(data["report_multiqc_command"], str):
        data["report_multiqc_command"] = _UNIX_ABS_PATH.sub(_scrub_absolute_path, data["report_multiqc_command"])
    return data


def _clean_multiqc_log_line(line: str) -> str:
    line = _UNIX_ABS_PATH.sub(_scrub_absolute_path, line)
    return line


def _clean_multiqc_html_line(line: str) -> str:
    pattern = r'(<code class="mqc_analysis_path">)(.*?)(</code>)'
    return re.sub(
        pattern,
        lambda m: f"{m.group(1)}{_basename_path(m.group(2))}{m.group(3)}",
        line,
    )


def strip_multiqc(processed_dir: Path, assay_suffix: str) -> list[str]:
    touched: list[str] = []
    pmultiqc_dir = processed_dir / "pmultiqc"
    if not pmultiqc_dir.is_dir():
        return touched

    for html in sorted(pmultiqc_dir.glob(f"multiqc*{assay_suffix}.html")):
        lines = html.read_text(encoding="utf-8", errors="replace").splitlines(keepends=True)
        cleaned = "".join(_clean_multiqc_html_line(line) for line in lines)
        if cleaned != "".join(lines):
            html.write_text(cleaned, encoding="utf-8")
            touched.append(str(html))

    for zip_path in sorted(pmultiqc_dir.glob(f"multiqc*{assay_suffix}_data.zip")):
        tmp_zip = zip_path.with_suffix(".zip.tmp")
        with zipfile.ZipFile(zip_path, "r") as zf_in:
            with zipfile.ZipFile(tmp_zip, "w", compression=zipfile.ZIP_DEFLATED) as zf_out:
                for member in zf_in.infolist():
                    if member.is_dir():
                        continue
                    name = member.filename
                    payload = zf_in.read(name)
                    if name.endswith("multiqc_data.json"):
                        obj = json.loads(payload.decode("utf-8"))
                        payload = (json.dumps(_clean_multiqc_json(obj), indent=4) + "\n").encode("utf-8")
                    elif name.endswith("multiqc.log"):
                        text = payload.decode("utf-8", errors="replace")
                        payload = "".join(
                            _clean_multiqc_log_line(line)
                            for line in text.splitlines(keepends=True)
                        ).encode("utf-8")
                    elif name.endswith("multiqc_sources.txt"):
                        lines = payload.decode("utf-8", errors="replace").splitlines(keepends=True)
                        out = []
                        for line in lines:
                            if line.strip():
                                parts = line.split("\t")
                                if len(parts) >= 4:
                                    parts[3] = _basename_path(parts[3])
                                out.append("\t".join(parts) + ("\n" if line.endswith("\n") else ""))
                            else:
                                out.append(line)
                        payload = "".join(out).encode("utf-8")
                    elif name.endswith(".html"):
                        lines = payload.decode("utf-8", errors="replace").splitlines(keepends=True)
                        payload = "".join(_clean_multiqc_html_line(line) for line in lines).encode("utf-8")
                    zf_out.writestr(member, payload)
        tmp_zip.replace(zip_path)
        touched.append(str(zip_path))

    return touched


def strip_fragpipe_zip(processed_dir: Path, assay_suffix: str) -> list[str]:
    touched: list[str] = []
    fragpipe_dir = processed_dir / "FragPipe"
    if not fragpipe_dir.is_dir():
        return touched

    for zip_path in sorted(fragpipe_dir.glob(f"fragpipe{assay_suffix}.zip")):
        with tempfile.TemporaryDirectory(prefix="strip_fragpipe_") as tmp:
            tmp_dir = Path(tmp)
            with zipfile.ZipFile(zip_path, "r") as zf:
                for member in zf.infolist():
                    if member.is_dir():
                        continue
                    data = zf.read(member.filename)
                    out_path = tmp_dir / member.filename
                    out_path.parent.mkdir(parents=True, exist_ok=True)
                    if _should_scrub_fragpipe_member(member.filename):
                        text = data.decode("utf-8", errors="replace")
                        out_path.write_text(_scrub_fragpipe_text(member.filename, text), encoding="utf-8")
                    else:
                        out_path.write_bytes(data)
            zip_path.unlink()
            with zipfile.ZipFile(zip_path, "w", compression=zipfile.ZIP_DEFLATED) as zf:
                for file_path in sorted(tmp_dir.rglob("*")):
                    if file_path.is_file():
                        zf.write(file_path, arcname=file_path.relative_to(tmp_dir).as_posix())
        touched.append(str(zip_path))

    return touched


def strip_all(processed_dir: Path, assay_suffix: str) -> dict[str, list[str]]:
    return {
        "processing_info": strip_processing_info(processed_dir, assay_suffix),
        "multiqc": strip_multiqc(processed_dir, assay_suffix),
        "fragpipe_zip": strip_fragpipe_zip(processed_dir, assay_suffix),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--processed-dir", required=True, type=Path)
    parser.add_argument("--assay-suffix", default="_GLProteomics")
    parser.add_argument(
        "--clean-paths",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Scrub paths from published artifacts (default: true). Use --no-clean-paths to no-op.",
    )
    parser.add_argument(
        "--targets",
        default="all",
        help="Comma-separated: all, processing_info, multiqc, fragpipe_zip",
    )
    args = parser.parse_args()

    if not args.clean_paths:
        print("clean_paths disabled; skipping path cleanup", file=sys.stderr)
        return

    if not args.processed_dir.is_dir():
        raise SystemExit(f"Processed directory not found: {args.processed_dir}")

    targets = {t.strip() for t in args.targets.split(",") if t.strip()}
    if "all" in targets:
        targets = {"processing_info", "multiqc", "fragpipe_zip"}

    results: dict[str, list[str]] = {}
    if "processing_info" in targets:
        results["processing_info"] = strip_processing_info(args.processed_dir, args.assay_suffix)
    if "multiqc" in targets:
        results["multiqc"] = strip_multiqc(args.processed_dir, args.assay_suffix)
    if "fragpipe_zip" in targets:
        results["fragpipe_zip"] = strip_fragpipe_zip(args.processed_dir, args.assay_suffix)

    for name, paths in results.items():
        print(f"{name}: {len(paths)} file(s) updated", file=sys.stderr)
        for path in paths:
            print(f"  {path}", file=sys.stderr)


if __name__ == "__main__":
    main()
