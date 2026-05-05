#!/usr/bin/env python3
"""MD5 manifests for NF_Proteomics post-processing.

- raw: files under <outdir>/RawData/.
- processed: selected published outputs under <outdir>.
"""

import os
import sys
import hashlib
import argparse
import fnmatch


def calculate_md5(filepath):
    md5_hash = hashlib.md5()
    actual_path = os.path.realpath(filepath) if os.path.islink(filepath) else filepath
    try:
        with open(actual_path, "rb") as f:
            for chunk in iter(lambda: f.read(4096), b""):
                md5_hash.update(chunk)
        return md5_hash.hexdigest()
    except OSError as e:
        sys.stderr.write(f"Error calculating MD5 for {filepath}: {e}\n")
        return "ERROR"


def collect_raw_files(outdir):
    """Published staging outputs: all regular files under RawData/."""
    raw_dir = os.path.join(outdir, "RawData")
    if not os.path.isdir(raw_dir):
        return []
    out = []
    for root, _, files in os.walk(raw_dir):
        for fn in files:
            fp = os.path.join(root, fn)
            if os.path.isfile(fp):
                out.append(fp)
    return sorted(out, key=lambda p: os.path.relpath(p, raw_dir).lower())


def collect_matches(base_dir, patterns):
    if not os.path.isdir(base_dir):
        return []
    out = []
    for root, _, files in os.walk(base_dir):
        for fn in files:
            if any(fnmatch.fnmatch(fn, pattern) for pattern in patterns):
                fp = os.path.join(root, fn)
                if os.path.isfile(fp):
                    out.append(fp)
    return sorted(out, key=lambda p: os.path.relpath(p, base_dir).lower())


def collect_first_match(base_dir, pattern):
    matches = collect_matches(base_dir, [pattern])
    return matches[:1]


def collect_processed_files(outdir, assay_suffix):
    out = []
    out.extend(collect_matches(os.path.join(outdir, "RawBeans"), ["*.zip"]))
    out.extend(collect_matches(os.path.join(outdir, "Metadata"), [
        "*.csv",
        "*.workflow",
        "manifest*.tsv",
        "experiment_annotation*.tsv",
    ]))
    out.extend(collect_matches(os.path.join(outdir, "Proteome"), ["*.fas"]))
    out.extend(collect_matches(os.path.join(outdir, "pmultiqc"), [
        "multiqc*.html",
        "multiqc*_data.zip",
    ]))
    out.extend(collect_matches(os.path.join(outdir, "MSstats"), [
        "msstats_comparison*.csv",
        "msstats_contrasts*.csv",
    ]))
    out.extend(collect_matches(os.path.join(outdir, "FragPipe"), [
        "combined_*.tsv",
        "msstats*.csv",
        "abundance_*.tsv",
        "ratio_*.tsv",
    ]))

    fpa_dir = os.path.join(outdir, "FragPipeAnalystR")
    out.extend(collect_matches(fpa_dir, [
        f"FragPipeAnalystR_parameters_*{assay_suffix}.txt",
        f"*matrix_*{assay_suffix}.csv",
        f"QC_plots_*{assay_suffix}.zip",
        f"comparison_plots_*{assay_suffix}.zip",
        f"pathway_analysis_plots_*{assay_suffix}.zip",
        f"DE_plots_*{assay_suffix}.zip",
        f"DE_results_*{assay_suffix}.csv",
    ]))
    out.extend(collect_first_match(fpa_dir, f"SampleTable{assay_suffix}.csv"))
    out.extend(collect_first_match(fpa_dir, f"contrasts{assay_suffix}.csv"))

    out.extend(collect_matches(os.path.join(outdir, "GeneLab"), [
        f"software_versions{assay_suffix}.md",
        f"processing_info{assay_suffix}.zip",
    ]))

    seen = set()
    unique = []
    for fp in sorted(out, key=lambda p: os.path.relpath(p, outdir).lower()):
        key = os.path.abspath(fp)
        if key in seen:
            continue
        seen.add(key)
        unique.append(fp)
    return unique


def dedup_by_basename(filename):
    seen = set()
    lines = []
    try:
        with open(filename, "r") as f:
            for line in f:
                if not line.strip():
                    continue
                key = line.split("\t", 1)[0]
                if key not in seen:
                    seen.add(key)
                    lines.append(line)
    except FileNotFoundError:
        return
    with open(filename, "w") as f:
        f.writelines(lines)


def main():
    parser = argparse.ArgumentParser(description="Generate MD5 TSVs for GeneLab proteomics outputs.")
    parser.add_argument("--outdir", required=True, help="Completed workflow output root (contains RawData/, etc.)")
    parser.add_argument("--assay_suffix", default="", help="e.g. _GLProteomics")
    args = parser.parse_args()

    outdir = os.path.abspath(args.outdir)
    raw_md5_file = f"raw_md5sum{args.assay_suffix}.tsv"
    processed_md5_file = f"processed_md5sum{args.assay_suffix}.tsv"

    raw_paths = collect_raw_files(outdir)
    raw_count = 0
    with open(raw_md5_file, "w") as f:
        for filepath in raw_paths:
            md5sum = calculate_md5(filepath)
            f.write(f"{os.path.basename(filepath)}\t{md5sum}\n")
            raw_count += 1

    print(f"Scanning (processed): {outdir}")
    processed_paths = collect_processed_files(outdir, args.assay_suffix)
    processed_lines = []
    for filepath in processed_paths:
        md5sum = calculate_md5(filepath)
        processed_lines.append(f"{os.path.basename(filepath)}\t{md5sum}\n")

    processed_count = len(processed_lines)
    with open(processed_md5_file, "w") as f:
        f.writelines(processed_lines)

    print(f"Raw (RawData/*): {raw_count} rows -> {raw_md5_file}")
    print(f"Processed: {processed_count} rows -> {processed_md5_file}")

    dedup_by_basename(raw_md5_file)
    dedup_by_basename(processed_md5_file)


if __name__ == "__main__":
    main()
