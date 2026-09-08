#!/usr/bin/env python

"""
Validate GeneLab Proteomics processed datasets.
"""

import argparse
import glob
import os
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from decoy_contam import DECOY_PREFIXES, scan_table


def fail(message, log):
    log.write(f"FAIL: {message}\n")
    sys.stderr.write(f"Validation failed: {message}\n")
    sys.exit(1)


def warn(message, log):
    log.write(f"WARN: {message}\n")


def check_dir(path, log, required=True):
    if not os.path.isdir(path):
        if required:
            fail(f"missing directory: {path}", log)
        warn(f"missing directory: {path}", log)
        return False
    log.write(f"OK: directory found: {path}\n")
    return True


def check_nonempty_file(path, log, required=True):
    if not os.path.isfile(path):
        if required:
            fail(f"missing file: {path}", log)
        warn(f"missing file: {path}", log)
        return False
    if os.path.getsize(path) == 0:
        if required:
            fail(f"empty file: {path}", log)
        warn(f"empty file: {path}", log)
        return False
    log.write(f"OK: file found: {path}\n")
    return True


def check_nonempty_glob(pattern, log, required=True):
    matches = sorted(glob.glob(pattern, recursive=True))
    if not matches:
        if required:
            fail(f"missing file matching: {pattern}", log)
        warn(f"missing file matching: {pattern}", log)
        return []
    for path in matches:
        check_nonempty_file(path, log, required=required)
    return matches


def check_decoy_free(paths, log, decoy_prefix="rev_"):
    prefixes = tuple(dict.fromkeys((decoy_prefix, *DECOY_PREFIXES)))
    for path in paths:
        hits = scan_table(Path(path), decoy_prefixes=prefixes)
        if hits:
            fail(
                f"{len(hits)} decoy/contam row(s) remain in {path} "
                f"(first: line {hits[0]['line']} {hits[0]['ids']})",
                log,
            )
        log.write(f"OK: no decoy/contam rows in {path}\n")


def main():
    parser = argparse.ArgumentParser(description="Validate NF_Proteomics processed outputs.")
    parser.add_argument("--outdir", required=True, help="Completed workflow output root")
    parser.add_argument("--assay_suffix", default="", help="e.g. _GLProteomics")
    parser.add_argument("--output", default=None, help="Validation log filename")
    parser.add_argument(
        "--require-post-processing",
        default="true",
        help="Require GeneLab/ md5s + processing_info (default true for NF post-processing).",
    )
    parser.add_argument(
        "--check-decoys",
        default="true",
        help="Fail if published DE / MSstats tables still contain decoys or cRAP.",
    )
    parser.add_argument(
        "--decoy-prefix",
        default="rev_",
        help="Philosopher decoy prefix; VALIDATE also scans the leftover default prefixes.",
    )
    args = parser.parse_args()

    require_pp = str(args.require_post_processing).strip().lower() in {"1", "true", "yes"}
    check_decoys = str(args.check_decoys).strip().lower() in {"1", "true", "yes"}

    outdir = os.path.abspath(args.outdir)
    output = args.output or f"validate_processed_proteomics{args.assay_suffix}.log"
    suffix = args.assay_suffix

    with open(output, "w") as log:
        log.write("NF_Proteomics processed output validation\n")
        log.write(f"outdir: {outdir}\n")
        log.write(f"assay_suffix: {suffix}\n\n")

        check_dir(outdir, log)
        for dirname in ("RawData", "RawBeans", "Metadata", "FragPipe"):
            check_dir(os.path.join(outdir, dirname), log, required=False)

        fragpipe = os.path.join(outdir, "FragPipe")
        if os.path.isdir(fragpipe):
            fp_zip = glob.glob(os.path.join(fragpipe, f"fragpipe{suffix}.zip"))
            tables = []
            for pat in (
                "combined_protein.tsv",
                "combined_peptide.tsv",
                "abundance_protein*.tsv",
                "msstats.csv",
            ):
                tables.extend(glob.glob(os.path.join(fragpipe, pat)))
            if not fp_zip and not tables:
                fail("FragPipe/ exists but no fragpipe zip or combined/tmt/msstats tables", log)
            for path in fp_zip + tables:
                check_nonempty_file(path, log)

        # pmultiqc / MultiQC if search ran
        for qc_dir in ("pmultiqc", "MultiQC"):
            d = os.path.join(outdir, qc_dir)
            if os.path.isdir(d):
                check_nonempty_glob(os.path.join(d, f"multiqc*{suffix}*"), log, required=False)

        # MSstats / MSstatsTMT
        comparison_files = []
        for dirname, pat in (
            ("MSstats", f"msstats_comparison{suffix}.csv"),
            ("MSstatsTMT", f"msstatstmt_comparison{suffix}.csv"),
        ):
            d = os.path.join(outdir, dirname)
            if os.path.isdir(d):
                comparison_files.extend(check_nonempty_glob(os.path.join(d, pat), log, required=True))

        # FPAR
        fpar = os.path.join(outdir, "FragPipeAnalystR")
        de_files = []
        if os.path.isdir(fpar):
            de_files = check_nonempty_glob(
                os.path.join(fpar, "**", f"DE_results_*{suffix}.csv"),
                log,
                required=True,
            )
            if not de_files:
                fail("FragPipeAnalystR/ exists but no DE_results_*.csv", log)

        if check_decoys:
            check_decoy_free(comparison_files + de_files, log, decoy_prefix=args.decoy_prefix)

        if require_pp:
            for dirname in ("GeneLab", "processing_info"):
                check_dir(os.path.join(outdir, dirname), log, required=True)

            processing_info_dir = os.path.join(outdir, "processing_info")
            check_nonempty_file(os.path.join(processing_info_dir, "samples.txt"), log)
            logs = check_nonempty_glob(
                os.path.join(processing_info_dir, f"nextflow_log*{suffix}.txt"),
                log,
                required=False,
            )
            if not logs:
                check_nonempty_glob(
                    os.path.join(processing_info_dir, f"nextflow_processing_info*{suffix}.txt"),
                    log,
                    required=True,
                )
            check_nonempty_glob(
                os.path.join(processing_info_dir, f"nextflow_run_command*{suffix}.txt"), log
            )

            check_nonempty_file(os.path.join(outdir, "GeneLab", f"raw_md5sum{suffix}.tsv"), log)
            check_nonempty_file(os.path.join(outdir, "GeneLab", f"processed_md5sum{suffix}.tsv"), log)
            check_nonempty_file(
                os.path.join(outdir, "GeneLab", f"processing_info{suffix}.zip"), log
            )

        log.write("\nValidation completed successfully.\n")

    print(f"Validation log written: {output}")


if __name__ == "__main__":
    main()
