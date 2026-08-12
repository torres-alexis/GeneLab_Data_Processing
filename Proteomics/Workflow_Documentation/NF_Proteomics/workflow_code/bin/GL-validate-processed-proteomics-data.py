#!/usr/bin/env python

"""
Validate GeneLab Proteomics processed datasets.
"""

import argparse
import glob
import os
import sys


def fail(message, log):
    log.write(f"FAIL: {message}\n")
    sys.stderr.write(f"Validation failed: {message}\n")
    sys.exit(1)


def check_dir(path, log):
    if not os.path.isdir(path):
        fail(f"missing directory: {path}", log)
    log.write(f"OK: directory found: {path}\n")


def check_nonempty_file(path, log):
    if not os.path.isfile(path):
        fail(f"missing file: {path}", log)
    if os.path.getsize(path) == 0:
        fail(f"empty file: {path}", log)
    log.write(f"OK: file found: {path}\n")


def check_nonempty_glob(pattern, log):
    matches = sorted(glob.glob(pattern))
    if not matches:
        fail(f"missing file matching: {pattern}", log)
    for path in matches:
        check_nonempty_file(path, log)


def main():
    parser = argparse.ArgumentParser(description="Validate NF_Proteomics processed outputs.")
    parser.add_argument("--outdir", required=True, help="Completed workflow output root")
    parser.add_argument("--assay_suffix", default="", help="e.g. _GLProteomics")
    parser.add_argument("--output", default=None, help="Validation log filename")
    args = parser.parse_args()

    outdir = os.path.abspath(args.outdir)
    output = args.output or f"validate_processed_proteomics{args.assay_suffix}.log"

    with open(output, "w") as log:
        log.write("NF_Proteomics processed output validation\n")
        log.write(f"outdir: {outdir}\n")
        log.write(f"assay_suffix: {args.assay_suffix}\n\n")

        check_dir(outdir, log)
        for dirname in ("RawData", "RawBeans", "Metadata", "FragPipe", "pmultiqc", "GeneLab", "processing_info"):
            check_dir(os.path.join(outdir, dirname), log)

        processing_info_dir = os.path.join(outdir, "processing_info")
        check_nonempty_file(os.path.join(processing_info_dir, "samples.txt"), log)
        check_nonempty_glob(os.path.join(processing_info_dir, f"nextflow_processing_info*{args.assay_suffix}.txt"), log)
        check_nonempty_glob(os.path.join(processing_info_dir, f"nextflow_run_command*{args.assay_suffix}.txt"), log)

        check_nonempty_file(os.path.join(outdir, "GeneLab", f"raw_md5sum{args.assay_suffix}.tsv"), log)
        check_nonempty_file(os.path.join(outdir, "GeneLab", f"processed_md5sum{args.assay_suffix}.tsv"), log)
        check_nonempty_file(os.path.join(outdir, "GeneLab", f"processing_info{args.assay_suffix}.zip"), log)

        log.write("\nValidation completed successfully.\n")

    print(f"Validation log written: {output}")


if __name__ == "__main__":
    main()
