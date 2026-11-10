# Opt-in pick-best technical-replicate selection

This helper is a **reviewable, opt-in two-pass procedure**, not a new workflow
policy. `--tech_rep first` remains the Nextflow default and `all` remains the
other supported workflow value. Promoting pick-best to a default still requires
the AWG decision. No automatic second FragPipe run is added here.

## What is ranked

The proposed metric is the number of **distinct unmodified peptide sequences
identified in each run**, after the existing FragPipe filtering. Multiple spectra,
charges, modified forms, and repeated rows for one sequence do not inflate its
count. This is a peptide-yield metric, not a comprehensive best-QC score. The
choice of metric, its suitability for each assay, and default-policy adoption
still need scientific review.

Only rows explicitly marked `Has Tech Reps=TRUE` can compete. Groups preserve
existing metadata boundaries:

- LFQ: `Source Name` + `Factor Value[...]` columns + optional `fraction`.
- TMT: `plex` + `TechRepMixture` + `fraction`, all required on technical-run rows.

Unflagged rows are always retained, even when identifiers repeat. Single-row
groups are retained without ranking. Highest count wins; a tie uses original
sheet order. An observed run with zero eligible peptides is a valid zero. A run
absent from the evidence is **missing**, not zero, and prevents selection for its
group. The log distinguishes ties, all-zero groups, and unranked rows.

## Why the previous combined-table approach was insufficient

FragPipe `combined_peptide.tsv` columns can summarize an entire
experiment/bioreplicate folder, for example `TMTa_1`. That aggregate cannot
separate technical runs or fractions within the folder. `combined_protein.tsv`
rows count proteins, not peptides, and are now rejected.

Prefer filtered `psm.tsv` output when several technical runs share a FragPipe
folder. The helper extracts the run from `Spectrum` using the FragPipe pattern
`<run>.<scan>.<scan>.<charge>`, and uses `Peptide` for sequence identity. It excludes
rows explicitly marked `Is Decoy` or `Is Contaminant`. It does not apply a new FDR
threshold: supply **final filtered FragPipe outputs from the same first-pass
search and filtering settings** for all candidates, not raw search hits. The
helper cannot independently certify how those files were produced.

## Two-pass example

Run commands from `workflow_code/`. First perform the otherwise configured
workflow with `--tech_rep all`, preserving the original sheet and the filtered
FragPipe output. Then select technical runs explicitly:

```bash
python3 bin/select_tech_reps_by_peptides.py \
  --mode tmt \
  --input original_data_sheet.csv \
  --psm first_pass/FragPipe/*/psm.tsv \
  --output review/selected_data_sheet.csv \
  --log review/tech_reps_selected.tsv
```

For LFQ use `--mode lfq` and the original LFQ runsheet. PSM files can contain
multiple runs, and multiple files are combined using sequence sets. The glob
must include evidence for every competing run. Do not mix independent searches
or omit low-yield runs.

Read the selected sheet and audit log before running the selected inputs in a
**new output directory** with `--tech_rep first` and all other scientific settings
unchanged. For TMT, replace only the data sheet and preserve the corresponding
sample sheet. Keep the first-pass files and original sheet for provenance.

A combined peptide table is still accepted when each candidate has one unique,
run-level column:

```bash
python3 bin/select_tech_reps_by_peptides.py \
  --mode lfq \
  --input original_runsheet.csv \
  --quant first_pass/combined_peptide.tsv \
  --output review/selected_runsheet.csv \
  --log review/tech_reps_selected.tsv
```

Use a decoy/contaminant-cleaned combined table. Recognizable `rev_`/`contam_`
primary protein identifiers are also excluded as a safeguard, using the existing
workflow helper. Positive spectral counts are
preferred. Intensity-only tables are rejected unless
`--allow-intensity-fallback` is explicitly supplied. That fallback counts
**detected** sequences, may include match-between-runs (MBR) detections, and is
logged as a different metric; it is not an identified-peptide ranking.

## Mapping and audit safeguards

Evidence labels may match `Sample Name`, `run`, or the `data_file` basename
without `.mzML`/`.raw` (extension matching is case-insensitive). `Original Sample
Name` is not a run identifier and is never used. An alias shared by multiple
sheet rows is not accepted as run-specific evidence. If several distinct
aliases for one row match, the helper fails rather than guessing. Align run
labels with the actual staged filenames before trying again. Do not rename an
aggregated column to make it look run-specific.

Missing or ambiguous evidence for any competing row aborts before writing either
output. LFQ fraction protection is reused from the keep-first helper; invalid
technical-replicate flags and incomplete TMT group identifiers also fail. Output
paths must differ from inputs and each other.

The TSV preserves original sheet order and records every row's action, group,
peptide count, matched evidence label, metric, reason, original row number and
kept row number, evidence type and resolved source file paths. Row numbers are one-based data-row positions, excluding the CSV
header. Original columns and selected-row order are preserved in the output CSV.

## Verification and remaining work

Run the standard-library regression suite without FragPipe, Nextflow, or data
downloads:

```bash
PYTHONDONTWRITEBYTECODE=1 python3 -m unittest discover \
  -s tests -p 'test_select_tech_reps_by_peptides.py' -v
```

Tests cover LFQ/TMT grouping, fraction protection, unflagged rows, ties, zeros,
missing and ambiguous labels, peptide counting, PSM flags, malformed input,
explicit intensity fallback, input preservation, audit output and unchanged
workflow default. These are helper-level regressions, not end-to-end validation
of the proteomics pipeline.

Before adopting this as production policy, review the metric with AWG and run a
controlled two-pass comparison for representative LFQ, TMT and fractionated
studies. Check downstream annotations, peptide/protein yield and differential
abundance outputs. No production data is modified by developing or testing this
helper.
