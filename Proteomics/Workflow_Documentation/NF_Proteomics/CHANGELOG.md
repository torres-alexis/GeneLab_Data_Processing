## [1.0.0](https://github.com/nasa/GeneLab_Data_Processing/tree/NF_Proteomics_1.0.0/Proteomics/Workflow_Documentation/NF_Proteomics) - 2025-XX-XX

### Added

- First release version of GeneLab Proteomics Nextflow workflow
- `DROP_DECOYS_CONTAMS`: `decoy_contam.py` drops `rev_` / `contam_` after FragPipe, before MSstats and FragPipeAnalystR. Leading protein-group token (or all tokens). Params: `drop_decoys_contams`, `fp_analyst_keep_contaminants` (FPAR-only).
- Scrub MSstats / MSstatsTMT comparison rows with Issue, non-finite log2FC / logFC, missing p, or DF≤0 so adj.pvalue cannot stay 0.
- Runsheet preflight: per-row `require_bioreplicate` (default true) and error on reused Bioreplicate across Source Names.
- `--tech_rep` `first` (default) / `all`. `tech_reps_dropped.tsv` from keep-first.
- `VV_STEP` after FragPipe / MSstats / FPAR (`--skip_vv` to bypass). Logs under `VV_Logs/VV_log_<step>_GLProteomics.log`. VALIDATE_PROCESSING still checks published products + residual decoys/cRAP.
- Processed protocol notes decoy drop, FPAR 50% filter + seed 40, dual-stats, human pathway libs, Docker plot quality, MBR compatibility.
- README: new params, interpretation box; PPP commands drop `--enrichment_direction`.

### Changed

- Default `philosopher_contaminants_prefix` is `contam_` (Philosopher `--contamprefix` on UniProt pull and custom FASTA add).
- `drop_decoys_contams` now gates MSstats as well as FPAR.
- Public `s3://` uses `aws.client.anonymous`. Docker profile no longer enables conda.
