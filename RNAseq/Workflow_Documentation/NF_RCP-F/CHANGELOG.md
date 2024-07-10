# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.0.5](https://github.com/nasa/GeneLab_Data_Processing/tree/NF_RCP-F_1.0.5/RNAseq/Workflow_Documentation/NF_RCP-F) - 2024-07-08

### Added

- Added organism support for bacillus subtilis 
- Added the '--osdAccession OSD-###' parameter, which is required for processing OSDR datasets and optional for workflows using a runsheet.
- Added the '--isaArchivePath "/path/to/*-ISA.zip*" parameter to run the workflow using an OSD dataset ISA archive
- Added parallel processing to Perform_DGE.rmd execution using BiocParallel
- Added '--technicalReplicates "/path/to/techReps.csv" for collapsing technical replicates with DESeq2::collapseReplicates
- STAR alignment now outputs unmapped and partially mapped (discordant) reads into separate FASTQ files (Unmapped.out.mate1/2). The reads are also still present in the output BAMs.
- Added dp_tools software reporting to VV_CONCAT_FILTER

### Fixed

- Fixed deprecated syntax issues and related warnings in VV-related code.
- Fixed issue where input runsheet would not be copied to gldsAccession/Metadata/ by adding MOVE_RUNSHEET
- Output folder is now the gldsAccession and not the OSD number.
- UPDATE ISA TABLE now converts ", " to ","
- Fixed issue where RSeQC Inner Distance plots would be cutoff for datasets with reads of length 151+

### Changed

- Updated software aside from ERCC notebook, Singularity.
  - Nextflow 23.10.1 -> 24.04.2 
  - FastQC 0.11.9 -> 0.12.1
  - MultiQC 1.12 -> 1.22.3
  - Cutadapt 3.7 -> 4.9
  - TrimGalore! 0.6.7 -> 0.6.10
  - STAR 2.7.10a -> 2.7.11b
  - RSEM 1.3.1 -> 1.3.3
  - Samtools 1.15 -> 1.20
  - MultiQC 1.12 -> 1.22.3
  - gtfToGenePred,genePredToBed 377 = 377
  - infer_experiment,geneBody_coverage,inner_distance,read_distribution (rseqc) 4.0.0 -> 5.0.3
  - R 4.1.3 -> 4.4.0
  - Bioconductor 3.14.0 -> 3.19
  - DESeq2 1.34 -> 1.44.0
  - tximport 1.27.1 -> 1.32.0
  - tidyverse 1.3.1 -> 2.0.0
  - stringr 1.4.1 -> 1.5.1
  - Removed / Commented out ERCC-normalization DGE steps

## [1.0.4](https://github.com/nasa/GeneLab_Data_Processing/tree/NF_RCP-F_1.0.4/RNAseq/Workflow_Documentation/NF_RCP-F) - 2024-02-08

### Fixed

- Workflow usage files will all follow output directory set by workflow user
- ERCC Notebook:
  - Moved gene prefix definition to start of notebook
  - Added fallback for scenarios where every gene has zeros: use "poscounts" estimator to calculate a modified geometric mean
  - Reordered box-whisker plots from descending to ascending reference concentration order, ordered bar plots similarly
  
### Changed

- TrimGalore! will now use autodetect for adaptor type
- V&V migrated from dp_tools version 1.1.8 to 1.3.4 including:
  - Migration of V&V protocol code to this codebase instead of dp_tools
  - Fix for sample wise checks reusing same sample
- Added '_GLbulkRNAseq' to output file names

## [1.0.3](https://github.com/nasa/GeneLab_Data_Processing/tree/NF_RCP-F_1.0.3/RNAseq/Workflow_Documentation/NF_RCP-F) - 2023-01-25

### Added

- Test coverage using [nf-test](https://github.com/askimed/nf-test) approach

### Changed

- Updated software versions (via container update)
  - tximport == 1.27.1

### Fixed

- 'ERCC Non detection causes non-silent error' #65
- 'This function is not compatible with certain updated ISA archive metadata filenaming' #56
- 'Groups can become misassigned during group statistic calculation' #55
- 'sample to filename mapping fails when sample names are prefix substrings of other sample names' #60
- Fixed Singularity specific container issue related to DESeq2 steps

## [1.0.2](https://github.com/nasa/GeneLab_Data_Processing/tree/NF_RCP-F_1.0.2/RNAseq/Workflow_Documentation/NF_RCP-F) - 2022-11-30

### Added

- Manual tool version reporting functionality for [script](https://github.com/nasa/GeneLab_Data_Processing/tree/NF_RCP-F_1.0.2/RNAseq/Workflow_Documentation/NF_RCP-F/workflow_code/bin/format_software_versions.py) that consolidates tool versions for full workflow.
  - Currently includes manual version reporting for gtfToGenePred and genePredToBed

### Fixed

- Updated Cutadapt version in workflow from 3.4 to 3.7 in accordance with pipeline [specification](https://github.com/nasa/GeneLab_Data_Processing/tree/NF_RCP-F_1.0.2/RNAseq/Pipeline_GL-DPPD-7101_Versions/GL-DPPD-7101-F.md)

## [1.0.1](https://github.com/nasa/GeneLab_Data_Processing/tree/NF_RCP-F_1.0.1/RNAseq/Workflow_Documentation/NF_RCP-F) - 2022-11-17

### Changed

- Updated to dp_tools version 1.1.8 from 1.1.7: This addresses api changes from the release of the [OSDR](https://osdr.nasa.gov/bio/)

### Removed

- Docs: Recommendation to use Nextflow Version 21.10.6 removed as newer stable releases address original issue that had merited the recommendation

## [1.0.0](https://github.com/nasa/GeneLab_Data_Processing/tree/NF_RCP-F_1.0.0/RNAseq/Workflow_Documentation/NF_RCP-F) - 2022-11-04

### Added

- First internal production ready release of the RNASeq Consensus Pipeline Nextflow Workflow