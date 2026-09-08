# Runsheet Specification

## Description

The runsheet is a CSV file that contains the metadata required for processing mass spectrometry-based proteomics datasets through GeneLab's Proteomics processing pipeline.

- **LFQ-MBR workflow** uses a single runsheet: one row per mzML file, with sample and experiment metadata.
- **TMT workflows** use two CSV files: a **data sheet** with file paths and run identifiers; and a **sample sheet** with channel-to-sample mapping and experiment metadata.

## Examples

1. **LFQ-MBR**:
   - [Runsheet for OSD-209](OSD-209_proteomics_LFQ_v1_runsheet.csv)
2. **TMT**:
   - [Data sheet for OSD-514](OSD-514_proteomics_TMT_v1_data_sheet.csv)
   - [Sample sheet for OSD-514](OSD-514_proteomics_TMT_v1_sample_sheet.csv)

## Runsheet

### Required columns

| Column Name | Type | Description | Example |
|:------------|:-----|:------------|:--------|
| Sample Name | string | Sample Name, added as a prefix to sample-specific processed data output files. Should not include spaces or weird characters. | RR10_KDN_WT_BSL_B1 |
| organism | string | Species name used to map to the appropriate gene annotations file. Supported species can be found in the `species` column of the [GL-DPPD-7110-A_annotations.csv](https://github.com/nasa/GeneLab_Data_Processing/blob/master/GeneLab_Reference_Annotations/Pipeline_GL-DPPD-7110_Versions/GL-DPPD-7110-A/GL-DPPD-7110-A_annotations.csv) file. | Mus musculus |
| data_file | string (url or local path) | Location of the mass spectrometry data file in mzML format. | /path/to/plex1/RR10_KDN_WT_BSL_B1.mzML |
| data_type | string | Mass spectrometry acquisition method. Written to FragPipe manifest. (Options: DDA) | DDA |
| Factor Value[<name, e.g. Spaceflight>] | string | A set of one or more columns specifying the experimental group the sample belongs to. In the simplest form, a column named 'Factor Value[group]' is sufficient. | Basal Control |

### Optional columns

| Column Name | Type | Description | Example |
|:------------|:-----|:------------|:--------|
| Original Sample Name | string | Used to map the sample name that will be used for processing to the original sample name. This is often identical except in cases where the original name includes spaces or weird characters. | 150615ZV_P01 |
| Source Name | string | Identifier linking samples. Used for handling technical replicates during processing. Multiple samples with the same Source Name may be collapsed during analysis depending on the Has Tech Reps setting. | RR3_BSL_B7 |
| Has Tech Reps | bool | Indicates whether this sample is a technical replicate that should be collapsed with other samples sharing the same Source Name and Factor Value columns. Set to True for technical replicates that should be collapsed, False for distinct samples that should remain separate even if they share a Source Name. When not provided, treated as False. | True |
| Bioreplicate | string | Numeric biological replicate ID for the FragPipe manifest. Required when `require_bioreplicate` is true (default). If that flag is false: from this column if set; else one ID per Source Name within each condition when Has Tech Reps is present; else sheet order. See [FragPipe tutorial](https://fragpipe.nesvilab.org/docs/tutorial_fragpipe.html). | 1 |
| fraction | string | Optional LC fraction ID. Only needed when the dataset has both fractionated runs and technical replicates. | 1 |

Note: If the same Sample Name appears on more than one row, staging/manifest file names use each `data_file` basename instead of `Sample Name`.

---

## Data sheet

### Required columns

| Column Name | Type | Description | Example |
|:------------|:-----|:------------|:--------|
| run | string | Unique identifier for each mzML file (MS run). | NASA_Flies_TMTA_Fr00 |
| plex | string | Plex / experiment identifier (e.g. TMTa, TMTb). | TMTa |
| TechRepMixture | string | Mixture technical replicate identifier for a plex. Sets FragPipe manifest `Bioreplicate` and MSstatsTMT `TechRepMixture`. | 1 |
| fraction | string | LC fraction identifier within a plex (MSstatsTMT `Fraction`). | 1 |
| data_file | string | Path to mzML file. | /path/to/NASA_Flies_TMTA_Fr00.mzML |
| data_type | string | Mass spectrometry acquisition method. Written to FragPipe manifest. (Options: DDA) | DDA |

### Optional columns (data sheet)

| Column Name | Type | Description | Example |
|:------------|:-----|:------------|:--------|
| Has Tech Reps | bool | Indicates whether this run is a technical replicate that should be collapsed with other runs sharing the same plex, TechRepMixture, and fraction. Set to True for technical replicates that should be collapsed, False for distinct runs that should remain separate. | TRUE |

## Sample sheet

### Required columns

| Column Name | Type | Description | Example |
|:------------|:-----|:------------|:--------|
| Sample Name | string | Sample Name, added as a prefix to sample-specific processed data output files. Should not include spaces or weird characters. | SFug_M1 |
| organism | string | Species name used to map to the appropriate gene annotations file. Supported species can be found in the `species` column of the [GL-DPPD-7110-A_annotations.csv](https://github.com/nasa/GeneLab_Data_Processing/blob/master/GeneLab_Reference_Annotations/Pipeline_GL-DPPD-7110_Versions/GL-DPPD-7110-A/GL-DPPD-7110-A_annotations.csv) file. | Mus musculus |
| plex | string | Plex identifier. Must match data sheet. | TMTa |
| channel | string | TMT channel. | 127N |
| Factor Value[<name, e.g. Spaceflight>] | string | A set of one or more columns specifying the experimental group the sample belongs to. In the simplest form, a column named 'Factor Value[group]' is sufficient. | male |
| Original Sample Name | string | Used to map the sample name that will be used for processing to the original sample name. This is often identical except in cases where the original name includes spaces or weird characters. | Earth M1 |

### Required for production (biological independence)

TMT sample sheets must set **`Bioreplicate`** so MSstatsTMT / FPAR do not invent independence from `Sample Name` or row order (`--require_bioreplicate false` for smoke). Same `Bioreplicate` on different `Source Name` values in one condition is an error. ISA `Parameter Value[Biological Replicate]` maps here when curated.

### Optional columns

| Column Name | Type | Description | Example |
|:------------|:-----|:------------|:--------|
| Source Name | string | Identifier linking samples. Used to assign BioReplicate across channels and plexes when Bioreplicate is not set. | Spaceflight microgravity Male 1 |
| Bioreplicate | string | Biological replicate ID for MSstatsTMT `BioReplicate` and FragPipeAnalystR `replicate`. Required when `require_bioreplicate` is true (default). If that flag is false: from this column if set; else `Source Name`; else `Sample Name`. | 1 |

Pool / bridge channel `Sample Name` values should include the tag set by `tmtintegrator.ref_d_tag` in the FragPipe TMT-Integrator config (default: `Pool`). Only needed when a pool channel is present in the experiment.
