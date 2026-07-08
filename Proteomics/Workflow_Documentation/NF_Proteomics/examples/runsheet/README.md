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
| Factor Value[<name, e.g. Spaceflight>] | string | A set of one or more columns specifying the experimental group the sample belongs to. In the simplest form, a column named 'Factor Value[group]' is sufficient. Used to create the Experiment field in the FragPipe manifest. | Basal Control |

### Optional columns

| Column Name | Type | Description | Example |
|:------------|:-----|:------------|:--------|
| Source Name | string | Identifier linking samples from the same biological subject. Used for handling technical replicates during processing. Multiple samples with the same Source Name may be collapsed during analysis depending on the Has Tech Reps setting. | RR3_BSL_B7 |
| Has Tech Reps | bool | True: collapse with other True rows sharing the same Source Name and Factor Value columns (first row in the sheet wins). False or blank: keep this row. Omit the column to keep every row. | True |
| Bioreplicate | string | Numeric biological replicate ID for the FragPipe manifest. From runsheet value if set; else one ID per Source Name within each condition when Has Tech Reps column is present; else assigned in sheet order within each condition. See [FragPipe tutorial](https://fragpipe.nesvilab.org/docs/tutorial_fragpipe.html). | 1 |

---

## Data sheet

### Required columns

| Column Name | Type | Description | Example |
|:------------|:-----|:------------|:--------|
| run | string | Unique identifier for each mzML file (MS run). | NASA_Flies_TMTA_Fr00 |
| plex | string | Plex / experiment identifier (e.g. TMTa, TMTb). | TMTa |
| TechRepMixture | string | Mixture technical replicate identifier for a plex. Use `1` when the mixture was run once. Sets FragPipe manifest `Bioreplicate` and MSstatsTMT `TechRepMixture`; staging folder `{plex}_{TechRepMixture}`. | 1 |
| fraction | string | LC fraction identifier within a plex (MSstatsTMT `Fraction`). Use `1` when samples were not fractionated. | 1 |
| data_file | string | Path to mzML file. | /path/to/NASA_Flies_TMTA_Fr00.mzML |
| data_type | string | Mass spectrometry acquisition method. Written to FragPipe manifest. (Options: DDA) | DDA |

### Optional columns (data sheet)

| Column Name | Type | Description | Example |
|:------------|:-----|:------------|:--------|
| Has Tech Reps | bool | True: collapse with other True rows sharing the same plex, TechRepMixture, and fraction (first row wins). False or blank: keep this run. Omit the column to keep every row. | TRUE |
| Source Name | string | Biological subject identifier. | RR-10_BL-01 |

## Sample sheet

### Required columns

| Column Name | Type | Description | Example |
|:------------|:-----|:------------|:--------|
| Sample Name | string | Sample Name, added as a prefix to sample-specific processed data output files. Should not include spaces or weird characters. | SFug_M1 |
| organism | string | Species name used to map to the appropriate gene annotations file. Supported species can be found in the `species` column of the [GL-DPPD-7110-A_annotations.csv](https://github.com/nasa/GeneLab_Data_Processing/blob/master/GeneLab_Reference_Annotations/Pipeline_GL-DPPD-7110_Versions/GL-DPPD-7110-A/GL-DPPD-7110-A_annotations.csv) file. | Mus musculus |
| plex | string | Plex identifier. Must match data sheet. | TMTa |
| channel | string | TMT channel. | 127N |
| Factor Value[...] | string | A set of one or more columns specifying the experimental group the sample belongs to. In the simplest form, a column named 'Factor Value[group]' is sufficient. | male |


### Optional columns

| Column Name | Type | Description | Example |
|:------------|:-----|:------------|:--------|
| Source Name | string | Identifier linking TMT channels from the same biological subject. | Spaceflight microgravity Male 1 |
| Bioreplicate | string | Biological replicate ID for MSstatsTMT `BioReplicate` and FragPipeAnalystR `replicate`. From sample sheet value if set; else by `Source Name` if present; else from `Sample Name`. | 1 |
