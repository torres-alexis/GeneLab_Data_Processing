# Runsheet Specification

## Description

The runsheet is a CSV file that contains the metadata required for processing mass spectrometry-based proteomics datasets through GeneLab's Proteomics processing pipeline.

- **LFQ-MBR workflow** uses a single runsheet: one row per mzML file, with sample and experiment metadata.
- **TMT workflows** use two CSV files: a **data sheet** with file paths and run identifiers; and a **sample sheet** with channel-to-sample mapping and experiment metadata.


## Examples

1. [Runsheet for OSD-581](OSD-209_proteomics_LFQ_v1_runsheet.csv)
2. **TMT**:
   - [Data sheet for OSD-514](OSD-514_proteomics_TMT_v1_data_sheet.csv)
   - [Sample sheet for OSD-514](OSD-514_proteomics_TMT_v1_sample_sheet.csv)

## Runsheet

### Required columns

| Column Name | Type | Description | Example |
|:------------|:-----|:------------|:--------|
| Sample Name | string | Sample Name, added as a prefix to sample-specific processed data output files. Should not include spaces or weird characters. | RR10_KDN_WT_BSL_B1 |
| organism | string | Species name used to map to the appropriate gene annotations file. Supported species can be found in the `species` column of the [GL-DPPD-7110-A_annotations.csv](https://github.com/nasa/GeneLab_Data_Processing/blob/GL_RefAnnotTable-A_1.1.0/GeneLab_Reference_Annotations/Pipeline_GL-DPPD-7110_Versions/GL-DPPD-7110-A/GL-DPPD-7110-A_annotations.csv) file. | Mus musculus |
| data_file | string (url or local path) | Location of the mass spectrometry data file in mzML format. | /path/to/plex1/RR10_KDN_WT_BSL_B1.mzML |
| Bioreplicate | string | Alphanumeric biological replicate identifier. If omitted, assigned sequentially from runsheet row order (1, 2, 3...). | 1 |
| Factor Value[<name, e.g. Spaceflight>] | string | A set of one or more columns specifying the experimental group the sample belongs to. In the simplest form, a column named 'Factor Value[group]' is sufficient. Used to create the Experiment field in the FragPipe manifest. | Basal Control |
<!--| data_type | string | Mass spectrometry acquisition method. Options: DDA | DDA | -->
<!-- | data_type | string | Mass spectrometry acquisition method. Options: DDA, DIA, GPF-DIA, DIA-Quant, DIA-Lib. | DDA | -->

### Optional columns

| Column Name | Type | Description | Example |
|:------------|:-----|:------------|:--------|
| Source Name | string | Identifier linking samples. | RR-10_BL-01 |

---

## Data sheet

### Required columns

| Column Name | Type | Description | Example |
|:------------|:-----|:------------|:--------|
| run | string | Unique identifier for each mzML file (MS run). | NASA_Flies_TMTA_Fr00 |
| plex | string | Plex identifier (e.g. TMTa, TMTb). | TMTa |
| TechRepMixture | string | Technical replicate of same mixture. Also maps to FragPipe manifest Bioreplicate. (Default: 1) | 1 |
| data_file | string | Path to mzML file. | /path/to/NASA_Flies_TMTA_Fr00.mzML |
<!--| data_type | string | Mass spectrometry acquisition method. Options: DDA | DDA | -->
<!-- | data_type | string | Mass spectrometry acquisition method. Options: DDA, DIA, GPF-DIA, DIA-Quant, DIA-Lib. | DDA | -->

## Sample sheet

### Required columns

| Column Name | Type | Description | Example |
|:------------|:-----|:------------|:--------|
| Sample Name | string | Sample Name, added as a prefix to sample-specific processed data output files. Should not include spaces or weird characters. | SFug_M1 |
| organism | string | Species name used to map to the appropriate gene annotations file. Supported species can be found in the `species` column of the [GL-DPPD-7110-A_annotations.csv](https://github.com/nasa/GeneLab_Data_Processing/blob/GL_RefAnnotTable-A_1.1.0/GeneLab_Reference_Annotations/Pipeline_GL-DPPD-7110_Versions/GL-DPPD-7110-A/GL-DPPD-7110-A_annotations.csv) file. | Mus musculus |
| plex | string | Plex identifier. Must match data sheet. | TMTa |
| channel | string | TMT channel. | 127N |
| Bioreplicate | string | Alphanumeric biological replicate identifier. If omitted, assigned sequentially from runsheet row order (1, 2, 3...). | 1 |
| Factor Value[...] | string | A set of one or more columns specifying the experimental group the sample belongs to. In the simplest form, a column named 'Factor Value[group]' is sufficient. | male |


### Optional columns

| Column Name | Type | Description | Example |
|:------------|:-----|:------------|:--------|
| Source Name | string | Identifier linking samples. | Spaceflight microgravity Male 1 |
