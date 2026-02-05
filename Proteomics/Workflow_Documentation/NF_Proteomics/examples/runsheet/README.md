# Runsheet Specification

## Description

* The Runsheet is a csv file that contains the metadata required for processing mass spectrometry-based proteomics datasets through GeneLab's Proteomics processing pipeline.


## Examples

1. [LFQ-MBR runsheet for OSD-581](LFQ-MBR_runsheet/OSD-581_LFQ-MBR_v1_runsheet.csv)
2. [TMT10 runsheet for OSD-514](TMT10_runsheet/OSD-514_TMT10_v1_runsheet.csv)
3. [TMT16 runsheet for OSD-462](TMT16_runsheet/OSD-462_TMT16_v1_runsheet.csv)
4. [TMT16-phospho runsheet for OSD-462](TMT16_phospho_runsheet/OSD-462_TMT16-phospho_v1_runsheet.csv)

## Required columns

| Column Name | Type | Description | Example |
|:------------|:-----|:------------|:--------|
| Sample Name | string | Sample Name, added as a prefix to sample-specific processed data output files. Should not include spaces or weird characters. | RR10_KDN_WT_BSL_B1 |
| data_file | string (url or local path) | Location of the mass spectrometry data file in mzML format. | /path/to/plex1/RR10_KDN_WT_BSL_B1.mzML |
| data_type | string | Mass spectrometry data type. Options: DDA, DIA, GPF-DIA, DIA-Quant, DIA-Lib. | DDA |
| organism | string | Species name used to map to the appropriate UniProt proteome ID if a reference proteome or UniProt ID are not provided as a Nextflow parameter. Also used to extend output tables with annotations corresponding to the species. | Mus musculus |
| Source Name | string | Biological replicate identifier. Multiple fractions from the same biological sample should share the same Source Name. | RR-10_BL-01 |
| Factor Value[<name, e.g. Spaceflight>] | string | A set of one or more columns specifying the experimental group the sample belongs to. Used to create the Experiment field in the FragPipe manifest. In the simplest form, a column named 'Factor Value[group]' is sufficient. | Basal Control |

## TMT-specific columns (required for TMT workflows)

| Column Name | Type | Description | Example |
|:------------|:-----|:------------|:--------|
| plex | string | Alphanumeric identifier indicating which TMT multiplexed run the sample belongs to. | TMT1 |
| label | string | TMT channel label assigned to the sample (e.g., 126, 127N, 127C, 128N, 128C, 129N, 129C, 130N, 130C, 131N for TMT-10; 126, 127N, 127C, 128N, 128C, 129N, 129C, 130N, 130C, 131N, 131C, 132N, 132C, 133N, 133C, 134N for TMT-16). | 126 |
