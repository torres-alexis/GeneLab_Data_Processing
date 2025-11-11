# Runsheet Specification

## Description

* The Runsheet is a csv file that contains the metadata required for processing bulk RNA sequence datasets through GeneLab's RNAseq consensus processing pipeline (RCP).


## Examples

1. [Runsheet for GLDS-574](GLDS-574_proteomics_v1_runsheet.csv) 



## Required columns

| Column Name | Type | Description | Example |
|:------------|:-----|:------------|:--------|
| Sample Name | string | Sample Name, added as a prefix to sample-specific processed data output files. Should not include spaces or weird characters. | Mmus_BAL-TAL_LRTN_BSL_Rep1_B7 |
| Source Name | string | Biological replicate identifier. Used by MSstats to set BioReplicate for statistical analysis. | Mouse_1 |
| input_file | string (url or local path) | Location of the mzML file. | /my/data/sample_1.mzML |
| Factor Value[<name, e.g. Spaceflight>] | string | A set of one or more columns specifying the experimental group the sample belongs to. In the simplest form, a column named 'Factor Value[group]' is sufficient. | Space Flight |

