# GeneLab bioinformatics processing pipeline for Mass Spectrometry-based Proteomics Data (TMT10 Workflow)

> **This page holds an overview and instructions for how GeneLab processes mass spectrometry-based proteomics data using the TMT10 (Tandem Mass Tag 10-plex) workflow. Exact processing commands, GL-DPPD-[STUB] version used, and processed data output files for specific datasets are provided in the [Open Science Data Repository (OSDR)](https://osdr.nasa.gov/bio/repo/).**  

---

**Date:** March [STUB], 2026  
**Revision:** A  
**Document Number:** GL-DPPD-[STUB]-A  

**Submitted by:**  
Alexis Torres (GeneLab Data Processing Team)  

**Approved by:**  
[STUB]

---

# Table of contents  

- [**Software used**](#software-used)
- [**General processing overview with example commands**](#general-processing-overview-with-example-commands)
  - [**1. Raw Data QC**](#1-raw-data-qc)
    - [1a. RawBeans QC (Samplewise)](#1a-rawbeans-qc-samplewise)
    - [1b. RawBeans QC (All Samples)](#1b-rawbeans-qc-all-samples)
  - [**2. Download Reference Proteome, Add Decoys and Contaminants to FASTA**](#2-download-reference-proteome-add-decoys-and-contaminants-to-fasta)
  - [**3. Configure Metadata**](#3-configure-metadata)
    - [3a. Create data sheet and sample sheet](#3a-create-data-sheet-and-sample-sheet)
    - [3b. Create manifest and experiment annotation](#3b-create-manifest-and-experiment-annotation)
    - [3c. Get organism-specific gene annotations table](#3c-get-organism-specific-gene-annotations-table)
    - [3d. Stage mzML by plex and create channel-sample annotation](#3d-stage-mzml-by-plex-and-create-channel-sample-annotation)
  - [**4. FragPipe Processing Pipeline**](#4-fragpipe-processing-pipeline)
    - [4a. Launch FragPipe](#4a-launch-fragpipe)
    - [4b. Check Spectral Files Centroid Status](#4b-check-spectral-files-centroid-status)
    - [4c. Initialize Workspace](#4c-initialize-workspace)
    - [4d. MSFragger Database Search](#4d-msfragger-database-search)
    - [4e. MSBooster Deep Learning Feature Addition](#4e-msbooster-deep-learning-feature-addition)
    - [4f. Percolator PSM Rescoring and Statistical Validation](#4f-percolator-psm-rescoring-and-statistical-validation)
        - [4f1. Perform Percolator PSM Rescoring and Statistical Validation](#4f1-perform-percolator-psm-rescoring-and-statistical-validation)
        - [4f2. Add Percolator Validation Information to pepXML](#4f2-add-percolator-validation-information-to-pepxml)
    - [4g. ProteinProphet Protein Inference and Statistical Validation](#4g-proteinprophet-protein-inference-and-statistical-validation)
    - [4h. Database Annotation](#4h-database-annotation)
    - [4i. Filter Results by FDR](#4i-filter-results-by-fdr)
    - [4j. Generate Reports](#4j-generate-reports)
    - [4k. IonQuant TMT Reporter Ion Extraction](#4k-ionquant-tmt-reporter-ion-extraction)
    - [4l. TMTIntegrator TMT Quantification](#4l-tmtintegrator-tmt-quantification)
  - [**5. Compile FragPipe QC Reports**](#5-compile-fragpipe-qc-reports)
  - [**6. FragPipeAnalystR Downstream Analysis**](#6-fragpipeanalystr-downstream-analysis)

---

# Software used  

|Program|Version|Relevant Links|
|:------|:------:|:-------------|
|dp_tools|1.3.8|[https://github.com/J-81/dp_tools](https://github.com/J-81/dp_tools)|
|RawBeans|1.6.4|[https://bitbucket.org/incpm/prot-qc/src/master/protqc/](https://bitbucket.org/incpm/prot-qc/src/master/protqc/)|
|Philosopher|5.1.3|[https://github.com/Nesvilab/philosopher/releases/latest](https://github.com/Nesvilab/philosopher/releases/latest)|
|FragPipe|24.0|[https://fragpipe.nesvilab.org/](https://fragpipe.nesvilab.org/)|
|MultiQC|1.32|[https://multiqc.info/](https://multiqc.info/)|
|pmultiqc|0.0.40|[https://github.com/bigbio/pmultiqc](https://github.com/bigbio/pmultiqc)|
|R|4.5.2|[https://www.r-project.org/](https://www.r-project.org/)|
|FragPipeAnalystR|1.1.1|[https://github.com/Nesvilab/FragPipeAnalystR](https://github.com/Nesvilab/FragPipeAnalystR)|


---

# General processing overview with example commands  

<img src="../Workflow_Documentation/NF_Proteomics/images/draft_pipeline.png" align="center" alt="Proteomics TMT-10 processing workflow [STUB]" />

> Exact processing commands and output files listed in **bold** below are included with each relevant mass spectrometry-based proteomics processed dataset in the [Open Science Data Repository (OSDR)](https://osdr.nasa.gov/bio/repo/). 

---

## 1. Raw Data QC  

<br>

### 1a. RawBeans QC (Samplewise)

```bash
create-qc-report.py \
  --input *.mzML \
  --output-dir . \
  --batch \
  --cores 1

cd *
zip -r ../*_GLProteomics_qc-report.zip qc-report.html resources/
cd ..
```

**Parameter Definitions:**

- `--input` – input mzML file path
- `--output-dir` – the output directory to store results
- `--batch` – creates a report for each sample
- `--cores` – number of CPU cores to use for processing

**Input Data:**

- \*.mzML (input mass spectrometry raw data in mzML format)

**Output Data:**

- qc-report.html (RawBeans QC report HTML file)
- resources/ (directory containing supporting files for the QC report HTML)
- **\*_GLProteomics_qc-report.zip** (zip archive containing RawBeans QC report HTML file and supporting files)

<br>

### 1b. RawBeans QC (All Samples)

```bash
create-qc-report.py \
  --input *.mzML \
  --output-dir . \
  --cores 1

zip -r All_GLProteomics_qc-report.zip qc-report.html resources/
```

**Parameter Definitions:**

- `--input` – multiple mzML files provided as individual paths separated by spaces
- `--output-dir` – the output directory to store results
- `--cores` – number of CPU cores to use for processing

**Input Data:**

- \*.mzML (all input mass spectrometry raw data files in mzML format)

**Output Data:**

- qc-report.html (RawBeans QC report HTML file for all samples)
- resources/ (directory containing supporting files for the QC report HTML)
- **All_GLProteomics_qc-report.zip** (zip archive containing qc-report.html and resources/ folder for all samples)

<br>

---

## 2. Download Reference Proteome, Add Decoys and Contaminants to FASTA

```bash
  philosopher database \
    --id UPXXXXXXXXX \
    --reviewed \
    --contam
```
**Parameter Definitions:**

- `--id` – UniProt proteome ID (e.g., UP000059680)
- `--reviewed` – restrict to reviewed (Swiss-Prot) proteome entries
- `--contam` – add 116 common contaminant proteins to the FASTA database (see [Philosopher Database Wiki](https://github.com/Nesvilab/philosopher/wiki/Database))

**Output Data:**

- \*-decoys-reviewed-contam-*.fas (FASTA database containing the proteome with reversed decoy sequences and common contaminants added)

<br>

---

## 3. Configure Metadata

<br>

### 3a. Create data sheet and sample sheet

> Note: The data sheet and sample sheet may be created manually by following the [data sheet](../Workflow_Documentation/NF_Proteomics/examples/runsheet/README.md#data-sheet) and [sample sheet](../Workflow_Documentation/NF_Proteomics/examples/runsheet/README.md#sample-sheet) specifications.

```bash
### Download the *ISA.zip file from the Open Science Data Repository ###

dpt-get-isa-archive \
 --accession OSD-###

### Parse the metadata from the *ISA.zip file to create data sheet and sample sheet ###

dpt-isa-to-runsheet --accession OSD-# \
  --isa-archive *ISA.zip \
  --plugin-dir dp_tools__NF_Proteomics_TMT/
```

**Parameter Definitions:**

- `--accession` – OSD accession ID or GLDS accession ID (`GLDS-#`), used to retrieve the URLs for the ISA archive and raw data hosted in OSDR
- `--isa-archive` – Specifies the *ISA.zip file for the respective OSD dataset, downloaded in the `dpt-get-isa-archive` command
- `--plugin-dir` – Directory containing the `dp_tools` plugin used to extract sample sheet and data sheet fields from ISA metadata

**Input Data:**

- No input data required other than the OSD (or GLDS) accession ID, which is used to download the respective ISA archive

**Output Data:**

- \*ISA.zip (compressed ISA directory containing Investigation, Study, and Assay (ISA) metadata files for the respective OSD dataset, used to define sample groups — the *ISA.zip file is located in the [OSDR repository](https://osdr.nasa.gov/bio/repo/) under 'Files' → 'Study Metadata Files')

- **{OSD-Accession-ID}_Proteomics_TMT_v{version}_data_sheet.csv** (table containing mass spectrometry data file information required for processing; version denotes the dp_tools schema used to specify the metadata to extract from the ISA archive)

- **{OSD-Accession-ID}_Proteomics_TMT_v{version}_sample_sheet.csv** (table containing sample metadata required for processing; version denotes the dp_tools schema used to specify the metadata to extract from the ISA archive)

<br>

### 3b. Create manifest and experiment annotation

```bash
runsheet_to_fp_metadata.py \
  --data_sheet {OSD-Accession-ID}_Proteomics_TMT_v{version}_data_sheet.csv \
  --sample_sheet {OSD-Accession-ID}_Proteomics_TMT_v{version}_sample_sheet.csv \
  --assay_suffix _GLProteomics
```

**Parameter Definitions:**

- `--data_sheet` – path to data sheet CSV (see [data sheet specification](../Workflow_Documentation/NF_Proteomics/examples/runsheet/README.md#data-sheet))
- `--sample_sheet` – path to sample sheet CSV (see [sample sheet specification](../Workflow_Documentation/NF_Proteomics/examples/runsheet/README.md#sample-sheet))
- `--assay_suffix` – assay suffix for output filenames

**Input Data:**

- {OSD-Accession-ID}_Proteomics_TMT_v{version}_data_sheet.csv (table containing mass spectrometry data file information required for processing, output from [Step 3a](#3a-create-data-sheet-and-sample-sheet), or created manually)
- {OSD-Accession-ID}_Proteomics_TMT_v{version}_sample_sheet.csv (table containing sample metadata required for processing, output from [Step 3a](#3a-create-data-sheet-and-sample-sheet), or created manually)

**Output Data:**

- **manifest_GLProteomics.tsv** (FragPipe input table; headerless columns in order:
  - Path (mzML basename from data sheet `run` (`*.mzML`))
  - Experiment (plex identifier from data sheet `plex`)
  - Bioreplicate (mixture technical replicate alphanumeric identifier (from data sheet `TechRepMixture`, or `1` if empty))
  - Data type (preset value (`DDA`)))

- **experiment_annotation_GLProteomics.tsv** (FragPipeAnalystR input table with additional `condition_name` column; columns in order:
  - plex (`{Experiment}_{Bioreplicate}` (matches manifest `Experiment` and `Bioreplicate`))
  - channel (TMT reporter channel from sample sheet `channel`)
  - sample (R-safe sample identifier from sample sheet `Sample Name`)
  - sample_name (same as `sample`)
  - condition (R-safe condition symbol from joined sample sheet `Factor Value[...]` values)
  - condition_name (human-readable condition)
  - replicate (biological replicate replicate alphanumeric identifier (from sample sheet `Bioreplicate` column if present; else `1`)))

<br>

### 3c. Get organism-specific gene annotations table

```r
### Sample sheet from Step 3a; organism must match the value in the species column of GL-DPPD-7110-A_annotations.csv ###
sample_sheet_path <- "{OSD-Accession-ID}_Proteomics_TMT_v{version}_sample_sheet.csv"
sample_sheet <- read.csv(sample_sheet_path, stringsAsFactors = FALSE, check.names = FALSE)
organism <- trimws(as.character(sample_sheet[["organism"]][1]))

### Pull in the GeneLab annotation table (GL-DPPD-7110-A_annotations.csv) ###
org_table_link <- "https://raw.githubusercontent.com/nasa/GeneLab_Data_Processing/master/GeneLab_Reference_Annotations/Pipeline_GL-DPPD-7110_Versions/GL-DPPD-7110-A/GL-DPPD-7110-A_annotations.csv"

org_table <- read.table(org_table_link, sep = ",", header = TRUE)

### URL of the organism-specific GeneLab gene annotation table ###
annotations_link <- org_table[org_table$species == organism, "genelab_annots_link"]
```

**Input Data:**

- {OSD-Accession-ID}_Proteomics_TMT_v{version}_sample_sheet.csv (output from [Step 3a](#3a-create-data-sheet-and-sample-sheet); `organism` column value must match a value in the `species` column of [GL-DPPD-7110-A_annotations.csv](../../GeneLab_Reference_Annotations/Pipeline_GL-DPPD-7110_Versions/GL-DPPD-7110-A/GL-DPPD-7110-A_annotations.csv))

**Output Data:**

- annotations_link (variable containing URL of organism-specific GeneLab gene annotation table)

<br>

### 3d. Stage mzML by plex and create channel-sample annotation

```bash
# TMT: reorganize mzML into plex-specific (manifest.tsv Experiment_Bioreplicate) folders + create annotation.txt before FragPipe
bash tmt_stage_by_plex.sh manifest_GLProteomics.tsv experiment_annotation_GLProteomics.tsv TMT10
```

**Parameter Definitions:**

- `manifest_GLProteomics.tsv` – path to FragPipe manifest table
- `experiment_annotation_GLProteomics.tsv` – path to FragPipe experiment_annotation table 
- `TMT10` – isobaric multiplex type; sets reporter channel order and row count for each annotation table created

**Input Data:**

- manifest_GLProteomics.tsv (FragPipe input table, output from [Step 3b](#3b-create-manifest-and-experiment-annotation))
- experiment_annotation_GLProteomics.tsv (FragPipeAnalystR input table with additional `condition_name` column, output from [Step 3b](#3b-create-manifest-and-experiment-annotation))
- \*.mzML (input mass spectrometry raw data in mzML format)

**Output Data:**

- {Experiment}_{Bioreplicate}/\*.mzML (input mass spectrometry raw data in mzML format organized into subdirectories by plex)
- {Experiment}_{Bioreplicate}/annotation.txt (channel-sample annotation table for each plex)
- manifest_GLProteomics.tsv (FragPipe input table with mzML file paths updated with plex directories)

<br>

---

## 4. FragPipe Processing Pipeline

<br>

### 4a. Launch FragPipe

```bash
/fragpipe_bin/fragpipe-24.0/fragpipe-24.0/bin/fragpipe \
  --headless \
  --workflow TMT10.workflow \
  --manifest manifest_GLProteomics.tsv \
  --workdir . \
  --ram 64 \
  --threads 16 \
  --config-tools-folder tools
```

**Parameter Definitions:**

- `--headless` – run FragPipe in headless mode (no GUI)
- `--workflow` – path to FragPipe workflow configuration file
- `--manifest` – path to manifest TSV file containing plex information and file paths
- `--workdir` – working directory for FragPipe execution
- `--ram` – Memory (GB) allocated to FragPipe
- `--threads` – number of CPU threads allocated to FragPipe
- `--config-tools-folder` – path to folder containing FragPipe tools not included in the Docker image (MSFragger JAR, IonQuant JAR, diaTracer JAR, ext/bruker/, ext/thermo/)

**Input Data:**

- TMT10.workflow (FragPipe workflow configuration file for TMT-10 workflow)
- manifest_GLProteomics.tsv (FragPipe input table, output from [Step 3d](#3d-stage-mzml-by-plex-and-create-channel-sample-annotation))
- \*-decoys-reviewed-contam-*.fas (proteome FASTA database, output from [Step 2](#2-download-reference-proteome-add-decoys-and-contaminants-to-fasta))
- \*.mzML (input mass spectrometry raw data in mzML format)
- \*_annotation.txt (channel-sample annotation table for each plex, output from [Step 3d](#3d-stage-mzml-by-plex-and-create-channel-sample-annotation))

**Output Data:**

- fragger.params (MSFragger parameter configuration file)
- msbooster_params.txt (MSBooster parameter configuration file)
- tmt-integrator-conf.yml (TMTIntegrator configuration file)
- filelist_proteinprophet.txt (list of interact.pep.xml files to be passed to ProteinProphet)
- filelist_ionquant.txt (file list for IonQuant)
- modmasses_ionquant.txt (modification masses file for IonQuant)
- experiment_annotation.tsv (experiment annotation file mapping TMT channels to samples)
- fragpipe.workflow (FragPipe output workflow configuration file)
- fragpipe-files.fp-manifest (FragPipe output manifest)
- fragpipe.job (FragPipe job configuration file)
- log_\*.txt (FragPipe execution log file with timestamp)
- sdrf.tsv (Sample and Data Relationship Format file)

<br>

### 4b. Check Spectral Files Centroid Status

```bash
java -Xmx64G -cp /fragpipe_bin/fragpipe-24.0/fragpipe-24.0/lib/fragpipe-24.0.jar:/fragpipe_bin/fragpipe-24.0/fragpipe-24.0/tools/batmass-io-1.35.4.jar org.nesvilab.fragpipe.util.CheckCentroid *.mzML 16
```

**Parameter Definitions:**

- `-Xmx64G` – Java memory limit (e.g., `-Xmx64G` for 64 GB RAM)
- `-cp` – Java classpath to FragPipe and BatMass libraries
- `org.nesvilab.fragpipe.util.CheckCentroid` – CheckCentroid main class
- `*.mzML` – input mzML file(s) to check
- `16` – number of CPU threads to use

**Input Data:**

- \*.mzML (input mass spectrometry raw data in mzML format)

**Output Data:**

- (No output files; checks if mzML files are centroided or profile mode; FragPipe exits if files are not centroided)

<br>

### 4c. Initialize Workspace

```bash
/fragpipe_bin/fragpipe-24.0/fragpipe-24.0/tools/Philosopher/philosopher-v5.1.3-RC9 workspace --clean --nocheck
/fragpipe_bin/fragpipe-24.0/fragpipe-24.0/tools/Philosopher/philosopher-v5.1.3-RC9 workspace --init --nocheck --temp /tmp/temp_directory
```

**Parameter Definitions:**

- `workspace` – Philosopher subcommand for managing workspace
- `--clean` – removes any existing workspace files
- `--init` – initializes a new Philosopher workspace
- `--nocheck` – skips workspace validation checks
- `--temp` – specifies temporary directory for workspace initialization

**Output Data:**

- .meta/ (Philosopher workspace metadata directory containing binary database files)

<br>

### 4d. MSFragger Database Search

```bash
java -jar -Dfile.encoding=UTF-8 -Xmx64G MSFragger-4.3.jar fragger.params plexA_1/sample1.mzML plexA_1/sample2.mzML
```

**Parameter Definitions:**

- `-jar` – executes JAR file
- `-Dfile.encoding=UTF-8` – sets file encoding to UTF-8
- `-Xmx64G` – Java memory limit (e.g., `-Xmx64G` for 64 GB RAM)
- `fragger.params` – MSFragger parameter configuration file
- `*.mzML` – multiple mzML files provided as individual paths separated by spaces

**Input Data:**

- fragger.params (MSFragger parameter configuration file, output from [Step 4a](#4a-launch-fragpipe))
- \*.mzML (input mass spectrometry raw data in mzML format)
- \*-decoys-reviewed-contam-\*.fas (proteome FASTA database with decoys and contaminants, output from [Step 2](#2-download-reference-proteome-add-decoys-and-contaminants-to-fasta))

**Output Data:**

- **\*.pepXML** (peptide-spectrum matches from the MSFragger database search)
- **\*.pin** (peptide-spectrum matches from the MSFragger database search in Percolator input format (PIN) for statistical validation)
- **\*.pepindex** (peptide index files for the FASTA database)
- **\*.tsv** (MSFragger results in tab-separated format)

<br>

### 4e. MSBooster Deep Learning Feature Addition

```bash
java -Djava.awt.headless=true -Xmx64G -cp MSBooster-1.3.17.jar:batmass-io-1.35.4.jar mainsteps.MainClass --paramsList msbooster_params.txt
```

**Parameter Definitions:**

- `-Djava.awt.headless=true` – runs in headless mode (no GUI)
- `-Xmx64G` – Java memory limit (e.g., `-Xmx64G` for 64 GB RAM)
- `-cp` – Java classpath to MSBooster and BatMass libraries
- `mainsteps.MainClass` – MSBooster main class
- `--paramsList` – path to MSBooster parameter configuration file

**Input Data:**

- msbooster_params.txt (MSBooster parameter configuration file, output from [Step 4a](#4a-launch-fragpipe))
- \*.pin (Percolator input files from MSFragger, output from [Step 4d](#4d-msfragger-database-search))
- \*.mzML (input mass spectrometry raw data in mzML format)

**Output Data:**

- \*_edited.pin (Percolator input files with added deep learning features from MSBooster: unweighted spectral entropy, weighted spectral entropy, hypergeometric probability, intersection, predicted RT real units, and delta RT LOESS)
- spectraRT_full.tsv (full spectra retention time data)
- spectraRT.predicted.bin (binary file containing predicted spectra, retention times, and ion mobilities from DIA-NN)
- spectraRT.tsv (spectra retention time data)
- MSBooster_plots/ (Directory containing MSBooster calibration and diagnostic plots)

<br>

### 4f. Percolator PSM Rescoring and Statistical Validation

<br>

#### 4f1. Perform Percolator PSM Rescoring and Statistical Validation

```bash
/fragpipe_bin/fragpipe-24.0/fragpipe-24.0/tools/percolator_3_7_1/linux/percolator \
  --only-psms \
  --no-terminate \
  --post-processing-tdc \
  --num-threads 16 \
  --results-psms *_percolator_target_psms.tsv \
  --decoy-results-psms *_percolator_decoy_psms.tsv \
  --protein-decoy-pattern rev_ \
  *_edited.pin
```

**Parameter Definitions:**

- `--only-psms` – do not remove redundant peptides, keep PSMs, exclude peptide level probabilities
- `--no-terminate` – do not terminate execution when encountering issues with SVM inputs or results
- `--post-processing-tdc` – replace mix-max method with target-decoy competition for assigning q-values and PEPs. For input PSMs from separate target/decoy searches, Percolator SVM scores eliminate lower-scoring target or decoy PSMs for each scan+expMass combination. Automatically enabled for concatenated searches
- `--num-threads` – number of CPU threads to use
- `--results-psms` – output file path for target PSM results
- `--decoy-results-psms` – output file path for decoy PSM results
- `--protein-decoy-pattern` – text pattern used to identify decoy proteins in the database
- `*_edited.pin` – input Percolator input files with MSBooster features

**Input Data:**

- \*_edited.pin (Percolator input files with MSBooster features, output from [Step 4e](#4e-msbooster-deep-learning-feature-addition))

**Output Data:**

- \*_percolator_target_psms.tsv (Percolator target PSM results in TSV format)
- \*_percolator_decoy_psms.tsv (Percolator decoy PSM results in TSV format)

<br>

#### 4f2. Add Percolator Validation Information to pepXML

```bash
java -cp /fragpipe_bin/fragpipe-24.0/fragpipe-24.0/lib/* \
  org.nesvilab.fragpipe.tools.percolator.PercolatorOutputToPepXML \
  *.pin \
  * \
  *_percolator_target_psms.tsv \
  *_percolator_decoy_psms.tsv \
  interact-* \
  DDA \
  0.5 \
  *.mzML
```

**Parameter Definitions:**

- `-cp` – Java classpath to FragPipe libraries
- `org.nesvilab.fragpipe.tools.percolator.PercolatorOutputToPepXML` – FragPipe utility class for adding Percolator validation information to pepXML files
- `*.pin` – original Percolator input PIN file
- `*` – sample name
- `*_percolator_target_psms.tsv` – Percolator target PSM results
- `*_percolator_decoy_psms.tsv` – Percolator decoy PSM results
- `interact-*` – output pepXML file prefix
- `DDA` – data acquisition type (DDA|DIA|GPF-DIA|DIA-Quant|DIA-Lib)
- `0.5` – minimum probability threshold (1 - PEP); filters out PSMs with PEP > 0.5
- `*.mzML` – input mzML file path

**Input Data:**

- \*.pin (original Percolator input files from MSFragger, output from [Step 4d](#4d-msfragger-database-search))
- \*_percolator_target_psms.tsv (Percolator target PSM results, output from [Step 4f1](#4f1-perform-percolator-psm-rescoring-and-statistical-validation))
- \*_percolator_decoy_psms.tsv (Percolator decoy PSM results, output from [Step 4f1](#4f1-perform-percolator-psm-rescoring-and-statistical-validation))
- \*.mzML (input mass spectrometry raw data in mzML format)

**Output Data:**

- interact-\*.pep.xml (peptide-spectrum matches with validation information generated by Percolator)

<br>

### 4g. ProteinProphet Protein Inference and Statistical Validation

```bash
/fragpipe_bin/fragpipe-24.0/fragpipe-24.0/tools/Philosopher/philosopher-v5.1.3-RC9 proteinprophet --maxppmdiff 2000000 --minprob 0.5 --output combined filelist_proteinprophet.txt
```

**Parameter Definitions:**

- `proteinprophet` – run ProteinProphet to generate probabilities for protein identifications based on MS/MS data
- `--maxppmdiff 2000000` – maximum peptide mass difference in ppm 
- `--minprob 0.5` – PeptideProphet minimum probability threshold 
- `--output combined` – output file name; results in `combined.prot.xml`
- `filelist_proteinprophet.txt` – list of interact.pep.xml files to be passed to ProteinProphet

**Input Data:**

- filelist_proteinprophet.txt (list of interact.pep.xml files to be passed to ProteinProphet, output from [Step 4a](#4a-launch-fragpipe))
- interact-\*.pep.xml (pepXML files listed in filelist_proteinprophet.txt, output from [Step 4f2](#4f2-add-percolator-validation-information-to-pepxml))

**Output Data:**

- combined.prot.xml (protein identifications with validation information generated by ProteinProphet via Philosopher)

<br>

### 4h. Database Annotation

```bash
/fragpipe_bin/fragpipe-24.0/fragpipe-24.0/tools/Philosopher/philosopher-v5.1.3-RC9 database --annotate *.fas --prefix rev_
```

**Parameter Definitions:**

- `database --annotate` – annotate FASTA database file (creates binary database files for Philosopher tools)
- `*.fas` – path to FASTA database file
- `--prefix rev_` – decoy prefix used in the database

**Input Data:**

- \*-decoys-reviewed-contam-\*.fas (proteome FASTA database with decoys and contaminants, output from [Step 2](#2-download-reference-proteome-add-decoys-and-contaminants-to-fasta))

**Output Data:**

- .meta/ (Philosopher workspace metadata directory containing binary database files)

<br>

### 4i. Filter Results by FDR

```bash
# First plex (initializes database annotation)
/fragpipe_bin/fragpipe-24.0/fragpipe-24.0/tools/Philosopher/philosopher-v5.1.3-RC9 filter \
  --sequential \
  --prot 0.01 \
  --picked \
  --tag rev_ \
  --pepxml plex_directory \
  --protxml combined.prot.xml \
  --razor

# Subsequent plexes (reuse database annotation from first plex)
/fragpipe_bin/fragpipe-24.0/fragpipe-24.0/tools/Philosopher/philosopher-v5.1.3-RC9 filter \
  --sequential \
  --prot 0.01 \
  --picked \
  --tag rev_ \
  --pepxml plex_directory \
  --dbbin first_plex_directory \
  --protxml combined.prot.xml \
  --probin first_plex_directory \
  --razor
```

**Parameter Definitions:**

- `filter` – filter PSMs, peptides, and proteins by FDR threshold
- `--sequential` – apply sequential FDR filtering at PSM, peptide, and ion levels in addition to protein level FDR
- `--prot 0.01` – protein-level FDR threshold
- `--picked` – apply picked FDR algorithm prior to protein scoring
- `--tag rev_` – decoy sequence prefix
- `--pepxml` – path to plex directory containing run-specific pepXML files
- `--protxml combined.prot.xml` – path to protXML file
- `--dbbin` – (for subsequent plexes) path to first plex directory containing database annotation
- `--probin` – (for subsequent plexes) path to first plex directory containing protein annotation
- `--razor` – use razor peptides for protein-level FDR scoring

**Input Data:**

- interact-\*.pep.xml (peptide-spectrum matches with validation information generated by Percolator, output from [Step 4f2](#4f2-add-percolator-validation-information-to-pepxml))
- combined.prot.xml (protein identifications with validation information generated by ProteinProphet via Philosopher, output from [Step 4g](#4g-proteinprophet-protein-inference-and-statistical-validation))
- .meta/ (Philosopher workspace metadata, output from [Step 4h](#4h-database-annotation))

**Output Data:**

- filter.log (Philosopher filter execution log file)
- Filtered data stored in Philosopher workspace as binary files (db.bin, ion.bin, pep.bin, pro.bin, protxml.bin, psm.bin, razor.bin, etc.; filtered PSM, peptide, and protein data ready for report generation)

<br>

### 4j. Generate Reports

```bash
/fragpipe_bin/fragpipe-24.0/fragpipe-24.0/tools/Philosopher/philosopher-v5.1.3-RC9 report
```

**Input Data:**

- Philosopher workspace containing filtered data (output from [Step 4i](#4i-filter-results-by-fdr))

**Output Data:**

- protein.fas (FASTA file containing FDR-filtered protein sequences identified)
- protein.tsv (plex-specific protein report)
- peptide.tsv (plex-specific peptide report)
- psm.tsv (plex-specific PSM report)
- ion.tsv (plex-specific ion report)

<br>

### 4k. IonQuant TMT Reporter Ion Extraction

```bash
# MS1 — initial
java -Djava.awt.headless=true -Xmx64G \
  -Dlibs.bruker.dir=tools/ext/bruker \
  -Dlibs.thermo.dir=tools/ext/thermo \
  -cp /fragpipe_bin/fragpipe-24.0/fragpipe-24.0/tools/jfreechart-1.5.3.jar:tools/IonQuant-1.11.20.jar \
  ionquant.IonQuant \
  --threads 16 \
  --perform-ms1quant 1 \
  --perform-isoquant 0 \
  --isotol 20.0 \
  --isolevel 2 \
  --isotype tmt10 \
  --ionmobility 0 \
  --site-reports 0 \
  --msstats 0 \
  --minexps 1 \
  --mbr 0 \
  --maxlfq 0 \
  --requantify 0 \
  --mztol 10 \
  --imtol 0.05 \
  --rttol 1 \
  --normalization 0 \
  --minisotopes 1 \
  --minscans 1 \
  --writeindex 0 \
  --tp 0 \
  --minfreq 0 \
  --minions 1 \
  --locprob 0 \
  --uniqueness 0 \
  --multidir . \
  --filelist filelist_ionquant.txt \
  --modlist modmasses_ionquant.txt

# MS1 — full parameters (second FragPipe MS1 call)
java -Djava.awt.headless=true -Xmx64G \
  -Dlibs.bruker.dir=tools/ext/bruker \
  -Dlibs.thermo.dir=tools/ext/thermo \
  -cp /fragpipe_bin/fragpipe-24.0/fragpipe-24.0/tools/jfreechart-1.5.3.jar:tools/IonQuant-1.11.20.jar \
  ionquant.IonQuant \
  --threads 16 \
  --perform-ms1quant 1 \
  --perform-isoquant 0 \
  --isotol 20.0 \
  --isolevel 2 \
  --isotype tmt10 \
  --ionmobility 0 \
  --site-reports 1 \
  --msstats 1 \
  --minexps 1 \
  --mbr 0 \
  --maxlfq 1 \
  --requantify 1 \
  --mztol 10 \
  --imtol 0.05 \
  --rttol 0.4 \
  --mbrmincorr 0 \
  --mbrrttol 1 \
  --mbrimtol 0.05 \
  --mbrtoprun 10 \
  --ionfdr 0.01 \
  --proteinfdr 1 \
  --peptidefdr 1 \
  --normalization 1 \
  --minisotopes 2 \
  --intensitymode 2 \
  --minscans 3 \
  --writeindex 0 \
  --tp 0 \
  --minfreq 0 \
  --minions 1 \
  --locprob 0.75 \
  --uniqueness 0 \
  --multidir . \
  --filelist filelist_ionquant.txt \
  --modlist modmasses_ionquant.txt

# TMT-10 reporters
java -Djava.awt.headless=true -Xmx64G \
  -Dlibs.bruker.dir=tools/ext/bruker \
  -Dlibs.thermo.dir=tools/ext/thermo \
  -cp /fragpipe_bin/fragpipe-24.0/fragpipe-24.0/tools/jfreechart-1.5.3.jar:tools/IonQuant-1.11.20.jar \
  ionquant.IonQuant \
  --threads 16 \
  --perform-ms1quant 0 \
  --perform-isoquant 1 \
  --isotol 20.0 \
  --isolevel 2 \
  --isotype TMT-10 \
  --ionmobility 0 \
  --site-reports 0 \
  --msstats 0 \
  --annotation plexA_1/psm.tsv=plexA_1/annotation.txt \
  --annotation plexB_1/psm.tsv=plexB_1/annotation.txt \
  --minexps 1 \
  --mbr 0 \
  --maxlfq 0 \
  --requantify 0 \
  --mztol 10 \
  --imtol 0.05 \
  --rttol 1 \
  --normalization 0 \
  --minisotopes 1 \
  --minscans 1 \
  --writeindex 0 \
  --tp 0 \
  --minfreq 0 \
  --minions 1 \
  --locprob 0 \
  --uniqueness 0 \
  --multidir . \
  --filelist filelist_ionquant.txt \
  --modlist modmasses_ionquant.txt
```

**Parameter Definitions:**

- `-Djava.awt.headless=true` – run in headless mode (no GUI)
- `-Xmx64G` – Java memory limit (e.g., `-Xmx64G` for 64 GB RAM)
- `-Dlibs.bruker.dir` – directory for Bruker libraries
- `-Dlibs.thermo.dir` – directory for Thermo libraries
- `-cp` – Java classpath to jfreechart and IonQuant JAR files
- `ionquant.IonQuant` – IonQuant main class
- `--threads` – number of CPU threads to use
- `--perform-ms1quant 1` – perform MS1 quantification (first two IonQuant steps for TMT10; 0 = no, 1 = yes)
- `--perform-isoquant 1` – perform isobaric labeling quantification (third IonQuant step only; 0 = no, 1 = yes)
- `--isotol 20.0` – MS2 tolerance in ppm for isobaric quantification
- `--isolevel 2` – isobaric quantification level (2 = MS2, 3 = MS3, 4 = ZOOM-HR)
- `--isotype tmt10` / `--isotype TMT-10` – isobaric quantification type (first two MS1 steps: `tmt10`; reporter step: `TMT-10`)
- `--ionmobility 0` – data has ion mobility information (0 = no, 1 = yes)
- `--site-reports 0` – generate site reports (0 = no, 1 = yes)
- `--msstats 0` – generate MSstats input files (0 = no, 1 = yes)
- `--annotation` – `{plex}_{TechRepMixture}/psm.tsv={plex}_{TechRepMixture}/annotation.txt` (repeatable)
- `--minexps 1` – minimum experiments in picking an ion for quantifying proteins (only for intensity, not MaxLFQ)
- `--mbr 0` – perform match-between-runs (0 = no, 1 = yes)
- `--maxlfq 0` – calculate MaxLFQ intensity (0 = no, 1 = yes)
- `--requantify 0` – re-quantify unidentified features based on identified features (0 = no, 1 = yes)
- `--mztol 10` – MS1 tolerance in ppm
- `--imtol 0.05` – 1/K0 tolerance
- `--rttol 1` – retention time tolerance in minutes
- `--normalization 0` – normalize intensities across all runs (0 = no, 1 = yes)
- `--minisotopes 1` – minimum isotopes required in feature extraction
- `--minscans 1` – minimum MS1 scans required in feature extraction
- `--writeindex 0` – write indexed file on disk for further usage (0 = no, 1 = yes)
- `--tp 0` – number of ions used in quantifying each protein (0 = use all ions; only for intensity, not MaxLFQ)
- `--minfreq 0` – minimum required frequency of an ion being selected for protein quantification (only for intensity, not MaxLFQ)
- `--minions 1` – minimum ions required for quantifying proteins (only for MaxLFQ intensity)
- `--locprob 0` – localization probability threshold
- `--uniqueness 0` – peptide-protein uniqueness (0 = unique+razor, 1 = unique only, 2 = all)
- `--multidir .` – output directory for multi-experimental results
- `--filelist` – file list for IonQuant listing, for each plex, a `--psm` entry pointing to the `psm.tsv` file and a `--specdir` entry pointing to the directory containing mzML files
- `--modlist` – file listing modification masses (used to remove mass discrepancy due to rounding errors)

**Input Data:**

- filelist_ionquant.txt (file list for IonQuant from [Step 4a](#4a-launch-fragpipe); listing for each plex, a `--psm` entry pointing to the `psm.tsv` file and a `--specdir` entry pointing to the directory containing mzML files)
- modmasses_ionquant.txt (modification masses file for IonQuant, output from [Step 4a](#4a-launch-fragpipe))
- \*_annotation.txt (plex-specific annotation files mapping TMT channels to sample names, output from [Step 3d](#3d-stage-mzml-by-plex-and-create-channel-sample-annotation))
- \*.mzML (input mass spectrometry raw data in mzML format; per plex, specified with `--specdir` in filelist_ionquant.txt)

**Output Data:**

- protein.tsv (plex-specific protein report with TMT reporter ion intensities and additional data added from IonQuant)
- peptide.tsv (plex-specific peptide report with TMT reporter ion intensities and additional data added from IonQuant)
- ion.tsv (plex-specific ion report with TMT reporter ion intensities and additional data added from IonQuant)
- psm.tsv (plex-specific PSM report with TMT reporter ion intensities and additional data added from IonQuant)
- **combined_protein.tsv** (combined protein report with TMT reporter ion intensities across all plexes)
- **combined_peptide.tsv** (combined peptide report with TMT reporter ion intensities and additional data across all plexes)
- **combined_ion.tsv** (combined ion report with TMT reporter ion intensities and additional data across all plexes)
- **combined_modified_peptide.tsv** (combined modified peptide report with TMT reporter ion intensities and additional data across all plexes)
- **msstats.csv** (input file for MSstats downstream differential analysis)
- **msstats_ptm.csv** (input file for MSstatsPTM PTM (post-translational modification) analysis)
- reprint.int.tsv (input file for the Resource for Evaluation of Protein Interaction Networks (REPRINT) containing protein intensities)
- reprint.spc.tsv (input file for the Resource for Evaluation of Protein Interaction Networks (REPRINT) containing protein spectral counts)

<br>

### 4l. TMTIntegrator TMT Quantification

```bash
java -Xmx64G -jar TMT-Integrator-6.1.1.jar \
  tmt-integrator-conf.yml \
  sample1/psm.tsv \
  sample2/psm.tsv
```

**Parameter Definitions:**

- `-Xmx64G` – Java memory limit (e.g., `-Xmx64G` for 64 GB RAM)
- `-jar` – executes JAR file
- `tmt-integrator-conf.yml` – TMTIntegrator configuration file
- `*/psm.tsv` – PSM files with TMT reporter ion intensities

**Input Data:**

- tmt-integrator-conf.yml (TMTIntegrator configuration file, output from [Step 4a](#4a-launch-fragpipe))
- psm.tsv (PSM reports with TMT reporter ion intensities, output from [Step 4k](#4k-ionquant-tmt-reporter-ion-extraction))

**Output Data:**

- abundance_protein_MD.tsv (TMTIntegrator protein-level table; log2 abundance per channel)
- abundance_peptide_MD.tsv (TMTIntegrator peptide-level table; log2 abundance per channel)
- abundance_gene_MD.tsv (TMTIntegrator gene-level table; log2 abundance per channel)
- ratio_protein_MD.tsv (TMTIntegrator protein-level table; log2(channel/reference) per channel)
- ratio_peptide_MD.tsv (TMTIntegrator peptide-level table; log2(channel/reference) per channel)
- ratio_gene_MD.tsv (TMTIntegrator gene-level table; log2(channel/reference) per channel)

<br>

---

## 5. Compile FragPipe QC Reports

```bash
multiqc --fragpipe-plugin \
  -o /path/to/pmultiqc/output/directory \
  -n multiqc_GLProteomics \
  /path/to/FragPipe/output/directory

clean_multiqc_paths.py multiqc_GLProteomics_data /path/to/pmultiqc/output/directory
```

**Parameter Definitions:**

- `--fragpipe-plugin` – enable FragPipe plugin for MultiQC to process FragPipe output files
- `-o` – the output directory to store results
- `-n` – prefix name for output files
- `/path/to/FragPipe/output/directory` – the directory containing FragPipe output files, provided as a positional argument
- `clean_multiqc_paths.py` – Python script to clean absolute paths (if present) from MultiQC data files and create a zip archive
- `multiqc_GLProteomics_data` – name of the MultiQC data directory to process
- `/path/to/pmultiqc/output/directory` – output directory where the zip file will be created

**Input Data:**

- psm.tsv (plex-specific PSM reports, output from [Step 4k](#4k-ionquant-tmt-reporter-ion-extraction))
- ion.tsv (plex-specific ion reports, output from [Step 4k](#4k-ionquant-tmt-reporter-ion-extraction))
- combined_protein.tsv (combined protein report, output from [Step 4k](#4k-ionquant-tmt-reporter-ion-extraction))
- combined_peptide.tsv (combined peptide report, output from [Step 4k](#4k-ionquant-tmt-reporter-ion-extraction))
- combined_ion.tsv (combined ion report, output from [Step 4k](#4k-ionquant-tmt-reporter-ion-extraction))
- fragpipe.workflow (FragPipe workflow configuration file, output from [Step 4a](#4a-launch-fragpipe))
- fragger.params (MSFragger parameters file, output from [Step 4a](#4a-launch-fragpipe))

**Output Data:**

- **multiqc_GLProteomics.html** (MultiQC output html summary)
- **multiqc_GLProteomics_data.zip** (zipped directory containing MultiQC output data with cleaned paths)

<br>

---

## 6. FragPipeAnalystR Downstream Analysis

The FragPipeAnalystR downstream analysis script is executed three times: once using the **protein**-level quantification file (abundance_protein_MD.tsv), once using the **gene**-level quantification file (abundance_gene_MD.tsv), and once using the **peptide**-level quantification file (abundance_peptide_MD.tsv).

**Protein run:**

```bash
Rscript FragPipeAnalystR_main.R \
  --experiment_annotation "experiment_annotation_GLProteomics.tsv" \
  --quantification_file "abundance_protein_MD.tsv" \
  --mode "TMT" \
  --level "protein" \
  --feature_list_protein "" \
  --feature_list_gene "" \
  --top_n_protein 10 \
  --top_n_gene 10 \
  --enrichment_database "Hallmark,GO_Biological_Process_2021" \
  --enrichment_direction "Up,Down" \
  --gsea_database "Hallmark,GO_Biological_Process_2021" \
  --normalization_method "none" \
  --de_alpha 0.05 \
  --de_lfc 1.0 \
  --de_fdr "Benjamini Hochberg" \
  --imputation_type "Perseus-type" \
  --imputation_shift 1.8 \
  --imputation_scale 0.3 \
  --qc_plot_data "nonimputed" \
  --sample_cvs_full_range "false" \
  --volcano_display_names "true" \
  --volcano_show_gene "true" \
  --gene_annotations $annotations_link \
  --output_dir "output/"
```

**Gene run:**

```bash
Rscript FragPipeAnalystR_main.R \
  --experiment_annotation "experiment_annotation_GLProteomics.tsv" \
  --quantification_file "abundance_gene_MD.tsv" \
  --mode "TMT" \
  --level "gene" \
  --feature_list_gene "" \
  --top_n_gene 10 \
  --enrichment_database "Hallmark,GO_Biological_Process_2021" \
  --enrichment_direction "Up,Down" \
  --gsea_database "Hallmark,GO_Biological_Process_2021" \
  --normalization_method "none" \
  --de_alpha 0.05 \
  --de_lfc 1.0 \
  --de_fdr "Benjamini Hochberg" \
  --imputation_type "Perseus-type" \
  --imputation_shift 1.8 \
  --imputation_scale 0.3 \
  --qc_plot_data "nonimputed" \
  --sample_cvs_full_range "false" \
  --volcano_display_names "true" \
  --volcano_show_gene "true" \
  --gene_annotations $annotations_link \
  --output_dir "output/"
```

**Peptide run:**

```bash
Rscript FragPipeAnalystR_main.R \
  --experiment_annotation "experiment_annotation_GLProteomics.tsv" \
  --quantification_file "abundance_peptide_MD.tsv" \
  --mode "TMT" \
  --level "peptide" \
  --feature_list_peptide "" \
  --top_n_peptide 10 \
  --enrichment_database "Hallmark,GO_Biological_Process_2021" \
  --enrichment_direction "Up,Down" \
  --normalization_method "none" \
  --de_alpha 0.05 \
  --de_lfc 1.0 \
  --de_fdr "Benjamini Hochberg" \
  --imputation_type "Perseus-type" \
  --imputation_shift 1.8 \
  --imputation_scale 0.3 \
  --qc_plot_data "nonimputed" \
  --sample_cvs_full_range "false" \
  --volcano_display_names "true" \
  --volcano_show_gene "true" \
  --gene_annotations $annotations_link \
  --output_dir "output/"
```

**Parameter Definitions:**

- `--experiment_annotation` – path to experiment annotation TSV file (table mapping TMT channels to samples)
- `--quantification_file` – path to TMTIntegrator abundance file (abundance_protein_MD.tsv, abundance_gene_MD.tsv, or abundance_peptide_MD.tsv)
- `--mode` – quantification mode: `LFQ`, `TMT`, or `DIA`
- `--level` – analysis level: `protein`, `gene`, or `peptide`
- `--normalization_method` – normalization method: `none`, `vsn` (Variance Stabilizing Normalization), `MD` (median subtraction), or `GN` (global median + MAD scaling)
- `--de_alpha` – adjusted p-value threshold for DE significance
- `--de_lfc` – log2 fold change threshold for DE significance
- `--de_fdr` – FDR correction: `Benjamini Hochberg` or `Local and tail area-based` 
- `--imputation_type` – imputation method: `none`, `Perseus-type`, `knn`, `MLE`, `min`, `zero`, `bpca`, `QRILC`, `MinDet`, `MinProb`, `nbavg`, `mixed` 
- `--imputation_shift` – Perseus-type: manual_impute shift in SD units 
- `--imputation_scale` – Perseus-type: manual_impute scale factor 
- `--feature_list_protein` – comma-separated protein IDs for feature plots (protein level). Empty = use `--top_n_protein` 
- `--feature_list_gene` – comma-separated gene names for feature plots. Empty = use `--top_n_gene`. Protein level only
- `--feature_list_peptide` – comma-separated peptide IDs for feature plots (peptide level). Empty = use `--top_n_peptide` 
- `--top_n_protein` – when feature_list_protein empty, plot top N most variable by protein ID 
- `--top_n_gene` – when feature_list_gene empty, plot top N most variable by gene 
- `--top_n_peptide` – when feature_list_peptide empty, plot top N most variable by peptide ID 
- `--qc_plot_data` – data for PCA, correlation, feature plots, sample CVs: `imputed` or `nonimputed`. If nonimputed has <2 complete features, PCA falls back to imputed with a warning
- `--sample_cvs_full_range` – sample CVs: `true` = full range, `false` = 0–1 
- `--volcano_display_names` – display names on significant volcano points 
- `--volcano_show_gene` – show gene names (`true`) or protein/peptide ID (`false`) on volcano. Peptide level uses Index; set `false` for peptide
- `--enrichment_database` – Enrichr database(s): `GO_Biological_Process_2021`, `Hallmark`, `KEGG_2021_Human`, `Reactome_2022`, etc. Comma-separated for multiple. Empty = skip 
- `--enrichment_direction` – enrichment direction(s): `Up`, `Down`, or comma-separated (e.g. `Up,Down`) 
- `--gsea_database` – GSEA database(s): `Hallmark`, `GO_Biological_Process_2021`, `GO_Cellular_Component_2021`, `GO_Molecular_Function_2021`, `KEGG_2021_Human`. Comma-separated. Protein/gene/site only. Empty = skip
- `--gene_annotations` – path or URL of gene annotations TSV/CSV; merges into DE_results on Gene. Empty = skip
- `--assay_suffix` – assay suffix for output filenames; empty = no suffix
- `--output_dir` – output directory for results

**Input Data:**

- experiment_annotation_GLProteomics.tsv (experiment annotation file, output from [Step 3b](#3b-create-manifest-and-experiment-annotation))
- abundance_protein_MD.tsv (TMTIntegrator protein-level abundance table, output from [Step 4l](#4l-tmtintegrator-tmt-quantification))
- abundance_gene_MD.tsv (TMTIntegrator gene-level abundance table, output from [Step 4l](#4l-tmtintegrator-tmt-quantification))
- abundance_peptide_MD.tsv (TMTIntegrator peptide-level abundance table, output from [Step 4l](#4l-tmtintegrator-tmt-quantification))
- annotations_link (variable containing URL of GeneLab gene annotation table for the organism; output from [Step 3c](#3c-get-organism-specific-gene-annotations-table))

**Output Data:**

- **FragPipeAnalystR_parameters.txt** (run parameters)
- **nonimputed_matrix.csv** (from abundance_protein_MD.tsv / abundance_gene_MD.tsv / abundance_peptide_MD.tsv: contaminants removed. NAs where feature not detected.)
- **imputed_matrix.csv** (same structure as nonimputed_matrix; NAs filled by Perseus-type imputation: missing values replaced with random numbers sampled from a normal distribution with mean shifted 1.8 standard deviations below and a width (SD) of 0.3, per sample.)
- **QC_plots.zip** (QC plots folder)
  - pca.pdf, .png (PCA plot)
  - missing_value_heatmap.pdf, .png (missing value pattern heatmap)
  - feature_numbers.pdf, .png (feature count per sample)
  - sample_cvs.pdf, .png (sample coefficient of variation)
  - density.pdf, .png (intensity distribution)
- **comparison_plots.zip** (comparison plots folder)
  - correlation_heatmap.pdf, .png (sample correlation heatmap)
  - feature/protein/boxplot/, feature/protein/violinplot/, feature/gene/boxplot/, feature/gene/violinplot/ (protein run: top 10 by protein ID and gene; boxplot_\*.pdf, .png and violinplot_\*.pdf, .png)
  - feature/gene/boxplot/, feature/gene/violinplot/ (gene run: top 10 by gene; boxplot_\*.pdf, .png and violinplot_\*.pdf, .png)
  - feature/peptide/boxplot/, feature/peptide/violinplot/ (peptide run: top 10 by peptide ID; boxplot_\*.pdf, .png and violinplot_\*.pdf, .png)
- **pathway_analysis_plots.zip** (pathway analysis plots folder)
  - or/ (over-representation analysis: or_database_direction.csv, .pdf, .png per database and direction)
  - gsea/ (GSEA: gsea_database_contrast.csv, .pdf, .png per database and contrast)
- **DE_plots.zip** (DE plots folder)
  - DE_heatmap.pdf, .png (DE heatmap)
  - volcano/ (volcano plots per contrast: contrast_volcano.pdf, .png)
- **SampleTable.csv** (table specifying the group or set of factor levels for each sample)
- **contrasts.csv** (table listing all pairwise group comparisons )
- **DE_results.csv** (differential expression results table; columns in order:
    - Organism-specific gene annotations
    - Protein level:
      - Index (protein-group key from the quantification table; same string as `ProteinID` at protein level)
      - NumberPSM (PSM count for the protein group)
      - MaxPepProb (maximum peptide probability among PSMs used in quantification)
      - ReferenceIntensity (log2 reference / bridge-channel intensity)
      - ProteinID (protein group identifier; duplicate of `Index` at protein level)
      - name (FragPipeAnalystR: plot/table label from `ProteinID` / `Index`)
      - ID (FragPipeAnalystR: copy of `Index`)
    - Gene level:
      - Index (gene name)
      - NumberPSM (PSMs mapping to the gene that are used in quantification)
      - ProteinID (protein identifier mapped to the gene)
      - MaxPepProb (highest PeptideProphet probability among PSMs mapping to the gene that are used in quantification)
      - ReferenceIntensity (log2 reference-channel abundance; real reference if provided, otherwise virtual reference from mean abundance across channels in the plex; global minimum reference for imputation; multi-plex: averaged across plexes)
      - name (FragPipeAnalystR: plot/table label from `ProteinID` — protein identifier mapped to the gene)
      - ID (FragPipeAnalystR: copy of `Index`)
    - Peptide level:
      - Index (FASTA protein sequence header with start and end positions of the peptide within the protein)
      - Gene (originating gene name)
      - Peptide (stripped peptide sequence)
      - NumberPSM (PSMs mapping to the peptide that are used in quantification)
      - ProteinID (protein identifier)
      - SequenceWindow (sequence window in the peptide report)
      - MaxPepProb (highest PeptideProphet probability among PSMs for this peptide sequence used in quantification)
      - ReferenceIntensity (log2 reference-channel abundance; real reference if provided, otherwise virtual reference from mean abundance across channels in the plex; global minimum reference for imputation; multi-plex: averaged across plexes)
      - name (FragPipeAnalystR: plot/table label from `ProteinID`)
      - ID (FragPipeAnalystR: copy of `Index`)
    - \* (sample / channel columns; log2 reporter abundance or ratio from the quantification matrix)
    - For each pairwise group comparison (B)v(A):
      - CI.L_(B)v(A) (lower bound of log2 fold-change confidence interval)
      - CI.R_(B)v(A) (upper bound of log2 fold-change confidence interval)
      - Log2fc_(B)v(A) (log2 fold change)
      - P.value_(B)v(A) (unadjusted p-value)
      - Adj.p.value_(B)v(A) (Benjamini-Hochberg adjusted p-value)
      - Significant_(B)v(A) (boolean at chosen FDR and fold-change thresholds)
    - significant (global; TRUE if significant in any contrast)
    - All.mean (mean across all samples)
    - All.stdev (standard deviation across all samples)
    - For each group:
      - Group.Mean_(group) (mean within group)
      - Group.Stdev_(group) (standard deviation within group))
