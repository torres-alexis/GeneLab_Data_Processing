# GeneLab bioinformatics processing pipeline for Mass Spectrometry-based Proteomics Data (TMT16 Workflow)

> **This page holds an overview and instructions for how GeneLab processes mass spectrometry-based proteomics data using the TMT16 (Tandem Mass Tag 16-plex) workflow. Exact processing commands, GL-DPPD-[STUB] version used, and processed data output files for specific datasets are provided in the [Open Science Data Repository (OSDR)](https://osdr.nasa.gov/bio/repo/).**  

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
    - [4e. Percolator PSM Rescoring and Statistical Validation](#4e-percolator-psm-rescoring-and-statistical-validation)
        - [4e1. Perform Percolator PSM Rescoring and Statistical Validation](#4e1-perform-percolator-psm-rescoring-and-statistical-validation)
        - [4e2. Add Percolator Validation Information to pepXML](#4e2-add-percolator-validation-information-to-pepxml)
    - [4f. ProteinProphet Protein Inference and Statistical Validation](#4f-proteinprophet-protein-inference-and-statistical-validation)
    - [4g. Database Annotation](#4g-database-annotation)
    - [4h. Filter Results by FDR](#4h-filter-results-by-fdr)
    - [4i. Generate Reports](#4i-generate-reports)
    - [4j. IonQuant TMT Reporter Ion Extraction](#4j-ionquant-tmt-reporter-ion-extraction)
    - [4k. TMTIntegrator TMT Quantification](#4k-tmtintegrator-tmt-quantification)
  - [**5. Compile FragPipe QC Reports**](#5-compile-fragpipe-qc-reports)

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

<img src="../Workflow_Documentation/NF_Proteomics/images/draft_pipeline.png" align="center" alt="Proteomics TMT-16 processing workflow [STUB]"/>

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

- \*-decoys-reviewed-contam-\*.fas (FASTA database containing the proteome with reversed decoy sequences and common contaminants added)

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

- **manifest_GLProteomics.tsv** (FragPipe input table; headerless columns in order:)
  - Path (mzML basename from data sheet `run` (`*.mzML`))
  - Experiment (plex identifier from data sheet `plex`)
  - Bioreplicate (mixture technical replicate alphanumeric identifier (from data sheet `TechRepMixture`, or `1` if empty))
  - Data type (preset value (`DDA`))

- **experiment_annotation_GLProteomics.tsv** (FragPipeAnalystR input table with additional `condition_name` column; columns in order:)
  - plex (`{Experiment}_{Bioreplicate}` (matches manifest `Experiment` and `Bioreplicate`))
  - channel (TMT reporter channel from sample sheet `channel`)
  - sample (R-safe sample identifier from sample sheet `Sample Name`)
  - sample_name (same as `sample`)
  - condition (R-safe condition symbol from joined sample sheet `Factor Value[...]` values)
  - condition_name (human-readable condition)
  - replicate (biological replicate replicate alphanumeric identifier (from sample sheet `Bioreplicate` column if present; else `1`))

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
bash tmt_stage_by_plex.sh manifest_GLProteomics.tsv experiment_annotation_GLProteomics.tsv TMT16
```

**Parameter Definitions:**

- `manifest_GLProteomics.tsv` – path to FragPipe manifest table
- `experiment_annotation_GLProteomics.tsv` – path to FragPipe experiment_annotation table 
- `TMT16` – isobaric multiplex type; sets reporter channel order and row count for each annotation table created

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
  --workflow TMT16.workflow \
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

- TMT16.workflow (FragPipe workflow configuration file for TMT-16 workflow)
- manifest_GLProteomics.tsv (FragPipe input table, output from [Step 3d](#3d-stage-mzml-by-plex-and-create-channel-sample-annotation))
- \*-decoys-reviewed-contam-\*.fas (proteome FASTA database, output from [Step 2](#2-download-reference-proteome-add-decoys-and-contaminants-to-fasta))
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

### 4e. Percolator PSM Rescoring and Statistical Validation

<br>

#### 4e1. Perform Percolator PSM Rescoring and Statistical Validation

```bash
/fragpipe_bin/fragpipe-24.0/fragpipe-24.0/tools/percolator_3_7_1/linux/percolator \
  --only-psms \
  --no-terminate \
  --post-processing-tdc \
  --num-threads 16 \
  --results-psms *_percolator_target_psms.tsv \
  --decoy-results-psms *_percolator_decoy_psms.tsv \
  --protein-decoy-pattern rev_ \
  *.pin
```

**Parameter Definitions:**

- `--only-psms` – do not remove redundant peptides, keep PSMs, exclude peptide level probabilities
- `--no-terminate` – do not terminate execution when encountering issues with SVM inputs or results
- `--post-processing-tdc` – replace mix-max method with target-decoy competition for assigning q-values and PEPs. For input PSMs from separate target/decoy searches, Percolator SVM scores eliminate lower-scoring target or decoy PSMs for each scan+expMass combination. Automatically enabled for concatenated searches
- `--num-threads` – number of CPU threads to use
- `--results-psms` – output file path for target PSM results
- `--decoy-results-psms` – output file path for decoy PSM results
- `--protein-decoy-pattern` – text pattern used to identify decoy proteins in the database
- `*.pin` – input Percolator input files from MSFragger

**Input Data:**

- \*.pin (Percolator input files from MSFragger, output from [Step 4d](#4d-msfragger-database-search))

**Output Data:**

- \*_percolator_target_psms.tsv (Percolator target PSM results in TSV format)
- \*_percolator_decoy_psms.tsv (Percolator decoy PSM results in TSV format)

<br>

#### 4e2. Add Percolator Validation Information to pepXML

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
- \*_percolator_target_psms.tsv (Percolator target PSM results, output from [Step 4e1](#4e1-perform-percolator-psm-rescoring-and-statistical-validation))
- \*_percolator_decoy_psms.tsv (Percolator decoy PSM results, output from [Step 4e1](#4e1-perform-percolator-psm-rescoring-and-statistical-validation))
- \*.mzML (input mass spectrometry raw data in mzML format)

**Output Data:**

- interact-\*.pep.xml (peptide-spectrum matches with validation information generated by Percolator)

<br>

### 4f. ProteinProphet Protein Inference and Statistical Validation

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
- interact-\*.pep.xml (pepXML files listed in filelist_proteinprophet.txt, output from [Step 4e2](#4e2-add-percolator-validation-information-to-pepxml))

**Output Data:**

- combined.prot.xml (protein identifications with validation information generated by ProteinProphet via Philosopher)

<br>

### 4g. Database Annotation

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

### 4h. Filter Results by FDR

```bash
# First plex (initializes database annotation)
/fragpipe_bin/fragpipe-24.0/fragpipe-24.0/tools/Philosopher/philosopher-v5.1.3-RC9 filter \
  --sequential \
  --picked \
  --prot 0.01 \
  --tag rev_ \
  --pepxml plex_directory \
  --protxml combined.prot.xml \
  --razor

# Subsequent plexes (reuse database annotation from first plex)
/fragpipe_bin/fragpipe-24.0/fragpipe-24.0/tools/Philosopher/philosopher-v5.1.3-RC9 filter \
  --sequential \
  --picked \
  --prot 0.01 \
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

- interact-\*.pep.xml (peptide-spectrum matches with validation information generated by Percolator, output from [Step 4e2](#4e2-add-percolator-validation-information-to-pepxml))
- combined.prot.xml (protein identifications with validation information generated by ProteinProphet via Philosopher, output from [Step 4f](#4f-proteinprophet-protein-inference-and-statistical-validation))
- .meta/ (Philosopher workspace metadata, output from [Step 4g](#4g-database-annotation))

**Output Data:**

- filter.log (Philosopher filter execution log file)
- Filtered data stored in Philosopher workspace as binary files (db.bin, ion.bin, pep.bin, pro.bin, protxml.bin, psm.bin, razor.bin, etc.; filtered PSM, peptide, and protein data ready for report generation)

<br>

### 4i. Generate Reports

```bash
/fragpipe_bin/fragpipe-24.0/fragpipe-24.0/tools/Philosopher/philosopher-v5.1.3-RC9 report
```

**Input Data:**

- Philosopher workspace containing filtered data (output from [Step 4h](#4h-filter-results-by-fdr))

**Output Data:**

- protein.fas (FASTA file containing FDR-filtered protein sequences identified)
- protein.tsv (plex-specific protein report)
- peptide.tsv (plex-specific peptide report)
- psm.tsv (plex-specific PSM report)
- ion.tsv (plex-specific ion report)

<br>

### 4j. IonQuant TMT Reporter Ion Extraction

```bash
# MS1
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

# Second pass: Isobaric TMT reporter ion extraction
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
  --isotype TMT-16 \
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
- `--perform-ms1quant 1` – perform MS1 quantification (first IonQuant step; 0 = no, 1 = yes)
- `--perform-isoquant 1` – perform isobaric labeling quantification (second IonQuant step; 0 = no, 1 = yes)
- `--isotol 20.0` – MS2 tolerance in ppm for isobaric quantification
- `--isolevel 2` – isobaric quantification level (2 = MS2, 3 = MS3, 4 = ZOOM-HR)
- `--isotype tmt10` / `--isotype TMT-16` – isobaric quantification type (MS1 step: `tmt10`; reporter step: `TMT-16`)
- `--ionmobility 0` – data has ion mobility information (0 = no, 1 = yes)
- `--site-reports 0` – generate site reports (0 = no, 1 = yes)
- `--msstats 0` – generate MSstats input files (0 = no, 1 = yes)
- `--annotation` – annotation file for isobaric quantification (format: `{plex}_{TechRepMixture}/psm.tsv={plex}_{TechRepMixture}/annotation.txt`; can specify multiple)
- `--minexps 1` – minimum experiments in picking an ion for quantifying proteins (only for intensity, not MaxLFQ)
- `--mbr 0` – perform match-between-runs (0 = no, 1 = yes)
- `--maxlfq 0` – calculate MaxLFQ intensity (0 = no, 1 = yes)
- `--requantify 0` – re-quantify unidentified features based on identified features (0 = no, 1 = yes;)
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
- `--multidir .` – output directory for multi-experimental results (optional)
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

### 4k. TMTIntegrator TMT Quantification

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
- psm.tsv (PSM reports with TMT reporter ion intensities, output from [Step 4j](#4j-ionquant-tmt-reporter-ion-extraction))

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

- psm.tsv (plex-specific PSM reports, output from [Step 4j](#4j-ionquant-tmt-reporter-ion-extraction))
- ion.tsv (plex-specific ion reports, output from [Step 4j](#4j-ionquant-tmt-reporter-ion-extraction))
- combined_protein.tsv (combined protein report, output from [Step 4j](#4j-ionquant-tmt-reporter-ion-extraction))
- combined_peptide.tsv (combined peptide report, output from [Step 4j](#4j-ionquant-tmt-reporter-ion-extraction))
- combined_ion.tsv (combined ion report, output from [Step 4j](#4j-ionquant-tmt-reporter-ion-extraction))
- fragpipe.workflow (FragPipe workflow configuration file, output from [Step 4a](#4a-launch-fragpipe))
- fragger.params (MSFragger parameters file, output from [Step 4a](#4a-launch-fragpipe))

**Output Data:**

- **multiqc_GLProteomics.html** (MultiQC output html summary)
- **multiqc_GLProteomics_data.zip** (zipped directory containing MultiQC output data with cleaned paths)

<br>
