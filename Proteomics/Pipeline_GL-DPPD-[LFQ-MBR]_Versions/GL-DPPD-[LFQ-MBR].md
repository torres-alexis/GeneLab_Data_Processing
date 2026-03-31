# GeneLab bioinformatics processing pipeline for Mass Spectrometry-based Proteomics Data (LFQ-MBR Workflow)

> **This page holds an overview and instructions for how GeneLab processes mass spectrometry-based proteomics data using the LFQ-MBR (Label-Free Quantification with Match-Between-Runs) workflow. Exact processing commands, GL-DPPD-[STUB] version used, and processed data output files for specific datasets are provided in the [Open Science Data Repository (OSDR)](https://osdr.nasa.gov/bio/repo/).**  

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
  - [**2. Create Proteome FASTA Database**](#2-create-proteome-fasta-database)
    - [2a. Download Proteome from UniProt](#2a-download-proteome-from-uniprot)
    - [2b. Add Decoys and Contaminants to FASTA](#2b-add-decoys-and-contaminants-to-fasta)
  - [**3. Configure Metadata**](#3-configure-metadata)
    - [3a. Create Sample Runsheet](#3a-create-sample-runsheet)
    - [3b. Create Manifest and Experiment Annotation from Runsheet](#3b-create-manifest-and-experiment-annotation-from-runsheet)
    - [3c. Get organism-specific gene annotations table](#3c-get-organism-specific-gene-annotations-table)
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
    - [4k. IonQuant Label-Free Quantification](#4k-ionquant-label-free-quantification)
  - [**5. Compile FragPipe QC Reports**](#5-compile-fragpipe-qc-reports)
  - [**6. MSstats Differential Abundance Analysis**](#6-msstats-differential-abundance-analysis)
  - [**7. FragPipeAnalystR Downstream Analysis**](#7-fragpipeanalystr-downstream-analysis)

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
|MSstats|4.18.0|[https://github.com/Vitek-Lab/MSstats](https://github.com/Vitek-Lab/MSstats)|
|FragPipeAnalystR|1.1.1|[https://github.com/Nesvilab/FragPipeAnalystR](https://github.com/Nesvilab/FragPipeAnalystR)|


---

# General processing overview with example commands  

<img src="../Workflow_Documentation/NF_Proteomics/images/draft_pipeline.png" align="center" alt="Proteomics LFQ-MBR processing workflow [STUB]"/>

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
  --input sample1.mzML sample2.mzML \
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

### 3a. Create Sample Runsheet

> Note: Rather than running the command below to create the runsheet needed for processing, the runsheet may also be created manually by following the [file specification](../Workflow_Documentation/NF_Proteomics/examples/runsheet/README.md).

```bash
### Download the *ISA.zip file from the Open Science Data Repository ###

dpt-get-isa-archive \
 --accession OSD-###

### Parse the metadata from the *ISA.zip file to create a sample runsheet ###

dpt-isa-to-runsheet --accession OSD-# \
  --isa-archive *ISA.zip \
  --plugin-dir dp_tools__NF_Proteomics_LFQ/
```

**Parameter Definitions:**

- `--accession` – OSD accession ID or GLDS accession ID (`GLDS-#`), used to retrieve the urls for the ISA archive and raw data hosted in OSDR
- `--isa-archive` – Specifies the *ISA.zip file for the respective OSD dataset, downloaded in the `dpt-get-isa-archive` command
- `--plugin-dir` – Directory containing the `dp_tools` plugin used to extract runsheet fields from ISA metadata

**Input Data:**

- No input data required but the OSD (or GLDS) accession ID needs to be indicated, which is used to download the respective ISA archive

**Output Data:**

- \*ISA.zip (compressed ISA directory containing Investigation, Study, and Assay (ISA) metadata files for the respective OSD dataset, used to define sample groups — the *ISA.zip file is located in the [OSDR repository](https://osdr.nasa.gov/bio/repo/) under 'Files' → 'Study Metadata Files')

- **{OSD-Accession-ID}_Proteomics_LFQ_v{version}_runsheet.csv** (table containing metadata required for processing; version denotes the dp_tools schema used to specify the metadata to extract from the ISA archive)

<br>

### 3b. Create Manifest and Experiment Annotation from Runsheet

```bash
runsheet_to_fp_metadata.py \
  --runsheet {OSD-Accession-ID}_Proteomics_LFQ_v{version}_runsheet.csv \
  --assay_suffix _GLProteomics
```

**Parameter Definitions:**

- `--runsheet` – path to runsheet CSV (see [Runsheet Specification](../Workflow_Documentation/NF_Proteomics/examples/runsheet/README.md))
- `--assay_suffix` – assay suffix for output filenames

**Input Data:**

- {OSD-Accession-ID}_Proteomics_LFQ_v{version}_runsheet.csv (table containing file paths and metadata required for processing, output from [Step 3a](#3a-create-sample-runsheet) or created manually)

**Output Data:**

- **manifest_GLProteomics.tsv** (FragPipe input table; headerless columns in order:
  - Path (mzML basename, from runsheet `Sample Name` (`*.mzML`))
  - Experiment (FragPipe experiment string (from `Factor Value[...]` columns))
  - Bioreplicate (biological replicate replicate alphanumeric identifier (from runsheet `Bioreplicate` column if present; else sequential by `condition`))
  - Data type (data acquisition type; preset value (`DDA`)))

- **experiment_annotation_GLProteomics.tsv** (FragPipeAnalystR input table with additional `condition_name` column; columns in order:
  - file (mzML basename (`*.mzML`))
  - sample (`{Experiment}_{Bioreplicate}` (matches manifest `Experiment` and `Bioreplicate`))
  - sample_name (sample name from runsheet `Sample Name`)
  - condition (R-safe condition symbol from joined `Factor Value[...]` values)
  - condition_name (human-readable condition)
  - replicate (biological replicate replicate alphanumeric identifier (from runsheet `Bioreplicate` column if present; else sequential by `condition`)))

<br>

### 3c. Get organism-specific gene annotations table

```r
### Runsheet from Step 3b input; organism must match the value in the species column of GL-DPPD-7110-A_annotations.csv ###
runsheet_path <- "{OSD-Accession-ID}_Proteomics_LFQ_v{version}_runsheet.csv"
runsheet <- read.csv(runsheet_path, stringsAsFactors = FALSE, check.names = FALSE)
organism <- trimws(as.character(runsheet[["organism"]][1]))

### Pull in the GeneLab annotation table (GL-DPPD-7110-A_annotations.csv) ###
org_table_link <- "https://raw.githubusercontent.com/nasa/GeneLab_Data_Processing/master/GeneLab_Reference_Annotations/Pipeline_GL-DPPD-7110_Versions/GL-DPPD-7110-A/GL-DPPD-7110-A_annotations.csv"

org_table <- read.table(org_table_link, sep = ",", header = TRUE)

### URL of the organism-specific GeneLab gene annotation table ###
annotations_link <- org_table[org_table$species == organism, "genelab_annots_link"]
```

**Input Data:**

- {OSD-Accession-ID}_Proteomics_LFQ_v{version}_runsheet.csv (output from [Step 3a](#3a-create-sample-runsheet); `organism` column value must match a value in the `species` column of [GL-DPPD-7110-A_annotations.csv](../../GeneLab_Reference_Annotations/Pipeline_GL-DPPD-7110_Versions/GL-DPPD-7110-A/GL-DPPD-7110-A_annotations.csv))

**Output Data:**

- annotations_link (variable containing URL of organism-specific GeneLab gene annotation table)

<br>

---

## 4. FragPipe Processing Pipeline

<br>

### 4a. Launch FragPipe

```bash
/fragpipe_bin/fragpipe-24.0/fragpipe-24.0/bin/fragpipe \
  --headless \
  --workflow LFQ-MBR.workflow \
  --manifest manifest_GLProteomics.tsv \
  --workdir . \
  --ram 64 \
  --threads 16 \
  --config-tools-folder tools_folder
```

**Parameter Definitions:**

- `--headless` – run FragPipe in headless mode (no GUI)
- `--workflow` – path to FragPipe workflow configuration file
- `--manifest` – path to manifest TSV file containing sample information and file paths
- `--workdir` – working directory for FragPipe execution
- `--ram` – Memory (GB) allocated to FragPipe
- `--threads` – number of CPU threads allocated to FragPipe
- `--config-tools-folder` – path to folder containing FragPipe tools not included in the Docker image (MSFragger JAR, IonQuant JAR, diaTracer JAR, ext/bruker/, ext/thermo/)

**Input Data:**

- LFQ-MBR.workflow (FragPipe LFQ-MBR workflow configuration file)
- manifest_GLProteomics.tsv (manifest file with sample information and file paths, output from [Step 3b](#3b-create-manifest-and-experiment-annotation-from-runsheet))
- tools_folder/ (directory containing FragPipe tools not included in the Docker image)
- \*.mzML (input mass spectrometry raw data in mzML format)
- \*-decoys-reviewed-contam-\*.fas (proteome FASTA database with decoys and contaminants, output from [Step 2](#2-create-proteome-fasta-database))

**Output Data:**

- fragger.params (MSFragger parameter configuration file)
- msbooster_params.txt (MSBooster parameter configuration file)
- filelist_proteinprophet.txt (list of interact.pep.xml files to be passed to ProteinProphet)
- filelist_ionquant.txt (file list for IonQuant)
- modmasses_ionquant.txt (modification masses file for IonQuant)
- experiment_annotation.tsv (experiment annotation file)
- fragpipe.workflow (FragPipe output workflow configuration file)
- fragpipe-files.fp-manifest (FragPipe output manifest)
- fragpipe.job (FragPipe job configuration file)
- log_\*.txt (FragPipe execution log file with timestamp)
- sdrf.tsv (Sample and Data Relationship Format file)

<br>

### 4b. Check Spectral Files Centroid Status

```bash
java -Xmx64G -cp /fragpipe_bin/fragpipe-24.0/fragpipe-24.0/lib/fragpipe-24.0.jar:/fragpipe_bin/fragpipe-24.0/fragpipe-24.0/tools/batmass-io-1.36.5.jar org.nesvilab.fragpipe.util.CheckCentroid *.mzML 16
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
philosopher-v5.1.3-RC9 workspace --clean --nocheck
philosopher-v5.1.3-RC9 workspace --init --nocheck --temp /tmp/temp_directory
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
java -jar -Dfile.encoding=UTF-8 -Xmx64G MSFragger-4.4.1.jar fragger.params sample1.mzML sample2.mzML
```
<!-- CLI mode (backup) - same command, no changes needed for headless mode -->

**Parameter Definitions:**

- `-jar` – executes JAR file
- `-Dfile.encoding=UTF-8` – sets file encoding to UTF-8
- `-Xmx64G` – Java memory limit (e.g., `-Xmx64G` for 64 GB RAM)
- `MSFragger-4.4.1.jar` – MSFragger JAR file
- `fragger.params` – MSFragger parameter configuration file
- `*.mzML` – multiple mzML files provided as individual paths separated by spaces

**Input Data:**

- fragger.params (MSFragger parameter configuration file, output from [Step 4a](#4a-launch-fragpipe))
- *.mzML (input mass spectrometry raw data in mzML format)
- \*-decoys-reviewed-contam-\*.fas (proteome FASTA database with decoys and contaminants, output from [Step 2](#2-create-proteome-fasta-database))

**Output Data:**

- ***.pepXML** (peptide-spectrum matches from the MSFragger database search)
- ***.pin** (peptide-spectrum matches from the MSFragger database search in Percolator input format (PIN) for statistical validation)
- ***.pepindex** (peptide index files for the FASTA database)
- ***.tsv** (MSFragger results in tab-separated format)

<!-- > **Note:** MSFragger performs the database search and reports PSMs and associated search scores in pin files. See [MSFragger GitHub](https://github.com/Nesvilab/MSFragger) for details. -->

<!-- > **Note:** **PIN (Percolator Input)** files are tab-delimited files containing peptide-spectrum matches (PSMs) with features and scores. MSFragger generates PIN files with basic features (e.g., hyperscore, delta score, retention time, charge). MSBooster adds additional deep learning-based features (e.g., spectral entropy, hypergeometric probability, intersection, predicted RT, delta RT LOESS) to these PIN files. **pepXML (Peptide XML)** is an open data format developed at the SPC/Institute for Systems Biology for the storage, exchange, and processing of peptide sequence assignments of MS/MS scans. It provides a common data output format for many different MS/MS search engines and subsequent peptide-level analyses. See [pepXML format documentation](http://tools.proteomecenter.org/wiki/index.php?title=Formats:pepXML) for details. -->

<br>

### 4e. MSBooster Deep Learning Feature Addition

```bash
java -Djava.awt.headless=true -Xmx64G -cp /fragpipe_bin/fragpipe-24.0/fragpipe-24.0/tools/MSBooster-1.4.14.jar:/fragpipe_bin/fragpipe-24.0/fragpipe-24.0/tools/batmass-io-1.36.5.jar mainsteps.MainClass --paramsList msbooster_params.txt
```
<!-- CLI mode (backup):
```bash
java -Xmx64G -cp MSBooster-1.3.17.jar:batmass-io-1.35.4.jar mainsteps.MainClass --paramsList msbooster_params.txt
``` -->

**Parameter Definitions:**

- `-Djava.awt.headless=true` – runs in headless mode (no GUI)
- `-Xmx64G` – Java memory limit (e.g., `-Xmx64G` for 64 GB RAM)
- `-cp` – Java classpath to MSBooster and BatMass libraries
- `mainsteps.MainClass` – MSBooster main class
- `--paramsList` – path to MSBooster parameter configuration file

**Input Data:**

- msbooster_params.txt (MSBooster parameter configuration file, output from [Step 4a](#4a-launch-fragpipe))
- *.pin (Percolator input files from MSFragger, output from [Step 4d](#4d-msfragger-database-search))
- *.mzML (original mass spectrometry raw data in mzML format)

**Output Data:**

- *_edited.pin (Percolator input files with added deep learning features from MSBooster: unweighted spectral entropy, weighted spectral entropy, hypergeometric probability, intersection, predicted RT real units, and delta RT LOESS)
- spectraRT_full.tsv (full spectra retention time data)
- spectraRT.predicted.bin (binary file containing predicted spectra, retention times, and ion mobilities from DIA-NN)
- spectraRT.tsv (spectra retention time data)
- MSBooster_plots/ (Directory containing MSBooster calibration and diagnostic plots)

<!-- > **Note:** MSBooster extracts peptides from pin files and creates input for a deep learning model (DIA-NN in FragPipe) to predict physicochemical properties (RT, IM, and/or MS/MS spectra). Predictions are performed for candidate peptides reported by MSFragger. MSBooster generates features based on agreement between experimental and predicted values and adds them to the pin files, which are then passed to Percolator. See [MSBooster GitHub](https://github.com/Nesvilab/MSBooster) and [Yang et al. (2023) Nature Communications](https://pmc.ncbi.nlm.nih.gov/articles/PMC10374903/). -->

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
<!-- CLI mode (backup):
```bash
percolator \
  --only-psms \
  --no-terminate \
  --post-processing-tdc \
  --num-threads 16 \
  --results-psms *_percolator_target_psms.tsv \
  --decoy-results-psms *_percolator_decoy_psms.tsv \
  --protein-decoy-pattern rev_ \
  *_edited.pin
``` -->

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

- *_edited.pin (Percolator input files with MSBooster features, output from [Step 4e](#4e-msbooster-deep-learning-feature-addition))

**Output Data:**

- *_percolator_target_psms.tsv (Percolator target PSM results in TSV format)
- *_percolator_decoy_psms.tsv (Percolator decoy PSM results in TSV format)

<!-- > **Note:** Percolator learns a linear support vector machine (SVM) to differentiate true target PSMs from decoy PSMs using features in the pin files (including deep learning-based features added by MSBooster). Percolator assigns an SVM score and posterior error probability to each PSM. -->

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
<!-- CLI mode (backup) - N/A -->

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
- `*.mzML` – original mzML file path

**Input Data:**

- *.pin (original Percolator input files from MSFragger, output from [Step 4d](#4d-msfragger-database-search))
- *_percolator_target_psms.tsv (Percolator target PSM results, output from [Step 4f1](#4f1-perform-percolator-psm-rescoring-and-statistical-validation))
- *_percolator_decoy_psms.tsv (Percolator decoy PSM results, output from [Step 4f1](#4f1-perform-percolator-psm-rescoring-and-statistical-validation))
- *.mzML (original mass spectrometry raw data in mzML format)

**Output Data:**

- interact-*.pep.xml (peptide-spectrum matches with validation information generated by Percolator)

<!-- > **Note:** The temporary `*_percolator_target_psms.tsv` and `*_percolator_decoy_psms.tsv` files are deleted after conversion to pepXML format. -->

<br>

### 4g. ProteinProphet Protein Inference and Statistical Validation

```bash
/fragpipe_bin/fragpipe-24.0/fragpipe-24.0/tools/Philosopher/philosopher-v5.1.3-RC9 proteinprophet --maxppmdiff 2000000 --output combined filelist_proteinprophet.txt
```
<!-- CLI mode (backup):
```bash
philosopher proteinprophet --maxppmdiff 2000000 --output combined filelist_proteinprophet.txt
``` -->

**Parameter Definitions:**

- `proteinprophet` – run ProteinProphet to generate probabilities for protein identifications based on MS/MS data
- `--maxppmdiff` – maximum peptide mass difference in ppm
- `--output combined` – output file name
- `filelist_proteinprophet.txt` – list of interact.pep.xml files to be passed to ProteinProphet

**Input Data:**

- filelist_proteinprophet.txt (list of interact.pep.xml files to be passed to ProteinProphet, output from [Step 4a](#4a-launch-fragpipe))
- interact-*.pep.xml (pepXML files listed in filelist_proteinprophet.txt, output from [Step 4f2](#4f2-add-percolator-validation-information-to-pepxml))

**Output Data:**

- combined.prot.xml (protein identifications with validation information generated by ProteinProphet via Philosopher)

<!-- > **Note:** ProteinProphet generates probabilities for protein identifications by combining peptide identifications corresponding to the same protein and using peptide probabilities. It addresses peptide degeneracy (when one peptide corresponds to several different proteins) and groups proteins into clusters within the protXML `<protein group>` element. Proteins sharing identified peptides are grouped together, and Occam's Razor is applied to assign probabilities (often assigning probability of zero to unneeded proteins in a group to present the shortest list of proteins needed to explain the data). See [ProteinProphet documentation](http://tools.proteomecenter.org/wiki/index.php?title=Software:ProteinProphet) for details. -->

<br>

### 4h. Database Annotation

```bash
/fragpipe_bin/fragpipe-24.0/fragpipe-24.0/tools/Philosopher/philosopher-v5.1.3-RC9 database --annotate *.fas --prefix rev_
```
<!-- CLI mode (backup):
```bash
philosopher database --annotate *.fas --prefix rev_
``` -->

**Parameter Definitions:**

- `database --annotate` – annotate FASTA database file (creates binary database files for Philosopher tools)
- `*.fas` – path to FASTA database file
- `--prefix rev_` – decoy prefix used in the database

**Input Data:**

- \*-decoys-reviewed-contam-\*.fas (proteome FASTA database with decoys and contaminants, output from [Step 2](#2-create-proteome-fasta-database))

**Output Data:**

- .meta/ (Philosopher workspace metadata directory containing binary database files)

<br>

### 4i. Filter Results by FDR

```bash
# First sample (initializes database annotation)
/fragpipe_bin/fragpipe-24.0/fragpipe-24.0/tools/Philosopher/philosopher-v5.1.3-RC9 filter \
  --sequential \
  --prot 0.01 \
  --picked \
  --tag rev_ \
  --pepxml sample_directory \
  --protxml combined.prot.xml \
  --razor

# Subsequent samples (reuse database annotation from first sample)
/fragpipe_bin/fragpipe-24.0/fragpipe-24.0/tools/Philosopher/philosopher-v5.1.3-RC9 filter \
  --sequential \
  --prot 0.01 \
  --picked \
  --tag rev_ \
  --pepxml sample_directory \
  --dbbin first_sample_directory \
  --protxml combined.prot.xml \
  --probin first_sample_directory \
  --razor
```
<!-- CLI mode (backup):
```bash
# First sample (initializes database annotation)
philosopher filter \
  --sequential \
  --prot 0.01 \
  --picked \
  --tag rev_ \
  --pepxml sample_directory \
  --protxml combined.prot.xml \
  --razor

# Subsequent samples (reuse database annotation from first sample)
philosopher filter \
  --sequential \
  --prot 0.01 \
  --picked \
  --tag rev_ \
  --pepxml sample_directory \
  --dbbin first_sample_directory \
  --protxml combined.prot.xml \
  --probin first_sample_directory \
  --razor
``` -->

**Parameter Definitions:**

- `filter` – filter PSMs, peptides, and proteins by FDR threshold
- `--sequential` – apply sequential FDR filtering at PSM, peptide, and ion levels in addition to protein level FDR
- `--prot 0.01` – protein-level FDR threshold
- `--picked` – apply picked FDR algorithm prior to protein scoring
- `--tag rev_` – decoy sequence prefix
- `--pepxml` – path to pepXML file(s) or directory containing pepXML files
- `--protxml combined.prot.xml` – path to protXML file
- `--dbbin` – (for subsequent samples) path to first sample directory containing database annotation
- `--probin` – (for subsequent samples) path to first sample directory containing protein annotation
- `--razor` – use razor peptides for protein-level FDR scoring

**Input Data:**

- interact-*.pep.xml (peptide-spectrum matches with validation information generated by Percolator, output from [Step 4f2](#4f2-add-percolator-validation-information-to-pepxml))
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
<!-- CLI mode (backup):
```bash
philosopher report
``` -->

**Input Data:**

- Philosopher workspace containing filtered data (output from [Step 4i](#4i-filter-results-by-fdr))

**Output Data:**

- protein.fas (FASTA file containing FDR-filtered protein sequences identified)
- protein.tsv (sample-specific protein report)
- peptide.tsv (sample-specific peptide report)
- psm.tsv (sample-specific PSM report)
- ion.tsv (sample-specific ion report)

<br>

### 4k. IonQuant Label-Free Quantification

```bash
java -Djava.awt.headless=true -Xmx64G \
  -Dlibs.bruker.dir=tools/ext/bruker \
  -Dlibs.thermo.dir=tools/ext/thermo \
  -cp /fragpipe_bin/fragpipe-24.0/fragpipe-24.0/tools/jfreechart-1.5.3.jar \
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
  --mbr 1 \
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
  --intensitymode 0 \
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
```
<!-- CLI mode (backup):
```bash
java -Xmx64G \
  -Dlibs.bruker.dir=tools/ext/bruker \
  -Dlibs.thermo.dir=tools/ext/thermo \
  -cp jfreechart-1.5.3.jar:IonQuant-1.11.11.jar \
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
  --mbr 1 \
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
  --intensitymode 0 \
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
``` -->

**Parameter Definitions:**

- `-Djava.awt.headless=true` – run in headless mode (no GUI)
- `-Xmx64G` – Java memory limit (e.g., `-Xmx64G` for 64 GB RAM)
- `-Dlibs.bruker.dir` – directory for Bruker libraries
- `-Dlibs.thermo.dir` – directory for Thermo libraries
- `-cp` – Java classpath to jfreechart and IonQuant JAR files
- `ionquant.IonQuant` – IonQuant main class
- `--threads` – number of CPU threads to use (0 = all logical cores)
- `--perform-ms1quant 1` – perform MS1 quantification (0 = no, 1 = yes)
- `--perform-isoquant 0` – perform isobaric labeling quantification (0 = no, 1 = yes)
- `--mbr 1` – perform match-between-runs (0 = no, 1 = yes)
- `--maxlfq 1` – calculate MaxLFQ intensity (0 = no, 1 = yes)
- `--msstats 1` – generate MSstats input files (0 = no, 1 = yes)
- `--site-reports 1` – generate site reports (0 = no, 1 = yes; requires modification localization columns in psm.tsv)
- `--multidir .` – output directory for multi-experimental results (optional)
- `--filelist` – file containing flags (tab-delimited file with `--psm` entries pointing to sample-specific psm.tsv files and `--specdir` entry pointing to the directory containing mzML files)
- `--modlist` – file listing modification masses (used to remove mass discrepancy due to rounding errors)
- `--specdir` – directory containing spectral files (mzML/mzXML/raw/quantindex); can specify multiple

**Input Data:**

- filelist_ionquant.txt (file list for IonQuant, output from [Step 4a](#4a-launch-fragpipe))
- modmasses_ionquant.txt (modification masses file for IonQuant, output from [Step 4a](#4a-launch-fragpipe))
- protein.tsv (sample-specific protein report, output from [Step 4j](#4j-generate-reports))
- peptide.tsv (sample-specific peptide report, output from [Step 4j](#4j-generate-reports))
- psm.tsv (sample-specific PSM report, output from [Step 4j](#4j-generate-reports))
- ion.tsv (sample-specific ion report, output from [Step 4j](#4j-generate-reports))
- *.mzML (original mass spectrometry raw data in mzML format; accessed via `--specdir` parameter specified in filelist_ionquant.txt to extract intensity data for MS1 quantification and match-between-runs feature matching)

**Output Data:**

- protein.tsv (sample-specific protein report with MS1 quantification data added from IonQuant)
- peptide.tsv (sample-specific peptide report with MS1 quantification data and additional data added from IonQuant)
- ion.tsv (sample-specific ion report with MS1 quantification data and additional data added from IonQuant)
- psm.tsv (sample-specific PSM report with MS1 quantification data and additional data added from IonQuant)
- *_model.png (sample-specific IonQuant model visualization plot showing quantification model fits)
- **combined_protein.tsv** (combined protein report with MS1 quantification data across all samples)
- **combined_peptide.tsv** (combined peptide report with MS1 quantification data and additional data across all samples)
- **combined_modified_peptide.tsv** (combined modified peptide report with MS1 quantification data and additional data across all samples)
- **combined_ion.tsv** (combined ion report with MS1 quantification data and additional data across all samples)
- **combined_site_*.tsv** (site-specific modification reports, e.g., combined_site_C_57.0215.tsv for carbamidomethylation, combined_site_M_15.9949.tsv for oxidation)
- reprint.int.tsv (input file for the Resource for Evaluation of Protein Interaction Networks (REPRINT) containing protein intensities)
- reprint.spc.tsv (input file for the Resource for Evaluation of Protein Interaction Networks (REPRINT) containing protein spectral counts)
- **msstats.csv** (MSstats input file for downstream differential analysis)
- **msstats_ptm.csv** (MSstatsPTM input file for PTM (post-translational modification) analysis)
<!-- - *.mbrbin (match-between-runs binary data file for MBR feature matching) -->
<!-- - *.quantbin2 (quantification binary cache file) -->

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

- psm.tsv (sample-specific PSM reports, output from [Step 4k](#4k-ionquant-label-free-quantification))
- ion.tsv (sample-specific ion reports, output from [Step 4k](#4k-ionquant-label-free-quantification))
- combined_protein.tsv (combined protein report, output from [Step 4k](#4k-ionquant-label-free-quantification))
- combined_peptide.tsv (combined peptide report, output from [Step 4k](#4k-ionquant-label-free-quantification))
- combined_ion.tsv (combined ion report, output from [Step 4k](#4k-ionquant-label-free-quantification))
- *.workflow (FragPipe workflow file, output from [Step 4a](#4a-launch-fragpipe))
- fragger.params (MSFragger parameters file, output from [Step 4a](#4a-launch-fragpipe))

**Output Data:**

- **multiqc_GLProteomics.html** (MultiQC output html summary)
- **multiqc_GLProteomics_data.zip** (zipped directory containing MultiQC output data with cleaned paths)

<br>

---

## 6. MSstats Differential Abundance Analysis

```bash
msstats_analysis.R . experiment_annotation_GLProteomics.tsv msstats.csv _GLProteomics
```

**Parameter Definitions:**

- `msstats_analysis.R` – R script for MSstats differential abundance analysis
- `.` – root directory for output
- `experiment_annotation_GLProteomics.tsv` – experiment annotation (sample metadata, condition assignments)
- `msstats.csv` – MSstats input file from IonQuant
- `_GLProteomics` – assay suffix: stripped from Run column for matching; appended to output filenames. 

**Input Data:**

- msstats.csv (MSstats input file, output from [Step 4k](#4k-ionquant-label-free-quantification))
- experiment_annotation_GLProteomics.tsv (sample metadata and condition assignments)

**Output Data:**


- **msstats_comparison_GLProteomics.csv** (all MSstats pairwise comparisons)
- **msstats_contrasts_GLProteomics.csv** (contrast definitions)

<br>

---

## 7. FragPipeAnalystR Downstream Analysis

The FragPipeAnalystR downstream analysis script is executed twice: once using the **protein**-level quantification file (combined_protein.tsv) and once using the **peptide**-level quantification file (combined_peptide.tsv).

**Protein run:**

```bash
Rscript FragPipeAnalystR_main.R \
  --experiment_annotation "experiment_annotation_GLProteomics.tsv" \
  --quantification_file "combined_protein.tsv" \
  --mode "LFQ" \
  --level "protein" \
  --feature_list_protein "" \
  --feature_list_gene "" \
  --top_n_protein 10 \
  --top_n_gene 10 \
  --enrichment_database "Hallmark,GO_Biological_Process_2021" \
  --enrichment_direction "Up,Down" \
  --gsea_database "Hallmark,GO_Biological_Process_2021" \
  --lfq_type "Intensity" \
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
  --quantification_file "combined_peptide.tsv" \
  --mode "LFQ" \
  --level "peptide" \
  --feature_list_peptide "" \
  --top_n_peptide 10 \
  --enrichment_database "Hallmark,GO_Biological_Process_2021" \
  --enrichment_direction "Up,Down" \
  --lfq_type "Intensity" \
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

- `--experiment_annotation` – path to experiment annotation TSV file (sample metadata and condition assignments)
- `--quantification_file` – path to combined quantification file (combined_protein.tsv or combined_peptide.tsv, output from [Step 4k](#4k-ionquant-label-free-quantification))
- `--mode` – quantification mode: `LFQ`, `TMT`, or `DIA`
- `--level` – analysis level: `protein` or `peptide`
- `--lfq_type` – LFQ column type: `Intensity`, `MaxLFQ`, or `Spectral Count`
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

- experiment_annotation_GLProteomics.tsv (experiment annotation file, output from [Step 3b](#3b-create-manifest-and-experiment-annotation-from-runsheet))
- combined_protein.tsv (combined protein report, output from [Step 4k](#4k-ionquant-label-free-quantification))
- combined_peptide.tsv (combined peptide report, output from [Step 4k](#4k-ionquant-label-free-quantification))
- annotations_link (variable containing URL of GeneLab gene annotation table for the organism)

**Output Data:**

- **FragPipeAnalystR_parameters.txt** (run parameters)
- **nonimputed_matrix.csv** (from combined_protein/peptide: contaminants removed; selected `--lfq_type` quantification columns reappended at the end of the table as either log2 intensity (Intensity/MaxLFQ) or raw counts (Spectral Count). NAs where feature not detected.)
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
  - feature/peptide/boxplot/, feature/peptide/violinplot/ (peptide run: top 10 by peptide ID; boxplot_\*.pdf, .png and violinplot_\*.pdf, .png)
<!--   - feature/site/boxplot/, feature/site/violinplot/ (site run: top N by site ID; boxplot_*.pdf, .png and violinplot_*.pdf, .png)-->
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
      - Protein (protein sequence header from search database FASTA; razor protein when peptide maps to multiple)
      - Protein.ID (UniProt primary accession; second pipe-delimited field of Protein)
      - Entry.Name (UniProt entry name; third pipe-delimited field of Protein)
      - Protein.Length (number of amino acids in protein)
      - Organism (species)
      - Protein.Existence (UniProt evidence type)
      - Description (protein name)
      - Protein.Probability (ProteinProphet confidence)
      - Top.Peptide.Probability (highest PeptideProphet score among mapped peptides)
      - Combined.Total.Peptides (number of peptides mappable to protein)
      - Combined.Spectral.Count (PSMs for razor peptides)
      - Combined.Unique.Spectral.Count (PSMs for unique peptides)
      - Combined.Total.Spectral.Count (PSMs for total peptides)
      - *.Spectral.Count (per-sample PSM count)
      - *.Unique.Spectral.Count (per-sample PSMs for unique peptides only)
      - *.Total.Spectral.Count (per-sample PSMs for unique and razor peptides)
      - *.MaxLFQ.Intensity (per-sample MaxLFQ-normalized intensity)
      - Indistinguishable.Proteins (proteins not distinguishable given evidence)
      - name (Gene)
      - ID (Protein ID)
    - Peptide level:
      - Peptide.Sequence (stripped sequence, no modifications)
      - Prev.AA (residue preceding peptide in protein)
      - Next.AA (residue following peptide in protein)
      - Start (position of peptide start in protein)
      - End (position of peptide end in protein)
      - Peptide.Length (number of residues)
      - Charges (observed charge states)
      - Protein (protein sequence header from search database FASTA; razor protein when peptide maps to multiple)
      - Protein.ID (UniProt primary accession; second pipe-delimited field of Protein)
      - Entry.Name (UniProt entry name; third pipe-delimited field of Protein)
      - Description (protein name of parent protein)
      - Mapped.Genes (additional genes peptide may originate from)
      - Mapped.Proteins (additional proteins peptide maps to)
      - *.Spectral.Count (per-sample PSM count)
      - *.MaxLFQ.Intensity (per-sample MaxLFQ-normalized intensity)
      - *.Match.Type (per-sample; direct = observed in run, transferred = matched between runs)
      - Index (Protein ID + Peptide Sequence, unique peptide identifier)
      - name (Protein ID)
      - ID (Index, Protein ID + Peptide Sequence)
    - \* (sample ID; assay columns; log2 precursor intensity when lfq_type=Intensity)
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
