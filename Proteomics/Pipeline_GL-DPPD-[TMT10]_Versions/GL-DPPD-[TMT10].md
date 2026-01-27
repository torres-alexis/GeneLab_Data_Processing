# GeneLab bioinformatics processing pipeline for Mass Spectrometry-based Proteomics Data (TMT-10 Workflow)

> **This page holds an overview and instructions for how GeneLab processes mass spectrometry-based proteomics data using the TMT-10 (Tandem Mass Tag 10-plex) workflow. Exact processing commands, GL-DPPD-[TMT10] version used, and processed data output files for specific datasets are provided in the [Open Science Data Repository (OSDR)](https://osdr.nasa.gov/bio/repo/).**  

---

**Date:** January X, 2026  
**Revision:** A  
**Document Number:** GL-DPPD-[TMT10]-A  

**Submitted by:**  
Alexis Torres (GeneLab Data Processing Team)  

**Approved by:**  
TBD

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
  - [**3. FragPipe Processing Pipeline**](#3-fragpipe-processing-pipeline)
    - [3a. Launch FragPipe](#3a-launch-fragpipe)
    - [3b. Check Spectral Files Centroid Status](#3b-check-spectral-files-centroid-status)
    - [3c. Initialize Workspace](#3c-initialize-workspace)
    - [3d. MSFragger Database Search](#3d-msfragger-database-search)
    - [3e. MSBooster Deep Learning Feature Addition](#3e-msbooster-deep-learning-feature-addition)
    - [3f. Percolator PSM Rescoring and Statistical Validation](#3f-percolator-psm-rescoring-and-statistical-validation)
        - [3f.1. Add Percolator Validation Information to pepXML](#3f1-add-percolator-validation-information-to-pepxml)
    - [3g. ProteinProphet Protein Inference and Statistical Validation](#3g-proteinprophet-protein-inference-and-statistical-validation)
    - [3h. Database Annotation](#3h-database-annotation)
    - [3i. Filter Results by FDR](#3i-filter-results-by-fdr)
    - [3j. Generate Reports](#3j-generate-reports)
    - [3k. IonQuant TMT Reporter Ion Extraction](#3k-ionquant-tmt-reporter-ion-extraction)
    - [3l. TMTIntegrator TMT Quantification](#3l-tmtintegrator-tmt-quantification)
  - [**4. Compile FragPipe QC Reports**](#4-compile-fragpipe-qc-reports)

---

# Software used  

|Program|Version|Relevant Links|
|:------|:------:|:-------------|
|dp_tools|1.3.9|[https://github.com/torres-alexis/dp_tools](https://github.com/torres-alexis/dp_tools)|
|rawBeans|1.6.4|[https://github.com/torres-alexis/rawBeans](https://github.com/torres-alexis/rawBeans)|
|FragPipe|23.1|[https://fragpipe.nesvilab.org/](https://fragpipe.nesvilab.org/)|
|BatMass|1.35.4|[https://batmass.org/](https://batmass.org/)|
|MSFragger|4.3|[http://msfragger-upgrader.nesvilab.org/upgrader/](http://msfragger-upgrader.nesvilab.org/upgrader/)|
|MSBooster|1.3.17|[https://github.com/Nesvilab/MSBooster](https://github.com/Nesvilab/MSBooster)|
|DIA-NN|1.8.2 Beta 8|[https://github.com/vdemichev/DiaNN](https://github.com/vdemichev/DiaNN)|
|Percolator|3.7.1|[https://github.com/percolator/percolator](https://github.com/percolator/percolator)|
|Philosopher|5.1.2|[https://github.com/Nesvilab/philosopher/releases/latest](https://github.com/Nesvilab/philosopher/releases/latest)|
|IonQuant|1.11.11|[https://github.com/Nesvilab/IonQuant/releases/latest](https://github.com/Nesvilab/IonQuant/releases/latest)|
|TMTIntegrator|6.1.1|[https://github.com/Nesvilab/TMTIntegrator](https://github.com/Nesvilab/TMTIntegrator)|
|MultiQC|1.32|[https://multiqc.info/](https://multiqc.info/)|
|pmultiqc|0.0.40|[https://github.com/bigbio/pmultiqc](https://github.com/bigbio/pmultiqc)|


---

# General processing overview with example commands  

<img src="../Workflow_Documentation/NF_Proteomics/images/draft_pipeline.png" align="center" alt="Proteomics TMT-10 processing workflow"/>

> Exact processing commands and output files listed in **bold** below are included with each relevant mass spectrometry-based proteomics processed dataset in the [Open Science Data Repository (OSDR)](https://osdr.nasa.gov/bio/repo/). 

---

## 1. Raw Data QC  

### 1a. RawBeans QC (Samplewise)

```bash
create-qc-report.py \
  --input *.mzML \
  --output-dir . \
  --batch \
  --cores NumberOfThreads

cd *
zip -r ../*-report.zip qc-report.html resources/
cd ..
```

**Parameter Definitions:**

- `--input` – one mzML file path
- `--output-dir` – the output directory to store results
- `--batch` – process file in batch mode (creates subdirectory for output)
- `--cores` – number of CPU cores to use for processing

**Input Data:**

- `*.mzML` (input mass spectrometry raw data in mzML format)

**Output Data:**

- `qc-report.html` (interactive HTML QC report with metrics and visualizations)
- `resources/` (directory containing QC report assets: CSS, JavaScript, and data files)

<br>

### 1b. RawBeans QC (All Samples)

```bash
create-qc-report.py \
  --input *.mzML \
  --output-dir . \
  --cores NumberOfThreads

zip -r *-report.zip qc-report.html resources/
```

**Parameter Definitions:**

- `--input` – multiple mzML file paths (all samples)
- `--output-dir` – the output directory to store results
- `--cores` – number of CPU cores to use for processing

**Input Data:**

- `*.mzML` (all input mass spectrometry raw data files in mzML format)

**Output Data:**

- `qc-report.html` (interactive HTML QC report with metrics and visualizations across all samples)
- `resources/` (directory containing QC report assets: CSS, JavaScript, and data files)

<br>

## 2. Create Proteome FASTA Database  

### 2a. Download Proteome from UniProt

```bash
philosopher-v5.1.2 workspace --init --nocheck
philosopher-v5.1.2 database --custom UP000005640 --reviewed --isoforms --contam
```

**Parameter Definitions:**

- `workspace --init` – initializes a Philosopher workspace
- `--nocheck` – skips workspace validation checks
- `database --custom` – downloads proteome from UniProt using UniProt ID
- `--reviewed` – includes only reviewed (Swiss-Prot) entries
- `--isoforms` – includes protein isoforms
- `--contam` – adds common contaminants to the database

**Input Data:**

- UniProt proteome ID (e.g., `UP000005640` for *Homo sapiens*)

**Output Data:**

- `*-reviewed-isoforms-contam-*.fas` (proteome FASTA file with reviewed entries, isoforms, and contaminants)

<br>

### 2b. Add Decoys and Contaminants to FASTA

```bash
philosopher-v5.1.2 database --annotate --custom UP000005640 --reviewed --isoforms --contam
```

**Parameter Definitions:**

- `database --annotate` – adds decoy sequences and contaminants to existing FASTA
- `--custom` – UniProt proteome ID
- `--reviewed` – includes only reviewed entries
- `--isoforms` – includes protein isoforms
- `--contam` – adds common contaminants

**Input Data:**

- `*-reviewed-isoforms-contam-*.fas` (proteome FASTA file, output from [Step 2a](#2a-download-proteome-from-uniprot))

**Output Data:**

- `*-decoys-reviewed-isoforms-contam-*.fas` (proteome FASTA database with decoy sequences (prefixed with `rev_`), reviewed entries, isoforms, and contaminants)

<br>

## 3. FragPipe Processing Pipeline

### 3a. Launch FragPipe

```bash
/fragpipe_bin/fragpipe-23.1/fragpipe-23.1/bin/fragpipe \
  --headless \
  --workflow TMT10.workflow \
  --manifest manifest_GLProteomics.tsv \
  --workdir . \
  --config-tools-folder tools \
  --config-python /usr/bin/python3.11
```

**Parameter Definitions:**

- `--headless` – runs FragPipe in headless mode (no GUI)
- `--workflow` – path to FragPipe workflow configuration file
- `--manifest` – path to FragPipe manifest TSV file (contains sample information and file paths)
- `--workdir` – working directory for FragPipe execution
- `--config-tools-folder` – directory containing FragPipe tools
- `--config-python` – path to Python executable

**Input Data:**

- `TMT10.workflow` (FragPipe workflow configuration file for TMT-10 workflow)
- `manifest_GLProteomics.tsv` (FragPipe manifest TSV file with columns: input_file, joined_factor_values, bioreplicate, data_type)
- `*-decoys-reviewed-isoforms-contam-*.fas` (proteome FASTA database, output from [Step 2b](#2b-add-decoys-and-contaminants-to-fasta))
- `*.mzML` (input mass spectrometry raw data in mzML format)

**Output Data:**

- `fragger.params` (MSFragger parameter configuration file)
- `msbooster_params.txt` (MSBooster parameter configuration file)
- `tmt-integrator-conf.yml` (TMTIntegrator configuration file)
- `filelist_proteinprophet.txt` (list of interact.pep.xml files to be passed to ProteinProphet)
- `filelist_ionquant.txt` (file list for IonQuant)
- `modmasses_ionquant.txt` (modification masses file for IonQuant)
- **`experiment_annotation.tsv`** (experiment annotation file mapping TMT channels to samples)
- `*_annotation.txt` (plex-specific annotation files mapping TMT channels to sample names)
- `fragpipe.workflow` (updated FragPipe workflow configuration file)
- `fragpipe-files.fp-manifest` (FragPipe files manifest)
- `fragpipe.job` (FragPipe job configuration file)
- `log_*.txt` (FragPipe execution log file with timestamp)
- `sdrf.tsv` (Sample and Data Relationship Format file)

<br>

### 3b. Check Spectral Files Centroid Status

```bash
java -Xmx55G -cp /fragpipe_bin/fragpipe-23.1/fragpipe-23.1/lib/fragpipe-23.1.jar:/fragpipe_bin/fragpipe-23.1/fragpipe-23.1/tools/batmass-io-1.35.4.jar org.nesvilab.fragpipe.util.CheckCentroid *.mzML 31
```
<!-- CLI mode (backup) - same command, no changes needed for headless mode -->

**Parameter Definitions:**

- `-Xmx55G` – Java memory limit (e.g., `-Xmx55G` for 55 GB RAM)
- `-cp` – Java classpath to FragPipe and BatMass libraries
- `org.nesvilab.fragpipe.util.CheckCentroid` – CheckCentroid main class
- `*.mzML` – input mzML file(s) to check
- `31` – number of CPU threads to use

**Input Data:**

- `*.mzML` (input mass spectrometry raw data in mzML format)

**Output Data:**

- (No output files; checks if mzML files are centroided or profile mode; FragPipe exits if files are not centroided)

<br>

### 3c. Initialize Workspace

```bash
philosopher-v5.1.2 workspace --clean --nocheck
philosopher-v5.1.2 workspace --init --nocheck --temp /tmp/temp_directory
```
<!-- ```bash
philosopher workspace --clean --nocheck
philosopher workspace --init --nocheck --temp /tmp/temp_directory
``` -->

**Parameter Definitions:**

- `workspace` – Philosopher subcommand for managing workspace
- `--clean` – removes any existing workspace files
- `--init` – initializes a new Philosopher workspace
- `--nocheck` – skips workspace validation checks
- `--temp` – specifies temporary directory for workspace initialization

**Output Data:**

- .meta/ (Philosopher workspace metadata directory containing binary database files)

<br>

### 3d. MSFragger Database Search

```bash
java -jar -Dfile.encoding=UTF-8 -Xmx55G MSFragger-4.3.jar fragger.params sample1.mzML sample2.mzML
```
<!-- CLI mode (backup) - same command, no changes needed for headless mode -->

**Parameter Definitions:**

- `-jar` – executes JAR file
- `-Dfile.encoding=UTF-8` – sets file encoding to UTF-8
- `-Xmx55G` – Java memory limit (e.g., `-Xmx55G` for 55 GB RAM)
- `MSFragger-4.3.jar` – MSFragger JAR file
- `fragger.params` – MSFragger parameter configuration file
- `*.mzML` – multiple mzML files provided as individual paths separated by spaces

**Input Data:**

- `fragger.params` (MSFragger parameter configuration file, output from [Step 3a](#3a-launch-fragpipe))
- `*.mzML` (input mass spectrometry raw data in mzML format)
- `*-decoys-reviewed-contam-*.fas` (proteome FASTA database with decoys and contaminants, output from [Step 2](#2-create-proteome-fasta-database))

**Output Data:**

- **\*.pepXML** (peptide-spectrum matches from the MSFragger database search)
- **\*.pin** (peptide-spectrum matches from the MSFragger database search in Percolator input format (PIN) for statistical validation)
- **\*.pepindex** (peptide index files for the FASTA database)
- **\*.tsv** (MSFragger results in tab-separated format)

<br>

### 3e. MSBooster Deep Learning Feature Addition

```bash
java -Djava.awt.headless=true -Xmx55G -cp MSBooster-1.3.17.jar:batmass-io-1.35.4.jar mainsteps.MainClass --paramsList msbooster_params.txt
```
<!-- CLI mode (backup):
```bash
java -Xmx55G -cp MSBooster-1.3.17.jar:batmass-io-1.35.4.jar mainsteps.MainClass --paramsList msbooster_params.txt
``` -->

**Parameter Definitions:**

- `-Djava.awt.headless=true` – runs in headless mode (no GUI)
- `-Xmx55G` – Java memory limit (e.g., `-Xmx55G` for 55 GB RAM)
- `-cp` – Java classpath to MSBooster and BatMass libraries
- `mainsteps.MainClass` – MSBooster main class
- `--paramsList` – path to MSBooster parameter configuration file

**Input Data:**

- `msbooster_params.txt` (MSBooster parameter configuration file, output from [Step 3a](#3a-launch-fragpipe))
- `*.pin` (Percolator input files from MSFragger, output from [Step 3d](#3d-msfragger-database-search))
- `*.mzML` (original mass spectrometry raw data in mzML format)

**Output Data:**

- \*_edited.pin (Percolator input files with added deep learning features from MSBooster: unweighted spectral entropy, weighted spectral entropy, hypergeometric probability, intersection, predicted RT real units, and delta RT LOESS)
- spectraRT_full.tsv (full spectra retention time data)
- spectraRT.predicted.bin (binary file containing predicted spectra, retention times, and ion mobilities from DIA-NN)
- spectraRT.tsv (spectra retention time data)
- MSBooster_plots/ (directory containing diagnostic plots: RT_calibration_curves/ with retention time calibration plots per pin file (up to top 5000 PSMs (peptide-spectrum matches) used for calibration), IM_calibration_curves/ with ion mobility calibration plots per charge state (up to top 1000 PSMs (peptide-spectrum matches) used for calibration; if IM features enabled), and score_histograms/ with overlayed histograms of all target and decoy PSMs for each pin file for all deep learning features (some features plotted on log scale for visualization, but original values are used in pin files))

<br>

### 3f. Percolator PSM Rescoring and Statistical Validation

```bash
/fragpipe_bin/fragpipe-23.1/fragpipe-23.1/tools/percolator_3_7_1/linux/percolator \
  --only-psms \
  --no-terminate \
  --post-processing-tdc \
  --num-threads 31 \
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
  --num-threads 31 \
  --results-psms *_percolator_target_psms.tsv \
  --decoy-results-psms *_percolator_decoy_psms.tsv \
  --protein-decoy-pattern rev_ \
  *_edited.pin
``` -->

**Parameter Definitions:**

- `--only-psms` – do not remove redundant peptides, keep PSMs, exclude peptide level probabilities
- `--no-terminate` – do not terminate execution when encountering issues with SVM inputs or results (default: false)
- `--post-processing-tdc` – replace mix-max method with target-decoy competition for assigning q-values and PEPs. For input PSMs from separate target/decoy searches, Percolator SVM scores eliminate lower-scoring target or decoy PSMs for each scan+expMass combination. Automatically enabled for concatenated searches
- `--num-threads` – number of CPU threads to use
- `--results-psms` – output file path for target PSM results
- `--decoy-results-psms` – output file path for decoy PSM results
- `--protein-decoy-pattern` – text pattern used to identify decoy proteins in the database
- `*_edited.pin` – input Percolator input files with MSBooster features

**Input Data:**

- `*_edited.pin` (Percolator input files with MSBooster features, output from [Step 3e](#3e-msbooster-deep-learning-feature-addition))

**Output Data:**

- *_percolator_target_psms.tsv (Percolator target PSM results in TSV format)
- *_percolator_decoy_psms.tsv (Percolator decoy PSM results in TSV format)

<br>

#### 3f.1. Add Percolator Validation Information to pepXML

```bash
java -cp /fragpipe_bin/fragpipe-23.1/fragpipe-23.1/lib/* \
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
- `0.5` – minimum probability threshold (1 - PEP); filters out PSMs with PEP > 0.5 (default: 0.5)
- `*.mzML` – original mzML file path

**Input Data:**

- `*.pin` (original Percolator input files from MSFragger, output from [Step 3d](#3d-msfragger-database-search))
- `*_percolator_target_psms.tsv` (Percolator target PSM results, output from [Step 3f](#3f-percolator-psm-rescoring-and-statistical-validation))
- `*_percolator_decoy_psms.tsv` (Percolator decoy PSM results, output from [Step 3f](#3f-percolator-psm-rescoring-and-statistical-validation))
- `*.mzML` (original mass spectrometry raw data in mzML format)

**Output Data:**

- interact-*.pep.xml (peptide-spectrum matches with validation information generated by Percolator)

<br>

### 3g. ProteinProphet Protein Inference and Statistical Validation

```bash
philosopher-v5.1.2 proteinprophet --maxppmdiff 2000000 --minprob 0.5 --output combined filelist_proteinprophet.txt
```
<!-- CLI mode (backup):
```bash
philosopher proteinprophet --maxppmdiff 2000000 --minprob 0.5 --output combined filelist_proteinprophet.txt
``` -->

**Parameter Definitions:**

- `proteinprophet` – run ProteinProphet to generate probabilities for protein identifications based on MS/MS data
- `--maxppmdiff 2000000` – maximum peptide mass difference in ppm (default: 20)
- `--minprob 0.5` – PeptideProphet minimum probability threshold (default: 0.05)
- `--output combined` – output file name (default: "interact.prot.xml"); results in `combined.prot.xml`
- `filelist_proteinprophet.txt` – list of interact.pep.xml files to be passed to ProteinProphet

**Input Data:**

- `filelist_proteinprophet.txt` (list of interact.pep.xml files to be passed to ProteinProphet, output from [Step 3a](#3a-launch-fragpipe))
- `interact-*.pep.xml` (pepXML files listed in filelist_proteinprophet.txt, output from [Step 3f.1](#3f1-add-percolator-validation-information-to-pepxml))

**Output Data:**

- combined.prot.xml (protein identifications with validation information generated by ProteinProphet via Philosopher)

<br>

### 3h. Database Annotation

```bash
philosopher-v5.1.2 database --annotate *.fas --prefix rev_
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

- `*-decoys-reviewed-contam-*.fas` (proteome FASTA database with decoys and contaminants, output from [Step 2](#2-create-proteome-fasta-database))

**Output Data:**

- .meta/ (Philosopher workspace metadata directory containing binary database files)

<br>

### 3i. Filter Results by FDR

```bash
# First sample (initializes database annotation)
philosopher-v5.1.2 filter \
  --sequential \
  --prot 0.01 \
  --picked \
  --tag rev_ \
  --pepxml sample_directory \
  --protxml combined.prot.xml \
  --razor

# Subsequent samples (reuse database annotation from first sample)
philosopher-v5.1.2 filter \
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

- `interact-*.pep.xml` (peptide-spectrum matches with validation information generated by Percolator, output from [Step 3f.1](#3f1-add-percolator-validation-information-to-pepxml))
- `combined.prot.xml` (protein identifications with validation information generated by ProteinProphet via Philosopher, output from [Step 3g](#3g-proteinprophet-protein-inference-and-statistical-validation))
- .meta/ (Philosopher workspace metadata, output from [Step 3h](#3h-database-annotation))

**Output Data:**

- filter.log (Philosopher filter execution log file)
- Filtered data stored in Philosopher workspace as binary files (db.bin, ion.bin, pep.bin, pro.bin, protxml.bin, psm.bin, razor.bin, etc.; filtered PSM, peptide, and protein data ready for report generation)

<br>

### 3j. Generate Reports

```bash
philosopher-v5.1.2 report
```
<!-- CLI mode (backup):
```bash
philosopher report
``` -->

**Input Data:**

- Philosopher workspace containing filtered data (output from [Step 3i](#3i-filter-results-by-fdr))

**Output Data:**

- protein.fas (FASTA file containing FDR-filtered protein sequences identified in the analysis, generated by Philosopher)
- protein.tsv (protein report)
- peptide.tsv (peptide report)
- psm.tsv (PSM report)
- ion.tsv (ion report)

<br>

### 3k. IonQuant TMT Reporter Ion Extraction

```bash
# First pass: MS1 quantification
java -Djava.awt.headless=true -Xmx55G \
  -Dlibs.bruker.dir=tools/ext/bruker \
  -Dlibs.thermo.dir=tools/ext/thermo \
  -cp jfreechart-1.5.3.jar:IonQuant-1.11.11.jar \
  ionquant.IonQuant \
  --threads 31 \
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
java -Djava.awt.headless=true -Xmx55G \
  -Dlibs.bruker.dir=tools/ext/bruker \
  -Dlibs.thermo.dir=tools/ext/thermo \
  -cp jfreechart-1.5.3.jar:IonQuant-1.11.11.jar \
  ionquant.IonQuant \
  --threads 31 \
  --perform-ms1quant 0 \
  --perform-isoquant 1 \
  --isotol 20.0 \
  --isolevel 2 \
  --isotype TMT-10 \
  --ionmobility 0 \
  --site-reports 0 \
  --msstats 0 \
  --annotation sample1/psm.tsv=sample1/sample1_annotation.txt \
  --annotation sample2/psm.tsv=sample2/sample2_annotation.txt \
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
- `-Xmx55G` – Java memory limit (e.g., `-Xmx55G` for 55 GB RAM)
- `-Dlibs.bruker.dir` – directory for Bruker libraries
- `-Dlibs.thermo.dir` – directory for Thermo libraries
- `-cp` – Java classpath to jfreechart and IonQuant JAR files
- `ionquant.IonQuant` – IonQuant main class
- `--threads` – number of CPU threads to use (0 = all logical cores)
- `--perform-ms1quant 1` – perform MS1 quantification (first pass; 0 = no, 1 = yes)
- `--perform-isoquant 1` – perform isobaric labeling quantification (second pass; 0 = no, 1 = yes)
- `--isotol 20.0` – MS2 tolerance in ppm for isobaric quantification (default: 10)
- `--isolevel 2` – isobaric quantification level (2 = MS2, 3 = MS3, 4 = ZOOM-HR; default: 2)
- `--isotype tmt10` / `--isotype TMT-10` – isobaric quantification type (first pass: `tmt10` for MS1 quantification; second pass: `TMT-10` for isobaric quantification)
- `--ionmobility 0` – data has ion mobility information (0 = no, 1 = yes; default: 0)
- `--site-reports 0` – generate site reports (0 = no, 1 = yes; requires modification localization columns in psm.tsv; default: 1)
- `--msstats 0` – generate MSstats input files (0 = no, 1 = yes; default: 0)
- `--annotation` – annotation file for isobaric quantification (format: `psm.tsv=annotation.txt`; can specify multiple; used in second pass only)
- `--minexps 1` – minimum experiments in picking an ion for quantifying proteins (only for intensity, not MaxLFQ; default: 1)
- `--mbr 0` – perform match-between-runs (0 = no, 1 = yes; default: 0)
- `--maxlfq 0` – calculate MaxLFQ intensity (0 = no, 1 = yes; default: 1)
- `--requantify 0` – re-quantify unidentified features based on identified features (0 = no, 1 = yes; default: 1)
- `--mztol 10` – MS1 tolerance in ppm (default: 10.0)
- `--imtol 0.05` – 1/K0 tolerance (default: 0.05)
- `--rttol 1` – retention time tolerance in minutes (default: 0.4)
- `--normalization 0` – normalize intensities across all runs (0 = no, 1 = yes; default: 1)
- `--minisotopes 1` – minimum isotopes required in feature extraction (default: 2)
- `--minscans 1` – minimum MS1 scans required in feature extraction (default: 3)
- `--writeindex 0` – write indexed file on disk for further usage (0 = no, 1 = yes; default: 0)
- `--tp 0` – number of ions used in quantifying each protein (0 = use all ions; only for intensity, not MaxLFQ; default: 0)
- `--minfreq 0` – minimum required frequency of an ion being selected for protein quantification (only for intensity, not MaxLFQ; default: 0)
- `--minions 1` – minimum ions required for quantifying proteins (only for MaxLFQ intensity; default: 1)
- `--locprob 0` – localization probability threshold (default: 0)
- `--uniqueness 0` – peptide-protein uniqueness (0 = unique+razor, 1 = unique only, 2 = all; default: 0)
- `--multidir .` – output directory for multi-experimental results (optional)
- `--filelist` – file containing flags (tab-delimited file with `--psm` entries pointing to `psm.tsv` files and `--specdir` entry pointing to the directory containing mzML files)
- `--modlist` – file listing modification masses (used to remove mass discrepancy due to rounding errors)

**Input Data:**

- `filelist_ionquant.txt` (file list for IonQuant containing `--psm` entries pointing to `psm.tsv` files and `--specdir` entry pointing to mzML directory, output from [Step 3a](#3a-launch-fragpipe))
- `modmasses_ionquant.txt` (modification masses file for IonQuant, output from [Step 3a](#3a-launch-fragpipe))
- `protein.tsv` (protein report, output from [Step 3j](#3j-generate-reports))
- `peptide.tsv` (peptide report, output from [Step 3j](#3j-generate-reports))
- `psm.tsv` (PSM report, output from [Step 3j](#3j-generate-reports))
- `ion.tsv` (ion report, output from [Step 3j](#3j-generate-reports))
- `*_annotation.txt` (plex-specific annotation files mapping TMT channels to sample names, used in second pass only, output from [Step 3a](#3a-launch-fragpipe))
- `*.mzML` (original mass spectrometry raw data in mzML format; accessed via `--specdir` parameter specified in `filelist_ionquant.txt` to extract TMT reporter ion intensities)

**Output Data:**

- protein.tsv (protein report with TMT reporter ion intensities and additional data added from IonQuant)
- peptide.tsv (peptide report with TMT reporter ion intensities and additional data added from IonQuant)
- ion.tsv (ion report with TMT reporter ion intensities and additional data added from IonQuant)
- psm.tsv (PSM report with TMT reporter ion intensities (individual TMT channel columns) and additional data added from IonQuant)
- combined_protein.tsv (combined protein report with TMT reporter ion intensities across all samples)
- combined_peptide.tsv (combined peptide report with TMT reporter ion intensities and additional data across all samples)
- combined_ion.tsv (combined ion report with TMT reporter ion intensities and additional data across all samples)
- combined_modified_peptide.tsv (combined modified peptide report with TMT reporter ion intensities and additional data across all samples)
- reprint.int.tsv (input file for the Resource for Evaluation of Protein Interaction Networks (REPRINT) containing protein intensities, generated by Philosopher)
- reprint.spc.tsv (input file for the Resource for Evaluation of Protein Interaction Networks (REPRINT) containing protein spectral counts, generated by Philosopher)

<br>

### 3l. TMTIntegrator TMT Quantification

```bash
java -Xmx55G -jar TMT-Integrator-6.1.1.jar \
  tmt-integrator-conf.yml \
  sample1/psm.tsv \
  sample2/psm.tsv
```

**Parameter Definitions:**

- `-Xmx55G` – Java memory limit (e.g., `-Xmx55G` for 55 GB RAM)
- `-jar` – executes JAR file
- `TMT-Integrator-6.1.1.jar` – TMTIntegrator JAR file
- `tmt-integrator-conf.yml` – TMTIntegrator configuration file
- `sample*/psm.tsv` – PSM files with TMT reporter ion intensities

**Input Data:**

- `tmt-integrator-conf.yml` (TMTIntegrator configuration file, output from [Step 3a](#3a-launch-fragpipe))
- `psm.tsv` (PSM reports with TMT reporter ion intensities, output from [Step 3k](#3k-ionquant-tmt-reporter-ion-extraction))

**Output Data:**

- tmt-report/abundance_protein_MD.tsv (protein-level TMT abundance values with median-centering normalization, log2-transformed)
- tmt-report/abundance_peptide_MD.tsv (peptide-level TMT abundance values with median-centering normalization, log2-transformed)
- tmt-report/abundance_gene_MD.tsv (gene-level TMT abundance values with median-centering normalization, log2-transformed)
- tmt-report/ratio_protein_MD.tsv (protein-level TMT ratio values; ratios are channel abundance divided by reference channel abundance, log2-transformed)
- tmt-report/ratio_peptide_MD.tsv (peptide-level TMT ratio values; ratios are channel abundance divided by reference channel abundance, log2-transformed)
- tmt-report/ratio_gene_MD.tsv (gene-level TMT ratio values; ratios are channel abundance divided by reference channel abundance, log2-transformed)

<br>

---

## 4. Compile FragPipe QC Reports

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

- `psm.tsv` (plex-specific PSM reports, output from [Step 3k](#3k-ionquant-tmt-reporter-ion-extraction))
- `ion.tsv` (plex-specific ion reports, output from [Step 3k](#3k-ionquant-tmt-reporter-ion-extraction))
- `combined_protein.tsv` (combined protein report, output from [Step 3k](#3k-ionquant-tmt-reporter-ion-extraction))
- `combined_peptide.tsv` (combined peptide report, output from [Step 3k](#3k-ionquant-tmt-reporter-ion-extraction))
- `combined_ion.tsv` (combined ion report, output from [Step 3k](#3k-ionquant-tmt-reporter-ion-extraction))
- `*.workflow` (FragPipe workflow file, output from [Step 3a](#3a-launch-fragpipe))
- `fragger.params` (MSFragger parameters file, output from [Step 3a](#3a-launch-fragpipe))

**Output Data:**

- **multiqc_GLProteomics.html** (MultiQC output html summary)
- **multiqc_GLProteomics_data.zip** (zipped directory containing MultiQC output data with cleaned paths)

<br>
