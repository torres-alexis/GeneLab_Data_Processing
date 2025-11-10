# GeneLab bioinformatics processing pipeline for Mass Spectrometry-based Proteomics Data

> **This page holds an overview and instructions for how GeneLab processes mass spectrometry-based proteomics data. Exact processing commands, GL-DPPD-X version used, and processed data output files for specific datasets are provided in the [Open Science Data Repository (OSDR)](https://osdr.nasa.gov/bio/repo/).**  

---

**Date:** November X, 2025  
**Revision:** A  
**Document Number:** GL-DPPD-X-A  

**Submitted by:**  
Alexis Torres (GeneLab Data Processing Team)  

**Approved by:**  
X (X)
---

# Table of contents  

- [**Software used**](#software-used)
- [**General processing overview with example commands**](#general-processing-overview-with-example-commands)
  - [**1. Raw Data QC**](#1-raw-data-qc)
    - [1a. Convert Raw Data to mzML](#1a-convert-raw-data-to-mzml)
    - [1b. Raw Data QC](#1b-raw-data-qc)
    - [1c. Compile Raw Data QC](#1c-compile-raw-data-qc)
  - [**2. Create Proteome FASTA Database**](#2-create-proteome-fasta-database)
    - [2a. Download Proteome from UniProt](#2a-download-proteome-from-uniprot)
    - [2b. Add Decoys and Contaminants to FASTA](#2b-add-decoys-and-contaminants-to-fasta)
  - [**3. Run MSFragger Database Search**](#3-run-msfragger-database-search)
  - [**4. Mass Recalibration with Crystal-C (Open Searches Only)**](#4-mass-recalibration-with-crystal-c-open-searches-only)
    - [4a. Run Crystal-C](#4a-run-crystal-c)
  - [**5. Peptide and Protein Validation with Philosopher**](#5-peptide-and-protein-validation-with-philosopher)
    - [5a. Initialize Philosopher Workspace](#5a-initialize-philosopher-workspace)
    - [5b. Annotate Database](#5b-annotate-database)
    - [5c. Validate Peptides with PeptideProphet](#5c-validate-peptides-with-peptideprophet)
    - [5d. Infer Proteins with ProteinProphet](#5d-infer-proteins-with-proteinprophet)
    - [5e. Filter Results by FDR](#5e-filter-results-by-fdr)
    - [5f. Generate Reports](#5f-generate-reports)
  - [**6. Label-Free Quantification with IonQuant**](#6-label-free-quantification-with-ionquant)
    

---

# Software used  

|Program|Version|Relevant Links|
|:------|:------:|:-------------|
|dp_tools|1.3.8|[https://github.com/torres-alexis/dp_tools](https://github.com/torres-alexis/dp_tools)|
|msconvert|3.0.9992 (ProteoWizard)|[http://proteowizard.sourceforge.net/](http://proteowizard.sourceforge.net/)|
|rawBeans|1.6.4|[https://github.com/torres-alexis/rawBeans](https://github.com/torres-alexis/rawBeans)|
|MSFragger|4.2|[http://msfragger-upgrader.nesvilab.org/upgrader/](http://msfragger-upgrader.nesvilab.org/upgrader/)|
|Crystal-C|TBD|[https://github.com/Nesvilab/Crystal-C/releases/latest](https://github.com/Nesvilab/Crystal-C/releases/latest)|
|Philosopher|5.1.2|[https://github.com/Nesvilab/philosopher/releases/latest](https://github.com/Nesvilab/philosopher/releases/latest)|
|IonQuant|1.9.9|[https://github.com/Nesvilab/IonQuant/releases/latest](https://github.com/Nesvilab/IonQuant/releases/latest)|


---

# General processing overview with example commands  

> Exact processing commands and output files listed in **bold** below are included with each mass spectrometry-based proteomics processed dataset in the [Open Science Data Repository (OSDR)](https://osdr.nasa.gov/bio/repo/). 

---

### 1. Raw Data QC  

```bash
create-qc-report.py \
  --input sample1.mzML sample2.mzML sample3.mzML \
  --output-dir . \
  --batch \
  --cores NumberOfThreads
```

**Parameter Definitions:**

- `--input` – one or multiple mzML files provided as individual paths separated by spaces
- `--output-dir` – the output directory to store results
- `--batch` – process files in batch mode
- `--cores` – number of CPU cores to use for processing

**Input Data:**

- `*.mzML` (input mass spectrometry raw data in mzML format)

**Output Data:**

- GLProteomics-rawbeans_qc-report.html(rawBeans QC output html summary)
- resources/ (directory containing mass deviation matrices and supporting files for the html summary)

<br>

---

### 2. Download Reference Proteome, Add Decoys and Contaminants to FASTA

```bash
  # Download proteome using philosopher database command with --id
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

## 3. Run MSFragger Database Search

```bash
  # MSFragger run for batch of mzML files
  java -XmxN -jar MSFragger.jar \
    /path/to/MSFragger/params/file \
    sample1.mzML sample2.mzML sample3.mzML
```

**Parameter Definitions:**

- `-XmxN` – Java memory limit (e.g., `-Xmx32G` for 32 GB RAM)

**Input Data:**

- `/path/to/MSFragger/params/file` – path to MSFragger params file
- `*.mzML` (input mass spectrometry raw data in mzML format, output from [Step 1](#1-raw-data-qc))
- `*-decoys-reviewed-contam-*.fas` (proteome FASTA database with decoys and contaminants, output from [Step 2](#2-download-reference-proteome-add-decoys-and-contaminants-to-fasta))

**Output Data:**

- **\*.pepXML** (peptide-spectrum matches in pepXML format; used as input for Crystal-C in open searches and for Philosopher validation)
- \*.pin (Percolator input format for statistical validation)
- \*.pepindex (peptide index files for the FASTA database)
- \*.tsv (MSFragger results in tab-separated format; not used by Crystal-C, moved for organizational purposes)

<br>

---

## 4. Mass Recalibration with Crystal-C (Open Search Only)

<br>

### 4a. Run Crystal-C

> Note: Crystal-C performs mass recalibration for open search results. This step is only required when using open search mode in MSFragger. Crystal-C recalibrates the mass measurements in the pepXML files based on the original spectral data.

```bash
  # Crystal-C mass recalibration for open searches (for mzML files)
  java -Xmx53G -cp "CrystalC-1.2.1.jar" crystalc.Run \
    crystalc.params \
    *.pepXML
```

**Parameter Definitions:**

- `-XmxN` – Java memory limit (e.g., `-Xmx53G` for 53 GB RAM)
- `-cp "CrystalC-1.2.1.jar"` – classpath to Crystal-C JAR file
- `crystalc.Run` – Crystal-C main class
- `crystalc.params` – Crystal-C parameter configuration file
- `*.pepXML` – input pepXML files from MSFragger open search

**Input Data:**

- `*.pepXML` (peptide search results from MSFragger open search, output from [Step 3](#3-run-msfragger-database-search))
- `*.mzML` (original mass spectrometry raw data in mzML format; location specified via `raw_file_location` parameter in `crystalc.params`)
- Crystal-C parameter configuration file (`crystalc.params`)
- `*-decoys-reviewed-contam-*.fas` (proteome FASTA database, same as used for MSFragger search; specified via `fasta` parameter in `crystalc.params`, output from [Step 2](#2-download-reference-proteome-add-decoys-and-contaminants-to-fasta))

**Output Data:**

- **\*.pepXML** (mass-recalibrated peptide search results)

<br>

---

## 5. Peptide and Protein Validation with Philosopher

<br>

### 5a. Initialize Philosopher Workspace

```bash
  # Clean and initialize Philosopher workspace
  philosopher workspace --clean
  philosopher workspace --init
```

**Parameter Definitions:**

- `workspace --clean` – removes any existing workspace files
- `workspace --init` – initializes a new Philosopher workspace

<br>

### 5b. Annotate Database

```bash
  # Annotate database with decoy prefix
  philosopher database --annotate database.fas --prefix rev_
```

**Parameter Definitions:**

- `--annotate` – path to FASTA database file
- `--prefix` – decoy prefix used in the database (default: `rev_`)

**Input Data:**

- `*-decoys-reviewed-contam-*.fas` (proteome FASTA database with decoys and contaminants, output from [Step 2](#2-download-reference-proteome-add-decoys-and-contaminants-to-fasta))

**Output Data:**

- .meta/ (Philosopher workspace metadata directory containing binary database files; internal workspace state required for subsequent Philosopher steps)

<br>

### 5c. Validate Peptides with PeptideProphet

```bash
  # Closed search
  philosopher peptideprophet \
    --nonparam \
    --expectscore \
    --decoyprobs \
    --ppm \
    --accmass \
    --decoy rev_ \
    --database database.fas \
    ./*.pepXML

  # Open search if you ran Crystal-C (uses *_c.pepXML files)
  philosopher peptideprophet \
    --nonparam \
    --expectscore \
    --decoyprobs \
    --masswidth 1000.0 \
    --clevel -2 \
    --decoy rev_ \
    --combine \
    --database database.fas \
    ./*_c.pepXML

  # Open search if you did NOT run Crystal-C (uses *.pepXML files)
  philosopher peptideprophet \
    --nonparam \
    --expectscore \
    --decoyprobs \
    --masswidth 1000.0 \
    --clevel -2 \
    --decoy rev_ \
    --combine \
    --database database.fas \
    ./*.pepXML

  # Non-specific closed search
  philosopher peptideprophet \
    --nonparam \
    --expectscore \
    --decoyprobs \
    --ppm \
    --accmass \
    --nontt \
    --decoy rev_ \
    --database database.fas \
    ./*.pepXML
```

**Parameter Definitions:**

- `--nonparam` – use semi-parametric modeling (must be used with `--decoy`)
- `--expectscore` – use expectation value as the only contributor to the f-value for modeling
- `--decoyprobs` – compute possible non-zero probabilities for decoy entries on the last iteration
- `--ppm` – use parts-per-million mass error instead of Dalton for mass modeling
- `--accmass` – use accurate mass model binning
- `--masswidth` – model mass width in Da (default: 5.0)
- `--clevel` – conservative level in neg_stdev from the neg_mean; low numbers are less conservative, high numbers are more conservative
- `--combine` – combine the results from PeptideProphet into a single result file
- `--nontt` – disable NTT enzymatic termini model
- `--decoy` – semi-supervised mode, protein name prefix to identify decoy entries (default: `rev_`)
- `--database` – path to FASTA database file
- `*.pepXML` or `*_c.pepXML` – input pepXML files from MSFragger (or Crystal-C for open searches)

**Input Data:**

- `*.pepXML` (peptide search results from MSFragger, output from [Step 3](#3-run-msfragger-database-search)) or `*.pepXML` (mass-recalibrated peptide search results, output from [Step 4](#4-mass-recalibration-with-crystal-c-open-search-only))
- `*-decoys-reviewed-contam-*.fas` (proteome FASTA database, output from [Step 2](#2-download-reference-proteome-add-decoys-and-contaminants-to-fasta))

**Output Data:**

- **interact-*.pep.xml** (PeptideProphet-validated peptide identifications; default output prefix is "interact", creates one file per input pepXML file) or **interact.pep.xml** (when using `--combine` flag, creates a single combined file)

<br>

### 5d. Infer Proteins with ProteinProphet

```bash
  # Run ProteinProphet to infer proteins from peptides
  philosopher proteinprophet \
    --maxppmdiff 2000000 \
    --output combined \
    ./*.pep.xml
```

**Parameter Definitions:**

- `--maxppmdiff` – maximum parts-per-million difference for protein grouping (default: 2000000)
- `--output` – output filename prefix (e.g., `combined` produces `combined.prot.xml`)
- `interact-*.pep.xml` or `interact.pep.xml` – input PeptideProphet-validated pepXML files

**Input Data:**

- `interact-*.pep.xml` or `interact.pep.xml` (PeptideProphet-validated peptide identifications, output from [Step 5c](#5c-validate-peptides-with-peptideprophet))

**Output Data:**

- **combined.prot.xml** (ProteinProphet-inferred protein identifications)

<br>

### 5e. Filter Results by FDR

```bash
  # Filter results by FDR (closed or non-specific closed search)
  philosopher filter \
    --sequential \
    --razor \
    --mapmods \
    --tag rev_ \
    --pepxml ./ \
    --protxml ./combined.prot.xml

  # Filter results by FDR (open search)
  philosopher filter \
    --sequential \
    --razor \
    --mapmods \
    --tag rev_ \
    --pepxml ./interact.pep.xml \
    --protxml ./combined.prot.xml
```

**Parameter Definitions:**

- `--sequential` – use sequential FDR filtering
- `--razor` – use razor peptide assignment (assign peptides to minimal set of proteins)
- `--mapmods` – map modifications
- `--tag` – decoy prefix (e.g., `rev_`)
- `--pepxml` – directory containing pepXML files (closed/non-specific search) or path to `interact.pep.xml` file (open search)
- `--protxml` – path to ProteinProphet protXML file

**Input Data:**

- `interact-*.pep.xml` or `interact.pep.xml` (PeptideProphet-validated peptide identifications, output from [Step 5c](#5c-validate-peptides-with-peptideprophet))
- `combined.prot.xml` (ProteinProphet-inferred protein identifications, output from [Step 5d](#5d-infer-proteins-with-proteinprophet))

**Output Data:**

- **interact-*.pep.xml** or **interact.pep.xml** (FDR-filtered peptide identifications; modifies existing `interact-*.pep.xml` or `interact.pep.xml` files in place by removing identifications that do not meet FDR thresholds)
- **combined.prot.xml** (FDR-filtered protein identifications; modifies existing `combined.prot.xml` file in place by removing proteins that do not meet FDR thresholds)

> Note: Filter modifies the input files in place. If all identifications already meet the FDR thresholds, the files may remain unchanged. The filtered files are the same filenames as the inputs.

<br>

### 5f. Generate Reports

```bash
  # Generate summary reports
  philosopher report
```

**Output Data:**

- **psm.tsv** (peptide-spectrum match results with probabilities, q-values, and protein assignments)
- **peptide.tsv** (peptide-level identification results with probabilities and q-values; aggregated across all samples)
- **protein.tsv** (protein-level identification results with probabilities, q-values, coverage, and spectral counts; aggregated across all samples)
- **ion.tsv** (ion-level quantification data for label-free quantification)
- **protein.fas** (FASTA file containing sequences of identified proteins)

> Note: When PeptideProphet uses the `--combine` flag (open searches), Philosopher report generates `combined_peptide.tsv` and `combined_protein.tsv` instead of `peptide.tsv` and `protein.tsv`. These combined reports contain per-sample columns (e.g., `<sample> Total Spectral Count`, `<sample> Unique Intensity`) rather than aggregated values. For closed searches without `--combine`, the standard `peptide.tsv` and `protein.tsv` files are generated with aggregated values across all samples.

<br>

---

## 6. Label-Free Quantification with IonQuant

<br>

```bash
java -Xmx32G -jar IonQuant.jar \
  --specdir . \
  --psm psm.tsv \
  --perform-ms1quant 1 \
  --site-reports 1 \
  --msstats 1 \
  --multidir . \
  --normalization 1 \
  --maxlfq 1 \
  --threads 8
```

**Input Data:**

- `psm.tsv` (peptide-spectrum match results from Philosopher, output from [Step 5f](#5f-generate-reports))
- `protein.tsv` (protein identification results from Philosopher, output from [Step 5f](#5f-generate-reports))
- `*.mzML` (original mass spectrometry raw data in mzML format; located in directories specified via `--specdir`)

**Output Data:**

IonQuant performs label-free quantification and updates the Philosopher report files with intensity measurements. IonQuant generates the following outputs:

- protein.tsv (updated protein report with quantification columns added: `Total Intensity`, `Unique Intensity`, `Razor Intensity`, `Coverage`, and other quantification metrics; modifies existing `protein.tsv` file from Philosopher)
- peptide.tsv (peptide-level quantification report with `Intensity` and `Spectral Count` columns)
- ion.tsv (ion-level quantification report with detailed intensity measurements, retention times, and charge states for each identified ion)
- psm.tsv (updated PSM report with quantification data; modifies existing `psm.tsv` file from Philosopher)
- combined_protein.tsv (combined protein quantification across all samples with per-sample intensity columns)
- combined_peptide.tsv (combined peptide quantification across all samples)
- combined_modified_peptide.tsv (combined modified peptide quantification across all samples)
- combined_ion.tsv (combined ion-level quantification across all samples)
- combined_site_*.tsv (site-specific modification reports, e.g., `combined_site_C_57.0214.tsv` for carbamidomethylation, `combined_site_M_15.9949.tsv` for oxidation)
- reprint.int.tsv (reprint intensity file)
- reprint.spc.tsv (reprint spectral count file)
- **msstats.csv** (formatted input file for the MSstats R package for downstream differential analysis; contains columns: `ProteinName`, `PeptideSequence`, `PrecursorCharge`, `FragmentIon`, `Condition`, `BioReplicate`, `Run`, `Intensity`; main input for downstream statistical analysis)
- **msstats_ptm.csv** (MSstats input file for PTM (post-translational modification) analysis; includes PTM site columns)

> Note: IonQuant automatically matches `psm.tsv` files to their corresponding mzML files by parsing the "Spectrum File" column in `psm.tsv`, which references the pepXML files that internally reference the mzML files. The `--specdir` directories are searched to locate the matching mzML files. IonQuant updates the `protein.tsv` file by adding quantification columns (e.g., `Total Intensity`, `Unique Intensity`, `Razor Intensity`) based on MS1 precursor intensities extracted from the mzML files.

<br>

---