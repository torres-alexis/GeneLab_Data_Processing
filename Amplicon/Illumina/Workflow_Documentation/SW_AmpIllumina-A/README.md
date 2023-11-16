# SW_AmpIllumina-A Workflow Information and Usage Instructions <!-- omit in toc -->


## General workflow info <!-- omit in toc -->
The current GeneLab Illumina amplicon sequencing data processing pipeline (AmpIllumina), [GL-DPPD-7104-A.md](../../Pipeline_GL-DPPD-7104_Versions/GL-DPPD-7104-A.md), is implemented as a [Snakemake](https://snakemake.readthedocs.io/en/stable/) workflow and utilizes [conda](https://docs.conda.io/en/latest/) environments to install/run all tools. This workflow (SW_AmpIllumina-A) is run using the command line interface (CLI) of any unix-based system. The workflow can be used even if you are unfamiliar with Snakemake and conda, but if you want to learn more about those, [this Snakemake tutorial](https://snakemake.readthedocs.io/en/stable/tutorial/tutorial.html) within [Snakemake's documentation](https://snakemake.readthedocs.io/en/stable/) is a good place to start for that, and an introduction to conda with installation help and links to other resources can be found [here at Happy Belly Bioinformatics](https://astrobiomike.github.io/unix/conda-intro).  

<br>

---

## Utilizing the workflow <!-- omit in toc -->

- [1. Install conda, mamba, and `genelab-utils` package](#1-install-conda-mamba-and-genelab-utils-package)
- [2. Download the workflow template files](#2-download-the-workflow-template-files)
- [3. Run the workflow using `run_workflow.py`](#3-run-the-workflow-using-run_workflowpy)
  - [3a. Approach 1: Run the workflow on GeneLab Illumina amplicon sequencing dataset with automatic retrieval of raw read files and metadata](#3a-approach-1-run-the-workflow-on-genelab-illumina-amplicon-sequencing-dataset-with-automatic-retrieval-of-raw-read-files-and-metadata)
  - [3b. Approach 2: Run the workflow using a local or user-created runsheet](#3b-approach-2-run-the-workflow-using-a-local-or-user-created-runsheet)
  - [3c. Approach 3: Run the workflow using pre-configured settings](#3c-approach-3-run-the-workflow-using-pre-configured-settings)


<br>

___

### 1. Install conda, mamba, and `genelab-utils` package
We recommend installing a Miniconda, Python3 version appropriate for your system, as exemplified in [the above link](https://astrobiomike.github.io/unix/conda-intro#getting-and-installing-conda).  

Once conda is installed on your system, we recommend installing [mamba](https://github.com/mamba-org/mamba#mamba), as it generally allows for much faster conda installations:

```bash
conda install -n base -c conda-forge mamba
```

> You can read a quick intro to mamba [here](https://astrobiomike.github.io/unix/conda-intro#bonus-mamba-no-5) if wanted.

Once mamba is installed, you can install the genelab-utils conda package in a new environment with the following command:

```bash
mamba create -n genelab-utils -c conda-forge -c bioconda -c defaults -c astrobiomike 'genelab-utils>=1.1.02'
```

The environment then needs to be activated:

```bash
conda activate genelab-utils
```

<br>

___

### 2. Download the workflow template files
All files required for utilizing the GeneLab workflow for processing Illumina amplicon sequencing data are in the [workflow_code](workflow_code) directory. To get a copy of latest SW_AmpIllumina-A version on to your system, the code can be downloaded as a zip file from the release page then unzipped after downloading by running the following commands:

```bash
wget https://github.com/nasa/GeneLab_Data_Processing/releases/download/SW_AmpIllumina-A_1.2.0/SW_AmpIllumina-A_1.2.0.zip

unzip SW_AmpIllumina-A_1.2.0.zip
```

This downloaded the workflow into a directory called `SW_AmpIllumina-*/`, with the workflow version number at the end.

<br>

___

### 3. Run the workflow using `run_workflow.py`

The `run_workflow.py` script in the `workflow_code/scripts` directory can be used to set up configuration files needed for the workflow and execute the workflow.

Navigate to the `workflow_code` directory to use this script.

<br>

___

#### 3a. Approach 1: Run the workflow on GeneLab Illumina amplicon sequencing dataset with automatic retrieval of raw read files and metadata

```bash
python ./scripts/run_workflow.py --OSD ### --run
```

This approach processes data from the NASA GeneLab Open Science Data Repository (OSDR). Upon execution, the command downloads the OSD/GLDS ISA file, which contains the necessary metadata for creating a runsheet. It will then create the runsheet, download the raw read files listed in the runsheet into `./raw_reads/` and prepare the necessary configuration files before executing the workflow using the default Snakemake run command.

<br>

___

### 3b. Approach 2: Run the workflow using a local or user-created runsheet

If you are working with a GeneLab Illumina amplicon sequencing dataset or a non-OSD dataset, the process for setting up and running the workflow using an existing runsheet is the same.

If you are working with a non-OSD dataset, you must manually create the runsheet for your dataset to run the workflow.

> Note: Specifications for creating a runsheet manually are described [here](examples/runsheet/README.md).

```bash
python ./scripts/run_workflow.py --runsheetPath </path/to/runsheet> --run
```

- For a GeneLab dataset, this command will download the raw read files listed in the runsheet to `./raw_reads/` and set up the configuration files before executing the workflow.
- For a non-OSD dataset, it will set up and run the workflow on local raw read data using a user-created runsheet.

### 3c. Approach 3: Run the workflow using pre-configured settings

Use this approach if you have already set up or manually modified `config.yaml` and other necessary files. To execute the Snakemake workflow with these pre-configured settings, simply use the `--run` argument.

```bash
python ./scripts/run_workflow.py --run
```

This command will start the Snakemake workflow using the existing `config.yaml` file, skipping the setup steps.

<br>

___

**Parameters for `run_workflow.py`:**


- **Set up for Approach 1**

  - `-o`, `--OSD` - Sets up the Snakemake workflow for a GeneLab OSD dataset. Pulls necessary metadata and read files. Acceptable Formats: `###,` `OSD-###`, `GLDS-###`.

- **Set up for Approach 2**

  - `-r`, `--runsheetPath` - Sets up the Snakemake workflow using a specified runsheet file. If the runsheet contains URLs, the script downloads the read files to `./raw_reads/`.

- **Optional parameters**

  - `-d`, `--outputDir` - Specifies the output directory for the workflow. Used only if `-o` or `-r` is specified. Default: `./workflow_output/`

<br>

- **Executing the workflow**
  - `-x`, `--run` - Command used to execute the snakemake workflow. Default: `snakemake --use-conda --conda-prefix ${CONDA_PREFIX}/envs -j 2 -p` 

___

**Below are instructions to prepare the configuration files and run the workflow manually if necessary.**

**Configuration files needed for the workflow**

`config.yaml` contains variables which the Snakemake workflow uses to perform its processing. `unique-sample-IDS.txt` contains the unique identifiers for each sample's read file(s).



You can modify the variables in the [config.yaml](workflow_code/config.yaml) file as needed. For example, you will have to provide a text file containing a single-column list of unique sample identifiers (see an example of how to set this up below). You will also need to indicate the paths to your input data (raw reads) and, if necessary, modify each variable to be consistent with the study you want to process. 

> Note: If you are unfamiliar with how to specify paths, one place you can learn more is [here](https://astrobiomike.github.io/unix/getting-started#the-unix-file-system-structure).  

**Example for how to create a single-column list of unique sample identifiers from your raw data file names**

For example, if you have paired-end read data for 2 samples located in `../Raw_Data/` relative to your workflow directory, that would look like this:

```bash
ls ../Raw_Data/
```

```bash
Sample-1_R1_raw.fastq.gz
Sample-1_R2_raw.fastq.gz
Sample-2_R1_raw.fastq.gz
Sample-2_R2_raw.fastq.gz
```

You would set up your `unique-sample-IDs.txt` file as follows:

```bash
cat unique-sample-IDs.txt
```

```
Sample-1
Sample-2
```

<br>

___


**Running the workflow from the command line**

While in the `workflow_code/` directory containing the `Snakefile`, `config.yaml`, and other workflow files that you downloaded in [step 2](#2-download-the-workflow-template-files), here is one example command of how to run the workflow from the command line:

```bash
snakemake --use-conda --conda-prefix ${CONDA_PREFIX}/envs -j 2 -p
```

**Parameter Definitions:**

* `--use-conda` – specifies to use the conda environments included in the workflow (these are specified in the [envs](workflow_code/envs) directory)
* `--conda-prefix` – indicates where the needed conda environments will be stored. Adding this option will also allow the same conda environments to be re-used when processing additional datasets, rather than making new environments each time you run the workflow. The value listed for this option, `${CONDA_PREFIX}/envs`, points to the default location for conda environments (note: the variable `${CONDA_PREFIX}` will be expanded to the appropriate location on whichever system it is run on).
* `-j` – assigns the number of jobs Snakemake should run concurrently
* `-p` – specifies to print out each command being run to the screen

See `snakemake -h` and [Snakemake's documentation](https://snakemake.readthedocs.io/en/stable/) for more options and details.

---
