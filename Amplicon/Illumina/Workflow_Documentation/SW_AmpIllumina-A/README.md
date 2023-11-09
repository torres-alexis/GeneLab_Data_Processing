# SW_AmpIllumina-A Workflow Information and Usage Instructions


## General workflow info
The current GeneLab Illumina amplicon sequencing data processing pipeline (AmpIllumina), [GL-DPPD-7104-A.md](../../Pipeline_GL-DPPD-7104_Versions/GL-DPPD-7104-A.md), is implemented as a [Snakemake](https://snakemake.readthedocs.io/en/stable/) workflow and utilizes [conda](https://docs.conda.io/en/latest/) environments to install/run all tools. This workflow (SW_AmpIllumina-A) is run using the command line interface (CLI) of any unix-based system. The workflow can be used even if you are unfamiliar with Snakemake and conda, but if you want to learn more about those, [this Snakemake tutorial](https://snakemake.readthedocs.io/en/stable/tutorial/tutorial.html) within [Snakemake's documentation](https://snakemake.readthedocs.io/en/stable/) is a good place to start for that, and an introduction to conda with installation help and links to other resources can be found [here at Happy Belly Bioinformatics](https://astrobiomike.github.io/unix/conda-intro).  

## Utilizing the workflow

1. [Install conda, mamba, and `genelab-utils` package](#1-install-conda-mamba-and-genelab-utils-package)  
2. [Download the workflow template files](#2-download-the-workflow-template-files)  
3. [Set up the runsheet](#3-set-up-the-runsheet)
  3a. [Approach 1: Running the workflow on a GeneLab Illumina amplicon sequencing dataset with automatic generation of the runsheet and retrieval of the raw sequencing data](#3a-approach-1-running-the-workflow-on-a-genelab-illumina-amplicon-sequencing-dataset-with-automatic-generation-of-the-runsheet-and-retrieval-of-the-raw-sequencing-data)
  3b. [Approach 2: Running the workflow on Non-GLDS Datasets using a user-generated runsheet](#3b-approach-2-running-the-workflow-on-non-glds-datasets-using-a-user-generated-runsheet)
4. [Configure config.yaml and unique-sample-IDs.txt](#4-configure-configyaml-and-unique-sample-idstxt)
5. [Run the workflow](#5-run-the-workflow)
 
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

### 2. Download the workflow template files
All files required for utilizing the GeneLab workflow for processing Illumina amplicon sequencing data are in the [workflow_code](workflow_code) directory. 

<!-- To get a copy of the latest SW_AmpIllumina-A version on to your system, run the following command:

```bash
GL-get-workflow Amplicon-Illumina
```

This downloaded the workflow into a directory called `SW_AmpIllumina-*/`, with the workflow version number at the end.

> Note: If wanting an earlier version, the wanted version can be provided as an optional argument like so:
> ```bash
> GL-get-workflow Amplicon-Illumina --wanted-version 1.0.0
> ``` -->

### 3. Set up the runsheet

To process your dataset using the Snakemake workflow, you will need to prepare a runsheet containing the required metadata for your dataset.

#### 3.a Approach 1: Running the workflow on a GeneLab Illumina amplicon sequencing dataset with automatic generation of the runsheet and retrieval of the raw sequencing data

If you are working with a GeneLab OSDR dataset, you can use dp_tools scripts to generate this runsheet.

You will need to install a development branch of dp_tools that contains updates for processing amplicon sequencing data. This branch can be found at [https://github.com/torres-alexis/dp_tools/tree/amplicon_updates](https://github.com/torres-alexis/dp_tools/tree/amplicon_updates).

You can install this development branch to your active conda environment using pip with the following command:
```bash
pip install git+https://github.com/torres-alexis/dp_tools.git@amplicon_updates
```

**Download the ISA files:** Use the dp_tools script `dpt-get-isa-archive` to download a study's ISA files.
```bash
dpt-get-isa-archive --accession GLDS-###
```
**Generate the runsheet:** Use the dp_tools script `dpt-isa-to-archive` to generate a runsheet from the study's ISA files.

```bash
dpt-isa-to-runsheet --accession GLDS-### --config-type amplicon --config-version Latest --isa-archive /path/to/GLDS-###_GLDS-###-ISA.zip
```

When using the [runsheet-to-config.sh](workflow_code/scripts/runsheet-to-config.sh) as detailed in section [4](#4-configure-configyaml-and-unique-sample-idstxt), raw reads will be automatically retrieved and downloaded to the `workflow_code/raw_reads` directory.

#### 3.b Approach 2: Running the workflow on Non-GLDS Datasets using a user-generated runsheet

If you are working with a non-GLDS dataset, you must manually create the runsheet for your dataset to run the workflow.

> Note: Specifications for creating a runsheet manually are described [here](examples/runsheet/README.md).

### 4. Configure config.yaml and unique-sample-IDs.txt

`config.yaml` contains variables which the Snakemake workflow uses to perform its processing. `unique-sample-IDS.txt` contains the unique identifiers for each sample's read file(s).

These two files can be generated from a runsheet using the [runsheet-to-config.sh](workflow_code/scripts/runsheet-to-config.sh) script.

This script should be executed from the [workflow_code/](workflow_code/) folder. When running the script, provide either absolute paths or paths relative to the Snakefile's location. For consistency and to avoid path-related errors during execution of the Snakemake workflow, it is recommended to use absolute paths.

> Note: You may have to make this script executable with the following command: 

```bash
chmod +x /PATH/TO/workflow_code/scripts/runsheet_to_config.sh
```

Here is the usage syntax for using the `runsheet_to_config.sh` script from the [workflow_code/](workflow_code/) directory:

```bash
./scripts/runsheet_to_config.sh -r <path_to_runsheet.csv> [-o <path_to_output_directory>] [-m <minimum_trimmed_read_length>]
```

**Parameter Definitions:**

* `-r` or `--runsheet` – Path to the runsheet CSV file, which must be either absolute or relative to the Snakefile. The runsheet should be formatted correctly with all necessary columns.
* `-o` or `--output` (optional) – Path to the output directory where the Snakefile outputs will be created. Default: `workflow_output/`
* `-m `or `--min_length` (optional) – The minimum length allowed for trimmed reads. If a read is shorter than this length after trimming, it will be discarded. For paired-end reads, if one read of the pair is discarded, the other will be as well. Default: `1`

Example execution from the [workflow_code/](workflow_code/) directory:
```bash
./scripts/runsheet_to_config.sh -r runsheet.csv -o /path/to/custom_outputs_directory/ -m 1
```



Below are instructions to prepare the two files manually if needed.

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

### 5. Run the workflow

While in the directory holding the Snakefile, config.yaml, and other workflow files that you downloaded in [step 2](#2-download-the-workflow-template-files), here is one example command of how to run the workflow:

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