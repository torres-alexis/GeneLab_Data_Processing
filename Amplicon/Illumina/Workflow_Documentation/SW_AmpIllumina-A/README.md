# SW_AmpIllumina-A Workflow Information and Usage Instructions


## General workflow info
The current GeneLab Illumina amplicon sequencing data processing pipeline (AmpIllumina), [GL-DPPD-7104-A.md](../../Pipeline_GL-DPPD-7104_Versions/GL-DPPD-7104-A.md), is implemented as a [Snakemake](https://snakemake.readthedocs.io/en/stable/) workflow and utilizes [conda](https://docs.conda.io/en/latest/) environments to install/run all tools. This workflow (SW_AmpIllumina-A) is run using the command line interface (CLI) of any unix-based system. The workflow can be used even if you are unfamiliar with Snakemake and conda, but if you want to learn more about those, [this Snakemake tutorial](https://snakemake.readthedocs.io/en/stable/tutorial/tutorial.html) within [Snakemake's documentation](https://snakemake.readthedocs.io/en/stable/) is a good place to start for that, and an introduction to conda with installation help and links to other resources can be found [here at Happy Belly Bioinformatics](https://astrobiomike.github.io/unix/conda-intro).  

## Utilizing the workflow

- [SW\_AmpIllumina-A Workflow Information and Usage Instructions](#sw_ampillumina-a-workflow-information-and-usage-instructions)
  - [General workflow info](#general-workflow-info)
  - [Utilizing the workflow](#utilizing-the-workflow)
    - [1. Install conda, mamba, and `genelab-utils` package](#1-install-conda-mamba-and-genelab-utils-package)
    - [2. Download the workflow template files](#2-download-the-workflow-template-files)
    - [3. Configure the input files](#3-configure-the-input-files)
      - [3a. Automate Configuration with Scripts](#3a-automate-configuration-with-scripts)
      - [3b. Prepare the runsheet](#3b-prepare-the-runsheet)
      - [3c. Set up unique-sample-IDs.txt](#3c-set-up-unique-sample-idstxt)
      - [3d. Set up config.yaml](#3d-set-up-configyaml)
    - [4. Run the workflow](#4-run-the-workflow)

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
All files required for utilizing the GeneLab workflow for processing Illumina amplicon sequencing data are in the [workflow_code](workflow_code) directory. To get a copy of the latest SW_AmpIllumina-A version on to your system, run the following command:

```bash
GL-get-workflow Amplicon-Illumina
```

This downloaded the workflow into a directory called `SW_AmpIllumina-*/`, with the workflow version number at the end.

> Note: If wanting an earlier version, the wanted version can be provided as an optional argument like so:
> ```bash
> GL-get-workflow Amplicon-Illumina --wanted-version 1.0.0
> ```


### 3. Configure the input files

#### 3a. Automate Configuration with Scripts

If you are working with GeneLab OSDR datasets, the entire process of generating configuration files can be automated using scripts from dp_tools and the [runsheet-to-config.sh](workflow_code/scripts/runsheet-to-config.sh) script. 

You will need to install a specific development branch of dp_tools that contains updates for amplicon sequencing data. This branch can be found at [https://github.com/torres-alexis/dp_tools/tree/amplicon_updates](https://github.com/torres-alexis/dp_tools/tree/amplicon_updates).

You can install this development branch to your active conda environment using pip with the following command:
```bash
pip install git+https://github.com/torres-alexis/dp_tools.git@amplicon_updates
```
**Download the ISA files:** Use the dp_tools script `dpt-get-isa-archive`.
```sh
dpt-get-isa-archive GLDS-###
```
**Generate the runsheet:** Use the dp_tools script `dpt-isa-to-archive` to generate a runsheet from the study's ISA files.
```bash
dpt-isa-to-runsheet --accession GLDS-### --config-type amplicon --config-version Latest --isa-archive /path/to/GLDS-###_GLDS-###-ISA.zip
```

**Create unique-sample-IDs.txt and config.yaml:** Once the runsheet is ready, the [runsheet-to-config.sh](workflow_code/scripts/runsheet-to-config.sh) script located in the scripts directory to generate the unique-sample-IDs.txt and config.yaml files. When running the script, provide either absolute paths or paths relative to the Snakefile's location. For consistency and to avoid path-related errors during execution of the Snakemake workflow, it is recommended to use absolute paths. This ensures that config.yaml is the only input file required in the workflow_code directory when initiating the workflow.


#### 3b. Prepare the runsheet
To process your dataset with the workflow, you should first prepare a runsheet CSV file that contains the required information for your samples. If you are working with a GeneLab OSDR dataset, the `dpt-isa-to-runsheet` script can be used to generate a runsheet from the associated study's ISA files. If you are using other datasets, the runsheet will need to be prepared manually.

The runsheet requires the following columns, in no particular order after the first:

- **Sample Name:** A unique identifier for each sample.
- **Parameter Value[Library Selection]:** Library type, either '16S' or 'ITS'
- **paired_end:** 'TRUE' or 'FALSE' indicating whether the sequencing is paired-end (TRUE) or single-end (FALSE).
- **F_Primer and R_Primer:** The sequences of the forward and reverse primers. Only one primer sequence should be entered in each respective column. For single-end data, only **F_Primer** is required.
- **read1_path** and **read2_path:** For paired-end data, the filenames or pathnames for the forward and reverse sequence files, respectively. For single-end data, only **read1_path** is required.
- **raw_R1_suffix** and **raw_R2_suffix:** The suffixes that identify the raw read files. These columns should contain only one value each. For single-end data, only **raw_R1_suffix** is required.
- **groups:** Define the experimental groups for each sample, using an '&' to differentiate between factors. For example, 'Ground Control & 2 weeks' indicates the sample's conditions related to the Spaceflight and Time factors.


| Sample Name | Parameter Value[Library Selection] | paired_end | F_Primer                     | R_Primer                     | read1_path                                        | read2_path                                        | raw_R1_suffix   | raw_R2_suffix   | Factor Value[Spaceflight] | Factor Value[Time] | groups                   |
|-------------|------------------------------------|------------|------------------------------|------------------------------|---------------------------------------------------|---------------------------------------------------|-----------------|-----------------|----------------------------|-------------------|--------------------------|
| Sample-1    | 16S                                | TRUE       | AGAGTTTGATCCTGGCTCAG         | TTACCGCGGCTGCTGGCAC          | Sample1_R1_raw.fastq.gz filename or pathname     | Sample1_R2_raw.fastq.gz filename or pathname     | _R1_raw.fastq.gz | _R2_raw.fastq.gz | Ground Control             | 2 Weeks           | Ground Control & 2 weeks |
| Sample-2    | 16S                                | TRUE       | AGAGTTTGATCCTGGCTCAG         | TTACCGCGGCTGCTGGCAC          | Sample2_R1_raw.fastq.gz filename or pathname     | Sample2_R2_raw.fastq.gz filename or pathname     | _R1_raw.fastq.gz | _R2_raw.fastq.gz | Space Flight               | 3 Weeks           | Space Flight & 3 weeks   |

After the runsheet is prepared, you can manually configure [unique-sample-IDs.txt](workflow_code/unique-sample-IDs.txt)  and [config.yaml](workflow_code/config.yaml) or use the [runsheet-to-config.sh](workflow_code/scripts/runsheet-to-config.sh) script to generate these two files based on your runsheet.

#### 3c. Set up unique-sample-IDs.txt

You will have to provide a text file containing a single-column list of unique sample identifiers (see an example of how to set this up below). You will also need to indicate the paths to your input data (raw reads) and, if necessary, modify each variable to be consistent with the study you want to process. 

> Note: If you are unfamiliar with how to specify paths, one place you can learn more is [here](https://astrobiomike.github.io/unix/getting-started#the-unix-file-system-structure).  

**Example for how to create a single-column list of unique sample identifiers from your raw data file names**

For example, if you have paired-end read data for 2 samples located in `../Raw_Data/` relative to your workflow directory, that would look like this:

```bash
ls ../Raw_Data/
```

```
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

#### 3d. Set up config.yaml
Once the runsheet CSV file and unique-samples-IDs.txt file are set up, you can modify the variables in the [config.yaml](workflow_code/config.yaml) file as needed. All paths in the configuration file should be absolute paths or paths relative to the Snakefile.

### 4. Run the workflow

While in the directory holding the Snakefile, config.yaml, and other workflow files that you downloaded in [step 2](#2-download-the-workflow-template-files), here is one example command of how to run the workflow:

```bash
snakemake --use-conda --conda-prefix ${CONDA_PREFIX}/envs -j 2 -p
```

**Parameter Definitions:**

* `--use-conda` – specifies to use the conda environments included in the workflow (these are specified in the [envs](workflow_code/envs) directory)
* `--conda-prefix` – indicates where the needed conda environments will be stored. Adding this option will also allow the same conda environments to be re-used when processing additional datasets, rather than making new environments each time you run the workflow. The value listed for this option, `${CONDA_PREFIX}/envs`, points to the default location for conda environments (note: the variable `${CONDA_PREFIX}` will be expanded to the appropriate location on whichever system it is run on).
* `-j` – assigns the number of jobs Snakemake should run concurrently
* `-p` – specifies to print out each command being run to the screen

See `snakemake -h` and [Snakemake's documentation](https://snakemake.readthedocs.io/en/stable/) for more options and details.

---
