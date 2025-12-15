# GeneLab Proteomics Consensus Processing Workflow

> GeneLab, part of [NASA's Open Science Data Repository (OSDR)](https://www.nasa.gov/osdr), has wrapped each step of the Proteomics processing pipeline ([PPP](https://github.com/nasa/GeneLab_Data_Processing/tree/master/Proteomics)) into a Nextflow workflow with validation and verification of output files built in after each step. This repository contains the Nextflow workflow code (NF_Proteomics) along with instructions for installation and usage. Exact workflow run info and PPP version used to process specific datasets that have been released are available in the \*nextflow_processing_info.txt file on the [Open Science Data Repository (OSDR)](https://osdr.nasa.gov/bio/repo/), which can be found under 'Files' -> 'GeneLab Processed Proteomics Files' -> 'Supplemental Materials'.

## General Workflow Information

### Implementation Tools

The current GeneLab Proteomics processing pipelines (PPP) for label-free quantification with match-between-runs ([GL-DPPD-[LFQ-MBR]](../../Pipeline_GL-DPPD-[LFQ-MBR]_Versions/GL-DPPD-[LFQ-MBR].md)), TMT10-plex ([GL-DPPD-[TMT10]](../../Pipeline_GL-DPPD-[TMT10]_Versions/GL-DPPD-[TMT10].md)), TMT16-plex ([GL-DPPD-[TMT16]](../../Pipeline_GL-DPPD-[TMT16]_Versions/GL-DPPD-[TMT16].md)), and TMT16-plex phosphoproteomics ([GL-DPPD-[TMT16-phospho]](../../Pipeline_GL-DPPD-[TMT16-phospho]_Versions/GL-DPPD-[TMT16-phospho].md)) are implemented as a single [Nextflow](https://nextflow.io/) DSL2 workflow that utilizes [Singularity](https://docs.sylabs.io/guides/3.10/user-guide/introduction.html) to run all tools in containers. This workflow (NF_Proteomics) is run using the command line interface (CLI) of any unix-based system. While knowledge of creating workflows in Nextflow is not required to run the workflow as is, [the Nextflow documentation](https://nextflow.io/docs/latest/index.html) is a useful resource for users who want to modify and/or extend this workflow. See the [NF_Proteomics Workflow & Subworkflows](#nf_proteomics-workflow--subworkflows) section below for more details on the NF_Proteomics workflow, including installation and execution information.

<br>

# NF_Proteomics Workflow & Subworkflows

### NF_Proteomics Resource Requirements

The table below details the default maximum resource allocations for individual Nextflow processes.

| Workflow Type | Default CPU Cores | Default Memory |
|---------------|-------------------|----------------|
| All workflows | 8                 | 64 GB          |

> **Note:** These per-process resource allocations are defaults. They can be adjusted by modifying `cpus` and `memory`  directives in the configuration files: [`local.config`](workflow_code/conf/local.config) (local execution) and [`slurm.config`](workflow_code/conf/slurm.config) (SLURM clusters).

> **Click links below to show/hide workflow diagrams**

<details open>
<summary>NF_Proteomics workflow</summary>
<p align="center">
<a href="images/draft_pipeline.png"><img src="images/draft_pipeline.png"></a>
</p>
</details>

<!-- <details>
<summary>NF_Proteomics workflow for GL-DPPD-[TMT10]</summary>
<p align="center">
<a href="images/draft_pipeline.png"><img src="images/draft_pipeline.png"></a>
</p>
</details>

<details>
<summary>NF_Proteomics workflow for GL-DPPD-[TMT16]</summary>
<p align="center">
<a href="images/draft_pipeline.png"><img src="images/draft_pipeline.png"></a>
</p>
</details>

<details>
<summary>NF_Proteomics workflow for GL-DPPD-[TMT16-phospho]</summary>
<p align="center">
<a href="images/draft_pipeline.png"><img src="images/draft_pipeline.png"></a>
</p>
</details> -->

---
The NF_Proteomics workflow uses a generic workflow that runs the same steps for all quantification methods, but calls FragPipe with different commands and configurations based on the specified workflow mode (LFQ-MBR, TMT10, TMT16, or TMT16-phospho). The workflow mode is determined from the runsheet data acquisition type and labeling method, or can be specified via a workflow configuration file. The pipeline documents linked below describe what happens during each workflow run for each specific quantification method.

Below is a description of each subworkflow and the additional output files generated that are not already indicated in the [GL-DPPD-[LFQ-MBR]](../../Pipeline_GL-DPPD-[LFQ-MBR]_Versions/GL-DPPD-[LFQ-MBR].md), [GL-DPPD-[TMT10]](../../Pipeline_GL-DPPD-[TMT10]_Versions/GL-DPPD-[TMT10].md), [GL-DPPD-[TMT16]](../../Pipeline_GL-DPPD-[TMT16]_Versions/GL-DPPD-[TMT16].md), and [GL-DPPD-[TMT16-phospho]](../../Pipeline_GL-DPPD-[TMT16-phospho]_Versions/GL-DPPD-[TMT16-phospho].md) pipeline documents:

1. **Analysis Staging Subworkflow**

   - Description:
     - This subworkflow extracts the metadata parameters (e.g. data_type, sample information) needed for processing from the OSD/GLDS ISA archive and retrieves the raw mass spectrometry files hosted on the [Open Science Data Repository (OSDR)](https://osdr.nasa.gov/bio/repo/).
       > *OSD/GLDS ISA archive*: ISA directory containing Investigation, Study, and Assay (ISA) metadata files for a respective GLDS dataset - the *ISA.zip file is located under 'Files' -> 'Study Metadata Files' for any GeneLab Data Set (GLDS) in the [OSDR](https://osdr.nasa.gov/bio/repo/).

2. **Proteomics Consensus Pipeline Subworkflow**

   - Description:
     - This subworkflow processes the staged raw data and metadata parameters from the Analysis Staging Subworkflow. The workflow automatically selects the appropriate FragPipe sequential tools and configurations based on the runsheet `data_type` column and labeling information, or the workflow mode can be explicitly specified via a workflow configuration file. The selected mode determines which of the following pipelines is executed:
       - [The GeneLab LFQ-MBR Pipeline](../../Pipeline_GL-DPPD-[LFQ-MBR]_Versions/GL-DPPD-[LFQ-MBR].md) for label-free quantification with match-between-runs
       - [The GeneLab TMT10 Pipeline](../../Pipeline_GL-DPPD-[TMT10]_Versions/GL-DPPD-[TMT10].md) for TMT10-plex labeling
       - [The GeneLab TMT16 Pipeline](../../Pipeline_GL-DPPD-[TMT16]_Versions/GL-DPPD-[TMT16].md) for TMT16-plex labeling
       - [The GeneLab TMT16-phospho Pipeline](../../Pipeline_GL-DPPD-[TMT16-phospho]_Versions/GL-DPPD-[TMT16-phospho].md) for TMT16-plex phosphoproteomics labeling

<br>

---
## Utilizing the Workflow

1. [Install Nextflow and Singularity](#1-install-nextflow-and-singularity)  
   1a. [Install Nextflow](#1a-install-nextflow)  
   1b. [Install Singularity](#1b-install-singularity)
2. [Download the Workflow Files](#2-download-the-workflow-files)  
3. [Fetch Singularity Images](#3-fetch-singularity-images)  
4. [Run the Workflow](#4-run-the-workflow)  
   4a. [Approach 1: Run the workflow on a GeneLab Proteomics dataset](#4a-approach-1-run-the-workflow-on-a-genelab-proteomics-dataset)   
   4b. [Approach 2: Run the workflow on a GeneLab Proteomics dataset with a custom reference proteome fasta file](#4b-approach-2-run-the-workflow-on-a-genelab-dataset-with-a-custom-reference-proteome-fasta-file)  
   4c. [Approach 3: Run the workflow on a custom dataset](#4c-approach-3-run-the-workflow-on-a-custom-dataset)  
5. [Additional Output Files](#5-additional-output-files)  

<br>

---

### 1. Install Nextflow and Singularity 

#### 1a. Install Nextflow

Nextflow can be installed either through [Anaconda](https://anaconda.org/bioconda/nextflow) or as documented on the [Nextflow documentation page](https://www.nextflow.io/docs/latest/getstarted.html).

> Note: If you want to install Anaconda, we recommend installing a Miniforge version appropriate for your system, as documented on the [conda-forge website](https://conda-forge.org/download/), where you can find basic binaries for most systems. More detailed miniforge documentation is available in the [miniforge github repository](https://github.com/conda-forge/miniforge).
> 
> Once conda is installed on your system, you can install the latest version of Nextflow by running the following commands:
> 
> ```bash
> conda install -c bioconda nextflow
> nextflow self-update
> ```

<br>

#### 1b. Install Singularity

Singularity is a container platform that allows usage of containerized software. This enables the GeneLab PPP workflow to retrieve and use all software required for processing without the need to install the software directly on the user's system.

We recommend installing Singularity on a system wide level as per the associated [documentation](https://docs.sylabs.io/guides/3.10/admin-guide/admin_quickstart.html).

> Note: Singularity is also available through [Anaconda](https://anaconda.org/conda-forge/singularity).

> Note: Alternatively, Docker can be used in place of Singularity. See the [Docker CE installation documentation](https://docs.docker.com/engine/install/).

<br>

---

### 2. Download the Workflow Files

All files required for utilizing the NF_Proteomics GeneLab workflow for processing Proteomics data are in the [workflow_code](workflow_code) directory. To get a 
copy of latest NF_Proteomics version on to your system, the code can be downloaded as a zip file from the release page then unzipped after downloading by running the following commands: 

```bash
wget https://github.com/nasa/GeneLab_Proteomics_Workflow/releases/download/NF_PPP_1.0.0/NF_PPP_1.0.0.zip

unzip NF_PPP_1.0.0.zip
```

<br>

---

### 3. Fetch Singularity Images

Although Nextflow can fetch Singularity images from a url, doing so may cause issues as detailed [here](https://github.com/nextflow-io/nextflow/issues/1210).

To avoid this issue, run the following command to fetch the Singularity images prior to running the NF_PPP workflow:
> Note: This command should be run in the location containing the `NF_PPP_1.0.0` directory that was downloaded in [step 2](#2-download-the-workflow-files) above. Depending on your network speed, fetching the images will take ~20 minutes. Approximately 8GB of RAM is needed to download and build the Singularity images.

```bash
bash NF_PPP_1.0.0/bin/prepull_singularity.sh NF_PPP_1.0.0/config/by_docker_image.config
```


Once complete, a `singularity` folder containing the Singularity images will be created. Run the following command to export this folder as a Nextflow configuration environment variable to ensure Nextflow can locate the fetched images:

```bash
export NXF_SINGULARITY_CACHEDIR=$(pwd)/singularity
```

<br>

---

### 4. Run the Workflow

While in the location containing the `NF_PPP_1.0.0` directory that was downloaded in [step 2](#2-download-the-workflow-files), you are now able to run the workflow.

 Below are examples of how to run the NF_Proteomics workflow:
> Note: Nextflow commands use both single hyphen arguments (e.g. -help) that denote general nextflow arguments and double hyphen arguments (e.g. --reference_version) that denote workflow specific parameters.  Take care to use the proper number of hyphens for each argument.

> Note: To use Docker instead of Singularity, use `-profile docker` in the Nextflow run command. Nextflow will automatically pull images as needed.

> Note: The `-resume` parameter can be used to resume a previously interrupted workflow from where it left off (see [Nextflow documentation](https://www.nextflow.io/docs/latest/getstarted.html#modify-and-resume)) or to restart the workflow from a specific point by changing relevant parameters, which will re-execute that process and all downstream affected processes.

<br>

#### 4a. Approach 1: Run the workflow on a GeneLab Proteomics dataset with Uniprot refence proteome ID(s)

```bash
nextflow run NF_PPP_1.0.0/main.nf \ 
   -profile singularity,local \
   --accession OSD-581 \ 
   --msfragger_fragment_mass_tolerance 300 \
   --uniprot_id UP001231189
```

<br>

#### 4b. Approach 2: Run the workflow on a GeneLab Proteomics dataset with a custom reference proteome fasta file

```bash
nextflow run NF_PPP_1.0.0/main.nf \ 
   -profile singularity,local \
   --runsheet </path/to/runsheet> \ 
   --msfragger_fragment_mass_tolerance 300 \
   --reference_proteome </path/to/fasta>
```

> Note: Specifications for creating a runsheet manually are described [here](examples/runsheet/README.md).

<br>

#### 4c. Approach 3: Run the workflow on a custom dataset

```bash
nextflow run NF_PPP_1.0.0/main.nf \ 
   -profile singularity,local \
   --runsheet_path </path/to/runsheet> \ 
   --msfragger_fragment_mass_tolerance 300 \ 
   --uniprot_id UP001231189 
```

> Note: Specifications for creating a runsheet manually are described [here](examples/runsheet/README.md).

<br>


#### Required Parameters For All Approaches:

* `NF_PPP_1.0.0/main.nf` - Instructs Nextflow to run the NF_Proteomics workflow 

* `-profile` - Specifies the configuration profile(s) to load, `singularity` instructs Nextflow to setup and use singularity for all software called in the workflow; use `local` for local execution ([local.config](workflow_code/conf/local.config)) or `slurm` for SLURM cluster execution ([slurm.config](workflow_code/conf/slurm.config))
  > Note: The output directory will be named `GLDS-#` when using a OSD or GLDS accession as input, or `results` when running the workflow with only a runsheet as input.

* `--msfragger_fragment_mass_tolerance` - Specifies the fragment mass tolerance for MSFragger search (in ppm). This parameter is required for all approaches.


<br>

**Additional Required Parameters For [Approach 1](#4a-approach-1-run-the-workflow-on-a-genelab-proteomics-dataset):**

* `--accession` - The OSD or GLDS ID for the dataset to be processed, eg. `GLDS-194` or `OSD-194`

* `--uniprot_id` - UniProt proteome ID(s) (e.g., `UP001231189`, ). The workflow will download the proteome FASTA from UniProt.

<br>

**Additional Required Parameters For [Approach 2](#4b-approach-2-run-the-workflow-on-a-genelab-dataset-with-a-custom-reference-proteome-fasta-file):**

* `--runsheet` - Path to the runsheet file containing sample metadata and input file paths

* `--reference_proteome` - Path to a custom reference proteome FASTA file

<br>

**Additional Required Parameters For [Approach 3](#4c-approach-3-run-the-workflow-on-a-custom-dataset):**

* `--runsheet_path` - Path to a local runsheet file containing sample metadata and input file paths

* `--uniprot_id` - UniProt proteome ID (e.g., `UP001231189`). The workflow will download the proteome FASTA from UniProt.

<br>

#### Optional Parameters:

* `--stub` - stub

* **Stub Parameters** - Options for stb:

  * `--stub` - stub

<br>

**Additional Optional Parameters:**

All parameters listed above and additional optional arguments for the Proteomics workflow, including debug related options that may not be immediately useful for most users, can be viewed by running the following command:

```bash
nextflow run NF_PPP_1.0.0/main.nf --help
```

See `nextflow run -h` and [Nextflow's CLI run command documentation](https://nextflow.io/docs/latest/cli.html#run) for more options and details common to all nextflow workflows.

<br>

---

### 5. Additional Output Files

The outputs from the Analysis Staging Subworkflow are described below:
> Note: The outputs from the Proteomics Consensus Pipeline Subworkflow are documented in the [GL-DPPD-[LFQ-MBR]](../../Pipeline_GL-DPPD-[LFQ-MBR]_Versions/GL-DPPD-[LFQ-MBR].md), [GL-DPPD-[TMT10]](../../Pipeline_GL-DPPD-[TMT10]_Versions/GL-DPPD-[TMT10].md), [GL-DPPD-[TMT16]](../../Pipeline_GL-DPPD-[TMT16]_Versions/GL-DPPD-[TMT16].md), and [GL-DPPD-[TMT16-phospho]](../../Pipeline_GL-DPPD-[TMT16-phospho]_Versions/GL-DPPD-[TMT16-phospho].md) processing protocols.

**Analysis Staging Subworkflow**

   - Output:
     - Metadata/\*_proteomics_v1_runsheet.csv (table containing metadata required for processing, including the raw data files location)
     - Metadata/\*-ISA.zip (the ISA archive of the OSD datasets to be processed, downloaded from the OSDR)
   
**Processing Information Archive**

   - Output:
     - GeneLab/processing_info_GLProteomics.zip (Archive containing workflow execution metadata)
       - processing_info/samples.txt (single column list of all sample names in the dataset)
       - processing_info/nextflow_log_GLProteomics.txt (Nextflow execution logs captured via `nextflow log`)
       - processing_info/nextflow_run_command_GLProteomics.txt (Exact command line used to initiate the workflow)

<br>

Standard Nextflow resource usage logs are also produced as follows:
> Further details about these logs can also found within [this Nextflow documentation page](https://www.nextflow.io/docs/latest/tracing.html#execution-report).

**Nextflow Resource Usage Logs**

   - Output:
     - nextflow_info/execution_report_{timestamp}.html (an html report that includes metrics about the workflow execution including computational resources and exact workflow process commands)
     - nextflow_info/execution_timeline_{timestamp}.html (an html timeline for all processes executed in the workflow)
     - nextflow_info/execution_trace_{timestamp}.txt (an execution tracing file that contains information about each process executed in the workflow, including: submission time, start time, completion time, cpu and memory used, machine-readable output)
     - nextflow_info/pipeline_dag_{timestamp}.html (a visualization of the workflow process DAG)

<br>

---

# Licenses

The software for the Proteomics pipeline and workflow is released under the [NASA Open Source Agreement (NOSA) Version 1.3](License/RNA_Sequencing_NOSA_License.pdf).


### 3rd Party Software Licenses

Licenses for the 3rd party open source software utilized in the Proteomics pipeline and workflow can be found in the [3rd_Party_Licenses sub-directory](License/3rd_Party_Licenses). 

<br>

---

## Notices

Copyright © 2025 United States Government as represented by the Administrator of the National Aeronautics and Space Administration.  All Rights Reserved. 

### Disclaimers

No Warranty: THE SUBJECT SOFTWARE IS PROVIDED "AS IS" WITHOUT ANY WARRANTY OF ANY KIND, EITHER EXPRESSED, IMPLIED, OR STATUTORY, INCLUDING, BUT NOT LIMITED TO, ANY WARRANTY THAT THE SUBJECT SOFTWARE WILL CONFORM TO SPECIFICATIONS, ANY IMPLIED WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, OR FREEDOM FROM INFRINGEMENT, ANY WARRANTY THAT THE SUBJECT SOFTWARE WILL BE ERROR FREE, OR ANY WARRANTY THAT DOCUMENTATION, IF PROVIDED, WILL CONFORM TO THE SUBJECT SOFTWARE. THIS AGREEMENT DOES NOT, IN ANY MANNER, CONSTITUTE AN ENDORSEMENT BY GOVERNMENT AGENCY OR ANY PRIOR RECIPIENT OF ANY RESULTS, RESULTING DESIGNS, HARDWARE, SOFTWARE PRODUCTS OR ANY OTHER APPLICATIONS RESULTING FROM USE OF THE SUBJECT SOFTWARE.  FURTHER, GOVERNMENT AGENCY DISCLAIMS ALL WARRANTIES AND LIABILITIES REGARDING THIRD-PARTY SOFTWARE, IF PRESENT IN THE ORIGINAL SOFTWARE, AND DISTRIBUTES IT "AS IS."

Waiver and Indemnity:  RECIPIENT AGREES TO WAIVE ANY AND ALL CLAIMS AGAINST THE UNITED STATES GOVERNMENT, ITS CONTRACTORS AND SUBCONTRACTORS, AS WELL AS ANY PRIOR RECIPIENT.  IF RECIPIENT'S USE OF THE SUBJECT SOFTWARE RESULTS IN ANY LIABILITIES, DEMANDS, DAMAGES, EXPENSES OR LOSSES ARISING FROM SUCH USE, INCLUDING ANY DAMAGES FROM PRODUCTS BASED ON, OR RESULTING FROM, RECIPIENT'S USE OF THE SUBJECT SOFTWARE, RECIPIENT SHALL INDEMNIFY AND HOLD HARMLESS THE UNITED STATES GOVERNMENT, ITS CONTRACTORS AND SUBCONTRACTORS, AS WELL AS ANY PRIOR RECIPIENT, TO THE EXTENT PERMITTED BY LAW.  RECIPIENT'S SOLE REMEDY FOR ANY SUCH MATTER SHALL BE THE IMMEDIATE, UNILATERAL TERMINATION OF THIS AGREEMENT. 

The "GeneLab RNA Sequencing Processing Pipeline and Workflow" software also makes use of 3rd party Open Source software, released under the licenses indicated above.  A complete listing of 3rd Party software notices and licenses made use of in "GeneLab RNA Sequencing Processing Pipeline and Workflow" can be found in the [3rd Party Licenses README.md](License/3rd_Party_Licenses/README.md) file. 

<br>

---
**Developed by:**  
A

**Maintained by:**  
B

**Contributors:**
C
