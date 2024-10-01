# GL_RefAnnotTable Workflow Information and Usage Instructions

## General workflow info
The current GeneLab Reference Annotation Table (GL_RefAnnotTable-A) pipeline is implemented as an R workflow that can be run from a command line interface (CLI) using bash. The workflow can be used even if you are unfamiliar with R, but if you want to learn more about R, visit the [R-project about page here](https://www.r-project.org/about.html). Additionally, an introduction to R along with installation help and information about using R for bioinformatics can be found [here at Happy Belly Bioinformatics](https://astrobiomike.github.io/R/basics).  

## Utilizing the workflow

1. [Install R and R packages](#1-install-r-and-r-packages)  
2. [Download the workflow files](#2-download-the-workflow-files)  
3. [Setup Execution Permission for Workflow Scripts](#3-setup-execution-permission-for-workflow-scripts)
4. [Run the workflow](#4-run-the-workflow)  
5. [Run the annotations database creation function as a stand-alone script](#5-run-the-annotations-database-creation-function-as-a-stand-alone-script)
6. [Run the Workflow Using Docker or Singularity](#6-run-the-workflow-using-docker-or-singularity)
<br>

### 1. Install R and R packages

We recommend installing R via the [Comprehensive R Archive Network (CRAN)](https://cran.r-project.org/) as follows: 

1. Select the [CRAN Mirror](https://cran.r-project.org/mirrors.html) closest to your location.
2. Click the link under the "Download and Install R" section that's consistent with your machine.
3. Click on the R-4.4.0 package consistent with your machine to download.
4. Double click on the R-4.4.0.pkg downloaded in step 3 and follow the installation instructions.

Once R is installed, open a CLI terminal and run the following command to activate R:

```bash
R
```
`
Within an active R environment, run the following commands to install the required R packages:

```R
install.packages("tidyverse")

install.packages("BiocManager")

BiocManager::install("STRINGdb")
BiocManager::install("PANTHER.db")
BiocManager::install("rtracklayer")
BiocManager::install("AnnotationForge")
BiocManager::install("biomaRt")
BiocManager::install("GO.db")
```

<br>

### 2. Download the Workflow Files

All files required for utilizing the GL_RefAnnotTable-A workflow for generating reference annotation tables are in the [workflow_code](workflow_code) directory. To get a copy of latest GL_RefAnnotTable version on to your system, run the following command:

```bash
curl -LO https://github.com/nasa/GeneLab_Data_Processing/releases/download/GL_RefAnnotTable-A_1.1.0/GL_RefAnnotTable-A_1.1.0.zip
``` 

<br>

### 3. Setup Execution Permission for Workflow Scripts

Once you've downloaded the GL_RefAnnotTable-A workflow directory as a zip file, unzip the workflow then `cd` into the GL_RefAnnotTable-A_1.1.0 directory on the CLI. Next, run the following command to set the execution permissions for the R script:

```bash
unzip GL_RefAnnotTable-A_1.1.0.zip
cd GL_RefAnnotTable-A_1.1.0
chmod -R u+x *R
```

<br>

### 4. Run the Workflow

While in the GL_RefAnnotTable workflow directory, you are now able to run the workflow. Below is an example of how to run the workflow to build an annotation table for Mus musculus (mouse):

```bash
Rscript GL-DPPD-7110-A_build-genome-annots-tab.R 'Mus musculus'
```

**Input data:**

- No input files are required. Specify the target organism using a positional command line argument. `Mus musculus` is used in the example above. To see a list of all available organisms, run `Rscript GL-DPPD-7110-A_build-genome-annots-tab.R` without positional arguments. The correct argument for each organism can also be found in the 'species' column of the [GL-DPPD-7110-A_annotations.csv](../../Pipeline_GL-DPPD-7110_Versions/GL-DPPD-7110-A/GL-DPPD-7110-A_annotations.csv)

- Optional: a reference table CSV can be supplied as a second positional argument instead of using the default [GL-DPPD-7110-A_annotations.csv](../../Pipeline_GL-DPPD-7110_Versions/GL-DPPD-7110-A/GL-DPPD-7110-A_annotations.csv)

**Output data:**

- *-GL-annotations.tsv (Tab delineated table of gene annotations)
- *-GL-build-info.txt (Text file containing information used to create the annotation table, including tool and tool versions and date of creation)

### 5. Run the annotations database creation function as a stand-alone script

When the workflow is run, if the reference table does not specify an annotations database for the target_organism in the `annotations` column, the `install_annotations` function, defined in the `install-org-db.R` script, will be executed. This script will locally create and install an annotations database R package using AnnotationForge. This function can also be run as a stand-alone script from the command line:

```bash
Rscript install-org-db.R 'Bacillus subtilis' /path/to/GL-DPPD-7110-A_annotations.csv
```

**Input data:**

- The target organism must be specified as the first positional command line argument, `Bacillus subtilis` is used in the example above. The correct argument for each organism can be found in the 'species' column of the [GL-DPPD-7110-A_annotations.csv](../../Pipeline_GL-DPPD-7110_Versions/GL-DPPD-7110-A/GL-DPPD-7110-A_annotations.csv)

- The path to a local reference table must also be supplied as the second positional argument

**Output data:**

- org.*.eg.db/ (species-specific annotation database, as a local R package)

### 6. Run the Workflow Using Docker or Singularity

Rather than running the workflow in your local environment, you can use a Docker or Singularity container. This method ensures that all dependencies are correctly installed.

1. **Pull the container image:**

   Docker:
   ```bash
   docker pull quay.io/nasa_genelab/gl-refannottable:v1.0.0
   ```

   Singularity:
   ```bash
   singularity pull docker://quay.io/nasa_genelab/gl-refannottable:v1.0.0
   ```

2. **Download the workflow files:**

   ```bash
   curl -LO https://github.com/nasa/GeneLab_Data_Processing/releases/download/GL_RefAnnotTable-A_1.1.0/GL_RefAnnotTable-A_1.1.0.zip
   unzip GL_RefAnnotTable-A_1.1.0.zip
   ```

3. **Run the workflow:**

   Docker:
   ```bash
   docker run -it -v $(pwd)/GL_RefAnnotTable-A_1.1.0:/work \
     quay.io/nasa_genelab/gl-refannottable:v1.0.0 \
     bash -c "cd /work && Rscript GL-DPPD-7110-A_build-genome-annots-tab.R 'Mus musculus'"
   ```

   Singularity:
   ```bash
   singularity exec -B $(pwd)/GL_RefAnnotTable-A_1.1.0:/work \
     gl-refannottable_v1.0.0.sif \
     bash -c "cd /work && Rscript GL-DPPD-7110-A_build-genome-annots-tab.R 'Mus musculus'"
   ```

**Input data:**

- No input files are required. Specify the target organism using a positional command line argument. `Mus musculus` is used in the example above. To see a list of all available organisms, run `Rscript GL-DPPD-7110-A_build-genome-annots-tab.R` without positional arguments. The correct argument for each organism can also be found in the 'species' column of the [GL-DPPD-7110-A_annotations.csv](../../Pipeline_GL-DPPD-7110_Versions/GL-DPPD-7110-A/GL-DPPD-7110-A_annotations.csv)

- Optional: a reference table CSV can be supplied as a second positional argument instead of using the default [GL-DPPD-7110-A_annotations.csv](../../Pipeline_GL-DPPD-7110_Versions/GL-DPPD-7110-A/GL-DPPD-7110-A_annotations.csv)

**Output data:**

- *-GL-annotations.tsv (Tab delineated table of gene annotations)
- *-GL-build-info.txt (Text file containing information used to create the annotation table, including tool and tool versions and date of creation)