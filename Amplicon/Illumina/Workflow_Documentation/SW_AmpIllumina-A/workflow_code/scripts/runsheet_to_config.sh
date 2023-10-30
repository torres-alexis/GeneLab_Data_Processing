#!/bin/bash

# Creates a config.yaml file https://github.com/nasa/GeneLab_Data_Processing/blob/master/Amplicon/Illumina/Workflow_Documentation/SW_AmpIllumina-A/workflow_code/config.yaml
# Based on values in the runsheet

# Check number of arguments passed to the script
if [ "$#" -ne 3 ]; then
    echo "Error: Incorrect number of arguments.\nUsage: $0 <Path to Runsheet CSV> <Path of Raw Reads Directory>
Use full paths or paths relative to the location of the Snakefile.
Ex: $0 ./runsheet.csv ../raw_reads/ ./workflow/"
    exit 1
fi

csv_file=$1
raw_reads_directory=$2
output_dir=$3


sample_ids_file="unique-sample-IDs.txt"

primers_linked="TRUE"

# Make the output directory if it doesn't exist
mkdir -p "$output_dir"

# Check if the runsheet file exists
if [ ! -f "$csv_file" ]; then
    echo "Error: "$csv_file" not found."
    exit 1
fi

# Function to get column number by its name
get_col_num() {
    header=$1
    file=$2
    awk -v var="$header" -F, 'NR==1 {for (i=1; i<=NF; i++) if ($i == var) print i; exit}' "$file"
}

# Create unique-sample-IDs.txt
sample_name_col_num=$(get_col_num "Sample Name" "$csv_file")



# Extract unique sample names and write to sample_names.txt
awk -F, -v col="$sample_name_col_num" 'NR > 1 && !seen[$col]++ {print $col}' $csv_file > $sample_ids_file
echo "Unique sample names saved to $sample_ids_file"


# Get the column number for "paired_end"
paired_col=$(get_col_num "paired_end" "$csv_file")
# Extract paired_end values, convert to uppercase, and take the unique value
paired=$(awk -v col="$paired_col" -F, 'NR > 1 {print toupper($col)}' $csv_file | uniq)
# # Count the number of lines in paired_end
# line_count=$(echo "$paired_end" | wc -l)

# if [ "$line_count" -gt 1 ]; then
#     echo "Error: Multiple values detected for paired_end. Please ensure consistency in the runsheet."
#     exit 1
# fi
# Check the value and set data_type accordingly
if [ "$paired" == "TRUE" ]; then
    data_type="PE"
else
    data_type="SE"
fi

# Get the column number for "raw_R1_suffix"
raw_r1_suffix_col=$(get_col_num "raw_R1_suffix" "$csv_file")
# Extract raw_R1_suffix values and take the unique value
raw_r1_suffix=$(awk -v col="$raw_r1_suffix_col" -F, 'NR > 1 {print $col}' $csv_file | uniq)
if [ "$paired" == "TRUE" ]; then
    # Get the column number for "raw_R2_suffix"
    raw_r2_suffix_col=$(get_col_num "raw_R2_suffix" "$csv_file")
    # Extract raw_R2_suffix values and take the unique value
    raw_r2_suffix=$(awk -v col="$raw_r2_suffix_col" -F, 'NR > 1 {print $col}' $csv_file | uniq)
else
    raw_r2_suffix=""
fi

# Get the column number for "F_primer"
f_primer_col=$(get_col_num "F_Primer" "$csv_file")
# Extract F_primer values and take the unique value
f_primer=$(awk -v col="$f_primer_col" -F, 'NR > 1 {print $col}' $csv_file | uniq)
if [ "$paired" == "TRUE" ]; then
    # Get the column number for "R_primer"
    r_primer_col=$(get_col_num "R_Primer" "$csv_file")
    # Extract R_primer values and take the unique value
    r_primer=$(awk -v col="$r_primer_col" -F, 'NR > 1 {print $col}' $csv_file | uniq)

    # Also get linked primers if paired, using reverse complement of other primer:
    f_linked_primer="^${f_primer}...$(echo "$r_primer" | rev | tr 'ACGTMRWSYKVHDBN' 'TGCAMKWSRMBDHVNY')"
    r_linked_primer="^${r_primer}...$(echo "$f_primer" | rev | tr 'ACGTMRWSYKVHDBN' 'TGCAMKWSRMBDHVNY')"

else
    r_primer=""
    f_linked_primer=""
    r_linked_primer=""
fi

# Get the column number for "Parameter Value[Library Selection]"
target_region_col=$(get_col_num "Parameter Value[Library Selection]" "$csv_file")
# Extract target_region values and take the unique value
target_region=$(awk -v col="$target_region_col" -F, 'NR > 1 {print $col}' $csv_file | uniq)


# Write header to config.yaml
echo "############################################################################################
## Configuration file for GeneLab Illumina amplicon processing workflow                   ##
## Developed by Michael D. Lee (Mike.Lee@nasa.gov)                                        ##
############################################################################################

############################################################
##################### VARIABLES TO SET #####################
############################################################

###########################################################################
##### These need to match what is specific to our system and our data #####
###########################################################################
" > config.yaml

# Write variables to config.yaml


echo "## Path to runsheet:
runsheet:
    \"$csv_file\"
" >> config.yaml

echo "## Set to "PE" for paired-end, "SE" for single-end. (if single-end, enter appropriate info below for the "_R1_" variables and ignore the "_R2_" ones)
data_type:
    \"$data_type\"
" >> config.yaml

echo "## single-column file with unique sample identifiers:
sample_info_file:
    \"$sample_ids_file\"
" >> config.yaml

echo "## input reads directory (can be relative to workflow directory, or needs to be full path):
raw_reads_dir:
    \"$raw_reads_directory\"
" >> config.yaml

echo "## raw read suffixes (region following the unique part of the sample names)
  # e.g. for paired-end data, "Sample-1_R1_raw.fastq.gz" would be "_R1_raw.fastq.gz" for 'raw_R1_suffix' below
  # e.g. if single-end, "Sample-1.fastq.gz" would be ".fastq.gz" for 'raw_R1_suffix' below, and 'raw_R2_suffix' won't be used
raw_R1_suffix:
    \"$raw_r1_suffix\"
raw_R2_suffix:
    \"$raw_r2_suffix\"
" >> config.yaml

echo "## if we are trimming primers or not ("TRUE", or "FALSE")
trim_primers:
    \"TRUE\"
" >> config.yaml

echo "## primer sequences if we are trimming them (include anchoring symbols, e.g. '^', as needed, see: https://cutadapt.readthedocs.io/en/stable/guide.html#adapter-types)
F_primer:
    \"^$f_primer\"
R_primer:
    \"^$r_primer\"
" >> config.yaml

echo "## should cutadapt treat these as linked primers? (https://cutadapt.readthedocs.io/en/stable/recipes.html#trimming-amplicon-primers-from-paired-end-reads)
primers_linked:
    \"$primers_linked\"
" >> config.yaml

echo "## if primers are linked, we need to provide them as below, where the second half, following three periods, is the other primer reverse-complemented (see https://cutadapt.readthedocs.io/en/stable/recipes.html#trimming-amplicon-primers-from-paired-end-reads)
  # (can reverse complement while retaining ambiguous bases at this site: http://arep.med.harvard.edu/labgc/adnan/projects/Utilities/revcomp.html)
  # include anchoring symbols, e.g. '^', as needed, see: https://cutadapt.readthedocs.io/en/stable/guide.html#adapter-types
F_linked_primer:
    \"$f_linked_primer\"
R_linked_primer:
    \"$r_linked_primer\"
" >> config.yaml

echo "## discard untrimmed, sets the "--discard-untrimmed" option if TRUE
discard_untrimmed:
    \"TRUE\"
" >> config.yaml

echo "## target region (16S or ITS acceptable; determines which reference database is used for taxonomic classification)
target_region:
    \"$target_region\"
" >> config.yaml

echo "## values to be passed to dada2's filterAndTrim() function:
left_trunc:
    0
right_trunc:
    0
left_maxEE:
    1
right_maxEE:
    1

## minimum length threshold for cutadapt
min_cutadapt_len:
    150


######################################################################
##### The rest only need to be altered if we want to change them #####
######################################################################

## filename suffixes
primer_trimmed_R1_suffix:
    \"_R1_trimmed.fastq.gz\"
primer_trimmed_R2_suffix:
    \"_R2_trimmed.fastq.gz\"

filtered_R1_suffix:
    \"_R1_filtered.fastq.gz\"
filtered_R2_suffix:
    \"_R2_filtered.fastq.gz\"


## output prefix (if needed to distinguish from multiple primer sets, leave as empty string if not, include connecting symbol if adding, e.g. "ITS-")
output_prefix:
    \"\"

## output directories (all relative to processing directory, they will be created if needed)
fastqc_out_dir:
    \"${output_dir}FastQC_Outputs/\"
trimmed_reads_dir:
    \"${output_dir}Trimmed_Sequence_Data/\"
filtered_reads_dir:
    \"${output_dir}Filtered_Sequence_Data/\"
final_outputs_dir:
    \"${output_dir}Final_Outputs/\"

" >> config.yaml


# Append closing notes
echo '############################################################
###################### GENERAL INFO ########################
############################################################
# Workflow is currently equipped to work with paired-end data only, and reads are expected to be gzipped

## example usage command ##
# snakemake --use-conda --conda-prefix ${CONDA_PREFIX}/envs -j 2 -p

# `--use-conda` – this specifies to use the conda environments included in the workflow
# `--conda-prefix` – this allows us to point to where the needed conda environments should be stored. Including this means if we use the workflow on a different dataset somewhere else in the future, it will re-use the same conda environments rather than make new ones. The value listed here, `${CONDA_PREFIX}/envs`, is the default location for conda environments (the variable `${CONDA_PREFIX}` will be expanded to the appropriate location on whichever system it is run on).
# `-j` – this lets us set how many jobs Snakemake should run concurrently (keep in mind that many of the thread and cpu parameters set in the config.yaml file will be multiplied by this)
# `-p` – specifies to print out each command being run to the screen

# See `snakemake -h` for more options and details.' >> config.yaml
# Print message to the user
echo "Config file created with data_type set to $data_type."