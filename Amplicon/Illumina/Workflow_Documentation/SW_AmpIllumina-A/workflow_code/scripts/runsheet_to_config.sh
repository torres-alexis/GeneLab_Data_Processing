#!/bin/bash

# config.yaml format is based on https://github.com/nasa/GeneLab_Data_Processing/blob/master/Amplicon/Illumina/Workflow_Documentation/SW_AmpIllumina-A/workflow_code/config.yaml

# Initialize default values
output_dir="workflow/"
min_trimmed_read_length=1
    echo ""
# Function to display usage
print_usage() {
    echo "Usage: $0"
    echo "       -r|--runsheet <Path to runsheet CSV,"
    echo "        relative to the Snakefile or as an absolute path."
    echo "       -d|--raw_reads <Path to directory containing the raw reads,"
    echo "        relative to the Snakefile or as an absolute path."
    echo "       [-o|--output <Path to output directory, relative to the Snakefile"
    echo "         or as an absolute path. (default: workflow/)>]"
    echo "       [-m|--min_length <Minimum trimmed read length (default: 1)>]"
    echo "       [-h|--help Display the help menu]"

    echo "Note:"
    echo "       -m specifies the minimum post-trimming read length."
    echo "       For paired end: if one read gets filtered, both reads"
    echo "       in the pair are discarded."

    echo "Example:"
    echo "       $0 -r runsheet.csv -d ../raw_reads/ -o workflow/ -m 1"
}

# Parse named arguments
while [ $# -gt 0 ]; do
	key="$1"
	shift

	case $key in
	-r|--runsheet)
		runsheet_csv="$1"
		shift
		;;
	-d|--raw_reads)
		raw_reads_directory="$1"
		shift
		;;
	-o|--output)
		output_dir="$1"
		shift
		;;
	-m|--minimum-length)
		min_trimmed_read_length="$1"
		shift
		;;
    -h|--help)
        print_usage
        exit 0
        ;;
	*)
        echo "Error: Unknown argument $key"
        echo ""
        print_usage
        exit 1
		;;
	esac
done

# Check if required named arguments are provided
if [[ -z "$runsheet_csv" ]] || [[ -z "$raw_reads_directory" ]]; then
    echo "Error: Missing required arguments."
    echo ""
    print_usage
    exit 1
fi

echo "Runsheet: $runsheet_csv"
echo "Raw Reads Directory: $raw_reads_directory"
echo "Output Directory: $output_dir"
echo "Minimum Trimmed Read Length: $min_trimmed_read_length"

sample_ids_file="unique-sample-IDs.txt"

primers_linked="TRUE"



# Check if the runsheet file exists
if [ ! -f "$runsheet_csv" ]; then
    echo "Error: "$runsheet_csv" not found."
    exit 1
fi

# Function to get column number by its name
get_col_num() {
    header=$1
    file=$2
    awk -v var="$header" -F, 'NR==1 {for (i=1; i<=NF; i++) if ($i == var) print i; exit}' "$file"
}

# Get the column number for "paired_end"
paired_col=$(get_col_num "paired_end" "$runsheet_csv")
if [ -z "$paired_col" ]; then
    echo "Error: 'paired_end' column not found in the runsheet."
    exit 1
fi
# Extract paired_end values, convert to uppercase, and take the unique value
paired=$(awk -v col="$paired_col" -F, 'NR > 1 {print toupper($col)}' $runsheet_csv | uniq)
# # Count the number of lines in paired_end
# line_count=$(echo "$paired_end" | wc -l)

# if [ "$line_count" -gt 1 ]; then
#     echo "Error: Multiple values detected for paired_end. Please ensure consistency in the runsheet."
#     exit 1
# fi
# Check the value and set data_type accordingly
if [ "$paired" == "TRUE" ]; then
    data_type="PE"
    echo "Data Type: Paired-End (PE)"
elif [ "$paired" == "FALSE" ]; then
    data_type="SE"
    echo "Data Type: Single-End (SE)"
else
    echo "Error: Expected 'TRUE' or 'FALSE' in 'paired_end' column."
    exit 1
fi

echo ""
# Get the column number for "raw_R1_suffix"
raw_r1_suffix_col=$(get_col_num "raw_R1_suffix" "$runsheet_csv")
if [ -z "$raw_r1_suffix_col" ]; then
    echo "Error: 'raw_R1_suffix' column not found in the runsheet."
    exit 1
fi
# Extract raw_R1_suffix values and take the unique value
raw_r1_suffix=$(awk -v col="$raw_r1_suffix_col" -F, 'NR > 1 {print $col}' $runsheet_csv | uniq)
if [ "$paired" == "TRUE" ]; then
    # Get the column number for "raw_R2_suffix"
    raw_r2_suffix_col=$(get_col_num "raw_R2_suffix" "$runsheet_csv")
    if [ -z "$raw_r2_suffix_col" ]; then
        echo "Error: 'raw_R2_suffix' column not found in the runsheet."
        exit 1
    fi
    # Extract raw_R2_suffix values and take the unique value
    raw_r2_suffix=$(awk -v col="$raw_r2_suffix_col" -F, 'NR > 1 {print $col}' $runsheet_csv | uniq)
else
    raw_r2_suffix=""
fi


# Create unique-sample-IDs.txt
# sample_name_col_num=$(get_col_num "Sample Name" "$runsheet_csv")
#awk -F, -v col="$sample_name_col_num" 'NR > 1 && !seen[$col]++ {print $col}' $runsheet_csv > $sample_ids_file
# echo "Unique sample names saved to $sample_ids_file"

# Refactored to get the sample IDs from the file names, for cases where sample name col does not match read{1,2}_path names
# Extract sample IDs using raw_r1_suffix
# Get the columns for "read1_path" and "read2_path"
# Extract unique sample names and write to sample_names.txt

read1_path_col=$(get_col_num "read1_path" "$runsheet_csv")
if [ -z "$read1_path_col" ]; then
    echo "Error: 'read1_path' column not found in the runsheet."
    exit 1
fi
if [ "$paired" == "TRUE" ]; then
    # Get the column number for "read2_path"
    read2_path_col=$(get_col_num "read2_path" "$runsheet_csv")
    if [ -z "$read2_path_col" ]; then
        echo "Error: 'read2_path' column not found in the runsheet."
        exit 1
    fi
else
    read2_path_col=""
fi

# Extract sample IDs from read1_path column using raw_r1_suffix
# Logic: get basenames, replace suffix with ""
sample_ids_r1=$(awk -F, 'NR > 1 {split($'"$read1_path_col"', arr, "/"); name=arr[length(arr)]; gsub("'"$raw_r1_suffix"'","",name); print name}' $runsheet_csv | sort -u)

# If paired end, extract sample IDs from read2_path and compare them with those extracted from read1_path
if [ -n "$raw_r2_suffix" ]; then
    # If raw_r2_suffix is not empty, extract sample IDs using it
    sample_ids_r2=$(awk -v r2="$raw_r2_suffix" -v col="$read2_path_col" -F, 'NR > 1 {split($col, arr, "/"); name=arr[length(arr)]; gsub(r2,"",name); print name}' $runsheet_csv | sort -u)

    # Compare the two sample ID lists directly
    diff_lines=$(paste <(echo "$sample_ids_r1") <(echo "$sample_ids_r2") | awk '$1!=$2{print NR}')

    if [ -n "$diff_lines" ]; then
        echo "$diff_lines" | while read -r line_num; do
            echo "Error: Line $line_num of the runsheet has mismatching read filenames."
            echo "Ensure forward and reverse read names differ only in their suffixes."
        done
        exit 1
    fi
fi

# Write the sample IDs to the output file
echo "$sample_ids_r1" > $sample_ids_file

if [ -z "$sample_ids_r1" ]; then
    echo "Warning: $sample_ids_file is empty."
else
    echo "$sample_ids_file was successfully created."
fi



f_primer=""

# Get the column number for "F_primer"
f_primer_col=$(get_col_num "F_Primer" "$runsheet_csv")
if [ -z "$f_primer_col" ]; then
    echo "Error: 'F_Primer' column not found in the runsheet."
    exit 1
fi
# Extract F_primer values and take the unique value
f_primer_value=$(awk -v col="$f_primer_col" -F, 'NR > 1 {print $col}' $runsheet_csv | uniq)
# Check if the f_primer_value is a csv, tsv, or text file
if [[ "$f_primer_value" == *".csv"* ]] || [[ "$f_primer_value" == *".tsv"* ]] || [[ "$f_primer_value" == *".txt"* ]]; then
    f_primer="SET THIS MANUALLY FROM THE PRIMERS FILE"
    echo ""
    echo "Warning: The F_Primer column references a file."
    echo "Set 'f_primer' in config.yaml using the primers file."
else
    f_primer=$f_primer_value
fi
# Extract reverse primers and create the linked primers, also echo CSV primers info
if [ "$paired" == "TRUE" ]; then
    # Get the column number for "R_primer"
    r_primer_col=$(get_col_num "R_Primer" "$runsheet_csv")
    if [ -z "$r_primer_col" ]; then
        echo "Error: 'R_Primer' column not found in the runsheet."
        exit 1
    fi
    # Extract R_primer values and take the unique value
    r_primer_value=$(awk -v col="$r_primer_col" -F, 'NR > 1 {print $col}' $runsheet_csv | uniq)

    # Check if the r_primer_value is a csv, tsv, or text file
    if [[ "$r_primer_value" == *".csv"* ]] || [[ "$r_primer_value" == *".tsv"* ]] || [[ "$r_primer_value" == *".txt"* ]]; then
        r_primer="SET THIS MANUALLY FROM THE PRIMERS FILE"
        echo "Warning: The R_Primer column references a file."
        echo "Set 'r_primer' in config.yaml using the primers file."
        echo "For linked primers in the YAML:"
        echo "Use sequence & its reverse complement:"
        echo "F_linked_primer: \"^F_SEQ...REVERSE_OF_R_SEQ\""
        echo "R_linked_primer: \"^R_SEQ...REVERSE_OF_F_SEQ\""
        echo ""
    else
        r_primer=$r_primer_value

        # Also get linked primers if paired, using reverse complement of other primer:
        f_linked_primer="^${f_primer}...$(echo "$r_primer" | rev | tr 'ACGTMRWSYKVHDBN' 'TGCAMKWSRMBDHVNY')"
        r_linked_primer="^${r_primer}...$(echo "$f_primer" | rev | tr 'ACGTMRWSYKVHDBN' 'TGCAMKWSRMBDHVNY')"
    fi
else
    r_primer=""
    f_linked_primer=""
    r_linked_primer=""
fi

# Get the column number for "Parameter Value[Library Selection]"
target_region_col=$(get_col_num "Parameter Value[Library Selection]" "$runsheet_csv")
if [ -z "$target_region_col" ]; then
    echo "Error: 'Parameter Value[Library Selection]' column not found in the runsheet."
    exit 1
fi
# Extract target_region values and take the unique value
target_region=$(awk -v col="$target_region_col" -F, 'NR > 1 {print $col}' $runsheet_csv | uniq)
# Check if target_region is neither 16S nor ITS
if [[ "$target_region" != "16S" && "$target_region" != "ITS" ]]; then
    echo "Warning: Neither '16S' nor 'ITS' found in the runsheet. Defaulting to '16S'."
    target_region="16S"
    echo "Set target_region in config.yaml to '16S' or 'ITS' to change this default"
    echo "value before running the Snakemake workflow."
    echo ""
fi

# Make the output directory if it doesn't exist
mkdir -p "$output_dir"

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
    \"$runsheet_csv\"
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
    $min_trimmed_read_length


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
echo "config.yaml was successfully created."