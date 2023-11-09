#!/bin/bash

# config.yaml format is based on https://github.com/nasa/GeneLab_Data_Processing/blob/master/Amplicon/Illumina/Workflow_Documentation/SW_AmpIllumina-A/workflow_code/config.yaml

# Initialize default values
raw_reads_directory="raw_reads/" 
output_dir="workflow_output/"
min_trimmed_read_length=1
download_counter_file=$(mktemp)

    echo ""
# Function to display usage
print_usage() {
    echo "Usage: $0"
    echo "This script must be executed from the 'workflow_code' directory."
    echo ""
    echo "Required arguments:"
    echo "  -r, --runsheet <path>    Path to runsheet CSV file, relative to the"
    echo "                           'workflow_code' directory or as an absolute path."
    
    echo "Optional arguments:"
    echo "  -o, --output <path>      Path to output directory, relative to the"
    echo "                           'workflow_code' directory or as an absolute path."
    echo "                           Default is 'workflow_output/'"
    echo "  -m, --min_length <num>   Minimum trimmed read length."
    echo "                           For paired end data: if one read gets filtered,"
    echo "                           both reads are discarded."
    echo "                           Default is 1"
    echo "  -h, --help               Display this help menu and exit."
    echo ""
    echo "Example:"
    echo "    $0 -r runsheet.csv -o workflow_output/ -m 1"
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
	-o|--output)
		output_dir="$1"
		shift
		;;
	-m|--min_length)
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
if [[ -z "$runsheet_csv" ]]; then
    echo "Error: Missing required arguments."
    echo ""
    print_usage
    exit 1
fi

echo "Runsheet: $runsheet_csv"
#echo "Raw Reads Directory: $raw_reads_directory"
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


check_for_links() {
    local col_num=$1
    if awk -v col="$col_num" -F, 'NR > 1 {print $col}' "$runsheet_csv" | grep -E -q '(http|https)://|genelab-data.ndc.nasa.gov'; then
        return 0 # True, URLs are detected
    else
        return 1 # False, URLs are not detected
    fi
}

# Function to set the boolean uses_links based on column checks
set_uses_links_flag() {
    local read1_col_num=$1
    local uses_links=false
    
    # Check for links in "read1_path" column
    if check_for_links "$read1_col_num"; then
        uses_links=true
    fi
    echo "$uses_links"
}

uses_links=$(set_uses_links_flag "$read1_path_col")

if [ "$uses_links" == "true" ]; then
    echo "Links are being used in the runsheet."
else
    echo "Local paths are being used in the runsheet."
fi

# Prompt function definition
prompt_user_for_download() {
    read -p "The $raw_reads_directory directory is not empty. Download read files anyway? [y/N]: " yn
    [[ $yn =~ ^[Yy]$ ]]
}

# Function to download a file and ensure a proper filename
download_file() {
    local url=$1
    local out_dir=$2
    local filename=$(echo "$url" | sed -n 's/.*file=\([^&]*\).*/\1/p')

    if [ -n "$filename" ]; then
        if wget -O "$out_dir/$filename" "$url" > /dev/null 2>&1; then
            echo "Downloaded $filename"
            echo "$filename" >> "$download_counter_file"
        else
            # If wget fails, use curl as a backup
            if curl -o "$out_dir/$filename" "$url" > /dev/null 2>&1; then
                echo "Downloaded $filename"
                echo "$filename" >> "$download_counter_file"
            else
                echo "Failed to download: $filename"
            fi
        fi
    else
        echo "Not a valid URL: $url"
    fi
}

# Function to process the runsheet and initiate downloads
process_and_download_reads() {
    local runsheet=$1
    local r1_col=$2
    local r2_col=$3
    local out_dir=$4

    # Read the file line by line starting from the second line
    tail -n +2 "$runsheet" | while IFS=, read -r -a fields; do
        local r1_path="${fields[r1_col-1]}"
        local r2_path="${fields[r2_col-1]}"

        # Download read 1 path
        download_file "$r1_path" "$out_dir" 

        # If paired is true, download read 2 path
        if [ ! -z "$r2_path" ]; then
            download_file "$r2_path" "$out_dir" 
        fi
    done
}

# Check if $raw_reads_directory is empty
if [ "$uses_links" == "true" ]; then
    mkdir -p "$raw_reads_directory"
    if [ -d "$raw_reads_directory" ] && [ "$(ls -A "$raw_reads_directory")" ]; then
        # Directory exists and is not empty, prompt the user
        if prompt_user_for_download; then
            echo "Downloading read files to $raw_reads_directory"
            process_and_download_reads "$runsheet_csv" "$read1_path_col" "$read2_path_col" "$raw_reads_directory"
            number_of_downloaded_files=$(awk 'END {print NR}' "$download_counter_file")
            echo "Downloaded $number_of_downloaded_files files."
            rm "$download_counter_file"
        else
            echo "Skipping download."
        fi
    else
        echo "Downloading files to $raw_reads_directory"
        process_and_download_reads "$runsheet_csv" "$read1_path_col" "$read2_path_col" "$raw_reads_directory"
        number_of_downloaded_files=$(awk 'END {print NR}' "$download_counter_file")
        echo "Downloaded $number_of_downloaded_files files."
        rm "$download_counter_file"
    fi
else
    echo "Local paths are being used."
fi


# Set raw reads directory based on the path in read1_path if the runsheet is not using download links for paths
if [ "$uses_links" == "true" ]; then
    :
else
    first_read1_path=$(awk -v col="$read1_path_col" -F, 'NR==2 {print $col; exit}' "$runsheet_csv")
    raw_reads_directory="$(dirname "$first_read1_path")/"
fi


if [ "$uses_links" == "true" ]; then
    # Logic for URLs using sed to extract the filename
    sample_ids_r1=$(awk -F, -v col="$read1_path_col" 'NR > 1 {
        print $col
      }' "$runsheet_csv" | sed -n 's/.*[?&]file=\([^&]*\).*/\1/p' | sed "s/$raw_r1_suffix\$//" | sort -u)
else
    # Logic for local paths
    sample_ids_r1=$(awk -F, -v col="$read1_path_col" -v r1_suffix="$raw_r1_suffix" 'NR > 1 {
        split($col, arr, "/")
        name=arr[length(arr)]
        gsub(r1_suffix, "", name)
        print name
      }' "$runsheet_csv" | sort -u)
fi

# If paired end, extract sample IDs from read2_path and compare them with those extracted from read1_path
if [ -n "$raw_r2_suffix" ]; then
    if [ "$uses_links" == "true" ]; then
        # Logic for URLs using sed to extract the filename
        sample_ids_r2=$(awk -F, -v col="$read2_path_col" 'NR > 1 {
            print $col
          }' "$runsheet_csv" | sed -n 's/.*[?&]file=\([^&]*\).*/\1/p' | sed "s/$raw_r2_suffix\$//" | sort -u)
    else
        # Logic for local paths
        sample_ids_r2=$(awk -F, -v col="$read2_path_col" -v r2_suffix="$raw_r2_suffix" 'NR > 1 {
            split($col, arr, "/")
            name=arr[length(arr)]
            gsub(r2_suffix, "", name)
            print name
          }' "$runsheet_csv" | sort -u)
    fi
    
   
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