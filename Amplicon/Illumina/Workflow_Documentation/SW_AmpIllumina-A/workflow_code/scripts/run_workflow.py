import argparse
import subprocess
import os
import sys
import tempfile
import re 
import shutil
import pandas as pd
#import pandera as pa
import requests 

####################
## 1.  For OSD ARG #
####################
# 1. Process the OSD arg to proper format
# 2. Download the ISA file
# 3. Convert to runsheet(s)
# 4. Select which runsheet to use

########################
## 1. For runsheet arg #
########################
# 1. Select which runsheet to use

##########################
## 2. Neutral flow after #
##########################
# 1. Validate schema of runsheet
# 2. Check if read_paths are URLs, prompt for download
# 3. Create config.yaml and unique-sample-IDs.txt
# 4. If --run is used: run the workflow

# Process OSD arg: if numeric, append OSD-, if OSD-# or GLDS-#, leave it
def process_osd_argument(osd_arg):
    # Check if the argument is just numeric
    if osd_arg.isdigit():
        return f"OSD-{osd_arg}"
    # Check if it's already in the correct format (OSD-numeric or GLDS-numeric)
    elif re.match(r'^(OSD|GLDS)-\d+$', osd_arg):
        return osd_arg
    else:
        print("Invalid format for --OSD argument. Use 'numeric', 'OSD-numeric', or 'GLDS-numeric'.")
        sys.exit(1)

# Run dpt-get-isa-archive in a temp folder, move it back to cd, return the filename
def download_isa_archive(accession_number):
    with tempfile.TemporaryDirectory() as temp_dir:
        try:
            # Run the command in the temporary directory
            subprocess.run(
                ["dpt-get-isa-archive", "--accession", str(accession_number)],
                check=True,
                text=True,
                cwd=temp_dir
            )

            # Find the downloaded zip file in the temp directory
            downloaded_files = [f for f in os.listdir(temp_dir) if f.endswith('.zip')]
            if not downloaded_files:
                print("No ISA archive file was downloaded.", file=sys.stderr)
                return None

            # Assuming there's only one file, get its name
            downloaded_file = downloaded_files[0]

            # Move the file back to the current directory
            shutil.move(os.path.join(temp_dir, downloaded_file), downloaded_file)

            return downloaded_file

        except subprocess.CalledProcessError as e:
            print("An error occurred while downloading ISA archive.", file=sys.stderr)
            sys.exit(1)

# Run dpt-isa-to-runsheet in a temp folder, move runsheet(s) back to cd, return list of runsheet(s)
def convert_isa_to_runsheet(accession_number, isa_zip):
    with tempfile.TemporaryDirectory() as temp_dir:
        # Copy the ISA archive to the temporary directory
        temp_isa_zip_path = shutil.copy(isa_zip, temp_dir)

        try:
            # Run the dpt-isa-to-runsheet command in the temporary directory
            subprocess.run(
                ["dpt-isa-to-runsheet", "--accession", accession_number, "--config-type", "amplicon", "--config-version", "Latest", "--isa-archive", os.path.basename(temp_isa_zip_path)],
                check=True,
                cwd=temp_dir,
                stdout=sys.stdout,
                stderr=sys.stderr
            )

            # Get the list of created files in the temp directory
            created_files = [f for f in os.listdir(temp_dir) if os.path.isfile(os.path.join(temp_dir, f)) and f != os.path.basename(temp_isa_zip_path)]

            # Move the created files back to the current directory
            moved_files = []
            for file in created_files:
                shutil.move(os.path.join(temp_dir, file), file)
                moved_files.append(file)

            return moved_files

        except subprocess.CalledProcessError as e:
            print("An error occurred while converting ISA archive to runsheet.", file=sys.stderr)
            sys.exit(1)

# Prompt the user to select a runsheet if more than 1 were created, else just return the single runsheet
def handle_runsheet_selection(runsheet_files):
    if len(runsheet_files) == 1:
        # Automatically use the single runsheet file
        selected_runsheet = runsheet_files[0]
        print(f"Using runsheet: {selected_runsheet}")
        return selected_runsheet
    elif len(runsheet_files) > 1:
        return None





        # Logic for prompting the user which runsheet to use
        # Prompt the user to select a runsheet file
        # print("Select a runsheet to use:")
        # for idx, file in enumerate(runsheet_files, start=1):
        #     print(f"[{idx}] {file}")

        # # Wait for user input and validate it
        # while True:
        #     try:
        #         selection = int(input("Enter the number of the runsheet: "))
        #         if 1 <= selection <= len(runsheet_files):
        #             selected_runsheet = runsheet_files[selection - 1]
        #             print(f"Selected {selected_runsheet}.")
        #             return selected_runsheet
        #         else:
        #             print("Invalid selection. Please enter a number from the list.")
        #     except ValueError:
        #         print("Invalid input. Please enter a number.")


# Neutral functions

# Verify that the runsheet follows the schema
# def validate_runsheet_schema(csv_file):
#     schema = pa.DataFrameSchema({
#         "Sample Name": pa.Column(str),
#         "paired_end": pa.Column(bool),
#         "read1_path": pa.Column(str),
#         "read2_path": pa.Column(str, nullable=True),  # Optional if paired_end is False
#         "F_Primer": pa.Column(str),
#         "R_Primer": pa.Column(str, nullable=True),  # Optional if paired_end is False
#         "raw_R1_suffix": pa.Column(str),
#         "raw_R2_suffix": pa.Column(str, nullable=True),  # Optional if paired_end is False
#         "groups": pa.Column(str)
#     })

#     # Load CSV file into a pandas DataFrame
#     df = pd.read_csv(csv_file)

#     # Validate DataFrame against the schema
#     try:
#         validated_df = schema.validate(df)
#         print(f"CSV file '{csv_file}' successfully validated.")
#         return validated_df
#     except pa.errors.SchemaError as e:
#         print(f"Schema validation error: {e}", file=sys.stderr)
#         return None

def check_runsheet_read_paths(runsheet_df):
    # Check if a string is a URL / genelab URL
    def is_url(s):
        return "http://" in s or "https://" in s or "genelab-data.ndc.nasa.gov" in s


    # Check if 'read2_path' column exists 
    paired_end = runsheet_df['paired_end'].eq(True).all()

    # Check the first row to determine if the paths are URLs or local paths
    first_row = runsheet_df.iloc[0]

    uses_url = is_url(first_row['read1_path'])
    if uses_url:
        print("Runsheet references URLs.")
    else:
        print("Runsheet references local read files.")

    return uses_url

def sample_IDs_from_local(runsheet_df, output_file='unique-sample-IDs.txt'):
    # Check if the DataFrame is paired-end
    paired_end = runsheet_df['paired_end'].eq(True).all()

    with open(output_file, 'w') as file:
        for index, row in runsheet_df.iterrows():
            # Extract base names minus the suffixes
            base_read1 = os.path.basename(row['read1_path']).replace(row['raw_R1_suffix'], '')

            if paired_end:
                base_read2 = os.path.basename(row['read2_path']).replace(row['raw_R2_suffix'], '')
                # Check if base names match for paired-end data, necessary for snakemake arg expansion
                if base_read1 != base_read2:
                    print(f"Mismatch in sample IDs in row {index}: {base_read1} vs {base_read2}")
                    sys.exit(1)
            
            # Write the base name to the file
            file.write(f"{base_read1}\n")
    
    print(f"Unique sample IDs written to {output_file}")

def handle_url_downloads(runsheet_df, output_file='unique-sample-IDs.txt'):
    print("Downloading read files...")
    # Check if the DataFrame is paired-end
    paired_end = runsheet_df['paired_end'].eq(True).all()
    # Write 'Sample Name' into unique-sample-IDs.txt
    with open(output_file, 'w') as file:
        for sample_name in runsheet_df['Sample Name']:
            file.write(sample_name + '\n')

    # Create ./raw_reads/ directory if it does not exist
    raw_reads_dir = os.path.abspath('./raw_reads/')
    if not os.path.exists(raw_reads_dir):
        os.makedirs(raw_reads_dir)

    # Initialize count for skipped downloads
    skipped_downloads_count = 0
    # Iterate over each row and download files if they don't exist
    for _, row in runsheet_df.iterrows():
        sample_id = row['Sample Name']
        read1_path = os.path.join(raw_reads_dir, sample_id + row['raw_R1_suffix'])
        read2_path = os.path.join(raw_reads_dir, sample_id + row['raw_R2_suffix']) if paired_end else None

        # Download Read 1 if it doesn't exist
        if not os.path.exists(read1_path):
            download_url_to_file(row['read1_path'], read1_path)
        else:
            skipped_downloads_count += 1

        # Download Read 2 if it doesn't exist and if paired_end
        if paired_end and read2_path and not os.path.exists(read2_path):
            download_url_to_file(row['read2_path'], read2_path)
        elif paired_end and read2_path:
            skipped_downloads_count += 1

    # Print the number of skipped downloads
    if skipped_downloads_count > 0:
        print(f"{skipped_downloads_count} read files were present. Skipped downloads for those.")

def download_url_to_file(url, file_path):
    response = requests.get(url, stream=True)
    if response.status_code == 200:
        with open(file_path, 'wb') as file:
            shutil.copyfileobj(response.raw, file)
    else:
        print(f"Failed to download file from {url}")


def reverse_complement(seq):
    complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 
                  'R': 'Y', 'Y': 'R', 'S': 'S', 'W': 'W', 
                  'K': 'M', 'M': 'K', 'B': 'V', 'V': 'B', 
                  'D': 'H', 'H': 'D', 'N': 'N'}
    return ''.join(complement.get(base, base) for base in reversed(seq))

def create_config_yaml(runsheet_file, runsheet_df, min_trimmed_length, uses_urls, output_dir):
    # Extract necessary variables from runsheet_df
    data_type = "PE" if runsheet_df['paired_end'].eq(True).all() else "SE"
    raw_R1_suffix = runsheet_df['raw_R1_suffix'].unique()[0]
    raw_R2_suffix = runsheet_df['raw_R2_suffix'].unique()[0] if data_type == "PE" else ""
    f_primer = runsheet_df['F_Primer'].unique()[0]
    r_primer = runsheet_df['R_Primer'].unique()[0] if data_type == "PE" else ""
    target_region = runsheet_df['Parameter Value[Library Selection]'].unique()[0]

    # Determine raw_reads_directory
    if uses_urls:
        raw_reads_directory = os.path.abspath('./raw_reads/') + '/'
    else:
        read1_path_dir = os.path.dirname(runsheet_df['read1_path'].iloc[0])
        raw_reads_directory = os.path.abspath(read1_path_dir) + '/' if read1_path_dir else "./"


    # Other default values
    output_dir = os.path.abspath(output_dir) + '/'
    trim_primers = "TRUE"
    primers_linked = "TRUE"

    f_linked_primer = f"^{f_primer}...{reverse_complement(r_primer)}"
    r_linked_primer = f"^{r_primer}...{reverse_complement(f_primer)}"

    # Make output_dir if it doesn't exist 
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Write to config.yaml
    with open('config.yaml', 'w') as file:
        file.write("############################################################################################\n")
        file.write("## Configuration file for GeneLab Illumina amplicon processing workflow                   ##\n")
        file.write("## Developed by Michael D. Lee (Mike.Lee@nasa.gov)                                        ##\n")
        file.write("############################################################################################\n\n")

        file.write("############################################################\n")
        file.write("##################### VARIABLES TO SET #####################\n")
        file.write("############################################################\n\n")

        file.write("###########################################################################\n")
        file.write("##### These need to match what is specific to our system and our data #####\n")
        file.write("###########################################################################\n\n")

        file.write("## Path to runsheet:\n")
        file.write(f"runsheet:\n    \"{os.path.abspath(runsheet_file)}\"\n\n")

        file.write("## Set to \"PE\" for paired-end, \"SE\" for single-end.\n")
        file.write(f"data_type:\n    \"{data_type}\"\n\n")

        file.write("## single-column file with unique sample identifiers:\n")
        file.write("sample_info_file:\n    \"unique-sample-IDs.txt\"\n\n")

        file.write("## input reads directory (can be relative to workflow directory, or needs to be full path):\n")
        file.write(f"raw_reads_dir:\n    \"{raw_reads_directory}\"\n\n")

        file.write("## raw read suffixes:\n")
        file.write("  # e.g. for paired-end data, Sample-1_R1_raw.fastq.gz would be _R1_raw.fastq.gz for 'raw_R1_suffix' below\n")
        file.write("  # e.g. if single-end, Sample-1.fastq.gz would be .fastq.gz for 'raw_R1_suffix' below, and 'raw_R2_suffix' won't be used\n")
        file.write(f"raw_R1_suffix:\n    \"{raw_R1_suffix}\"\n")
        file.write(f"raw_R2_suffix:\n    \"{raw_R2_suffix}\"\n\n")

        file.write("## if we are trimming primers or not (\"TRUE\", or \"FALSE\")\n")
        file.write(f"trim_primers:\n    \"{trim_primers}\"\n\n")

        file.write("## primer sequences if we are trimming them (include anchoring symbols, e.g. '^', as needed, see: https://cutadapt.readthedocs.io/en/stable/guide.html#adapter-types)\n")
        file.write(f"F_primer:\n    \"^{f_primer}\"\n")
        file.write(f"R_primer:\n    \"^{r_primer}\"\n\n")

        # For linked primers
        file.write("## should cutadapt treat these as linked primers? (https://cutadapt.readthedocs.io/en/stable/recipes.html#trimming-amplicon-primers-from-paired-end-reads)\n")
        file.write(f"primers_linked:\n    \"{primers_linked}\"\n\n")
        file.write("## if primers are linked, we need to provide them as below, where the second half, following three periods, is the other primer reverse-complemented\n")
        file.write(f"  # (can reverse complement while retaining ambiguous bases at this site: http://arep.med.harvard.edu/labgc/adnan/projects/Utilities/revcomp.html)\n")
        file.write(f"  # include anchoring symbols, e.g. '^', as needed, see: https://cutadapt.readthedocs.io/en/stable/guide.html#adapter-types\n")
        file.write(f"F_linked_primer:\n    \"{f_linked_primer}\"\n")
        file.write(f"R_linked_primer:\n    \"{r_linked_primer}\"\n\n")

        file.write("## discard untrimmed, sets the \"--discard-untrimmed\" option if TRUE\n")
        file.write("discard_untrimmed:\n    \"TRUE\"\n\n")

        file.write("## target region (16S or ITS acceptable; determines which reference database is used for taxonomic classification)\n")
        file.write(f"target_region:\n    \"{target_region}\"\n\n")

        file.write("## concatenate only with dada2 instead of merging paired reads if TRUE\n")
        file.write("  # this is typically used with primers like 515-926, that captured 18S fragments that are typically too long to merge\n")
        file.write("  # note that 16S and 18S should have been separated already prior to running this workflow\n")
        file.write("  # this should likely be left as FALSE for any option other than \"18S\" above\n\n")
        file.write("concatenate_reads_only:\n    \"FALSE\"\n\n")

        file.write("## values to be passed to dada2's filterAndTrim() function:\n")
        file.write("left_trunc:\n    0\n")
        file.write("right_trunc:\n    0\n")
        file.write("left_maxEE:\n    1\n")
        file.write("right_maxEE:\n    1\n\n")

        file.write("## minimum length threshold for cutadapt\n")
        file.write(f"pcutadapt_len:\n    {min_trimmed_length}\n\n")

        file.write("######################################################################\n")
        file.write("##### The rest only need to be altered if we want to change them #####\n")
        file.write("######################################################################\n\n")

        file.write("## filename suffixes\n")
        file.write("primer_trimmed_R1_suffix:\n    \"_R1_trimmed.fastq.gz\"\n")
        file.write("primer_trimmed_R2_suffix:\n    \"_R2_trimmed.fastq.gz\"\n\n")

        file.write("filtered_R1_suffix:\n    \"_R1_filtered.fastq.gz\"\n")
        file.write("filtered_R2_suffix:\n    \"_R2_filtered.fastq.gz\"\n\n")

        file.write("## output prefix (if needed to distinguish from multiple primer sets, leave as empty string if not, include connecting symbol if adding, e.g. \"ITS-\")\n")
        file.write("output_prefix:\n    \"\"\n\n")

        file.write("## output directories (all relative to processing directory, they will be created if needed)\n")
        file.write(f"fastqc_out_dir:\n    \"{output_dir}FastQC_Outputs/\"\n")
        file.write(f"trimmed_reads_dir:\n    \"{output_dir}Trimmed_Sequence_Data/\"\n")
        file.write(f"filtered_reads_dir:\n    \"{output_dir}Filtered_Sequence_Data/\"\n")
        file.write(f"final_outputs_dir:\n    \"{output_dir}Final_Outputs/\"\n\n")

        # For general info and example usage command
        file.write("############################################################\n")
        file.write("###################### GENERAL INFO ########################\n")
        file.write("############################################################\n")
        file.write("# Workflow is currently equipped to work with paired-end data only, and reads are expected to be gzipped\n\n")
        file.write("## example usage command ##\n")
        file.write("# snakemake --use-conda --conda-prefix ${CONDA_PREFIX}/envs -j 2 -p\n")
        file.write("# `--use-conda` – this specifies to use the conda environments included in the workflow\n")
        file.write("# `--conda-prefix` – this allows us to point to where the needed conda environments should be stored...\n")
        file.write("# `-j` – this lets us set how many jobs Snakemake should run concurrently...\n")
        file.write("# `-p` – specifies to print out each command being run to the screen\n\n")
        file.write("# See `snakemake -h` for more options and details.\n")
    print("config.yaml was successfully created.")

# Example usage
# create_config_yaml(runsheet_df, uses_urls)


def main():
    # Argument parser setup with short argument names and an automatic help option
    parser = argparse.ArgumentParser(
        description='Run workflow for GeneLab data processing.',
        add_help=True,
        usage='%(prog)s [options]'  # Custom usage message
    )
    
    parser.add_argument('-o', '--OSD',
                        help='Set up the Snakemake workflow for a GeneLab OSD dataset and pull necessary read files and metadata. Acceptable formats: ###, OSD-###, GLDS-###',
                        type=str)
    
    parser.add_argument('-r', '--runsheetPath',
                        help='Set up the Snakemake workflow using a specified runsheet file.',
                        type=str)
    
    parser.add_argument('-m', '--min_trimmed_length',
                    default=130,  # Default value
                    help='Minimum length of trimmed reads. For paired-end data: if one read gets filtered, both reads are discarded. Default: 130',
                    type=int)

    parser.add_argument('-d', '--outputDir',
                        default='./workflow_output/',  # Default value
                        help='Specify the output directory for the workflow. This argument is only used if -o or -r is specified. Default: ./workflow/',
                        type=str)
    
    parser.add_argument('-x', '--run',
                        nargs='?',
                        const="snakemake --use-conda --conda-prefix ${CONDA_PREFIX}/envs -j 2 -p",
                        type=str,
                        help='Execute the Snakemake workflow. Optionally provide a custom Snakemake command. Default: "snakemake --use-conda --conda-prefix ${CONDA_PREFIX}/envs -j 2 -p"')
    args = parser.parse_args()
    output_dir = args.outputDir
    min_trimmed_length = args.min_trimmed_length

    # If OSD is used, pull ISA metadata for the study, create and select the runsheet
    if args.OSD:
        accession_number = process_osd_argument(args.OSD)
        isa_zip = download_isa_archive(accession_number)
        if isa_zip:
            runsheet_files = convert_isa_to_runsheet(accession_number, isa_zip)
            if runsheet_files:
                runsheet_file = handle_runsheet_selection(runsheet_files)
                if runsheet_file is None:
                    sys.exit("This OSD dataset contains multiple assays for the indicated amplicon. This situation is not handled in the current workflow version. To process this dataset using the workflow, please use approach 2 and set up a manual runsheet containing the information from the metadata on OSDR, and point to the raw reads links for the specific assay you want to process.")
            else:
                print("No runsheet files were created.")
        else:
            print("No ISA archive was downloaded. Cannot proceed to runsheet conversion.", file=sys.stderr)
            sys.exit(1)
    
    # If a runsheet is specified, use that runsheet
    elif args.runsheetPath:
        runsheet_file = args.runsheetPath

    # Validate and load the runsheet if a file is specified
    # Create unique-sample-IDs.txt based on filenames or 'Sample Name' if URLs
    # Download files if necessary
    if args.OSD or args.runsheetPath:
        if runsheet_file:
            #runsheet_df = validate_runsheet_schema(runsheet_file)
            runsheet_df = pd.read_csv(runsheet_file)
            if runsheet_df is not None:
                uses_urls = check_runsheet_read_paths(runsheet_df)

                # Create the 'unique-sample-IDs.txt' file and download read files if necessary
                if uses_urls:
                    handle_url_downloads(runsheet_df, output_file='unique-sample-IDs.txt')
                else:
                    sample_IDs_from_local(runsheet_df, output_file='unique-sample-IDs.txt')

                # Create the config.yaml file
                create_config_yaml(runsheet_file, runsheet_df, min_trimmed_length, uses_urls, output_dir)
                print("Snakemake workflow setup is complete.")
            else:
                print("Failed to validate the runsheet file.", file=sys.stderr)
                sys.exit(1)
        else:
            print("No runsheet file specified.", file=sys.stderr)
            sys.exit(1)
    
    # Run the snakemake workflow if --run is used
    if args.run:
        snakemake_command = args.run if args.run is not None else "snakemake --use-conda --conda-prefix ${CONDA_PREFIX}/envs -j 2 -p"
        print(f"Running Snakemake command: {snakemake_command}")
        subprocess.run(snakemake_command, shell=True, check=True)
    




if __name__ == "__main__":
    main()