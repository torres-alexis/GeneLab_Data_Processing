#!/usr/bin/env python
import pandas as pd
import sys
import csv

# Input and output file paths from command line arguments
input_file = sys.argv[1]
output_file = sys.argv[2]

# Load the data
data = pd.read_csv(input_file, header=1, sep='\t')

# Set 'Geneid' as the index (row names) for the DataFrame
data.set_index('Geneid', inplace=True)

# Remove columns not needed for counting
data.drop(columns=["Chr", "Start", "End", "Strand", "Length"], inplace=True)

# Remove '.bam' from column names
data.columns = [col.replace('.bam', '') for col in data.columns]

# Convert the remaining data to numeric, assuming these are the count columns
data = data.apply(pd.to_numeric, errors='coerce')

# Count non-zero entries for each sample
non_zero_counts = (data > 0).sum().reset_index()
non_zero_counts.columns = ["", "Number of genes with non-zero counts"]

# Save the results to a CSV
non_zero_counts.to_csv(output_file, index=False, quotechar='"', quoting=csv.QUOTE_NONNUMERIC)