#!/usr/bin/env Rscript

rm(list = ls())
library(stringr)
library(readr)
library(MSstats)

# Get root directory, assay suffix, runsheet path, fragpipe manifest path, and msstats_csv path from command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 0) {
    rootDir <- args[1]
} else {
    rootDir <- getwd()
}

assay_suffix <- if (length(args) > 1) args[2] else ""
runsheet_path <- args[3]
fragpipe_manifest_path <- args[4]
msstats_csv_path <- args[5]

# Ensure rootDir ends with /
if (!grepl("/$", rootDir)) {
    rootDir <- str_c(rootDir, "/")
}

print(str_c("Using IonQuant's result from ", rootDir))

# Read MSstats.csv file.
raw <- read_csv(msstats_csv_path, na = c("", "NA", "0"))
raw$ProteinName <- factor(raw$ProteinName)
raw$PeptideSequence <- factor(raw$PeptideSequence)

# Remove assay suffix from Run column entries
if (assay_suffix != "") {
    raw$Run <- gsub(assay_suffix, "", raw$Run, fixed = TRUE)
}

# Read runsheet and update Condition column based on Factor Value columns
runsheet <- read.csv(runsheet_path, stringsAsFactors = FALSE, check.names = FALSE)

# Debug: Print runsheet column names
print("Runsheet column names:")
print(names(runsheet))

# Identify Factor Value columns (assume at least 1 exists)
factor_cols <- names(runsheet)[grepl("^Factor Value\\[", names(runsheet))]
print(str_c("Found Factor Value columns: ", paste(factor_cols, collapse = ", ")))

# Create condition string by joining Factor Value columns with " & "
if (length(factor_cols) > 0) {
    condition_raw <- apply(runsheet[, factor_cols, drop = FALSE], 1, 
                           function(x) paste(x[!is.na(x) & x != ""], collapse = " & "))
} else {
    stop("No Factor Value columns found in runsheet!")
}

# Apply safe name handling: use make.names with BLOCKER_ prefix to handle names starting with numbers
# Then strip BLOCKER_ prefix to preserve original format
condition_safe <- sub("^BLOCKER_", "", make.names(str_c("BLOCKER_", condition_raw)))

# Create lookup vector: Sample Name -> Condition (safe name)
condition_lookup <- setNames(condition_safe, runsheet$`Sample Name`)

# Debug: Print condition mapping
print("Condition mapping from runsheet:")
print(data.frame(Sample_Name = runsheet$`Sample Name`, Condition_Raw = condition_raw, Condition_Safe = condition_safe, stringsAsFactors = FALSE))

# Match Run to Sample Name and replace Condition
raw$Condition <- condition_lookup[raw$Run]

# Debug: Print unique Condition values in raw data
print("Unique Condition values in raw data after matching:")
print(table(raw$Condition, useNA = "always"))

# Set BioReplicate from Source Name in runsheet
# Unique Source Names are naturally sorted (numeric-alphabetic) and assigned sequential integers
if ("Source Name" %in% names(runsheet)) {
    # Create mapping from unique Source Names using natural sort
    # Natural sort: handles numeric ordering (e.g., "Mouse_1", "Mouse_2", "Mouse_10" -> 1, 2, 3)
    source_names <- runsheet$`Source Name`
    unique_sources <- unique(source_names)
    
    # Split strings into parts (text and numbers) and pad numbers for proper numeric sorting
    natural_sort_key <- function(x) {
        parts <- strsplit(x, "(?<=\\D)(?=\\d)|(?<=\\d)(?=\\D)", perl=TRUE)[[1]]
        key <- sapply(parts, function(p) {
            num <- suppressWarnings(as.numeric(p))
            if (!is.na(num)) {
                sprintf("%020d", num)  # Pad numbers to 20 digits for sorting
            } else {
                p
            }
        })
        paste(key, collapse="")
    }
    unique_sources_sorted <- unique_sources[order(sapply(unique_sources, natural_sort_key))]
    source_to_num <- setNames(seq_along(unique_sources_sorted), unique_sources_sorted)
    bioreplicate_nums <- source_to_num[source_names]
    
    bioreplicate_lookup <- setNames(bioreplicate_nums, runsheet$`Sample Name`)
    raw$BioReplicate <- bioreplicate_lookup[raw$Run]
} else {
    stop("ERROR: 'Source Name' column required in runsheet for MSstats BioReplicate")
}

# Write processed CSV as intermediate file
write.csv(raw, "msstats_input.csv", row.names = FALSE)

# Change root directory for MSstats
print(str_c("Root DIR: ", rootDir))
setwd(rootDir)

#https://fragpipe.nesvilab.org/docs/tutorial_msstats.html
# Processing the data using MSstats
processedData <- dataProcess(raw, logTrans = 10)

# Generate all pairwise comparisons between conditions (using combn like amplicon workflow)
# Sort conditions first to ensure consistent ordering regardless of runsheet order
conditions <- unique(raw$Condition)
conditions <- conditions[!is.na(conditions)]
conditions <- sort(conditions)  # Sort to ensure consistent comparison order
if (length(conditions) > 1) {
    # Use combn to generate pairwise combinations (sorted order ensures consistency)
    contrast.names <- combn(conditions, 2)
    
    # Create contrast matrix for MSstats
    n_comparisons <- ncol(contrast.names)
    comparison <- matrix(0, nrow = n_comparisons, ncol = length(conditions))
    colnames(comparison) <- conditions
    rownames(comparison) <- paste(contrast.names[1,], contrast.names[2,], sep = "_v_")
    
    # Fill contrast matrix: +1 for first condition, -1 for second condition
    for (i in 1:n_comparisons) {
        comparison[i, contrast.names[1, i]] <- 1
        comparison[i, contrast.names[2, i]] <- -1
    }
    
    # Perform group comparison
    comparisonResults <- groupComparison(contrast.matrix = comparison, data = processedData)
    
    # Write comparison results - separate file for each comparison
    comparison_df <- comparisonResults$ComparisonResult
    for (comp_name in rownames(comparison)) {
        # Convert comparison name to filename-safe format
        # Replace spaces and special chars (except _v_) with single underscore, then collapse multiple underscores
        comp_safe <- gsub(" & ", "_", comp_name, fixed = TRUE)  # Replace " & " with single underscore
        comp_safe <- gsub("[^A-Za-z0-9_]", "_", comp_safe)  # Replace other special chars
        comp_safe <- gsub("_{2,}", "_", comp_safe)  # Collapse multiple underscores to single
        comp_safe <- gsub("^_|_$", "", comp_safe)  # Remove leading/trailing underscores
        comp_safe <- tolower(comp_safe)
        filename <- str_c("msstats_comparison_", comp_safe, ".csv")
        comp_data <- comparison_df[comparison_df$Label == comp_name, ]
        write.csv(comp_data, filename, row.names = FALSE)
    }
    
    # Write all comparisons combined
    write.csv(comparison_df, "msstats_comparison_all.csv", row.names = FALSE)
    
    # Write contrast matrix for reference
    write.csv(comparison, "msstats_contrasts.csv", row.names = TRUE)
}

