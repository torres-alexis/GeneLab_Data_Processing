#!/usr/bin/env Rscript

rm(list = ls())
library(stringr)
library(MSstats)

# Get root directory, experiment_annotation, msstats_csv, and optional assay_suffix 
args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 0) {
    rootDir <- args[1]
} else {
    rootDir <- getwd()
}

experiment_annotation_path <- args[2]
msstats_csv_path <- args[3]
assay_suffix <- if (length(args) > 3) args[4] else ""

# Ensure rootDir ends with /
if (!grepl("/$", rootDir)) {
    rootDir <- str_c(rootDir, "/")
}

print(str_c("Using IonQuant's result from ", rootDir))

# Read MSstats.csv file.
raw <- read.csv(msstats_csv_path, na.strings = c("", "NA", "0"), stringsAsFactors = FALSE)
raw$ProteinName <- factor(raw$ProteinName)
raw$PeptideSequence <- factor(raw$PeptideSequence)

# Remove assay suffix from Run column entries
if (assay_suffix != "") {
    raw$Run <- gsub(assay_suffix, "", raw$Run, fixed = TRUE)
}

# Read experiment_annotation (same source as fp_analyst; condition/condition_label from runsheet_to_fp_metadata)
anno <- read.table(experiment_annotation_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
colnames(anno) <- tolower(colnames(anno))
if (!"condition" %in% colnames(anno)) {
    stop("experiment_annotation must have condition column")
}
# MSstats Run column: file basename (from FragPipe manifest). Map via file column.
if (!"file" %in% colnames(anno)) stop("experiment_annotation must have file column")
Run_for_lookup <- gsub("\\.mzML$", "", raw$Run, ignore.case = TRUE)
file_basename <- gsub("\\.mzML$", "", anno$file, ignore.case = TRUE)
condition_lookup <- setNames(anno$condition, file_basename)
raw$Condition <- condition_lookup[Run_for_lookup]
if (all(is.na(raw$Condition))) {
    stop("MSstats Run values do not match experiment_annotation file. Run sample: ", paste(head(unique(raw$Run), 5), collapse = ", "))
}
safe_to_label <- if ("condition_label" %in% colnames(anno) && all(nzchar(trimws(anno$condition_label)))) {
    u <- unique(anno[, c("condition", "condition_label")])
    setNames(u$condition_label, u$condition)
} else {
    setNames(anno$condition, anno$condition)
}

print("Condition mapping from experiment_annotation:")
print(anno[, intersect(c("file", "condition", "condition_label"), colnames(anno)), drop = FALSE])

# Debug: Print unique Condition values in raw data
print("Unique Condition values in raw data after matching:")
print(table(raw$Condition, useNA = "always"))

# Write processed CSV as intermediate file
write.csv(raw, "msstats_input.csv", row.names = FALSE)

# Change root directory for MSstats
print(str_c("Root DIR: ", rootDir))
setwd(rootDir)

#https://fragpipe.nesvilab.org/docs/tutorial_msstats.html
# Processing the data using MSstats
processedData <- dataProcess(raw, logTrans = 10)

# Generate all pairwise comparisons between conditions
# Sort conditions first to ensure consistent ordering regardless of runsheet order
conditions <- unique(raw$Condition)
conditions <- conditions[!is.na(conditions)]
conditions <- sort(conditions)  # Sort to ensure consistent comparison order

if (length(conditions) > 1) {
    # Pairwise combinations: numerator = second, denominator = first
    contrast.names <- combn(conditions, 2)
    
    # Contrast matrix: numerator = second, denominator = first
    n_comparisons <- ncol(contrast.names)
    comparison <- matrix(0, nrow = n_comparisons, ncol = length(conditions))
    colnames(comparison) <- conditions
    rownames(comparison) <- paste(contrast.names[2,], contrast.names[1,], sep = "_v_")
    
    for (i in 1:n_comparisons) {
        comparison[i, contrast.names[2, i]] <- 1   # numerator
        comparison[i, contrast.names[1, i]] <- -1  # denominator
    }
    
    # Perform group comparison
    comparisonResults <- groupComparison(contrast.matrix = comparison, data = processedData)
    
    format_label <- function(cond1, cond2) paste0("(", cond1, ")v(", cond2, ")")
    comparison_df <- comparisonResults$ComparisonResult
    for (i in 1:n_comparisons) {
        c1 <- contrast.names[1, i]
        c2 <- contrast.names[2, i]
        r1 <- safe_to_label[c1]
        r2 <- safe_to_label[c2]
        if (is.na(r1)) r1 <- c1
        if (is.na(r2)) r2 <- c2
        lbl_new <- format_label(r2, r1)  # (numerator)v(denominator) = (c2)v(c1)
        comparison_df$Label[comparison_df$Label == rownames(comparison)[i]] <- lbl_new
    }
    
    # All MSstats pairwise comparisons
    write.csv(comparison_df, str_c("msstats_comparison", assay_suffix, ".csv"), row.names = FALSE)

    ## per-contrast files, redundant with rows in msstats_comparison.csv
    # for (i in 1:n_comparisons) {
    #     c1 <- contrast.names[1, i]
    #     c2 <- contrast.names[2, i]
    #     r1 <- safe_to_label[c1]
    #     r2 <- safe_to_label[c2]
    #     if (is.na(r1)) r1 <- c1
    #     if (is.na(r2)) r2 <- c2
    #     comp_name_raw <- format_label(r2, r1)
    #     comp_safe <- gsub(" & ", "_", comp_name_raw, fixed = TRUE)
    #     comp_safe <- gsub("[^A-Za-z0-9_]", "_", comp_safe)
    #     comp_safe <- gsub("_{2,}", "_", comp_safe)
    #     comp_safe <- gsub("^_|_$", "", comp_safe)
    #     comp_safe <- tolower(comp_safe)
    #     filename <- str_c("msstats_comparison_", comp_safe, assay_suffix, ".csv")
    #     comp_data <- comparison_df[comparison_df$Label == comp_name_raw, ]
    #     write.csv(comp_data, filename, row.names = FALSE)
    # }
    # write.csv(comparison_df, str_c("msstats_comparison_all", assay_suffix, ".csv"), row.names = FALSE)
    
    contrasts_df <- data.frame(row.names = c("1", "2"))
    for (i in 1:n_comparisons) {
        c1 <- contrast.names[1, i]
        c2 <- contrast.names[2, i]
        r1 <- if (is.na(safe_to_label[c1])) c1 else safe_to_label[c1]
        r2 <- if (is.na(safe_to_label[c2])) c2 else safe_to_label[c2]
        col_name <- format_label(r2, r1)
        contrasts_df[[col_name]] <- c(c2, c1)
    }
    write.csv(contrasts_df, str_c("msstats_contrasts", assay_suffix, ".csv"), row.names = TRUE)
}

