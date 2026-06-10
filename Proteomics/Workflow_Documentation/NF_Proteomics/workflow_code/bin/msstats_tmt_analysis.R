#!/usr/bin/env Rscript
# MSstatsTMT for FragPipe TMT (Philosopher msstats.csv + pipeline annotation).

rm(list = ls())
library(MSstatsTMT)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: msstats_tmt_analysis.R <rootDir> <MSstatsTMT_annotation.csv> <msstats_dir_or_file> [assay_suffix]")
}
rootDir <- args[1]
annotation_path <- args[2]
msstats_path <- args[3]
assay_suffix <- if (length(args) > 3) args[4] else ""

if (!grepl("/$", rootDir)) rootDir <- paste0(rootDir, "/")

annotation <- read.csv(annotation_path, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
required <- c("Run", "Fraction", "TechRepMixture", "Mixture", "Channel", "BioReplicate", "Condition")
missing <- setdiff(required, colnames(annotation))
if (length(missing) > 0) {
  stop("MSstatsTMT annotation missing columns: ", paste(missing, collapse = ", "))
}

collect_msstats_files <- function(path) {
  if (dir.exists(path)) {
    files <- list.files(path, pattern = "^msstats\\.csv$", recursive = TRUE, full.names = TRUE)
    files <- files[!grepl("msstats_input|msstats_comparison|msstats_contrasts", files, ignore.case = TRUE)]
    return(files)
  }
  if (file.exists(path)) return(path)
  character(0)
}

read_msstats_table <- function(files) {
  tables <- lapply(
    files,
    function(f) read.csv(f, check.names = FALSE, stringsAsFactors = FALSE)
  )
  if (length(tables) == 1) return(tables[[1]])
  do.call(rbind, tables)
}

msstats_files <- collect_msstats_files(msstats_path)
if (length(msstats_files) == 0) {
  stop("No msstats.csv found at ", msstats_path)
}

msstats_data <- read_msstats_table(msstats_files)

input_tmt <- PhilosophertoMSstatsTMTFormat(
  input = msstats_data,
  annotation = annotation
)
write.csv(input_tmt, "msstats_input.csv", row.names = FALSE)

setwd(rootDir)

has_norm <- any(annotation$Condition == "Norm", na.rm = TRUE)
use_reference_norm <- has_norm && length(unique(annotation$Run)) > 1

quant <- proteinSummarization(
  input_tmt,
  method = "msstats",
  global_norm = TRUE,
  reference_norm = use_reference_norm,
  remove_norm_channel = TRUE,
  remove_empty_channel = TRUE
)

all_conditions <- sort(unique(annotation$Condition))
all_conditions <- all_conditions[
  !is.na(all_conditions) & all_conditions != "" & !all_conditions %in% c("Empty", "Norm")
]
conditions <- all_conditions

if (length(conditions) > 1) {
  contrast.names <- combn(conditions, 2)
  n_comp <- ncol(contrast.names)
  comparison <- matrix(0, nrow = n_comp, ncol = length(all_conditions))
  colnames(comparison) <- all_conditions
  rownames(comparison) <- paste(contrast.names[2, ], contrast.names[1, ], sep = "_v_")
  for (i in seq_len(n_comp)) {
    comparison[i, contrast.names[2, i]] <- 1
    comparison[i, contrast.names[1, i]] <- -1
  }

  comparisonResults <- groupComparisonTMT(contrast.matrix = comparison, data = quant)
  comparison_df <- comparisonResults$ComparisonResult

  format_label <- function(num, den) paste0("(", num, ")v(", den, ")")
  for (i in seq_len(n_comp)) {
    c1 <- contrast.names[1, i]
    c2 <- contrast.names[2, i]
    lbl_new <- format_label(c2, c1)
    comparison_df$Label[comparison_df$Label == rownames(comparison)[i]] <- lbl_new
  }

  write.csv(comparison_df, paste0("msstats_comparison", assay_suffix, ".csv"), row.names = FALSE)

  contrasts_df <- data.frame(row.names = c("1", "2"))
  for (i in seq_len(n_comp)) {
    c1 <- contrast.names[1, i]
    c2 <- contrast.names[2, i]
    col_name <- format_label(c2, c1)
    contrasts_df[[col_name]] <- c(c2, c1)
  }
  write.csv(contrasts_df, paste0("msstats_contrasts", assay_suffix, ".csv"), row.names = TRUE)
}
