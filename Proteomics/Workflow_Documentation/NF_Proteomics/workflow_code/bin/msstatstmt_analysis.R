#!/usr/bin/env Rscript
# MSstatsTMT for FragPipe TMT (Philosopher msstats.csv + pipeline annotation).

rm(list = ls())
library(MSstatsTMT)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: msstatstmt_analysis.R <rootDir> <MSstatsTMT_annotation.csv> <msstats_dir_or_file> [assay_suffix]")
}
rootDir <- args[1]
msstats_tmt_annotation_path <- args[2]
msstats_csv_path <- args[3]
assay_suffix <- if (length(args) > 3) args[4] else ""

if (!grepl("/$", rootDir)) rootDir <- paste0(rootDir, "/")

annotation <- read.csv(msstats_tmt_annotation_path, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
required <- c("Run", "Fraction", "TechRepMixture", "Mixture", "Channel", "BioReplicate", "Condition")
missing <- setdiff(required, colnames(annotation))
if (length(missing) > 0) {
  stop("MSstatsTMT annotation missing columns: ", paste(missing, collapse = ", "))
}

anno <- annotation
colnames(anno) <- tolower(colnames(anno))
if (!"condition" %in% colnames(anno)) {
  stop("MSstatsTMT annotation must have condition column")
}
safe_to_label <- if ("condition_name" %in% colnames(anno) && all(nzchar(trimws(anno$condition_name)))) {
  u <- unique(anno[, c("condition", "condition_name")])
  setNames(u$condition_name, u$condition)
} else {
  setNames(anno$condition, anno$condition)
}

print("Condition mapping from MSstatsTMT annotation:")
print(anno[, intersect(c("condition", "condition_name"), colnames(anno)), drop = FALSE])

annotation_msstats <- annotation[, required, drop = FALSE]

collect_msstats_files <- function(path) {
  if (dir.exists(path)) {
    files <- list.files(path, pattern = "^msstats\\.csv$", recursive = TRUE, full.names = TRUE)
    files <- files[!grepl("msstats_input|msstatstmt_comparison|msstatstmt_contrasts|msstats_comparison|msstats_contrasts", files, ignore.case = TRUE)]
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

msstats_files <- collect_msstats_files(msstats_csv_path)
if (length(msstats_files) == 0) {
  stop("No msstats.csv found at ", msstats_csv_path)
}

msstats_data <- read_msstats_table(msstats_files)

input_tmt <- PhilosophertoMSstatsTMTFormat(
  input = msstats_data,
  annotation = annotation_msstats
)

setwd(rootDir)

has_norm <- any(annotation_msstats$Condition == "Norm", na.rm = TRUE)
use_reference_norm <- has_norm && length(unique(annotation_msstats$Run)) > 1

quant <- proteinSummarization(
  input_tmt,
  method = "msstats",
  global_norm = TRUE,
  reference_norm = use_reference_norm,
  remove_norm_channel = TRUE,
  remove_empty_channel = TRUE
)

annotation_conditions <- function(annot) {
  conds <- sort(unique(annot$Condition))
  conds[!is.na(conds) & conds != "" & !conds %in% c("Empty", "Norm")]
}

summarized_conditions <- function(quant_obj) {
  pld <- quant_obj$ProteinLevelData
  if (is.null(pld)) {
    stop("MSstatsTMT proteinSummarization output missing ProteinLevelData")
  }
  col <- if ("Condition" %in% names(pld)) {
    "Condition"
  } else if ("Group" %in% names(pld)) {
    "Group"
  } else {
    stop("MSstatsTMT ProteinLevelData missing Condition/Group column")
  }
  conds <- sort(unique(pld[[col]]))
  conds[!is.na(conds) & conds != "" & !conds %in% c("Empty", "Norm")]
}

write_dropped_conditions_notice <- function(annotated, retained, dropped, notice_path) {
  lines <- c(
    "GeneLab MSstatsTMT dropped conditions notice",
    paste("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
    "",
    "Some Condition labels from MSstatsTMT_annotation were removed during",
    "proteinSummarization (empty channel / insufficient data after filtering).",
    "MSstatsTMT group comparisons use only retained conditions below.",
    "",
    paste0("Annotation conditions (n=", length(annotated), "):"),
    if (length(annotated)) paste0("  - ", annotated) else "  (none)",
    "",
    paste0("Retained for groupComparisonTMT (n=", length(retained), "):"),
    if (length(retained)) paste0("  - ", retained) else "  (none)"
  )
  if (length(dropped)) {
    lines <- c(
      lines,
      "",
      paste0("Dropped conditions (n=", length(dropped), "):"),
      paste0("  - ", dropped)
    )
  }
  writeLines(lines, notice_path)
}

annotated_conditions <- annotation_conditions(annotation_msstats)
conditions <- summarized_conditions(quant)
dropped_conditions <- setdiff(annotated_conditions, conditions)

notice_suffix <- if (nzchar(assay_suffix)) assay_suffix else "_GLProteomics"
notice_path <- paste0("dropped-conditions-msstatstmt", notice_suffix, ".txt")
if (length(dropped_conditions)) {
  write_dropped_conditions_notice(annotated_conditions, conditions, dropped_conditions, notice_path)
  message("MSstatsTMT: some annotation conditions were dropped after summarization, see ", notice_path)
}

if (length(conditions) > 1) {
  contrast.names <- combn(conditions, 2)
  n_comparisons <- ncol(contrast.names)
  comparison <- matrix(0, nrow = n_comparisons, ncol = length(conditions))
  colnames(comparison) <- conditions
  rownames(comparison) <- paste(contrast.names[2, ], contrast.names[1, ], sep = "_v_")
  for (i in 1:n_comparisons) {
    comparison[i, contrast.names[2, i]] <- 1
    comparison[i, contrast.names[1, i]] <- -1
  }

  comparisonResults <- groupComparisonTMT(contrast.matrix = comparison, data = quant)
  comparison_df <- comparisonResults$ComparisonResult

  format_label <- function(cond1, cond2) paste0("(", cond1, ")v(", cond2, ")")
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

  write.csv(comparison_df, paste0("msstatstmt_comparison", assay_suffix, ".csv"), row.names = FALSE)

  contrasts_df <- data.frame(row.names = c("1", "2"))
  for (i in 1:n_comparisons) {
    c1 <- contrast.names[1, i]
    c2 <- contrast.names[2, i]
    r1 <- if (is.na(safe_to_label[c1])) c1 else safe_to_label[c1]
    r2 <- if (is.na(safe_to_label[c2])) c2 else safe_to_label[c2]
    col_name <- format_label(r2, r1)
    contrasts_df[[col_name]] <- c(c2, c1)
  }
  write.csv(contrasts_df, paste0("msstatstmt_contrasts", assay_suffix, ".csv"), row.names = TRUE)
}
