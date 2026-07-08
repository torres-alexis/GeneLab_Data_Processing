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

normalize_run_id <- function(x, assay_suffix = "") {
  x <- trimws(as.character(x))
  x <- gsub("\\.mzML$", "", x, ignore.case = TRUE)
  x <- gsub("\\.raw$", "", x, ignore.case = TRUE)
  if (nzchar(assay_suffix)) {
    x <- gsub(assay_suffix, "", x, fixed = TRUE)
  }
  x
}

msstats_run_ids <- function(data) {
  if ("Spectrum.File" %in% names(data)) {
    ids <- data[["Spectrum.File"]]
  } else if ("Run" %in% names(data)) {
    ids <- data[["Run"]]
  } else {
    stop("msstats.csv missing Spectrum.File and Run columns; cannot match MSstatsTMT annotation")
  }
  sort(unique(normalize_run_id(ids)))
}

filter_annotation_to_msstats <- function(annotation_msstats, msstats_runs, assay_suffix = "") {
  annotated_runs <- sort(unique(normalize_run_id(annotation_msstats$Run, assay_suffix)))
  retained_runs <- intersect(annotated_runs, msstats_runs)
  dropped_runs <- setdiff(annotated_runs, msstats_runs)
  extra_msstats_runs <- setdiff(msstats_runs, annotated_runs)
  keep <- normalize_run_id(annotation_msstats$Run, assay_suffix) %in% msstats_runs
  filtered <- annotation_msstats[keep, , drop = FALSE]
  list(
    annotation = filtered,
    annotated_runs = annotated_runs,
    retained_runs = sort(unique(normalize_run_id(filtered$Run, assay_suffix))),
    dropped_runs = dropped_runs,
    extra_msstats_runs = extra_msstats_runs
  )
}

format_msstatstmt_notice_section <- function(label, items) {
  c(
    paste0("  ", label),
    if (length(items)) paste0("    - ", items) else "    - (none)"
  )
}

write_msstatstmt_notice <- function(intro_lines, sections, notice_path) {
  cmd <- paste(commandArgs(trailingOnly = TRUE), collapse = " ")
  lines <- c(
    intro_lines,
    "",
    paste0("MSstatsTMT analysis executed as:\n    msstatstmt_analysis.R ", cmd),
    ""
  )
  for (i in seq_along(sections)) {
    lines <- c(lines, format_msstatstmt_notice_section(sections[[i]]$label, sections[[i]]$items), "")
  }
  writeLines(lines, notice_path)
}

write_dropped_runs_notice <- function(annotated, retained, dropped, extra_msstats, notice_path) {
  intro <- c(
    "MSstatsTMT annotation included runs that were absent from msstats.csv.",
    "Those annotation rows were removed before PhilosophertoMSstatsTMTFormat."
  )
  sections <- list(
    list(label = paste0("Annotation runs (n=", length(annotated), "):"), items = annotated),
    list(label = paste0("Retained for MSstatsTMT (n=", length(retained), "):"), items = retained)
  )
  if (length(dropped)) {
    sections <- c(sections, list(
      list(label = paste0("Dropped annotation runs (n=", length(dropped), "):"), items = dropped)
    ))
  }
  if (length(extra_msstats)) {
    sections <- c(sections, list(
      list(label = paste0("msstats.csv runs not in annotation (n=", length(extra_msstats), "):"), items = extra_msstats)
    ))
  }
  write_msstatstmt_notice(intro, sections, notice_path)
}

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

msstats_runs <- msstats_run_ids(msstats_data)
run_filter <- filter_annotation_to_msstats(annotation_msstats, msstats_runs, assay_suffix)
annotation_msstats <- run_filter$annotation

notice_suffix <- if (nzchar(assay_suffix)) assay_suffix else "_GLProteomics"
dropped_runs_path <- paste0("dropped-runs-msstatstmt", notice_suffix, ".txt")
if (length(run_filter$dropped_runs) || length(run_filter$extra_msstats_runs)) {
  write_dropped_runs_notice(
    run_filter$annotated_runs,
    run_filter$retained_runs,
    run_filter$dropped_runs,
    run_filter$extra_msstats_runs,
    dropped_runs_path
  )
  message(
    "MSstatsTMT: annotation/msstats run mismatch; removed ",
    length(run_filter$dropped_runs),
    " orphan run(s). See ",
    dropped_runs_path
  )
}

if (nrow(annotation_msstats) == 0) {
  stop(
    "MSstatsTMT annotation has no runs present in msstats.csv after filtering. ",
    "See ", dropped_runs_path
  )
}

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
  intro <- c(
    "Some Condition labels from MSstatsTMT_annotation were removed during",
    "proteinSummarization (empty channel / insufficient data after filtering).",
    "MSstatsTMT group comparisons use only retained conditions below."
  )
  sections <- list(
    list(label = paste0("Annotation conditions (n=", length(annotated), "):"), items = annotated),
    list(label = paste0("Retained for groupComparisonTMT (n=", length(retained), "):"), items = retained)
  )
  if (length(dropped)) {
    sections <- c(sections, list(
      list(label = paste0("Dropped conditions (n=", length(dropped), "):"), items = dropped)
    ))
  }
  write_msstatstmt_notice(intro, sections, notice_path)
}

annotated_conditions <- annotation_conditions(annotation_msstats)
conditions <- summarized_conditions(quant)
dropped_conditions <- setdiff(annotated_conditions, conditions)

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
