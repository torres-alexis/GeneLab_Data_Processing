#!/usr/bin/env Rscript
# MSstatsPTM for FragPipe TMT16-phospho (Philosopher msstats.csv + MSstatsTMT annotation).

rm(list = ls())
library(data.table)
library(MSstatsPTM)

.script_dir <- local({
  args <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", args[grepl("^--file=", args)])
  if (length(f)) dirname(normalizePath(f[[1]])) else getwd()
})
source(file.path(.script_dir, "decoy_contam.R"))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: msstatstmtptm_analysis.R <rootDir> <MSstatsTMT_annotation.csv> <msstats_dir_or_file> [assay_suffix] [mod_id_col] [msstats_protein.csv] [msstats_protein_annotation.csv]")
}
rootDir <- args[1]
msstats_tmt_annotation_path <- args[2]
msstats_csv_path <- args[3]
assay_suffix <- if (length(args) > 3) args[4] else ""
mod_id_col <- if (length(args) > 4 && nzchar(args[5])) args[5] else ""
msstats_protein_path <- if (length(args) > 5 && nzchar(args[6])) args[6] else ""
msstats_protein_annot_path <- if (length(args) > 6 && nzchar(args[7])) args[7] else ""

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

format_msstatstmtptm_notice_section <- function(label, items) {
  c(
    paste0("  ", label),
    if (length(items)) paste0("    - ", items) else "    - (none)"
  )
}

write_msstatstmtptm_notice <- function(intro_lines, sections, notice_path) {
  cmd <- paste(commandArgs(trailingOnly = TRUE), collapse = " ")
  lines <- c(
    intro_lines,
    "",
    paste0("MSstatsPTM analysis executed as:\n    msstatstmtptm_analysis.R ", cmd),
    ""
  )
  for (i in seq_along(sections)) {
    lines <- c(lines, format_msstatstmtptm_notice_section(sections[[i]]$label, sections[[i]]$items), "")
  }
  writeLines(lines, notice_path)
}

write_dropped_runs_notice <- function(annotated, retained, dropped, extra_msstats, notice_path) {
  intro <- c(
    "MSstatsTMT annotation included runs that were absent from msstats.csv.",
    "Those annotation rows were removed before FragPipetoMSstatsPTMFormat."
  )
  sections <- list(
    list(label = paste0("Annotation runs (n=", length(annotated), "):"), items = annotated),
    list(label = paste0("Retained for MSstatsPTM (n=", length(retained), "):"), items = retained)
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
  write_msstatstmtptm_notice(intro, sections, notice_path)
}

collect_msstats_files <- function(path) {
  if (dir.exists(path)) {
    files <- list.files(path, pattern = "^msstats\\.csv$", recursive = TRUE, full.names = TRUE)
    files <- files[!grepl("msstats_input|msstatstmtptm|msstatstmt_comparison|msstatstmt_contrasts|msstats_comparison|msstats_contrasts|msstatsptm", files, ignore.case = TRUE)]
    return(files)
  }
  if (file.exists(path)) return(path)
  character(0)
}

read_msstats_table <- function(files) {
  tables <- lapply(
    files,
    function(f) as.data.table(read.csv(f, check.names = FALSE, stringsAsFactors = FALSE))
  )
  if (length(tables) == 1) return(tables[[1]])
  rbindlist(tables, use.names = TRUE, fill = TRUE)
}

# Package grepl is unanchored, so "M" also hits Modified/Mapped. Pass the full column.
mod_column <- function(nms, token) {
  hits <- nms[grepl(paste0("^", token, "[:.][0-9]"), nms)]
  if (length(hits) == 1) return(hits)
  stop(
    "mod_id_col ", token, " matched ", length(hits), " localization column(s): ",
    paste(if (length(hits)) hits else nms, collapse = ", ")
  )
}

msstats_files <- collect_msstats_files(msstats_csv_path)
if (length(msstats_files) == 0) {
  stop("No msstats.csv found at ", msstats_csv_path)
}

msstats_data_raw <- read_msstats_table(msstats_files)

nms <- names(msstats_data_raw)
has_channel_cols <- "Channel" %in% nms || any(grepl("^Channel[ .]", nms))
has_is_unique <- "Is.Unique" %in% nms
if (!has_channel_cols || !has_is_unique) {
  notice_suffix <- if (nzchar(assay_suffix)) assay_suffix else "_GLProteomics"
  notice_path <- paste0("dropped-runs-msstatstmtptm", notice_suffix, ".txt")
  writeLines(
    c(
      "msstats.csv is not Philosopher TMT format (need Is.Unique and Channel / 'Channel *' columns).",
      "MSstatsPTM skipped. Use --tmt_extraction_tool Philosopher so FragPipe writes philosopher-msstats.",
      paste("columns:", paste(nms, collapse = ", "))
    ),
    notice_path
  )
  quit(save = "no", status = 0)
}

if (!nzchar(mod_id_col)) {
  mods <- unique(sub("^([A-Za-z]+)[:.].*$", "\\1", nms[grepl("^[A-Za-z]+[:.][0-9]", nms)]))
  if (length(mods) != 1) {
    stop("Pass mod_id_col. Localization columns: ", paste(mods, collapse = ", "))
  }
  mod_id_col <- mods
}
mod_col <- mod_column(nms, mod_id_col)
message("Using mod_id_col=", mod_id_col, " column=", mod_col)
mod_file_tag <- paste0("_", mod_id_col)
rm_oxidation_m <- !identical(mod_id_col, "M")

msstats_runs <- msstats_run_ids(msstats_data_raw)
run_filter <- filter_annotation_to_msstats(annotation_msstats, msstats_runs, assay_suffix)
annotation_msstats <- run_filter$annotation

notice_suffix <- paste0(mod_file_tag, if (nzchar(assay_suffix)) assay_suffix else "_GLProteomics")
dropped_runs_path <- paste0("dropped-runs-msstatstmtptm", notice_suffix, ".txt")
if (length(run_filter$dropped_runs) || length(run_filter$extra_msstats_runs)) {
  write_dropped_runs_notice(
    run_filter$annotated_runs,
    run_filter$retained_runs,
    run_filter$dropped_runs,
    run_filter$extra_msstats_runs,
    dropped_runs_path
  )
  message(
    "MSstatsPTM: annotation/msstats run mismatch; removed ",
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

raw_protein <- NULL
annotation_protein <- NULL
if (nzchar(msstats_protein_path)) {
  if (!nzchar(msstats_protein_annot_path)) {
    stop("msstats_protein.csv requires msstats_protein_annotation")
  }
  if (!file.exists(msstats_protein_path)) {
    stop("msstats_protein.csv not found: ", msstats_protein_path)
  }
  if (!file.exists(msstats_protein_annot_path)) {
    stop("msstats_protein_annotation.csv not found: ", msstats_protein_annot_path)
  }
  raw_protein <- read_msstats_table(msstats_protein_path)
  protein_annotation <- read.csv(msstats_protein_annot_path, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
  missing_prot <- setdiff(required, colnames(protein_annotation))
  if (length(missing_prot) > 0) {
    stop("Protein MSstatsTMT annotation missing columns: ", paste(missing_prot, collapse = ", "))
  }
  protein_annotation <- protein_annotation[, required, drop = FALSE]
  protein_runs <- msstats_run_ids(raw_protein)
  protein_filter <- filter_annotation_to_msstats(protein_annotation, protein_runs, assay_suffix)
  annotation_protein <- protein_filter$annotation
  protein_notice <- paste0("dropped-runs-msstatstmtptm", notice_suffix, "-protein.txt")
  if (length(protein_filter$dropped_runs) || length(protein_filter$extra_msstats_runs)) {
    write_dropped_runs_notice(
      protein_filter$annotated_runs,
      protein_filter$retained_runs,
      protein_filter$dropped_runs,
      protein_filter$extra_msstats_runs,
      protein_notice
    )
    message(
      "MSstatsPTM: annotation/msstats run mismatch; removed ",
      length(protein_filter$dropped_runs),
      " orphan run(s). See ",
      protein_notice
    )
  }
  if (nrow(annotation_protein) == 0) {
    stop(
      "MSstatsTMT annotation has no runs present in msstats.csv after filtering. ",
      "See ", protein_notice
    )
  }
}

msstats_data <- FragPipetoMSstatsPTMFormat(
  input = msstats_data_raw,
  annotation = annotation_msstats,
  input_protein = raw_protein,
  annotation_protein = annotation_protein,
  label_type = "TMT",
  mod_id_col = mod_col,
  localization_cutoff = 0.75,
  rmPeptide_OxidationM = rm_oxidation_m,
  remove_unlocalized_peptides = TRUE,
  protein_id_col = "Protein",
  peptide_id_col = "Peptide.Sequence",
  use_log_file = FALSE,
  append = FALSE,
  verbose = TRUE
)

if (is.null(msstats_data$PTM) || nrow(msstats_data$PTM) == 0) {
  stop("FragPipetoMSstatsPTMFormat returned empty PTM table (check mod_id_col / localization)")
}

setwd(rootDir)

has_norm <- any(annotation_msstats$Condition == "Norm", na.rm = TRUE)
use_reference_norm <- has_norm && length(unique(annotation_msstats$Run)) > 1

summarized <- dataSummarizationPTM_TMT(
  msstats_data,
  method = "msstats",
  global_norm = TRUE,
  global_norm.PTM = TRUE,
  reference_norm = use_reference_norm,
  reference_norm.PTM = use_reference_norm,
  remove_norm_channel = TRUE,
  remove_empty_channel = TRUE,
  use_log_file = FALSE,
  append = FALSE,
  verbose = TRUE
)

annotation_conditions <- function(annot) {
  conds <- sort(unique(annot$Condition))
  conds[!is.na(conds) & conds != "" & !conds %in% c("Empty", "Norm")]
}

summarized_conditions <- function(quant_obj) {
  pld <- quant_obj$PTM$ProteinLevelData
  if (is.null(pld)) {
    stop("MSstatsPTM dataSummarizationPTM_TMT output missing PTM$ProteinLevelData")
  }
  col <- if ("Condition" %in% names(pld)) {
    "Condition"
  } else if ("Group" %in% names(pld)) {
    "Group"
  } else {
    stop("PTM ProteinLevelData missing Condition/Group column")
  }
  conds <- sort(unique(pld[[col]]))
  conds[!is.na(conds) & conds != "" & !conds %in% c("Empty", "Norm")]
}

write_dropped_conditions_notice <- function(annotated, retained, dropped, notice_path) {
  intro <- c(
    "Some Condition labels from MSstatsTMT_annotation were removed during",
    "dataSummarizationPTM_TMT (empty channel / insufficient data after filtering).",
    "MSstatsPTM group comparisons use only retained conditions below."
  )
  sections <- list(
    list(label = paste0("Annotation conditions (n=", length(annotated), "):"), items = annotated),
    list(label = paste0("Retained for groupComparisonPTM (n=", length(retained), "):"), items = retained)
  )
  if (length(dropped)) {
    sections <- c(sections, list(
      list(label = paste0("Dropped conditions (n=", length(dropped), "):"), items = dropped)
    ))
  }
  write_msstatstmtptm_notice(intro, sections, notice_path)
}

annotated_conditions <- annotation_conditions(annotation_msstats)
conditions <- summarized_conditions(summarized)
dropped_conditions <- setdiff(annotated_conditions, conditions)

notice_path <- paste0("dropped-conditions-msstatstmtptm", notice_suffix, ".txt")
if (length(dropped_conditions)) {
  write_dropped_conditions_notice(annotated_conditions, conditions, dropped_conditions, notice_path)
  message("MSstatsPTM: some annotation conditions were dropped after summarization, see ", notice_path)
}

format_label <- function(cond1, cond2) paste0("(", cond1, ")v(", cond2, ")")

relabel_comparison <- function(comparison_df, contrast.names, comparison, safe_to_label) {
  n_comparisons <- ncol(contrast.names)
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
  scrub_msstats_comparison(comparison_df)
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

  model <- groupComparisonPTM(
    summarized,
    contrast.matrix = comparison,
    ptm_label_type = "TMT",
    protein_label_type = "TMT",
    use_log_file = FALSE,
    append = FALSE,
    verbose = TRUE
  )

  primary <- if (!is.null(model$ADJUSTED.Model) && nrow(model$ADJUSTED.Model) > 0) {
    model$ADJUSTED.Model
  } else {
    model$PTM.Model
  }
  primary <- relabel_comparison(primary, contrast.names, comparison, safe_to_label)
  write.csv(primary, paste0("msstatstmtptm_comparison", mod_file_tag, assay_suffix, ".csv"), row.names = FALSE)

  if (!is.null(model$PTM.Model)) {
    ptm_df <- relabel_comparison(model$PTM.Model, contrast.names, comparison, safe_to_label)
    write.csv(ptm_df, paste0("msstatstmtptm_ptm_comparison", mod_file_tag, assay_suffix, ".csv"), row.names = FALSE)
  }
  if (!is.null(model$PROTEIN.Model) && nrow(model$PROTEIN.Model) > 0) {
    prot_df <- relabel_comparison(model$PROTEIN.Model, contrast.names, comparison, safe_to_label)
    write.csv(prot_df, paste0("msstatstmtptm_protein_comparison", mod_file_tag, assay_suffix, ".csv"), row.names = FALSE)
  }
  if (!is.null(model$ADJUSTED.Model) && nrow(model$ADJUSTED.Model) > 0) {
    adj_df <- relabel_comparison(model$ADJUSTED.Model, contrast.names, comparison, safe_to_label)
    write.csv(adj_df, paste0("msstatstmtptm_adjusted_comparison", mod_file_tag, assay_suffix, ".csv"), row.names = FALSE)
  }

  contrasts_df <- data.frame(row.names = c("1", "2"))
  for (i in 1:n_comparisons) {
    c1 <- contrast.names[1, i]
    c2 <- contrast.names[2, i]
    r1 <- if (is.na(safe_to_label[c1])) c1 else safe_to_label[c1]
    r2 <- if (is.na(safe_to_label[c2])) c2 else safe_to_label[c2]
    contrasts_df[[format_label(r2, r1)]] <- c(c2, c1)
  }
  write.csv(contrasts_df, paste0("msstatstmtptm_contrasts", mod_file_tag, assay_suffix, ".csv"), row.names = TRUE)
}
