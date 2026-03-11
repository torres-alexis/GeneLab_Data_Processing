#!/usr/bin/Rscript
# FragPipe-Analyst downstream via FragPipeAnalystR package.
# Reference: Reference_Repos/FragPipeAnalystR
# CLI: optparse, types (LFQ/TMT/DIA), levels (protein/peptide/gene/site).

library(optparse)
library(ggplot2)  # ensure geom_violin available for FragPipeAnalystR::plot_feature

option_list <- list(
  # --- make_se_from_files params (FragPipeAnalystR) ---
  # mode: LFQ | TMT | DIA
  # level: protein | peptide | gene | site | glycan
  # lfq_type: Intensity | MaxLFQ | Spectral Count (LFQ only)
  make_option(c("--experiment_annotation"), type = "character", default = NULL,
    help = "Path to experiment annotation TSV", metavar = "FILE"),
  make_option(c("--quantification_file"), type = "character", default = NULL,
    help = "Path to quantification file", metavar = "FILE"),
  make_option(c("--mode"), type = "character", default = NULL,
    help = "Quantification mode: LFQ, TMT, or DIA", metavar = "MODE"),
  make_option(c("--level"), type = "character", default = NULL,
    help = "Analysis level: protein, peptide, gene, site, or glycan", metavar = "LEVEL"),
  make_option(c("--output_dir"), type = "character", default = "output/",
    help = "Output directory", metavar = "DIR"),
  make_option(c("--lfq_type"), type = "character", default = "Intensity",
    help = "LFQ column type: Intensity, MaxLFQ, or Spectral Count (LFQ mode only)", metavar = "STRING"),
  # --- filter ---
  make_option(c("--min_global_appearance"), type = "numeric", default = 0,
    help = "Min %% present across all samples (0-100). 0 = unfiltered", metavar = "NUMERIC"),
  make_option(c("--min_appearance_one_condition"), type = "numeric", default = 0,
    help = "Min %% present in at least one condition (0-100). 0 = unfiltered", metavar = "NUMERIC"),
  # --- norm ---
  make_option(c("--normalization_method"), type = "character", default = "none",
    help = "Normalization: none, vsn, MD, or GN (FragPipeAnalystR)", metavar = "STRING"),
  # --- impute (FragPipeAnalystR: manual_impute / impute) ---
  make_option(c("--imputation_type"), type = "character", default = "Perseus-type",
    help = "Imputation: none, Perseus-type (fun=man), knn, MLE, min, zero, bpca, QRILC, MinDet, MinProb, nbavg, mixed", metavar = "STRING"),
  make_option(c("--imputation_shift"), type = "numeric", default = 1.8,
    help = "Perseus-type: manual_impute shift (SD units)", metavar = "NUMERIC"),
  make_option(c("--imputation_scale"), type = "numeric", default = 0.3,
    help = "Perseus-type: manual_impute scale factor", metavar = "NUMERIC"),
  # --- DE ---
  make_option(c("--de_alpha"), type = "numeric", default = 0.05,
    help = "Adjusted p-value threshold for DE significance", metavar = "NUMERIC"),
  make_option(c("--de_lfc"), type = "numeric", default = 1.0,
    help = "Log2 fold change threshold for DE significance", metavar = "NUMERIC"),
  make_option(c("--de_fdr"), type = "character", default = "Benjamini Hochberg",
    help = "FDR correction: Benjamini Hochberg or Local and tail area-based", metavar = "STRING"),
  # --- feature plots ---
  make_option(c("--feature_list_protein"), type = "character", default = "",
    help = "Comma-separated protein IDs for feature plots (protein level). Empty = use top_n_protein", metavar = "STRING"),
  make_option(c("--feature_list_gene"), type = "character", default = "",
    help = "Comma-separated gene names for feature plots. Empty = use top_n_gene", metavar = "STRING"),
  make_option(c("--feature_list_peptide"), type = "character", default = "",
    help = "Comma-separated peptide IDs for feature plots (peptide level). Empty = use top_n_peptide", metavar = "STRING"),
  make_option(c("--feature_list_site"), type = "character", default = "",
    help = "Comma-separated site IDs for feature plots (site level). Empty = use top_n_site", metavar = "STRING"),
  make_option(c("--top_n_protein"), type = "integer", default = 10,
    help = "Top N variable by protein ID when list empty", metavar = "INTEGER"),
  make_option(c("--top_n_gene"), type = "integer", default = 10,
    help = "Top N variable by gene when list empty", metavar = "INTEGER"),
  make_option(c("--top_n_peptide"), type = "integer", default = 10,
    help = "Top N variable by peptide ID when list empty (peptide level)", metavar = "INTEGER"),
  make_option(c("--top_n_site"), type = "integer", default = 10,
    help = "Top N variable by site ID when list empty (site level)", metavar = "INTEGER"),
  # --- enrichment (FragPipeAnalystR or_test, Enrichr backend) ---
  make_option(c("--enrichment_database"), type = "character", default = "Hallmark,GO_Biological_Process_2021",
    help = "Comma-separated Enrichr DB(s): GO_Biological_Process_2021, GO_Cellular_Component_2021, GO_Molecular_Function_2021, MSigDB_Hallmark_2020, KEGG_2021_Human, Reactome_2022. Aliases: Hallmark, KEGG, Reactome. Empty = skip", metavar = "STRING"),
  make_option(c("--enrichment_direction"), type = "character", default = "Up,Down",
    help = "Enrichment direction(s): Up, Down, or comma-separated (e.g. Up,Down)", metavar = "STRING"),
  make_option(c("--gsea_database"), type = "character", default = "Hallmark,GO_Biological_Process_2021",
    help = "GSEA DB(s): Hallmark, GO_Biological_Process_2021, GO_Cellular_Component_2021, GO_Molecular_Function_2021, KEGG_2021_Human. Protein/gene/site only. Empty = skip", metavar = "STRING"),
  # --- legacy (not yet wired) ---
  make_option(c("--qc_plot_data"), type = "character", default = "nonimputed",
    help = "Data for PCA, correlation, feature, CVs: imputed or nonimputed", metavar = "STRING"),
  make_option(c("--sample_cvs_full_range"), type = "character", default = "false",
    help = "Sample CVs: true = full range, false = 0-1", metavar = "STRING"),
  make_option(c("--volcano_display_names"), type = "character", default = "true",
    help = "Volcano: add_names (label significant points): true/false", metavar = "STRING"),
  make_option(c("--volcano_show_gene"), type = "character", default = "true",
    help = "Volcano: name_col=Gene when true; else use row ID", metavar = "STRING"),
  make_option(c("--gene_annotations"), type = "character", default = "",
    help = "Path or URL to gene annotations TSV/CSV. Merges into DE_results on Gene using the annotation column with most matches.", metavar = "STRING"),
  make_option(c("--assay_suffix"), type = "character", default = "",
    help = "Assay suffix for volcano filenames, e.g. GLProteomics. Empty = no suffix.", metavar = "STRING")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

.or <- function(x, y) if (is.null(x) || is.na(x) || (is.character(x) && !nzchar(trimws(x)))) y else x
.oneof <- function(x, valid, param) {
  v <- tolower(trimws(as.character(x)))
  if (!v %in% tolower(valid)) stop("Invalid ", param, ": '", x, "'. Valid: ", paste(valid, collapse = ", "))
  valid[match(v, tolower(valid))]
}

# Validate required
if (is.null(opt$experiment_annotation) || is.null(opt$quantification_file) || is.null(opt$mode) || is.null(opt$level)) {
  stop("Required: --experiment_annotation, --quantification_file, --mode, --level")
}
mode <- .oneof(opt$mode, c("LFQ", "TMT", "DIA"), "mode")
level <- .oneof(opt$level, c("protein", "peptide", "gene", "site", "glycan"), "level")

output_dir <- normalizePath(opt$output_dir, mustWork = FALSE)
if (!dir.exists(output_dir)) dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
output_dir <- normalizePath(output_dir, mustWork = TRUE)
qc_dir <- file.path(output_dir, "qc")
dir.create(qc_dir, showWarnings = FALSE, recursive = TRUE)
comparison_dir <- file.path(output_dir, "comparison")
dir.create(comparison_dir, showWarnings = FALSE, recursive = TRUE)
de_dir <- file.path(output_dir, "de")
dir.create(de_dir, showWarnings = FALSE, recursive = TRUE)
volcano_dir <- file.path(de_dir, "volcano")
dir.create(volcano_dir, showWarnings = FALSE, recursive = TRUE)
enr_dir <- file.path(output_dir, "enrichment")
dir.create(enr_dir, showWarnings = FALSE, recursive = TRUE)
enr_dir_or <- file.path(enr_dir, "or")
enr_dir_gsea <- file.path(enr_dir, "gsea")
dir.create(enr_dir_or, showWarnings = FALSE, recursive = TRUE)
dir.create(enr_dir_gsea, showWarnings = FALSE, recursive = TRUE)
assay_suffix <- trimws(.or(opt$assay_suffix, ""))
fn_ <- function(base, ext) paste0(base, ".", ext)

# Log raw inputs (reproducibility)
param_path <- file.path(output_dir, fn_("fp_analyst_parameters", "txt"))
writeLines(c(
  "FragPipe-Analyst Parameters (FragPipeAnalystR package)",
  "=====================================================",
  paste("experiment_annotation:", opt$experiment_annotation),
  paste("quantification_file:", opt$quantification_file),
  paste("mode:", mode),
  paste("level:", level),
  paste("output_dir:", output_dir),
  paste("lfq_type:", opt$lfq_type),
  paste("normalization_method:", opt$normalization_method),
  paste("imputation_type:", opt$imputation_type),
  paste("imputation_shift:", opt$imputation_shift),
  paste("imputation_scale:", opt$imputation_scale),
  paste("min_global_appearance:", opt$min_global_appearance),
  paste("min_appearance_one_condition:", opt$min_appearance_one_condition),
  paste("de_alpha:", opt$de_alpha),
  paste("de_lfc:", opt$de_lfc),
  paste("enrichment_database:", opt$enrichment_database),
  paste("enrichment_direction:", opt$enrichment_direction),
  paste("gsea_database:", opt$gsea_database),
  paste("gene_annotations:", opt$gene_annotations)
), param_path)

# Parse typed params (for downstream use)
lfq_type <- .oneof(.or(opt$lfq_type, "Intensity"), c("Intensity", "MaxLFQ", "Spectral Count"), "lfq_type")
de_alpha <- as.numeric(.or(opt$de_alpha, 0.05))
de_lfc <- as.numeric(.or(opt$de_lfc, 1.0))
min_global <- as.numeric(.or(opt$min_global_appearance, 0))
min_cond <- as.numeric(.or(opt$min_appearance_one_condition, 0))
norm_method <- .oneof(.or(opt$normalization_method, "none"), c("none", "vsn", "MD", "GN"), "normalization_method")
imp_type_raw <- trimws(.or(opt$imputation_type, "Perseus-type"))
imp_valid <- c("none", "Perseus-type", "knn", "MLE", "min", "zero", "bpca", "QRILC", "MinDet", "MinProb", "nbavg", "mixed")
imp_type <- .oneof(imp_type_raw, imp_valid, "imputation_type")
imp_shift <- as.numeric(.or(opt$imputation_shift, 1.8))
imp_scale <- as.numeric(.or(opt$imputation_scale, 0.3))
enrichment_dbs <- if (nzchar(trimws(.or(opt$enrichment_database, ""))))
  unique(trimws(strsplit(trimws(opt$enrichment_database), "\\s*,\\s*")[[1]])) else character(0)
# Parse comma-separated enrichment directions; each must be Up or Down (Both -> Up,Down for backward compat)
enrichment_dir_raw <- trimws(.or(opt$enrichment_direction, "Up,Down"))
enrichment_dir_tokens <- unique(trimws(strsplit(enrichment_dir_raw, "\\s*,\\s*")[[1]]))
enrichment_dir_tokens <- enrichment_dir_tokens[nzchar(enrichment_dir_tokens)]
enrichment_dirs <- if (length(enrichment_dir_tokens) == 0) c("Up", "Down") else {
  expand <- function(t) if (tolower(trimws(t)) == "both") c("Up", "Down") else .oneof(t, c("Up", "Down"), "enrichment_direction")
  unique(unlist(lapply(enrichment_dir_tokens, expand)))
}
gsea_dbs <- if (nzchar(trimws(.or(opt$gsea_database, ""))))
  unique(trimws(strsplit(trimws(opt$gsea_database), "\\s*,\\s*")[[1]])) else character(0)
volcano_add_names <- .oneof(.or(opt$volcano_display_names, "true"), c("true", "false"), "volcano_display_names") == "true"
volcano_name_col <- .oneof(.or(opt$volcano_show_gene, "true"), c("true", "false"), "volcano_show_gene") == "true"

cat("Params logged to:", param_path, "\n")
cat("mode:", mode, "| level:", level, "| lfq_type:", lfq_type, "\n")

# --- make_se_from_files (FragPipeAnalystR) ---
# type:     LFQ | TMT | DIA
# level:    protein | peptide | gene | site | glycan  (TMT default: gene; LFQ/DIA default: protein)
# lfq_type: Intensity | MaxLFQ | Spectral Count  (LFQ mode only; maps to quant columns)
# Optional (not in CLI): exp_type (global|phospho|glyco|acetyl|ubiquit), log2transform, gencode, additional_cols
# LFQ: protein matrix has one column per sample; annotation may have multiple rows per sample (tech reps).
# Deduplicate by sample so rownames(expdesign) <- sample_name does not hit duplicate row.names.
exp_anno_path <- opt$experiment_annotation
if (mode == "LFQ") {
  anno_pre <- read.table(opt$experiment_annotation, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  colnames(anno_pre) <- tolower(colnames(anno_pre))
  if (!"sample" %in% colnames(anno_pre)) stop("experiment_annotation requires 'sample' column for LFQ")
  dup <- duplicated(anno_pre$sample)
  if (any(dup)) {
    anno_dedup <- anno_pre[!dup, , drop = FALSE]
    exp_anno_path <- tempfile(fileext = ".tsv")
    write.table(anno_dedup, exp_anno_path, sep = "\t", row.names = FALSE, quote = FALSE)
    cat("Deduplicated annotation: ", sum(dup), " duplicate sample rows removed, ", nrow(anno_dedup), " unique samples\n")
  }
}
library(FragPipeAnalystR)
data_se <- make_se_from_files(
  opt$quantification_file,
  exp_anno_path,
  type = mode,
  level = level,
  lfq_type = lfq_type
)
if (is.null(data_se)) stop("make_se_from_files failed")
cat("SummarizedExperiment:", nrow(data_se), "features,", ncol(data_se), "samples\n")

# --- Filter (our code; FragPipeAnalystR has no filter) ---
# global_filter: keep rows with <= (100 - min_global)% missing globally
# filter_by_condition: keep rows with >= min_cond% valid in at least one condition
global_filter <- function(se, pct_present) {
  pct_na_max <- (100 - pct_present) / 100
  ridx <- rowSums(is.na(assay(se))) / ncol(assay(se)) <= pct_na_max
  se[ridx, ]
}
filter_by_condition <- function(se, min_pct) {
  min_pct <- min_pct / 100
  conds <- unique(colData(se)$condition)
  keep <- rep(FALSE, nrow(se))
  for (c in conds) {
    se_c <- se[, colData(se)$condition == c]
    keep <- keep | (rowSums(!is.na(assay(se_c))) / ncol(se_c) >= min_pct)
  }
  se[keep, ]
}
filtered_se <- data_se
if (min_global > 0) {
  filtered_se <- global_filter(filtered_se, min_global)
  cat("global_filter: kept", nrow(filtered_se), "features (min", min_global, "% present globally)\n")
}
if (min_cond > 0) {
  filtered_se <- filter_by_condition(filtered_se, min_cond)
  cat("filter_by_condition: kept", nrow(filtered_se), "features (min", min_cond, "% in one condition)\n")
}

# --- Normalization (FragPipeAnalystR: MD_normalization, GN_normalization, VSN_normalization) ---
# none | MD (median subtraction) | GN (median + MAD scaling) | vsn (variance-stabilizing; LFQ/DIA intensity only)
if (norm_method != "none") {
  filtered_se <- switch(tolower(norm_method),
    "md" = MD_normalization(filtered_se),
    "gn" = GN_normalization(filtered_se),
    "vsn" = VSN_normalization(filtered_se),
    stop("Invalid normalization_method: ", norm_method)
  )
  cat("Normalization applied:", norm_method, "\n")
}
normalized_se <- filtered_se

# --- Imputation (FragPipeAnalystR: impute with fun=man for Perseus-type, else MSnbase) ---
# Perseus-type -> manual_impute(se, shift, scale); others -> impute(se, fun=...)
imputed_se <- normalized_se
if (imp_type != "none") {
  imp_fun <- if (imp_type == "Perseus-type") "man" else imp_type
  if (imp_fun == "man") {
    imputed_se <- manual_impute(normalized_se, shift = imp_shift, scale = imp_scale, seed = 123L)
  } else {
    imputed_se <- impute(normalized_se, fun = imp_fun, seed = 123L)
  }
  cat("Imputation applied:", imp_type, "\n")
}

# --- QC plots (FragPipeAnalystR) ---
# Data: imputed_se or normalized_se per qc_plot_data (imputed/nonimputed)
qc_use_imputed <- tolower(trimws(.or(opt$qc_plot_data, "nonimputed"))) == "imputed"
qc_se <- if (qc_use_imputed) imputed_se else normalized_se
# PCA (needs complete cases; plot_pca filters to complete.cases internally)
pca_se <- qc_se
n_pca <- sum(complete.cases(assay(pca_se)))
if (n_pca < 2 && !qc_use_imputed && imp_type != "none") {
  warning("Only ", n_pca, " complete features for PCA (nonimputed); using imputed data instead.")
  pca_se <- imputed_se
  n_pca <- sum(complete.cases(assay(pca_se)))
}
if (n_pca >= 2) {
  p_pca <- plot_pca(pca_se, indicate = "condition", n = min(500, n_pca), plot = TRUE)
  ggplot2::ggsave(file.path(qc_dir, fn_("pca", "pdf")), p_pca, width = 8, height = 6)
  ggplot2::ggsave(file.path(qc_dir, fn_("pca", "png")), p_pca, width = 8, height = 6, dpi = 150)
  cat("PCA plot saved\n")
} else if (n_pca < 2) {
  warning("Skipping PCA: only ", n_pca, " features (need >= 2). Use qc_plot_data=imputed for sparse data.")
}
# Correlation heatmap (cor() needs complete cases; plot_correlation_heatmap uses use="complete.obs")
cor_se <- qc_se
n_cor_complete <- sum(complete.cases(assay(cor_se)))
if (n_cor_complete < 2 && !qc_use_imputed && imp_type != "none") {
  cor_se <- imputed_se
  n_cor_complete <- sum(complete.cases(assay(cor_se)))
}
if (n_cor_complete >= 2) {
  ht_cor <- plot_correlation_heatmap(cor_se, indicate = "condition")
  pdf(file.path(comparison_dir, fn_("correlation_heatmap", "pdf")), width = 8, height = 7)
  ComplexHeatmap::draw(ht_cor, heatmap_legend_side = "top")
  dev.off()
  png(file.path(comparison_dir, fn_("correlation_heatmap", "png")), width = 8, height = 7, units = "in", res = 150)
  ComplexHeatmap::draw(ht_cor, heatmap_legend_side = "top")
  dev.off()
  cat("Correlation heatmap saved\n")
} else {
  warning("Skipping correlation heatmap: only ", n_cor_complete, " complete features (need >= 2). Imputation produced NAs.")
}
# Missing value heatmap (uses normalized_se; skip if no NAs)
if (any(is.na(assay(normalized_se)))) {
  pdf(file.path(qc_dir, fn_("missing_value_heatmap", "pdf")), width = 8, height = 6)
  plot_missval_heatmap(normalized_se)
  dev.off()
  png(file.path(qc_dir, fn_("missing_value_heatmap", "png")), width = 8, height = 6, units = "in", res = 150)
  plot_missval_heatmap(normalized_se)
  dev.off()
  cat("Missing value heatmap saved\n")
}
# plotCumulativeMissingPercent: not in FragPipeAnalystR README; commented out
# tryCatch({
#   p_cummiss <- plotCumulativeMissingPercent(normalized_se, title = paste0(mode, " ", level))
#   report_plots[[length(report_plots) + 1]] <- p_cummiss
#   ggplot2::ggsave(file.path(qc_dir, fn_("cumulative_missing_percent", "pdf")), p_cummiss, width = 8, height = 5)
#   ggplot2::ggsave(file.path(qc_dir, fn_("cumulative_missing_percent", "png")), p_cummiss, width = 8, height = 5, dpi = 150)
#   cat("Cumulative missing % plot saved\n")
# }, error = function(e) warning("Cumulative missing % plot failed: ", conditionMessage(e)))
# Feature numbers (barplot: features per sample)
p_fn <- plot_feature_numbers(normalized_se, fill = "condition")
ggplot2::ggsave(file.path(qc_dir, fn_("feature_numbers", "pdf")), p_fn, width = 8, height = 5)
ggplot2::ggsave(file.path(qc_dir, fn_("feature_numbers", "png")), p_fn, width = 8, height = 5, dpi = 150)
cat("Feature numbers plot saved\n")
# Sample CVs (can fail when sparse data yields all NA/Inf CVs)
cvs_scale <- !(tolower(trimws(.or(opt$sample_cvs_full_range, "false"))) == "true")
tryCatch({
  p_cvs <- plot_cvs(qc_se, id = "sample_name", scale = cvs_scale)
  ggplot2::ggsave(file.path(qc_dir, fn_("sample_cvs", "pdf")), p_cvs, width = 8, height = 5)
  ggplot2::ggsave(file.path(qc_dir, fn_("sample_cvs", "png")), p_cvs, width = 8, height = 5, dpi = 150)
  cat("Sample CVs plot saved\n")
}, error = function(e) {
  unlink(c(file.path(qc_dir, fn_("sample_cvs", "pdf")), file.path(qc_dir, fn_("sample_cvs", "png"))))
  warning("Skipping sample CVs plot: ", conditionMessage(e))
})
# Density (base R; can fail with sparse/all-NA data)
density_pdf <- file.path(qc_dir, fn_("density", "pdf"))
density_png <- file.path(qc_dir, fn_("density", "png"))
tryCatch({
  pdf(density_pdf, width = 8, height = 5)
  get_density(qc_se, tag = paste0(mode, " ", level))
  dev.off()
  png(density_png, width = 8, height = 5, units = "in", res = 150)
  get_density(qc_se, tag = paste0(mode, " ", level))
  dev.off()
  cat("Density plot saved\n")
}, error = function(e) {
  tryCatch(dev.off(), error = function(x) NULL)
  unlink(c(density_pdf, density_png))
  warning("Skipping density plot: ", conditionMessage(e))
})

# Feature boxplots/violins (top N variable or feature list)
# Helpers: top N variable by rowname or by rowData column
top_n_by_rownames <- function(se, n) {
  if (n <= 0 || nrow(se) == 0) return(character(0))
  v <- apply(assay(se), 1, function(x) var(x, na.rm = TRUE))
  v[is.na(v)] <- 0
  idx <- order(v, decreasing = TRUE)[seq_len(min(n, nrow(se)))]
  rownames(se)[idx]
}
top_n_by_col <- function(se, col, n) {
  if (n <= 0 || is.null(col) || !col %in% colnames(rowData(se)) || nrow(se) == 0) return(character(0))
  v <- apply(assay(se), 1, function(x) var(x, na.rm = TRUE))
  v[is.na(v)] <- 0
  idx <- order(v, decreasing = TRUE)[seq_len(min(n * 5L, nrow(se)))]
  seen <- character(0)
  for (i in idx) {
    id <- as.character(rowData(se)[[col]][i])
    if (is.na(id) || !nzchar(trimws(id)) || id %in% seen) next
    seen <- c(seen, id)
    if (length(seen) >= n) break
  }
  seen
}
cd <- as.data.frame(colData(qc_se))
assay_ids <- colnames(assay(qc_se))
feat_id <- if (all(assay_ids %in% cd$sample_name)) "sample_name" else if ("label" %in% colnames(cd) && all(assay_ids %in% cd$label)) "label" else "sample_name"
gene_col <- if (nrow(qc_se) > 0 && "Gene" %in% colnames(rowData(qc_se))) "Gene" else NULL

do_one_feature_plot <- function(se, feat, feat_dir, ptype, id, index = NULL, safe_name = NULL) {
  safe <- if (!is.null(safe_name)) gsub("[^a-zA-Z0-9_-]", "_", as.character(safe_name))[1] else gsub("[^a-zA-Z0-9_-]", "_", as.character(feat))[1]
  tryCatch({
    p_feat <- plot_feature(se, feat, index = index, type = ptype, id = id)
    ggplot2::ggsave(file.path(feat_dir, fn_(paste0(ptype, "_", safe), "pdf")), p_feat, width = 6, height = 4)
    ggplot2::ggsave(file.path(feat_dir, fn_(paste0(ptype, "_", safe), "png")), p_feat, width = 6, height = 4, dpi = 150)
  }, error = function(e) warning("Feature ", ptype, " failed for ", feat, ": ", conditionMessage(e)))
}

# Protein-level feature plots
feat_dir_prot <- file.path(comparison_dir, "feature", "protein")
dir.create(feat_dir_prot, showWarnings = FALSE, recursive = TRUE)
fl_protein <- trimws(strsplit(.or(opt$feature_list_protein, ""), "\\s*,\\s*")[[1]])
fl_protein <- fl_protein[nzchar(fl_protein)]
top_n_prot <- as.integer(.or(opt$top_n_protein, 10))
if (level == "protein" && (length(fl_protein) > 0 || top_n_prot > 0)) {
  features_prot <- if (length(fl_protein) > 0) intersect(fl_protein, rownames(qc_se)) else top_n_by_rownames(qc_se, top_n_prot)
  if (length(features_prot) > 0) {
    for (f in features_prot) {
      for (ptype in c("boxplot", "violin")) do_one_feature_plot(qc_se, f, feat_dir_prot, ptype, feat_id, index = NULL)
    }
    cat("Feature boxplots and violin plots (protein) saved (", length(features_prot), " features)\n")
  }
}

# Peptide-level feature plots
feat_dir_pep <- file.path(comparison_dir, "feature", "peptide")
fl_peptide <- trimws(strsplit(.or(opt$feature_list_peptide, ""), "\\s*,\\s*")[[1]])
fl_peptide <- fl_peptide[nzchar(fl_peptide)]
top_n_pep <- as.integer(.or(opt$top_n_peptide, 10))
if (level == "peptide" && (length(fl_peptide) > 0 || top_n_pep > 0)) {
  features_pep <- if (length(fl_peptide) > 0) intersect(fl_peptide, rownames(qc_se)) else top_n_by_rownames(qc_se, top_n_pep)
  if (length(features_pep) > 0) {
    dir.create(feat_dir_pep, showWarnings = FALSE, recursive = TRUE)
    for (f in features_pep) {
      for (ptype in c("boxplot", "violin")) do_one_feature_plot(qc_se, f, feat_dir_pep, ptype, feat_id, index = NULL)
    }
    cat("Feature boxplots and violin plots (peptide) saved (", length(features_pep), " features)\n")
  }
}

# Gene-level feature plots (when level is protein and Gene column exists)
# Pass protein IDs to plot_feature (index=NULL); use gene name for filename
fl_gene <- trimws(strsplit(.or(opt$feature_list_gene, ""), "\\s*,\\s*")[[1]])
fl_gene <- fl_gene[nzchar(fl_gene)]
top_n_gene <- as.integer(.or(opt$top_n_gene, 10))
if (level == "protein" && !is.null(gene_col) && (length(fl_gene) > 0 || top_n_gene > 0)) {
  features_gene <- if (length(fl_gene) > 0) {
    intersect(fl_gene, unique(as.character(rowData(qc_se)[[gene_col]])))
  } else {
    top_n_by_col(qc_se, gene_col, top_n_gene)
  }
  if (length(features_gene) > 0) {
    feat_dir_gene <- file.path(comparison_dir, "feature", "gene")
    dir.create(feat_dir_gene, showWarnings = FALSE, recursive = TRUE)
    for (g in features_gene) {
      prots <- rownames(qc_se)[rowData(qc_se)[[gene_col]] == g]
      if (length(prots) == 0) next
      for (ptype in c("boxplot", "violin")) do_one_feature_plot(qc_se, prots, feat_dir_gene, ptype, feat_id, index = NULL, safe_name = g)
    }
    cat("Feature boxplots and violin plots (gene) saved (", length(features_gene), " features)\n")
  }
}

# --- Sample table (LFQ: Sample Name | Experiment_Bioreplicate | condition) ---
anno <- read.table(opt$experiment_annotation, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
sample_col <- if ("sample_name" %in% colnames(anno)) "sample_name" else if ("sample" %in% colnames(anno)) "sample" else "label"
if (!sample_col %in% colnames(anno)) stop("experiment_annotation requires sample_name or sample column")
cond_col <- if ("condition" %in% colnames(anno)) "condition" else "group"
# LFQ: Sample Name = file basename (without extension); TMT: use sample_name/sample as before
sample_name_col <- if (mode == "LFQ" && "file" %in% colnames(anno)) {
  sub("\\.[^.]+$", "", basename(as.character(anno$file)))
} else {
  anno[[sample_col]]
}
sample_table <- data.frame(
  "Sample Name" = sample_name_col,
  "Experiment_Bioreplicate" = anno[[sample_col]],
  "condition" = anno[[cond_col]],
  stringsAsFactors = FALSE,
  check.names = FALSE
)
sample_fname <- paste0("SampleTable", if (nzchar(assay_suffix)) assay_suffix else "", ".csv")
write.csv(sample_table, file.path(de_dir, sample_fname), row.names = FALSE)
cat("Sample table saved:", sample_fname, "(", nrow(sample_table), "rows)\n")

# --- DE (FragPipeAnalystR: test_limma, add_rejections) ---
conditions <- sort(unique(as.character(colData(imputed_se)$condition)))
pairwise_pairs <- utils::combn(conditions, 2)
manual_contrasts <- apply(pairwise_pairs, 2, function(col) paste(col[2], col[1], sep = "_vs_"))
de_se <- test_limma(imputed_se, type = "manual", test = manual_contrasts)
de_se <- add_rejections(de_se, alpha = de_alpha, lfc = de_lfc)
contrast_cols <- grep("_diff$", colnames(rowData(de_se)), value = TRUE)
contrast_names <- gsub("_diff$", "", contrast_cols)
cat("DE: tested contrasts =", paste(contrast_names, collapse = ", "), "\n")
de_df <- cbind(as.data.frame(rowData(de_se)), as.data.frame(assay(de_se)))
de_df <- de_df[, !duplicated(colnames(de_df))]

# --- DE table additions ---
assay_mat <- as.matrix(assay(de_se))
de_df$All.mean <- rowMeans(assay_mat, na.rm = TRUE)
de_df$All.stdev <- matrixStats::rowSds(assay_mat, na.rm = TRUE)
# Per-condition mean and stdev
cd <- as.data.frame(colData(de_se))
conds <- cd$condition[match(colnames(assay_mat), rownames(cd))]
if (!any(is.na(conds))) {
  cond_to_label <- if ("condition_label" %in% colnames(anno)) {
    u <- unique(anno[, c(cond_col, "condition_label")])
    setNames(as.character(u$condition_label), as.character(u[[cond_col]]))
  } else NULL
  ucond <- unique(conds)
  for (c in ucond) {
    idx <- which(conds == c)
    lbl <- if (!is.null(cond_to_label) && c %in% names(cond_to_label)) cond_to_label[c] else c
    de_df[[paste0("Group.Mean_(", lbl, ")")]] <- rowMeans(assay_mat[, idx, drop = FALSE], na.rm = TRUE)
    if (length(idx) > 1L) {
      de_df[[paste0("Group.Stdev_(", lbl, ")")]] <- matrixStats::rowSds(assay_mat[, idx, drop = FALSE], na.rm = TRUE)
    } else {
      de_df[[paste0("Group.Stdev_(", lbl, ")")]] <- NA_real_
    }
  }
}

# --- Gene annotations merge (merge on Gene = best-matching annotation column) ---
gene_annotations <- .or(opt$gene_annotations, "")
if (nzchar(trimws(gene_annotations)) && gene_annotations != "null" && "Gene" %in% colnames(de_df)) {
  tryCatch({
    path <- trimws(gene_annotations)
    is_csv <- grepl("\\.csv$", path, ignore.case = TRUE)
    annot <- if (is_csv) read.csv(path, header = TRUE, quote = "\"", comment.char = "") else read.delim(path, sep = "\t", header = TRUE, quote = "", comment.char = "")
    # Try every annotation column; pick the one with most matches to de_df$Gene
    candidates <- colnames(annot)
    de_genes <- unique(na.omit(as.character(de_df$Gene)))
    best_col <- NULL
    best_n <- 0L
    for (col in candidates) {
      annot_vals <- unique(na.omit(as.character(annot[[col]])))
      n_match <- sum(de_genes %in% annot_vals)
      if (n_match > best_n) {
        best_n <- n_match
        best_col <- col
      }
    }
    if (!is.null(best_col) && best_n > 0L) {
      annot_merge <- annot
      annot_merge$Gene <- as.character(annot_merge[[best_col]])
      annot_merge <- annot_merge[!duplicated(annot_merge$Gene), , drop = FALSE]
      annot_cols <- c("Gene", setdiff(colnames(annot), best_col))
      de_df$..ord.. <- seq_len(nrow(de_df))
      de_df <- merge(annot_merge[, annot_cols, drop = FALSE], de_df, by = "Gene", all.y = TRUE)
      de_df <- de_df[order(de_df$..ord..), setdiff(colnames(de_df), "..ord..")]
      cat("GeneLab annotations merged into DE_results (", best_col, "→Gene, ", best_n, "/", length(de_genes), " matches)\n", sep = "")
    } else {
      warning("No annotation column matched Gene; skipping annotation merge")
    }
  }, error = function(e) warning("Could not merge annotations: ", conditionMessage(e)))
}

# --- Rename DE columns to pipeline doc expected format: Log2fc_, P.value_, Adj.p.value_, Significant_, CI.L_, CI.R_ ---
cond_label_col <- if ("condition_label" %in% colnames(anno)) "condition_label" else NULL
cond_to_label_de <- if (!is.null(cond_label_col)) {
  u <- unique(anno[, c(cond_col, cond_label_col)])
  setNames(as.character(u[[cond_label_col]]), as.character(u[[cond_col]]))
} else NULL
num_labels <- vapply(pairwise_pairs[2, ], function(cond) if (!is.null(cond_to_label_de) && cond %in% names(cond_to_label_de)) cond_to_label_de[cond] else cond, character(1))
denom_labels <- vapply(pairwise_pairs[1, ], function(cond) if (!is.null(cond_to_label_de) && cond %in% names(cond_to_label_de)) cond_to_label_de[cond] else cond, character(1))
comp_names <- paste0("(", num_labels, ")v(", denom_labels, ")")
suffix_map <- list(
  "_diff" = "Log2fc_", "_p.val" = "P.value_", "_p.adj" = "Adj.p.value_",
  "_CI.L" = "CI.L_", "_CI.R" = "CI.R_", "_significant" = "Significant_"
)
for (i in seq_along(contrast_names)) {
  idx <- match(contrast_names[i], manual_contrasts)
  if (is.na(idx)) next
  old_prefix <- contrast_names[i]
  new_prefix <- comp_names[idx]
  for (suffix in names(suffix_map)) {
    old_col <- paste0(old_prefix, suffix)
    new_col <- paste0(suffix_map[[suffix]], new_prefix)
    if (old_col %in% colnames(de_df)) {
      colnames(de_df)[colnames(de_df) == old_col] <- new_col
    }
  }
}

# --- Reorder columns to match pipeline doc: per-contrast CI.L, CI.R, Log2fc, P.value, Adj.p.value, Significant ---
stat_prefixes <- c("CI.L_", "CI.R_", "Log2fc_", "P.value_", "Adj.p.value_", "Significant_")
contrast_cols <- character()
for (cn in comp_names) {
  for (pfx in stat_prefixes) {
    ccol <- paste0(pfx, cn)
    if (ccol %in% colnames(de_df)) contrast_cols <- c(contrast_cols, ccol)
  }
}
fixed_order <- c("significant", "All.mean", "All.stdev")
group_mean_cols <- grep("^Group\\.Mean_", colnames(de_df), value = TRUE)
group_stdev_cols <- grep("^Group\\.Stdev_", colnames(de_df), value = TRUE)
other_cols <- setdiff(colnames(de_df), c(contrast_cols, fixed_order, group_mean_cols, group_stdev_cols))
de_df <- de_df[, c(other_cols, contrast_cols, fixed_order, group_mean_cols, group_stdev_cols)]

write.csv(de_df, file.path(de_dir, fn_("DE_results", "csv")), row.names = FALSE)
cat("DE_results.csv saved\n")

# --- Data table exports (only stages that changed the data); written to output_dir root ---
# Assay values: log2 for Intensity/MaxLFQ (FragPipeAnalystR make_se log2transform=T); raw for Spectral Count.
.write_table_stage <- function(se, fname, stage) {
  df <- cbind(as.data.frame(rowData(se)), as.data.frame(assay(se)))
  df <- df[, !duplicated(colnames(df))]
  write.csv(df, file.path(output_dir, fname), row.names = FALSE)
  cat(stage, "table saved:", fname, "\n")
}
sfx <- assay_suffix
.write_table_stage(data_se, paste0("raw_matrix", sfx, ".csv"), "Raw")
if (imp_type != "none") .write_table_stage(imputed_se, paste0("imputed_matrix", sfx, ".csv"), "Imputed")

# --- Contrasts table (row1=numerator, row2=denominator) ---
# Headers use condition_label when available; row entries use raw condition
cond_label_col <- if ("condition_label" %in% colnames(anno)) "condition_label" else NULL
cond_to_label <- if (!is.null(cond_label_col)) {
  u <- unique(anno[, c(cond_col, cond_label_col)])
  setNames(as.character(u[[cond_label_col]]), as.character(u[[cond_col]]))
} else NULL
num_labels <- vapply(pairwise_pairs[2, ], function(cond) if (!is.null(cond_to_label) && cond %in% names(cond_to_label)) cond_to_label[cond] else cond, character(1))
denom_labels <- vapply(pairwise_pairs[1, ], function(cond) if (!is.null(cond_to_label) && cond %in% names(cond_to_label)) cond_to_label[cond] else cond, character(1))
comp_names <- paste0("(", num_labels, ")v(", denom_labels, ")")
contrasts_df <- data.frame(row_index = c("1", "2"))
for (i in seq_along(comp_names)) {
  contrasts_df[[comp_names[i]]] <- c(pairwise_pairs[2, i], pairwise_pairs[1, i])
}
colnames(contrasts_df)[1] <- ""
contrasts_fname <- paste0("contrasts", if (nzchar(assay_suffix)) assay_suffix else "", ".csv")
write.csv(contrasts_df, file.path(de_dir, contrasts_fname), row.names = FALSE)
cat("Contrasts table saved:", contrasts_fname, "\n")

# --- Volcano plots (one per contrast; filenames use human-readable comp_names like ANCOMBC) ---
volcano_ncol <- if (volcano_name_col && "Gene" %in% colnames(rowData(de_se))) "Gene" else NULL
for (i in seq_along(contrast_names)) {
  tryCatch({
    p_v <- plot_volcano(de_se, contrast_names[i], name_col = volcano_ncol, add_names = volcano_add_names, alpha = de_alpha, lfc = de_lfc)
    fname_base <- gsub("[[:space:]]+", "_", comp_names[i])
    fname_base <- paste0(fname_base, if (nzchar(assay_suffix)) assay_suffix else "", "_volcano")
    ggplot2::ggsave(file.path(volcano_dir, paste0(fname_base, ".pdf")), p_v, width = 8, height = 6)
    ggplot2::ggsave(file.path(volcano_dir, paste0(fname_base, ".png")), p_v, width = 8, height = 6, dpi = 150)
  }, error = function(e) warning("Volcano for ", contrast_names[i], " failed: ", conditionMessage(e)))
}
cat("Volcano plots saved (", length(contrast_names), " contrasts)\n")

# --- DE heatmap (FragPipeAnalystR: get_cluster_heatmap, type=centered) ---
# get_cluster_heatmap maps colnames(df) via temp[colnames(df), "sample_name"] with rownames(temp)=label.
# make_se_from_files sets assay colnames=sample_name, so lookup fails (NA). Pass a copy with assay colnames=label.
tryCatch({
  de_se_hm <- de_se
  if ("label" %in% colnames(colData(de_se))) {
    cd <- as.data.frame(colData(de_se))
    rownames(cd) <- cd$label
    a <- assay(de_se)
    colnames(a) <- cd$label
    de_se_hm <- SummarizedExperiment::SummarizedExperiment(
      assays = list(a), colData = S4Vectors::DataFrame(cd), rowData = rowData(de_se), metadata = metadata(de_se))
  }
  ht_res <- get_cluster_heatmap(de_se_hm, type = "centered", indicate = "condition",
    alpha = de_alpha, lfc = de_lfc, col_limit = 6, plot = TRUE)
  if (inherits(ht_res, "list") && inherits(ht_res[[1]], "Heatmap")) {
    ht <- ht_res[[1]]
    pdf(file.path(de_dir, fn_("DE_heatmap", "pdf")), width = 10, height = 8)
    ComplexHeatmap::draw(ht, heatmap_legend_side = "top")
    dev.off()
    png(file.path(de_dir, fn_("DE_heatmap", "png")), width = 10, height = 8, units = "in", res = 150)
    ComplexHeatmap::draw(ht, heatmap_legend_side = "top")
    dev.off()
    cat("DE heatmap saved\n")
  } else     if (inherits(ht_res, "gg")) {
    ggplot2::ggsave(file.path(de_dir, fn_("DE_heatmap", "pdf")), ht_res, width = 8, height = 4)
    ggplot2::ggsave(file.path(de_dir, fn_("DE_heatmap", "png")), ht_res, width = 8, height = 4, dpi = 150)
    cat("DE heatmap: no significant features (empty plot saved)\n")
  }
}, error = function(e) warning("DE heatmap failed: ", conditionMessage(e)))

# --- Enrichment (FragPipeAnalystR: or_test, plot_or) ---
# or_test expects display names (e.g. "GO Biological Process"), not Enrichr lib names (e.g. "GO_Biological_Process_2021")
or_db_map <- c(
  GO_Biological_Process_2021 = "GO Biological Process",
  GO_Cellular_Component_2021 = "GO Cellular Component",
  GO_Molecular_Function_2021 = "GO Molecular Function",
  MSigDB_Hallmark_2020 = "Hallmark",
  KEGG_2021_Human = "KEGG",
  Reactome_2022 = "Reactome"
)
for (db in enrichment_dbs) {
  db_or <- if (db %in% names(or_db_map)) or_db_map[db] else db
  for (dir in enrichment_dirs) {
    tryCatch({
      or_res <- or_test(de_se, database = db_or, direction = toupper(dir), alpha = de_alpha, log2_threshold = de_lfc)
      if (!is.null(or_res) && nrow(or_res) > 0) {
        safe_name <- paste0("or_", gsub("[^A-Za-z0-9_-]", "_", db), "_", tolower(dir))
        write.csv(or_res, file.path(enr_dir_or, paste0(safe_name, ".csv")), row.names = FALSE)
        p_or <- plot_or(or_res, number = 15, alpha = de_alpha)
        if (!is.null(p_or)) {
          ggplot2::ggsave(file.path(enr_dir_or, paste0(safe_name, ".pdf")), p_or, width = 10, height = 6)
          ggplot2::ggsave(file.path(enr_dir_or, paste0(safe_name, ".png")), p_or, width = 10, height = 6, dpi = 150)
        }
        cat("Enrichment saved:", db, dir, "\n")
      } else {
        cat("No enrichment for", db, dir, "(0 genes)\n")
      }
    }, error = function(e) warning("Enrichment failed (", db, " ", dir, "): ", conditionMessage(e)))
  }
}
if (length(enrichment_dbs) > 0) cat("Enrichment complete\n")

# --- GSEA (FragPipeAnalystR: GSEA_test, plot_GSEA; protein/gene/site only, no peptide) ---
gsea_db_map <- c(
  GO_Biological_Process_2021 = "GO Biological Process",
  GO_Cellular_Component_2021 = "GO Cellular Component",
  GO_Molecular_Function_2021 = "GO Molecular Function",
  MSigDB_Hallmark_2020 = "Hallmark",
  Hallmark = "Hallmark",
  KEGG_2021_Human = "KEGG",
  KEGG = "KEGG"
)
gsea_dbs_valid <- setdiff(gsea_dbs, c("Reactome", "Reactome_2022"))  # GSEA_test: no Reactome
if (level != "peptide" && length(gsea_dbs_valid) > 0) {
  if (!requireNamespace("clusterProfiler", quietly = TRUE)) {
    warning("clusterProfiler not installed; skipping GSEA")
    gsea_dbs_valid <- character(0)
  } else {
    suppressPackageStartupMessages(library(clusterProfiler))
    need_org <- any(gsea_dbs_valid %in% c("GO_Biological_Process_2021", "GO_Cellular_Component_2021", "GO_Molecular_Function_2021", "KEGG", "KEGG_2021_Human", "GO Biological Process", "GO Cellular Component", "GO Molecular Function"))
    if (need_org && !requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
      warning("org.Hs.eg.db not installed; skipping GO/KEGG GSEA (Hallmark may still run)")
      gsea_dbs_valid <- setdiff(gsea_dbs_valid, c("GO_Biological_Process_2021", "GO_Cellular_Component_2021", "GO_Molecular_Function_2021", "KEGG", "KEGG_2021_Human"))
    } else if (need_org) {
      suppressPackageStartupMessages(library(org.Hs.eg.db, character.only = TRUE))
    }
  }
}
if (level != "peptide" && length(gsea_dbs_valid) > 0) {
  rd <- as.data.frame(rowData(de_se))
  has_gene <- "Gene" %in% colnames(rd)
  if (!has_gene) {
    warning("No Gene column; skipping GSEA")
    gsea_dbs_valid <- character(0)
  }
  for (db in gsea_dbs_valid) {
    db_gsea <- if (db %in% names(gsea_db_map)) gsea_db_map[db] else db
    for (i in seq_along(contrast_names)) {
      col_stat <- paste0(contrast_names[i], "_diff")
      if (!col_stat %in% colnames(rowData(de_se))) next
      rd <- as.data.frame(rowData(de_se))
      rd$ID <- sapply(strsplit(as.character(rd$Gene), ";"), function(x) trimws(x[1]))
      keep <- !is.na(rd$ID) & nzchar(trimws(rd$ID))
      rd <- rd[keep, , drop = FALSE]
      if (nrow(rd) == 0) next
      if (any(duplicated(rd$ID))) {
        ord <- order(rd$ID, -abs(rd[[col_stat]]), na.last = TRUE)
        rd <- rd[ord, ]
        rd <- rd[!duplicated(rd$ID), ]
      }
      idx <- intersect(rownames(rd), rownames(de_se))
      if (length(idx) == 0) next
      de_se_gsea <- de_se[idx, ]
      rowData(de_se_gsea)$ID <- rd$ID[match(idx, rownames(rd))]
      safe_contrast <- gsub("[^A-Za-z0-9_-]", "_", contrast_names[i])
      safe_name <- paste0("gsea_", gsub("[^A-Za-z0-9_-]", "_", db), "_", safe_contrast)
      tryCatch({
        gsea_res <- GSEA_test(de_se_gsea, col = col_stat, database = db_gsea, convert = TRUE)
        if (!is.null(gsea_res) && nrow(gsea_res) > 0) {
          write.csv(gsea_res, file.path(enr_dir_gsea, paste0(safe_name, ".csv")), row.names = FALSE)
          gsea_plot <- gsea_res[!duplicated(gsea_res$ID), , drop = FALSE]
          n_cat <- min(15L, max(1L, floor(nrow(gsea_plot) / 2)))
          p_gsea <- tryCatch(plot_GSEA(gsea_plot, categroies = n_cat), error = function(e) NULL)
          if (!is.null(p_gsea)) {
            ggplot2::ggsave(file.path(enr_dir_gsea, paste0(safe_name, ".pdf")), p_gsea, width = 10, height = 6)
            ggplot2::ggsave(file.path(enr_dir_gsea, paste0(safe_name, ".png")), p_gsea, width = 10, height = 6, dpi = 150)
          }
          cat("GSEA saved:", db, contrast_names[i], "\n")
        } else {
          cat("GSEA: no results for", db, contrast_names[i], "\n")
        }
      }, error = function(e) warning("GSEA failed (", db, " ", contrast_names[i], "): ", conditionMessage(e)))
    }
  }
  cat("GSEA complete\n")
}

