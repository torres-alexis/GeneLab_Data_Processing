#!/usr/bin/Rscript
# FragPipe-Analyst downstream analysis (refactored).
# Original full script in Reference_Repos/Old_implementation/
# Phase 1: params logger – parse all NF-passed args and write to fp_analyst_parameters.txt

library(optparse)

option_list <- list(
  make_option(c("--experiment_annotation"), type = "character", default = NULL,
    help = "Path to experiment annotation TSV", metavar = "FILE"),
  make_option(c("--quantification_file"), type = "character", default = NULL,
    help = "Path to quantification file", metavar = "FILE"),
  make_option(c("--mode"), type = "character", default = NULL,
    help = "Quantification mode: LFQ, TMT, or DIA", metavar = "MODE"),
  make_option(c("--level"), type = "character", default = NULL,
    help = "Analysis level: protein, peptide, gene, or site", metavar = "LEVEL"),
  make_option(c("--output_dir"), type = "character", default = "output/",
    help = "Output directory", metavar = "DIR"),
  make_option(c("--lfq_type"), type = "character", default = "Intensity",
    help = "LFQ column type: Intensity or MaxLFQ (LFQ mode only)", metavar = "STRING"),
  make_option(c("--normalization_method"), type = "character", default = "none",
    help = "Normalization: none, vsn, MD, or GN", metavar = "STRING"),
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
  make_option(c("--pathway_database"), type = "character", default = "",
    help = "Pathway DB(s): comma-separated. Hallmark, KEGG, KEGG Mouse, Reactome, WikiPathways Mouse, or any Enrichr libraryName (maayanlab.cloud/Enrichr/datasetStatistics). Empty = skip", metavar = "STRING"),
  make_option(c("--pathway_direction"), type = "character", default = "Both",
    help = "Pathway direction: Up, Down, Both", metavar = "STRING"),
  make_option(c("--go_database"), type = "character", default = "",
    help = "GO DB(s): comma-separated. GO Biological Process, GO Cellular Component, GO Molecular Function, or any Enrichr libraryName. Empty = skip", metavar = "STRING"),
  make_option(c("--go_direction"), type = "character", default = "Both",
    help = "GO direction: Up, Down, Both", metavar = "STRING"),
  make_option(c("--de_alpha"), type = "numeric", default = 0.05,
    help = "Adjusted p-value threshold for DE significance", metavar = "NUMERIC"),
  make_option(c("--de_lfc"), type = "numeric", default = 1.0,
    help = "Log2 fold change threshold for DE significance", metavar = "NUMERIC"),
  make_option(c("--de_fdr"), type = "character", default = "Benjamini Hochberg",
    help = "FDR correction: Benjamini Hochberg or Local and tail area-based", metavar = "STRING"),
  make_option(c("--imputation_type"), type = "character", default = "Perseus-type",
    help = "Imputation: none, Perseus-type, knn, MLE, min, zero, bpca, QRILC, MinDet, MinProb, RF, nbavg, mixed", metavar = "STRING"),
  make_option(c("--min_global_appearance"), type = "numeric", default = 0,
    help = "Min %% non-missing across all samples (0-100)", metavar = "NUMERIC"),
  make_option(c("--min_appearance_one_condition"), type = "numeric", default = 0,
    help = "Min %% non-missing in at least one condition (0-100)", metavar = "NUMERIC"),
  make_option(c("--qc_show_imputed"), type = "character", default = "true",
    help = "Use imputed data for QC plots: true/false", metavar = "STRING"),
  make_option(c("--qc_include_both"), type = "character", default = "false",
    help = "Generate both imputed and unimputed QC: true/false", metavar = "STRING"),
  make_option(c("--volcano_display_names"), type = "character", default = "true",
    help = "Volcano: display names on significant points: true/false", metavar = "STRING"),
  make_option(c("--volcano_show_gene"), type = "character", default = "true",
    help = "Volcano: show gene names: true/false", metavar = "STRING"),
  make_option(c("--volcano_highlight_feature"), type = "character", default = "",
    help = "Comma-delimited features to highlight on volcano", metavar = "STRING"),
  make_option(c("--volcano_show_other_peptides"), type = "character", default = "true",
    help = "Peptide/site volcano: show other peptides from same protein (blue) when highlighting. false = don't color", metavar = "STRING")
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
level <- .oneof(opt$level, c("protein", "peptide", "gene", "site"), "level")

output_dir <- opt$output_dir
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
qc_dir <- file.path(output_dir, "qc")
dir.create(qc_dir, showWarnings = FALSE, recursive = TRUE)
comparison_dir <- file.path(output_dir, "comparison")
dir.create(comparison_dir, showWarnings = FALSE, recursive = TRUE)

# Log raw inputs (reproducibility)
param_path <- file.path(output_dir, "fp_analyst_parameters.txt")
lines <- c(
  "FragPipe-Analyst Parameters (raw inputs)",
  "=======================================",
  paste("experiment_annotation:", opt$experiment_annotation),
  paste("quantification_file:", opt$quantification_file),
  paste("mode:", opt$mode),
  paste("level:", opt$level),
  paste("output_dir:", output_dir),
  paste("lfq_type:", opt$lfq_type),
  paste("normalization_method:", opt$normalization_method),
  paste("feature_list_protein:", opt$feature_list_protein),
  paste("feature_list_gene:", opt$feature_list_gene),
  paste("feature_list_peptide:", opt$feature_list_peptide),
  paste("feature_list_site:", opt$feature_list_site),
  paste("top_n_protein:", opt$top_n_protein),
  paste("top_n_gene:", opt$top_n_gene),
  paste("top_n_peptide:", opt$top_n_peptide),
  paste("top_n_site:", opt$top_n_site),
  paste("pathway_database:", opt$pathway_database),
  paste("pathway_direction:", opt$pathway_direction),
  paste("go_database:", opt$go_database),
  paste("go_direction:", opt$go_direction),
  paste("de_alpha:", opt$de_alpha),
  paste("de_lfc:", opt$de_lfc),
  paste("de_fdr:", opt$de_fdr),
  paste("imputation_type:", opt$imputation_type),
  paste("min_global_appearance:", opt$min_global_appearance),
  paste("min_appearance_one_condition:", opt$min_appearance_one_condition),
  paste("qc_show_imputed:", opt$qc_show_imputed),
  paste("qc_include_both:", opt$qc_include_both),
  paste("volcano_display_names:", opt$volcano_display_names),
  paste("volcano_show_gene:", opt$volcano_show_gene),
  paste("volcano_highlight_feature:", opt$volcano_highlight_feature),
  paste("volcano_show_other_peptides:", opt$volcano_show_other_peptides)
)
writeLines(lines, param_path)

# Normalize/coerce to typed variables (strict: only listed options, tolower for case)
feature_list_protein <- if (nzchar(trimws(.or(opt$feature_list_protein, "")))) {
  trimws(strsplit(trimws(opt$feature_list_protein), "\\s*,\\s*")[[1]])
} else character(0)
feature_list_gene <- if (nzchar(trimws(.or(opt$feature_list_gene, "")))) {
  trimws(strsplit(trimws(opt$feature_list_gene), "\\s*,\\s*")[[1]])
} else character(0)
feature_list_peptide <- if (nzchar(trimws(.or(opt$feature_list_peptide, "")))) {
  trimws(strsplit(trimws(opt$feature_list_peptide), "\\s*,\\s*")[[1]])
} else character(0)
feature_list_site <- if (nzchar(trimws(.or(opt$feature_list_site, "")))) {
  trimws(strsplit(trimws(opt$feature_list_site), "\\s*,\\s*")[[1]])
} else character(0)
top_n_protein <- max(0L, as.integer(.or(opt$top_n_protein, 10)))
top_n_gene <- max(0L, as.integer(.or(opt$top_n_gene, 10)))
top_n_peptide <- max(0L, as.integer(.or(opt$top_n_peptide, 10)))
top_n_site <- max(0L, as.integer(.or(opt$top_n_site, 10)))

pathway_database_raw <- trimws(.or(opt$pathway_database, ""))
pathway_database <- if (!nzchar(pathway_database_raw)) character(0) else {
  unique(trimws(strsplit(pathway_database_raw, "\\s*,\\s*")[[1]]))
}
pathway_database <- pathway_database[nzchar(pathway_database)]
pathway_direction <- .oneof(.or(opt$pathway_direction, "Both"), c("Up", "Down", "Both"), "pathway_direction")
go_database_raw <- trimws(.or(opt$go_database, ""))
go_database <- if (!nzchar(go_database_raw)) character(0) else {
  unique(trimws(strsplit(go_database_raw, "\\s*,\\s*")[[1]]))
}
go_database <- go_database[nzchar(go_database)]
go_direction <- .oneof(.or(opt$go_direction, "Both"), c("Up", "Down", "Both"), "go_direction")

lfq_type <- .oneof(.or(opt$lfq_type, "Intensity"), c("Intensity", "MaxLFQ", "Spectral Count"), "lfq_type")
de_alpha <- as.numeric(.or(opt$de_alpha, 0.05))
de_lfc <- as.numeric(.or(opt$de_lfc, 1.0))
de_fdr_raw <- tolower(trimws(.or(opt$de_fdr, "Benjamini Hochberg")))
de_fdr <- if (de_fdr_raw %in% c("benjamini hochberg", "bh", "benjamini-hochberg")) "Benjamini Hochberg" else
  if (de_fdr_raw %in% c("local and tail area-based", "fdrtool", "local_tail")) "Local and tail area-based" else
  stop("Invalid de_fdr: '", opt$de_fdr, "'. Valid: Benjamini Hochberg, Local and tail area-based")

imputation_valid <- c("Perseus-type", "knn", "MLE", "min", "zero", "bpca", "QRILC", "MinDet", "MinProb", "RF", "nbavg", "mixed")
imputation_type_raw <- trimws(.or(opt$imputation_type, "Perseus-type"))
imputation_type <- if (tolower(imputation_type_raw) == "none") "none" else
  if (tolower(imputation_type_raw) %in% c("man", "perseus", "perseus-type", "perseus_type")) "Perseus-type" else
  if (tolower(imputation_type_raw) == "mle") "MLE" else
  if (imputation_type_raw %in% imputation_valid) imputation_type_raw else
  stop("Invalid imputation_type: '", opt$imputation_type, "'. Valid: none, ", paste(imputation_valid, collapse = ", "))

normalization_valid <- c("none", "vsn", "MD", "GN")
normalization_method_raw <- tolower(trimws(.or(opt$normalization_method, "none")))
normalization_method <- if (normalization_method_raw %in% c("none", "")) "none" else
  if (normalization_method_raw == "vsn") "vsn" else
  if (normalization_method_raw == "md") "MD" else
  if (normalization_method_raw == "gn") "GN" else
  stop("Invalid normalization_method: '", opt$normalization_method, "'. Valid: ", paste(normalization_valid, collapse = ", "))

min_global_appearance <- as.numeric(.or(opt$min_global_appearance, 0))
min_appearance_one_condition <- as.numeric(.or(opt$min_appearance_one_condition, 0))

qc_show_imputed <- .oneof(.or(opt$qc_show_imputed, "true"), c("true", "false"), "qc_show_imputed") == "true"
qc_include_both <- .oneof(.or(opt$qc_include_both, "false"), c("true", "false"), "qc_include_both") == "true"
volcano_display_names <- .oneof(.or(opt$volcano_display_names, "true"), c("true", "false"), "volcano_display_names") == "true"
volcano_show_gene <- .oneof(.or(opt$volcano_show_gene, "true"), c("true", "false"), "volcano_show_gene") == "true"
volcano_highlight_feature <- if (!nzchar(trimws(.or(opt$volcano_highlight_feature, "")))) character(0) else
  trimws(strsplit(trimws(opt$volcano_highlight_feature), "\\s*,\\s*")[[1]])
volcano_show_other_peptides <- .oneof(.or(opt$volcano_show_other_peptides, "true"), c("true", "false"), "volcano_show_other_peptides") == "true"

# Print parsed (for testing)
cat("--- Parsed (normalized) ---\n")
cat("experiment_annotation:", opt$experiment_annotation, "\n")
cat("quantification_file:", opt$quantification_file, "\n")
cat("mode:", mode, "| level:", level, "| lfq_type:", lfq_type, "\n")
cat("normalization_method:", normalization_method, "| imputation_type:", imputation_type, "\n")
cat("de_alpha:", de_alpha, "| de_lfc:", de_lfc, "| de_fdr:", de_fdr, "\n")
cat("feature_list_protein:", paste(feature_list_protein, collapse = ", "), "\n")
cat("feature_list_gene:", paste(feature_list_gene, collapse = ", "), "\n")
cat("feature_list_peptide:", paste(feature_list_peptide, collapse = ", "), "\n")
cat("feature_list_site:", paste(feature_list_site, collapse = ", "), "\n")
cat("top_n_protein:", top_n_protein, "| top_n_gene:", top_n_gene, "| top_n_peptide:", top_n_peptide, "| top_n_site:", top_n_site, "\n")
cat("qc_show_imputed:", qc_show_imputed, "| qc_include_both:", qc_include_both, "\n")
cat("volcano_display_names:", volcano_display_names, "| volcano_show_gene:", volcano_show_gene, "\n")
cat("volcano_highlight_feature:", paste(volcano_highlight_feature, collapse = ", "), "| volcano_show_other_peptides:", volcano_show_other_peptides, "\n")

# Quick summary
anno <- read.table(opt$experiment_annotation, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
quant <- read.table(opt$quantification_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE,
  fill = TRUE, comment.char = "", blank.lines.skip = FALSE, check.names = FALSE)
cat("annotation:", nrow(anno), "rows,", ncol(anno), "cols | quantification:", nrow(quant), "rows,", ncol(quant), "cols\n")
cat("Params logged to:", param_path, "\n")

# Phase 2: Create SummarizedExperiment (fp_se_helper.R)
script_dir <- if (length(grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)) > 0) {
  dirname(sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1]))
} else "."
source(file.path(script_dir, "fp_se_helper.R"))
source(file.path(script_dir, "fp_plot_helper.R"))
source(file.path(script_dir, "fp_de_helper.R"))
source(file.path(script_dir, "fp_enrichment_helper.R"))
data_se <- make_se_from_files(
  quant_table_path = opt$quantification_file,
  exp_anno_path = opt$experiment_annotation,
  type = mode,
  level = level,
  lfq_type = lfq_type
)
if (is.null(data_se)) stop("SE creation returned NULL")
cat("SE created:", nrow(data_se), "features x", ncol(data_se), "samples\n")

# Phase 3: Export Raw_matrix (raw intensities before filter/norm/impute), apply filter, export Filtered_matrix
original_se <- data_se
original_df <- cbind(as.data.frame(rowData(data_se)), as.data.frame(assay(data_se)))
write.csv(original_df, file.path(output_dir, "Raw_matrix.csv"), row.names = FALSE)
cat("Raw_matrix.csv:", nrow(original_df), "features\n")

filtered_se <- data_se
if (min_global_appearance > 0) {
  filtered_se <- global_filter(filtered_se, 100 - min_global_appearance)
  cat("global_filter: kept", nrow(filtered_se), "features (min", min_global_appearance, "% non-missing globally)\n")
}
if (min_appearance_one_condition > 0) {
  filtered_se <- filter_by_condition(filtered_se, min_appearance_one_condition)
  cat("filter_by_condition: kept", nrow(filtered_se), "features (min", min_appearance_one_condition, "% in at least one condition)\n")
}

filtered_df <- cbind(as.data.frame(rowData(filtered_se)), as.data.frame(assay(filtered_se)))
write.csv(filtered_df, file.path(output_dir, "Filtered_matrix.csv"), row.names = FALSE)
cat("Filtered_matrix.csv:", nrow(filtered_df), "features\n")

# Phase 4: Normalize, export Normalized_matrix (only when normalization applied; else redundant with Filtered_matrix)
normalized_se <- normalize_se(filtered_se, normalization_method)
if (normalization_method != "none") {
  cat("normalize_se:", normalization_method, "applied\n")
  normalized_df <- cbind(as.data.frame(rowData(normalized_se)), as.data.frame(assay(normalized_se)))
  write.csv(normalized_df, file.path(output_dir, "Normalized_matrix.csv"), row.names = FALSE)
  cat("Normalized_matrix.csv:", nrow(normalized_df), "features\n")
}

# Phase 5: Impute (if requested), export Imputed_matrix
imputed_se <- NULL
if (imputation_type != "none") {
  imputed_se <- impute_se(normalized_se, fun = imputation_type)
  imputed_df <- cbind(as.data.frame(rowData(imputed_se)), as.data.frame(assay(imputed_se)))
  write.csv(imputed_df, file.path(output_dir, "Imputed_matrix.csv"), row.names = FALSE)
  cat("impute_se:", imputation_type, "applied | Imputed_matrix.csv:", nrow(imputed_df), "features\n")
  data_se <- imputed_se
} else {
  data_se <- normalized_se
}

# Phase 6: QC plots (PCA, correlation, missing heatmap, feature numbers, coverage, density)
has_missing <- any(is.na(assay(filtered_se)))
se_for_qc <- if (qc_show_imputed && !is.null(imputed_se)) imputed_se else normalized_se
qc_versions <- if (qc_include_both && !is.null(imputed_se)) {
  list(
    list(se = imputed_se, suffix = "_imputed", desc = "imputed"),
    list(se = normalized_se, suffix = "_unimputed", desc = "unimputed")
  )
} else {
  list(list(se = se_for_qc, suffix = "", desc = ""))
}

cat("QC plots using:", if (qc_show_imputed && !is.null(imputed_se)) "imputed" else "unimputed", "data\n")

for (vv in qc_versions) {
  tryCatch({
    n_rep <- if ("replicate" %in% colnames(colData(vv$se)))
      length(unique(as.character(colData(vv$se)$replicate))) else 0L
    pca_indicate <- if (n_rep >= 1L && n_rep <= 6L) c("condition", "replicate") else "condition"
    p_pca <- plot_pca_custom(vv$se, indicate = pca_indicate, plot = TRUE)
    base <- paste0("pca", vv$suffix)
    ggplot2::ggsave(file.path(qc_dir, paste0(base, ".pdf")), p_pca, width = 8, height = 6)
    ggplot2::ggsave(file.path(qc_dir, paste0(base, ".png")), p_pca, width = 8, height = 6, dpi = 150)
    cat("PCA saved:", base, vv$desc, "\n")
  }, error = function(e) warning("PCA failed: ", conditionMessage(e)))
}

for (vv in qc_versions) {
  tryCatch({
    n_samp <- ncol(vv$se)
    ht_corr <- plot_cor_customized(vv$se, indicate = "condition", plot = FALSE,
      font_size = if (n_samp > 10) 9 else 12)
    base <- paste0("correlation_heatmap", vv$suffix)
    fig_side <- min(14, max(7, 5 + n_samp * 0.35))
    pdf(file.path(comparison_dir, paste0(base, ".pdf")), width = fig_side, height = fig_side)
    ComplexHeatmap::draw(ht_corr, heatmap_legend_side = "top")
    dev.off()
    png(file.path(comparison_dir, paste0(base, ".png")), width = fig_side, height = fig_side, units = "in", res = 150)
    ComplexHeatmap::draw(ht_corr, heatmap_legend_side = "top")
    dev.off()
    cat("Correlation heatmap saved:", base, "\n")
  }, error = function(e) warning("Correlation heatmap failed: ", conditionMessage(e)))
}

if (has_missing) {
  tryCatch({
    pdf(file.path(qc_dir, "missing_value_heatmap.pdf"), width = 8, height = 6)
    plot_missval_customized(filtered_se)
    dev.off()
    png(file.path(qc_dir, "missing_value_heatmap.png"), width = 8, height = 6, units = "in", res = 150)
    plot_missval_customized(filtered_se)
    dev.off()
    cat("Missing value heatmap saved\n")
  }, error = function(e) warning("Missing value heatmap failed: ", conditionMessage(e)))
}

tryCatch({
  p_fn <- plot_feature_numbers_custom(filtered_se)
  ggplot2::ggsave(file.path(qc_dir, "feature_numbers.pdf"), p_fn, width = 8, height = 5)
  ggplot2::ggsave(file.path(qc_dir, "feature_numbers.png"), p_fn, width = 8, height = 5, dpi = 150)
  cat("Feature numbers plot saved\n")
}, error = function(e) warning("Feature numbers failed: ", conditionMessage(e)))

tryCatch({
  p_cov <- plot_coverage_customized(filtered_se, plot = TRUE)
  ggplot2::ggsave(file.path(qc_dir, "sample_coverage.pdf"), p_cov, width = 8, height = 5)
  ggplot2::ggsave(file.path(qc_dir, "sample_coverage.png"), p_cov, width = 8, height = 5, dpi = 150)
  cat("Sample coverage plot saved\n")
}, error = function(e) warning("Sample coverage failed: ", conditionMessage(e)))

for (vv in qc_versions) {
  tryCatch({
    p_cvs <- plot_cvs_custom(vv$se, id = "label", scale = TRUE)
    base <- paste0("sample_cvs", vv$suffix)
    ggplot2::ggsave(file.path(qc_dir, paste0(base, ".pdf")), p_cvs, width = 8, height = 6)
    ggplot2::ggsave(file.path(qc_dir, paste0(base, ".png")), p_cvs, width = 8, height = 6, dpi = 150)
    cat("Sample CVs saved:", base, "\n")
  }, error = function(e) warning("Sample CVs failed: ", conditionMessage(e)))
}

tryCatch({
  ses_dens <- list("original" = original_se, "filtered" = filtered_se)
  if (!is.null(imputed_se)) ses_dens[["imputed"]] <- imputed_se
  p_dens <- plot_density_custom(ses_dens)
  ggplot2::ggsave(file.path(qc_dir, "density.pdf"), p_dens, width = 8, height = 7)
  ggplot2::ggsave(file.path(qc_dir, "density.png"), p_dens, width = 8, height = 7, dpi = 150)
  cat("Density plot saved\n")
}, error = function(e) warning("Density plot failed: ", conditionMessage(e)))

# Phase 7: Comparison plots (Jaccard, Venn, UpSet, feature)
n_conditions <- length(unique(colData(filtered_se)$condition))
if (n_conditions >= 2 && mode %in% c("LFQ", "DIA")) {
  exp <- if (!is.null(metadata(filtered_se)$exp)) metadata(filtered_se)$exp else mode
  att_df <- data_attendance_custom(filtered_se, exp = exp, level = level)
  conditions <- unique(colData(filtered_se)$condition)
  tryCatch({
    pdf(file.path(comparison_dir, "jaccard.pdf"), width = 7, height = 6)
    plot_Jaccard_custom(filtered_se, plot = TRUE)
    dev.off()
    png(file.path(comparison_dir, "jaccard.png"), width = 7, height = 6, units = "in", res = 150)
    plot_Jaccard_custom(filtered_se, plot = TRUE)
    dev.off()
    cat("Jaccard saved\n")
  }, error = function(e) warning("Jaccard failed: ", conditionMessage(e)))
  vd_dir <- file.path(comparison_dir, "venndiagram")
  dir.create(vd_dir, showWarnings = FALSE, recursive = TRUE)
  for (pr in utils::combn(conditions, 2, simplify = FALSE)) {
    tryCatch({
      v <- plot_venn_custom(att_df, pr[1], pr[2], cond3 = NULL)
      if (!is.null(v)) {
        safe_name <- paste0(gsub("[^A-Za-z0-9_-]", "_", pr[1]), "_vs_", gsub("[^A-Za-z0-9_-]", "_", pr[2]))
        ggplot2::ggsave(file.path(vd_dir, paste0("venn_", safe_name, ".pdf")), v, width = 6, height = 6)
        ggplot2::ggsave(file.path(vd_dir, paste0("venn_", safe_name, ".png")), v, width = 6, height = 6, dpi = 150)
        cat("Venn saved:", paste(pr, collapse = " vs "), "\n")
      }
    }, error = function(e) warning("Venn failed: ", conditionMessage(e)))
  }
  if (length(grep("Occurences_", colnames(att_df))) >= 2) {
    tryCatch({
      pdf(file.path(comparison_dir, "upset.pdf"), width = 10, height = 6, onefile = FALSE)
      plot_upset_custom(att_df)
      dev.off()
      png(file.path(comparison_dir, "upset.png"), width = 10, height = 6, units = "in", res = 150)
      plot_upset_custom(att_df)
      dev.off()
      cat("UpSet saved\n")
    }, error = function(e) warning("UpSet failed: ", conditionMessage(e)))
  }
}

# Feature plots (boxplot/violin by protein, gene, peptide, site)
do_feature_plots <- function(features, feat_index, subdir) {
  if (length(features) == 0) return(invisible(NULL))
  for (feat in features) {
    if (is.na(feat) || !nzchar(trimws(as.character(feat)))) next
    feat <- as.character(feat)
    safe_name <- gsub("[^A-Za-z0-9_-]", "_", feat)
    for (vv in qc_versions) {
      prots <- if (is.null(feat_index)) feat else rownames(vv$se)[rowData(vv$se)[[feat_index]] == feat]
      if (length(prots) == 0) next
      for (ptype in c("boxplot", "violinplot")) {
        tryCatch({
          p_f <- plot_feature_custom(vv$se, prots, type = sub("plot$", "", ptype), show_gene = !is.null(feat_index))
          feat_dir <- file.path(comparison_dir, "feature", subdir, ptype)
          dir.create(feat_dir, showWarnings = FALSE, recursive = TRUE)
          base <- paste0(ptype, "_feature_", safe_name, vv$suffix)
          ggplot2::ggsave(file.path(feat_dir, paste0(base, ".pdf")), p_f, width = 6, height = 4)
          ggplot2::ggsave(file.path(feat_dir, paste0(base, ".png")), p_f, width = 6, height = 4, dpi = 150)
          cat("Feature plot saved:", file.path("feature", subdir, ptype, base), "\n")
        }, error = function(e) warning("Feature plot failed: ", conditionMessage(e)))
      }
    }
  }
}
lvl <- metadata(se_for_qc)$level
if (is.null(lvl)) lvl <- "protein"
gene_col <- if (nrow(se_for_qc) > 0 && "Gene" %in% colnames(rowData(se_for_qc))) "Gene" else NULL

# Helper: top N variable by rowname
top_n_by_rownames <- function(se, n, mult = 2L) {
  if (n <= 0 || nrow(se) == 0) return(character(0))
  vars <- apply(assay(se), 1, function(x) var(x, na.rm = TRUE))
  vars[is.na(vars)] <- 0
  top_idx <- order(vars, decreasing = TRUE)[seq_len(min(n * mult, length(vars)))]
  seen <- character(0)
  for (i in top_idx) {
    id <- rownames(se)[i]
    if (is.na(id) || !nzchar(id) || id %in% seen) next
    seen <- c(seen, id)
    if (length(seen) >= n) break
  }
  seen
}
# Helper: top N variable by column
top_n_by_col <- function(se, col, n, mult = 10L) {
  if (n <= 0 || is.null(col) || nrow(se) == 0) return(character(0))
  vars <- apply(assay(se), 1, function(x) var(x, na.rm = TRUE))
  vars[is.na(vars)] <- 0
  top_idx <- order(vars, decreasing = TRUE)[seq_len(min(n * mult, length(vars)))]
  seen <- character(0)
  for (i in top_idx) {
    id <- as.character(rowData(se)[[col]][i])
    if (is.na(id) || !nzchar(trimws(id)) || id %in% seen) next
    seen <- c(seen, id)
    if (length(seen) >= n) break
  }
  seen
}

if (lvl == "protein") {
  features_protein <- if (length(feature_list_protein) > 0) trimws(feature_list_protein[nzchar(trimws(feature_list_protein))]) else top_n_by_rownames(se_for_qc, top_n_protein)
  if (length(features_protein) > 0) do_feature_plots(features_protein, feat_index = NULL, subdir = "protein")
  features_gene <- if (length(feature_list_gene) > 0) trimws(feature_list_gene[nzchar(trimws(feature_list_gene))]) else top_n_by_col(se_for_qc, gene_col, top_n_gene)
  if (length(features_gene) > 0 && !is.null(gene_col)) do_feature_plots(features_gene, feat_index = gene_col, subdir = "gene")
} else if (lvl == "gene") {
  features_gene <- if (length(feature_list_gene) > 0) trimws(feature_list_gene[nzchar(trimws(feature_list_gene))]) else top_n_by_col(se_for_qc, gene_col, top_n_gene)
  if (length(features_gene) > 0 && !is.null(gene_col)) do_feature_plots(features_gene, feat_index = gene_col, subdir = "gene")
} else if (lvl == "peptide") {
  features_pep <- if (length(feature_list_peptide) > 0) trimws(feature_list_peptide[nzchar(trimws(feature_list_peptide))]) else if (top_n_peptide > 0) top_n_by_rownames(se_for_qc, top_n_peptide) else character(0)
  if (length(features_pep) > 0) do_feature_plots(features_pep, feat_index = NULL, subdir = "peptide")
} else if (lvl == "site") {
  features_site <- if (length(feature_list_site) > 0) trimws(feature_list_site[nzchar(trimws(feature_list_site))]) else if (top_n_site > 0) top_n_by_rownames(se_for_qc, top_n_site) else character(0)
  if (length(features_site) > 0) do_feature_plots(features_site, feat_index = NULL, subdir = "site")
}

# Phase 8: Differential expression (if >= 2 conditions, imputed data required)
if (n_conditions >= 2) {
  se_for_de <- if (!is.null(imputed_se)) imputed_se else normalized_se
  if (any(is.na(assay(se_for_de)))) {
    cat("Skipping DE: imputed data required but imputation was none or failed.\n")
  } else {
    if (!"label" %in% colnames(colData(se_for_de)) && "sample_name" %in% colnames(colData(se_for_de))) {
      colData(se_for_de)$label <- colData(se_for_de)$sample_name
    }
    de_fdr_norm <- tolower(trimws(as.character(de_fdr)))
    use_fdrtool <- de_fdr_norm %in% c("local and tail area-based", "fdrtool", "local_tail")
    tryCatch({
      if (use_fdrtool) {
        cat("Running DE (limma + fdrtool, Local and tail area-based FDR)...\n")
        de_result <- test_diff_customized(se_for_de, type = "all")
      } else {
        cat("Running DE (limma + Benjamini-Hochberg FDR)...\n")
        de_result <- test_limma_customized(se_for_de, type = "all")
      }
      de_result <- add_rejections_customized(de_result, alpha = de_alpha, lfc = de_lfc)
      de_df <- get_de_results_extended(de_result)
      de_dir <- file.path(output_dir, "de")
      dir.create(de_dir, showWarnings = FALSE, recursive = TRUE)
      write.csv(de_df, file.path(de_dir, "DE_results.csv"), row.names = FALSE)
      cat("DE_results.csv saved (extended: CI.L, CI.R, diff, p.val, p.adj, significant per contrast)\n")
      volcano_dir <- file.path(de_dir, "volcano")
      dir.create(volcano_dir, showWarnings = FALSE, recursive = TRUE)
      contrast_cols <- grep("_diff$", colnames(rowData(de_result)), value = TRUE)
      contrasts <- gsub("_diff$", "", contrast_cols)
      lvl <- metadata(de_result)$level
      for (cntrst in contrasts) {
        tryCatch({
          if (lvl %in% c("peptide", "site")) {
            p_v <- plot_peptide_volcano(de_result, cntrst, peptides = volcano_highlight_feature,
              show_other_peptides = volcano_show_other_peptides, show_gene = volcano_show_gene,
              add_names = volcano_display_names, alpha = de_alpha, lfc = de_lfc, adjusted = TRUE)
          } else {
            p_v <- plot_volcano_customized(de_result, cntrst, plot = TRUE, alpha = de_alpha, lfc = de_lfc,
              add_names = volcano_display_names, adjusted = TRUE, show_gene = volcano_show_gene,
              selected = if (length(volcano_highlight_feature) > 0) volcano_highlight_feature else NULL)
          }
          safe_name <- gsub("[^A-Za-z0-9_-]", "_", cntrst)
          ggplot2::ggsave(file.path(volcano_dir, paste0("volcano_", safe_name, ".pdf")), p_v, width = 8, height = 6)
          ggplot2::ggsave(file.path(volcano_dir, paste0("volcano_", safe_name, ".png")), p_v, width = 8, height = 6, dpi = 150)
          cat("Volcano saved:", cntrst, "\n")
        }, error = function(e) warning("Volcano ", cntrst, " failed: ", conditionMessage(e)))
      }
      tryCatch({
        hm_res <- get_cluster_heatmap_customized(de_result, type = "centered", alpha = de_alpha, lfc = de_lfc, indicate = "condition")
        if (!is.null(hm_res)) {
          pdf(file.path(de_dir, "de_heatmap.pdf"), width = 8, height = 8)
          get_cluster_heatmap_customized(de_result, type = "centered", alpha = de_alpha, lfc = de_lfc, indicate = "condition")
          dev.off()
          png(file.path(de_dir, "de_heatmap.png"), width = 8, height = 8, units = "in", res = 150)
          get_cluster_heatmap_customized(de_result, type = "centered", alpha = de_alpha, lfc = de_lfc, indicate = "condition")
          dev.off()
          cat("DE heatmap saved\n")
        }
      }, error = function(e) warning("DE heatmap failed: ", conditionMessage(e)))

      # ---------- Pathway and GO enrichment (from FragPipeAnalystR) ----------
      # pathway_database and go_database are vectors (comma-separated input).
      # Any Enrichr libraryName from maayanlab.cloud/Enrichr/datasetStatistics is supported.
      enr_dir <- file.path(output_dir, "enrichment")
      pathway_dirs <- if (tolower(pathway_direction) == "both") c("UP", "DOWN") else if (tolower(pathway_direction) == "up") "UP" else if (tolower(pathway_direction) == "down") "DOWN" else c("UP", "DOWN")
      go_dirs <- if (tolower(go_direction) == "both") c("UP", "DOWN") else if (tolower(go_direction) == "up") "UP" else if (tolower(go_direction) == "down") "DOWN" else c("UP", "DOWN")
      if (length(pathway_database) > 0) {
        dir.create(enr_dir, showWarnings = FALSE, recursive = TRUE)
        for (db in pathway_database) {
          for (dir in pathway_dirs) {
            tryCatch({
              or_result <- or_test(de_result, database = db, direction = dir, alpha = de_alpha, log2_threshold = de_lfc)
              if (!is.null(or_result) && nrow(or_result) > 0) {
                safe_name <- paste0(gsub("[^A-Za-z0-9_-]", "_", db), "_", dir)
                write.csv(or_result, file.path(enr_dir, paste0("pathway_", safe_name, ".csv")), row.names = FALSE)
                p_or <- plot_or(or_result, number = 15, alpha = de_alpha)
                ggplot2::ggsave(file.path(enr_dir, paste0("pathway_", safe_name, ".pdf")), p_or, width = 10, height = 6)
                ggplot2::ggsave(file.path(enr_dir, paste0("pathway_", safe_name, ".png")), p_or, width = 10, height = 6, dpi = 150)
                cat("Pathway enrichment saved:", db, dir, "\n")
              } else {
                cat("No pathway enrichment found for", db, dir, "\n")
              }
            }, error = function(e) warning("Pathway enrichment failed (", db, " ", dir, "): ", conditionMessage(e)))
          }
        }
      }
      if (length(go_database) > 0) {
        dir.create(enr_dir, showWarnings = FALSE, recursive = TRUE)
        for (db in go_database) {
          for (dir in go_dirs) {
            tryCatch({
              or_result <- or_test(de_result, database = db, direction = dir, alpha = de_alpha, log2_threshold = de_lfc)
              if (!is.null(or_result) && nrow(or_result) > 0) {
                safe_name <- paste0(gsub("[^A-Za-z0-9_-]", "_", db), "_", dir)
                write.csv(or_result, file.path(enr_dir, paste0("go_", safe_name, ".csv")), row.names = FALSE)
                p_or <- plot_or(or_result, number = 15, alpha = de_alpha)
                ggplot2::ggsave(file.path(enr_dir, paste0("go_", safe_name, ".pdf")), p_or, width = 10, height = 6)
                ggplot2::ggsave(file.path(enr_dir, paste0("go_", safe_name, ".png")), p_or, width = 10, height = 6, dpi = 150)
                cat("GO enrichment saved:", db, dir, "\n")
              } else {
                cat("No GO enrichment found for", db, dir, "\n")
              }
            }, error = function(e) warning("GO enrichment failed (", db, " ", dir, "): ", conditionMessage(e)))
          }
        }
      }

      # ---------- Report PDF ----------
      report_rmd_map <- c(
        "LFQ_protein" = "LFQ_report.Rmd",
        "LFQ_peptide" = "LFQ-peptide_report.Rmd",
        "LFQ_site" = "LFQ-site_report.Rmd",
        "TMT_protein" = "TMT_report.Rmd",
        "TMT_gene" = "TMT_report.Rmd",
        "TMT_peptide" = "TMT-peptide_report.Rmd",
        "TMT_site" = "TMT-site_report.Rmd"
      )
      report_rmd_key <- paste(mode, level, sep = "_")
      rmd_name <- report_rmd_map[report_rmd_key]
      if (!is.na(rmd_name)) {
        script_dir_abs <- if (dir.exists(script_dir)) normalizePath(script_dir, mustWork = FALSE) else getwd()
        rmd_path <- file.path(script_dir_abs, "reports", rmd_name)
        if (file.exists(rmd_path) && requireNamespace("rmarkdown", quietly = TRUE)) {
          tryCatch({
            de_df_rd <- as.data.frame(rowData(de_result))
            diff_cols <- grep("_diff$", colnames(de_df_rd), value = TRUE)
            valid_cntrsts <- gsub("_diff$", "", diff_cols)
            num_signif <- 0L
            if ("significant" %in% colnames(de_df_rd)) {
              num_signif <- sum(replace_na(de_df_rd$significant, FALSE), na.rm = TRUE)
            } else {
              sig_cols <- grep("_significant$", colnames(de_df_rd), value = TRUE)
              for (sc in sig_cols) num_signif <- num_signif + sum(de_df_rd[[sc]], na.rm = TRUE)
            }
            imputation_display <- if (imputation_type == "Perseus-type") "Perseus-type" else imputation_type
            report_params <- list(
              output_dir = normalizePath(output_dir, mustWork = TRUE),
              data = data_se,
              dep = function() de_result,
              alpha = de_alpha,
              lfc = de_lfc,
              normalization = normalization_method,
              imputation = imputation_display,
              fdr_correction = de_fdr,
              num_signif = num_signif,
              tested_contrasts = paste(valid_cntrsts, collapse = ", "),
              numbers_input = function() NULL,
              coverage_input = function() NULL,
              pca_input = function() NULL,
              correlation_input = function() NULL,
              missval_input = function() NULL,
              detect_input = function() NULL,
              density_input = function() NULL,
              p_hist_input = function() NULL,
              heatmap_input = function() NULL,
              cvs_input = function() NULL,
              volcano_input = function() NULL
            )
            owd <- getwd()
            on.exit(setwd(owd), add = TRUE)
            setwd(output_dir)
            rmarkdown::render(rmd_path, output_file = "report.pdf", output_dir = ".",
              params = report_params, envir = new.env(parent = parent.frame()), quiet = TRUE)
            cat("Report saved:", file.path(output_dir, "report.pdf"), "\n")
          }, error = function(e) warning("Report generation failed: ", conditionMessage(e)))
        }
      }
    }, error = function(e) {
      warning("DE failed: ", conditionMessage(e))
    })
  }
}

cat("Done.\n")
