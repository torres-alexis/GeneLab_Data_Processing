#!/usr/bin/Rscript
# FragPipe-Analyst downstream analysis: QC plots, DE (limma), feature plots, enrichment.
# Sources fp_analyst_de_fragpipe.R (filters and plot functions from MonashProteomics/FragPipe-Analyst).
print("STARTING ANALYSIS")
args = commandArgs(trailingOnly=TRUE)

# Load required libraries
library(optparse)
library(FragPipeAnalystR)
library(ggplot2)
library(fdrtool)
library(tidyr)
library(purrr)
# Source FragPipe-Analyst DE implementation (limma + fdrtool)
cmd_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", cmd_args, value = TRUE)
script_dir <- if (length(file_arg) > 0) dirname(sub("^--file=", "", file_arg)) else "."
source(file.path(script_dir, "fp_analyst_de_fragpipe.R"))
print("Libraries loaded")

# Parse command-line arguments
option_list = list(
    make_option(c("--experiment_annotation"), type="character", default=NULL, 
                help="Path to experiment annotation TSV file", metavar="FILE"),
    make_option(c("--quantification_file"), type="character", default=NULL,
                help="Path to quantification file", metavar="FILE"),
    make_option(c("--mode"), type="character", default=NULL,
                help="Quantification mode: LFQ, TMT, or DIA", metavar="MODE"),
    make_option(c("--level"), type="character", default=NULL,
                help="Analysis level: protein, peptide, gene, or site", metavar="LEVEL"),
    make_option(c("--feature_list_protein"), type="character", default="",
                help="Comma-separated protein IDs for feature plots. If empty, use top_n_protein.", metavar="STRING"),
    make_option(c("--feature_list_gene"), type="character", default="",
                help="Comma-separated gene names for feature plots. If empty, use top_n_gene.", metavar="STRING"),
    make_option(c("--top_n_protein"), type="integer", default=10,
                help="When feature_list_protein empty: plot top N variable by protein ID. 0 = skip.", metavar="INTEGER"),
    make_option(c("--top_n_gene"), type="integer", default=10,
                help="When feature_list_gene empty: plot top N variable by gene. 0 = skip.", metavar="INTEGER"),
    make_option(c("--pathway_database"), type="character", default="",
                help="Pathway database: Hallmark, KEGG, Reactome. Empty = skip.", metavar="STRING"),
    make_option(c("--pathway_direction"), type="character", default="Both",
                help="Pathway enrichment direction: Up, Down, or Both", metavar="STRING"),
    make_option(c("--go_database"), type="character", default="",
                help="Gene Ontology database: GO Biological Process, GO Cellular Component, GO Molecular Function. Empty = skip.", metavar="STRING"),
    make_option(c("--go_direction"), type="character", default="Both",
                help="GO enrichment direction: Up, Down, or Both", metavar="STRING"),
    make_option(c("--output_dir"), type="character", default="output/",
                help="Output directory", metavar="DIR"),
    make_option(c("--lfq_type"), type="character", default="Intensity",
                help="LFQ column type: Intensity or MaxLFQ (for LFQ mode only)", metavar="STRING"),
    make_option(c("--de_alpha"), type="numeric", default=0.05,
                help="Adjusted p-value threshold for DE significance", metavar="NUMERIC"),
    make_option(c("--de_lfc"), type="numeric", default=1.0,
                help="Log2 fold change threshold for DE significance", metavar="NUMERIC"),
    make_option(c("--de_fdr"), type="character", default="Benjamini Hochberg",
                help="Type of FDR correction: 'Benjamini Hochberg' or 'Local and tail area-based'", metavar="STRING"),
    make_option(c("--imputation_type"), type="character", default="man",
                help="Imputation type: none, man (Perseus-type), knn, MLE, min, zero", metavar="STRING"),
    make_option(c("--min_global_appearance"), type="numeric", default=0,
                help="At least X%% non-missing across all samples (0-100). 0 = no filter.", metavar="NUMERIC"),
    make_option(c("--min_appearance_one_condition"), type="numeric", default=0,
                help="At least X%% non-missing in at least one condition (0-100). 0 = no filter. On Ubiquitin LFQ: 10-33%% -> 2777 'reproducibly quantified'. Use --explore_filter to trace.", metavar="NUMERIC"),
    make_option(c("--explore_filter"), action="store_true", default=FALSE,
                help="Print filter value -> protein count table and exit (no DE). Use to trace source of 'reproducibly quantified' count."),
    make_option(c("--qc_show_imputed"), type="character", default="true",
                help="Use imputed data for QC plots when qc_include_both=false. true/false.", metavar="STRING"),
    make_option(c("--qc_include_both"), type="character", default="false",
                help="Generate both imputed and unimputed versions of QC plots (PCA, correlation, CVs, feature). true/false.", metavar="STRING"),
    make_option(c("--shiny_reference_de"), type="character", default=NULL,
                help="Path to DE_results.csv; if set, filter SE to same protein set before DE.", metavar="FILE"),
    make_option(c("--volcano_display_names"), type="character", default="true",
                help="Volcano plot: display names on significant points. true/false.", metavar="STRING"),
    make_option(c("--volcano_show_gene"), type="character", default="true",
                help="Volcano plot: show gene names (true) or protein/peptide ID (false).", metavar="STRING"),
    make_option(c("--volcano_highlight_feature"), type="character", default="",
                help="Comma-delimited list of features to highlight on volcano (e.g. Gene names or IDs).", metavar="STRING")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

# Assign arguments to variables
experiment_annotation = opt$experiment_annotation
quantification_file = opt$quantification_file
mode = opt$mode
level = opt$level
feature_list_protein_str = opt$feature_list_protein
feature_list_protein = if(feature_list_protein_str != "") unlist(strsplit(feature_list_protein_str, ",")) else character(0)
feature_list_gene_str = opt$feature_list_gene
feature_list_gene = if(feature_list_gene_str != "") unlist(strsplit(feature_list_gene_str, ",")) else character(0)
top_n_protein = max(0L, as.integer(opt$top_n_protein))
top_n_gene = max(0L, as.integer(opt$top_n_gene))
pathway_database = if (is.null(opt$pathway_database) || is.na(opt$pathway_database)) "" else opt$pathway_database
pathway_direction = if (is.null(opt$pathway_direction) || is.na(opt$pathway_direction)) "Both" else opt$pathway_direction
go_database = if (is.null(opt$go_database) || is.na(opt$go_database)) "" else opt$go_database
go_direction = if (is.null(opt$go_direction) || is.na(opt$go_direction)) "Both" else opt$go_direction
output_dir = opt$output_dir
lfq_type = opt$lfq_type
de_alpha = opt$de_alpha
de_lfc = opt$de_lfc
# Normalize FDR: accept "Benjamini Hochberg","BH","bh" -> BH; "Local and tail area-based","fdrtool" -> Local and tail area-based
de_fdr_raw = trimws(if (is.null(opt$de_fdr) || is.na(opt$de_fdr) || opt$de_fdr == "") "Benjamini Hochberg" else opt$de_fdr)
de_fdr = if (tolower(de_fdr_raw) %in% c("benjamini hochberg", "bh", "benjamini-hochberg")) "Benjamini Hochberg" else if (tolower(de_fdr_raw) %in% c("local and tail area-based", "fdrtool", "local_tail")) "Local and tail area-based" else de_fdr_raw
# Imputation: none, man (Perseus-type), knn, MLE, min, zero
imputation_type_raw = trimws(if (is.null(opt$imputation_type) || is.na(opt$imputation_type) || opt$imputation_type == "") "man" else opt$imputation_type)
imputation_type = if (tolower(imputation_type_raw) %in% c("none", "no")) "none" else if (tolower(imputation_type_raw) %in% c("man", "perseus", "perseus-type", "perseus_type")) "man" else if (tolower(imputation_type_raw) %in% c("knn", "min", "zero")) tolower(imputation_type_raw) else if (tolower(imputation_type_raw) == "mle") "MLE" else imputation_type_raw
min_global_appearance = opt$min_global_appearance
min_appearance_one_condition = opt$min_appearance_one_condition
explore_filter = isTRUE(opt$explore_filter)
qc_show_imputed = tolower(trimws(opt$qc_show_imputed)) %in% c("true", "1", "yes")
qc_include_both = tolower(trimws(opt$qc_include_both)) %in% c("true", "1", "yes")
shiny_reference_de = opt$shiny_reference_de
volcano_display_names <- tolower(trimws(opt$volcano_display_names)) %in% c("true", "1", "yes")
volcano_show_gene <- tolower(trimws(opt$volcano_show_gene)) %in% c("true", "1", "yes")
volcano_highlight_feature <- if (is.null(opt$volcano_highlight_feature) || trimws(opt$volcano_highlight_feature) == "") character(0) else
  trimws(strsplit(trimws(opt$volcano_highlight_feature), "\\s*,\\s*")[[1]])

print(paste("experiment_annotation:", experiment_annotation))
print(paste("quantification_file:", quantification_file))
print(paste("mode:", mode))
print(paste("level:", level))
print(paste("pathway_database:", pathway_database))
print(paste("pathway_direction:", pathway_direction))
print(paste("go_database:", go_database))
print(paste("go_direction:", go_direction))
print(paste("output_dir:", output_dir))

# Ensure output directory and subdirs exist
dir_out <- c(output_dir,
  file.path(output_dir, "enrichment"),
  file.path(output_dir, "feature", "protein"),
  file.path(output_dir, "feature", "gene"),
  file.path(output_dir, "venndiagram"),
  file.path(output_dir, "volcano"))
for (d in dir_out) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)
}

# Write variables to text file
var_file = file.path(output_dir, "fp_analyst_parameters.txt")
writeLines(c(
    "FragPipe-Analyst Parameters",
    "===========================",
    paste("experiment_annotation:", experiment_annotation),
    paste("quantification_file:", quantification_file),
    paste("mode:", mode),
    paste("level:", level),
    paste("lfq_type:", lfq_type),
    paste("de_alpha:", de_alpha),
    paste("de_lfc:", de_lfc),
    paste("de_fdr:", de_fdr),
    paste("imputation_type:", imputation_type),
    paste("feature_list_protein:", paste(feature_list_protein, collapse=", ")),
    paste("feature_list_gene:", paste(feature_list_gene, collapse=", ")),
    paste("top_n_protein:", top_n_protein),
    paste("top_n_gene:", top_n_gene),
    paste("pathway_database:", pathway_database),
    paste("pathway_direction:", pathway_direction),
    paste("go_database:", go_database),
    paste("go_direction:", go_direction),
    paste("min_global_appearance:", min_global_appearance),
    paste("min_appearance_one_condition:", min_appearance_one_condition),
    paste("qc_show_imputed:", qc_show_imputed),
    paste("qc_include_both:", qc_include_both),
    paste("volcano_display_names:", volcano_display_names),
    paste("volcano_show_gene:", volcano_show_gene),
    paste("volcano_highlight_feature:", paste(volcano_highlight_feature, collapse = ", ")),
    paste("output_dir:", output_dir)
), var_file)
print(paste("Parameters written to:", var_file))


# Create SummarizedExperiment (make_unique, make_se_customized, manual_impute, test_limma_customized, add_rejections_customized).
# From MonashProteomics/FragPipe-Analyst server.R (processed_data, make_se_customized, etc.).
print("Creating SummarizedExperiment object...")
make.unique.2 <- function(x, sep = ".") {
  ave(x, x, FUN = function(a) {
    if (length(a) > 1) paste(a, seq_along(a), sep = sep) else a
  })
}
if (mode == "LFQ" && level == "protein") {
  # LFQ protein
  temp_data <- read.table(quantification_file, header = TRUE, fill = TRUE, sep = "\t",
    quote = "", comment.char = "", blank.lines.skip = FALSE, check.names = FALSE)
  colnames(temp_data) <- make.unique.2(colnames(temp_data), "_")
  colnames(temp_data) <- gsub("-", ".", colnames(temp_data))
  temp_data <- temp_data[!grepl("contam", temp_data$Protein), ]
  data_unique <- FragPipeAnalystR::make_unique(temp_data, "Gene", "Protein ID")
  lfq_cols <- if (lfq_type == "MaxLFQ") grep("MaxLFQ", colnames(data_unique)) else
    if (lfq_type == "Spectral Count") {
      setdiff(grep("Spectral", colnames(data_unique)), grep("Total Spectral Count", colnames(data_unique)))
    } else {
      setdiff(grep("Intensity", colnames(data_unique)), grep("MaxLFQ", colnames(data_unique)))
    }
  if (length(lfq_cols) == 0) stop("No ", lfq_type, " columns found.")
  temp_df <- read.table(experiment_annotation, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  colnames(temp_df) <- tolower(colnames(temp_df))
  temp_df$condition <- make.names(temp_df$condition)
  temp_df$sample <- gsub("-", ".", temp_df$sample)
  temp_df$label <- temp_df$sample
  temp_df$label <- paste(temp_df$label,
    if (lfq_type == "MaxLFQ") "MaxLFQ.Intensity" else if (lfq_type == "Spectral Count") "Spectral.Count" else "Intensity",
    sep = " ")
  # sample_name must match assay colnames (label) for FragPipeAnalystR plot_pca, plot_cvs, plot_feature_numbers
  temp_df$sample_name <- temp_df$label
  data_se <- FragPipeAnalystR::make_se_customized(data_unique, lfq_cols, temp_df,
    log2transform = (lfq_type != "Spectral Count"), exp = "LFQ", lfq_type = lfq_type, level = "protein")
} else if (mode == "LFQ" && level == "peptide") {
  # LFQ peptide
  temp_data <- read.table(quantification_file, header = TRUE, fill = TRUE, sep = "\t",
    quote = "", comment.char = "", blank.lines.skip = FALSE, check.names = FALSE)
  colnames(temp_data) <- make.unique.2(colnames(temp_data), "_")
  colnames(temp_data) <- gsub("-", ".", colnames(temp_data))
  colnames(temp_data)[colnames(temp_data) == "Protein Description"] <- "Description"
  temp_data <- temp_data[!grepl("contam", temp_data$Protein), ]
  if (!"Modified Sequence" %in% colnames(temp_data)) {
    temp_data$Index <- paste0(temp_data$`Protein ID`, "_", temp_data$`Peptide Sequence`)
  } else {
    temp_data$Index <- paste0(temp_data$`Protein ID`, "_", temp_data$`Modified Sequence`)
  }
  data_unique <- FragPipeAnalystR::make_unique(temp_data, "Protein ID", "Index")
  lfq_cols <- if (lfq_type == "MaxLFQ") grep("MaxLFQ", colnames(data_unique)) else
    if (lfq_type == "Spectral Count") {
      setdiff(grep("Spectral", colnames(data_unique)), grep("Total Spectral Count", colnames(data_unique)))
    } else {
      setdiff(grep("Intensity", colnames(data_unique)), grep("MaxLFQ", colnames(data_unique)))
    }
  if (length(lfq_cols) == 0) stop("No ", lfq_type, " columns found in peptide table.")
  temp_df <- read.table(experiment_annotation, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  colnames(temp_df) <- tolower(colnames(temp_df))
  temp_df$condition <- make.names(temp_df$condition)
  temp_df$sample <- gsub("-", ".", temp_df$sample)
  temp_df$label <- temp_df$sample
  temp_df$label <- paste(temp_df$label,
    if (lfq_type == "MaxLFQ") "MaxLFQ.Intensity" else if (lfq_type == "Spectral Count") "Spectral.Count" else "Intensity",
    sep = " ")
  # sample_name must match assay colnames (label) for FragPipeAnalystR plot_pca, plot_cvs, plot_feature_numbers
  temp_df$sample_name <- temp_df$label
  data_se <- FragPipeAnalystR::make_se_customized(data_unique, lfq_cols, temp_df,
    log2transform = (lfq_type != "Spectral Count"), exp = "LFQ", lfq_type = lfq_type, level = "peptide")
} else {
  data_se <- make_se_from_files(
    quant_table_path = quantification_file,
    exp_anno_path = experiment_annotation,
    type = mode,
    level = level,
    lfq_type = lfq_type
  )
}

if (is.null(data_se)) {
    stop("Failed to create SummarizedExperiment object. Check input files and parameters.")
}

print(paste("Successfully created SummarizedExperiment object with", nrow(data_se), "features and", ncol(data_se), "samples"))

# ========== Filter exploration mode (trace source of "reproducibly quantified" count) ==========
# Report text like "2777 proteins were reproducibly quantified, 183 differ significantly" comes from
# FragPipe-Analyst's filter_by_condition (Min % non-missing in at least one condition).
# On Ubiquitin LFQ (4 conditions, 3 reps each): 10-33% -> 2777 proteins; 35%+ -> 2606.
if (explore_filter) {
  se <- data_se
  if (min_global_appearance > 0 && min_global_appearance <= 100) {
    pct_missing_allowed <- (100 - min_global_appearance) / 100
    ridx <- rowSums(is.na(assay(se))) / ncol(assay(se)) <= pct_missing_allowed
    se <- se[ridx, ]
  }
  cat("\n=== filter_by_condition exploration ===\n")
  cat("Source of 'reproducibly quantified' count in FragPipe-Analyst reports.\n")
  cat("Data: ", nrow(se), " proteins", if (min_global_appearance > 0) paste0(" (after global_filter ", min_global_appearance, "%)") else "", "\n\n")
  conditions <- unique(colData(se)$condition)
  for (pct in c(0, 10, 20, 25, 30, 33, 35, 40, 45, 50, 60, 70, 80, 90, 100)) {
    min_frac <- pct / 100
    row_ok <- rep(FALSE, nrow(se))
    for (c in conditions) {
      se_c <- se[, colData(se)$condition == c]
      ridx <- rowSums(!is.na(assay(se_c))) / ncol(assay(se_c)) >= min_frac
      row_ok <- row_ok | ridx
    }
    n <- sum(row_ok)
    marker <- if (n == 2777) "  <- matches Ubiquitin report" else ""
    cat(sprintf("  %3d%%: %d proteins%s\n", pct, n, marker))
  }
  cat("\n")
  quit(save = "no", status = 0)
}

# ========== Filtering (FragPipe-Analyst R/filter.R) ==========
# server.R: global_filter(processed_data(), 100 - min_global), filter_by_condition(..., min_appearance_each_condition)
processed_se <- data_se  # keep for Absence/Presence
filtered_se <- data_se
if (min_global_appearance > 0 && min_global_appearance <= 100) {
  filtered_se <- global_filter(filtered_se, 100 - min_global_appearance)
  print(paste("global_filter: kept", nrow(filtered_se), "features (min", min_global_appearance, "% non-missing globally)"))
}
if (min_appearance_one_condition > 0 && min_appearance_one_condition <= 100) {
  filtered_se <- filter_by_condition(filtered_se, min_appearance_one_condition)
  print(paste("filter_by_condition: kept", nrow(filtered_se), "features (min", min_appearance_one_condition, "% in at least one condition)"))
}
data_se <- filtered_se

# Optional: filter to same protein set as external DE_results.csv before imputation/DE (for variance moderation consistency)
if (!is.null(shiny_reference_de) && file.exists(shiny_reference_de)) {
  ref_df <- read.table(shiny_reference_de, header = TRUE, sep = ",", quote = "\"",
    stringsAsFactors = FALSE, check.names = FALSE)
  id_col <- if ("Protein ID" %in% colnames(ref_df)) "Protein ID" else
    if ("Protein.ID" %in% colnames(ref_df)) "Protein.ID" else "ID"
  ref_ids <- unique(ref_df[[id_col]])
  ref_ids <- ref_ids[nzchar(trimws(ref_ids))]
  our_ids <- rownames(data_se)
  # Match by ID (rowData$ID or rownames)
  if ("ID" %in% colnames(rowData(data_se))) {
    match_idx <- rowData(data_se)$ID %in% ref_ids
  } else {
    match_idx <- rownames(data_se) %in% ref_ids
  }
  n_before <- nrow(data_se)
  data_se <- data_se[match_idx, ]
  print(paste("shiny_reference_de: filtered to", nrow(data_se), "proteins (from", n_before, ")"))
}

# ========== PCA ==========
# PCA uses complete cases; for LFQ/DIA with many missing values, impute first
n_conditions <- length(unique(colData(data_se)$condition))
has_missing <- any(is.na(assay(data_se)))

# For LFQ/DIA with missing values: impute before PCA and DE
imputed_se <- NULL
if (has_missing && mode %in% c("LFQ", "DIA") && imputation_type != "none") {
    print(paste("Imputing missing values (", imputation_type, ", seed=123)...", sep = ""))
    imputed_se <- impute_customized(data_se, fun = imputation_type)
}

# QC plots: when qc_include_both=TRUE, generate both imputed and unimputed versions
# Otherwise use qc_show_imputed to pick one.
se_imputed <- if (!is.null(imputed_se)) imputed_se else data_se
se_unimputed <- data_se
qc_versions <- if (qc_include_both && !is.null(imputed_se)) {
  print("QC plots: generating both imputed and unimputed versions")
  list(
    list(se = se_imputed, suffix = "_imputed", desc = "imputed"),
    list(se = se_unimputed, suffix = "_unimputed", desc = "unimputed")
  )
} else {
  se_for_qc <- if (qc_show_imputed && !is.null(imputed_se)) se_imputed else se_unimputed
  print(paste("QC plots using:", if (qc_show_imputed && !is.null(imputed_se)) "imputed" else "unimputed", "data"))
  list(list(se = se_for_qc, suffix = "", desc = ""))
}

print("Plotting PCA...")
for (vv in qc_versions) {
  # Condition only (no replicate shapes)
  p_pca <- plot_pca(vv$se, indicate = "condition", plot = TRUE)
  base <- paste0("pca", vv$suffix)
  ggplot2::ggsave(file.path(output_dir, paste0(base, ".pdf")), p_pca, width = 8, height = 6)
  ggplot2::ggsave(file.path(output_dir, paste0(base, ".png")), p_pca, width = 8, height = 6, dpi = 150)
  print(paste("PCA saved:", base, vv$desc))
}

# ========== Correlation heatmap (FragPipe-Analyst plot_cor_customized) ==========
print("Plotting correlation heatmap...")
for (vv in qc_versions) {
  tryCatch({
    n_samp <- ncol(vv$se)
    ht_corr <- plot_cor_customized(vv$se, significant = FALSE, indicate = "condition", plot = FALSE,
      font_size = if (n_samp > 10) 9 else 12)
    base <- paste0("correlation_heatmap", vv$suffix)
    # Scale figure size so sample labels don't truncate in PDF
    fig_side <- min(14, max(7, 5 + n_samp * 0.35))
    pdf(file.path(output_dir, paste0(base, ".pdf")), width = fig_side, height = fig_side)
    ComplexHeatmap::draw(ht_corr, heatmap_legend_side = "top")
    dev.off()
    png(file.path(output_dir, paste0(base, ".png")), width = fig_side, height = fig_side, units = "in", res = 150)
    ComplexHeatmap::draw(ht_corr, heatmap_legend_side = "top")
    dev.off()
    print(paste("Correlation heatmap saved:", base, vv$desc))
  }, error = function(e) {
    warning("Correlation heatmap failed: ", conditionMessage(e))
  })
}

# ========== Missing value heatmap (only if data has NAs) ==========
if (has_missing) {
  print("Plotting missing value heatmap...")
  tryCatch({
    pdf(file.path(output_dir, "missing_value_heatmap.pdf"), width = 8, height = 6)
    plot_missval_customized(data_se)
    dev.off()
    png(file.path(output_dir, "missing_value_heatmap.png"), width = 8, height = 6, units = "in", res = 150)
    plot_missval_customized(data_se)
    dev.off()
    print(paste("Missing value heatmap saved to", output_dir))
  }, error = function(e) {
    warning("Missing value heatmap failed: ", conditionMessage(e))
  })
}

# ========== Feature numbers per sample (FragPipeAnalystR::plot_feature_numbers) ==========
print("Plotting feature numbers...")
tryCatch({
  p_fn <- FragPipeAnalystR::plot_feature_numbers(data_se, fill = "condition")
  ggplot2::ggsave(file.path(output_dir, "feature_numbers.pdf"), p_fn, width = 8, height = 5)
  ggplot2::ggsave(file.path(output_dir, "feature_numbers.png"), p_fn, width = 8, height = 5, dpi = 150)
  print(paste("Feature numbers plot saved to", output_dir))
}, error = function(e) {
  warning("Feature numbers plot failed: ", conditionMessage(e))
})

# Sample coverage (from FragPipe-Analyst plot_coverage_customized): barplot of features in 1, 2, ..., N samples.
print("Plotting sample coverage...")
tryCatch({
  p_cov <- plot_coverage_customized(data_se, plot = TRUE)
  ggplot2::ggsave(file.path(output_dir, "sample_coverage.pdf"), p_cov, width = 8, height = 5)
  ggplot2::ggsave(file.path(output_dir, "sample_coverage.png"), p_cov, width = 8, height = 5, dpi = 150)
  print(paste("Sample coverage plot saved to", output_dir))
}, error = function(e) {
  warning("Sample coverage plot failed: ", conditionMessage(e))
})

# ========== Density plot (QC) - Monash plot_density: original, filtered, imputed ==========
print("Plotting density plot...")
tryCatch({
  ses_dens <- list("original data" = processed_se, "filtered data" = data_se)
  if (!is.null(imputed_se)) ses_dens[["imputed data"]] <- imputed_se
  p_dens <- plot_density_custom(ses_dens)
  ggplot2::ggsave(file.path(output_dir, "density.pdf"), p_dens, width = 8, height = 7)
  ggplot2::ggsave(file.path(output_dir, "density.png"), p_dens, width = 8, height = 7, dpi = 150)
  print("Density plot saved")
}, error = function(e) { warning("Density plot failed: ", conditionMessage(e)) })

# Absence/presence (from MonashProteomics/FragPipe-Analyst): data_attendance, Venn, UpSet, Jaccard.
if (n_conditions >= 2 && mode %in% c("LFQ", "DIA")) {
  print("Plotting Absence/Presence analyses...")
  exp <- if (!is.null(metadata(processed_se)$exp)) metadata(processed_se)$exp else mode
  level <- if (!is.null(metadata(processed_se)$level)) metadata(processed_se)$level else "protein"
  att_df <- data_attendance_custom(processed_se, exp = exp, level = level)
  conditions <- unique(colData(processed_se)$condition)
  # Jaccard - Monash plot_Jaccard (sample-level); main output folder
  tryCatch({
    pdf(file.path(output_dir, "jaccard.pdf"), width = 7, height = 6)
    plot_Jaccard_custom(processed_se, plot = TRUE, exp = exp)
    dev.off()
    png(file.path(output_dir, "jaccard.png"), width = 7, height = 6, units = "in", res = 150)
    plot_Jaccard_custom(processed_se, plot = TRUE, exp = exp)
    dev.off()
    print("Jaccard similarity plot saved")
  }, error = function(e) { warning("Jaccard plot failed: ", conditionMessage(e)) })
  # Venn - Monash ggVennDiagram (2-way per pair); venndiagram/
  vd_dir <- file.path(output_dir, "venndiagram")
  for (pr in utils::combn(conditions, 2, simplify = FALSE)) {
    tryCatch({
      v <- plot_venn_custom(att_df, pr[1], pr[2], cond3 = NULL)
      if (!is.null(v)) {
        safe_name <- paste0(gsub("[^A-Za-z0-9_-]", "_", pr[1]), "_vs_", gsub("[^A-Za-z0-9_-]", "_", pr[2]))
        ggplot2::ggsave(file.path(vd_dir, paste0("venn_", safe_name, ".pdf")), v, width = 6, height = 6)
        ggplot2::ggsave(file.path(vd_dir, paste0("venn_", safe_name, ".png")), v, width = 6, height = 6, dpi = 150)
        print(paste("Venn plot saved:", paste(pr, collapse = " vs ")))
      }
    }, error = function(e) { warning("Venn plot failed: ", conditionMessage(e)) })
  }
  # UpSet - Monash upset_plot_input; main output folder
  tryCatch({
    if (length(grep("^#Occurences_", colnames(att_df))) >= 2) {
      pdf(file.path(output_dir, "upset.pdf"), width = 10, height = 6, onefile = FALSE)
      plot_upset_custom(att_df)
      dev.off()
      png(file.path(output_dir, "upset.png"), width = 10, height = 6, units = "in", res = 150, type = "cairo")
      plot_upset_custom(att_df)
      dev.off()
      print("UpSet plot saved")
    }
  }, error = function(e) { warning("UpSet plot failed: ", conditionMessage(e)) })
}

# ========== Sample CVs ==========
print("Plotting sample CVs...")
for (vv in qc_versions) {
  tryCatch({
    p_cvs <- plot_cvs(vv$se, id = "sample_name", scale = TRUE)
    base <- paste0("sample_cvs", vv$suffix)
    ggplot2::ggsave(file.path(output_dir, paste0(base, ".pdf")), p_cvs, width = 8, height = 5)
    ggplot2::ggsave(file.path(output_dir, paste0(base, ".png")), p_cvs, width = 8, height = 5, dpi = 150)
    print(paste("Sample CVs saved:", base, vv$desc))
  }, error = function(e) {
    warning("Sample CVs plot failed: ", conditionMessage(e))
  })
}

# Feature plots (plot_feature_custom from MonashProteomics/FragPipe-Analyst R/customized.R): protein/gene list or top N.
do_feature_plots <- function(features, feat_index, subdir) {
  if (length(features) == 0) return(invisible(NULL))
  # folder_name -> type for plot_feature_custom (violinplot uses type="violin")
  plot_specs <- list(boxplot = "boxplot", violinplot = "violin")
  for (i in seq_along(features)) {
    feat <- features[i]
    if (is.na(feat) || nchar(as.character(feat)) == 0) next
    feat <- as.character(feat)
    safe_name <- gsub("[^A-Za-z0-9_-]", "_", feat)
    for (vv in qc_versions) {
      prots <- if (is.null(feat_index)) feat else rownames(vv$se)[rowData(vv$se)[[feat_index]] == feat]
      if (length(prots) == 0) next
      for (folder_name in names(plot_specs)) {
        ptype <- plot_specs[[folder_name]]
        tryCatch({
          p_f <- plot_feature_custom(vv$se, prots, type = ptype, show_gene = !is.null(feat_index))
          feat_dir <- file.path(output_dir, "feature", subdir, folder_name)
          if (!dir.exists(feat_dir)) dir.create(feat_dir, recursive = TRUE)
          base <- paste0(folder_name, "_feature_", safe_name, vv$suffix)
          ggplot2::ggsave(file.path(feat_dir, paste0(base, ".pdf")), p_f, width = 6, height = 4)
          ggplot2::ggsave(file.path(feat_dir, paste0(base, ".png")), p_f, width = 6, height = 4, dpi = 150)
          print(paste("Feature plot saved:", file.path(folder_name, base), vv$desc))
        }, error = function(e) {
          warning("Feature plot for ", feat, " ", folder_name, " ", vv$desc, " failed: ", conditionMessage(e))
        })
      }
    }
  }
}

# Protein/peptide IDs (rownames). For peptide level, subdir = "peptide"
id_subdir <- if (!is.null(metadata(qc_versions[[1]]$se)$level) && metadata(qc_versions[[1]]$se)$level == "peptide") "peptide" else "protein"
features_protein <- character(0)
if (length(feature_list_protein) > 0) {
  features_protein <- trimws(feature_list_protein[nchar(trimws(feature_list_protein)) > 0])
} else if (top_n_protein > 0) {
  se_var <- qc_versions[[1]]$se
  if (nrow(se_var) > 0) {
    vars <- apply(assay(se_var), 1, function(x) var(x, na.rm = TRUE))
    vars[is.na(vars)] <- 0
    top_idx <- order(vars, decreasing = TRUE)[seq_len(min(top_n_protein * 2, length(vars)))]
    seen <- character(0)
    for (i in top_idx) {
      id <- rownames(se_var)[i]
      if (is.na(id) || !nzchar(id) || id %in% seen) next
      seen <- c(seen, id)
      features_protein <- c(features_protein, id)
      if (length(features_protein) >= top_n_protein) break
    }
    print(paste("Top", length(features_protein), "variable features by", id_subdir, "ID for plotting"))
  }
}
if (length(features_protein) > 0) {
  print(paste("Plotting feature(s) by", id_subdir, "ID..."))
  do_feature_plots(features_protein, feat_index = NULL, subdir = id_subdir)
}

# Gene names - map to protein IDs, plot_feature_custom with show_gene=TRUE
features_gene <- character(0)
gene_col <- if (nrow(qc_versions[[1]]$se) > 0 && "Gene" %in% colnames(rowData(qc_versions[[1]]$se))) "Gene" else NULL
if (length(feature_list_gene) > 0) {
  features_gene <- trimws(feature_list_gene[nchar(trimws(feature_list_gene)) > 0])
} else if (top_n_gene > 0 && !is.null(gene_col)) {
  se_var <- qc_versions[[1]]$se
  vars <- apply(assay(se_var), 1, function(x) var(x, na.rm = TRUE))
  vars[is.na(vars)] <- 0
  n_candidates <- min(top_n_gene * 10, length(vars))
  top_idx <- order(vars, decreasing = TRUE)[seq_len(n_candidates)]
  seen <- character(0)
  for (i in top_idx) {
    id <- as.character(rowData(se_var)[[gene_col]][i])
    if (is.na(id) || !nzchar(trimws(id)) || id %in% seen) next
    seen <- c(seen, id)
    features_gene <- c(features_gene, id)
    if (length(features_gene) >= top_n_gene) break
  }
  print(paste("Top", length(features_gene), "variable features by gene for plotting"))
}
if (is.null(gene_col) && (top_n_gene > 0 || length(feature_list_gene) > 0)) {
  print("Skipping gene feature plots: no 'Gene' column in rowData")
}
if (length(features_gene) > 0 && !is.null(gene_col)) {
  print("Plotting feature(s) by gene...")
  do_feature_plots(features_gene, feat_index = gene_col, subdir = "gene")
}

# ========== Export unimputed matrix (filtered, before imputation) ==========
unimputed_df <- cbind(
    as.data.frame(rowData(data_se)),
    as.data.frame(assay(data_se))
)
write.table(unimputed_df, file.path(output_dir, "unimputed_matrix.tsv"),
    sep = "\t", row.names = FALSE, quote = FALSE)
print(paste("Unimputed matrix saved to", output_dir))

# Export normalized matrix. From MonashProteomics/FragPipe-Analyst server.R normalized_data(); when norm="none", normalized = filtered.
normalized_df <- cbind(
    as.data.frame(rowData(data_se)),
    as.data.frame(assay(data_se))
)
write.table(normalized_df, file.path(output_dir, "normalized_matrix.tsv"),
    sep = "\t", row.names = FALSE, quote = FALSE)
print(paste("Normalized matrix saved to", output_dir))

# ========== Export imputed matrix (if imputation was done) ==========
if (!is.null(imputed_se)) {
    imputed_df <- cbind(
        as.data.frame(rowData(imputed_se)),
        as.data.frame(assay(imputed_se))
    )
    write.table(imputed_df, file.path(output_dir, "imputed_matrix.tsv"),
        sep = "\t", row.names = FALSE, quote = FALSE)
    print(paste("Imputed matrix saved to", output_dir))
}

# ========== Differential expression (if >= 2 conditions) ==========
if (n_conditions >= 2) {
    se_for_de <- if (!is.null(imputed_se)) imputed_se else data_se
    if (any(is.na(assay(se_for_de)))) {
        stop("Differential expression requires imputed data. Imputation failed or was skipped.")
    }

    # Ensure colData has "label" (FragPipe-Analyst DE expects it; FragPipeAnalystR may use sample_name)
    if (!"label" %in% colnames(colData(se_for_de)) && "sample_name" %in% colnames(colData(se_for_de))) {
      colData(se_for_de)$label <- colData(se_for_de)$sample_name
    }
    use_local_tail <- (de_fdr == "Local and tail area-based")
    if (use_local_tail) {
      print("Running differential expression (limma + fdrtool, Local and tail area-based FDR)...")
      de_result <- test_diff_customized(se_for_de, type = "all")
      de_result_updated <- add_rejections_customized(de_result, alpha = de_alpha, lfc = de_lfc)
    } else {
      print("Running differential expression (limma + Benjamini Hochberg FDR)...")
      de_result <- test_limma_customized(se_for_de, type = "all")
      de_result_updated <- add_rejections_customized(de_result, alpha = de_alpha, lfc = de_lfc)
    }

    # Export DE results
    de_df <- as.data.frame(rowData(de_result_updated))
    write.table(de_df, file.path(output_dir, "de_results.tsv"),
        sep = "\t", row.names = FALSE, quote = FALSE)
    print(paste("DE results saved to", output_dir))

    # Volcano plots per contrast (Display names, Show gene names, highlight selected features)
    volcano_dir <- file.path(output_dir, "volcano")
    valid_cntrsts <- gsub("_diff$", "", grep("_diff$", colnames(de_df), value = TRUE))
    v_name_col <- if (volcano_show_gene && "Gene" %in% colnames(de_df)) "Gene" else "ID"
    highlight_vec <- volcano_highlight_feature
    for (cntrst in valid_cntrsts) {
        p_v <- plot_volcano_customized(de_result_updated, cntrst, plot = TRUE, alpha = de_alpha, lfc = de_lfc,
          name_col = v_name_col, add_names = volcano_display_names, adjusted = TRUE)
        # Add highlight layer if features specified
        if (length(highlight_vec) > 0) {
            diff_col <- paste0(cntrst, "_diff")
            pval_col <- paste0(cntrst, "_p.adj")
            if (!pval_col %in% colnames(de_df)) pval_col <- paste0(cntrst, "_p.val")
            hl_vec_upper <- toupper(highlight_vec)
            match_id <- toupper(de_df$ID) %in% hl_vec_upper
            match_gene <- if ("Gene" %in% colnames(de_df)) toupper(de_df$Gene) %in% hl_vec_upper else FALSE
            match_name <- toupper(de_df$name) %in% hl_vec_upper
            hit <- which(match_id | match_gene | match_name)
            if (length(hit) > 0) {
                hl_df <- data.frame(
                    diff = de_df[[diff_col]][hit],
                    p_values = -log10(pmax(de_df[[pval_col]][hit], 1e-300)),
                    label = de_df[[v_name_col]][hit]
                )
                hl_df <- hl_df[is.finite(hl_df$p_values) & !is.na(hl_df$diff), , drop = FALSE]
                if (nrow(hl_df) > 0) {
                    p_v <- p_v +
                        ggplot2::geom_point(data = hl_df, ggplot2::aes(x = diff, y = p_values),
                            color = "maroon", size = 3.5, alpha = 0.9, inherit.aes = FALSE) +
                        ggrepel::geom_text_repel(data = hl_df, ggplot2::aes(x = diff, y = p_values, label = label),
                            color = "maroon", size = 3, inherit.aes = FALSE,
                            box.padding = unit(0.2, "lines"), point.padding = unit(0.2, "lines"), segment.size = 0.5)
                }
            }
        }
        safe_name <- gsub("[^A-Za-z0-9_-]", "_", cntrst)
        ggplot2::ggsave(file.path(volcano_dir, paste0("volcano_", safe_name, ".pdf")), p_v, width = 8, height = 6)
        ggplot2::ggsave(file.path(volcano_dir, paste0("volcano_", safe_name, ".png")), p_v, width = 8, height = 6, dpi = 150)
        print(paste("Volcano plot saved:", cntrst))
    }

    # DE heatmap (significant proteins, sample-level)
    # Uses get_cluster_heatmap_customized for correct col labels (label->sample_name mapping fallback)
    # When no significant features it returns NULL (report shows no figure; standalone gets a placeholder)
    print("Plotting DE heatmap...")
    tryCatch({
      hm_res <- get_cluster_heatmap_customized(de_result_updated, type = "centered", alpha = de_alpha, lfc = de_lfc,
        indicate = "condition")
      pdf(file.path(output_dir, "de_heatmap.pdf"), width = 8, height = 8)
      if (is.null(hm_res)) {
        plot.new()
        text(0.5, 0.5, "No differentially expressed features available for the heatmap", cex = 1.2)
      } else {
        get_cluster_heatmap_customized(de_result_updated, type = "centered", alpha = de_alpha, lfc = de_lfc,
          indicate = "condition")
      }
      dev.off()
      png(file.path(output_dir, "de_heatmap.png"), width = 8, height = 8, units = "in", res = 150)
      if (is.null(hm_res)) {
        plot.new()
        text(0.5, 0.5, "No differentially expressed features available for the heatmap", cex = 1.2)
      } else {
        get_cluster_heatmap_customized(de_result_updated, type = "centered", alpha = de_alpha, lfc = de_lfc,
          indicate = "condition")
      }
      dev.off()
      print(paste("DE heatmap saved to", output_dir))
    }, error = function(e) {
      warning("DE heatmap failed: ", conditionMessage(e))
    })

    # ---------- Pathway and GO enrichment ----------
    enr_dir <- file.path(output_dir, "enrichment")
    pathway_dbs <- c("Hallmark", "KEGG", "Reactome")
    pathway_dirs <- if (tolower(pathway_direction) == "both") c("UP", "DOWN") else if (tolower(pathway_direction) == "up") "UP" else if (tolower(pathway_direction) == "down") "DOWN" else c("UP", "DOWN")
    if (pathway_database != "" && pathway_database %in% pathway_dbs) {
      for (dir in pathway_dirs) {
        print(paste("Running pathway enrichment:", pathway_database, "direction:", dir))
        tryCatch({
          or_result <- or_test(de_result_updated, database = pathway_database, direction = dir, alpha = de_alpha, log2_threshold = de_lfc)
          if (!is.null(or_result) && nrow(or_result) > 0) {
            safe_name <- paste0(gsub(" ", "_", pathway_database), "_", dir)
            write.table(or_result, file.path(enr_dir, paste0("enrichment_pathway_", safe_name, ".tsv")),
              sep = "\t", row.names = FALSE, quote = FALSE)
            p_or <- plot_or(or_result, number = 15, alpha = de_alpha)
            ggplot2::ggsave(file.path(enr_dir, paste0("enrichment_pathway_", safe_name, ".pdf")), p_or, width = 10, height = 6)
            ggplot2::ggsave(file.path(enr_dir, paste0("enrichment_pathway_", safe_name, ".png")), p_or, width = 10, height = 6, dpi = 150)
            print(paste("Pathway enrichment saved:", pathway_database, dir))
          } else {
            print(paste("No pathway enrichment found for", pathway_database, dir))
          }
        }, error = function(e) {
          warning("Pathway enrichment failed (", pathway_database, " ", dir, "): ", conditionMessage(e))
        })
      }
    }

    # ---------- Gene Ontology enrichment (GO BP, CC, MF) ----------
    go_dbs <- c("GO Biological Process", "GO Cellular Component", "GO Molecular Function")
    go_dirs <- if (tolower(go_direction) == "both") c("UP", "DOWN") else if (tolower(go_direction) == "up") "UP" else if (tolower(go_direction) == "down") "DOWN" else c("UP", "DOWN")
    if (go_database != "" && go_database %in% go_dbs) {
      for (dir in go_dirs) {
        print(paste("Running GO enrichment:", go_database, "direction:", dir))
        tryCatch({
          or_result <- or_test(de_result_updated, database = go_database, direction = dir, alpha = de_alpha, log2_threshold = de_lfc)
          if (!is.null(or_result) && nrow(or_result) > 0) {
            safe_name <- paste0(gsub(" ", "_", go_database), "_", dir)
            write.table(or_result, file.path(enr_dir, paste0("enrichment_go_", safe_name, ".tsv")),
              sep = "\t", row.names = FALSE, quote = FALSE)
            p_or <- plot_or(or_result, number = 15, alpha = de_alpha)
            ggplot2::ggsave(file.path(enr_dir, paste0("enrichment_go_", safe_name, ".pdf")), p_or, width = 10, height = 6)
            ggplot2::ggsave(file.path(enr_dir, paste0("enrichment_go_", safe_name, ".png")), p_or, width = 10, height = 6, dpi = 150)
            print(paste("GO enrichment saved:", go_database, dir))
          } else {
            print(paste("No GO enrichment found for", go_database, dir))
          }
        }, error = function(e) {
          warning("GO enrichment failed (", go_database, " ", dir, "): ", conditionMessage(e))
        })
      }
    }
} else {
    print(paste("Skipping DE: only", n_conditions, "condition(s). Need >= 2 for differential expression."))
}

# Report PDF: map (mode, level) to report Rmd (from MonashProteomics/FragPipe-Analyst reports/).
report_rmd_map <- c(
  "LFQ_protein" = "LFQ_report.Rmd",
  "LFQ_peptide" = "LFQ-peptide_report.Rmd",
  "TMT_protein" = "TMT_report.Rmd",
  "TMT_gene" = "TMT_report.Rmd",
  "TMT_peptide" = "TMT-peptide_report.Rmd",
  "TMT_site" = "TMT-site_report.Rmd",
  "DIA_protein" = "DIA_report.Rmd",
  "DIA_peptide" = "DIA-peptide_report.Rmd",
  "DIA_site" = "DIA-site_report.Rmd"
)
report_rmd_key <- paste(mode, level, sep = "_")
rmd_name <- report_rmd_map[report_rmd_key]
if (n_conditions >= 2 && !is.na(rmd_name)) {
  rmd_path <- file.path(script_dir, rmd_name)
  if (file.exists(rmd_path) && requireNamespace("rmarkdown", quietly = TRUE)) {
    tryCatch({
      de_for_report <- if (exists("de_result_updated")) de_result_updated else NULL
      se_for_qc <- if (!is.null(imputed_se)) imputed_se else data_se
      valid_cntrsts <- character(0)
      num_signif <- 0L
      if (!is.null(de_for_report)) {
        de_df <- as.data.frame(rowData(de_for_report))
        diff_cols <- grep("_diff$", colnames(de_df), value = TRUE)
        valid_cntrsts <- gsub("_diff$", "", diff_cols)
        if ("significant" %in% colnames(de_df)) {
          num_signif <- sum(replace_na(de_df$significant, FALSE), na.rm = TRUE)
        } else {
          sig_cols <- grep("_significant$", colnames(de_df), value = TRUE)
          for (sc in sig_cols) num_signif <- num_signif + sum(de_df[[sc]], na.rm = TRUE)
        }
      }
      pca_input <- function() plot_pca(se_for_qc, indicate = "condition", plot = TRUE)
      correlation_input <- function() {
        ht <- plot_cor_customized(se_for_qc, significant = FALSE, indicate = "condition", plot = FALSE)
        ComplexHeatmap::draw(ht, heatmap_legend_side = "top")
      }
      cvs_input <- function() plot_cvs(se_for_qc, id = "sample_name", scale = TRUE)
      numbers_input <- function() FragPipeAnalystR::plot_feature_numbers(data_se, fill = "condition")
      coverage_input <- function() plot_coverage_customized(data_se, plot = TRUE)
      missval_input <- function() {
        if (any(is.na(assay(data_se)))) plot_missval_customized(data_se)
      }
      density_input <- function() {
        ses <- list("original data" = processed_se, "filtered data" = data_se)
        if (!is.null(imputed_se)) ses[["imputed data"]] <- imputed_se
        plot_density_custom(ses)
      }
      heatmap_input <- function() {
        if (is.null(de_for_report)) return(invisible(NULL))
        get_cluster_heatmap_customized(de_for_report, type = "centered", alpha = de_alpha, lfc = de_lfc,
          indicate = "condition")
      }
      dep_fn <- function() de_for_report
      # Volcano: always draw something (plot or placeholder)
      volcano_input <- function() {
        if (!is.null(de_for_report) && length(valid_cntrsts) > 0) {
          for (i in valid_cntrsts) {
            tryCatch(
              print(plot_volcano_customized(de_for_report, contrast = i, label_size = 2, add_names = FALSE,
                alpha = de_alpha, lfc = de_lfc, adjusted = TRUE)),
              error = function(e) NULL
            )
          }
        } else {
          print(ggplot2::ggplot() + ggplot2::annotate("text", x = 0.5, y = 0.5, label = "No contrasts tested.", size = 5) + ggplot2::theme_void())
        }
        invisible()
      }
      report_params <- list(
        data = data_se,
        dep = dep_fn,
        alpha = de_alpha,
        lfc = de_lfc,
        normalization = "none",
        imputation = if (imputation_type == "man") "Perseus-type" else imputation_type,
        fdr_correction = de_fdr,
        num_signif = num_signif,
        tested_contrasts = valid_cntrsts,
        numbers_input = numbers_input,
        coverage_input = coverage_input,
        pca_input = pca_input,
        correlation_input = correlation_input,
        missval_input = missval_input,
        detect_input = function() NULL,
        density_input = density_input,
        p_hist_input = function() NULL,
        heatmap_input = heatmap_input,
        cvs_input = cvs_input,
        volcano_input = volcano_input
      )
      rmarkdown::render(rmd_path, output_file = "report.pdf", output_dir = output_dir,
        params = report_params, envir = new.env(parent = parent.frame()), quiet = TRUE)
      print(paste("Report saved:", file.path(output_dir, "report.pdf")))
    }, error = function(e) {
      warning("Report generation failed: ", conditionMessage(e))
    })
  }
}

print("FragPipe-Analyst analysis finished successfully.")

