# fp_plot_helper.R
# QC plot functions for FragPipe-Analyst. From FragPipeAnalystR, FragPipe-Analyst.

library(ggplot2)
library(dplyr)
library(tibble)
library(tidyr)
library(purrr)
library(assertthat)

# plot_pca_custom: from FragPipeAnalystR.
# complete.cases, top n variable by SD, prcomp(t(df), scale=F), color by indicate.
plot_pca_custom <- function(dep, x = 1, y = 2, indicate = "condition",
  n = 500, point_size = 8, label_size = 3, plot = TRUE, ID_col = "sample_name", scale = FALSE) {
  assertthat::assert_that(
    inherits(dep, "SummarizedExperiment"),
    is.numeric(x), length(x) == 1,
    is.numeric(y), length(y) == 1,
    is.numeric(n), length(n) == 1,
    is.character(indicate), length(indicate) <= 2,
    is.logical(plot), length(plot) == 1
  )
  if (x > ncol(dep) || y > ncol(dep)) {
    stop("'x' and/or 'y' must be <= ", ncol(dep), call. = FALSE)
  }
  columns <- colnames(colData(dep))
  if (length(indicate) > 0 && any(!indicate %in% columns)) {
    stop("'", paste(indicate, collapse = "', '"), "' not in colData. Valid: ", paste(columns, collapse = ", "), call. = FALSE)
  }

  data <- assay(dep)
  data <- data[complete.cases(data), , drop = FALSE]
  if (nrow(data) == 0) stop("No complete cases for PCA. Impute missing values first.", call. = FALSE)

  var <- apply(data, 1, sd)
  if (n == 0) {
    df <- data
    n <- nrow(data)
  } else if (n > nrow(data)) {
    df <- data
    n <- nrow(data)
  } else {
    df <- data[order(var, decreasing = TRUE)[seq_len(n)], , drop = FALSE]
  }

  pca <- prcomp(t(df), scale = scale)
  pca_df <- pca$x %>%
    data.frame() %>%
    tibble::rownames_to_column() %>%
    dplyr::left_join(., data.frame(colData(dep)), by = c("rowname" = ID_col))

  percent <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)

  for (feat in indicate) {
    if (feat %in% colnames(pca_df)) pca_df[[feat]] <- as.factor(pca_df[[feat]])
  }

  p <- ggplot2::ggplot(pca_df, ggplot2::aes(get(paste0("PC", x)), get(paste0("PC", y)))) +
    ggplot2::labs(
      title = paste0("PCA plot - top ", n, " variable features"),
      x = paste0("PC", x, ": ", percent[x], "%"),
      y = paste0("PC", y, ": ", percent[y], "%")
    ) +
    ggplot2::coord_fixed() +
    theme_DEP1()

  if (length(indicate) == 0 || !indicate[1] %in% colnames(pca_df)) {
    p <- p + ggplot2::geom_point(size = point_size)
  } else if (length(indicate) >= 2 && indicate[2] %in% colnames(pca_df)) {
    p <- p + ggplot2::geom_point(ggplot2::aes(col = .data[[indicate[1]]], shape = .data[[indicate[2]]]),
      size = point_size) +
      ggplot2::labs(col = indicate[1], shape = indicate[2])
  } else {
    p <- p + ggplot2::geom_point(ggplot2::aes(col = .data[[indicate[1]]]), size = point_size) +
      ggplot2::labs(col = indicate[1])
  }

  if (plot) p else pca_df
}

# theme_DEP1 (DEP-style theme_bw)
theme_DEP1 <- function() {
  basesize <- 12
  theme <- ggplot2::theme_bw(base_size = basesize)
  theme$plot.title$face <- "bold"
  theme$plot.title$size <- basesize + 2
  theme$plot.title$hjust <- 0.5
  theme$axis.title.x$size <- basesize + 2
  theme$axis.title.y$size <- basesize + 2
  theme$axis.text$size <- basesize
  theme$axis.text$colour <- "black"
  theme$legend.title$size <- basesize + 2
  theme$legend.text$size <- basesize
  theme$strip.text$face <- "bold"
  theme$strip.text$size <- basesize + 2
  theme$strip.text$colour <- "black"
  theme
}

# plot_cor_customized: from FragPipe-Analyst.
plot_cor_customized <- function(dep, significant = FALSE, lower = -1, upper = 1,
  pal = "PRGn", pal_rev = FALSE, indicate = "condition", font_size = 12, plot = FALSE, ...) {
  assertthat::assert_that(inherits(dep, "SummarizedExperiment"),
    is.logical(significant), is.numeric(font_size), is.logical(plot))

  if (significant && "significant" %in% colnames(rowData(dep, use.names = FALSE))) {
    dep <- dep[tidyr::replace_na(rowData(dep, use.names = FALSE)$significant, FALSE), ]
  }

  data <- assay(dep)
  temp <- as.data.frame(colData(dep))
  if ("label" %in% colnames(temp) && "sample_name" %in% colnames(temp)) {
    rownames(temp) <- temp$label
    new_names <- temp[colnames(data), "sample_name"]
    if (!any(is.na(new_names))) colnames(data) <- new_names
  }
  cn <- colnames(data)
  cn <- gsub("_MaxLFQ\\.Intensity$| MaxLFQ\\.Intensity$", "", cn)
  cn <- gsub("_Intensity$| Intensity$", "", cn)
  cn <- gsub("_Spectral\\.Count$| Spectral\\.Count$", "", cn)
  cn <- gsub("_{2,}", "_", cn)
  colnames(data) <- make.unique(trimws(cn), sep = "_")

  cor_mat <- cor(data, use = "complete.obs")
  lower <- min(cor_mat)
  upper <- max(cor_mat)

  ha1 <- NULL
  if (!is.null(indicate) && indicate %in% colnames(temp)) {
    anno <- as.data.frame(colData(dep)) %>% dplyr::select(dplyr::all_of(indicate))
    var <- sort(unique(anno[[1]]))
    nv <- length(var)
    cols <- if (nv == 1) c("black") else if (nv == 2) c("orangered", "cornflowerblue") else
      if (nv <= 6) RColorBrewer::brewer.pal(max(3, nv), "Pastel1")[seq_len(nv)] else
      if (nv <= 12) RColorBrewer::brewer.pal(nv, "Set3") else
      colorRampPalette(RColorBrewer::brewer.pal(12, "Set3"))(nv)
    names(cols) <- var
    ha1 <- ComplexHeatmap::HeatmapAnnotation(df = anno, col = setNames(list(cols), indicate), show_annotation_name = TRUE)
  }

  ht1 <- ComplexHeatmap::Heatmap(cor_mat,
    col = circlize::colorRamp2(seq(lower, upper, (upper - lower) / 7),
      if (pal_rev) rev(RColorBrewer::brewer.pal(8, pal)) else RColorBrewer::brewer.pal(8, pal)),
    heatmap_legend_param = list(color_bar = "continuous", legend_direction = "horizontal",
      legend_width = grid::unit(5, "cm"), title_position = "topcenter"),
    name = "Pearson correlation",
    column_names_gp = grid::gpar(fontsize = font_size),
    row_names_gp = grid::gpar(fontsize = font_size),
    top_annotation = ha1, ...)
  if (plot) ComplexHeatmap::draw(ht1, heatmap_legend_side = "top")
  ht1
}

# plot_missval_customized: from FragPipe-Analyst.
plot_missval_customized <- function(se) {
  assertthat::assert_that(inherits(se, "SummarizedExperiment"))
  se_assay <- assay(se)
  if (!any(is.na(se_assay))) stop("No missing values", call. = FALSE)
  df <- se_assay %>% data.frame(check.names = FALSE)
  missval <- df[apply(df, 1, function(x) any(is.na(x))), , drop = FALSE]
  missval <- ifelse(is.na(missval), 0, 1)
  temp <- as.data.frame(colData(se))
  if ("label" %in% colnames(temp) && "sample_name" %in% colnames(temp)) {
    # Use rownames that match assay colnames (we use sample_name; Monash uses label)
    idx_col <- if (all(colnames(missval) %in% temp$sample_name)) "sample_name"
      else if (all(colnames(missval) %in% temp$label)) "label"
      else "label"
    rownames(temp) <- temp[[idx_col]]
    new_cn <- temp[colnames(missval), "sample_name"]
    if (!any(is.na(new_cn))) colnames(missval) <- new_cn
  }
  # Strip LFQ suffix for display (match FragPipe-Analyst / plot_cor_customized)
  cn <- colnames(missval)
  cn <- gsub("_MaxLFQ\\.Intensity$| MaxLFQ\\.Intensity$", "", cn)
  cn <- gsub("_Intensity$| Intensity$", "", cn)
  cn <- gsub("_Spectral\\.Count$| Spectral\\.Count$", "", cn)
  cn <- gsub("_{2,}", "_", cn)
  colnames(missval) <- make.unique(trimws(cn), sep = "_")
  nfeat <- dim(missval)[1]
  feat_label <- if (!is.null(metadata(se)$level)) {
    switch(metadata(se)$level, protein = "proteins", peptide = "peptides", site = "sites", "features")
  } else "features"
  if (nfeat >= 65536 && requireNamespace("factoextra", quietly = TRUE) && requireNamespace("fastcluster", quietly = TRUE)) {
    dist <- factoextra::get_dist(missval, "euclidean")
    mat.hc <- fastcluster::hclust(dist, method = "complete")
    ht2 <- ComplexHeatmap::Heatmap(missval, col = c("#FFFFFF", "#000000"),
      cluster_rows = as.dendrogram(mat.hc), column_names_side = "top",
      show_row_names = FALSE, show_column_names = TRUE, show_row_dend = FALSE,
      name = paste0("Missing values pattern (", nfeat, " ", feat_label, ")"),
      column_names_gp = grid::gpar(fontsize = 16),
      heatmap_legend_param = list(at = c(0, 1), labels = c("Missing value", "Valid value")))
  } else {
    ht2 <- ComplexHeatmap::Heatmap(missval, col = c("#FFFFFF", "#000000"),
      column_names_side = "top", show_row_names = FALSE, show_column_names = TRUE,
      name = paste0("Missing values pattern (", nfeat, " ", feat_label, ")"),
      column_names_gp = grid::gpar(fontsize = 16),
      heatmap_legend_param = list(at = c(0, 1), labels = c("Missing value", "Valid value")))
  }
  ComplexHeatmap::draw(ht2, heatmap_legend_side = "top")
}

# plot_cvs_custom: from FragPipeAnalystR. Sample CV distribution per condition.
coef_variation <- function(x) {
  m <- mean(x, na.rm = TRUE)
  if (is.na(m) || m == 0) return(NA)
  sd(x, na.rm = TRUE) / m
}

plot_cvs_custom <- function(se, id = "sample_name", scale = TRUE, check.names = FALSE) {
  assertthat::assert_that(inherits(se, "SummarizedExperiment"))
  untransformed_intensity <- 2^(assay(se))
  exp_design <- as.data.frame(colData(se))
  cvs_group <- untransformed_intensity %>%
    data.frame(check.names = check.names) %>%
    tibble::rownames_to_column() %>%
    tidyr::gather("ID", "Intensity", -rowname) %>%
    dplyr::left_join(exp_design, by = c("ID" = id)) %>%
    dplyr::group_by(rowname, condition) %>%
    dplyr::summarise(cvs = coef_variation(Intensity), .groups = "drop") %>%
    dplyr::group_by(condition) %>%
    dplyr::mutate(condition_median = median(cvs, na.rm = TRUE))
  cvs_group <- cvs_group %>% dplyr::filter(!is.na(cvs) & is.finite(cvs))
  if (nrow(cvs_group) == 0) stop("No valid CVs to plot")
  if (scale) {
    p1 <- ggplot2::ggplot(cvs_group, ggplot2::aes(cvs, color = condition, fill = condition)) +
      ggplot2::geom_histogram(alpha = 0.5, bins = 20, show.legend = FALSE) +
      ggplot2::facet_wrap(~condition) +
      ggplot2::geom_vline(ggplot2::aes(xintercept = condition_median, group = condition),
        color = "grey40", linetype = "dashed") +
      ggplot2::scale_x_continuous(labels = scales::percent, limits = c(0, 1)) +
      ggplot2::labs(title = "Sample Coefficient of Variation", x = "Coefficient of Variation", y = "Count") +
      ggplot2::theme_bw() +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
        axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5))
  } else {
    p1 <- ggplot2::ggplot(cvs_group, ggplot2::aes(cvs, color = condition, fill = condition)) +
      ggplot2::geom_histogram(alpha = 0.5, bins = 20, show.legend = FALSE) +
      ggplot2::facet_wrap(~condition) +
      ggplot2::geom_vline(ggplot2::aes(xintercept = condition_median, group = condition),
        color = "grey40", linetype = "dashed") +
      ggplot2::scale_x_continuous(labels = scales::percent) +
      ggplot2::labs(title = "Sample Coefficient of Variation", x = "Coefficient of Variation", y = "Count") +
      ggplot2::theme_bw() +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
        axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5))
  }
  ymax_vals <- ggplot2::ggplot_build(p1)$data[[1]]$ymax
  ymax <- if (length(ymax_vals) > 0) max(ymax_vals, na.rm = TRUE) * 1.1 else 1
  label_df <- cvs_group %>% dplyr::distinct(condition, condition_median)
  p1 + ggplot2::geom_text(ggplot2::aes(x = 0.9, y = ymax, color = condition,
    label = paste0("Median = ", round(condition_median, 2) * 100, "%")),
    show.legend = FALSE, size = 4, inherit.aes = FALSE,
    data = label_df)
}

# plot_feature_numbers_custom: from FragPipe-Analyst.
plot_feature_numbers_custom <- function(se, fill = "condition") {
  assertthat::assert_that(inherits(se, "SummarizedExperiment"))
  df <- assay(se) %>% data.frame(check.names = FALSE) %>% tibble::rownames_to_column() %>%
    tidyr::gather(ID, bin, -rowname) %>% dplyr::mutate(bin = ifelse(is.na(bin), 0, 1))
  stat <- df %>% dplyr::group_by(ID) %>% dplyr::summarize(n = dplyr::n(), sum = sum(bin))
  cd <- as.data.frame(colData(se))
  # Use column that matches assay colnames (we use sample_name for display; Monash uses label)
  id_col <- if ("sample_name" %in% colnames(cd) && all(stat$ID %in% cd$sample_name)) "sample_name"
    else if ("label" %in% colnames(cd) && all(stat$ID %in% cd$label)) "label"
    else if ("label" %in% colnames(cd)) "label"
    else "sample_name"
  stat <- dplyr::left_join(stat, cd, by = c("ID" = id_col))
  # Strip LFQ suffix for display (match FragPipe-Analyst)
  stat$display_name <- gsub("_MaxLFQ\\.Intensity$| MaxLFQ\\.Intensity$", "", stat$ID)
  stat$display_name <- gsub("_Intensity$| Intensity$", "", stat$display_name)
  stat$display_name <- gsub("_Spectral\\.Count$| Spectral\\.Count$", "", stat$display_name)
  stat$display_name <- make.unique(trimws(gsub("_{2,}", "_", stat$display_name)), sep = "_")
  ggplot2::ggplot(stat, ggplot2::aes(x = display_name, y = sum, fill = .data[[fill]])) +
    ggplot2::geom_col() + ggplot2::geom_hline(yintercept = unique(stat$n), linetype = "solid") +
    ggplot2::labs(title = "Features per sample", x = "", y = "Number of features") +
    theme_DEP1() + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1))
}

# plot_coverage_customized: barplot of features in 1, 2, ..., N samples
plot_coverage_customized <- function(se, plot = TRUE) {
  assertthat::assert_that(inherits(se, "SummarizedExperiment"), is.logical(plot))
  df <- assay(se) %>% data.frame() %>% tibble::rownames_to_column() %>%
    tidyr::gather(ID, bin, -rowname) %>% dplyr::mutate(bin = ifelse(is.na(bin), 0, 1))
  stat <- df %>% dplyr::group_by(rowname) %>% dplyr::summarize(sum = sum(bin))
  table <- data.frame(table(stat$sum))
  p <- ggplot2::ggplot(table, ggplot2::aes(x = "all", y = Freq, fill = Var1)) +
    ggplot2::geom_col(col = "white") + ggplot2::scale_fill_grey(start = 0.8, end = 0.2) +
    ggplot2::labs(title = "Feature coverage", x = "", y = "Number of features", fill = "Samples") + theme_DEP1()
  if (plot) p else { colnames(table) <- c("samples", "features"); table }
}

# plot_density_custom: overlaid density curves by condition
plot_density_custom <- function(ses) {
  gather_join <- function(se) {
    cd <- as.data.frame(colData(se))
    cd$..colkey.. <- rownames(cd)
    assay(se) %>% data.frame(check.names = FALSE) %>% tidyr::gather(ID, val, dplyr::everything()) %>%
      dplyr::left_join(cd, by = c("ID" = "..colkey.."))
  }
  df <- purrr::map_df(ses, gather_join, .id = "var") %>%
    dplyr::mutate(var = factor(var, levels = names(ses)))
  ggplot2::ggplot(df, ggplot2::aes(val, col = condition)) +
    ggplot2::geom_density(na.rm = TRUE) +
    ggplot2::facet_wrap(~var, ncol = 1, strip.position = "top") +
    ggplot2::labs(x = expression(log[2] ~ "Intensity"), y = "Density") + theme_DEP1()
}

# ---- Comparison plots (Jaccard, Venn, UpSet, feature) ----

# data_attendance_custom: occurrence matrix for Venn/UpSet. From FragPipe-Analyst.
data_attendance_custom <- function(se, exp = "LFQ", level = "protein") {
  assertthat::assert_that(inherits(se, "SummarizedExperiment"))
  df <- as.data.frame(assay(se), check.names = FALSE)
  col_data <- as.data.frame(colData(se))
  sample_cols <- colnames(df)
  conditions <- unique(col_data$condition)
  id_col <- if ("label" %in% colnames(col_data)) "label" else "sample_name"
  if (exp == "LFQ") {
    rd <- as.data.frame(rowData(se))
    df$Gene <- as.character(if ("Gene" %in% colnames(rd)) rd$Gene else if ("Genes" %in% colnames(rd)) rd$Genes else "NoGeneNameAvailable")
    df$Protein <- as.character(if ("Protein ID" %in% colnames(rd)) rd[["Protein ID"]] else if ("Protein" %in% colnames(rd)) rd$Protein else rd$ID)
    if (any(df$Gene == "" | is.na(df$Gene), na.rm = TRUE)) df$Gene[df$Gene == "" | is.na(df$Gene)] <- "NoGeneNameAvailable"
    df <- df[rowSums(!is.na(df[, sample_cols, drop = FALSE])) != 0, ]
    for (i in seq_along(conditions)) {
      cond <- conditions[i]
      temp <- col_data[col_data$condition == cond, , drop = FALSE]
      sel_cols <- intersect(rownames(temp), colnames(df))
      if (length(sel_cols) > 0) df[[paste0("#Occurences_", cond)]] <- rowSums(!is.na(df[, sel_cols, drop = FALSE]))
    }
    df <- dplyr::relocate(df, Protein, Gene, .before = 1)
  } else if (exp == "DIA" && level == "protein") {
    df$Gene <- as.character(rowData(se)$Genes)
    df$Protein <- as.character(if ("Protein.Ids" %in% colnames(rowData(se))) rowData(se)$Protein.Ids else rowData(se)$Protein.Group)
    if (any(df$Gene == "" | is.na(df$Gene), na.rm = TRUE)) df$Gene[df$Gene == "" | is.na(df$Gene)] <- "NoGeneNameAvailable"
    df <- df[rowSums(!is.na(df[, sample_cols, drop = FALSE])) != 0, ]
    for (i in seq_along(conditions)) {
      cond <- conditions[i]
      temp <- col_data[col_data$condition == cond, , drop = FALSE]
      sel_cols <- intersect(temp[[id_col]], colnames(df))
      if (length(sel_cols) > 0) df[[paste0("#Occurences_", cond)]] <- rowSums(!is.na(df[, sel_cols, drop = FALSE]))
    }
    df <- dplyr::relocate(df, Protein, Gene, .before = 1)
  } else {
    rd <- as.data.frame(rowData(se))
    df$Gene <- as.character(if ("Gene" %in% colnames(rd)) rd$Gene else if ("Genes" %in% colnames(rd)) rd$Genes else "NoGeneNameAvailable")
    if (any(df$Gene == "" | is.na(df$Gene), na.rm = TRUE)) df$Gene[df$Gene == "" | is.na(df$Gene)] <- "NoGeneNameAvailable"
    df <- df[rowSums(!is.na(df[, sample_cols, drop = FALSE])) != 0, ]
    for (i in seq_along(conditions)) {
      cond <- conditions[i]
      temp <- col_data[col_data$condition == cond, , drop = FALSE]
      sel_cols <- intersect(temp[[id_col]], colnames(df))
      if (length(sel_cols) > 0) df[[paste0("#Occurences_", cond)]] <- rowSums(!is.na(df[, sel_cols, drop = FALSE]))
    }
    df <- dplyr::relocate(df, Gene, .before = 1)
  }
  rownames(df) <- NULL
  df
}

# plot_Jaccard_custom: sample-level Jaccard similarity heatmap. From FragPipe-Analyst.
plot_Jaccard_custom <- function(dep, plot = TRUE, indicate = "condition") {
  assertthat::assert_that(inherits(dep, "SummarizedExperiment"))
  if (!requireNamespace("vegan", quietly = TRUE)) return(NULL)
  data <- assay(dep)
  cd <- as.data.frame(colData(dep))
  cn <- colnames(data)
  new_cn <- if ("label" %in% colnames(cd) && "sample_name" %in% colnames(cd)) {
    idx <- match(cn, cd$label)
    if (!any(is.na(idx))) cd$sample_name[idx] else cn
  } else if ("sample_name" %in% colnames(cd) && all(cn %in% rownames(cd))) {
    cd[cn, "sample_name"]
  } else cn
  new_cn <- gsub("_MaxLFQ\\.Intensity$| MaxLFQ\\.Intensity$", "", new_cn)
  new_cn <- gsub("_Intensity$| Intensity$", "", new_cn)
  new_cn <- gsub("_Spectral\\.Count$| Spectral\\.Count$", "", new_cn)
  colnames(data) <- make.unique(trimws(gsub("_{2,}", "_", new_cn)), sep = "_")
  cor_mat <- 1 - as.matrix(vegan::vegdist(t(data), method = "jaccard", na.rm = TRUE))
  lower <- min(cor_mat); upper <- max(cor_mat)
  ha1 <- NULL
  if (!is.null(indicate) && indicate %in% colnames(cd)) {
    anno <- cd[, indicate, drop = FALSE]; rownames(anno) <- colnames(cor_mat)
    var <- unique(anno[[1]])
    cols <- if (length(var) == 1) c("black") else if (length(var) == 2) c("orangered", "cornflowerblue") else
      if (length(var) < 7) RColorBrewer::brewer.pal(max(3, length(var)), "Pastel1")[seq_len(length(var))] else
      RColorBrewer::brewer.pal(length(var), "Set3")
    names(cols) <- var
    ha1 <- ComplexHeatmap::HeatmapAnnotation(df = anno, col = setNames(list(cols), indicate), show_annotation_name = TRUE)
  }
  ht1 <- ComplexHeatmap::Heatmap(cor_mat,
    col = circlize::colorRamp2(c(lower, (upper + lower) / 2, upper), c("blue", "lightyellow", "red")),
    heatmap_legend_param = list(color_bar = "continuous", legend_direction = "horizontal",
      legend_width = grid::unit(5, "cm"), title_position = "topcenter"),
    name = "Jaccard similarity", column_names_gp = grid::gpar(fontsize = 12), row_names_gp = grid::gpar(fontsize = 12),
    top_annotation = ha1)
  if (plot) ComplexHeatmap::draw(ht1, heatmap_legend_side = "top") else ht1
}

# plot_venn_custom: pairwise Venn (ggVennDiagram)
plot_venn_custom <- function(df, cond1, cond2, cond3 = NULL) {
  if (!requireNamespace("ggVennDiagram", quietly = TRUE)) return(NULL)
  occ1 <- paste0("#Occurences_", cond1); occ2 <- paste0("#Occurences_", cond2)
  if (!occ1 %in% colnames(df) || !occ2 %in% colnames(df)) return(NULL)
  set1 <- df[df[[occ1]] != 0, "Gene"]; set2 <- df[df[[occ2]] != 0, "Gene"]
  x <- list(set1, set2); names(x) <- c(cond1, cond2)
  if (!is.null(cond3) && cond3 != "NONE" && paste0("#Occurences_", cond3) %in% colnames(df)) {
    set3 <- df[df[[paste0("#Occurences_", cond3)]] != 0, "Gene"]
    x <- list(set1, set2, set3); names(x) <- c(cond1, cond2, cond3)
  }
  max_cond_len <- max(nchar(c(cond1, cond2, if (!is.null(cond3) && cond3 != "NONE") cond3 else character(0))), 0)
  set_size <- if (max_cond_len > 10) round(max(3, 12 - (max_cond_len - 10) / 3)) else NULL
  venn_args <- list(x = x, label_alpha = 0)
  if (!is.null(set_size)) venn_args$set_size <- set_size
  do.call(ggVennDiagram::ggVennDiagram, venn_args) +
    ggplot2::scale_fill_gradient(low = "#F4FAFE", high = "#4981BF") +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = 0.3)) +
    ggplot2::coord_flip() +
    ggplot2::theme(plot.margin = ggplot2::margin(12, 24, 12, 12, "pt"))
}

# plot_upset_custom: UpSetR (grid-based; must print() to render to device)
plot_upset_custom <- function(df) {
  if (!requireNamespace("UpSetR", quietly = TRUE)) return(invisible(NULL))
  df <- df[, grep("Occurences", colnames(df)), drop = FALSE]
  df <- ifelse(df != 0, 1, 0)
  df <- data.frame(df)
  colnames(df) <- gsub("X.Occurences_|#Occurences_", "", colnames(df))
  if (sum(colSums(df) != 0) <= 1) return(invisible(NULL))
  p <- UpSetR::upset(df, nsets = ncol(df), mb.ratio = c(0.6, 0.4), text.scale = 1.5, point.size = 3,
    order.by = "freq", decreasing = TRUE, nintersects = NA, mainbar.y.label = "#Features in intersection",
    sets.x.label = "#Features", set_size.scale_max = nrow(df) + 1000, set_size.show = TRUE)
  print(p)
  invisible(NULL)
}

# plot_feature_custom: boxplot or violin per feature. From FragPipe-Analyst.
plot_feature_custom <- function(dep, protein, type = "boxplot", id = NULL, show_gene = FALSE) {
  assertthat::assert_that(inherits(dep, "SummarizedExperiment"), is.character(protein), is.character(type))
  subset <- dep[protein, ]
  cd <- as.data.frame(colData(subset))
  assay_ids <- unique(colnames(assay(subset)))
  # Use column that matches assay colnames (make_se renames to sample_name; Monash uses label)
  id_col <- if (!is.null(id)) id
    else if ("sample_name" %in% colnames(cd) && all(assay_ids %in% cd$sample_name)) "sample_name"
    else if ("label" %in% colnames(cd) && all(assay_ids %in% cd$label)) "label"
    else if ("label" %in% colnames(cd)) "label"
    else "sample_name"
  df_reps <- data.frame(assay(subset), check.names = FALSE) %>%
    tibble::rownames_to_column() %>%
    tidyr::gather(ID, val, -rowname) %>%
    dplyr::left_join(cd, by = c("ID" = id_col))
  df_reps$rowname <- factor(as.character(df_reps$rowname), levels = protein)
  df_reps$condition <- as.factor(df_reps$condition)
  df_reps <- df_reps[!is.na(df_reps$val), ]
  if ("replicate" %in% colnames(df_reps)) {
    df_reps$replicate[is.na(df_reps$replicate)] <- 1L
    df_reps$replicate <- as.character(df_reps$replicate)
  }
  if (show_gene && nrow(df_reps) > 0) {
    rd <- rowData(subset)
    md <- metadata(dep)
    lvl <- if (!is.null(md$level)) md$level else "protein"
    if (!lvl %in% c("site", "peptide")) {
      df_reps$rowname <- rd[df_reps$rowname, "name"]
    } else {
      if (!is.null(md$exp) && md$exp == "DIA") {
        gene_col <- if ("Gene" %in% colnames(rd)) "Gene" else "Genes"
        if (lvl == "site") {
          df_reps$rowname <- paste0(rd[df_reps$rowname, gene_col], "_", gsub(".*_", "", df_reps$rowname))
        } else {
          df_reps$rowname <- paste0(rd[df_reps$rowname, gene_col], "_", gsub(".*_", "", df_reps$rowname))
        }
      } else if (!is.null(md$exp) && md$exp == "TMT") {
        gene_col <- if ("Gene" %in% colnames(rd)) "Gene" else "Genes"
        if (lvl == "site") {
          df_reps$rowname <- paste0(rd[df_reps$rowname, gene_col], "_", gsub(".*_", "", rd[df_reps$rowname, "ID"]))
        } else {
          df_reps$rowname <- paste0(rd[df_reps$rowname, gene_col], "_", gsub(".*_", "", rd[df_reps$rowname, "Peptide"]))
        }
      } else {
        gene_col <- if ("Gene" %in% colnames(rd)) "Gene" else "Genes"
        seq_col <- if ("Modified Sequence" %in% colnames(rd)) "Modified Sequence" else "Peptide.Sequence"
        df_reps$rowname <- paste0(rd[df_reps$rowname, gene_col], "_", rd[df_reps$rowname, seq_col])
      }
    }
  }
  nrep <- if ("replicate" %in% colnames(df_reps)) length(unique(as.character(df_reps$replicate))) else 1
  max_cond_len <- max(nchar(as.character(unique(df_reps$condition))), 0)
  th <- ggplot2::theme(axis.title.x = ggplot2::element_blank(), panel.border = ggplot2::element_blank(),
    panel.grid.major.x = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank(),
    panel.grid.major.y = ggplot2::element_line(colour = "grey90"),
    axis.line = ggplot2::element_line(colour = "black"))
  if (max_cond_len > 10) {
    x_size <- round(max(3, 12 - (max_cond_len - 10) / 3))
    th <- th + ggplot2::theme(axis.text.x = ggplot2::element_text(size = x_size))
  }
  if (type == "violin") {
    if (nrep <= 1) {
      p <- ggplot2::ggplot(df_reps, ggplot2::aes(condition, val)) +
        ggplot2::geom_violin(fill = "grey90", scale = "width", draw_quantiles = 0.5, trim = TRUE) +
        ggplot2::geom_jitter(size = 3, position = ggplot2::position_dodge(width = 0.3)) +
        ggplot2::labs(y = expression(log[2] ~ "Intensity")) + ggplot2::facet_wrap(~rowname) + ggplot2::theme_bw() + th
    } else {
      p <- ggplot2::ggplot(df_reps, ggplot2::aes(condition, val)) +
        ggplot2::geom_violin(fill = "grey90", scale = "width", draw_quantiles = 0.5, trim = TRUE) +
        ggplot2::geom_jitter(ggplot2::aes(color = factor(replicate)), size = 3, position = ggplot2::position_dodge(width = 0.3)) +
        ggplot2::labs(y = expression(log[2] ~ "Intensity"), col = "Replicate") + ggplot2::facet_wrap(~rowname) +
        ggplot2::scale_color_brewer(palette = "Dark2") + ggplot2::theme_bw() + th
    }
  } else {
    if (nrep <= 1) {
      p <- ggplot2::ggplot(df_reps, ggplot2::aes(condition, val)) +
        ggplot2::geom_boxplot() + ggplot2::geom_jitter(size = 3, position = ggplot2::position_dodge(width = 0.3)) +
        ggplot2::labs(y = expression(log[2] ~ "Intensity")) + ggplot2::facet_wrap(~rowname) + ggplot2::theme_bw() + th
    } else {
      p <- ggplot2::ggplot(df_reps, ggplot2::aes(condition, val)) +
        ggplot2::geom_boxplot() +
        ggplot2::geom_jitter(ggplot2::aes(color = factor(replicate)), size = 3, position = ggplot2::position_dodge(width = 0.3)) +
        ggplot2::labs(y = expression(log[2] ~ "Intensity"), col = "Replicate") + ggplot2::facet_wrap(~rowname) +
        ggplot2::scale_color_brewer(palette = "Dark2") + ggplot2::theme_bw() + th
    }
  }
  p
}
