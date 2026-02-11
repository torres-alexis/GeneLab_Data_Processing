# FragPipe-Analyst DE and plot functions (library)
# Source: MonashProteomics/FragPipe-Analyst (GPL-3.0) https://github.com/MonashProteomics/FragPipe-Analyst
# global_filter, filter_by_condition: from R/filter.R
# GUI calls: global_filter(se, 100 - min_global_appearance), filter_by_condition(se, min_appearance_each_condition)
global_filter <- function(se, percentage = 50) {
  percentage <- percentage / 100
  ridx <- rowSums(is.na(assay(se))) / ncol(assay(se)) <= percentage
  se <- se[ridx, ]
  return(se)
}

filter_by_condition <- function(se, min_percentage = 50) {
  min_percentage <- min_percentage / 100
  conditions <- unique(colData(se)$condition)
  row_ids <- rep(0, nrow(assay(se)))
  for (c in conditions) {
    se_c <- se[, colData(se)$condition == c]
    ridx <- rowSums(!is.na(assay(se_c))) / ncol(assay(se_c)) >= min_percentage
    row_ids <- row_ids + ridx
  }
  se <- se[row_ids > 0, ]
  return(se)
}

# plot_feature_monash: based on plot_feature from MonashProteomics/FragPipe-Analyst R/customized.R
# Always subset by protein IDs (rownames); show_gene replaces facet labels with name/Gene from rowData.
plot_feature_monash <- function(dep, protein, type = "boxplot", id = "sample_name", show_gene = FALSE) {
  assertthat::assert_that(inherits(dep, "SummarizedExperiment"),
                          is.character(protein),
                          is.character(type))
  subset <- dep[protein, ]
  df_reps <- data.frame(assay(subset), check.names = FALSE) %>%
    tibble::rownames_to_column() %>%
    tidyr::gather(ID, val, -rowname) %>%
    dplyr::left_join(., data.frame(colData(subset)), by = c("ID" = id))
  df_reps$rowname <- factor(as.character(df_reps$rowname), levels = protein)
  df_CI <- df_reps %>%
    dplyr::group_by(condition, rowname) %>%
    dplyr::summarize(mean = mean(val, na.rm = TRUE),
              sd = sd(val, na.rm = TRUE),
              n = dplyr::n()) %>%
    dplyr::mutate(error = qnorm(0.975) * sd / sqrt(n),
           CI.L = mean - error,
           CI.R = mean + error) %>%
    as.data.frame()
  df_CI$rowname <- factor(as.character(df_CI$rowname), levels = protein)
  df_reps$condition <- as.factor(df_reps$condition)
  df_reps <- df_reps[!is.na(df_reps$val), ]
  if ("replicate" %in% colnames(df_reps)) {
    df_reps$replicate[is.na(df_reps$replicate)] <- 1L
    df_reps$replicate <- as.character(df_reps$replicate)
  }
  if (show_gene) {
    md <- metadata(dep)
    level <- if (!is.null(md$level)) md$level else "protein"
    if (!level %in% c("site", "peptide")) {
      df_reps$rowname <- rowData(subset)[df_reps$rowname, "name"]
    } else {
      if (!is.null(md$exp) && md$exp == "DIA") {
        if (level == "site") {
          df_reps$rowname <- paste0(rowData(subset)[df_reps$rowname, "Gene"], "_", gsub(".*_", "", df_reps$rowname))
        } else {
          gene_col <- if ("Gene" %in% colnames(rowData(subset))) "Gene" else "Genes"
          df_reps$rowname <- paste0(rowData(subset)[df_reps$rowname, gene_col], "_", gsub(".*_", "", df_reps$rowname))
        }
      } else if (!is.null(md$exp) && md$exp == "TMT") {
        if (level == "site") {
          df_reps$rowname <- paste0(rowData(subset)[df_reps$rowname, "Gene"], "_", gsub(".*_", "", rowData(subset)[df_reps$rowname, "ID"]))
        } else {
          df_reps$rowname <- paste0(rowData(subset)[df_reps$rowname, "Gene"], "_", gsub(".*_", "", rowData(subset)[df_reps$rowname, "Peptide"]))
        }
      } else {
        seq_col <- if ("Modified Sequence" %in% colnames(rowData(subset))) "Modified Sequence" else "Peptide.Sequence"
        df_reps$rowname <- paste0(rowData(subset)[df_reps$rowname, "Gene"], "_", rowData(subset)[df_reps$rowname, seq_col])
      }
    }
  }
  nrep <- if ("replicate" %in% colnames(df_reps)) max(as.numeric(df_reps$replicate), na.rm = TRUE) else 1
  if (type == "violin") {
    if (nrep <= 1) {
      p <- ggplot2::ggplot(df_reps, ggplot2::aes(condition, val)) +
        ggplot2::geom_violin(fill = "grey90", scale = "width", draw_quantiles = 0.5, trim = TRUE) +
        ggplot2::geom_jitter(size = 3, position = ggplot2::position_dodge(width = 0.3)) +
        ggplot2::labs(y = expression(log[2]~"Intensity")) +
        ggplot2::facet_wrap(~rowname) +
        ggplot2::theme_bw() +
        ggplot2::theme(axis.title.x = ggplot2::element_blank(),
              panel.border = ggplot2::element_blank(), panel.grid.major = ggplot2::element_blank(),
              panel.grid.minor = ggplot2::element_blank(), axis.line = ggplot2::element_line(colour = "black"))
    } else {
      p <- ggplot2::ggplot(df_reps, ggplot2::aes(condition, val)) +
        ggplot2::geom_violin(fill = "grey90", scale = "width", draw_quantiles = 0.5, trim = TRUE) +
        ggplot2::geom_jitter(ggplot2::aes(color = factor(replicate)), size = 3, position = ggplot2::position_dodge(width = 0.3)) +
        ggplot2::labs(y = expression(log[2]~"Intensity"), col = "Replicates") +
        ggplot2::facet_wrap(~rowname) +
        ggplot2::scale_color_brewer(palette = "Dark2") +
        ggplot2::theme_bw() +
        ggplot2::theme(axis.title.x = ggplot2::element_blank(),
              panel.border = ggplot2::element_blank(), panel.grid.major = ggplot2::element_blank(),
              panel.grid.minor = ggplot2::element_blank(), axis.line = ggplot2::element_line(colour = "black"))
    }
  } else {
    if (nrep <= 1) {
      p <- ggplot2::ggplot(df_reps, ggplot2::aes(condition, val)) +
        ggplot2::geom_boxplot() +
        ggplot2::geom_jitter(size = 3, position = ggplot2::position_dodge(width = 0.3)) +
        ggplot2::labs(y = expression(log[2]~"Intensity")) +
        ggplot2::facet_wrap(~rowname) +
        ggplot2::theme_bw() +
        ggplot2::theme(axis.title.x = ggplot2::element_blank(),
              panel.border = ggplot2::element_blank(), panel.grid.major = ggplot2::element_blank(),
              panel.grid.minor = ggplot2::element_blank(), axis.line = ggplot2::element_line(colour = "black"))
    } else {
      p <- ggplot2::ggplot(df_reps, ggplot2::aes(condition, val)) +
        ggplot2::geom_boxplot() +
        ggplot2::geom_jitter(ggplot2::aes(color = factor(replicate)), size = 3, position = ggplot2::position_dodge(width = 0.3)) +
        ggplot2::labs(y = expression(log[2]~"Intensity"), col = "Replicates") +
        ggplot2::facet_wrap(~rowname) +
        ggplot2::scale_color_brewer(palette = "Dark2") +
        ggplot2::theme_bw() +
        ggplot2::theme(axis.title.x = ggplot2::element_blank(),
              panel.border = ggplot2::element_blank(), panel.grid.major = ggplot2::element_blank(),
              panel.grid.minor = ggplot2::element_blank(), axis.line = ggplot2::element_line(colour = "black"))
    }
  }
  return(p)
}

# from MonashProteomics/FragPipe-Analyst R/functions.R # 305 (test_limma, BH FDR)
# Modified: "others" type, left_join by ID. Matches Shiny limma path.
# ---- test_limma_customized: Benjamini Hochberg FDR (limma topTable adjust.method="BH") ----
test_limma_customized <- function(se, type = c("control", "all", "others", "manual"),
                                  control = NULL, test = NULL,
                                  design_formula = formula(~ 0 + condition),
                                  paired = FALSE) {
  assertthat::assert_that(inherits(se, "SummarizedExperiment"),
                          is.character(type),
                          class(design_formula) == "formula")
  if (!paired) {
    design_formula <- design_formula
  } else {
    design_formula <- formula(~ 0 + condition + replicate)
  }
  type <- match.arg(type)

  col_data <- colData(se)
  raw <- assay(se)

  if (any(!c("name", "ID") %in% colnames(rowData(se, use.names = FALSE)))) {
    stop("'name' and/or 'ID' columns are not present in '", deparse(substitute(se)),
         "'\nRun make_unique() and make_se() to obtain the required columns", call. = FALSE)
  }
  if (any(!c("label", "condition", "replicate") %in% colnames(col_data))) {
    stop("'label', 'condition' and/or 'replicate' columns are not present in '",
         deparse(substitute(se)), "'\nRun make_se() or make_se_parse() to obtain the required columns", call. = FALSE)
  }
  if (any(is.na(raw))) {
    warning("Missing values in '", deparse(substitute(se)), "'")
  }

  if (!is.null(control)) {
    if (any(!control %in% unique(col_data$condition))) {
      stop("run test_limma() with a valid control.\nValid controls are: '",
           paste0(unique(col_data$condition), collapse = "', '"), "'", call. = FALSE)
    }
  }

  variables <- terms.formula(design_formula) %>%
    attr(., "variables") %>%
    as.character() %>%
    .[-1]
  if (any(!variables %in% colnames(col_data))) {
    stop("run make_diff() with an appropriate 'design_formula'")
  }
  if (variables[1] != "condition") {
    stop("first factor of 'design_formula' should be 'condition'")
  }

  for (var in variables) {
    temp <- factor(col_data[[var]])
    assign(var, temp)
  }

  design <- model.matrix(design_formula, data = environment())
  colnames(design) <- gsub("condition", "", colnames(design))

  conditions <- as.character(unique(col_data$condition))
  if (type == "all") {
    cntrst <- apply(utils::combn(conditions, 2), 2, paste, collapse = " - ")
    if (!is.null(control)) {
      flip <- grep(paste("^", control, sep = ""), cntrst)
      if (length(flip) >= 1) {
        cntrst[flip] <- cntrst[flip] %>%
          gsub(paste(control, "- ", sep = " "), "", .) %>%
          paste(" - ", control, sep = "")
      }
    }
  } else if (type == "control") {
    if (is.null(control)) stop("Please specify 'control' condition")
    conditions_other_than_control <- conditions[!conditions %in% control]
    cntrst <- c()
    for (i in seq_along(control)) {
      cntrst <- c(cntrst, paste(conditions_other_than_control, control[i], sep = " - "))
    }
  } else if (type == "manual") {
    if (is.null(test)) stop("run test_diff(type = 'manual') with a 'test' argument")
    assertthat::assert_that(is.character(test))
    if (any(!unlist(strsplit(test, "_vs_")) %in% conditions)) {
      stop("run test_diff() with valid contrasts in 'test'",
           ".\nValid contrasts should contain combinations of: '",
           paste0(conditions, collapse = "', '"), "'.", call. = FALSE)
    }
    cntrst <- gsub("_vs_", " - ", test)
  } else if (type == "others") {
    for (i in seq_along(conditions)) {
      design <- cbind(design, ifelse(design[, conditions[i]], 0, 1))
      colnames(design)[ncol(design)] <- paste0("NOT_", conditions[i])
    }
  }

  if (type == "others") {
    message("Tested contrasts: ", paste(paste0(conditions, "_vs_others"), collapse = ", "))
    limma_res <- data.frame()
    for (c in conditions) {
      sub_design <- design[, c(c, paste0("NOT_", c))]
      fit <- limma::lmFit(raw, design = sub_design)
      made_contrasts <- limma::makeContrasts(contrasts = paste0(c, "-", "NOT_", c), levels = sub_design)
      contrast_fit <- limma::contrasts.fit(fit, made_contrasts)
      eB_fit <- limma::eBayes(contrast_fit)
      temp <- limma::topTable(eB_fit, sort.by = "t", adjust.method = "BH", coef = paste0(c, "-", "NOT_", c), number = Inf, confint = TRUE)
      temp <- temp[, c("logFC", "CI.L", "CI.R", "P.Value", "adj.P.Val")]
      colnames(temp) <- c("diff", "CI.L", "CI.R", "p.val", "p.adj")
      colnames(temp) <- paste0(c, "_vs_others_", colnames(temp))
      temp <- tibble::rownames_to_column(temp, "Row.names")
      if (nrow(limma_res) == 0) {
        limma_res <- temp
      } else {
        limma_res <- merge(limma_res, temp, by = "Row.names")
      }
    }
    rowData(se) <- as.data.frame(dplyr::left_join(as.data.frame(rowData(se)), limma_res, by = c("ID" = "Row.names")))
  } else {
    message("Tested contrasts: ", paste(gsub(" - ", "_vs_", cntrst), collapse = ", "))
    fit <- limma::lmFit(raw, design = design)
    made_contrasts <- limma::makeContrasts(contrasts = cntrst, levels = design)
    contrast_fit <- limma::contrasts.fit(fit, made_contrasts)

    if (!type %in% c("others")) {
      if (any(is.na(raw))) {
        for (i in cntrst) {
          covariates <- strsplit(i, " - ") %>% unlist()
          single_contrast <- limma::makeContrasts(contrasts = i, levels = design[, covariates])
          single_contrast_fit <- limma::contrasts.fit(fit[, covariates], single_contrast)
          contrast_fit$coefficients[, i] <- single_contrast_fit$coefficients[, 1]
          contrast_fit$stdev.unscaled[, i] <- single_contrast_fit$stdev.unscaled[, 1]
        }
      }
    }

    eB_fit <- limma::eBayes(contrast_fit)

    retrieve_fun <- function(comp, fit = eB_fit) {
      res <- limma::topTable(fit, sort.by = "t", adjust.method = "BH", coef = comp, number = Inf, confint = TRUE)
      res$comparison <- rep(comp, dim(res)[1])
      res <- tibble::rownames_to_column(res)
      return(res)
    }

    limma_res <- purrr::map_df(cntrst, retrieve_fun)

    table <- limma_res %>%
      dplyr::select(rowname, logFC, CI.L, CI.R, P.Value, adj.P.Val, comparison) %>%
      dplyr::mutate(comparison = gsub(" - ", "_vs_", comparison)) %>%
      tidyr::gather(variable, value, -c(rowname, comparison)) %>%
      dplyr::mutate(variable = dplyr::recode(variable, logFC = "diff", P.Value = "p.val", adj.P.Val = "p.adj")) %>%
      tidyr::unite(temp, comparison, variable) %>%
      tidyr::spread(temp, value)
    rowData(se) <- as.data.frame(dplyr::left_join(as.data.frame(rowData(se)), table, by = c("ID" = "rowname")))
  }
  return(se)
}

# from MonashProteomics/FragPipe-Analyst R/customized.R # 105
# ---- test_diff_customized: Local and tail area-based FDR (fdrtool on t-statistics) ----
test_diff_customized <- function(se, type = c("control", "all", "others", "manual"),
                                 control = NULL, test = NULL,
                                 design_formula = formula(~ 0 + condition)) {
  assertthat::assert_that(inherits(se, "SummarizedExperiment"),
                          is.character(type),
                          class(design_formula) == "formula")
  type <- match.arg(type)

  col_data <- colData(se)
  raw <- assay(se)

  if (any(!c("name", "ID") %in% colnames(rowData(se, use.names = FALSE)))) {
    stop("'name' and/or 'ID' columns are not present in '", deparse(substitute(se)),
         "'\nRun make_unique() and make_se() to obtain the required columns", call. = FALSE)
  }
  if (any(!c("label", "condition", "replicate") %in% colnames(col_data))) {
    stop("'label', 'condition' and/or 'replicate' columns are not present in '",
         deparse(substitute(se)), "'\nRun make_se() or make_se_parse() to obtain the required columns", call. = FALSE)
  }
  if (any(is.na(raw))) {
    warning("Missing values in '", deparse(substitute(se)), "'")
  }

  if (!is.null(control)) {
    if (any(!control %in% unique(col_data$condition))) {
      stop("run test_diff() with a valid control.\nValid controls are: '",
           paste0(unique(col_data$condition), collapse = "', '"), "'", call. = FALSE)
    }
  }

  variables <- terms.formula(design_formula) %>%
    attr(., "variables") %>%
    as.character() %>%
    .[-1]
  if (any(!variables %in% colnames(col_data))) {
    stop("run make_diff() with an appropriate 'design_formula'")
  }

  for (var in variables) {
    temp <- factor(col_data[[var]])
    assign(var, temp)
  }

  design <- model.matrix(design_formula, data = environment())
  colnames(design) <- gsub("condition", "", colnames(design))

  conditions <- as.character(unique(col_data$condition))
  if (type == "all") {
    cntrst <- apply(utils::combn(conditions, 2), 2, paste, collapse = " - ")
    if (!is.null(control)) {
      flip <- grep(paste("^", control, sep = ""), cntrst)
      if (length(flip) >= 1) {
        cntrst[flip] <- cntrst[flip] %>%
          gsub(paste(control, "- ", sep = " "), "", .) %>%
          paste(" - ", control, sep = "")
      }
    }
  } else if (type == "control") {
    if (is.null(control)) stop("run test_diff(type = 'control') with a 'control' argument")
    cntrst <- paste(conditions[!conditions %in% control], control, sep = " - ")
  } else if (type == "manual") {
    if (is.null(test)) stop("run test_diff(type = 'manual') with a 'test' argument")
    assertthat::assert_that(is.character(test))
    cntrst <- gsub("_vs_", " - ", test)
  } else if (type == "others") {
    for (i in seq_along(conditions)) {
      design <- cbind(design, ifelse(design[, conditions[i]], 0, 1))
      colnames(design)[ncol(design)] <- paste0("NOT_", conditions[i])
    }
  }

  if (type == "others") {
    message("Tested contrasts: ", paste(paste0(conditions, "_vs_others"), collapse = ", "))
    limma_res <- data.frame()
    for (c in conditions) {
      sub_design <- design[, c(c, paste0("NOT_", c))]
      fit <- limma::lmFit(raw, design = sub_design)
      made_contrasts <- limma::makeContrasts(contrasts = paste0(c, "-", "NOT_", c), levels = sub_design)
      contrast_fit <- limma::contrasts.fit(fit, made_contrasts)
      eB_fit <- limma::eBayes(contrast_fit)
      temp <- limma::topTable(eB_fit, sort.by = "t", coef = paste0(c, "-", "NOT_", c), number = Inf, confint = TRUE)
      temp <- temp[!is.na(temp$t), ]
      fdr_res <- fdrtool::fdrtool(temp$t, plot = FALSE, verbose = FALSE)
      temp$qval <- fdr_res$qval
      temp <- temp[, c("logFC", "CI.L", "CI.R", "P.Value", "qval")]
      colnames(temp) <- c("diff", "CI.L", "CI.R", "p.val", "p.adj")
      colnames(temp) <- paste0(c, "_vs_others_", colnames(temp))
      temp <- tibble::rownames_to_column(temp, "Row.names")
      if (nrow(limma_res) == 0) {
        limma_res <- temp
      } else {
        limma_res <- merge(limma_res, temp, by = "Row.names")
      }
    }
    rowData(se) <- as.data.frame(dplyr::left_join(as.data.frame(rowData(se)), limma_res, by = c("ID" = "Row.names")))
  } else {
    message("Tested contrasts: ", paste(gsub(" - ", "_vs_", cntrst), collapse = ", "))
    fit <- limma::lmFit(raw, design = design)
    made_contrasts <- limma::makeContrasts(contrasts = cntrst, levels = design)
    contrast_fit <- limma::contrasts.fit(fit, made_contrasts)

    if (type != "manual" && any(is.na(raw))) {
      for (i in cntrst) {
        covariates <- strsplit(i, " - ") %>% unlist()
        single_contrast <- limma::makeContrasts(contrasts = i, levels = design[, covariates])
        single_contrast_fit <- limma::contrasts.fit(fit[, covariates], single_contrast)
        contrast_fit$coefficients[, i] <- single_contrast_fit$coefficients[, 1]
        contrast_fit$stdev.unscaled[, i] <- single_contrast_fit$stdev.unscaled[, 1]
      }
    }
    eB_fit <- limma::eBayes(contrast_fit)

    retrieve_fun <- function(comp, fit = eB_fit) {
      res <- limma::topTable(fit, sort.by = "t", coef = comp, number = Inf, confint = TRUE)
      res <- res[!is.na(res$t), ]
      fdr_res <- fdrtool::fdrtool(res$t, plot = FALSE, verbose = FALSE)
      res$qval <- fdr_res$qval
      res$comparison <- rep(comp, dim(res)[1])
      res <- tibble::rownames_to_column(res)
      return(res)
    }

    limma_res <- purrr::map_df(cntrst, retrieve_fun)

    table <- limma_res %>%
      dplyr::select(rowname, logFC, CI.L, CI.R, P.Value, qval, comparison) %>%
      dplyr::mutate(comparison = gsub(" - ", "_vs_", comparison)) %>%
      tidyr::gather(variable, value, -c(rowname, comparison)) %>%
      dplyr::mutate(variable = dplyr::recode(variable, logFC = "diff", P.Value = "p.val", qval = "p.adj")) %>%
      tidyr::unite(temp, comparison, variable) %>%
      tidyr::spread(temp, value)
    rowData(se) <- as.data.frame(dplyr::left_join(as.data.frame(rowData(se)), table, by = c("ID" = "rowname")))
  }
  return(se)
}

# from MonashProteomics/FragPipe-Analyst R/customized.R # 324
# ---- add_rejections_customized: used by BOTH (marks significant from p.adj and diff) ----
add_rejections_customized <- function(diff, alpha = 0.05, lfc = 1) {
  if (is.integer(alpha)) alpha <- as.numeric(alpha)
  if (is.integer(lfc)) lfc <- as.numeric(lfc)
  assertthat::assert_that(inherits(diff, "SummarizedExperiment"),
                          is.numeric(alpha), length(alpha) == 1,
                          is.numeric(lfc), length(lfc) == 1)

  row_data <- rowData(diff, use.names = FALSE) %>% as.data.frame()
  if (any(!c("name", "ID") %in% colnames(row_data))) {
    stop("'name' and/or 'ID' columns are not present in '", deparse(substitute(diff)),
         "'\nRun make_unique() and make_se() to obtain the required columns", call. = FALSE)
  }
  if (length(grep("_p.adj|_diff", colnames(row_data))) < 1) {
    stop("'[contrast]_diff' and/or '[contrast]_p.adj' columns are not present in '",
         deparse(substitute(diff)), "'\nRun test_diff() to obtain the required columns", call. = FALSE)
  }

  cols_p <- grep("_p.adj", colnames(row_data))
  cols_diff <- grep("_diff", colnames(row_data))

  if (length(cols_p) == 1) {
    rowData(diff)$significant <-
      row_data[, cols_p] <= alpha & abs(row_data[, cols_diff]) >= lfc
    rowData(diff)$contrast_significant <- rowData(diff, use.names = FALSE)$significant
    colnames(rowData(diff))[ncol(rowData(diff, use.names = FALSE))] <-
      gsub("p.adj", "significant", colnames(row_data)[cols_p])
  }
  if (length(cols_p) > 1) {
    p_reject <- row_data[, cols_p] <= alpha
    p_reject[is.na(p_reject)] <- FALSE
    diff_reject <- abs(row_data[, cols_diff]) >= lfc
    diff_reject[is.na(diff_reject)] <- FALSE
    sign_df <- p_reject & diff_reject
    sign_df <- cbind(sign_df, significant = apply(sign_df, 1, function(x) any(x)))
    colnames(sign_df) <- gsub("_p.adj", "_significant", colnames(sign_df))
    sign_df <- cbind(ID = row_data$ID, as.data.frame(sign_df))
    rowData(diff) <- as.data.frame(dplyr::left_join(as.data.frame(rowData(diff)), sign_df, by = c("ID" = "ID")))
  }
  return(diff)
}


# from MonashProteomics/FragPipe-Analyst R/functions.R # 81 (get_cluster_heatmap)
# Modified: label->sample_name col mapping for heatmap when colnames=label.
# ========== DE heatmap (get_cluster_heatmap_customized) ==========
get_cluster_heatmap_customized <- function(dep, type = c("contrast", "centered"),
                                            kmeans = FALSE, k = 6, col_limit = 6, indicate = NULL,
                                            alpha = 0.01, lfc = 1,
                                            clustering_distance = c("euclidean", "maximum", "manhattan",
                                                "canberra", "binary", "minkowski", "pearson", "spearman", "kendall", "gower"),
                                            row_font_size = 6, col_font_size = 10, plot = TRUE, ...) {
  if (is.integer(k)) k <- as.numeric(k)
  if (is.integer(col_limit)) col_limit <- as.numeric(col_limit)
  if (is.integer(row_font_size)) row_font_size <- as.numeric(row_font_size)
  if (is.integer(col_font_size)) col_font_size <- as.numeric(col_font_size)
  assertthat::assert_that(inherits(dep, "SummarizedExperiment"), is.character(type),
                          is.logical(kmeans), is.numeric(k), length(k) == 1,
                          is.numeric(col_limit), length(col_limit) == 1,
                          is.numeric(row_font_size), length(row_font_size) == 1,
                          is.numeric(col_font_size), length(col_font_size) == 1,
                          is.logical(plot), length(plot) == 1)
  type <- match.arg(type)
  clustering_distance <- match.arg(clustering_distance)

  row_data <- rowData(dep)
  col_data <- as.data.frame(colData(dep))
  if (any(!c("label", "condition", "replicate") %in% colnames(col_data))) {
    stop("'label', 'condition' and/or 'replicate' columns are not present in colData.", call. = FALSE)
  }
  if (length(grep("_diff", colnames(row_data))) < 1) {
    stop("'[contrast]_diff' columns are not present. Run test_diff() first.", call. = FALSE)
  }
  if (!"significant" %in% colnames(row_data)) {
    stop("'significant' column is not present. Run add_rejections() first.", call. = FALSE)
  }

  ha1 <- NULL
  if (!is.null(indicate) && type == "centered") {
    ha1 <- FragPipeAnalystR:::get_annotation(dep, indicate)
  }

  conditions <- gsub("_diff", "", colnames(row_data)[grepl("_diff", colnames(row_data))])
  cols_p <- paste0(conditions, "_p.adj")
  cols_lfc <- paste0(conditions, "_diff")
  p <- as.matrix(row_data[, cols_p]) <= alpha
  lfc_mat <- abs(as.matrix(row_data[, cols_lfc])) >= lfc
  p[is.na(p)] <- FALSE
  lfc_mat[is.na(p)] <- FALSE
  combined_rejections <- p & lfc_mat
  filtered <- dep[apply(combined_rejections, 1, any), ]

  if (nrow(filtered) == 0) {
    return(ggplot2::ggplot() +
      ggplot2::annotate("text", x = 4, y = 25, size = 8, label = "No differential expressed genes available for the heatmap") +
      ggplot2::theme_void())
  }

  if (any(is.na(assay(filtered)))) {
    clustering_distance <- "gower"
    obs_NA <- TRUE
  } else {
    obs_NA <- FALSE
  }

  if (type == "centered") {
    rowData(filtered)$mean <- rowMeans(assay(filtered), na.rm = TRUE)
    df <- assay(filtered) - rowData(filtered)$mean
  } else {
    df <- as.data.frame(rowData(filtered)) %>%
      tibble::column_to_rownames(var = "name") %>%
      dplyr::select(dplyr::ends_with("_diff"))
    colnames(df) <- gsub("_vs_", " vs ", gsub("_diff", "", colnames(df)))
  }

  if (kmeans && obs_NA) kmeans <- FALSE
  if (kmeans && !obs_NA) {
    set.seed(1)
    df_kmeans <- kmeans(df, k)
  }

  col_clust <- ncol(df) > 1
  row_clust <- nrow(df) > 1
  if (clustering_distance == "gower") {
    clustering_distance <- function(x) {
      d <- cluster::daisy(x, metric = "gower")
      d[is.na(d)] <- max(d, na.rm = TRUE)
      d
    }
  }

  legend <- ifelse(type == "contrast", "log2 Fold change", "log2 Centered intensity")

  temp <- as.data.frame(colData(filtered))
  rownames(temp) <- temp$label
  new_names <- temp[colnames(df), "sample_name"]
  if (!any(is.na(new_names))) colnames(df) <- new_names

  ht1 <- ComplexHeatmap::Heatmap(df,
    col = circlize::colorRamp2(seq(-col_limit, col_limit, col_limit/5), rev(RColorBrewer::brewer.pal(11, "RdBu"))),
    split = if (kmeans) df_kmeans$cluster else NULL,
    cluster_rows = col_clust,
    cluster_columns = row_clust,
    row_names_side = "left",
    column_names_side = "top",
    clustering_distance_rows = clustering_distance,
    clustering_distance_columns = clustering_distance,
    heatmap_legend_param = list(color_bar = "continuous", legend_direction = "horizontal",
      legend_width = grid::unit(5, "cm"), title_position = "lefttop"),
    name = legend,
    row_names_gp = grid::gpar(fontsize = row_font_size),
    column_names_gp = grid::gpar(fontsize = col_font_size),
    top_annotation = ha1,
    ...)
  p <- ComplexHeatmap::draw(ht1, heatmap_legend_side = "top")
  list(ht1, ComplexHeatmap::row_order(ht1))
}


# from MonashProteomics/FragPipe-Analyst R/customized.R # 1621
plot_cor_customized <- function(dep, significant = FALSE, lower = -1, upper = 1,
                                pal = "PRGn", pal_rev = FALSE, indicate = NULL,
                                font_size = 12, plot = FALSE, ...) {
  assertthat::assert_that(inherits(dep, "SummarizedExperiment"),
                          is.logical(significant), length(significant) == 1,
                          is.numeric(lower), length(lower) == 1,
                          is.numeric(upper), length(upper) == 1,
                          is.character(pal), length(pal) == 1,
                          is.logical(pal_rev), length(pal_rev) == 1,
                          is.numeric(font_size), length(font_size) == 1,
                          is.logical(plot), length(plot) == 1)

  if (!(lower >= -1 & upper >= -1 & lower <= 1 & upper <= 1)) {
    stop("'lower' and/or 'upper' arguments are not valid. Use values between -1 and 1.", call. = FALSE)
  }

  pals <- RColorBrewer::brewer.pal.info %>%
    tibble::rownames_to_column() %>%
    dplyr::filter(category != "qual")
  if (!pal %in% pals$rowname) {
    stop("'", pal, "' is not a valid color panel. Try: ", paste(pals$rowname, collapse = ", "), call. = FALSE)
  }

  ha1 <- NULL
  if (!is.null(indicate)) {
    assertthat::assert_that(is.character(indicate))
    col_data <- as.data.frame(colData(dep))
    if (any(!indicate %in% colnames(col_data))) {
      stop("'", paste0(indicate, collapse = "' and/or '"), "' not in colData. Valid: ",
           paste(colnames(col_data), collapse = ", "), call. = FALSE)
    }
    anno <- as.data.frame(colData(dep)) %>% dplyr::select(dplyr::all_of(indicate))
    names <- colnames(anno)
    anno_col <- vector(mode = "list", length = length(names))
    names(anno_col) <- names
    for (i in names) {
      var <- sort(unique(anno[[i]]))
      nv <- length(var)
      cols <- if (nv == 1) c("black") else
              if (nv == 2) c("orangered", "cornflowerblue") else
              if (nv <= 6) RColorBrewer::brewer.pal(max(3, nv), "Pastel1")[1:nv] else
              if (nv <= 12) RColorBrewer::brewer.pal(nv, "Set3") else
              colorRampPalette(RColorBrewer::brewer.pal(12, "Set3"))(nv)
      names(cols) <- var
      anno_col[[i]] <- cols
    }
    ha1 <- ComplexHeatmap::HeatmapAnnotation(df = anno, col = anno_col, show_annotation_name = TRUE)
  }

  if (significant) {
    if (!"significant" %in% colnames(rowData(dep, use.names = FALSE))) {
      stop("'significant' column not present. Run add_rejections() first.", call. = FALSE)
    }
    dep <- dep[tidyr::replace_na(rowData(dep, use.names = FALSE)$significant, FALSE), ]
  }

  data <- assay(dep)
  temp <- as.data.frame(colData(dep))
  if ("label" %in% colnames(temp) && "sample_name" %in% colnames(temp)) {
    rownames(temp) <- temp$label
    new_names <- temp[colnames(data), "sample_name"]
    if (!any(is.na(new_names))) colnames(data) <- new_names
  }

  cor_mat <- cor(data, use = "complete.obs")
  lower <- min(cor_mat)
  upper <- max(cor_mat)

  ht1 <- ComplexHeatmap::Heatmap(cor_mat,
    col = circlize::colorRamp2(
      seq(lower, upper, (upper - lower) / 7),
      if (pal_rev) rev(RColorBrewer::brewer.pal(8, pal)) else RColorBrewer::brewer.pal(8, pal)
    ),
    heatmap_legend_param = list(
      color_bar = "continuous",
      legend_direction = "horizontal",
      legend_width = grid::unit(5, "cm"),
      title_position = "topcenter"
    ),
    name = "Pearson correlation",
    column_names_gp = grid::gpar(fontsize = font_size),
    row_names_gp = grid::gpar(fontsize = font_size),
    top_annotation = ha1,
    ...
  )
  if (plot) {
    ComplexHeatmap::draw(ht1, heatmap_legend_side = "top")
  }
  return(ht1)
}


# from arnesmits/DEP R/gg_theme.R # 21
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
  return(theme)
}

# plot_feature_numbers_monash: based on DEP plot_numbers (plot_functions_frequencies.R)
# Aligns with FragPipe-Analyst GUI: join by label (assay colnames = expdesign$label), fill by condition.
plot_feature_numbers_monash <- function(se, fill = "condition") {
  assertthat::assert_that(inherits(se, "SummarizedExperiment"))
  df <- assay(se) %>%
    data.frame(check.names = FALSE) %>%
    tibble::rownames_to_column() %>%
    tidyr::gather(ID, bin, -rowname) %>%
    dplyr::mutate(bin = ifelse(is.na(bin), 0, 1))
  stat <- df %>%
    dplyr::group_by(ID) %>%
    dplyr::summarize(n = dplyr::n(), sum = sum(bin))
  cd <- as.data.frame(colData(se))
  id_col <- if ("label" %in% colnames(cd)) "label" else "sample_name"
  stat <- dplyr::left_join(stat, cd, by = c("ID" = id_col))
  feature <- if (!is.null(metadata(se)$level)) {
    switch(metadata(se)$level, protein = "Proteins", peptide = "Peptides", gene = "Peptides", site = "Sites", "Features")
  } else "Features"
  p <- ggplot2::ggplot(stat, ggplot2::aes(x = ID, y = sum, fill = .data[[fill]])) +
    ggplot2::geom_col() +
    ggplot2::geom_hline(yintercept = unique(stat$n), linetype = "dashed") +
    ggplot2::labs(title = paste0("Number of ", feature, " per Sample (Total: ", nrow(se), ")"),
                  x = "", y = paste0("Number of ", feature)) +
    theme_DEP1() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1))
  return(p)
}

# from FragPipe-Analyst R/customized.R plot_coverage_customized (lines 356-393)
# Shows only sample counts that exist in data (e.g. 3-12 after global_filter; no 0,1,2 if filtered out)
plot_coverage_customized <- function(se, plot = TRUE) {
  assertthat::assert_that(inherits(se, "SummarizedExperiment"),
                          is.logical(plot), length(plot) == 1)
  df <- assay(se) %>%
    data.frame() %>%
    tibble::rownames_to_column() %>%
    tidyr::gather(ID, bin, -rowname) %>%
    dplyr::mutate(bin = ifelse(is.na(bin), 0, 1))
  stat <- df %>%
    dplyr::group_by(rowname) %>%
    dplyr::summarize(sum = sum(bin))
  table <- table(stat$sum) %>%
    data.frame()
  p <- ggplot2::ggplot(table, ggplot2::aes(x = "all", y = Freq, fill = Var1)) +
    ggplot2::geom_col(col = "white") +
    ggplot2::scale_fill_grey(start = 0.8, end = 0.2) +
    ggplot2::labs(title = "Feature coverage",
                  x = "",
                  y = "Number of features",
                  fill = "Samples") +
    theme_DEP1()
  if (plot) {
    return(p)
  } else {
    df_out <- as.data.frame(table)
    colnames(df_out) <- c("samples", "features")
    return(df_out)
  }
}


# ========== Absence/Presence plots - SOURCE from MonashProteomics/FragPipe-Analyst ==========
# server.R: data_attendance, venn_plot_input, upset_plot_input
# R/customized.R: plot_Jaccard, plot_density

# data_attendance - from Monash server.R lines 1825-1922
# Builds occurrence matrix (#Occurences per condition) for Venn/UpSet.
data_attendance_monash <- function(se, exp = "LFQ", level = "protein") {
  assertthat::assert_that(inherits(se, "SummarizedExperiment"))
  df <- as.data.frame(assay(se), check.names = FALSE)
  col_data <- as.data.frame(colData(se))
  sample_cols <- colnames(df)
  conditions <- unique(col_data$condition)

  if (exp == "LFQ") {
    df$Gene <- as.character(rowData(se)$Gene)
    rd <- as.data.frame(rowData(se))
    df$Protein <- as.character(if ("Protein ID" %in% colnames(rd)) rd[["Protein ID"]] else if ("Protein" %in% colnames(rd)) rd$Protein else rd$ID)
    if (any(df$Gene == "", na.rm = TRUE)) df$Gene[df$Gene == ""] <- "NoGeneNameAvailable"
    df <- df[rowSums(!is.na(df[, sample_cols])) != 0, ]
    for (i in seq_along(conditions)) {
      cond <- conditions[i]
      temp <- col_data[col_data$condition == cond, , drop = FALSE]
      sel_cols <- intersect(rownames(temp), colnames(df))
      df[[paste0("#Occurences_", cond)]] <- rowSums(!is.na(df[, sel_cols, drop = FALSE]))
    }
    df <- dplyr::relocate(df, Protein, Gene, .before = 1)
  } else if (exp == "DIA" && level == "protein") {
    df$Gene <- as.character(rowData(se)$Genes)
    df$Protein <- if ("Protein.Ids" %in% colnames(rowData(se))) as.character(rowData(se)$Protein.Ids) else as.character(rowData(se)$Protein.Group)
    if (any(df$Gene == "", na.rm = TRUE)) df$Gene[df$Gene == ""] <- "NoGeneNameAvailable"
    df <- df[rowSums(!is.na(df[, sample_cols])) != 0, ]
    for (i in seq_along(conditions)) {
      cond <- conditions[i]
      temp <- col_data[col_data$condition == cond, ]
      sel_cols <- intersect(temp$label, colnames(df))
      df[[paste0("#Occurences_", cond)]] <- rowSums(!is.na(df[, sel_cols, drop = FALSE]))
    }
    df <- dplyr::relocate(df, Protein, Gene, .before = 1)
  } else {
    df$Gene <- as.character(rowData(se)$Genes)
    if (any(df$Gene == "", na.rm = TRUE)) df$Gene[df$Gene == ""] <- "NoGeneNameAvailable"
    df <- df[rowSums(!is.na(df[, sample_cols])) != 0, ]
    for (i in seq_along(conditions)) {
      cond <- conditions[i]
      temp <- col_data[col_data$condition == cond, ]
      sel_cols <- intersect(temp$label, colnames(df))
      df[[paste0("#Occurences_", cond)]] <- rowSums(!is.na(df[, sel_cols, drop = FALSE]))
    }
    df <- dplyr::relocate(df, Gene, .before = 1)
  }
  rownames(df) <- NULL
  df
}

# Venn plot - from Monash server.R venn_plot_input, uses ggVennDiagram
plot_venn_monash <- function(df, cond1, cond2, cond3 = NULL) {
  if (!requireNamespace("ggVennDiagram", quietly = TRUE)) {
    warning("ggVennDiagram not installed. Skipping Venn. Install with: install.packages('ggVennDiagram')")
    return(NULL)
  }
  occ1 <- paste0("#Occurences_", cond1)
  occ2 <- paste0("#Occurences_", cond2)
  if (!occ1 %in% colnames(df) || !occ2 %in% colnames(df)) return(NULL)
  set1 <- df[df[[occ1]] != 0, "Gene"]
  set2 <- df[df[[occ2]] != 0, "Gene"]
  x <- list(set1, set2)
  names(x) <- c(cond1, cond2)
  if (!is.null(cond3) && cond3 != "NONE") {
    occ3 <- paste0("#Occurences_", cond3)
    if (occ3 %in% colnames(df)) {
      set3 <- df[df[[occ3]] != 0, "Gene"]
      x <- list(set1, set2, set3)
      names(x) <- c(cond1, cond2, cond3)
    }
  }
  ggVennDiagram::ggVennDiagram(x, label_alpha = 0) +
    ggplot2::scale_fill_gradient(low = "#F4FAFE", high = "#4981BF") +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = 0.3)) +
    ggplot2::coord_flip()
}

# UpSet plot - from Monash server.R upset_plot_input
plot_upset_monash <- function(df) {
  if (!requireNamespace("UpSetR", quietly = TRUE)) {
    warning("UpSetR not installed. Skipping UpSet.")
    return(invisible(NULL))
  }
  df <- df[, grep("#Occurences", colnames(df)), drop = FALSE]
  df <- ifelse(df != 0, 1, 0)
  df <- data.frame(df)
  colnames(df) <- gsub("X.Occurences_", "", colnames(df))
  if (sum(colSums(df) != 0) <= 1) return(invisible(NULL))
  p <- UpSetR::upset(df, nsets = ncol(df), mb.ratio = c(0.6, 0.4),
                     text.scale = 1.5, point.size = 3, order.by = "freq", decreasing = TRUE,
                     nintersects = NA, mainbar.y.label = "#Features in intersection",
                     sets.x.label = "#Features", set_size.scale_max = nrow(df) + 1000, set_size.show = TRUE)
  print(p)
  invisible(NULL)
}

# plot_Jaccard - from Monash R/customized.R (sample-level Jaccard via vegdist)
# Uses top_annotation for condition and sample_name for labels
plot_Jaccard_monash <- function(dep, plot = TRUE, exp = "LFQ", indicate = "condition") {
  assertthat::assert_that(inherits(dep, "SummarizedExperiment"))
  if (!requireNamespace("vegan", quietly = TRUE)) {
    warning("vegan not installed. Skipping Jaccard. Install with: install.packages('vegan')")
    return(NULL)
  }
  data <- assay(dep)
  cd <- as.data.frame(colData(dep))
  colnames(data) <- cd[match(colnames(data), rownames(cd)), "sample_name"]
  cor_mat <- 1 - as.matrix(vegan::vegdist(t(data), method = "jaccard", na.rm = TRUE))
  lower <- min(cor_mat)
  upper <- max(cor_mat)
  ha1 <- NULL
  if (!is.null(indicate) && indicate %in% colnames(cd)) {
    anno <- cd[, indicate, drop = FALSE]
    rownames(anno) <- colnames(cor_mat)
    var <- unique(anno[[1]])
    if (length(var) == 1) {
      cols <- c("black")
    } else if (length(var) == 2) {
      cols <- c("orangered", "cornflowerblue")
    } else if (length(var) < 7) {
      cols <- RColorBrewer::brewer.pal(max(3, length(var)), "Pastel1")[seq_len(length(var))]
    } else {
      cols <- RColorBrewer::brewer.pal(length(var), "Set3")
    }
    names(cols) <- var
    ha1 <- ComplexHeatmap::HeatmapAnnotation(df = anno, col = setNames(list(cols), indicate),
      show_annotation_name = TRUE)
  }
  ht1 <- ComplexHeatmap::Heatmap(cor_mat,
    col = circlize::colorRamp2(c(lower, (upper + lower) / 2, upper), c("blue", "lightyellow", "red")),
    heatmap_legend_param = list(color_bar = "continuous", legend_direction = "horizontal",
      legend_width = grid::unit(5, "cm"), title_position = "topcenter"),
    name = "Jaccard similarity", column_names_gp = grid::gpar(fontsize = 12), row_names_gp = grid::gpar(fontsize = 12),
    top_annotation = ha1
  )
  if (plot) ComplexHeatmap::draw(ht1, heatmap_legend_side = "top") else as.data.frame(cor_mat)
}

# plot_density - from Monash R/customized.R (list of SEs: original, filtered, imputed)
# Overlaid density curves colored by condition (Control, Ubi, etc.)
plot_density_monash <- function(ses) {
  gather_join <- function(se) {
    cd <- as.data.frame(colData(se))
    cd$..colkey.. <- rownames(cd)
    assay(se) %>%
      data.frame(check.names = FALSE) %>%
      tidyr::gather(ID, val, dplyr::everything()) %>%
      dplyr::left_join(cd, by = c("ID" = "..colkey.."))
  }
  df <- purrr::map_df(ses, gather_join, .id = "var") %>%
    dplyr::mutate(var = factor(var, levels = names(ses)))
  ggplot2::ggplot(df, ggplot2::aes(val, col = condition)) +
    ggplot2::geom_density(na.rm = TRUE) +
    ggplot2::facet_wrap(~var, ncol = 1, strip.position = "top") +
    ggplot2::labs(x = expression(log[2] ~ "Intensity"), y = "Density") +
    theme_DEP1()
}
