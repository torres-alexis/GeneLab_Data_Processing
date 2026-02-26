# fp_de_helper.R
# Differential expression: limma + BH or fdrtool. FragPipe-Analyst.
# Dependencies: limma, fdrtool (only for test_diff_customized), assertthat, dplyr, tidyr, purrr, tibble.

# test_limma_customized: Benjamini-Hochberg FDR via limma topTable.
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
      temp <- temp[, c("logFC", "CI.L", "CI.R", "P.Value", "adj.P.Val", "t")]
      colnames(temp) <- c("diff", "CI.L", "CI.R", "p.val", "p.adj", "t")
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
      dplyr::select(rowname, logFC, CI.L, CI.R, P.Value, adj.P.Val, t, comparison) %>%
      dplyr::mutate(comparison = gsub(" - ", "_vs_", comparison)) %>%
      tidyr::gather(variable, value, -c(rowname, comparison)) %>%
      dplyr::mutate(variable = dplyr::recode(variable, logFC = "diff", P.Value = "p.val", adj.P.Val = "p.adj", t = "t")) %>%
      tidyr::unite(temp, comparison, variable) %>%
      tidyr::spread(temp, value)
    rowData(se) <- as.data.frame(dplyr::left_join(as.data.frame(rowData(se)), table, by = c("ID" = "rowname")))
  }
  return(se)
}

# test_diff_customized: Local and tail area-based FDR (fdrtool on t-statistics). Requires fdrtool.
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
      temp <- temp[, c("logFC", "CI.L", "CI.R", "P.Value", "qval", "t")]
      colnames(temp) <- c("diff", "CI.L", "CI.R", "p.val", "p.adj", "t")
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
      dplyr::select(rowname, logFC, CI.L, CI.R, P.Value, qval, t, comparison) %>%
      dplyr::mutate(comparison = gsub(" - ", "_vs_", comparison)) %>%
      tidyr::gather(variable, value, -c(rowname, comparison)) %>%
      dplyr::mutate(variable = dplyr::recode(variable, logFC = "diff", P.Value = "p.val", qval = "p.adj", t = "t")) %>%
      tidyr::unite(temp, comparison, variable) %>%
      tidyr::spread(temp, value)
    rowData(se) <- as.data.frame(dplyr::left_join(as.data.frame(rowData(se)), table, by = c("ID" = "rowname")))
  }
  return(se)
}

# add_groupwise_stats: group means and stdevs from assay by condition (for DE_results).
# Returns data frame with Mean_<cond>, Stdev_<cond> columns; rownames = rownames(assay(se)).
add_groupwise_stats <- function(se) {
  assertthat::assert_that(inherits(se, "SummarizedExperiment"))
  cd <- as.data.frame(colData(se))
  if (!"condition" %in% colnames(cd)) return(NULL)
  assay_mat <- assay(se)
  conds <- cd$condition[match(colnames(assay_mat), rownames(cd))]
  if (any(is.na(conds))) return(NULL)
  ucond <- unique(conds)
  group_means <- matrix(NA_real_, nrow = nrow(assay_mat), ncol = length(ucond))
  group_stdev <- matrix(NA_real_, nrow = nrow(assay_mat), ncol = length(ucond))
  for (i in seq_along(ucond)) {
    idx <- which(conds == ucond[i])
    group_means[, i] <- rowMeans(assay_mat[, idx, drop = FALSE], na.rm = TRUE)
    if (length(idx) > 1L) {
      group_stdev[, i] <- matrixStats::rowSds(assay_mat[, idx, drop = FALSE], na.rm = TRUE)
    }
  }
  colnames(group_means) <- paste0("Mean_", make.names(ucond))
  colnames(group_stdev) <- paste0("Stdev_", make.names(ucond))
  out <- as.data.frame(cbind(group_means, group_stdev))
  rownames(out) <- rownames(assay_mat)
  out
}

# get_de_results_extended: full DE table (rowData + assay) with CI.L, CI.R, diff, p.val, p.adj, t (Stat), significant per contrast.
# Extended format (full rowData + assay).
# add_group_stats: if TRUE, append All.mean, All.stdev, Mean_<condition>, Stdev_<condition> (meeting spec).
get_de_results_extended <- function(se, add_group_stats = TRUE) {
  assertthat::assert_that(inherits(se, "SummarizedExperiment"))
  rd <- as.data.frame(rowData(se, use.names = FALSE))
  asy <- as.data.frame(assay(se))
  de_df <- cbind(rd, asy)
  # All.mean and All.stdev: mean and sd across all samples (RNA-seq parity)
  assay_mat <- as.matrix(assay(se))
  all_mean <- rowMeans(assay_mat, na.rm = TRUE)
  all_stdev <- matrixStats::rowSds(assay_mat, na.rm = TRUE)
  all_stats <- data.frame(All.mean = all_mean, All.stdev = all_stdev)
  all_stats_aligned <- all_stats[match(rd$ID, rownames(assay_mat)), , drop = FALSE]
  rownames(all_stats_aligned) <- NULL
  de_df <- cbind(de_df, all_stats_aligned)
  if (add_group_stats) {
    gs <- add_groupwise_stats(se)
    if (!is.null(gs)) {
      gs_aligned <- gs[match(rd$ID, rownames(gs)), , drop = FALSE]
      rownames(gs_aligned) <- NULL
      de_df <- cbind(de_df, gs_aligned)
    }
  }
  de_df
}

# add_sample_suffix_to_export: rename sample columns in df for export (FragPipe parity).
# LFQ: .Intensity, .MaxLFQ.Intensity, .Spectral.Count. TMT: suffix from metadata or channel.
# sample_cols = colnames(assay(se)), suffix = e.g. "Intensity" for LFQ Intensity mode.
add_sample_suffix_to_export <- function(df, sample_cols, suffix) {
  if (is.null(suffix) || !nzchar(suffix)) return(df)
  idx <- colnames(df) %in% sample_cols
  if (!any(idx)) return(df)
  colnames(df)[idx] <- paste0(colnames(df)[idx], ".", suffix)
  df
}

# sanitized_to_raw_map: build map from condition (design matrix) -> condition_raw (display).
# Uses (condition, condition_raw) pairs from colData so lookup matches contrast names exactly.
sanitized_to_raw_map <- function(se) {
  if (!inherits(se, "SummarizedExperiment") || !"condition" %in% colnames(colData(se))) return(NULL)
  cd <- as.data.frame(colData(se))
  if ("condition_raw" %in% colnames(cd) && all(nzchar(trimws(cd$condition_raw)))) {
    # Direct map: condition -> condition_raw from colData (matches contrast names exactly)
    pairs <- unique(cd[, c("condition", "condition_raw")])
    m <- setNames(as.character(pairs$condition_raw), as.character(pairs$condition))
  } else {
    m <- setNames(as.character(unique(cd$condition)), as.character(unique(cd$condition)))
  }
  if (length(m) == 0) return(NULL)
  function(sanitized) {
    if (is.na(sanitized) || !nzchar(sanitized)) return(sanitized)
    if (!is.na(m[sanitized])) return(m[sanitized])
    sanitized
  }
}

# contrast_to_display: convert sanitized contrast "A_vs_B" to display string "rawA vs rawB".
contrast_to_display <- function(contrast, se) {
  f <- sanitized_to_raw_map(se)
  if (is.null(f)) return(contrast)
  if (grepl("_vs_others$", contrast)) {
    a <- sub("_vs_others$", "", contrast)
    return(paste0(f(a), " vs Others"))
  }
  if (grepl("_vs_", contrast)) {
    parts <- strsplit(contrast, "_vs_")[[1]]
    if (length(parts) == 2) return(paste0(f(parts[1]), " vs ", f(parts[2])))
  }
  contrast
}

# apply_rnaseq_headers: rename DE and group columns to RNA-seq format for export.
# Contrast: [A_vs_B]_diff -> Log2fc_(A)v(B), _t -> Stat_, _p.val -> P.value_, _p.adj -> Adj.p.value_
# Group: Mean_X -> Group.Mean_(X), Stdev_X -> Group.Stdev_(X) (X = original condition name)
# Requires se (SummarizedExperiment) for condition names.
apply_rnaseq_headers <- function(de_df, se) {
  assertthat::assert_that(inherits(se, "SummarizedExperiment"))
  conds <- as.character(unique(colData(se)$condition))
  if (length(conds) == 0) return(de_df)
  orig_name <- sanitized_to_raw_map(se)
  if (is.null(orig_name)) orig_name <- function(x) x

  new_names <- colnames(de_df)
  # DE contrast columns: A_vs_B or A_vs_others
  for (suffix in c("_diff", "_t", "_p.val", "_p.adj", "_CI.L", "_CI.R", "_significant")) {
    cols <- grep(paste0(suffix, "$"), colnames(de_df), value = TRUE)
    for (c in cols) {
      base <- sub(suffix, "", c)
      if (grepl("_vs_others$", base)) {
        a <- sub("_vs_others$", "", base)
        rnaseq_contrast <- paste0("(", orig_name(a), ")v(Others)")
      } else if (grepl("_vs_", base)) {
        parts <- strsplit(base, "_vs_")[[1]]
        if (length(parts) == 2) {
          rnaseq_contrast <- paste0("(", orig_name(parts[1]), ")v(", orig_name(parts[2]), ")")
        } else {
          next
        }
      } else {
        next
      }
      prefix <- switch(suffix,
        `_diff` = "Log2fc_", `_t` = "Stat_", `_p.val` = "P.value_",
        `_p.adj` = "Adj.p.value_", `_CI.L` = "CI.L_", `_CI.R` = "CI.R_",
        `_significant` = "Significant_", suffix)
      new_names[colnames(de_df) == c] <- paste0(prefix, rnaseq_contrast)
    }
  }
  # Group stats: Mean_X -> Group.Mean_(X), Stdev_X -> Group.Stdev_(X)
  mean_cols <- grep("^Mean_", colnames(de_df), value = TRUE)
  for (c in mean_cols) {
    x <- sub("^Mean_", "", c)
    new_names[colnames(de_df) == c] <- paste0("Group.Mean_(", orig_name(x), ")")
  }
  stdev_cols <- grep("^Stdev_", colnames(de_df), value = TRUE)
  for (c in stdev_cols) {
    x <- sub("^Stdev_", "", c)
    new_names[colnames(de_df) == c] <- paste0("Group.Stdev_(", orig_name(x), ")")
  }
  colnames(de_df) <- new_names
  de_df
}

# add_rejections_customized: marks significant from p.adj and diff
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


# get_cluster_heatmap_customized: DE heatmap of significant features. FragPipeAnalystR.
get_cluster_heatmap_customized <- function(dep, type = c("contrast", "centered"),
                                           kmeans = FALSE, k = 6, col_limit = 6, indicate = NULL,
                                           alpha = 0.01, lfc = 1, plot = TRUE,
                                           clustering_distance = c("euclidean", "maximum", "manhattan",
                                               "canberra", "binary", "minkowski", "pearson", "spearman", "kendall", "gower"),
                                           row_font_size = 6, col_font_size = 10, ...) {
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
  if (!is.null(indicate) && indicate %in% colnames(col_data) && type == "centered") {
    anno <- as.data.frame(colData(dep)) %>% dplyr::select(dplyr::all_of(indicate))
    if (indicate == "condition" && "condition_raw" %in% colnames(col_data)) anno[[1]] <- col_data$condition_raw
    var <- sort(unique(anno[[1]]))
    nv <- length(var)
    cols <- if (nv == 1) c("black") else if (nv == 2) c("orangered", "cornflowerblue") else
      if (nv <= 6) RColorBrewer::brewer.pal(max(3, nv), "Pastel1")[seq_len(nv)] else
      if (nv <= 12) RColorBrewer::brewer.pal(nv, "Set3") else
      colorRampPalette(RColorBrewer::brewer.pal(12, "Set3"))(nv)
    names(cols) <- var
    ha1 <- ComplexHeatmap::HeatmapAnnotation(df = anno, col = setNames(list(cols), indicate), show_annotation_name = TRUE)
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
    return(NULL)
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
    raw_map <- sanitized_to_raw_map(dep)
    if (!is.null(raw_map)) {
      for (i in seq_along(colnames(df))) {
        cn <- colnames(df)[i]
        if (grepl(" vs ", cn)) {
          parts <- strsplit(cn, " vs ")[[1]]
          if (length(parts) == 2) {
            colnames(df)[i] <- paste0(raw_map(trimws(parts[1])), " vs ", raw_map(trimws(parts[2])))
          }
        }
      }
    }
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
  if (!plot) return(ht1)
  p <- ComplexHeatmap::draw(ht1, heatmap_legend_side = "top")
  list(ht1, ComplexHeatmap::row_order(p))
}


# volcano_font_sizes: scale title and corner-label font sizes for long condition names.
# Returns list(title_size, corner_size). Use when display/name1/name2 exceed thresholds.
volcano_font_sizes <- function(display, name1, name2,
                               title_thresh = 40, corner_thresh = 30,
                               title_base = 12, title_min = 7, corner_base = 5, corner_min = 2.5) {
  title_len <- nchar(display)
  corner_len <- max(nchar(name1), nchar(name2))
  title_size <- if (title_len > title_thresh) {
    extra <- title_len - title_thresh
    max(title_min, title_base - 0.15 * extra)
  } else title_base
  corner_size <- if (corner_len > corner_thresh) {
    extra <- corner_len - corner_thresh
    max(corner_min, corner_base - 0.1 * extra)
  } else corner_base
  list(title_size = title_size, corner_size = corner_size)
}

# volcano_plot_dims: suggest width/height for ggsave when condition names are long.
volcano_plot_dims <- function(contrast, dep, base_width = 8, base_height = 6) {
  display <- contrast_to_display(contrast, dep)
  name1 <- gsub("_vs_.*", "", contrast)
  name2 <- gsub(".*_vs_", "", contrast)
  raw_map <- sanitized_to_raw_map(dep)
  if (!is.null(raw_map)) { name1 <- raw_map(name1); name2 <- raw_map(name2) }
  max_len <- max(nchar(display), nchar(name1), nchar(name2))
  if (max_len <= 45) return(c(base_width, base_height))
  extra <- max_len - 45
  w <- min(14, base_width + 0.08 * extra)
  h <- min(9, base_height + 0.04 * extra)
  c(w, h)
}

# plot_peptide_volcano: peptide/site volcano with highlight + show_other_peptides. FragPipeAnalystR.
# peptides: IDs to highlight (maroon). show_other_peptides=T: also show other peptides from same protein (blue); ID prefix = protein.
# When peptides is NA or empty, falls back to plot_volcano_customized.
plot_peptide_volcano <- function(dep, contrast, peptides = NA, show_other_peptides = TRUE, show_gene = FALSE,
                                 label_size = 3, name_col = NULL, add_names = TRUE, adjusted = TRUE,
                                 alpha = 0.05, lfc = 1) {
  pep_vec <- if (is.null(peptides) || (length(peptides) == 1 && is.na(peptides)) || length(peptides) == 0) {
    character(0)
  } else {
    as.character(peptides)
  }
  if (length(pep_vec) == 0) {
    return(plot_volcano_customized(dep, contrast, label_size = label_size, name_col = name_col,
      add_names = add_names, adjusted = adjusted, lfc = lfc, alpha = alpha, plot = TRUE, show_gene = show_gene))
  }
  if (is.integer(label_size)) label_size <- as.numeric(label_size)
  assertthat::assert_that(inherits(dep, "SummarizedExperiment"), is.character(contrast), length(contrast) == 1,
    is.numeric(label_size), is.logical(add_names), is.logical(adjusted),
    !is.null(metadata(dep)$level), metadata(dep)$level %in% c("peptide", "site"))
  row_data <- SummarizedExperiment::rowData(dep, use.names = FALSE)
  if (is.null(name_col)) name_col <- "ID"
  gene_col <- if ("Gene" %in% colnames(row_data)) "Gene" else if ("Genes" %in% colnames(row_data)) "Genes" else NULL
  if (is.null(gene_col)) row_data$Gene <- as.character(row_data$ID) else row_data$Gene <- as.character(row_data[[gene_col]])
  if (any(!c("name", "ID", name_col) %in% colnames(row_data))) stop("'name'/'ID' not in rowData", call. = FALSE)
  if (length(grep("_p.adj|_diff", colnames(row_data))) < 1) stop("Run test_diff first", call. = FALSE)
  if (length(grep("_significant", colnames(row_data))) < 1) stop("Run add_rejections first", call. = FALSE)
  if (length(grep(paste0("^", contrast, "_diff"), colnames(row_data))) == 0) stop("Invalid contrast", call. = FALSE)
  diff_col <- grep(paste0("^", contrast, "_diff"), colnames(row_data))
  p_values_col <- grep(paste0("^", contrast, "_p.adj"), colnames(row_data))
  if (length(p_values_col) == 0) p_values_col <- grep(paste0("^", contrast, "_p.val"), colnames(row_data))
  signif <- abs(row_data[, diff_col]) >= lfc & row_data[, p_values_col] <= alpha
  df <- data.frame(diff = row_data[, diff_col], p_values = -log10(row_data[, p_values_col]),
    signif = signif, name = row_data$name, ID = row_data$ID, label = row_data[, name_col], Gene = row_data$Gene)
  df <- df %>% dplyr::filter(!is.na(signif)) %>% dplyr::arrange(signif)
  name1 <- gsub("_vs_.*", "", contrast)
  name2 <- gsub(".*_vs_", "", contrast)
  raw_map <- sanitized_to_raw_map(dep)
  if (!is.null(raw_map)) {
    name1 <- raw_map(name1)
    name2 <- raw_map(name2)
  }
  display <- contrast_to_display(contrast, dep)
  if (show_gene) df$ID_new <- paste0(df$Gene, gsub(".*_", "_", df$ID))
  label_col <- if (show_gene) "ID_new" else "ID"
  fs <- volcano_font_sizes(display, name1, name2)
  p <- ggplot2::ggplot(df, ggplot2::aes(diff, p_values)) +
    ggplot2::geom_vline(xintercept = 0) +
    ggplot2::geom_point(ggplot2::aes(col = signif)) +
    ggplot2::geom_text(data = data.frame(), ggplot2::aes(x = c(Inf, -Inf), y = c(-Inf, -Inf), hjust = c(1, 0), vjust = c(-1, -1),
      label = c(name1, name2), size = fs$corner_size, fontface = "bold")) +
    ggplot2::labs(title = display, x = expression(log[2] ~ "Fold change")) +
    ggplot2::theme_bw() +
    ggplot2::theme(plot.title = ggplot2::element_text(size = fs$title_size),
      panel.border = ggplot2::element_blank(), panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(), axis.line = ggplot2::element_line(colour = "black"), legend.position = "none") +
    ggplot2::scale_color_manual(values = c("TRUE" = "black", "FALSE" = "grey"))
  if (show_other_peptides && length(pep_vec) > 0) {
    gene_prefix <- unique(gsub("_.*", "", pep_vec))
    other_same_gene <- df$ID[gsub("_.*", "", df$ID) %in% gene_prefix & !df$ID %in% pep_vec]
    if (length(other_same_gene) > 0) {
      p <- p + ggplot2::geom_point(data = dplyr::filter(df, .data$ID %in% other_same_gene), color = "blue", size = 3)
    }
  }
  pep_in_df <- intersect(pep_vec, df$ID)
  if (length(pep_in_df) > 0) {
    p <- p + ggplot2::geom_point(data = dplyr::filter(df, .data$ID %in% pep_in_df), color = "maroon", size = 3) +
      ggrepel::geom_text_repel(data = dplyr::filter(df, .data$ID %in% pep_in_df), color = "maroon",
        ggplot2::aes(label = .data[[label_col]]), size = label_size, box.padding = grid::unit(0.1, "lines"),
        point.padding = grid::unit(0.1, "lines"), segment.size = 0.5, max.overlaps = 100)
  }
  if (add_names) {
    repel_df <- if (length(pep_in_df) > 0) dplyr::filter(df, signif, !.data$ID %in% pep_in_df) else dplyr::filter(df, signif)
    if (nrow(repel_df) > 0) {
      p <- p + ggrepel::geom_text_repel(data = repel_df, ggplot2::aes(label = .data[[label_col]]),
        size = label_size, box.padding = grid::unit(0.1, "lines"), point.padding = grid::unit(0.1, "lines"), segment.size = 0.5, max.overlaps = 100)
    }
  }
  p <- p + ggplot2::labs(y = if (adjusted) expression(-log[10] ~ "Adjusted p-value") else expression(-log[10] ~ "P-value"))
  p
}

# plot_volcano_customized: volcano plot per contrast. FragPipe-Analyst.
plot_volcano_customized <- function(dep, contrast, label_size = 3, name_col = NULL,
                                    add_names = TRUE, adjusted = TRUE, lfc = 1, alpha = 0.05,
                                    plot = TRUE, show_gene = FALSE, selected = NULL) {
  if (is.integer(label_size)) label_size <- as.numeric(label_size)
  assertthat::assert_that(
    inherits(dep, "SummarizedExperiment"),
    is.character(contrast), length(contrast) == 1,
    is.numeric(label_size), length(label_size) == 1,
    is.logical(add_names), length(add_names) == 1,
    is.logical(adjusted), length(adjusted) == 1,
    is.logical(plot), length(plot) == 1
  )
  row_data <- SummarizedExperiment::rowData(dep, use.names = FALSE)
  if (is.null(name_col)) name_col <- "ID"
  if (any(!c("name", "ID", name_col) %in% colnames(row_data))) {
    stop("'name' and/or 'ID' columns are not present. Run make_unique() first.", call. = FALSE)
  }
  if (length(grep("_p.adj|_diff", colnames(row_data))) < 1) {
    stop("'[contrast]_diff' and '[contrast]_p.adj' columns are not present. Run test_diff() and add_rejections().", call. = FALSE)
  }
  if (length(grep("_significant", colnames(row_data))) < 1) {
    stop("'[contrast]_significant' column not present. Run add_rejections() first.", call. = FALSE)
  }
  if (length(grep(paste("^", contrast, "_diff", sep = ""), colnames(row_data))) == 0) {
    valid <- row_data %>% data.frame() %>%
      dplyr::select(dplyr::ends_with("_diff")) %>% colnames() %>% gsub("_diff", "", .)
    stop("Not a valid contrast. Valid: ", paste0("'", valid, "'", collapse = ", "), call. = FALSE)
  }
  diff_col <- grep(paste("^", contrast, "_diff", sep = ""), colnames(row_data))
  p_values_col <- if (adjusted) {
    grep(paste("^", contrast, "_p.adj", sep = ""), colnames(row_data))
  } else {
    grep(paste("^", contrast, "_p.val", sep = ""), colnames(row_data))
  }
  signif <- abs(row_data[, diff_col]) >= lfc & row_data[, p_values_col] <= alpha
  exp <- if (!is.null(metadata(dep)$exp)) metadata(dep)$exp else "LFQ"
  lvl <- if (!is.null(metadata(dep)$level)) metadata(dep)$level else "protein"
  if (!show_gene) {
    if (exp == "LFQ") {
      if (lvl != "peptide") {
        df_tmp <- data.frame(diff = row_data[, diff_col], p_values = -log10(row_data[, p_values_col]),
                             signif = signif, name = row_data$name, ID = row_data$ID, label = row_data[, name_col])
      } else {
        df_tmp <- data.frame(diff = row_data[, diff_col], p_values = -log10(row_data[, p_values_col]),
                             signif = signif, name = if ("Index" %in% colnames(row_data)) row_data$Index else row_data$name,
                             ID = row_data$ID, label = row_data[, name_col])
      }
    } else if (exp == "TMT") {
      if (lvl == "protein") {
        df_tmp <- data.frame(diff = row_data[, diff_col], p_values = -log10(row_data[, p_values_col]),
                             signif = signif, name = row_data$ID, ID = row_data$ID, label = row_data[, name_col])
      } else if (lvl == "gene") {
        df_tmp <- data.frame(diff = row_data[, diff_col], p_values = -log10(row_data[, p_values_col]),
                             signif = signif, name = row_data$ID, ID = row_data$ID, label = row_data[, name_col])
      } else {
        df_tmp <- data.frame(diff = row_data[, diff_col], p_values = -log10(row_data[, p_values_col]),
                             signif = signif, name = if ("Index" %in% colnames(row_data)) row_data$Index else row_data$name,
                             ID = row_data$ID, label = row_data[, name_col])
      }
    } else if (exp == "DIA") {
      if (!lvl %in% c("site", "peptide")) {
        df_tmp <- data.frame(diff = row_data[, diff_col], p_values = -log10(row_data[, p_values_col]),
                             signif = signif, name = row_data$ID, ID = row_data$ID, label = row_data[, name_col])
      } else {
        df_tmp <- data.frame(diff = row_data[, diff_col], p_values = -log10(row_data[, p_values_col]),
                             signif = signif, name = if ("Index" %in% colnames(row_data)) row_data$Index else row_data$name,
                             ID = row_data$ID, label = row_data[, name_col])
      }
    } else {
      df_tmp <- data.frame(diff = row_data[, diff_col], p_values = -log10(row_data[, p_values_col]),
                           signif = signif, name = row_data$name, ID = row_data$ID, label = row_data[, name_col])
    }
  } else {
    if (exp == "LFQ") {
      if (lvl != "peptide") {
        df_tmp <- data.frame(diff = row_data[, diff_col], p_values = -log10(row_data[, p_values_col]),
                             signif = signif, name = if ("Gene" %in% colnames(row_data)) row_data$Gene else row_data$name,
                             ID = row_data$ID, label = if ("Gene" %in% colnames(row_data)) row_data$Gene else row_data[, name_col])
      } else {
        df_tmp <- data.frame(diff = row_data[, diff_col], p_values = -log10(row_data[, p_values_col]),
                             signif = signif, name = paste0(if ("Gene" %in% colnames(row_data)) row_data$Gene else row_data$ID, "_", if ("Peptide.Sequence" %in% colnames(row_data)) row_data$Peptide.Sequence else row_data$ID),
                             ID = row_data$ID, label = row_data[, name_col])
      }
    } else if (exp == "TMT") {
      if (lvl == "protein") {
        df_tmp <- data.frame(diff = row_data[, diff_col], p_values = -log10(row_data[, p_values_col]),
                             signif = signif, name = if ("Gene" %in% colnames(row_data)) row_data$Gene else row_data$ID,
                             ID = row_data$ID, label = row_data[, name_col])
      } else if (lvl == "gene") {
        df_tmp <- data.frame(diff = row_data[, diff_col], p_values = -log10(row_data[, p_values_col]),
                             signif = signif, name = row_data$ID, ID = row_data$ID, label = row_data[, name_col])
      } else {
        nm <- if (lvl == "site" && "ID" %in% colnames(row_data)) {
          paste0(if ("Gene" %in% colnames(row_data)) row_data$Gene else row_data$ID, "_", gsub(".*_", "", row_data$ID))
        } else {
          paste0(if ("Gene" %in% colnames(row_data)) row_data$Gene else row_data$ID, "_", if ("Peptide" %in% colnames(row_data)) row_data$Peptide else row_data$ID)
        }
        df_tmp <- data.frame(diff = row_data[, diff_col], p_values = -log10(row_data[, p_values_col]),
                             signif = signif, name = nm, ID = row_data$ID, label = row_data[, name_col])
      }
    } else if (exp == "DIA") {
      if (lvl == "peptide") {
        nm <- if ("Gene" %in% colnames(row_data)) paste0(row_data$Gene, "_", if ("Peptide" %in% colnames(row_data)) row_data$Peptide else row_data$ID) else paste0(row_data$ID, "_", row_data$ID)
        df_tmp <- data.frame(diff = row_data[, diff_col], p_values = -log10(row_data[, p_values_col]),
                             signif = signif, name = nm, ID = row_data$ID, label = row_data[, name_col])
      } else if (lvl == "site") {
        df_tmp <- data.frame(diff = row_data[, diff_col], p_values = -log10(row_data[, p_values_col]),
                             signif = signif, name = paste0(if ("Gene" %in% colnames(row_data)) row_data$Gene else row_data$ID, "_", gsub(".*_", "", row_data$ID)),
                             ID = row_data$ID, label = row_data[, name_col])
      } else {
        gene_col <- if ("Genes" %in% colnames(row_data)) "Genes" else "Gene"
        df_tmp <- data.frame(diff = row_data[, diff_col], p_values = -log10(row_data[, p_values_col]),
                             signif = signif, name = row_data[, gene_col], ID = row_data$ID, label = row_data[, name_col])
      }
    } else {
      df_tmp <- data.frame(diff = row_data[, diff_col], p_values = -log10(row_data[, p_values_col]),
                           signif = signif, name = row_data$name, ID = row_data$ID, label = row_data[, name_col])
    }
  }
  df <- df_tmp %>% data.frame() %>% dplyr::filter(!is.na(signif)) %>% dplyr::arrange(signif)
  display <- contrast_to_display(contrast, dep)
  name1 <- gsub("_vs_.*", "", contrast)
  name2 <- gsub(".*_vs_", "", contrast)
  raw_map <- sanitized_to_raw_map(dep)
  if (!is.null(raw_map)) {
    name1 <- raw_map(name1)
    name2 <- raw_map(name2)
  }
  fs <- volcano_font_sizes(display, name1, name2)
  p <- ggplot2::ggplot(df, ggplot2::aes(diff, p_values)) +
    ggplot2::geom_vline(xintercept = 0) +
    ggplot2::geom_point(ggplot2::aes(col = signif)) +
    ggplot2::geom_text(data = data.frame(), ggplot2::aes(x = c(Inf, -Inf), y = c(-Inf, -Inf),
                                                       hjust = c(1, 0), vjust = c(-1, -1),
                                                       label = c(name1, name2), size = fs$corner_size, fontface = "bold")) +
    ggplot2::labs(title = display, x = expression(log[2] ~ "Fold change")) +
    ggplot2::theme_bw() +
    ggplot2::theme(plot.title = ggplot2::element_text(size = fs$title_size), legend.position = "none") +
    ggplot2::scale_color_manual(values = c("TRUE" = "black", "FALSE" = "grey"))
  if (add_names) {
    if (!is.null(selected)) {
      p <- p + ggrepel::geom_text_repel(data = dplyr::filter(df, signif, !name %in% selected),
                                        ggplot2::aes(label = name), size = label_size,
                                        box.padding = grid::unit(0.1, "lines"), point.padding = grid::unit(0.1, "lines"), segment.size = 0.5, max.overlaps = 100)
    } else {
      p <- p + ggrepel::geom_text_repel(data = dplyr::filter(df, signif),
                                        ggplot2::aes(label = name), size = label_size,
                                        box.padding = grid::unit(0.1, "lines"), point.padding = grid::unit(0.1, "lines"), segment.size = 0.5, max.overlaps = 100)
    }
  }
  if (adjusted) {
    p <- p + ggplot2::labs(y = expression(-log[10] ~ "Adjusted p-value"))
  } else {
    p <- p + ggplot2::labs(y = expression(-log[10] ~ "P-value"))
  }
  if (plot) return(p)
  df_out <- df %>% dplyr::select(name, diff, p_values, signif)
  colnames(df_out)[c(1, 2, 3)] <- c("protein", "log2_fold_change", "p_value_-log10")
  if (adjusted) colnames(df_out)[3] <- "adjusted_p_value_-log10"
  return(df_out)
}
