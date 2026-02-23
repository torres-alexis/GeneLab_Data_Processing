# fp_se_helper.R
# SummarizedExperiment creation for FragPipe quant + annotation.
# From FragPipeAnalystR, DEP.
# Fixes: readExpDesign always lowercases colnames; sample_name <- label for plot joins.
# Filter functions from FragPipe-Analyst.

library(SummarizedExperiment)
library(dplyr)
library(tibble)

# Filter: keep rows with <= percentage (0-100) fraction of NAs globally.
global_filter <- function(se, percentage = 50) {
  percentage <- percentage / 100
  ridx <- rowSums(is.na(assay(se))) / ncol(assay(se)) <= percentage
  se <- se[ridx, ]
  return(se)
}

# Filter: keep rows with >= min_percentage (0-100) valid in at least one condition.
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

# Normalization. From FragPipeAnalystR.
# MD: median subtraction per sample. GN: median + MAD scaling. VSN: variance-stabilizing (LFQ/DIA intensity only).
MD_normalization <- function(se) {
  data <- assay(se)
  assay(se) <- sweep(data, 2, matrixStats::colMedians(data, na.rm = TRUE))
  se
}

GN_normalization <- function(se) {
  data <- assay(se)
  MD <- sweep(data, 2, matrixStats::colMedians(data, na.rm = TRUE))
  MAD <- apply(MD, 2, mad, na.rm = TRUE)
  MAD_0 <- median(MAD, na.rm = TRUE)
  MAD[MAD == 0 | is.na(MAD)] <- MAD_0
  assay(se) <- sweep(MD, 2, MAD, FUN = "/") * MAD_0
  se
}

VSN_normalization <- function(se) {
  if (!metadata(se)$exp %in% c("LFQ", "DIA")) {
    stop("VSN normalization is only for LFQ or DIA (not TMT)")
  }
  if (!is.null(metadata(se)$lfq_type) && metadata(se)$lfq_type == "Spectral Count") {
    stop("VSN normalization is not for Spectral Count data")
  }
  data <- assay(se)
  vsn.fit <- vsn::vsnMatrix(2^data)
  assay(se) <- vsn::predict(vsn.fit, 2^data)
  se
}

# Wrapper: apply normalization method. Returns se unchanged if method is "none".
# VSN only for LFQ/DIA intensity; MD and GN for all.
normalize_se <- function(se, method = "none") {
  if (method == "none") return(se)
  if (method == "MD") return(MD_normalization(se))
  if (method == "GN") return(GN_normalization(se))
  if (method == "vsn") return(VSN_normalization(se))
  stop("Invalid normalization method: ", method)
}

# Imputation. From FragPipeAnalystR, DEP, FragPipe-Analyst.
# manual_impute: Perseus-type, per-sample rnorm(median - shift*sd, sd*scale). No extra deps.
manual_impute <- function(se, scale = 0.3, shift = 1.8, seed = 123, ...) {
  if (is.integer(scale)) scale <- as.numeric(scale)
  if (is.integer(shift)) shift <- as.numeric(shift)
  se_assay <- assay(se)
  if (!any(is.na(se_assay))) {
    stop("No missing values in '", deparse(substitute(se)), "'", call. = FALSE)
  }
  set.seed(seed)
  for (i in seq_len(ncol(se_assay))) {
    vals <- se_assay[, i]
    n_infin <- sum(is.na(vals))
    if (n_infin == 0) next
    med <- median(vals, na.rm = TRUE)
    s <- sd(vals, na.rm = TRUE)
    if (is.na(s) || s == 0) s <- 1e-10
    assay(se)[is.na(assay(se)[, i]), i] <- rnorm(n_infin, mean = med - shift * s, sd = s * scale)
  }
  se
}

# impute_se: wrapper. Perseus-type = manual_impute; else MSnbase::impute (requires MSnbase).
# Valid fun: Perseus-type, knn, MLE, min, zero, bpca, QRILC, MinDet, MinProb, RF, nbavg, mixed.
impute_se <- function(se, fun = c("Perseus-type", "knn", "MLE", "min", "zero", "bpca", "QRILC",
    "MinDet", "MinProb", "RF", "nbavg", "mixed"), seed = 123, ...) {
  fun <- match.arg(fun)
  if (any(!c("name", "ID") %in% colnames(rowData(se, use.names = FALSE)))) {
    stop("'name' and 'ID' required in rowData. Run make_unique() and make_se first.", call. = FALSE)
  }
  rowData(se)$imputed <- apply(is.na(assay(se)), 1, any)
  rowData(se)$num_NAs <- rowSums(is.na(assay(se)))
  se <- se[!rowData(se)$num_NAs == ncol(se), ]
  if (!any(is.na(assay(se)))) {
    warning("No missing values. Returning unchanged object.", call. = FALSE)
    return(se)
  }
  if (fun == "Perseus-type") {
    return(manual_impute(se, seed = seed, ...))
  }
  if (!requireNamespace("MSnbase", quietly = TRUE)) {
    stop("MSnbase required for imputation method '", fun, "'. Install with BiocManager::install('MSnbase')", call. = FALSE)
  }
  MSnSet_data <- as(se, "MSnSet")
  set.seed(seed)
  MSnSet_imputed <- MSnbase::impute(MSnSet_data, method = fun, ...)
  assay(se) <- MSnbase::exprs(MSnSet_imputed)
  se
}

make.unique.2 <- function(x, sep = ".") {
  ave(x, x, FUN = function(a) {
    if (length(a) > 1) paste(a, 1:length(a), sep = sep) else a
  })
}

readQuantTable <- function(quant_table_path, type = "TMT", level = NULL, log2transform = FALSE,
  exp_type = NULL, additional_cols = NULL) {
  temp_data <- read.table(quant_table_path, header = TRUE, fill = TRUE, sep = "\t",
    quote = "", comment.char = "", blank.lines.skip = FALSE, check.names = FALSE)
  colnames(temp_data) <- make.unique.2(colnames(temp_data), "_")
  if (type == "TMT") {
    skip <- c("Index", "Gene", "Peptide", "NumberPSM", "ProteinID", "MaxPepProb", "SequenceWindow", "ReferenceIntensity")
    if (!is.null(additional_cols)) skip <- c(skip, additional_cols)
    mut.cols <- colnames(temp_data)[!colnames(temp_data) %in% skip]
    temp_data[mut.cols] <- sapply(temp_data[mut.cols], as.numeric)
  } else if (type == "LFQ") {
    if (level == "peptide") {
      colnames(temp_data) <- gsub("-", ".", colnames(temp_data))
      colnames(temp_data)[colnames(temp_data) == "Protein Description"] <- "Description"
      temp_data <- temp_data[!grepl("contam", temp_data$Protein), ]
      if (!"Modified Sequence" %in% colnames(temp_data)) {
        temp_data$Index <- paste0(temp_data$`Protein ID`, "_", temp_data$`Peptide Sequence`)
      } else {
        temp_data$Index <- paste0(temp_data$`Protein ID`, "_", temp_data$`Modified Sequence`)
      }
    } else {
      colnames(temp_data) <- gsub("-", ".", colnames(temp_data))
      temp_data <- temp_data[!grepl("contam", temp_data$Protein), ]
    }
  } else {
    library(data.table)
    if (level == "peptide") {
      if ("SequenceWindow" %in% colnames(temp_data)) {
        stop("Wrong format for level=peptide. Use level=site for single-site report.")
      }
      drop_cols <- c("Proteotypic", "Precursor.Charge", "Precursor.Id", "Modified.Sequence",
        "First.Protein.Description", "All Mapped Proteins", "All Mapped Genes")
      drop_cols <- drop_cols[drop_cols %in% colnames(temp_data)]
      temp <- data.table::melt.data.table(data.table::setDT(
        temp_data[, !colnames(temp_data) %in% drop_cols, drop = FALSE]),
        id.vars = c("Protein.Group", "Protein.Names", "Protein.Ids", "Genes", "Stripped.Sequence"),
        variable.name = "File.Name")
      temp_data <- as.data.frame(data.table::dcast.data.table(temp,
        Protein.Group + Protein.Names + Protein.Ids + Genes + Stripped.Sequence ~ File.Name,
        value.var = "value", fun.aggregate = function(x) max(x, na.rm = TRUE)))
      temp_data[sapply(temp_data, is.infinite)] <- NA
      temp_data$Index <- paste0(temp_data$Protein.Ids, "_", temp_data$Stripped.Sequence)
      temp_data <- temp_data[, c("Index", setdiff(colnames(temp_data), "Index"))]
    } else if (level == "site") {
      if (!"SequenceWindow" %in% colnames(temp_data)) {
        stop("No SequenceWindow column. Use single-site report for level=site.")
      }
      skip <- c("Index", "ProteinID", "Gene", "Peptide", "SequenceWindow")
      if (!is.null(additional_cols)) skip <- c(skip, additional_cols)
      mut.cols <- colnames(temp_data)[!colnames(temp_data) %in% skip]
      temp_data[mut.cols] <- sapply(temp_data[mut.cols], as.numeric)
    }
  }
  temp_data
}

readExpDesign <- function(exp_anno_path, type = "TMT", lfq_type = "Intensity") {
  temp_df <- read.table(exp_anno_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  colnames(temp_df) <- tolower(colnames(temp_df))
  if (type == "TMT") {
    if (ncol(temp_df) == 1) {
      temp_df <- tryCatch(read.table(exp_anno_path, header = TRUE, sep = " ", stringsAsFactors = FALSE),
        error = function(e) temp_df)
    }
    temp_df$condition <- make.names(temp_df$condition)
    temp_df$label <- temp_df$sample
    if (anyDuplicated(temp_df$label)) {
      temp_df$label <- paste(temp_df$label, temp_df$replicate, sep = "_")
      temp_df$label <- gsub("_1$", "", temp_df$label)
      samples_with_replicate <- unique(gsub("_\\d+$", "", temp_df$label[grepl("_", temp_df$label)]))
      idx <- temp_df$label %in% samples_with_replicate
      temp_df$label[idx] <- paste0(temp_df$label[idx], "_1")
    }
  } else if (type == "LFQ") {
    temp_df$condition <- make.names(temp_df$condition)
    if (!all(is.na(temp_df$replicate))) {
      temp_df$sample <- gsub("-", ".", temp_df$sample)
      temp_df$label <- temp_df$sample
      if (lfq_type == "Intensity") temp_df$label <- paste(temp_df$label, "Intensity", sep = " ")
      else if (lfq_type == "MaxLFQ") temp_df$label <- paste(temp_df$label, "MaxLFQ.Intensity", sep = " ")
      else if (lfq_type == "Spectral Count") temp_df$label <- paste(temp_df$label, "Spectral.Count", sep = " ")
    }
    # Use short sample name for display (plots); prefer sample_name from file (basename) when set
    if ("sample_name" %in% colnames(temp_df) && all(nzchar(trimws(temp_df$sample_name)))) {
      # Keep sample_name from file (e.g. runsheet_to_experiment_annotation basename)
    } else if ("sample" %in% colnames(temp_df) && all(nzchar(trimws(temp_df$sample)))) {
      temp_df$sample_name <- temp_df$sample
    } else {
      temp_df$sample_name <- temp_df$label
    }
  } else {
    temp_df$condition <- make.names(temp_df$condition)
    if (!all(is.na(temp_df$replicate))) temp_df$label <- temp_df$file
  }
  if (!"sample_name" %in% colnames(temp_df) || !all(nzchar(trimws(temp_df$sample_name)))) {
    temp_df$sample_name <- temp_df$label
  }
  temp_df
}

make_unique <- function(proteins, names, ids, delim = ";") {
  col_names <- colnames(proteins)
  if (!names %in% col_names) stop("'", names, "' not in ", deparse(substitute(proteins)))
  if (!ids %in% col_names) stop("'", ids, "' not in ", deparse(substitute(proteins)))
  if (tibble::is_tibble(proteins)) proteins <- as.data.frame(proteins)
  double_NAs <- apply(proteins[, c(names, ids)], 1, function(x) all(is.na(x)))
  if (any(double_NAs)) stop("NAs in both 'names' and 'ids' columns")
  proteins %>%
    dplyr::mutate(
      name = .data[[names]],
      ID = .data[[ids]],
      name = make.unique(ifelse(name == "" | is.na(name), ID, name))
    )
}

make_se_customized <- function(proteins_unique, columns, expdesign, log2transform = FALSE,
  exp = "LFQ", lfq_type = NULL, level = NULL, exp_type = NULL) {
  if (any(!c("name", "ID") %in% colnames(proteins_unique)))
    stop("Run make_unique() first")
  if (any(!c("label", "condition", "replicate", "sample_name") %in% colnames(expdesign)))
    stop("expdesign needs label, condition, replicate, sample_name")
  if (any(!apply(proteins_unique[, columns], 2, is.numeric)))
    stop("specified columns must be numeric")
  if (tibble::is_tibble(proteins_unique)) proteins_unique <- as.data.frame(proteins_unique)
  if (tibble::is_tibble(expdesign)) expdesign <- as.data.frame(expdesign)
  rownames(proteins_unique) <- proteins_unique$ID
  raw <- proteins_unique[, columns]
  raw[raw == 0] <- NA
  if (log2transform) raw <- log2(raw)
  rownames(expdesign) <- expdesign$label
  matched <- match(make.names(expdesign$label), make.names(colnames(raw)))
  if (any(is.na(matched))) {
    stop("Labels in experimental design do not match column names in quantification table")
  }
  rownames(expdesign) <- expdesign$sample_name
  colnames(raw)[matched] <- expdesign$sample_name
  raw <- raw[, rownames(expdesign), drop = FALSE]
  row_data <- proteins_unique[, -columns, drop = FALSE]
  rownames(row_data) <- row_data$ID
  SummarizedExperiment(
    assays = as.matrix(raw),
    colData = expdesign,
    rowData = row_data,
    metadata = list(log2transform = log2transform, exp = exp, lfq_type = lfq_type,
      exp_type = exp_type, level = level)
  )
}

make_se_from_files <- function(quant_table_path, exp_anno_path, type = "TMT", level = NULL,
  exp_type = NULL, log2transform = NULL, lfq_type = "Intensity", gencode = FALSE, additional_cols = NULL) {
  if (type == "TMT" && is.null(level)) level <- "gene"
  else if (is.null(level)) level <- "protein"
  if (type == "DIA" && is.null(log2transform)) log2transform <- TRUE
  else if (is.null(log2transform)) log2transform <- FALSE
  if (!level %in% c("gene", "protein", "peptide", "site", "glycan")) {
    stop("level must be gene, protein, peptide, site, or glycan")
  }
  quant_table <- readQuantTable(quant_table_path, type = type, level = level,
    exp_type = exp_type, additional_cols = additional_cols)
  if (is.null(quant_table)) return(NULL)
  exp_design <- readExpDesign(exp_anno_path, type = type, lfq_type = lfq_type)
  if (type == "LFQ") {
    quant_table <- quant_table[!grepl("contam", quant_table$Protein), ]
    if (level == "site") {
      data_unique <- make_unique(quant_table, "Protein ID", "Index")
      if (lfq_type == "Intensity") {
        lfq_columns <- setdiff(grep("Intensity", colnames(data_unique)), grep("MaxLFQ", colnames(data_unique)))
        lfq_columns <- setdiff(lfq_columns, grep("Total Intensity", colnames(data_unique)))
        lfq_columns <- setdiff(lfq_columns, grep("Unique Intensity", colnames(data_unique)))
      } else if (lfq_type == "MaxLFQ") {
        lfq_columns <- grep("MaxLFQ", colnames(data_unique))
        if (length(lfq_columns) == 0) stop("No MaxLFQ columns found")
      } else {
        lfq_columns <- setdiff(grep("Spectral", colnames(data_unique)), grep("Total Spectral Count", colnames(data_unique)))
        lfq_columns <- setdiff(lfq_columns, grep("Unique Spectral Count", colnames(data_unique)))
      }
      data_se <- make_se_customized(data_unique, lfq_columns, exp_design,
        log2transform = (lfq_type != "Spectral Count"), exp = "LFQ", lfq_type = lfq_type, level = "site")
    } else if (level != "peptide") {
      data_unique <- make_unique(quant_table, "Gene", "Protein ID")
      if (lfq_type == "Intensity") {
        lfq_columns <- setdiff(grep("Intensity", colnames(data_unique)), grep("MaxLFQ", colnames(data_unique)))
        lfq_columns <- setdiff(lfq_columns, grep("Total Intensity", colnames(data_unique)))
        lfq_columns <- setdiff(lfq_columns, grep("Unique Intensity", colnames(data_unique)))
      } else if (lfq_type == "MaxLFQ") {
        lfq_columns <- grep("MaxLFQ", colnames(data_unique))
        if (length(lfq_columns) == 0) stop("No MaxLFQ columns found")
      } else {
        lfq_columns <- setdiff(grep("Spectral", colnames(data_unique)), grep("Total Spectral Count", colnames(data_unique)))
        lfq_columns <- setdiff(lfq_columns, grep("Unique Spectral Count", colnames(data_unique)))
      }
      data_se <- make_se_customized(data_unique, lfq_columns, exp_design,
        log2transform = (lfq_type != "Spectral Count"), exp = "LFQ", lfq_type = lfq_type, level = level)
    } else {
      data_unique <- make_unique(quant_table, "Protein ID", "Index")
      if (lfq_type == "Intensity") {
        lfq_columns <- setdiff(grep("Intensity", colnames(data_unique)), grep("MaxLFQ", colnames(data_unique)))
        lfq_columns <- setdiff(lfq_columns, grep("Total Intensity", colnames(data_unique)))
        lfq_columns <- setdiff(lfq_columns, grep("Unique Intensity", colnames(data_unique)))
      } else if (lfq_type == "MaxLFQ") {
        lfq_columns <- grep("MaxLFQ", colnames(data_unique))
        if (length(lfq_columns) == 0) stop("No MaxLFQ columns found")
      } else {
        lfq_columns <- setdiff(grep("Spectral", colnames(data_unique)), grep("Total Spectral Count", colnames(data_unique)))
        lfq_columns <- setdiff(lfq_columns, grep("Unique Spectral Count", colnames(data_unique)))
      }
      data_se <- make_se_customized(data_unique, lfq_columns, exp_design,
        log2transform = (lfq_type != "Spectral Count"), exp = "LFQ", lfq_type = lfq_type, level = "peptide")
    }
  } else if (type == "DIA") {
    if (level == "protein") {
      if (gencode) quant_table <- quant_table[grepl("^ENS", quant_table$Protein.Group), ]
      data_unique <- make_unique(quant_table, "Genes", "Protein.Group")
      cols <- colnames(data_unique)
      selected_cols <- which(!(cols %in% c("Protein.Group", "Protein.Ids", "Protein.Names", "Genes",
        "First.Protein.Description", "ID", "name", additional_cols)))
      data_se <- make_se_customized(data_unique, selected_cols, exp_design,
        log2transform = log2transform, exp = "DIA", level = "protein")
      dimnames(data_se) <- list(dimnames(data_se)[[1]], colData(data_se)$sample_name)
      colData(data_se)$label <- colData(data_se)$sample_name
    } else if (level == "gene") {
      if (gencode) quant_table <- quant_table[grepl("^ENS", quant_table$Genes), ]
      quant_table$Index <- quant_table$Genes
      data_unique <- make_unique(quant_table, "Genes", "Index")
      cols <- colnames(data_unique)
      selected_cols <- which(!(cols %in% c("Genes", "Index", "ID", "name", additional_cols)))
      data_se <- make_se_customized(data_unique, selected_cols, exp_design,
        log2transform = log2transform, exp = "DIA", level = "gene")
      dimnames(data_se) <- list(dimnames(data_se)[[1]], colData(data_se)$sample_name)
      colData(data_se)$label <- colData(data_se)$sample_name
    } else {
      if (level == "site") data_unique <- make_unique(quant_table, "ProteinID", "Index")
      else data_unique <- make_unique(quant_table, "Protein.Group", "Index")
      cols <- colnames(data_unique)
      selected_cols <- which(!(cols %in% c("Index", "Protein.Group", "Protein.Ids", "Stripped.Sequence",
        "Protein.Names", "Genes", "First.Protein.Description", "ID", "name", "Gene", "ProteinID",
        "Peptide", "SequenceWindow", "All Mapped Proteins", "All Mapped Genes", additional_cols)))
      data_se <- make_se_customized(data_unique, selected_cols, exp_design,
        log2transform = log2transform, exp = "DIA", level = level)
      dimnames(data_se) <- list(dimnames(data_se)[[1]], colData(data_se)$sample_name)
      colData(data_se)$label <- colData(data_se)$sample_name
    }
    if (level == "peptide") rowData(data_se)$Gene <- rowData(data_se)$Genes
  } else {
    temp_exp_design <- exp_design[!is.na(exp_design$condition) & exp_design$condition != "", ]
    if (level == "protein") quant_table$ProteinID <- quant_table$Index
    data_unique <- make_unique(quant_table, "ProteinID", "Index")
    overlapped_samples <- intersect(colnames(data_unique), temp_exp_design$label)
    interest_cols <- switch(level,
      gene = c("Index", "NumberPSM", "ProteinID", "MaxPepProb", "ReferenceIntensity", "name", "ID"),
      protein = c("Index", "NumberPSM", "Gene", "ProteinID", "MaxPepProb", "ReferenceIntensity", "name", "ID"),
      peptide = , site = c("Index", "Gene", "Peptide", "NumberPSM", "ProteinID", "SequenceWindow",
        "MaxPepProb", "ReferenceIntensity", "name", "ID"),
      c("Index", "Gene", "ProteinID", "Peptide", "SequenceWindow", "Start", "End", "MaxPepProb",
        "ReferenceIntensity", "name", "ID", "Spectrum Number")
    )
    if (!is.null(additional_cols)) interest_cols <- c(interest_cols, additional_cols)
    data_unique <- data_unique[, colnames(data_unique) %in% c(interest_cols, overlapped_samples)]
    temp_exp_design <- temp_exp_design[temp_exp_design$label %in% overlapped_samples, ]
    cols <- colnames(data_unique)
    selected_cols <- which(!(cols %in% interest_cols))
    data_unique[selected_cols] <- apply(data_unique[selected_cols], 2, as.numeric)
    data_se <- make_se_customized(data_unique, selected_cols, temp_exp_design, exp = "TMT", level = level)
  }
  metadata(data_se)$exp_type <- exp_type
  metadata(data_se)$level <- level
  data_se
}
