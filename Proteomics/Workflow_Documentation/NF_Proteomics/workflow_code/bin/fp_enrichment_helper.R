# fp_enrichment_helper.R
# Enrichment: or_test, plot_or. From FragPipeAnalystR.
# Uses Enrichr API (human DBs). Deps: httr, dplyr, assertthat, ggplot2.

library(httr)
library(dplyr)
library(assertthat)
library(ggplot2)

# enrichr_mod: submit genes to Enrichr, query databases. Internal.
enrichr_mod <- function(genes, databases = NULL) {
  if (length(genes) != 0) {
    if (all(startsWith(genes, "ENSG"))) {
      if (requireNamespace("ensembldb", quietly = TRUE) && requireNamespace("EnsDb.Hsapiens.v86", quietly = TRUE)) {
        genes <- gsub("\\..*", "", genes)
        genes_map <- ensembldb::select(EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86,
          keys = genes, keytype = "GENEID", columns = c("SYMBOL", "GENEID"))
        genes <- genes_map$SYMBOL
      }
    }
    httr::set_config(httr::config(ssl_verifypeer = 0L))
    cat("Uploading data to Enrichr... ")
    if (is.vector(genes) && !all(genes == "") && length(genes) != 0) {
      temp <- httr::POST(url = "http://maayanlab.cloud/Enrichr/enrich",
        body = list(list = paste(genes, collapse = "\n")))
    } else if (is.data.frame(genes)) {
      temp <- httr::POST(url = "http://maayanlab.cloud/Enrichr/enrich",
        body = list(list = paste(paste(genes[, 1], genes[, 2], sep = ","), collapse = "\n")))
    } else {
      warning("genes must be a non-empty vector of gene names or a dataframe with genes and score.")
    }
    httr::GET(url = "http://maayanlab.cloud/Enrichr/share")
    cat("Done.\n")
    dbs <- as.list(databases)
    dfSAF <- options()$stringsAsFactors
    options(stringsAsFactors = FALSE)
    result <- lapply(dbs, function(x) {
      cat("  Querying ", x, "... ", sep = "")
      r <- httr::GET(url = "http://maayanlab.cloud/Enrichr/export",
        query = list(file = "API", backgroundType = x))
      r <- gsub("&#39;", "'", intToUtf8(r$content))
      tc <- textConnection(r)
      r <- read.table(tc, sep = "\t", header = TRUE, quote = "", comment.char = "")
      close(tc)
      cat("Done.\n")
      return(r)
    })
    options(stringsAsFactors = dfSAF)
    cat("Parsing results... ")
    names(result) <- dbs
    cat("Done.\n")
  } else {
    result <- data.frame(Term = character(), Overlap = character(), P.value = double(),
      Adjusted.P.value = double(), Old.P.value = double(), Old.Adjusted.P.value = double(),
      Odds.Ratio = double(), Combined.Score = double(), Genes = character())
  }
  return(result)
}

# test_ora_mod: over-representation via Enrichr. Internal.
test_ora_mod <- function(dep, databases, contrasts = TRUE, direction = "UP",
  log2_threshold = 0.7, alpha = 0.05) {
  assertthat::assert_that(inherits(dep, "SummarizedExperiment"),
    is.character(databases), is.logical(contrasts), length(contrasts) == 1)
  row_data <- SummarizedExperiment::rowData(dep, use.names = FALSE)
  if (any(!c("name", "ID") %in% colnames(row_data))) {
    stop("'name' and/or 'ID' columns are not present. Run make_unique() and make_se() first.", call. = FALSE)
  }
  if (length(grep("_p.adj|_diff", colnames(row_data))) < 1) {
    stop("'[contrast]_diff' and/or '[contrast]_p.adj' columns are not present. Run test_diff() first.", call. = FALSE)
  }
  message("Background")
  if (!is.null(metadata(dep)$exp) && metadata(dep)$exp == "TMT" && !is.null(metadata(dep)$level) && metadata(dep)$level == "protein") {
    background <- unique(row_data$Gene)
  } else if (!is.null(metadata(dep)$exp) && metadata(dep)$exp == "TMT" && !is.null(metadata(dep)$level) && metadata(dep)$level == "gene") {
    background <- unique(row_data$ID)
  } else if (!is.null(metadata(dep)$level) && metadata(dep)$level == "protein") {
    background <- unique(gsub("[.].*", "", row_data$name))
  } else if (!is.null(metadata(dep)$level) && metadata(dep)$level %in% c("peptide", "site", "glycan")) {
    background <- unique(row_data$Gene)
  } else {
    background <- unique(gsub("[.].*", "", row_data$name))
  }
  background_enriched <- enrichr_mod(background, databases)
  df_background <- NULL
  for (db in databases) {
    temp <- background_enriched[db][[1]] %>% dplyr::mutate(var = db)
    df_background <- rbind(df_background, temp)
  }
  df_background$contrast <- "background"
  df_background$n <- length(background)
  OUT <- df_background %>%
    dplyr::mutate(bg_IN = as.numeric(gsub("/.*", "", Overlap)), bg_OUT = n - bg_IN) %>%
    dplyr::select(Term, bg_IN, bg_OUT)
  if (contrasts) {
    df <- row_data %>% as.data.frame() %>% dplyr::mutate(name = gsub("[.].*", "", name))
    contrast_cols <- df %>% dplyr::select(dplyr::ends_with("_significant")) %>% colnames()
    df_enrich <- NULL
    for (contrast in contrast_cols) {
      message(gsub("_significant", "", contrast))
      df[is.na(df[[contrast]]), contrast] <- FALSE
      significant <- df %>% dplyr::filter(!is.na(.[[gsub("_significant", "_diff", contrast)]]))
      if (direction == "UP") {
        significant <- significant %>% dplyr::filter(.[[gsub("_significant", "_diff", contrast)]] > log2_threshold)
      } else if (direction == "DOWN") {
        significant <- significant %>% dplyr::filter(.[[gsub("_significant", "_diff", contrast)]] < -log2_threshold)
      }
      significant <- significant %>% dplyr::filter(.[[gsub("_significant", "_p.adj", contrast)]] < alpha)
      if (!is.null(metadata(dep)$exp) && metadata(dep)$exp == "TMT" && !is.null(metadata(dep)$level) && metadata(dep)$level == "protein") {
        genes <- unique(significant$Gene)
      } else if (!is.null(metadata(dep)$exp) && metadata(dep)$exp == "TMT" && !is.null(metadata(dep)$level) && metadata(dep)$level == "gene") {
        genes <- unique(significant$ID)
      } else if (!is.null(metadata(dep)$level) && metadata(dep)$level == "protein") {
        genes <- significant$name
      } else if (!is.null(metadata(dep)$level) && metadata(dep)$level %in% c("peptide", "site", "glycan")) {
        genes <- unique(significant$Gene)
      } else {
        genes <- significant$name
      }
      message(paste0(length(genes), " genes are submitted"))
      if (length(genes) != 0) {
        enriched <- enrichr_mod(genes, databases)
        contrast_enrich <- NULL
        for (db in databases) {
          temp <- enriched[db][[1]] %>% dplyr::mutate(var = db)
          contrast_enrich <- rbind(contrast_enrich, temp)
        }
        if (nrow(contrast_enrich) != 0) {
          contrast_enrich$contrast <- contrast
          contrast_enrich$n <- length(genes)
          cat("Background correction... ")
          contrast_enrich <- contrast_enrich %>%
            dplyr::mutate(IN = as.numeric(gsub("/.*", "", Overlap)), OUT = n - IN) %>%
            dplyr::select(-n) %>%
            dplyr::left_join(OUT, by = "Term") %>%
            dplyr::mutate(log_odds = log2((IN * bg_OUT) / (OUT * bg_IN)))
          cat("Done.")
          contrast_enrich$contrast <- gsub("_significant", "", contrast_enrich$contrast)
          df_enrich <- rbind(df_enrich, contrast_enrich)
        }
      } else {
        cat("No significant genes for enrichment analysis")
      }
    }
  } else {
    significant <- row_data %>% as.data.frame() %>%
      dplyr::select(name, significant) %>% dplyr::filter(significant) %>%
      dplyr::mutate(name = gsub("[.].*", "", name))
    if (!is.null(metadata(dep)$exp) && metadata(dep)$exp == "TMT" && !is.null(metadata(dep)$level) && metadata(dep)$level == "protein") {
      genes <- unique(significant$Gene)
    } else if (!is.null(metadata(dep)$level) && metadata(dep)$level %in% c("peptide", "site", "glycan")) {
      genes <- unique(significant$Gene)
    } else {
      genes <- significant$name
    }
    enriched <- enrichr_mod(genes, databases)
    df_enrich <- NULL
    for (db in databases) {
      temp <- enriched[db][[1]] %>% dplyr::mutate(var = db)
      df_enrich <- rbind(df_enrich, temp)
    }
    df_enrich$contrast <- "significant"
    df_enrich$n <- length(genes)
    cat("Background correction... ")
    df_enrich <- df_enrich %>%
      dplyr::mutate(IN = as.numeric(gsub("/.*", "", Overlap)), OUT = n - IN) %>%
      dplyr::select(-n) %>%
      dplyr::left_join(OUT, by = "Term") %>%
      dplyr::mutate(log_odds = log2((IN * bg_OUT) / (OUT * bg_IN)))
    cat("Done.")
  }
  if (is.null(df_enrich) || nrow(df_enrich) == 0) return(NULL)
  df_enrich$p_hyper <- phyper(q = df_enrich$IN - 1, m = df_enrich$bg_IN, n = df_enrich$bg_OUT,
    k = df_enrich$IN + df_enrich$OUT, lower.tail = FALSE)
  df_enrich$p.adjust_hyper <- p.adjust(df_enrich$p_hyper, method = "BH")
  df_enrich
}

# or_test: over-representation test. Enrichr backend. From FragPipeAnalystR.
# database: friendly name (mapped below) or any Enrichr libraryName from
# https://maayanlab.cloud/Enrichr/datasetStatistics (passthrough).
or_test <- function(se, database = "GO Biological Process", backend = "enrichr",
  direction = "UP", log2_threshold = 0.7, alpha = 0.05) {
  if (backend == "enrichr") {
    database_mappings <- c(
      "GO Biological Process" = "GO_Biological_Process_2021",
      "GO Cellular Component" = "GO_Cellular_Component_2021",
      "GO Molecular Function" = "GO_Molecular_Function_2021",
      "Hallmark" = "MSigDB_Hallmark_2020",
      "KEGG" = "KEGG_2021_Human",
      "KEGG Mouse" = "KEGG_2019_Mouse",
      "Reactome" = "Reactome_2022",
      "WikiPathways Mouse" = "WikiPathways_2024_Mouse"
    )
    reverse_mappings <- c(
      "GO_Biological_Process_2021" = "GO Biological Process",
      "GO_Cellular_Component_2021" = "GO Cellular Component",
      "GO_Molecular_Function_2021" = "GO Molecular Function",
      "MSigDB_Hallmark_2020" = "Hallmark",
      "KEGG_2021_Human" = "KEGG",
      "KEGG_2019_Mouse" = "KEGG Mouse",
      "Reactome_2022" = "Reactome",
      "WikiPathways_2024_Mouse" = "WikiPathways Mouse"
    )
    enrichr_name <- if (database %in% names(database_mappings)) database_mappings[database] else database
    result <- test_ora_mod(se, databases = enrichr_name, contrasts = TRUE,
      direction = direction, log2_threshold = log2_threshold, alpha = alpha)
    if (!is.null(result)) {
      idx <- result$var %in% names(reverse_mappings)
      result$var[idx] <- reverse_mappings[result$var[idx]]
    }
    return(result)
  }
  cat("Only enrichr backend supported. Use backend='enrichr'.\n")
  NULL
}

# plot_or: bar plot of enrichment results. From FragPipeAnalystR.
plot_or <- function(or_result, number = 10, alpha = 0.05, contrasts = NULL, databases = NULL,
  adjust = FALSE, use_whole_proteome = FALSE, nrow = 1, term_size = 8) {
  assertthat::assert_that(is.data.frame(or_result), is.numeric(number), length(number) == 1,
    is.numeric(alpha), length(alpha) == 1, is.numeric(term_size), length(term_size) == 1,
    is.numeric(nrow), length(nrow) == 1)
  if (any(!c("Term", "var", "contrast", "Adjusted.P.value") %in% colnames(or_result))) {
    stop("or_result must contain Term, var, contrast, Adjusted.P.value. Ensure HGNC gene symbols in Gene column.", call. = FALSE)
  }
  no_enrichment_text <- "\n   No enrichment found.\n   You can still download enrichment result table.\n"
  if (!is.null(contrasts)) {
    valid_contrasts <- unique(or_result$contrast)
    if (!all(contrasts %in% valid_contrasts)) {
      return(ggplot2::ggplot() + ggplot2::annotate("text", x = 4, y = 25, size = 8, label = no_enrichment_text) + ggplot2::theme_void())
    }
    or_result <- dplyr::filter(or_result, contrast %in% contrasts)
  }
  if (!is.null(databases)) {
    valid_dbs <- unique(or_result$var)
    if (all(!databases %in% valid_dbs)) stop("Invalid database(s). Valid: ", paste(valid_dbs, collapse = ", "), call. = FALSE)
    or_result <- dplyr::filter(or_result, var %in% databases)
  }
  if (!use_whole_proteome) {
    if (adjust) {
      terms <- or_result %>% dplyr::group_by(contrast, var) %>%
        dplyr::filter(p.adjust_hyper <= alpha) %>% dplyr::arrange(p.adjust_hyper) %>%
        dplyr::slice(seq_len(number)) %>% dplyr::pull(Term)
      subset <- or_result %>% dplyr::filter(Term %in% terms) %>% dplyr::arrange(var, p.adjust_hyper)
    } else {
      terms <- or_result %>% dplyr::group_by(contrast, var) %>%
        dplyr::filter(p_hyper <= alpha) %>% dplyr::arrange(p_hyper) %>%
        dplyr::slice(seq_len(number)) %>% dplyr::pull(Term)
      subset <- or_result %>% dplyr::filter(Term %in% terms) %>% dplyr::arrange(var, p_hyper)
    }
  } else {
    if (adjust) {
      terms <- or_result %>% dplyr::group_by(contrast, var) %>%
        dplyr::filter(Adjusted.P.value <= alpha) %>% dplyr::arrange(Adjusted.P.value) %>%
        dplyr::slice(seq_len(number)) %>% dplyr::pull(Term)
      subset <- or_result %>% dplyr::filter(Term %in% terms) %>% dplyr::arrange(var, Adjusted.P.value)
    } else {
      terms <- or_result %>% dplyr::group_by(contrast, var) %>%
        dplyr::filter(P.value <= alpha) %>% dplyr::arrange(P.value) %>%
        dplyr::slice(seq_len(number)) %>% dplyr::pull(Term)
      subset <- or_result %>% dplyr::filter(Term %in% terms) %>% dplyr::arrange(var, P.value)
    }
  }
  subset$Term <- factor(subset$Term, levels = unique(subset$Term))
  subset$var <- factor(subset$var, levels = unique(subset$var))
  if (nrow(subset) == 0) {
    return(ggplot2::ggplot() + ggplot2::annotate("text", x = 4, y = 25, size = 8, label = no_enrichment_text) + ggplot2::theme_void())
  }
  if (!use_whole_proteome) {
    if (adjust) {
      ggplot2::ggplot(subset, ggplot2::aes(y = reorder(Term, log_odds), x = log_odds, size = IN, color = p.adjust_hyper)) +
        ggplot2::geom_point() + ggplot2::facet_wrap(~contrast, nrow = nrow) +
        ggplot2::scale_color_continuous(low = "red", high = "blue", name = "p.adjust",
          guide = ggplot2::guide_colorbar(reverse = TRUE, label.theme = ggplot2::element_text(angle = 90), label.vjust = 0.5)) +
        ggplot2::labs(y = "Term", x = "log2 Odds ratio", size = "size") +
        ggplot2::theme_bw() + ggplot2::theme(legend.position = "top", legend.text = ggplot2::element_text(size = 9))
    } else {
      ggplot2::ggplot(subset, ggplot2::aes(y = reorder(Term, log_odds), x = log_odds, size = IN, color = p_hyper)) +
        ggplot2::geom_point() + ggplot2::facet_wrap(~contrast, nrow = nrow) +
        ggplot2::scale_color_continuous(low = "red", high = "blue", name = "p",
          guide = ggplot2::guide_colorbar(reverse = TRUE, label.theme = ggplot2::element_text(angle = 90), label.vjust = 0.5)) +
        ggplot2::labs(y = "Term", x = "log2 Odds ratio", size = "size") +
        ggplot2::theme_bw() + ggplot2::theme(legend.position = "top", legend.text = ggplot2::element_text(size = 9))
    }
  } else {
    if (adjust) {
      ggplot2::ggplot(subset, ggplot2::aes(y = reorder(Term, Odds.Ratio), x = log2(Odds.Ratio), size = IN, color = Adjusted.P.value)) +
        ggplot2::geom_point() + ggplot2::facet_wrap(~contrast, nrow = nrow) +
        ggplot2::scale_color_continuous(low = "red", high = "blue", name = "Adjusted.P.value",
          guide = ggplot2::guide_colorbar(reverse = TRUE, label.theme = ggplot2::element_text(angle = 90), label.vjust = 0.5)) +
        ggplot2::labs(y = "Term", x = "log2 Odds ratio", size = "size") +
        ggplot2::theme_bw() + ggplot2::theme(legend.position = "top", legend.text = ggplot2::element_text(size = 9))
    } else {
      ggplot2::ggplot(subset, ggplot2::aes(y = reorder(Term, Odds.Ratio), x = log2(Odds.Ratio), size = IN, color = P.value)) +
        ggplot2::geom_point() + ggplot2::facet_wrap(~contrast, nrow = nrow) +
        ggplot2::scale_color_continuous(low = "red", high = "blue", name = "P.value",
          guide = ggplot2::guide_colorbar(reverse = TRUE, label.theme = ggplot2::element_text(angle = 90), label.vjust = 0.5)) +
        ggplot2::labs(y = "Term", x = "log2 Odds ratio", size = "size") +
        ggplot2::theme_bw() + ggplot2::theme(legend.position = "top", legend.text = ggplot2::element_text(size = 9))
    }
  }
}
