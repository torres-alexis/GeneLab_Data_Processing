# Shared decoy / contaminant filters + MSstats comparison scrub.

DECOY_PREFIXES_DEFAULT <- c("rev_", "REV_", "decoy_", "DECOY_")
CONTAM_PREFIXES_DEFAULT <- c("contam_", "Cont_", "CONTAM_", "CRAP_", "crap_")
ID_COLUMNS_DEFAULT <- c(
  "ProteinName", "Protein ID", "Protein.ID", "ProteinID", "Protein",
  "Entry Name", "Entry.Name", "Index", "name", "ID"
)

.dc_env_flag <- function(name, default = TRUE) {
  raw <- Sys.getenv(name, unset = if (isTRUE(default)) "true" else "false")
  tolower(trimws(raw)) %in% c("1", "true", "yes", "y")
}

.dc_env_str <- function(name, default = "rev_") {
  raw <- trimws(Sys.getenv(name, unset = default))
  if (!nzchar(raw)) default else raw
}

.dc_group_parts <- function(value) {
  if (is.null(value) || length(value) == 0 || is.na(value)) return(character(0))
  text <- trimws(as.character(value))
  if (!nzchar(text) || tolower(text) %in% c("na", "nan", "none")) return(character(0))
  parts <- unlist(strsplit(text, "[;,]+"))
  parts <- trimws(parts)
  parts[nzchar(parts)]
}

token_is_junk <- function(token,
                          decoy_prefixes = DECOY_PREFIXES_DEFAULT,
                          contam_prefixes = CONTAM_PREFIXES_DEFAULT) {
  t <- trimws(as.character(token))
  if (!nzchar(t)) return(FALSE)
  for (p in decoy_prefixes) {
    if (startsWith(t, p) || startsWith(tolower(t), tolower(p))) return(TRUE)
  }
  for (p in contam_prefixes) {
    if (startsWith(t, p) || startsWith(tolower(t), tolower(p))) return(TRUE)
  }
  FALSE
}

value_is_junk <- function(value, ...) {
  parts <- .dc_group_parts(value)
  if (!length(parts)) return(FALSE)
  if (all(vapply(parts, token_is_junk, logical(1), ...))) return(TRUE)
  token_is_junk(parts[[1]], ...)
}

drop_decoy_contam_rows <- function(df,
                                   id_cols = NULL,
                                   decoy_prefix = "rev_") {
  if (is.null(df) || !nrow(df)) return(df)
  cols <- id_cols
  if (is.null(cols)) {
    cols <- intersect(ID_COLUMNS_DEFAULT, names(df))
  }
  cols <- intersect(cols, names(df))
  if (!length(cols)) return(df)
  junk <- rep(FALSE, nrow(df))
  prefix <- trimws(as.character(decoy_prefix)[[1]])
  if (!nzchar(prefix)) prefix <- "rev_"
  kwargs <- list(decoy_prefixes = prefix)
  for (col in cols) {
    junk <- junk | vapply(df[[col]], function(v) do.call(value_is_junk, c(list(v), kwargs)), logical(1))
  }
  n_drop <- sum(junk, na.rm = TRUE)
  if (n_drop) {
    message("decoy_contam: dropped ", n_drop, " / ", nrow(df), " rows")
  }
  df[!junk, , drop = FALSE]
}

scrub_msstats_comparison <- function(df) {
  if (is.null(df) || !nrow(df)) return(df)
  n <- nrow(df)
  issue_col <- intersect(c("issue", "Issue"), names(df))
  bad <- rep(FALSE, n)
  if (length(issue_col)) {
    iss <- df[[issue_col[1]]]
    bad <- bad | (!is.na(iss) & nzchar(trimws(as.character(iss))))
  }
  fc_col <- intersect(c("log2FC", "logFC"), names(df))
  if (length(fc_col)) {
    bad <- bad | !is.finite(df[[fc_col[1]]])
  }
  if ("pvalue" %in% names(df)) {
    bad <- bad | is.na(df$pvalue)
  }
  if ("DF" %in% names(df)) {
    bad <- bad | is.na(df$DF) | df$DF <= 0
  }
  if ("pvalue" %in% names(df)) df$pvalue[bad] <- NA_real_
  if ("adj.pvalue" %in% names(df)) df$adj.pvalue[bad] <- NA_real_
  n_bad <- sum(bad, na.rm = TRUE)
  if (n_bad) {
    message("scrub_msstats_comparison: NA'd p/adj.p on ", n_bad, " failed-test row(s)")
  }
  df
}
