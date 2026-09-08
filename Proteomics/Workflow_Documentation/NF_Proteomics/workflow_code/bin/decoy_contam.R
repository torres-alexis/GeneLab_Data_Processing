# MSstats comparison scrub (Issue / non-finite FC / NA p / DF≤0 → NA p/adj.p).

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
