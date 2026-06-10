#!/usr/bin/env Rscript
# GeneLab script for building organism-specific UniProt proteome FASTA (targets only).
#
# Usage:
#   Rscript GL-DPPD-7110-A_build-reference-proteome.R 'Organism' [annotations.csv] philosopher_bin
#
# Requires `uniprot_id` in GL-DPPD-7110-A_annotations.csv.
# Runs: philosopher database --id UP... --reviewed --nodecoys
#   (Swiss-Prot canonical / reference proteome; no isoforms; no decoys/contaminants — add at search time)
#
# Outputs (cwd):
#   YYYY-MM-DD-reviewed-UP....-{N}.fas  (Philosopher name + entry count; --isoform off)
#   YYYY-MM-DD-reviewed-UP....-{N}-GL-build-info.txt
options(timeout = 3600)

GL_DPPD_ID <- "GL-DPPD-7110-A"
ref_tab_path <- "https://raw.githubusercontent.com/nasa/GeneLab_Data_Processing/master/GeneLab_Reference_Annotations/Pipeline_GL-DPPD-7110_Versions/GL-DPPD-7110-A/GL-DPPD-7110-A_annotations.csv"

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 1) {
  stop(
    "Usage: Rscript GL-DPPD-7110-A_build-reference-proteome.R 'Organism' [annotations.csv] philosopher_bin"
  )
}

target_organism <- args[1]
ref_tab_path <- if (length(args) >= 2 && nzchar(args[2])) args[2] else ref_tab_path
philosopher_bin <- if (length(args) >= 3 && nzchar(args[3])) {
  args[3]
} else if (nzchar(Sys.getenv("PHILOSOPHER"))) {
  Sys.getenv("PHILOSOPHER")
} else {
  stop("Philosopher path required as 3rd argument or PHILOSOPHER env var.")
}
if (!file.exists(philosopher_bin)) {
  stop("Philosopher not found: ", philosopher_bin)
}

ref_table <- read.csv(ref_tab_path, stringsAsFactors = FALSE)
if (!"uniprot_id" %in% names(ref_table)) {
  stop("Reference table missing 'uniprot_id' column: ", ref_tab_path)
}

target_info <- ref_table[ref_table$species == target_organism, , drop = FALSE]
if (nrow(target_info) != 1) {
  stop("Expected one row for '", target_organism, "' in reference table.")
}

uniprot_id <- trimws(target_info$uniprot_id)
if (is.na(uniprot_id) || !nzchar(uniprot_id)) {
  cat("\n  No uniprot_id for '", target_organism, "' in reference table. Exiting.\n\n", sep = "")
  quit(status = 0)
}

run_philosopher <- function(bin, args) {
  status <- system2(bin, args, stdout = "", stderr = "")
  if (!identical(status, 0L)) {
    stop("Philosopher failed: ", bin, " ", paste(args, collapse = " "))
  }
}

build_philosopher_fasta <- function(bin, up_id, use_reviewed) {
  run_philosopher(bin, c("workspace", "--clean"))
  run_philosopher(bin, c("workspace", "--init"))
  db_args <- c("database", "--id", up_id)
  if (use_reviewed) {
    db_args <- c(db_args, "--reviewed")
  }
  db_args <- c(db_args, "--nodecoys")
  run_philosopher(bin, db_args)
  fasta_files <- list.files(getwd(), pattern = "\\.fas$", full.names = TRUE)
  fasta_files <- fasta_files[grepl(up_id, basename(fasta_files), fixed = TRUE)]
  if (length(fasta_files) == 0) {
    stop("Philosopher did not write a .fas for ", up_id)
  }
  list(
    path = fasta_files[[which.max(file.info(fasta_files)$mtime)]],
    command = paste("philosopher", paste(db_args, collapse = " "))
  )
}

use_reviewed <- TRUE
build <- tryCatch(
  build_philosopher_fasta(philosopher_bin, uniprot_id, use_reviewed),
  error = function(e) {
    if (!use_reviewed) stop(e)
    use_reviewed <<- FALSE
    build_philosopher_fasta(philosopher_bin, uniprot_id, FALSE)
  }
)
philosopher_fasta <- build$path
philosopher_command <- build$command

headers <- readLines(philosopher_fasta, warn = FALSE)
headers <- headers[startsWith(headers, ">")]
n_total <- length(headers)
if (use_reviewed && n_total < 100) {
  file.remove(philosopher_fasta)
  use_reviewed <- FALSE
  build <- build_philosopher_fasta(philosopher_bin, uniprot_id, FALSE)
  philosopher_fasta <- build$path
  philosopher_command <- build$command
  headers <- readLines(philosopher_fasta, warn = FALSE)
  headers <- headers[startsWith(headers, ">")]
  n_total <- length(headers)
}
n_sp <- sum(grepl("^>sp\\|", headers))
n_tr <- sum(grepl("^>tr\\|", headers))

out_fasta_basename <- sub(
  "\\.fas$",
  sprintf("-%d.fas", n_total),
  basename(philosopher_fasta)
)
out_fasta_filename <- file.path(dirname(philosopher_fasta), out_fasta_basename)
if (!identical(philosopher_fasta, out_fasta_filename)) {
  file.rename(philosopher_fasta, out_fasta_filename)
}
out_log_filename <- sub("\\.fas$", "-GL-build-info.txt", out_fasta_basename)

if (file.exists(out_log_filename)) {
  cat("\n  '", out_log_filename, "' exists already. Move it to rebuild.\n", sep = "")
  quit(status = 0)
}

philosopher_version <- paste(system2(philosopher_bin, "version", stdout = TRUE, stderr = TRUE), collapse = " ")
date_generated <- format(Sys.time(), "%d-%B-%Y")

writeLines(paste("Based on:\n    ", GL_DPPD_ID), out_log_filename)
write(paste("\nBuild done on:\n    ", date_generated, sep = ""), out_log_filename, append = TRUE)
write(paste("\nTarget organism:\n    ", target_organism, sep = ""), out_log_filename, append = TRUE)
write(paste("\nUniProt proteome ID:\n    ", uniprot_id, sep = ""), out_log_filename, append = TRUE)
write(paste("\nPhilosopher command:\n    ", philosopher_command, sep = ""), out_log_filename, append = TRUE)
write(paste("\nOutput FASTA:\n    ", out_fasta_basename, sep = ""), out_log_filename, append = TRUE)
write(paste("\nSwiss-Prot (reviewed):\n    ", if (use_reviewed) "yes" else "no (auto: <100 reviewed entries;  dropped `--reviewed` flag)", sep = ""), out_log_filename, append = TRUE)
write("\nIsoforms included:\n    no", out_log_filename, append = TRUE)
write(paste("\nTotal FASTA entries ({N} in filename):\n    ", n_total, sep = ""), out_log_filename, append = TRUE)
write(paste("\nSwiss-Prot (sp|) entries:\n    ", n_sp, sep = ""), out_log_filename, append = TRUE)
write(paste("\nTrEMBL (tr|) entries:\n    ", n_tr, sep = ""), out_log_filename, append = TRUE)
write(paste("\nUsed Philosopher version:\n    ", philosopher_version, sep = ""), out_log_filename, append = TRUE)
write("\n\nAll session info:\n", out_log_filename, append = TRUE)
write(capture.output(sessionInfo()), out_log_filename, append = TRUE)

cat("\nFASTA:     ", out_fasta_basename, "\n", sep = "")
cat("Build log: ", out_log_filename, "\n", sep = "")
