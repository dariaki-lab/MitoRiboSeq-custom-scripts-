#!/usr/bin/env Rscript
# Processed HEK and N2a: Rscript install_R_packages.R
# Also install optional unprocessed dependencies:
#   Rscript install_R_packages.R --with-unprocessed
cran <- c("jsonlite", "dplyr", "ggplot2")
missing <- cran[!vapply(cran, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) install.packages(missing, repos = "https://cloud.r-project.org")
if ("--with-unprocessed" %in% commandArgs(trailingOnly = TRUE)) {
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager", repos = "https://cloud.r-project.org")
  bio <- c("Biostrings", "rtracklayer")
  missing <- bio[!vapply(bio, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) BiocManager::install(missing, ask = FALSE, update = FALSE)
}
