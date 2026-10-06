#!/usr/bin/env Rscript
# Run from the repository root. Optional argument: --with-unprocessed
args <- commandArgs(trailingOnly = TRUE)
cran <- c("tidyverse", "jsonlite", "rmarkdown", "knitr")
missing <- cran[!vapply(cran, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) install.packages(missing, repos = "https://cloud.r-project.org")
if ("--with-unprocessed" %in% args) {
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager", repos = "https://cloud.r-project.org")
  bio <- c("Biostrings", "rtracklayer")
  missing <- bio[!vapply(bio, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) BiocManager::install(missing, ask = FALSE, update = FALSE)
}
