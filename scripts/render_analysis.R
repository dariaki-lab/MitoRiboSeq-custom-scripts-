#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
folders <- c(hek_processed = "hek", n2a_processed = "n2a", hek_unprocessed = "hek_unprocessed")
if (length(args) != 1L || !args[1] %in% names(folders))
  stop("Run from repository root: Rscript scripts/render_analysis.R hek_processed | n2a_processed | hek_unprocessed")
root <- normalizePath(getwd(), mustWork = TRUE)
if (!file.exists(file.path(root, "metadata", "sample_info.csv"))) stop("Run from the repository root.")
needed <- c("rmarkdown", "jsonlite", "tidyverse", "knitr")
if (args[1] == "hek_unprocessed") needed <- c(needed, "Biostrings", "rtracklayer")
missing <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Missing R packages: ", paste(missing, collapse = ", "))
if (!rmarkdown::pandoc_available()) stop("Pandoc is required; use RStudio or install Pandoc.")
Sys.setenv(MITORIBO_ROOT = root)
config <- jsonlite::fromJSON(file.path(root, "config", paste0(args[1], ".json")))
report_dir <- file.path(root, config$output_dir, "R")
dir.create(report_dir, recursive = TRUE, showWarnings = FALSE)
rmarkdown::render(file.path(root, "analysis", folders[[args[1]]], "analysis.Rmd"),
                  output_file = paste0(args[1], "_analysis.html"), output_dir = report_dir,
                  knit_root_dir = root, envir = new.env(parent = globalenv()))
