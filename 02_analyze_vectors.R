#!/usr/bin/env Rscript
# Run from the directory containing these files:
#   Rscript 02_analyze_vectors.R HEK
#   Rscript 02_analyze_vectors.R N2a
# For Source in RStudio, change the following choice if needed:
ANALYSIS <- "HEK"
args <- commandArgs(trailingOnly = TRUE)
if (length(args)) ANALYSIS <- args[1]
if (length(args) > 1L || !ANALYSIS %in% c("HEK", "N2a", "HEK_unprocessed"))
  stop("Choose HEK, N2a or HEK_unprocessed: Rscript 02_analyze_vectors.R HEK")
settings <- paste0(ANALYSIS, "_settings.json")
if (!file.exists(settings) || !file.exists("R_functions.R"))
  stop("Run from the directory containing the scripts and settings files.")
needed <- c("jsonlite", "dplyr", "ggplot2")
if (ANALYSIS == "HEK_unprocessed") needed <- c(needed, "Biostrings", "rtracklayer")
missing <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Install missing packages with install_R_packages.R: ", paste(missing, collapse = ", "))
library(dplyr)
library(ggplot2)
Sys.setenv(MITORIBO_ROOT = normalizePath(getwd(), mustWork = TRUE))
source("R_functions.R")
inputs <- load_analysis_inputs(settings)
print(inputs$samples[c("sample_id", "Sample", "Condition", "Replicate")])
run_vector_analysis(inputs)
