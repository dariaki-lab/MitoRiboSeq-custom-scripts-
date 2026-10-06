#!/usr/bin/env Rscript
# Optional inventory only. The public analysis does not depend on a saved workspace.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) stop("Usage: Rscript scripts/inspect_r_workspace.R path/to/RData inventory.tsv")
workspace <- new.env(parent = emptyenv())
objects <- load(args[1], envir = workspace)
inventory <- data.frame(
  Object = objects,
  Class = vapply(objects, function(n) paste(class(get(n, workspace)), collapse = "/"), character(1)),
  Bytes = vapply(objects, function(n) as.numeric(object.size(get(n, workspace))), numeric(1)),
  stringsAsFactors = FALSE
)
write.table(inventory, args[2], sep = "\t", quote = FALSE, row.names = FALSE)
