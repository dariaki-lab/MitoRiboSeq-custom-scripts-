# Parse all active R code and all Rmd chunks. No analysis dependencies needed.
files <- c(list.files("scripts", pattern = "\\.R$", full.names = TRUE),
           list.files("analysis", pattern = "\\.R$", recursive = TRUE, full.names = TRUE),
           list.files("tests", pattern = "\\.R$", full.names = TRUE))
for (file in files) parse(file = file)
notebooks <- list.files("analysis", pattern = "\\.Rmd$", recursive = TRUE, full.names = TRUE)
for (file in notebooks) {
  lines <- readLines(file, warn = FALSE)
  inside <- FALSE
  code <- character()
  for (line in lines) {
    if (grepl("^```\\{r", line)) {
      if (inside) stop("Nested code fence in ", file)
      inside <- TRUE
    } else if (inside && grepl("^```[[:space:]]*$", line)) {
      inside <- FALSE
    } else if (inside) {
      code <- c(code, line)
    }
  }
  if (inside) stop("Unclosed R fence in ", file)
  parse(text = code)
}
message("All active R files and notebook chunks parse.")
