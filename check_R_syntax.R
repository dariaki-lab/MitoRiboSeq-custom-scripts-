# Run from the script directory: Rscript check_R_syntax.R
for (file in list.files(pattern = "\\.R$")) parse(file = file)
message("All R files parse.")
