# Meaningful synthetic checks for vector coordinates and codon normalization.
# Run from repository root with: Rscript check_R_calculations.R
source("R_functions.R")
work <- tempfile("r_helper_test_", tmpdir = ".")
dir.create(work)
on.exit_cleanup <- function() unlink(work, recursive = TRUE)

run_tests <- function() {
  on.exit(on.exit_cleanup(), add = TRUE)
  sample_dir <- file.path(work, "test_sample")
  dir.create(sample_dir)
  file <- file.path(sample_dir, "gene.txt")
  write.table(c(100, 1, 2, 0, 3, 0, 0, 200), file,
              row.names = FALSE, col.names = FALSE, quote = FALSE)
  inputs <- list(
    vector_dir = work,
    config = list(r_analysis = list(candidate_stop_codons = c("UAA","UAG","AGA","AGG"))),
    samples = data.frame(folder = "test_sample", sample_id = "test_sample",
      Condition = "WT", Replicate = "1", Sample = "WT1"),
    mapping = data.frame(Gene = "test_gene", vector_file = "gene.txt",
      sequence_gene = "test_gene", vector_start_1based = 2L)
  )
  result <- collect_occupancy(inputs, list(test_gene = "ATGTAA"))
  stopifnot(nrow(result$data) == 2L,
            identical(result$data$codon, c("ATG","TAA")),
            all(result$data$Percent == c(50,50)),
            identical(result$data$is_stop, c(FALSE,TRUE)),
            result$audit$UpstreamPositionsExcluded == 1,
            result$audit$DownstreamPositionsExcluded == 1,
            result$gene_qc$TotalReads == 6,
            nrow(result$nucleotide_data) == 6,
            sum(result$nucleotide_data$Count) == 6,
            abs(sum(result$nucleotide_data$Percent) - 100) < 1e-8)
  # Incomplete codons must not recycle bases into a fabricated codon.
  partial <- collect_occupancy(inputs, list(test_gene = "ATGTA"))
  stopifnot(nrow(partial$data) == 1,
            partial$audit$PartialCodonPositions == 2,
            partial$audit$PartialCodonCounts == 3)
  write.table(rep(0, 8), file, row.names = FALSE, col.names = FALSE, quote = FALSE)
  zero <- collect_occupancy(inputs, list(test_gene = "ATGTAA"))
  stopifnot(all(is.na(zero$data$Percent)))
  # An out-of-range sequence region must fail rather than pad the vector.
  bad <- tryCatch(collect_occupancy(inputs, list(test_gene = "ATGTAAATGTAA")),
                  error = function(e) e)
  stopifnot(inherits(bad, "error"))
  set.seed(17)
  noise <- runif(180, -0.1, 0.1)
  three <- 2 + sin(2 * pi * (0:179) / 3) + noise
  four <- 2 + sin(2 * pi * (0:179) / 4) + noise
  stopifnot(fourier_3nt_score(three) > 100 * fourier_3nt_score(four))
  message("R helper tests passed.")
}
run_tests()
