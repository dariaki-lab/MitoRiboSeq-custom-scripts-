# Input handling and calculations for 02_analyze_vectors.R.
# No R workspace is loaded. No vectors are silently padded, resized, or skipped.

analysis_root <- function() {
  root <- Sys.getenv("MITORIBO_ROOT", unset = getwd())
  normalizePath(root, mustWork = TRUE)
}

repository_path <- function(root, path) {
  if (!is.character(path) || length(path) != 1L || !nzchar(path))
    stop("A required input path is empty.")
  if (grepl("^(/|[A-Za-z]:)", path)) path else file.path(root, path)
}

load_analysis_inputs <- function(settings_file) {
  root <- analysis_root()
  cfg <- jsonlite::fromJSON(repository_path(root, settings_file))
  if (!isTRUE(cfg$settings_confirmed) || !isTRUE(cfg$r_analysis$coordinates_confirmed)) {
    stop("Review the analysis settings and vector coordinates, then set settings_confirmed ",
         "and r_analysis.coordinates_confirmed to true in the profile JSON.")
  }
  if (!cfg$mapping_site %in% c("P", "A", "custom"))
    stop("Record mapping_site as P, A or custom for the supplied vectors.")
  sheet <- repository_path(root, cfg$sample_sheet)
  read_sheet <- if (grepl("\\.csv$", sheet, ignore.case = TRUE)) read.csv else read.delim
  samples <- read_sheet(sheet, stringsAsFactors = FALSE, check.names = FALSE,
                        colClasses = "character", na.strings = character())
  needed <- c("sample_id", "cell_line", "condition", "replicate", "sample_label", "vector_folder")
  if (!all(needed %in% names(samples))) stop("Missing columns in sample_info.csv")
  samples <- samples[samples$cell_line == cfg$cell_line, , drop = FALSE]
  if (!nrow(samples)) stop("No samples for configured cell line: ", cfg$cell_line)
  if (anyDuplicated(samples$sample_id)) stop("Duplicate sample_id in sample sheet")
  samples$folder <- ifelse(nzchar(samples$vector_folder), samples$vector_folder, samples$sample_id)
  samples$Condition <- factor(samples$condition, levels = cfg$r_analysis$condition_order)
  if (anyNA(samples$Condition)) stop("A condition is missing from r_analysis.condition_order")
  samples$Clone <- if ("clone" %in% names(samples)) samples$clone else ""
  samples$Replicate <- samples$replicate
  samples$Batch <- if ("batch" %in% names(samples)) samples$batch else ""
  # Numeric labels and provenance batches are separate; independence is not inferred.
  samples$Sample <- samples$sample_label
  if (anyDuplicated(samples$Sample)) stop("Sample plotting labels are not unique")
  mapping <- read.delim(repository_path(root, cfg$r_analysis$vector_mapping),
                        stringsAsFactors = FALSE, check.names = FALSE)
  if (anyDuplicated(mapping$Gene) || anyDuplicated(mapping$vector_file))
    stop("Duplicate gene or vector filename in transcript mapping")
  starts <- suppressWarnings(as.numeric(mapping$vector_start_1based))
  if (anyNA(starts) || any(!is.finite(starts)) || any(starts < 1 | starts != floor(starts)))
    stop("Confirm every gene-to-vector mapping and fill vector_start_1based in the transcript TSV.")
  mapping$vector_start_1based <- as.integer(starts)
  list(root = root, config = cfg, samples = samples, mapping = mapping,
       vector_dir = if (is.null(cfg$r_analysis$vector_dir)) repository_path(root, file.path(cfg$output_dir, "vectors"))
                    else repository_path(root, cfg$r_analysis$vector_dir),
       phasing_dir = repository_path(root, file.path(cfg$output_dir, "phasing")))
}

read_fasta_list <- function(path) {
  lines <- readLines(path, warn = FALSE)
  headers <- which(substr(lines, 1, 1) == ">")
  if (!length(headers) || headers[1] != 1) stop("Invalid FASTA: ", path)
  keys <- sub("^>([^[:space:]]+).*", "\\1", lines[headers])
  if (anyDuplicated(keys)) stop("Duplicate FASTA identifiers")
  ends <- c(headers[-1] - 1L, length(lines))
  values <- lapply(seq_along(headers), function(i) {
    if (ends[i] <= headers[i]) stop("Empty FASTA sequence")
    toupper(paste0(lines[(headers[i] + 1L):ends[i]], collapse = ""))
  })
  stats::setNames(values, keys)
}

extract_cds_sequences <- function(fasta, annotation, gene_column, mapping) {
  genome <- Biostrings::readDNAStringSet(fasta)
  names(genome) <- sub("[[:space:]].*$", "", names(genome))
  if (anyDuplicated(names(genome))) stop("Duplicate genome FASTA identifiers")
  gff <- rtracklayer::import(annotation)
  cds <- gff[gff$type == "CDS"]
  if (!gene_column %in% colnames(S4Vectors::mcols(cds)))
    stop("Configured gene column missing from GFF3: ", gene_column)
  gene_ids <- as.character(S4Vectors::mcols(cds)[[gene_column]])
  sequences <- list()
  for (i in seq_len(nrow(mapping))) {
    key <- mapping$sequence_gene[i]
    features <- cds[gene_ids == key & !is.na(gene_ids)]
    if (!length(features)) stop("Missing CDS feature for sequence gene: ", key)
    coordinates <- as.data.frame(features)
    chromosomes <- unique(as.character(coordinates$seqnames))
    strands <- unique(as.character(coordinates$strand))
    if (length(chromosomes) != 1 || length(strands) != 1 || !strands %in% c("+", "-"))
      stop("Mixed chromosome or strand for ", key)
    if (!chromosomes %in% names(genome)) stop("FASTA/GFF3 chromosome names disagree for ", key)
    # The supplied workflow targets mitochondrial CDSs. Additional CDS segments
    # need an explicit phase-aware extraction rule rather than implicit concatenation.
    if (length(features) != 1) stop("Multiple CDS segments need an explicit phase-aware rule: ", key)
    if ("phase" %in% colnames(S4Vectors::mcols(features))) {
      phase <- as.character(S4Vectors::mcols(features)$phase)
      if (!is.na(phase) && !phase %in% c("0", ".")) stop("Nonzero CDS phase: ", key)
    }
    piece <- Biostrings::subseq(genome[[chromosomes]],
                               start = coordinates$start,
                               end = coordinates$end)
    if (strands == "-") piece <- Biostrings::reverseComplement(piece)
    sequences[[mapping$Gene[i]]] <- as.character(piece)
  }
  sequences
}

read_count_vector <- function(path) {
  if (!file.exists(path)) stop("Missing count vector: ", path)
  data <- read.delim(path, header = FALSE, comment.char = "#")
  if (ncol(data) != 1L || !nrow(data)) stop("Expected one nonempty count column: ", path)
  counts <- suppressWarnings(as.numeric(data[[1]]))
  if (anyNA(counts) || any(!is.finite(counts)) || any(counts < 0))
    stop("Counts must be finite and nonnegative: ", path)
  counts
}

fourier_3nt_score <- function(counts) {
  if (length(counts) < 60L) return(NA_real_)
  x <- counts - mean(counts)
  power <- Mod(fft(x))^2
  # R indexes from 1; index 1 is the zero-frequency component.
  target <- round(length(power) / 3) + 1L
  background <- stats::median(power[-1])
  if (!is.finite(background) || background <= 0) return(NA_real_)
  power[target] / background
}

collect_occupancy <- function(inputs, sequences) {
  tables <- list()
  audits <- list()
  qc <- list()
  nucleotides <- list()
  candidate_stops <- inputs$config$r_analysis$candidate_stop_codons
  k <- 0L
  for (j in seq_len(nrow(inputs$mapping))) {
    mapping <- inputs$mapping[j, ]
    gene <- mapping$Gene
    sequence <- sequences[[gene]]
    if (is.null(sequence) || length(sequence) != 1L || !grepl("^[ACGTN]+$", sequence))
      stop("Missing or invalid sequence for gene: ", gene)
    bases <- strsplit(sequence, "")[[1]]
    n_sequence <- length(bases)
    n_codons <- n_sequence %/% 3L
    if (n_codons < 1) stop("Sequence shorter than one codon: ", gene)
    for (i in seq_len(nrow(inputs$samples))) {
      sample <- inputs$samples[i, ]
      path <- file.path(inputs$vector_dir, sample$folder, mapping$vector_file)
      vector <- read_count_vector(path)
      first <- mapping$vector_start_1based
      last <- first + n_sequence - 1L
      if (last > length(vector)) stop("Confirmed sequence region exceeds vector length: ", path,
                                      ". Review the reference and vector coordinates.")
      counts <- vector[seq.int(first, last)]
      total <- sum(counts)
      use <- seq_len(n_codons * 3L)
      codons <- matrix(bases[use], ncol = 3L, byrow = TRUE)
      sums <- rowSums(matrix(counts[use], ncol = 3L, byrow = TRUE))
      dna <- apply(codons, 1L, paste0, collapse = "")
      rna <- chartr("T", "U", dna)
      k <- k + 1L
      tables[[k]] <- data.frame(
        codon_pos = seq_len(n_codons), Sum = sums,
        A = codons[, 1], B = codons[, 2], C = codons[, 3],
        codon = dna, codon_rna = rna,
        Percent = if (total > 0) sums / total * 100 else rep(NA_real_, n_codons),
        Condition = sample$Condition, Replicate = sample$Replicate,
        Sample = sample$Sample, sample_id = sample$sample_id, Gene = gene,
        is_stop = rna %in% candidate_stops, stringsAsFactors = FALSE)
      nucleotides[[k]] <- data.frame(
        nucleotide_pos = seq_len(n_sequence), nucleotide = bases, Count = counts,
        Percent = if (total > 0) counts / total * 100 else rep(NA_real_, n_sequence),
        Condition = sample$Condition, Replicate = sample$Replicate,
        Sample = sample$Sample, sample_id = sample$sample_id, Gene = gene,
        stringsAsFactors = FALSE)
      audits[[k]] <- data.frame(sample_id = sample$sample_id, Gene = gene,
        VectorLength = length(vector), SequenceLength = n_sequence,
        VectorStart = first, VectorEnd = last, UpstreamPositionsExcluded = first - 1L,
        DownstreamPositionsExcluded = length(vector) - last,
        PartialCodonPositions = n_sequence %% 3L,
        PartialCodonCounts = if (n_sequence %% 3L) sum(tail(counts, n_sequence %% 3L)) else 0,
        TotalCountsInNormalizationRegion = total)
      frames <- vapply(0:2, function(frame) sum(counts[(seq_along(counts) - 1L) %% 3L == frame]), numeric(1))
      qc[[k]] <- data.frame(Sample = sample$Sample, sample_id = sample$sample_id,
        Condition = sample$Condition, Gene = gene, TotalReads = total,
        Frame0 = if (total > 0) 100 * frames[1] / total else NA_real_,
        Frame1 = if (total > 0) 100 * frames[2] / total else NA_real_,
        Frame2 = if (total > 0) 100 * frames[3] / total else NA_real_,
        FrameBias = if (total > 0) frames[1] / max((frames[2] + frames[3]) / 2, 1) else NA_real_,
        Fourier3nt = fourier_3nt_score(counts))
    }
  }
  list(data = do.call(rbind, tables), nucleotide_data = do.call(rbind, nucleotides), audit = do.call(rbind, audits),
       gene_qc = do.call(rbind, qc))
}

read_phasing <- function(inputs) {
  samples <- inputs$samples
  batches <- inputs$config$r_analysis$phasing_batches
  if (length(batches)) samples <- samples[samples$batch %in% batches, , drop = FALSE]
  if (!nrow(samples)) stop("No samples selected for phasing QC.")
  files <- file.path(inputs$phasing_dir, paste0(samples$sample_id, "_phasing.txt"))
  if (!any(file.exists(files))) {
    message("No phase_by_size tables found; phasing report skipped.")
    return(NULL)
  }
  if (!all(file.exists(files))) stop("Some expected phasing files are missing.")
  tables <- lapply(seq_along(files), function(i) {
    table <- read.delim(files[i], comment.char = "#")
    needed <- c("read_length", "reads_counted", "phase0", "phase1", "phase2")
    if (!all(needed %in% names(table))) stop("Unexpected phase_by_size columns: ", files[i])
    table$Sample <- samples$Sample[i]
    table$sample_id <- samples$sample_id[i]
    table$Condition <- samples$Condition[i]
    table
  })
  do.call(rbind, tables)
}


run_vector_analysis <- function(inputs) {
  OUTPUT_DIR <- repository_path(inputs$root, file.path(inputs$config$output_dir, "analysis"))
  dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
  # Export occupancy tables, QC and plots from the selected count vectors.
  # Tables always retain sample_id and the requested sample_label.
  set.seed(1)
  cfg <- inputs$config
  ra <- cfg$r_analysis
  if (nzchar(ra$sequence_fasta)) {
    raw_sequences <- read_fasta_list(repository_path(inputs$root, ra$sequence_fasta))
    seqs <- setNames(lapply(inputs$mapping$sequence_gene, function(g) raw_sequences[[g]]),
                     inputs$mapping$Gene)
  } else {
    seqs <- extract_cds_sequences(repository_path(inputs$root, ra$genome_fasta),
                                  repository_path(inputs$root, ra$cds_gff3),
                                  ra$gene_column, inputs$mapping)
  }
  result <- collect_occupancy(inputs, seqs)
  all_data <- result$data
  all_data$Condition <- factor(all_data$Condition, levels = ra$condition_order)
  all_data$Replicate <- factor(all_data$Replicate, levels = unique(inputs$samples$replicate))
  write.csv(inputs$samples, file.path(OUTPUT_DIR, "samples_used.csv"), row.names = FALSE)
  write.csv(all_data, file.path(OUTPUT_DIR, "codon_occupancy.csv"), row.names = FALSE)
  write.csv(result$nucleotide_data, file.path(OUTPUT_DIR, "nucleotide_counts.csv"), row.names = FALSE)
  write.csv(result$audit, file.path(OUTPUT_DIR, "vector_coordinate_audit.csv"), row.names = FALSE)
  write.csv(result$gene_qc, file.path(OUTPUT_DIR, "vector_QC_by_gene.csv"), row.names = FALSE)

  # Transcript end profiles; each numbered sample is visible in the legend.
  last_codons <- all_data %>%
    group_by(Gene) %>%
    filter(codon_pos > max(codon_pos) - ra$end_window_codons) %>%
    ungroup()
  for (gene in unique(last_codons$Gene)) {
    df <- last_codons %>% filter(Gene == gene)
    unique_codons <- df %>% distinct(codon_pos, codon_rna, is_stop) %>% arrange(codon_pos)
    tick_labels <- paste0(unique_codons$codon_pos, "_",
                         ifelse(unique_codons$is_stop, paste0("STOP_", unique_codons$codon_rna),
                                unique_codons$codon_rna))
    p <- ggplot(df, aes(codon_pos, Percent, color = Sample, group = sample_id)) +
      geom_line(linewidth = 0.7, na.rm = TRUE) + facet_wrap(~Condition, nrow = 1) +
      scale_x_continuous(breaks = unique_codons$codon_pos, labels = tick_labels) +
      theme_bw() + theme(axis.text.x = element_text(angle = 90, size = 6, hjust = 1)) +
      labs(title = paste(cfg$cell_line, "last", ra$end_window_codons, "codons:", gene),
           x = "Codon in supplied analysis sequence", y = "% occupancy", color = "Sample")
    ggsave(file.path(OUTPUT_DIR, paste0("last_codons_", gene, ".png")), p,
           width = 16, height = 5, dpi = 300)
  }

  # The first candidate is selected within each sample and gene.
  # This is not a validated termination-site annotation.
  stop_data <- all_data %>% filter(is_stop) %>%
    group_by(Gene, Condition, sample_id, Sample, Replicate) %>%
    slice_min(codon_pos, n = 1, with_ties = FALSE) %>% ungroup()
  write.csv(stop_data, file.path(OUTPUT_DIR, "candidate_stop_occupancy.csv"), row.names = FALSE)
  stop_mean <- stop_data %>% group_by(Gene, Condition) %>%
    summarise(mean_percent = if (all(is.na(Percent))) NA_real_ else mean(Percent, na.rm = TRUE),
              samples_with_counts = sum(!is.na(Percent)), .groups = "drop")
  write.csv(stop_mean, file.path(OUTPUT_DIR, "candidate_stop_condition_means.csv"), row.names = FALSE)
  if (nrow(stop_data)) {
    p <- ggplot() +
      geom_point(data = stop_mean, aes(Gene, mean_percent, color = Condition),
                 size = 4, shape = 18, position = position_dodge(width = 0.7), na.rm = TRUE) +
      geom_point(data = stop_data, aes(Gene, Percent, color = Condition, shape = Replicate),
                 size = 2.5, position = position_jitterdodge(dodge.width = 0.7, jitter.width = 0.15),
                 na.rm = TRUE) +
      theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      labs(title = paste(cfg$cell_line, "occupancy at first candidate stop"),
           x = "Transcript", y = "% occupancy", shape = "Replicate number", color = "Condition")
    ggsave(file.path(OUTPUT_DIR, "candidate_stop_occupancy.png"), p, width = 13, height = 5, dpi = 300)
  }

  # Descriptive vector quality; thresholds are inherited source heuristics.
  qc_summary <- result$gene_qc %>% group_by(Sample, sample_id, Condition) %>%
    summarise(TotalCounts = sum(TotalReads), VectorsLoaded = n_distinct(Gene),
              TranscriptsWithCounts = sum(TotalReads > 0),
              Frame0 = mean(Frame0, na.rm = TRUE), Frame1 = mean(Frame1, na.rm = TRUE),
              Frame2 = mean(Frame2, na.rm = TRUE), FrameBias = mean(FrameBias, na.rm = TRUE),
              Fourier3nt = mean(Fourier3nt, na.rm = TRUE), .groups = "drop") %>%
    mutate(QC_Flag = case_when(TotalCounts == 0 ~ "NoCounts",
                              Frame0 >= 70 & FrameBias >= 3 ~ "Excellent",
                              Frame0 >= 50 & FrameBias >= 2 ~ "Good", TRUE ~ "Poor"))
  write.csv(qc_summary, file.path(OUTPUT_DIR, "vector_QC_summary.csv"), row.names = FALSE)
  print(qc_summary)

  # Read length QC is optional when working only from existing count vectors.
  phasing_data <- read_phasing(inputs)
  if (!is.null(phasing_data)) {
    write.csv(phasing_data, file.path(OUTPUT_DIR, "phasing_by_read_length.csv"), row.names = FALSE)
    qc_phasing <- phasing_data %>% group_by(Sample, sample_id, Condition) %>%
      summarise(TotalReads = sum(reads_counted),
                Phase0 = if (sum(reads_counted) > 0) weighted.mean(phase0, reads_counted) else NA_real_,
                Phase1 = if (sum(reads_counted) > 0) weighted.mean(phase1, reads_counted) else NA_real_,
                Phase2 = if (sum(reads_counted) > 0) weighted.mean(phase2, reads_counted) else NA_real_,
                .groups = "drop")
    write.csv(qc_phasing, file.path(OUTPUT_DIR, "phasing_QC_summary.csv"), row.names = FALSE)
    p <- ggplot(phasing_data, aes(read_length, phase0, color = Condition, group = sample_id)) +
      geom_line(alpha = 0.5) + geom_point(size = 2) + theme_bw() +
      labs(title = paste(cfg$cell_line, "phase 0 by read length"), x = "Read length (nt)", y = "Phase 0 fraction")
    ggsave(file.path(OUTPUT_DIR, "phase0_by_read_length.png"), p, width = 10, height = 5, dpi = 300)
  }
  writeLines(capture.output(sessionInfo()), file.path(OUTPUT_DIR, "sessionInfo.txt"))
  writeLines(jsonlite::toJSON(cfg, pretty = TRUE, auto_unbox = TRUE), file.path(OUTPUT_DIR, "analysis_settings.json"))
  message("Tables and plots saved to: ", OUTPUT_DIR)
  invisible(OUTPUT_DIR)
}
