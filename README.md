# MitoRiboSeq custom scripts

Reusable scripts for mitochondrial ribosome profiling: adapter trimming, alignment, offset review, count-vector export and R analysis. HEK and N2a profiles provide examples that readers can adapt to their own samples, references and protocols.

Repository: [dariaki-lab/MitoRiboSeq-custom-scripts-](https://github.com/dariaki-lab/MitoRiboSeq-custom-scripts-)

## Start here

Keep the scripts and supporting files in the same directory. Input data may be elsewhere if the paths are set correctly. Run commands from the script directory.

| File | Purpose |
| --- | --- |
| `01_process_reads.py` | Process reads through named stages |
| `02_analyze_vectors.R`, `R_functions.R` | Analyze vectors and export tables and plots |
| `sample_info.csv` | Sample names, conditions, replicate metadata and input paths |
| `HEK_settings.json`, `N2a_settings.json` | Editable processing and analysis profiles |
| `HEK_transcripts.tsv`, `N2a_transcripts.tsv` | Gene-to-vector and sequence-coordinate mapping |
| [WORKFLOW.md](WORKFLOW.md) | Setup, commands and adaptation instructions |

## Workflow

| Step | Script stage |
| --- | --- |
| Build a reference index; trim and align reads | Python: `index`, `trim`, `align` |
| Estimate candidate offsets and review them | Python: `psite` |
| Check phasing; export counts and vectors | Python: `phase`, `counts`, `vectors` |
| Calculate occupancy and create tables and plots | R: `02_analyze_vectors.R` |

1. Fill the sample table with your data paths and actual replicate metadata.
2. Review the selected JSON profile: adapter, read lengths, reference, annotations, calibration regions and output paths.
3. Supply reviewed offsets and record the mapped ribosome site. For standard plastid P-site calibration at start codons, use `offset_roi_key: "roi_start"`; the supplied profiles currently select `roi_stop` and require review. See [WORKFLOW.md](WORKFLOW.md#3-review-processing-settings).
4. Before R analysis, verify the transcript TSV, analysis sequences and vector coordinates. Set the confirmation flags only after these checks.

The supplied settings are examples. Raw reads or count vectors, the matching alignment reference and final offsets must be supplied for your analysis. The current alignment wrapper uses single-end reads.

Preview a stage without running sequencing tools:

```bash
python3 01_process_reads.py --analysis HEK --stage align --dry-run
```

After preparing the vectors and confirming the R inputs:

```bash
Rscript install_R_packages.R
Rscript 02_analyze_vectors.R HEK
```

Use `N2a` for the mouse example. Outputs go to the configured `output_dir`; R tables and plots are in its `analysis/` subdirectory. R can also analyze existing compatible vectors. See [WORKFLOW.md](WORKFLOW.md) for installation and all processing commands.

## Sample-name templates

Use these display labels for each cell line. They are editable template slots.

| Cell line | Condition | Replicate 1 | Replicate 2 | Replicate 3 |
| --- | --- | --- | --- | --- |
| HEK | WT | WT1 | WT2 | WT3 |
| HEK | MTRF1 KO | MTRF1 KO1 | MTRF1 KO2 | MTRF1 KO3 |
| HEK | MTRF1A KO | MTRF1A KO1 | MTRF1A KO2 | MTRF1A KO3 |
| HEK | MTRF1/MTRF1A dKO | MTRF1/MTRF1A dKO1 | MTRF1/MTRF1A dKO2 | MTRF1/MTRF1A dKO3 |
| N2a | WT | WT1 | WT2 | WT3 |
| N2a | MTRF1 KO | MTRF1 KO1 | MTRF1 KO2 | MTRF1 KO3 |
| N2a | MTRF1A KO | MTRF1A KO1 | MTRF1A KO2 | MTRF1A KO3 |
| N2a | MTRF1/MTRF1A dKO | MTRF1/MTRF1A dKO1 | MTRF1/MTRF1A dKO2 | MTRF1/MTRF1A dKO3 |

File IDs include the cell line and use underscores, for example `HEK_WT1`, `HEK_MTRF1A_KO2` and `N2a_MTRF1_MTRF1A_dKO3`. The scripts support other conditions and replicate counts. Numbering alone does not establish replicate independence.

For a WT/KO comparison, use a working sample CSV containing the intended samples and point `sample_sheet` to it. Python `--samples` limits processing only; R analyzes every row for the configured cell line in that CSV.

## Optional Word guides

These concise guides can be uploaded individually alongside the scripts:

- [HEK_processing_workflow.docx](HEK_processing_workflow.docx)
- [HEK_MTRF1_KO_processing_notes.docx](HEK_MTRF1_KO_processing_notes.docx)
- [HEK_MTRF1A_KO_processing_notes.docx](HEK_MTRF1A_KO_processing_notes.docx)
- [HEK_MTRF1_MTRF1A_dKO_processing_notes.docx](HEK_MTRF1_MTRF1A_dKO_processing_notes.docx)
- [N2a_processing_workflow.docx](N2a_processing_workflow.docx)

## Credit and citation

Original script/workflow credit: **Annika Krüger**, recorded as creator or editor of the supplied source documents. When reporting an analysis, cite the associated manuscript when available and record the repository release or commit used. Complete the author-approved citation and license information before a manuscript release.



# Reusing the MitoRiboSeq workflow

This guide describes the shared workflow. HEK and N2a are example profiles; sample names, conditions, references and scientific parameters are editable. Processing depends on inputs and settings appropriate to your experiment.

## 1. Install software and keep the files together

Run commands from the directory containing `01_process_reads.py`, `02_analyze_vectors.R`, `R_functions.R` and the supporting configuration files.

Requirements: Python 3.9 or later, Cutadapt, Bowtie2, SAMtools and plastid for processing; R with `jsonlite`, `dplyr` and `ggplot2` for vector analysis. The optional genomic-CDS extraction route also needs `Biostrings` and `rtracklayer`.

```bash
conda env create -f preprocess_environment.yml
conda env create -f plastid_environment.yml
Rscript install_R_packages.R
```

These environment files are installation examples, rather than frozen software-version records. Record the versions used for your analysis. The supplied Bowtie2 wrapper uses single-end input and local alignment; other read layouts or alignment strategies require adapting the alignment commands.

## 2. Prepare a sample table and reference inputs

Edit a working copy of `sample_info.csv`. Use one row per sequencing library.

| Column | Meaning |
| --- | --- |
| `sample_id` | Unique file-safe ID; use letters, numbers, underscores, dots or hyphens |
| `cell_line` | Group selected by the settings profile, such as `HEK` or `N2a` |
| `condition` | Experimental group; must appear in `r_analysis.condition_order` |
| `replicate` | Actual replicate identifier or number |
| `sample_label` | Unique display label for plots |
| `fastq_path` | Path to the input FASTQ; needed for trimming |
| `vector_folder` | Sample folder below the vector root; an empty value falls back to `sample_id` in R |
| `offset_file` | Reviewed sample-specific offset table; an empty value uses the JSON `offset_file` |

The [README.md](README.md#sample-name-templates) lists the WT1–3, MTRF1 KO1–3, MTRF1A KO1–3 and MTRF1/MTRF1A dKO1–3 templates for both cell lines. For example, `HEK_MTRF1_KO2` is the file ID for display label `MTRF1 KO2`. Templates can be expanded or replaced. Record the actual experimental replicate structure.

For one WT/KO comparison, use a CSV containing only the intended samples, point `sample_sheet` to it and use a separate output directory. Set `r_analysis.condition_order` to the conditions present. Python `--samples` selects processing IDs only; it does not restrict the R analysis, which reads all rows for the configured `cell_line`.

Supply the following inputs:

| Input | Requirement |
| --- | --- |
| Alignment FASTA | The actual reference to which reads will be aligned |
| GTF and BED | Features used for transcript counts and nucleotide vectors |
| Start-region ROI file | plastid-compatible start-codon windows for the standard calibration and phasing route |
| Final offset tables | Reviewed read-length-to-position offsets |
| Transcript TSV and analysis FASTA | Gene-to-vector filenames, sequence identifiers and coordinate origins for R |

Reference sequence names, annotations and ROI coordinates must agree. Use inputs for the correct species and reference build or custom transcript construction. Per-gene analysis FASTAs are R sequence lookups; they do not replace the alignment FASTA.

## 3. Review processing settings

Select `HEK_settings.json` or `N2a_settings.json` and edit the following fields. Existing adapter sequences, read-length ranges and reference choices are source-workflow examples, not universal recommendations.

| Field | What to review |
| --- | --- |
| `cell_line`, `sample_sheet` | The samples included in this run |
| `adapter`, `max_n`, stage read-length ranges | Library preparation and read-quality choices |
| `reference_fasta`, `bowtie2_index` | Alignment reference and index location |
| `annotation_gtf`, `annotation_bed`, ROI paths | Features and coordinates matching that reference |
| `offset_roi_key` | Calibration landmark; use `roi_start` for standard start-codon P-site calibration |
| `offset_file`, `mapping_site` | Reviewed offsets and the ribosome position they represent |
| `output_dir`, `trimmed_dir` | Locations for this run; keep distinct analyses separate |
| `threads`, `codon_buffer` | Computing resources and the phasing window |
| `settings_confirmed` | Set to `true` after reviewing processing settings |

**Calibration needs explicit review.** The supplied profiles currently use `offset_roi_key: "roi_stop"`. Standard plastid `psite` estimates P-site offsets from start-codon windows. For that mode, change the key to `roi_start`, supply the appropriate ROI file and verify that the method fits your footprint protocol. Inspect calibration plots rather than automatically accepting the estimated offsets. See the [official plastid psite documentation](https://plastid.readthedocs.io/en/latest/generated/plastid.bin.psite.html).

Select the profile and preview commands:

```bash
PROFILE=HEK
# Set PROFILE=N2a for the mouse example.
python3 01_process_reads.py --analysis "$PROFILE" --stage align --dry-run
```

`--dry-run` prints commands without executing tools or creating results. It is a command preview, not a check that the supplied data or references are valid.

## 4. Process reads

After confirming the settings and supplying the inputs, build the index, trim and align:

```bash
conda activate mitoribo-preprocessing
python3 01_process_reads.py --analysis "$PROFILE" --stage index
python3 01_process_reads.py --analysis "$PROFILE" --stage trim
python3 01_process_reads.py --analysis "$PROFILE" --stage align
```

Review trimming and alignment logs. Alignment produces mapped, sorted and indexed BAM files. Then estimate candidate offsets:

```bash
conda activate mitoribo-plastid
python3 01_process_reads.py --analysis "$PROFILE" --stage psite
```

Review the plots and select the retained read lengths. Supply final offset tables with an offset for every length used downstream. Set each sample's `offset_file`, or the shared JSON path; the JSON path supports `{sample_id}` for separate tables. Record `mapping_site` as `P`, `A` or `custom` according to the mapped position. This field does not convert the offsets.

Check phasing, then generate counts and vectors:

```bash
python3 01_process_reads.py --analysis "$PROFILE" --stage phase
python3 01_process_reads.py --analysis "$PROFILE" --stage counts
python3 01_process_reads.py --analysis "$PROFILE" --stage vectors
```

Inspect phasing by read length and resolve poor or ambiguous calibration before interpretation. The wrapper requires reviewed offsets for `phase`, `counts` and `vectors`. Use `--force` only when intentionally rerunning a stage and replacing its outputs.

## 5. Confirm coordinates and run R

Review the transcript TSV selected by `r_analysis.vector_mapping`:

| Column | Meaning |
| --- | --- |
| `Gene` | Gene label in output tables |
| `vector_file` | Filename within each sample's vector folder |
| `sequence_gene` | Matching sequence identifier in the analysis FASTA |
| `vector_start_1based` | Vector position of the first base of that analysis sequence, using 1-based indexing |

Confirm that vectors and sequences have the same orientation and coordinate origin, that the selected sequence starts in the intended reading frame, and that its complete length fits inside the vector. Set `r_analysis.coordinates_confirmed` to `true` after review. Check `r_analysis.condition_order` and any candidate-codon or sequence-window choices for the biological question being asked.

```bash
Rscript 02_analyze_vectors.R "$PROFILE"
```

For **existing vectors**, skip read processing and use this R step once settings and coordinates are confirmed. Set `r_analysis.vector_dir` to the existing vector root and `vector_folder` to each sample's subfolder. Files should contain one numeric column of nonnegative nucleotide counts, with names matching `vector_file` in the TSV.

The built-in R entry point also accepts `HEK_unprocessed` when its matching inputs are prepared. That is an optional reference choice; the required input-review steps remain the same.

## 6. Outputs and interpretation

Processing creates `trimmed/`, `bam/`, `offset_estimates/`, `phasing/`, `counts/`, `vectors/` and `logs/` below the configured output locations. R writes an `analysis/` directory containing sample metadata, nucleotide counts, codon occupancy, coordinate/QC tables and PNG plots. With the supplied profiles, output roots are `results/HEK/` and `results/N2a/`.

Codon occupancy is calculated as the sum of three adjacent nucleotide counts divided by total counts in the selected sequence region, multiplied by 100. Confirm the normalization region and the mapped ribosome site before interpretation. Candidate-codon summaries are exploratory; condition means are descriptive and the scripts do not perform differential hypothesis tests.

Keep the final sample table, settings, reference versions, offset tables and software versions with the analysis. Record the repository release or commit used.

## 7. Adapt to another cell line or reference

Copy an existing JSON profile to `custom_settings.json` and update its `analysis_id`, `cell_line`, sample table, references, mappings, output paths and scientific parameters.

For Python, use the custom-config option in place of `--analysis`:

```bash
python3 01_process_reads.py --config custom_settings.json --stage align --dry-run
# After review, use --config custom_settings.json for each processing stage.
```

The numbered R entry point accepts the three built-in profile names. For a custom profile, save the following as `custom_analysis.R` in the script directory and run `Rscript custom_analysis.R` after confirming inputs:

```r
library(dplyr)
library(ggplot2)
Sys.setenv(MITORIBO_ROOT = normalizePath(getwd(), mustWork = TRUE))
source("R_functions.R")
inputs <- load_analysis_inputs("custom_settings.json")
run_vector_analysis(inputs)
```

This uses the same shared functions as `02_analyze_vectors.R`. Protocols requiring different trimming, paired-end alignment, reference construction or coordinate handling need corresponding changes to the processing logic.

## Credit

Original script/workflow credit: **Annika Krüger**, recorded as creator or editor of the supplied source documents.
