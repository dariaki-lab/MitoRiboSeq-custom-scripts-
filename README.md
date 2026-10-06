[README.md](https://github.com/user-attachments/files/33101193/README.md)
# MitoRiboSeq custom scripts

Custom scripts and processing notes for mitochondrial ribosome profiling in **HEK (human)** and **N2a (mouse)** cells, comparing WT, MTRF1 KO, MTRF1A KO and MTRF1/MTRF1A double KO samples.

The repository provides a general workflow from sequencing reads to mitochondrial ribosome occupancy tables and plots. The scripts can be adapted by editing the sample table, reference paths and analysis settings.

Repository: [dariaki-lab/MitoRiboSeq-custom-scripts-](https://github.com/dariaki-lab/MitoRiboSeq-custom-scripts-)

## 1. Workflow overview

| Step | What the step does | Tool or script |
| --- | --- | --- |
| 1 | Remove adapters and filter read lengths | Cutadapt, through `01_process_reads.py` |
| 2 | Align reads to the selected mitochondrial reference | Bowtie2, through `01_process_reads.py` |
| 3 | Keep mapped reads, sort BAM files and create BAM indexes | SAMtools, through `01_process_reads.py` |
| 4 | Estimate read-length offsets and inspect reading-frame phasing | plastid: `psite` and `phase_by_size` |
| 5 | Generate transcript counts and nucleotide count vectors | plastid: `counts_in_region` and `get_count_vectors` |
| 6 | Calculate codon occupancy and export tables and plots | `02_analyze_vectors.R` |

The main workflow uses references representing processed mitochondrial transcripts with poly(A) extensions. An optional HEK workflow analyzes an unprocessed reference. Reference sequences, annotations and offsets must correspond to the selected workflow.

## 2. Sample and replicate names

Use the following labels for each cell line. Map each label to its actual sequencing library in [sample_info.csv](sample_info.csv); the template numbering must match the final experimental sample metadata.

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

For filenames and script input, add the cell line and replace spaces and `/` with underscores. Examples: `HEK_WT1`, `HEK_MTRF1A_KO2` and `N2a_MTRF1_MTRF1A_dKO3`. Use each sample ID consistently in the sample table and generated results; FASTQ paths may retain the original library filenames.

Original clone and batch names are retained in [original_sample_names.csv](original_sample_names.csv). Assigning a numbered label does not establish that a library is an independent biological replicate.

## 3. Main repository files

Keep the supplied source files together in one directory. Input data can remain elsewhere if the paths are entered in the sample table and settings files.

| File | Purpose |
| --- | --- |
| `01_process_reads.py` | Run processing stages from FASTQ to count vectors |
| `02_analyze_vectors.R` | Analyze vectors and produce CSV tables and PNG plots |
| `R_functions.R` | Shared functions used by the R analysis |
| `sample_info.csv` | Sample IDs, conditions, replicate numbers and input paths |
| `HEK_settings.json`, `N2a_settings.json` | Cell-line-specific processing and analysis settings |
| `HEK_transcripts.tsv`, `N2a_transcripts.tsv` | Gene-to-vector mapping and sequence coordinates |
| `HEK_*`, `N2a_*` reference files | Annotations, regions of interest and analysis sequences |
| `install_R_packages.R` | Install the required R packages |
| `preprocess_environment.yml`, `plastid_environment.yml` | Example Conda environments |
| [WORKFLOW.md](WORKFLOW.md) | Additional instructions and scientific notes |
| [original_workflows.zip](original_workflows.zip) | Historical notebooks and processing notes |

`HEK_unprocessed_settings.json` and `HEK_unprocessed_transcripts.tsv` support the optional unprocessed comparison. `summarize_processing_logs.py` summarizes trimming and alignment logs; the `check_*` files provide optional software checks.

## 4. Requirements and input data

- Python 3.9 or later, Cutadapt, Bowtie2, SAMtools and plastid for read processing.
- R with `jsonlite`, `dplyr` and `ggplot2` for processed-vector analysis.
- FASTQ files, or existing count vectors for R-only analysis.
- The exact alignment FASTA, matching annotations and region-of-interest files, and reviewed offset tables.

The optional unprocessed R analysis also requires `Biostrings` and `rtracklayer`.

Supply raw reads or count vectors, and the exact custom alignment FASTAs, separately. N2a processing also needs the final mouse region-of-interest files. The included per-gene analysis FASTAs are sequence lookups for R; supply the corresponding custom alignment references separately.

## 5. Set up your analysis

1. Edit `sample_info.csv` with your FASTQ paths or existing vector-folder names. Raw sequencing files do not need to be renamed.
2. Select `HEK_settings.json` or `N2a_settings.json`. Review the adapter, read-length ranges, reference paths and output paths, then set `settings_confirmed` to `true`.
3. Estimate and review offsets. Supply the final offset-table paths and record `mapping_site` as `P`, `A` or `custom`, according to the position represented by your count vectors.
4. Before running R, confirm gene-to-vector mapping, strand, sequence origin and `vector_start_1based` in the matching transcript table. Set `r_analysis.coordinates_confirmed` to `true` after this review.

The public settings use the adapter `TGGAATTCTCGGGTGCCAAGG`, trimming lengths of 25–45 nt and count/vector lengths of 30–45 nt. Historical notes include alternative filters; select settings that match your experiment and manuscript Methods.

Preview a stage without executing sequencing tools or writing outputs:

```bash
python3 01_process_reads.py --analysis HEK --stage align --dry-run
python3 01_process_reads.py --analysis N2a --stage vectors --dry-run
```

## 6. Run the workflow

Run commands from the directory containing the scripts. The environment files are editable examples; record the versions used for your analysis.

Create the processing environments once:

```bash
conda env create -f preprocess_environment.yml
conda env create -f plastid_environment.yml
```

Build the reference index, trim reads and align:

```bash
conda activate mitoribo-preprocessing
python3 01_process_reads.py --analysis HEK --stage index
python3 01_process_reads.py --analysis HEK --stage trim
python3 01_process_reads.py --analysis HEK --stage align
```

Estimate candidate offsets:

```bash
conda activate mitoribo-plastid
python3 01_process_reads.py --analysis HEK --stage psite
```

**Review the offset plots and supply final offset tables before continuing.** The source notes use a COX1 stop landmark for calibration. Confirm the landmark and mapped ribosome site for your analysis; the wrapper does not automatically convert P-site offsets to A-site offsets.

Then check phasing and export counts and vectors:

```bash
python3 01_process_reads.py --analysis HEK --stage phase
python3 01_process_reads.py --analysis HEK --stage counts
python3 01_process_reads.py --analysis HEK --stage vectors
```

For mouse samples, replace `HEK` with `N2a` and use the mouse inputs. Add `--samples HEK_WT1` or another sample ID to process selected samples.

Analyze the completed vectors:

```bash
Rscript install_R_packages.R
Rscript 02_analyze_vectors.R HEK
```

After preparing N2a vectors, run `Rscript 02_analyze_vectors.R N2a`. For R-only use, the default vector locations are `results/HEK/vectors/<sample_id>/` and `results/N2a/vectors/<sample_id>/`. Each sample folder contains the vector files specified in its transcript table. [WORKFLOW.md](WORKFLOW.md) explains how to use existing folders.

## 7. Outputs and interpretation

Results are written under `results/HEK/` or `results/N2a/`. The `analysis/` directory contains codon occupancy, nucleotide counts, candidate-stop summaries, coordinate/QC tables and PNG plots. Processing logs are saved in the corresponding `logs/` directory.

Codon occupancy is calculated as the sum of three adjacent nucleotide counts divided by total counts in the selected sequence region, multiplied by 100. Confirm the normalization region and vector coordinates before interpreting results. Candidate-stop summaries and condition means are descriptive outputs.

The Python wrapper passed 14 tests during package preparation. The sequencing tools, R analysis and manuscript results still require validation on the study inputs.
