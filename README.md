[README (3).md](https://github.com/user-attachments/files/33103515/README.3.md)
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
