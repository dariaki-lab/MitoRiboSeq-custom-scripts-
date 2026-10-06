# Reusable R Markdown analysis templates

These notebooks replace study-specific folder lists, clone-name parsing, embedded reference sequences and fixed plot descriptions with editable metadata and configuration files.

| Notebook | Default profile |
| --- | --- |
| [HEK_analysis_template.Rmd](HEK_analysis_template.Rmd) | `HEK_settings.json` |
| [N2a_analysis_template.Rmd](N2a_analysis_template.Rmd) | `N2a_settings.json` |
| [HEK_unprocessed_analysis_template.Rmd](HEK_unprocessed_analysis_template.Rmd) | `HEK_unprocessed_settings.json` |

## Files to keep together

Upload the notebooks and `R_functions.R` individually alongside the sample table, settings JSON files, transcript TSVs and reference inputs already provided with the reusable script package. `R_functions.R` is supplied here as the same shared helper used by `02_analyze_vectors.R`; its calculations have not been changed for these notebooks.

The notebooks are R Markdown entry points for that workflow. They require the helper and configuration inputs; they are not independent replacements for those files. They start from count vectors, so raw-read processing is separate.

## Install once

Run in the R console:

```r
install.packages(c("rmarkdown", "knitr", "jsonlite", "dplyr", "ggplot2"))
```

For genomic-CDS extraction, also install:

```r
install.packages("BiocManager")
BiocManager::install(c("Biostrings", "rtracklayer"))
```

HTML rendering also requires Pandoc. Use RStudio's Knit workflow or an R installation with Pandoc available.

## Edit three notebook parameters

| YAML parameter | Meaning |
| --- | --- |
| `project_dir` | Directory containing the helper and configuration files; default `.` means the notebook directory when knitting normally |
| `settings_file` | Selected JSON profile, or a custom settings file |
| `run_analysis` | `false` shows the guide and code; `true` runs the configured analysis |

Start by opening a notebook in RStudio and clicking **Knit**. The default guide mode requires no experimental data and does not calculate results. To run the analysis, review the inputs below and set `run_analysis: true`, or use:

```r
rmarkdown::render("HEK_analysis_template.Rmd", params = list(run_analysis = TRUE))
rmarkdown::render("N2a_analysis_template.Rmd", params = list(run_analysis = TRUE))
rmarkdown::render("HEK_unprocessed_analysis_template.Rmd", params = list(run_analysis = TRUE))
```

Choose the command for the prepared profile. Each HTML report includes the workflow, visible code, selected sample metadata, table previews and sequence-end plots.

## Review the inputs

1. Fill `sample_info.csv` with actual sample IDs, conditions, replicate metadata and vector-folder names. The WT1–3, MTRF1 KO1–3, MTRF1A KO1–3 and MTRF1/MTRF1A dKO1–3 labels remain as editable examples for each cell line.
2. Set the JSON `cell_line`, `sample_sheet`, `output_dir` and `r_analysis.condition_order`. For a comparison, use a working CSV containing only the intended samples. R reads every row for the configured cell line; Python `--samples` does not restrict it.
3. Supply the vector root through `r_analysis.vector_dir`, or use the default `<output_dir>/vectors`. Each sample's `vector_folder` contains the files named in the transcript TSV.
4. Review gene-to-vector and sequence identifiers, orientation, reading frame, sequence extent and `vector_start_1based`. The latter is the 1-based vector position corresponding to the analysis sequence's first base.
5. Provide matching analysis sequences. A nonempty `sequence_fasta` uses a FASTA lookup; an empty value uses `genome_fasta`, `cds_gff3` and `gene_column` to extract CDS sequences. The supplied extraction helper supports one CDS segment per gene with phase 0 or unspecified phase; more complex annotations require an explicit extraction rule.
6. Confirm `mapping_site` as `P`, `A` or `custom`, and set `settings_confirmed` and `r_analysis.coordinates_confirmed` to `true` after review. These values document reviewed inputs; they do not convert count-vector positions.

## Use another cell line or reference

Copy a settings profile, update its inputs and use it with any notebook:

```r
rmarkdown::render(
  "HEK_analysis_template.Rmd",
  params = list(
    project_dir = ".",
    settings_file = "custom_settings.json",
    run_analysis = TRUE
  ),
  output_file = "custom_analysis_report.html"
)
```

Change the example subtitle if you want the HTML heading to describe your custom analysis. Set a distinct `output_dir` for each reference or sample selection. Reusing an output location replaces matching outputs.

## What is retained and what to interpret carefully

The original occupancy formula is retained: three nucleotide counts per codon divided by total counts in the selected sequence region, multiplied by 100. All notebooks use the shared functions rather than three copies of the calculation and plotting code.

The shared workflow uses confirmed vector coordinates and stops for missing or out-of-range inputs instead of silently skipping libraries or padding counts. Partial trailing codons and zero-count regions are explicitly recorded. These templates are an adaptation for reuse, not a claim to reproduce every historical plot or output exactly.

Candidate-codon flags and `STOP_` labels indicate configured sequence matches. They do not establish termination sites or stalling. Condition means and QC summaries are descriptive; source-derived QC flags are not universal acceptance thresholds.

## Credit and documentation

Original script/workflow credit: **Annika Krüger**, recorded as creator or editor of the supplied workflow documents. Cite the associated manuscript when available and record the repository release or commit used.

R Markdown parameter usage follows the [R Markdown Cookbook](https://bookdown.org/yihui/rmarkdown-cookbook/parameterized-reports.html). Figure inclusion uses [knitr's include_graphics](https://search.r-project.org/CRAN/refmans/knitr/html/include_graphics.html).
