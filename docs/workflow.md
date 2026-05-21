# Workflow

This repository presents the small non-coding RNA pipeline as sectioned code blocks in the main `README.md`, with matching R templates and one Bash template in `scripts/`.

Use the README sections in this order:

1. Conventional Mapping and DE
2. SPORTS Bash workflow, then SPORTS Summary to Count Matrix and R DE
3. tDRMapper Bash workflow, then tDRMapper R DE
4. proTRAC Bash cluster preparation, then proTRAC Cluster R DE
5. General Count Matrix DE
6. proTRAC Non-DE Unique Cluster Workflow
7. Volcano Plot
8. MA Plot
9. Heatmap

## Script Templates

The same workflows are available as editable script templates:

```text
scripts/
|-- conventional_de.R
|-- sports_run_workflow.sh
|-- sports_summary_to_counts.R
|-- tdrmapper_run_workflow.sh
|-- tdrmapper_de.R
|-- protrac_run_prepare_clusters.sh
|-- protrac_cluster_de.R
|-- general_count_matrix_de.R
|-- protrac_non_de_unique_clusters.sh
`-- de_plots.R
```

## Tool-Specific Order

SPORTS:

```text
scripts/sports_run_workflow.sh
scripts/sports_summary_to_counts.R
scripts/general_count_matrix_de.R
```

tDRMapper:

```text
scripts/tdrmapper_run_workflow.sh
scripts/tdrmapper_de.R
```

proTRAC cluster DE:

```text
scripts/protrac_run_prepare_clusters.sh
scripts/protrac_cluster_de.R
```

proTRAC unique/no-overlap clusters:

```text
scripts/protrac_non_de_unique_clusters.sh
```

This unique-cluster workflow is not a DESeq2 step.

## How to Use the Code Sections

Copy the section you need into R, RStudio, or Bash, then replace placeholders such as:

```text
your_working_directory
your_output_directory
yourfile.fastq.gz
your_reference.fa
your_count_matrix.txt
your_metadata.csv
your_TREATMENT
CONTROL
```

Keep sample names consistent between count matrices and metadata files. DESeq2 requires count matrix columns and metadata row names to match exactly.

Section 6 is intentionally not a DESeq2 step. It follows the proTRAC workflow for merging treatment and control GTFs separately, finding no-overlap unique clusters with BEDtools, and intersecting those unique clusters with genes for functional study.

## Suggested Output Organization

```text
results/
|-- conventional/
|-- sports/
|   |-- miRNA/
|   |-- tsRNA/
|   `-- rsRNA/
|-- tDRMapper/
`-- proTRAC/
```

Generated outputs are ignored by Git because they can be large and project-specific.
