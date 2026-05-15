# Workflow

This repository presents the small non-coding RNA pipeline as sectioned code blocks in the main `README.md` and as matching R script templates in `scripts/`.

Use the README sections in this order:

1. Conventional Mapping and DE
2. SPORTS Summary to Count Matrix
3. tDRMapper DE
4. proTRAC Cluster DE
5. General Count Matrix DE
6. Non-overlapping Cluster Comparison
7. Volcano Plot
8. MA Plot
9. Heatmap

## Script Templates

The same workflows are available as editable R scripts:

```text
scripts/
|-- conventional_de.R
|-- sports_summary_to_counts.R
|-- tdrmapper_de.R
|-- protrac_cluster_de.R
|-- general_count_matrix_de.R
|-- non_overlapping_clusters.R
`-- de_plots.R
```

## How to Use the Code Sections

Copy the section you need into R or RStudio, then replace placeholders such as:

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
