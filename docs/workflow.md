# Workflow

This pipeline organizes differential expression analysis for small non-coding RNA datasets. The current scripts are templates, so each project should update paths, sample labels, treatment names, and filtering thresholds before execution.

## 1. Prepare Project Metadata

Create a metadata table with one row per sample. At minimum, keep these fields consistent across scripts:

- `sample`
- `condition`
- `file`

Use `CONTROL` as the reference condition unless your study design requires a different baseline.

## 2. Choose an Analysis Module

### Conventional Reference-based Workflow

Use [scripts/Conventional_DE.R](../scripts/Conventional_DE.R) when you have a custom reference FASTA containing small RNA transcripts such as miRNAs, tsRNAs, piRNAs, or rsRNAs.

Main steps:

1. Build an Rsubread index from the reference FASTA.
2. Align FASTQ files with `subjunc`.
3. Build a feature annotation table from the FASTA.
4. Count reads with `featureCounts`.
5. Run DESeq2.
6. Export tables and plots.

### SPORTS Workflow

Use [scripts/SPORTS_DE.R](../scripts/SPORTS_DE.R) when SPORTS has already produced per-sample summary files.

Main steps:

1. Read all `*_summary.txt` files.
2. Classify entries into miRNA, tsRNA, rsRNA, or other.
3. Build count matrices for each group.
4. Run DESeq2 separately for miRNA, tsRNA, and rsRNA.
5. Export tables and plots.

### tDRMapper Workflow

Use [scripts/tDRMapper_DE.R](../scripts/tDRMapper_DE.R) when tDRMapper has produced a tRNA-derived RNA count matrix.

Main steps:

1. Read `tRNA_count_matrix.txt`.
2. Define control and treatment metadata.
3. Filter low-count tRNA clusters.
4. Run DESeq2.
5. Export significant and total tRNA expression tables.
6. Generate volcano plot and heatmap.

### proTRAC Workflow

Use [scripts/proTRAC_DE.R](../scripts/proTRAC_DE.R) when proTRAC has produced piRNA cluster count matrices.

Main steps:

1. Read `piRNA_cluster_count_matrix_renamed.txt`.
2. Define control and treatment metadata.
3. Filter low-count piRNA clusters.
4. Run DESeq2.
5. Export significant and total piRNA cluster expression tables.
6. Generate MA plot, volcano plot, and heatmap.

## 3. Update Placeholders

Before running a script, search for these placeholders and replace them with study-specific values:

```text
your_directory_here
your_output_path
your_summary_directory_path
your_TREATMENT
your_treatment
CONTROL
your_fasta.fa
your_treatment_file.fastqsanger.gz
your_control_file.fastqsanger.gz
```

## 4. Check Sample Ordering

DESeq2 requires the count matrix columns and metadata row names to be in the same order. Each script includes a `stopifnot()` check or matching logic. If this check fails, fix the sample names before continuing.

## 5. Review Thresholds

The scripts currently use common defaults:

- adjusted p-value cutoff: `0.05` or `0.10`, depending on module
- low-count filtering: examples include `rowSums(counts(dds) >= 10) >= 2`
- volcano plot fold-change label threshold: `abs(log2FoldChange) > 1`

Record any threshold changes in your project notes so results remain reproducible.

## 6. Save Results

Recommended output organization:

```text
results/
├── conventional/
├── sports/
│   ├── miRNA/
│   ├── tsRNA/
│   └── rsRNA/
├── tDRMapper/
└── proTRAC/
```

Generated outputs are ignored by Git by default because they can be large and project-specific.
