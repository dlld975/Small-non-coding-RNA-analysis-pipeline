# Small Non-coding RNA Analysis Pipeline

This repository contains an analysis workflow for differential expression of small non-coding RNAs from multiple upstream pipelines:

- conventional reference-based small RNA mapping
- SPORTS summary outputs for miRNA, tsRNA, and rsRNA groups
- tDRMapper tRNA-derived RNA count matrices
- proTRAC piRNA cluster count matrices

The scripts are written as R analysis templates. Replace the placeholder paths, treatment names, sample names, and count-matrix filenames with values from your project before running.

## Repository Layout

```text
.
├── scripts/
│   ├── Conventional_DE.R
│   ├── SPORTS_DE.R
│   ├── proTRAC_DE.R
│   └── tDRMapper_DE.R
├── docs/
│   ├── workflow.md
│   └── tool_documentation/
├── examples/
│   ├── metadata_template.csv
│   ├── sports_samples_template.tsv
│   └── targets_template.tsv
├── config/
│   └── config.example.yml
└── README.md
```

## Pipeline Overview

1. Prepare small RNA FASTQ files and sample metadata.
2. Run the appropriate upstream small RNA tool:
   - `Conventional_DE.R` uses Rsubread alignment and featureCounts.
   - `SPORTS_DE.R` summarizes SPORTS class/sub-class outputs.
   - `tDRMapper_DE.R` analyzes tDRMapper tRNA count matrices.
   - `proTRAC_DE.R` analyzes proTRAC piRNA cluster count matrices.
3. Run DESeq2 differential expression analysis.
4. Export total and significant differential expression tables.
5. Generate volcano plots, MA plots, and heatmaps for publication or QC.

## Requirements

Install R and the required R/Bioconductor packages before running the scripts.

Core R packages:

- `data.table`
- `dplyr`
- `ggplot2`
- `ggrepel`
- `pheatmap`
- `stringr`
- `writexl`

Bioconductor packages:

- `DESeq2`
- `Rsubread`
- `Biostrings`

Example install commands:

```r
install.packages(c("data.table", "dplyr", "ggplot2", "ggrepel", "pheatmap", "stringr", "writexl"))

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

BiocManager::install(c("DESeq2", "Rsubread", "Biostrings"))
```

## Quick Start

Clone the repository:

```bash
git clone https://github.com/dlld975/Small-non-coding-RNA-analysis-pipeline.git
cd Small-non-coding-RNA-analysis-pipeline
```

Copy the example config and metadata templates:

```bash
cp config/config.example.yml config/config.yml
cp examples/metadata_template.csv metadata.csv
```

Edit the placeholders in the script you plan to run:

- `your_directory_here`
- `your_TREATMENT`
- `CONTROL`
- input FASTQ, FASTA, count matrix, summary, or target filenames
- output directory paths

Run one analysis module from R or RStudio:

```r
source("scripts/tDRMapper_DE.R")
```

## Inputs

The expected input depends on the module.

| Module | Main input | Main output |
| --- | --- | --- |
| Conventional | FASTQ files, reference FASTA, feature-count targets table | DE tables, MA plot, volcano plot, heatmap |
| SPORTS | `*_summary.txt` files from SPORTS | miRNA/tsRNA/rsRNA count matrices and DE tables |
| tDRMapper | `tRNA_count_matrix.txt` | tRNA DE tables, volcano plot, heatmap |
| proTRAC | `piRNA_cluster_count_matrix_renamed.txt` | piRNA cluster DE tables, MA plot, volcano plot, heatmap |

See [docs/workflow.md](docs/workflow.md) for a more detailed workflow.

## Output Files

Common outputs include:

- `Total_expressed_*.csv`
- `differentially_expressed_*.csv`
- volcano plots
- MA plots
- heatmaps
- saved R objects such as `dds.RData` and `results.RData`

Large raw data files and generated analysis outputs are ignored by Git. Keep FASTQ, BAM, large count outputs, and generated figures outside version control unless you intentionally add small example files.

## Documentation

Tool documentation PDFs are stored in [docs/tool_documentation](docs/tool_documentation):

- tDRMapper documentation
- proTRAC documentation
- SPORTS documentation

## Citation

Please cite the upstream tools used in your analysis, including DESeq2, Rsubread, SPORTS, tDRMapper, and proTRAC where applicable.

## License

No license has been selected yet. Add a license before public reuse if you want others to be able to copy, modify, or redistribute this workflow.
