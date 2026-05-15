# Small Non-coding RNA Analysis Pipeline

A bioinformatics workflow for differential expression analysis of small non-coding RNAs, including piRNAs, tsRNAs, tRNA-derived fragments, miRNAs, rsRNAs, and related small RNA classes.

This repository includes R analysis templates for four complementary workflows:

- conventional reference-based small RNA mapping
- SPORTS summary-based profiling for miRNA, tsRNA, and rsRNA groups
- tDRMapper tRNA-derived RNA count matrices
- proTRAC piRNA cluster count matrices

The scripts are templates. Replace placeholder paths, sample names, treatment names, reference files, and count-matrix filenames with values from your own project before running.

## Repository Layout

```text
.
|-- scripts/
|   |-- Conventional_DE.R
|   |-- SPORTS_DE.R
|   |-- proTRAC_DE.R
|   `-- tDRMapper_DE.R
|-- docs/
|   |-- workflow.md
|   `-- tool_documentation/
|-- examples/
|   |-- metadata_template.csv
|   |-- sports_samples_template.tsv
|   `-- targets_template.tsv
|-- config/
|   `-- config.example.yml
|-- LICENSE
`-- README.md
```

## Pipeline Overview

1. Prepare small RNA FASTQ files and sample metadata.
2. Run the appropriate upstream tool or mapping workflow.
3. Build count matrices for each small RNA feature class.
4. Run DESeq2 differential expression analysis.
5. Export total and significant result tables.
6. Generate volcano plots, MA plots, and heatmaps.

## Workflow Modules

| Module | Main input | Main output |
| --- | --- | --- |
| Conventional | FASTQ files, reference FASTA, feature-count targets table | DE tables, MA plot, volcano plot, heatmap |
| SPORTS | `*_summary.txt` files from SPORTS | miRNA/tsRNA/rsRNA count matrices and DE tables |
| tDRMapper | `tRNA_count_matrix.txt` | tRNA DE tables, volcano plot, heatmap |
| proTRAC | `piRNA_cluster_count_matrix_renamed.txt` | piRNA cluster DE tables, MA plot, volcano plot, heatmap |

Detailed workflow notes are in [docs/workflow.md](docs/workflow.md).

## Requirements

Core command-line tools used by upstream workflows:

- Perl 5.x or later
- Python 3.7 or later
- R 4.0 or later
- Bowtie 1.x
- SAMtools 1.10 or later
- BEDtools 2.30 or later

R packages:

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

Example R install commands:

```r
install.packages(c("data.table", "dplyr", "ggplot2", "ggrepel", "pheatmap", "stringr", "writexl"))

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

BiocManager::install(c("DESeq2", "Rsubread", "Biostrings"))
```

Upstream tools:

- proTRAC v2.4.3
- tDRMapper
- SPORTS v1.1

See the PDFs in [docs/tool_documentation](docs/tool_documentation) for tool-specific installation and run instructions.

## Quick Start

Clone the repository:

```bash
git clone https://github.com/dlld975/Small-non-coding-RNA-analysis-pipeline.git
cd Small-non-coding-RNA-analysis-pipeline
```

Copy the example configuration and metadata templates:

```bash
cp config/config.example.yml config/config.yml
cp examples/metadata_template.csv metadata.csv
```

Edit the placeholders in the script you plan to run:

- `your_directory_here`
- `your_output_path`
- `your_summary_directory_path`
- `your_TREATMENT`
- `your_treatment`
- `CONTROL`
- `your_fasta.fa`
- FASTQ, FASTA, count matrix, summary, or target filenames

Run an analysis module from R or RStudio:

```r
source("scripts/tDRMapper_DE.R")
```

## Typical Inputs

This pipeline is designed for small RNA sequencing data, usually starting from `.fastq`, `.fastq.gz`, or `.fastqsanger.gz` files. The organism and reference databases can be changed by updating the reference paths and annotation files in the relevant script.

Recommended project metadata fields:

- `sample`
- `condition`
- `file`

Keep count matrix column names and metadata row names in the same order before running DESeq2.

## Output Files

Common outputs include:

- `Total_expressed_*.csv`
- `differentially_expressed_*.csv`
- count matrices
- volcano plots
- MA plots
- heatmaps
- saved R objects such as `dds.RData` and `results.RData`

Large raw data files and generated analysis outputs are ignored by Git. Keep FASTQ, BAM, count-output folders, and generated figures outside version control unless you intentionally add small example files.

## Documentation

Tool documentation PDFs are stored in [docs/tool_documentation](docs/tool_documentation):

- `proTRAC_documentation.pdf`
- `SPORTS_documentation.pdf`
- `tDRMapper_documentation.pdf`

## Citation

If you use this pipeline, cite the upstream tools used in your analysis, including DESeq2, Rsubread, SPORTS, tDRMapper, and proTRAC where applicable.

## Authors

Da Lu, Anthony Hannan Group, Hannan Laboratory, The Florey Institute of Neuroscience and Mental Health.

## License

This repository is released under the MIT License. See [LICENSE](LICENSE).
