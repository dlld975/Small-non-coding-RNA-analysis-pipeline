# Small Non-Coding RNA Analysis Pipeline

A comprehensive bioinformatics pipeline for analyzing small non-coding RNAs (sncRNAs) including piRNAs, tsRNAs, tRNA-derived fragments (tDRs), miRNAs and other small RNA species from high-throughput sequencing data.

## Overview

This repository contains R scripts and documentation for processing and analyzing small RNA sequencing data using a conventional method as well as using a method with three complementary tools:
- **proTRAC** - piRNA cluster identification and analysis
- **tDRMapper** - tRNA-derived fragment mapping and quantification
- **SPORTS** - Small RNA annotation and profiling

The pipeline performs differential expression analysis and functional characterization of small non-coding RNAs between experimental conditions (e.g., treatment vs. control).

## Data Requirements

This pipeline is designed for:
- **Input data**: RNA-seq FASTQ files
- **Recommended format**: `.fastqsanger` or `.fastq.gz`
- **Organism**: Mouse (*Mus musculus*) - scripts use GRCm39/mm10 reference genome
  - Can be adapted for other organisms by modifying reference paths

### Typical Applications
- tRNA-derived fragment (tDR) profiling
- Comprehensive small RNA characterization
- Comparative analysis between experimental conditions

## Repository Structure

```
.
├── src/                    # R source code for analysis
│   ├── Conventional_DE.R   # Conventional DE analysis
│   └── ...                 # Additional analysis scripts
├── doc/                    # Documentation
│   ├── proTRAC_documentation.pdf
│   ├── SPORTS_documentation.pdf
│   └── tDRMapper_documentation.pdf
└── README.md
```

## System Requirements

### Platform
This pipeline was developed and tested on:
- **Windows Subsystem for Linux (WSL)** or **Linux-based system**
- Windows users should use WSL2 with Ubuntu 20.04 or later

### Dependencies

#### Core Tools
- **Perl** (v5.x or later) - for proTRAC and tDRMapper
- **Python** (v3.7+) - for PILFER collapse script
- **R** (v4.0+) - for statistical analysis
- **Bowtie** (v1.x) - for sequence alignment
- **SAMtools** (v1.10+) - for BAM file processing
- **BEDtools** (v2.30+) - for genomic interval operations

#### Python Packages
```bash
pip install cutadapt
```

#### R Packages
```R
# Install required R packages
install.packages(c("DESeq2", "edgeR", "ggplot2", "dplyr"))
# Bioconductor packages
BiocManager::install(c("DESeq2", "GenomicFeatures", "Rsubread"))
```

#### Additional Software
- **proTRAC** v2.4.3 - Download from [SourceForge](https://sourceforge.net/projects/protrac/)
- **tDRMapper** - Clone from [GitHub](https://github.com/sararselitsky/tDRmapper)
- **SPORTS** v1.1 - Download from [SPORTS repository](https://github.com/junchaoshi/SPORTS1.1)

## Installation

### 1. Clone this repository
```bash
git clone https://github.com/yourusername/sncRNA-analysis-pipeline.git
cd sncRNA-analysis-pipeline
```

### 2. Set up Conda environment
```bash
# Create environment for PILFER
conda create -n pilfer_env python=3.9
conda activate pilfer_env

# Create environment for bioinformatics tools
conda create -n bioinfo
conda activate bioinfo
conda install -c bioconda bowtie samtools bedtools
```

### 3. Install external tools
Follow the detailed installation instructions in the documentation files:
- `doc/proTRAC_documentation.pdf`
- `doc/SPORTS_documentation.pdf`
- `doc/tDRMapper_documentation.pdf`

### 4. Download reference databases
Required reference files for mouse (GRCm39/mm10):
- Genome FASTA
- Gene annotations (GTF/GFF)
- miRBase, piRBase, GtRNAdb, Rfam databases
- See documentation for download links

## Usage

### Preprocessing (Recommended)

Before running the main pipeline, perform quality control and adapter trimming:

```bash
# Quality control with FastQC
fastqc your_sample.fastq.gz -o qc_reports/

# Adapter trimming with Cutadapt
cutadapt -a ADAPTER_SEQUENCE \
  -m 15 -M 35 \
  -o trimmed_sample.fastq.gz \
  your_sample.fastq.gz
```

**Note**: Adapter sequences vary by library preparation kit. Common adapters:
- Illumina TruSeq: `TGGAATTCTCGGGTGCCAAGG`
- NEBNext: `AGATCGGAAGAGCACACGTCT`

### Pipeline Workflow

#### 1. Run proTRAC (piRNA analysis)
```bash
# See doc/proTRAC_documentation.pdf for detailed steps
# Key steps:
# - Collapse reads with PILFER
# - Map to genome with Bowtie
# - Identify piRNA clusters
# - Differential expression in R
```

#### 2. Run SPORTS (comprehensive small RNA profiling)
```bash
# See doc/SPORTS_documentation.pdf for detailed steps
sports.pl -i input.fastq -g genome -m miRNA -r rRNA -t tRNA -w piRNA -o output/
```

#### 3. Run tDRMapper (tRNA-derived fragments)
```bash
# See doc/tDRMapper_documentation.pdf for detailed steps
perl TdrMappingScripts.pl reference.fa input.fastq.gz
```

#### 4. Differential Expression Analysis
```R
# Run R scripts in src/ directory
source("src/Conventional_DE.R")
# Follow script-specific parameters
```

## Output Files

The pipeline generates:
- **Count matrices** - Read counts per feature (piRNA cluster, tRNA, etc.)
- **Differential expression results** - Statistical comparisons with FDR correction
- **Annotation files** - BED/GTF files of identified features
- **Visualizations** - Plots and Heatmaps reports

## Documentation

Detailed step-by-step instructions for each tool are available in the `doc/` folder:
- **proTRAC_documentation.pdf** - piRNA cluster identification workflow
- **SPORTS_documentation.pdf** - Small RNA profiling workflow
- **tDRMapper_documentation.pdf** - tRNA fragment analysis workflow

Each documentation file includes:
- Installation instructions
- Reference database setup
- Complete command-line examples

## Citation

If you use this pipeline, please cite the original tools:

- **proTRAC**: Rosenkranz D, et al. (2015) *Bioinformatics*
- **SPORTS**: Shi J, et al. (2018) *Bioinformatics*
- **tDRMapper**: Selitsky SR, et al. (2015) *RNA*
- and Da Lu for the development of this pipeline/Hannan Lab

## Support

For issues or questions:
- Check the documentation in `doc/` folder
- Review tool-specific manuals and publications
- Open an issue on this GitHub repository

## License

This pipeline is provided as-is for research purposes. Individual tools retain their original licenses.

## Authors

Hannan Group- Da Lu

## Acknowledgments

This pipeline integrates several published bioinformatics tools developed by the scientific community.
