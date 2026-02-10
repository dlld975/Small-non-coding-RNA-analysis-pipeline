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

---

## Tool-Specific Workflows

### proTRAC - piRNA Analysis

#### Overview
proTRAC identifies and characterizes piRNA clusters from small RNA sequencing data. This module handles piRNA-specific analysis including cluster detection, quantification, and differential expression.

#### What It Does
1. **Collapses** redundant reads to reduce computational load
2. **Maps** collapsed reads to the reference genome
3. **Identifies** piRNA clusters based on mapping patterns
4. **Quantifies** piRNA expression across samples
5. **Performs** differential expression analysis between conditions

#### Key Steps

**1. Environment Setup**
```bash
conda activate pilfer_env  # For read collapse
conda activate bioinfo     # For mapping and analysis
```

**2. Read Collapse**
Uses PILFER's collapse.py to reduce identical sequences:
- Input: Raw FASTQ files
- Output: Collapsed FASTA files
- Reduces file size and processing time

**3. Genome Mapping**
Maps collapsed reads using Bowtie1:
- Build genome index once (time-consuming)
- Map each sample with mismatch tolerance (-v 2)
- Outputs SAM/BAM alignment files

**4. piRNA Cluster Identification**
Runs proTRAC.pl with key parameters:
- `-pimin 21 -pimax 33`: piRNA size range (21-33 nucleotides)
- `-1Tor10A 0.75`: 1T or 10A bias (typical piRNA signature)
- `-clstrand 0.75`: Strand clustering threshold
- Generates cluster GTF files for each sample

**5. Cluster Merging & Quantification**
- Combines cluster GTF files from all samples
- Creates unified cluster coordinates
- Generates count matrix using bedtools intersect
- Each row = one piRNA cluster, each column = one sample

**6. Differential Expression**
- Uses DESeq2 or edgeR in R (see `Conventional_DE.R`)
- Identifies upregulated/downregulated clusters
- Outputs FDR-corrected statistics

#### Key Outputs
- `clusters.gtf`: Identified piRNA cluster coordinates
- `piRNA_cluster_count_matrix.txt`: Read counts per cluster per sample
- `differentially_expressed_clusters.csv`: DE results with statistics
- HTML reports with cluster visualizations

#### Alternative: Condition-Specific Clusters
Instead of differential expression, you can:
1. Merge treatment samples separately from controls
2. Identify clusters unique to each condition
3. Compare overlapping vs. condition-specific clusters

This approach finds clusters present in one condition but absent in the other.

#### Target Gene Prediction
Intersect DE clusters with gene annotations to identify potential piRNA targets:
- Convert cluster coordinates to BED format
- Intersect with gene GTF using bedtools
- Generate list of genes overlapping/near piRNA clusters

---

### SPORTS - Comprehensive Small RNA Profiling

#### Overview
SPORTS (Small RNA-seq Portal for Organizing Reference TRanscriptome Sequences) provides comprehensive annotation and profiling of small RNAs across multiple categories: miRNA, piRNA, tRNA, rRNA, and other non-coding RNAs.

#### What It Does
1. **Annotates** small RNA reads against multiple reference databases
2. **Classifies** reads into different small RNA categories
3. **Quantifies** expression levels for each RNA type
4. **Generates** summary statistics and length distributions
5. **Creates** publication-ready visualizations

#### Key Steps

**1. Reference Database Setup**
Build Bowtie1 indexes for multiple databases:
- **Genome**: Full reference genome (mm10/GRCm39)
- **miRNA**: miRBase database
- **rRNA**: Ribosomal RNA sequences
- **tRNA**: Transfer RNA sequences (mature + mt)
- **piRNA**: piRBase database
- **Ensembl ncRNA**: Other non-coding RNAs
- **Rfam**: RNA families database

**Note**: Build indexes once, reuse for all samples

**2. Input Preparation**
- Convert `.fastqsanger` to `.fastq` (hardlink)
- Decompress if needed (`.gz` files)
- Verify file integrity (non-zero size)

**3. Run SPORTS**
```bash
sports.pl \
  -i input.fastq \
  -p 4                    # threads
  -M 1                    # allow 1 mismatch
  -k                      # keep intermediate files
  -g genome_index \
  -m miRNA_index \
  -r rRNA_index \
  -t tRNA_index \
  -w piRNA_index \
  -e ensembl_index \
  -f rfam_index \
  -o output_directory
```

**4. Batch Processing**
Process multiple samples in loop:
- Creates separate output folder per sample
- Generates individual logs for troubleshooting
- Automatically checks for successful completion

**5. Results Collection**
Collect summary files from all samples:
- `{sample}_summary.txt`: Overall mapping statistics
- `{sample}_output.txt`: Detailed annotations
- `{sample}_length_distribution.txt`: Size profiles

#### Key Outputs

**Summary Statistics**
- Total reads processed
- % mapped to each RNA category
- % unmapped reads
- Read count per category

**Detailed Files**
- **Annotated reads**: Which database each read mapped to
- **Length distributions**: Size profiles for each RNA type
- **Category counts**: Expression levels per RNA class
- **HTML reports**: Visual summaries (if generated)

#### Workflow Logic
SPORTS uses hierarchical mapping:
1. Reads map to miRNA first (highest priority)
2. Unmapped reads → rRNA
3. Still unmapped → tRNA
4. Continue through piRNA, ncRNA, etc.
5. Final unmapped reads reported separately

This prevents ambiguous assignments across categories.

#### Differential Expression
After SPORTS completes:
1. Extract counts for your RNA type of interest
2. Create count matrix (rows = features, columns = samples)
3. Run DE analysis in R using DESeq2/edgeR
4. See `SPORTSDE.R` script for implementation

#### Tips
- **Storage**: Keep intermediate files (`-k`) for debugging, but they consume disk space
- **Threads**: Use `-p 4` or more for faster processing
- **Memory**: ~4-8 GB RAM per sample typical
- **Mismatches**: `-M 1` balances sensitivity and specificity
- **Database order**: Hierarchy matters—miRNA before piRNA prevents misclassification

#### Quality Checks
Before running SPORTS:
- Verify all index files exist (`.ebwt` extensions)
- Check input FASTQ is non-empty
- Ensure sufficient disk space (outputs ~2-5x input size)

After SPORTS:
- Check `run.log` for errors
- Verify `*_summary.txt` was generated
- Confirm mapping percentages are reasonable:
  - rRNA typically 10-40% (varies by sample prep)
  - miRNA typically 5-30%
  - High unmapped (>70%) suggests issues

#### Typical Runtime
- Index building: 1-2 hours (one-time)
- Per sample processing: 30-90 minutes
- Batch of 10 samples: 5-15 hours

#### Common Issues
- **"Index not found"**: Check paths to `.ebwt` files
- **Empty output**: Verify input `.fastq` hardlink exists
- **Low mapping**: May need adapter trimming first
- **Memory errors**: Reduce thread count or process fewer samples simultaneously

---

### tDRMapper - tRNA-Derived Fragment Analysis

#### Overview
tDRMapper (tRNA-Derived Fragment Mapper) specializes in identifying and quantifying tRNA-derived small RNAs (tDRs), also called tRNA fragments or tRFs. These fragments are increasingly recognized as important regulatory molecules.

#### What It Does
1. **Maps** reads to tRNA reference sequences
2. **Classifies** fragments by tRNA type and fragment position
3. **Quantifies** tDR expression across samples
4. **Distinguishes** between different tDR subtypes (5'-tRFs, 3'-tRFs, etc.)

#### Key Steps

**1. Setup**
Clone tDRMapper from GitHub:
```bash
git clone https://github.com/sararselitsky/tDRmapper
```

Prepare tRNA reference:
- Use mature + precursor tRNA sequences
- Format: `mm10_mature_pre_for_tdrMapper.fa`
- Available from GtRNAdb

**2. Run Mapping**
```bash
conda activate bioinfo

perl TdrMappingScripts.pl \
  /path/to/mm10_mature_pre_for_tdrMapper.fa \
  your_sample.fastqsanger.gz
```

**Note**: Can run directly on compressed `.gz` files

**3. Process All Samples**
Loop through all samples in your experiment:
- Treatment samples
- Control samples
- Each generates its own output file

**4. Build Count Matrix**
Combine all samples into a single count matrix:
- Rows: tRNA IDs (which tRNA the fragment derives from)
- Columns: Sample names
- Values: Read counts

Uses AWK script to:
- Parse `.hq_cs.mapped.speciesInfo.txt` output files
- Aggregate counts by tRNA ID
- Handle missing values (fill with 0)

**5. Differential Expression**
- Import count matrix into R
- Run DESeq2 or edgeR analysis
- See `tDRMapper_DE.R` script

#### Key Outputs

**Per Sample**
- `{sample}.hq_cs.mapped.speciesInfo.txt`: Main results file
  - Fragment sequences
  - tRNA IDs
  - Fragment counts
  - Fragment positions

**Combined**
- `tRNA_count_matrix.txt`: Count table for all samples
  - Ready for statistical analysis
  - Tab-delimited format
  - Header row with sample names

#### Understanding tDR Types
tDRMapper can distinguish:
- **5'-tRFs**: Fragments from 5' end of tRNA
- **3'-tRFs**: Fragments from 3' end of tRNA
- **i-tRFs**: Internal fragments
- **tRNA halves**: ~30-35 nt fragments from stress-induced cleavage

#### Count Matrix Building Notes
The AWK script for matrix building:
1. Reads all `.speciesInfo.txt` files
2. Extracts tRNA ID (column 2) and count (column 4)
3. Matches filename to sample name
4. Aggregates counts per tRNA per sample
5. Fills missing values with 0

**Important**: Update sample name mappings in the AWK script to match your filenames.

#### Quality Checks
After mapping:
- Check output file sizes (non-zero)
- Verify reasonable mapping rates (typically 5-20% for tDRs)
- Examine length distributions (should peak at tDR sizes: 15-35 nt)
- Look for tRNA bias (some tRNAs produce more fragments)

Before DE analysis:
- Verify count matrix dimensions (tRNAs × samples)
- Check for zero-sum rows (tRNAs with no reads)
- Confirm sample order matches metadata

#### Typical Runtime
- Per sample: 5-15 minutes
- Batch of 10 samples: 1-2 hours
- Count matrix building: <1 minute

#### Common Issues
- **Low mapping**: Check adapter trimming quality
- **Missing output**: Verify Perl dependencies installed
- **Count matrix errors**: Check sample name matching in AWK script
- **Reference errors**: Ensure tRNA reference is properly formatted

#### Integration with Other Tools
tDRMapper complements:
- **SPORTS**: SPORTS quantifies tRNAs; tDRMapper quantifies tRNA *fragments*
- **proTRAC**: Different fragment classes—piRNAs vs tDRs
- Use all three for comprehensive small RNA profiling

---

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

Hannan Group - Da Lu

## Acknowledgments

This pipeline integrates several published bioinformatics tools developed by the scientific community.
