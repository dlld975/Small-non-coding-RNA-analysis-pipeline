# Small Non-coding RNA Analysis Pipeline

A section-based workflow for small non-coding RNA analysis, including conventional mapping, SPORTS, tDRMapper, proTRAC, cluster differential expression, non-overlapping cluster comparison, and shared R differential expression visualization.

Replace every placeholder such as `yourfile.fastq.gz`, `your_reference.fa`, `your_count_matrix.txt`, `your_output_directory`, `your_TREATMENT`, and `CONTROL` with your own project values.

## Repository Layout

```text
.
|-- README.md
|-- LICENSE
|-- scripts/
|   |-- conventional_de.R
|   |-- sports_run_workflow.sh
|   |-- sports_summary_to_counts.R
|   |-- tdrmapper_run_workflow.sh
|   |-- tdrmapper_de.R
|   |-- protrac_run_prepare_clusters.sh
|   |-- protrac_cluster_de.R
|   |-- general_count_matrix_de.R
|   |-- protrac_non_de_unique_clusters.sh
|   `-- de_plots.R
|-- config/
|   `-- config.example.yml
|-- examples/
|   |-- metadata_template.csv
|   |-- sports_samples_template.tsv
|   `-- targets_template.tsv
`-- docs/
    |-- workflow.md
    `-- tool_documentation/
```

## Requirements

R packages:

```r
install.packages(c(
  "data.table", "dplyr", "ggplot2", "ggrepel",
  "pheatmap", "stringr", "writexl", "seqinr"
))

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

BiocManager::install(c("DESeq2", "Rsubread", "Biostrings"))
```

External tools used by the upstream workflows may include Perl, Python, Bowtie, SAMtools, BEDtools, proTRAC, SPORTS, and tDRMapper. See the PDFs in `docs/tool_documentation/`.

Upstream tool links:

- SPORTS1.1 GitHub: https://github.com/junchaoshi/sports1.1
- tDRmapper GitHub: https://github.com/sararselitsky/tDRmapper
- proTRAC official source page: https://sourceforge.net/projects/protrac/

## R Script Templates

The full workflow code is shown in sections below, and matching R script templates are also provided:

| Section | Script |
| --- | --- |
| Conventional Mapping and DE | `scripts/conventional_de.R` |
| SPORTS Bash workflow | `scripts/sports_run_workflow.sh` |
| SPORTS Summary to Count Matrix | `scripts/sports_summary_to_counts.R` |
| tDRMapper Bash workflow | `scripts/tdrmapper_run_workflow.sh` |
| tDRMapper DE | `scripts/tdrmapper_de.R` |
| proTRAC Bash cluster preparation | `scripts/protrac_run_prepare_clusters.sh` |
| proTRAC Cluster DE | `scripts/protrac_cluster_de.R` |
| General Count Matrix DE | `scripts/general_count_matrix_de.R` |
| proTRAC non-DE unique/no-overlap clusters | `scripts/protrac_non_de_unique_clusters.sh` |
| Volcano, MA, and heatmap plots | `scripts/de_plots.R` |

Edit the user settings at the top of each script before running.

## Bash Then R DE Workflow Order

For the tool-specific modules, run the upstream command-line workflow first, then run the R DE script:

| Tool | Step 1: Bash/upstream workflow | Step 2: R DE |
| --- | --- | --- |
| SPORTS | `scripts/sports_run_workflow.sh` creates per-sample SPORTS outputs and `*_summary.txt` files | `scripts/sports_summary_to_counts.R`, then `scripts/general_count_matrix_de.R` for miRNA, tsRNA, or rsRNA matrices |
| tDRMapper | `scripts/tdrmapper_run_workflow.sh` runs tDRMapper and builds `your_tRNA_count_matrix.txt` | `scripts/tdrmapper_de.R` |
| proTRAC | `scripts/protrac_run_prepare_clusters.sh` runs proTRAC and prepares a piRNA cluster count matrix | `scripts/protrac_cluster_de.R` |
| proTRAC unique clusters | `scripts/protrac_non_de_unique_clusters.sh` merges treatment/control GTFs separately and finds no-overlap clusters | No DESeq2 step; use unique BED files for gene intersection and functional study |

## 1. Conventional Mapping and DE

Use this section when mapping FASTQ reads directly to a custom small RNA reference FASTA.

```r
setwd("your_working_directory")

library(Rsubread)
library(seqinr)
library(Biostrings)
library(DESeq2)
library(writexl)

# Build a reference index.
buildindex(
  basename = "your_index_name",
  reference = "your_reference.fa"
)

# Map treatment and control FASTQ files.
subjunc(index = "your_index_name", readfile1 = "your_treatment_file_1.fastq.gz", output_file = "your_treatment_1.bam")
subjunc(index = "your_index_name", readfile1 = "your_treatment_file_2.fastq.gz", output_file = "your_treatment_2.bam")
subjunc(index = "your_index_name", readfile1 = "your_control_file_1.fastq.gz", output_file = "your_control_1.bam")
subjunc(index = "your_index_name", readfile1 = "your_control_file_2.fastq.gz", output_file = "your_control_2.bam")

# Build a simple annotation table from the FASTA.
reference_sequences <- readDNAStringSet("your_reference.fa")
annotation <- data.frame(
  GeneID = names(reference_sequences),
  Chr = names(reference_sequences),
  start = 0,
  end = width(reference_sequences),
  strand = "*"
)

# Target table columns should include File_name and Treatment.
targets <- read.delim("your_targets_table.txt", check.names = FALSE)

FC <- featureCounts(
  files = targets$File_name,
  annot.ext = annotation,
  countMultiMappingReads = FALSE,
  useMetaFeatures = TRUE
)

dds <- DESeqDataSetFromMatrix(
  countData = FC$counts,
  colData = targets,
  design = ~ Treatment
)

# Filter low-count features before DESeq2.
# min_count_per_feature = minimum raw count required in a sample.
# min_samples_with_count = minimum number of samples that must reach that count.
min_count_per_feature <- your_min_count_threshold
min_samples_with_count <- your_min_sample_threshold
dds <- dds[rowSums(counts(dds) >= min_count_per_feature) >= min_samples_with_count, ]
dds <- DESeq(dds)

res <- results(dds, contrast = c("Treatment", "your_TREATMENT", "CONTROL"))
sig_res <- res[which(res$padj < 0.05 & !is.na(res$padj)), ]

colnames(res)[colnames(res) == "padj"] <- "adj.P.Val"
colnames(sig_res)[colnames(sig_res) == "padj"] <- "adj.P.Val"
res$Gene <- rownames(res)
sig_res$Gene <- rownames(sig_res)

res_sorted <- res[order(res$adj.P.Val), ]
sig_res_sorted <- sig_res[order(sig_res$adj.P.Val), ]

write_xlsx(as.data.frame(res_sorted), path = "your_total_DE_results.xlsx")
write_xlsx(as.data.frame(sig_res_sorted), path = "your_significant_DE_results.xlsx")
save(dds, res, sig_res, res_sorted, sig_res_sorted, file = "your_DE_results.RData")
```

## 2. SPORTS Summary to Count Matrix

Use this section after SPORTS has generated one `*_summary.txt` file per sample. SPORTS1.1 is available at https://github.com/junchaoshi/sports1.1.

Upstream Bash step:

```bash
bash scripts/sports_run_workflow.sh
```

Then build count matrices with R:

```r
setwd("your_working_directory")

library(data.table)
library(stringr)

summary_dir <- "your_sports_summary_directory"
files <- list.files(summary_dir, pattern = "_summary\\.txt$", full.names = TRUE)

sample_names <- str_replace(basename(files), "_summary\\.txt$", "")
sample_table <- data.frame(
  sample = sample_names,
  file = files,
  condition = ifelse(str_detect(sample_names, "^your_treatment_prefix"), "your_TREATMENT", "CONTROL"),
  stringsAsFactors = FALSE
)

sample_table <- sample_table[order(sample_table$condition, sample_table$sample), ]

read_summary <- function(file) {
  dt <- fread(file, sep = "\t", header = TRUE, data.table = TRUE)
  setnames(dt, tolower(names(dt)))
  dt[, reads := as.numeric(reads)]
  dt[is.na(sub_class) | sub_class == "", sub_class := "-"]
  dt
}

assign_group <- function(class_str) {
  x <- tolower(class_str)
  if (str_detect(x, "mirbase.*mirna")) return("miRNA")
  if (str_detect(x, "trna")) return("tsRNA")
  if (str_detect(x, "rrna|yrna|rny|\\b12s\\b|\\b16s\\b|\\b18s\\b|\\b28s\\b|5\\.8s|\\b5s\\b")) return("rsRNA")
  "other"
}

all_dt <- rbindlist(lapply(seq_len(nrow(sample_table)), function(i) {
  dt <- read_summary(sample_table$file[i])
  dt[, sample := sample_table$sample[i]]
  dt[, group := vapply(class, assign_group, character(1))]
  dt[, feature := sub_class]
  dt[, counts := as.integer(round(reads))]
  dt
}), fill = TRUE)

make_matrix <- function(group_name) {
  dtg <- all_dt[group == group_name]
  dtg <- dtg[feature != "-" & !is.na(feature)]
  dtg <- dtg[, .(counts = sum(counts, na.rm = TRUE)), by = .(feature, sample)]

  mat <- dcast(dtg, feature ~ sample, value.var = "counts", fill = 0)
  rn <- mat$feature
  mat$feature <- NULL

  m <- as.matrix(mat)
  rownames(m) <- rn
  m <- m[, sample_table$sample, drop = FALSE]
  storage.mode(m) <- "integer"
  m
}

out_dir <- "your_output_directory"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

write.table(make_matrix("miRNA"), file.path(out_dir, "your_miRNA_counts.txt"), sep = "\t", quote = FALSE, col.names = NA)
write.table(make_matrix("tsRNA"), file.path(out_dir, "your_tsRNA_counts.txt"), sep = "\t", quote = FALSE, col.names = NA)
write.table(make_matrix("rsRNA"), file.path(out_dir, "your_rsRNA_counts.txt"), sep = "\t", quote = FALSE, col.names = NA)
```

## 3. tDRMapper DE

Use this section for tDRMapper tRNA-derived RNA count matrices. tDRmapper is available at https://github.com/sararselitsky/tDRmapper.

Upstream Bash step:

```bash
bash scripts/tdrmapper_run_workflow.sh
```

Then run R DE on the generated tRNA count matrix:

```r
setwd("your_tdrmapper_results_directory")

library(DESeq2)

count_data <- read.csv(
  "your_tRNA_count_matrix.txt",
  header = TRUE,
  sep = "\t",
  row.names = 1,
  check.names = FALSE
)

col_data <- data.frame(
  Treatment = factor(c(rep("CONTROL", 5), rep("your_TREATMENT", 5)), levels = c("CONTROL", "your_TREATMENT")),
  row.names = colnames(count_data)
)

stopifnot(all(colnames(count_data) == rownames(col_data)))

dds <- DESeqDataSetFromMatrix(
  countData = count_data,
  colData = col_data,
  design = ~ Treatment
)

# Filter low-count tDR/tRNA features before DESeq2.
# Choose these thresholds based on sequencing depth and replicate number.
min_count_per_feature <- your_min_count_threshold
min_samples_with_count <- your_min_sample_threshold
dds <- dds[rowSums(counts(dds) >= min_count_per_feature) >= min_samples_with_count, ]
dds <- DESeq(dds)

res <- results(dds, contrast = c("Treatment", "your_TREATMENT", "CONTROL"))
sig_res <- res[which(res$padj < 0.05 & !is.na(res$padj)), ]

colnames(res)[colnames(res) == "padj"] <- "adj.P.Val"
colnames(sig_res)[colnames(sig_res) == "padj"] <- "adj.P.Val"
res$Gene <- rownames(res)
sig_res$Gene <- rownames(sig_res)

write.csv(as.data.frame(res), "your_total_tDRMapper_DE_results.csv", row.names = TRUE)
write.csv(as.data.frame(sig_res), "your_significant_tDRMapper_DE_results.csv", row.names = TRUE)
```

## 4. proTRAC Cluster DE

Use this section for proTRAC piRNA cluster count matrices. proTRAC is officially distributed at https://sourceforge.net/projects/protrac/.

Upstream Bash step:

```bash
bash scripts/protrac_run_prepare_clusters.sh
```

Then run R DE on the generated piRNA cluster count matrix:

```r
setwd("your_protrac_results_directory")

library(DESeq2)

count_data <- read.csv(
  "your_piRNA_cluster_count_matrix.txt",
  header = TRUE,
  sep = "\t",
  row.names = 1,
  check.names = FALSE
)

col_data <- data.frame(
  Treatment = factor(c(rep("your_TREATMENT", 5), rep("CONTROL", 5)), levels = c("CONTROL", "your_TREATMENT")),
  row.names = colnames(count_data)
)

stopifnot(all(colnames(count_data) == rownames(col_data)))

dds <- DESeqDataSetFromMatrix(
  countData = count_data,
  colData = col_data,
  design = ~ Treatment
)

# Filter low-count piRNA clusters before DESeq2.
# Choose these thresholds based on sequencing depth and replicate number.
min_count_per_cluster <- your_min_count_threshold
min_samples_with_count <- your_min_sample_threshold
dds <- dds[rowSums(counts(dds) >= min_count_per_cluster) >= min_samples_with_count, ]
dds <- DESeq(dds)

res <- results(dds, contrast = c("Treatment", "your_TREATMENT", "CONTROL"))
sig_res <- res[which(res$padj < 0.10 & !is.na(res$padj)), ]

colnames(res)[colnames(res) == "padj"] <- "adj.P.Val"
colnames(sig_res)[colnames(sig_res) == "padj"] <- "adj.P.Val"
res$Gene <- rownames(res)
sig_res$Gene <- rownames(sig_res)

write.csv(as.data.frame(res), "your_total_cluster_DE_results.csv", row.names = TRUE)
write.csv(as.data.frame(sig_res), "your_significant_cluster_DE_results.csv", row.names = TRUE)
```

## 5. General Count Matrix DE

Use this reusable function for any count matrix, including SPORTS miRNA, tsRNA, rsRNA, tDRMapper, conventional mapping, or cluster-level matrices.

```r
library(DESeq2)

run_deseq2 <- function(
  count_file,
  metadata_file,
  output_prefix,
  treatment_label = "your_TREATMENT",
  control_label = "CONTROL",
  padj_cutoff = 0.05,
  min_count = your_min_count_threshold,
  min_samples = your_min_sample_threshold
) {
  count_data <- read.csv(
    count_file,
    sep = "\t",
    header = TRUE,
    row.names = 1,
    check.names = FALSE
  )

  metadata <- read.csv(metadata_file, row.names = 1, check.names = FALSE)
  metadata$Treatment <- factor(metadata$Treatment, levels = c(control_label, treatment_label))

  count_data <- count_data[, rownames(metadata), drop = FALSE]
  stopifnot(all(colnames(count_data) == rownames(metadata)))

  dds <- DESeqDataSetFromMatrix(
    countData = round(count_data),
    colData = metadata,
    design = ~ Treatment
  )

  # Filter low-count features before DESeq2.
  # min_count = minimum raw count required in a sample.
  # min_samples = minimum number of samples that must reach min_count.
  dds <- dds[rowSums(counts(dds) >= min_count) >= min_samples, ]
  dds <- DESeq(dds)

  res <- results(dds, contrast = c("Treatment", treatment_label, control_label))
  sig_res <- res[which(res$padj < padj_cutoff & !is.na(res$padj)), ]

  colnames(res)[colnames(res) == "padj"] <- "adj.P.Val"
  colnames(sig_res)[colnames(sig_res) == "padj"] <- "adj.P.Val"
  res$Gene <- rownames(res)
  sig_res$Gene <- rownames(sig_res)

  write.csv(as.data.frame(res), paste0(output_prefix, "_total_DE.csv"), row.names = TRUE)
  write.csv(as.data.frame(sig_res), paste0(output_prefix, "_significant_DE.csv"), row.names = TRUE)

  list(dds = dds, res = res, sig_res = sig_res)
}

de <- run_deseq2(
  count_file = "your_count_matrix.txt",
  metadata_file = "your_metadata.csv",
  output_prefix = "your_output_name"
)
```

## 6. proTRAC Non-DE Unique Cluster Workflow

This is not a DESeq2 differential-expression step. Use this workflow when you want to merge treatment proTRAC clusters into one combined GTF, merge control proTRAC clusters into one combined GTF, and then identify clusters that have no overlap in the opposite group.

The key outputs are:

- `your_TREATMENT_unique_noOverlap.bed`: treatment-specific clusters with no overlap in control. Use these for gene intersection and functional study as treatment-specific/upregulated-like clusters.
- `CONTROL_unique_noOverlap.bed`: control-specific clusters with no overlap in treatment. Use these for gene intersection and functional study as control-specific/downregulated-like clusters.
- `your_TREATMENT_clusters_withGenes.tsv` and `CONTROL_clusters_withGenes.tsv`: unique clusters annotated with overlapping genes.

The full reusable script is `scripts/protrac_non_de_unique_clusters.sh`.

```bash
conda activate bioinfo
cd "your_protrac_project_directory"

treatment_label="your_TREATMENT"
control_label="CONTROL"
gene_gtf="your_gene_annotation.gtf"

# 1. Merge treatment and control proTRAC cluster GTF files separately.
cat proTRAC_your_treatment_sample*/clusters.gtf > "${treatment_label}combined.gtf"
cat proTRAC_your_control_sample*/clusters.gtf > "${control_label}combined.gtf"

# 2. Add gene_id to proTRAC cluster GTF attributes.
awk 'BEGIN{FS=OFS="\t"} {
  gsub(/.*piRNA cluster no: ([^;]+);.*/, "gene_id \"cluster_"$9"\";", $9)
  print
}' "${treatment_label}combined.gtf" > "${treatment_label}modified.gtf"

awk 'BEGIN{FS=OFS="\t"} {
  gsub(/.*piRNA cluster no: ([^;]+);.*/, "gene_id \"cluster_"$9"\";", $9)
  print
}' "${control_label}combined.gtf" > "${control_label}modified.gtf"

# 3. Convert merged GTF files to BED6.
awk 'BEGIN{FS=OFS="\t"} $3=="piRNA_cluster" {print $1, $4-1, $5, $9, ".", $7}' \
  "${treatment_label}modified.gtf" > "${treatment_label}modified.bed"

awk 'BEGIN{FS=OFS="\t"} $3=="piRNA_cluster" {print $1, $4-1, $5, $9, ".", $7}' \
  "${control_label}modified.gtf" > "${control_label}modified.bed"

sort -k1,1 -k2,2n "${treatment_label}modified.bed" > "${treatment_label}modified.sorted.bed"
sort -k1,1 -k2,2n "${control_label}modified.bed" > "${control_label}modified.sorted.bed"

# 4. Check chromosome naming before overlap.
echo "${treatment_label} chroms:"
cut -f1 "${treatment_label}modified.sorted.bed" | sort -u

echo "${control_label} chroms:"
cut -f1 "${control_label}modified.sorted.bed" | sort -u

echo "Chroms in ${treatment_label} but not in ${control_label}:"
comm -23 \
  <(cut -f1 "${treatment_label}modified.sorted.bed" | sort -u) \
  <(cut -f1 "${control_label}modified.sorted.bed" | sort -u)

echo "Chroms in ${control_label} but not in ${treatment_label}:"
comm -13 \
  <(cut -f1 "${treatment_label}modified.sorted.bed" | sort -u) \
  <(cut -f1 "${control_label}modified.sorted.bed" | sort -u)

# 5. Classify any overlap as WITHIN / INCLUDES / OVERLAP_ONLY.
bedtools intersect \
  -a "${treatment_label}modified.sorted.bed" \
  -b "${control_label}modified.sorted.bed" \
  -wa -wb \
  > "${treatment_label}_vs_${control_label}.overlap.tsv"

awk -v treatment="$treatment_label" -v control="$control_label" 'BEGIN{FS=OFS="\t"}
{
  aS=$2; aE=$3
  bS=$8; bE=$9
  rel="OVERLAP_ONLY"
  if (aS>=bS && aE<=bE) rel=treatment"_WITHIN_"control
  else if (aS<=bS && aE>=bE) rel=treatment"_INCLUDES_"control
  print $1,aS,aE,rel,$4,bS,bE,$10
}' "${treatment_label}_vs_${control_label}.overlap.tsv" \
  > "${treatment_label}_vs_${control_label}.classified.tsv"

# 6. Main no-overlap outputs for downstream gene intersection.
bedtools intersect \
  -a "${treatment_label}modified.sorted.bed" \
  -b "${control_label}modified.sorted.bed" \
  -v \
  > "${treatment_label}_unique_noOverlap.bed"

bedtools intersect \
  -a "${control_label}modified.sorted.bed" \
  -b "${treatment_label}modified.sorted.bed" \
  -v \
  > "${control_label}_unique_noOverlap.bed"

awk -v rel="${treatment_label}_WITHIN_${control_label}" '$4==rel' \
  "${treatment_label}_vs_${control_label}.classified.tsv" \
  > "${treatment_label}_within_${control_label}.tsv"

awk -v rel="${treatment_label}_INCLUDES_${control_label}" '$4==rel' \
  "${treatment_label}_vs_${control_label}.classified.tsv" \
  > "${treatment_label}_includes_${control_label}.tsv"

# 7. Convert gene annotation GTF to BED6.
awk 'BEGIN{FS=OFS="\t"} $3=="gene" {
  g=$9
  sub(/.*gene_name "/, "", g)
  sub(/".*/, "", g)
  print $1, $4-1, $5, g, ".", $7
}' "$gene_gtf" | sort -k1,1 -k2,2n > genes.bed

# 8. Intersect treatment unique clusters with genes.
# If cluster strand is ".", do not use -s. Add -s only when strand is reliable.
bedtools intersect \
  -a "${treatment_label}_unique_noOverlap.bed" \
  -b genes.bed \
  -wa -wb \
  > "${treatment_label}_intersected_genes.bed"

awk 'BEGIN{FS=OFS="\t"}
{
  key=$1":"$2"-"$3":"$6
  gene=$10
  if (gene=="") gene="NA"
  if (!seen[key,gene]++) {
    genes[key] = (genes[key]=="" ? gene : genes[key]","gene)
  }
}
END{
  for (k in genes) print k, genes[k]
}' "${treatment_label}_intersected_genes.bed" > "${treatment_label}_genes_collapsed.tsv"

awk 'BEGIN{FS=OFS="\t";
  print "Cluster_ID","Chromosome","Start","End","Genes","Score","Strand"
}
ARGIND==1 {g[$1]=$2; next}
ARGIND==2 {
  key=$1":"$2"-"$3":"$6
  cid="cluster_" FNR
  gene_list=(key in g ? g[key] : "NA")
  print cid,$1,$2,$3,gene_list,$5,$6
}' "${treatment_label}_genes_collapsed.tsv" "${treatment_label}_unique_noOverlap.bed" \
  > "${treatment_label}_clusters_withGenes.tsv"

# 9. Intersect control unique clusters with genes.
bedtools intersect \
  -a "${control_label}_unique_noOverlap.bed" \
  -b genes.bed \
  -wa -wb \
  > "${control_label}_intersected_genes.bed"

awk 'BEGIN{FS=OFS="\t"}
{
  key=$1":"$2"-"$3":"$6
  gene=$10
  if (gene=="") gene="NA"
  if (!seen[key,gene]++) {
    genes[key] = (genes[key]=="" ? gene : genes[key]","gene)
  }
}
END{
  for (k in genes) print k, genes[k]
}' "${control_label}_intersected_genes.bed" > "${control_label}_genes_collapsed.tsv"

awk 'BEGIN{FS=OFS="\t";
  print "Cluster_ID","Chromosome","Start","End","Genes","Score","Strand"
}
ARGIND==1 {g[$1]=$2; next}
ARGIND==2 {
  key=$1":"$2"-"$3":"$6
  cid="cluster_" FNR
  gene_list=(key in g ? g[key] : "NA")
  print cid,$1,$2,$3,gene_list,$5,$6
}' "${control_label}_genes_collapsed.tsv" "${control_label}_unique_noOverlap.bed" \
  > "${control_label}_clusters_withGenes.tsv"
```

Notes:

- Use `*_unique_noOverlap.bed` files for downstream gene intersection and functional study.
- Treatment unique clusters can be interpreted as treatment-specific/upregulated-like clusters, and control unique clusters as control-specific/downregulated-like clusters, but this is a presence/absence cluster comparison rather than a statistical DE test.
- If the cluster strand is `.`, avoid `bedtools intersect -s`; strand-specific matching can return no gene overlaps.
- In the gene-intersection output, gene names are expected in column 10 because `genes.bed` contributes `chr/start/end/gene` after `-wa -wb`.

## 7. Volcano Plot

```r
library(ggplot2)
library(ggrepel)

plot_volcano <- function(res, output_file, title = "your_TREATMENT vs CONTROL", padj_cutoff = 0.05) {
  res_df <- as.data.frame(res)
  res_df$Gene <- rownames(res_df)
  res_df$log10FDR <- -log10(res_df$adj.P.Val)
  res_df$Regulation <- ifelse(
    res_df$adj.P.Val < padj_cutoff & res_df$log2FoldChange > 0,
    "Upregulated",
    ifelse(res_df$adj.P.Val < padj_cutoff & res_df$log2FoldChange < 0, "Downregulated", "Not Significant")
  )

  significance_level <- -log10(padj_cutoff)
  max_abs_logFC <- max(abs(res_df$log2FoldChange), na.rm = TRUE)

  volcano <- ggplot(res_df, aes(x = log2FoldChange, y = log10FDR)) +
    geom_point(aes(color = Regulation), alpha = 0.8) +
    geom_text_repel(
      data = subset(res_df, adj.P.Val < padj_cutoff & abs(log2FoldChange) > 1),
      aes(label = Gene),
      size = 4,
      max.overlaps = 10
    ) +
    geom_hline(yintercept = significance_level, linetype = "dashed", color = "black") +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +
    geom_vline(xintercept = 0, color = "black") +
    xlim(c(-max_abs_logFC, max_abs_logFC)) +
    scale_color_manual(values = c("Upregulated" = "firebrick", "Downregulated" = "royalblue", "Not Significant" = "grey70")) +
    theme_minimal() +
    labs(title = title, x = "Log2 Fold Change", y = "-Log10 FDR", color = "Regulation")

  ggsave(output_file, plot = volcano, dpi = 600, width = 8, height = 6, bg = "white")
  volcano
}

plot_volcano(de$res, "your_volcano_plot.png")
```

## 8. MA Plot

```r
library(ggplot2)
library(ggrepel)

plot_ma <- function(res, output_file, padj_cutoff = 0.05) {
  res_df <- as.data.frame(res)
  res_df$Gene <- rownames(res_df)
  res_df$regulation <- ifelse(
    res_df$adj.P.Val < padj_cutoff & res_df$log2FoldChange > 0,
    "Upregulated",
    ifelse(res_df$adj.P.Val < padj_cutoff & res_df$log2FoldChange < 0, "Downregulated", "Not Significant")
  )

  max_log2fc <- max(abs(res_df$log2FoldChange), na.rm = TRUE)

  ma_plot <- ggplot(res_df, aes(x = log2(baseMean + 1), y = log2FoldChange)) +
    geom_point(aes(color = regulation), alpha = 0.4) +
    geom_text_repel(
      data = subset(res_df, regulation != "Not Significant"),
      aes(label = Gene, color = regulation),
      size = 3,
      max.overlaps = 10,
      show.legend = FALSE
    ) +
    scale_color_manual(values = c("Upregulated" = "firebrick", "Downregulated" = "royalblue", "Not Significant" = "grey70")) +
    theme_minimal() +
    theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1)) +
    labs(title = "MA plot", x = "Log2(Average Expression + 1)", y = "Log2 Fold Change") +
    ylim(-max_log2fc, max_log2fc)

  ggsave(output_file, plot = ma_plot, dpi = 600, width = 8, height = 6, bg = "white")
  ma_plot
}

plot_ma(de$res, "your_MA_plot.png")
```

## 9. Heatmap

```r
library(pheatmap)

plot_de_heatmap <- function(dds, sig_res, output_file, top_n = 20) {
  selected_res <- head(sig_res[order(sig_res$adj.P.Val), ], top_n)
  normalized_counts <- counts(dds, normalized = TRUE)
  selected_counts <- normalized_counts[rownames(selected_res), , drop = FALSE]
  selected_counts <- selected_counts[complete.cases(selected_counts), , drop = FALSE]

  log_selected_counts <- log1p(selected_counts)
  my_palette <- colorRampPalette(c("royalblue", "white", "firebrick"))(100)

  hm <- pheatmap(
    log_selected_counts,
    scale = "row",
    clustering_distance_rows = "euclidean",
    clustering_distance_cols = "euclidean",
    clustering_method = "complete",
    cluster_cols = FALSE,
    show_rownames = TRUE,
    show_colnames = TRUE,
    annotation_col = NULL,
    annotation_legend = FALSE,
    fontsize_row = 10,
    fontsize_col = 10,
    angle_col = 45,
    color = my_palette,
    legend = TRUE
  )

  ggsave(output_file, plot = hm$gtable, dpi = 600, width = 10, height = 8, bg = "white")
  hm
}

plot_de_heatmap(de$dds, de$sig_res, "your_heatmap.png")
```

## Output Files

Common outputs include:

- `your_total_DE_results.csv`
- `your_significant_DE_results.csv`
- `your_volcano_plot.png`
- `your_MA_plot.png`
- `your_heatmap.png`
- `your_DE_results.RData`

## Documentation

Tool documentation PDFs are stored in `docs/tool_documentation/`.

## Citation

Please cite the upstream tools used in your analysis, including SPORTS, tDRMapper, proTRAC, DESeq2, Rsubread, and featureCounts where applicable.

## References

- Shi J, Ko EA, Sanders KM, Chen Q, Zhou T. SPORTS1.0: a tool for annotating and profiling non-coding RNAs optimized for rRNA- and tRNA-derived small RNAs. *Genomics, Proteomics & Bioinformatics*. 2018;16(2):144-151. doi:10.1016/j.gpb.2018.04.004. GitHub: https://github.com/junchaoshi/sports1.1.
- Selitsky SR, Sethupathy P. tDRmapper: challenges and solutions to mapping, naming, and quantifying tRNA-derived RNAs from human small RNA-sequencing data. *BMC Bioinformatics*. 2015;16:354. doi:10.1186/s12859-015-0800-0. GitHub: https://github.com/sararselitsky/tDRmapper.
- Rosenkranz D, Zischler H. proTRAC: a software for probabilistic piRNA cluster detection, visualization and analysis. *BMC Bioinformatics*. 2012;13:5. doi:10.1186/1471-2105-13-5. Official source: https://sourceforge.net/projects/protrac/.
- Love MI, Huber W, Anders S. Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. *Genome Biology*. 2014;15:550. doi:10.1186/s13059-014-0550-8.
- Liao Y, Smyth GK, Shi W. featureCounts: an efficient general purpose program for assigning sequence reads to genomic features. *Bioinformatics*. 2014;30(7):923-930. doi:10.1093/bioinformatics/btt656.
- Liao Y, Smyth GK, Shi W. The R package Rsubread is easier, faster, cheaper and better for alignment and quantification of RNA sequencing reads. *Nucleic Acids Research*. 2019;47(8):e47. doi:10.1093/nar/gkz114.

## Authors

Da Lu, Anthony Hannan Group, Hannan Laboratory, The Florey Institute of Neuroscience and Mental Health.

## License

This repository is released under the MIT License. See `LICENSE`.
