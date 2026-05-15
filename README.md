# Small Non-coding RNA Analysis Pipeline

A section-based workflow for small non-coding RNA analysis, including conventional mapping, SPORTS, tDRMapper, proTRAC, cluster differential expression, non-overlapping cluster comparison, and shared R differential expression visualization.

Replace every placeholder such as `yourfile.fastq.gz`, `your_reference.fa`, `your_count_matrix.txt`, `your_output_directory`, `your_TREATMENT`, and `CONTROL` with your own project values.

## Repository Layout

```text
.
|-- README.md
|-- LICENSE
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

dds <- dds[rowSums(counts(dds) > 10) > 2, ]
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

Use this section after SPORTS has generated one `*_summary.txt` file per sample.

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

Use this section for tDRMapper tRNA-derived RNA count matrices.

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

dds <- dds[rowSums(counts(dds) >= 100) >= 2, ]
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

Use this section for proTRAC piRNA cluster count matrices.

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

dds <- dds[rowSums(counts(dds) >= 10) >= 2, ]
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
  min_count = 10,
  min_samples = 2
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

## 6. Non-overlapping Cluster Comparison

Use this section when comparing clusters that are unique to treatment, unique to control, or shared between groups.

```r
library(dplyr)

treatment_clusters <- read.delim("your_treatment_clusters.bed", header = FALSE)
control_clusters <- read.delim("your_control_clusters.bed", header = FALSE)

colnames(treatment_clusters)[1:3] <- c("chr", "start", "end")
colnames(control_clusters)[1:3] <- c("chr", "start", "end")

cluster_key <- function(df) {
  paste(df$chr, df$start, df$end, sep = ":")
}

treatment_clusters$cluster_id <- cluster_key(treatment_clusters)
control_clusters$cluster_id <- cluster_key(control_clusters)

treatment_only <- treatment_clusters[!treatment_clusters$cluster_id %in% control_clusters$cluster_id, ]
control_only <- control_clusters[!control_clusters$cluster_id %in% treatment_clusters$cluster_id, ]
shared_clusters <- treatment_clusters[treatment_clusters$cluster_id %in% control_clusters$cluster_id, ]

write.table(treatment_only, "your_treatment_only_clusters.bed", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(control_only, "your_control_only_clusters.bed", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(shared_clusters, "your_shared_clusters.bed", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
```

For coordinate-aware overlap, use BEDtools before importing into R:

```bash
bedtools intersect -v -a your_treatment_clusters.bed -b your_control_clusters.bed > your_treatment_only_clusters.bed
bedtools intersect -v -a your_control_clusters.bed -b your_treatment_clusters.bed > your_control_only_clusters.bed
bedtools intersect -u -a your_treatment_clusters.bed -b your_control_clusters.bed > your_shared_clusters.bed
```

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

Please cite the upstream tools used in your analysis, including DESeq2, Rsubread, SPORTS, tDRMapper, and proTRAC where applicable.

## Authors

Da Lu, Anthony Hannan Group, Hannan Laboratory, The Florey Institute of Neuroscience and Mental Health.

## License

This repository is released under the MIT License. See `LICENSE`.
