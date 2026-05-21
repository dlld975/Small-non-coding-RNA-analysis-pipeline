# proTRAC piRNA cluster count matrix DESeq2 analysis.
# Edit the user settings below before running.

suppressPackageStartupMessages({
  library(DESeq2)
})

working_directory <- "your_protrac_results_directory"
count_file <- "your_piRNA_cluster_count_matrix.txt"
treatment_label <- "your_TREATMENT"
control_label <- "CONTROL"
output_prefix <- "your_proTRAC_cluster"
padj_cutoff <- 0.10
min_count_per_cluster <- your_min_count_threshold
min_samples_with_count <- your_min_sample_threshold

setwd(working_directory)

count_data <- read.csv(
  count_file,
  header = TRUE,
  sep = "\t",
  row.names = 1,
  check.names = FALSE
)

col_data <- data.frame(
  Treatment = factor(c(rep(treatment_label, 5), rep(control_label, 5)), levels = c(control_label, treatment_label)),
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
dds <- dds[rowSums(counts(dds) >= min_count_per_cluster) >= min_samples_with_count, ]
dds <- DESeq(dds)

res <- results(dds, contrast = c("Treatment", treatment_label, control_label))
sig_res <- res[which(res$padj < padj_cutoff & !is.na(res$padj)), ]

colnames(res)[colnames(res) == "padj"] <- "adj.P.Val"
colnames(sig_res)[colnames(sig_res) == "padj"] <- "adj.P.Val"
res$Gene <- rownames(res)
sig_res$Gene <- rownames(sig_res)

write.csv(as.data.frame(res), paste0(output_prefix, "_total_DE_results.csv"), row.names = TRUE)
write.csv(as.data.frame(sig_res), paste0(output_prefix, "_significant_DE_results.csv"), row.names = TRUE)
save(dds, res, sig_res, file = paste0(output_prefix, "_DE_results.RData"))
