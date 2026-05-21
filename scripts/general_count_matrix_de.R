# General DESeq2 differential expression for any small RNA count matrix.
# Edit the user settings below before running.

suppressPackageStartupMessages({
  library(DESeq2)
})

count_file <- "your_count_matrix.txt"
metadata_file <- "your_metadata.csv"
output_prefix <- "your_output_name"
treatment_label <- "your_TREATMENT"
control_label <- "CONTROL"
adjusted_p_value_cutoff <- your_adjusted_p_value_cutoff
min_count <- your_min_count_threshold
min_samples <- your_min_sample_threshold

run_deseq2 <- function(count_file,
                       metadata_file,
                       output_prefix,
                       treatment_label = "your_TREATMENT",
                       control_label = "CONTROL",
                       adjusted_p_value_cutoff = your_adjusted_p_value_cutoff,
                       min_count = your_min_count_threshold,
                       min_samples = your_min_sample_threshold) {
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
  sig_res <- res[which(res$padj < adjusted_p_value_cutoff & !is.na(res$padj)), ]

  colnames(res)[colnames(res) == "padj"] <- "adj.P.Val"
  colnames(sig_res)[colnames(sig_res) == "padj"] <- "adj.P.Val"
  res$Gene <- rownames(res)
  sig_res$Gene <- rownames(sig_res)

  write.csv(as.data.frame(res), paste0(output_prefix, "_total_DE.csv"), row.names = TRUE)
  write.csv(as.data.frame(sig_res), paste0(output_prefix, "_significant_DE.csv"), row.names = TRUE)
  save(dds, res, sig_res, file = paste0(output_prefix, "_DESeq2_results.RData"))

  list(dds = dds, res = res, sig_res = sig_res)
}

de <- run_deseq2(
  count_file = count_file,
  metadata_file = metadata_file,
  output_prefix = output_prefix,
  treatment_label = treatment_label,
  control_label = control_label,
  adjusted_p_value_cutoff = adjusted_p_value_cutoff,
  min_count = min_count,
  min_samples = min_samples
)
