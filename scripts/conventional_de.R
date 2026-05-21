# Conventional reference-based small RNA mapping and DESeq2 analysis.
# Edit the user settings below before running.

suppressPackageStartupMessages({
  library(Rsubread)
  library(Biostrings)
  library(DESeq2)
  library(writexl)
})

working_directory <- "your_working_directory"
reference_fasta <- "your_reference.fa"
index_basename <- "your_index_name"
targets_file <- "your_targets_table.txt"
treatment_label <- "your_TREATMENT"
control_label <- "CONTROL"
output_prefix <- "your_conventional"
min_count_per_feature <- your_min_count_threshold
min_samples_with_count <- your_min_sample_threshold

setwd(working_directory)

buildindex(
  basename = index_basename,
  reference = reference_fasta
)

# Optional mapping examples. Add one line per sample or map externally.
# subjunc(index = index_basename, readfile1 = "your_treatment_file_1.fastq.gz", output_file = "your_treatment_1.bam")
# subjunc(index = index_basename, readfile1 = "your_control_file_1.fastq.gz", output_file = "your_control_1.bam")

reference_sequences <- readDNAStringSet(reference_fasta)
annotation <- data.frame(
  GeneID = names(reference_sequences),
  Chr = names(reference_sequences),
  start = 0,
  end = width(reference_sequences),
  strand = "*"
)

targets <- read.delim(targets_file, check.names = FALSE)

fc <- featureCounts(
  files = targets$File_name,
  annot.ext = annotation,
  countMultiMappingReads = FALSE,
  useMetaFeatures = TRUE
)

dds <- DESeqDataSetFromMatrix(
  countData = fc$counts,
  colData = targets,
  design = ~ Treatment
)

# Filter low-count features before DESeq2.
# min_count_per_feature = minimum raw count required in a sample.
# min_samples_with_count = minimum number of samples that must reach that count.
dds <- dds[rowSums(counts(dds) >= min_count_per_feature) >= min_samples_with_count, ]
dds <- DESeq(dds)

res <- results(dds, contrast = c("Treatment", treatment_label, control_label))
sig_res <- res[which(res$padj < 0.05 & !is.na(res$padj)), ]

colnames(res)[colnames(res) == "padj"] <- "adj.P.Val"
colnames(sig_res)[colnames(sig_res) == "padj"] <- "adj.P.Val"
res$Gene <- rownames(res)
sig_res$Gene <- rownames(sig_res)

res_sorted <- res[order(res$adj.P.Val), ]
sig_res_sorted <- sig_res[order(sig_res$adj.P.Val), ]

write_xlsx(as.data.frame(res_sorted), path = paste0(output_prefix, "_total_DE_results.xlsx"))
write_xlsx(as.data.frame(sig_res_sorted), path = paste0(output_prefix, "_significant_DE_results.xlsx"))
save(dds, res, sig_res, res_sorted, sig_res_sorted, file = paste0(output_prefix, "_DE_results.RData"))
