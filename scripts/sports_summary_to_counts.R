# Convert SPORTS summary files into miRNA, tsRNA, and rsRNA count matrices.
# Edit the user settings below before running.

suppressPackageStartupMessages({
  library(data.table)
  library(stringr)
})

summary_dir <- "your_sports_summary_directory"
out_dir <- "your_output_directory"
treatment_prefix <- "^your_treatment_prefix"
treatment_label <- "your_TREATMENT"
control_label <- "CONTROL"

files <- list.files(summary_dir, pattern = "_summary\\.txt$", full.names = TRUE)

sample_names <- str_replace(basename(files), "_summary\\.txt$", "")
sample_table <- data.frame(
  sample = sample_names,
  file = files,
  condition = ifelse(str_detect(sample_names, treatment_prefix), treatment_label, control_label),
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

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

write.table(make_matrix("miRNA"), file.path(out_dir, "your_miRNA_counts.txt"), sep = "\t", quote = FALSE, col.names = NA)
write.table(make_matrix("tsRNA"), file.path(out_dir, "your_tsRNA_counts.txt"), sep = "\t", quote = FALSE, col.names = NA)
write.table(make_matrix("rsRNA"), file.path(out_dir, "your_rsRNA_counts.txt"), sep = "\t", quote = FALSE, col.names = NA)
