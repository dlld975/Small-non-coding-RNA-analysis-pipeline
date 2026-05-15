# Compare treatment-only, control-only, and shared clusters by exact coordinates.
# For coordinate-aware overlap, prefer the BEDtools commands in README.md.

treatment_clusters_file <- "your_treatment_clusters.bed"
control_clusters_file <- "your_control_clusters.bed"

treatment_clusters <- read.delim(treatment_clusters_file, header = FALSE)
control_clusters <- read.delim(control_clusters_file, header = FALSE)

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
