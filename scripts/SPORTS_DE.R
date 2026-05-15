setwd("your_directory_here")



# sports_DE_from_summary.R
suppressPackageStartupMessages({
  library(data.table)
  library(stringr)
  library(DESeq2)
})

# -----------------------------
# 1)PART A path to summaries
# -----------------------------
summary_dir <- "your_summary_directory_path"
files <- list.files(summary_dir, pattern = "_summary\\.txt$", full.names = TRUE)

sample_names <- str_replace(basename(files), "_summary\\.txt$", "")
sample_table <- data.frame(
  sample = sample_names,
  file = files,
  condition = ifelse(str_detect(sample_names, "^yourtreatment_prefix"), "your_treatment", "Control"),
  stringsAsFactors = FALSE
)


# keep order: Control first, then treatment (or swap if prefer)
sample_table <- sample_table[order(sample_table$condition, sample_table$sample), ]

read_summary <- function(f){
  dt <- data.table::fread(f, sep="\t", header=TRUE, data.table = TRUE)
  setnames(dt, tolower(names(dt)))
  dt[, reads := as.numeric(reads)]
  dt[is.na(sub_class) | sub_class=="", sub_class := "-"]
  dt
}

# Extract Clean_Reads per sample
get_clean_reads <- function(dt){
  dt[tolower(class) == "clean_reads" & sub_class == "-", reads][1]
}

lib_sizes <- sapply(sample_table$file, function(f){
  dt <- read_summary(f)
  get_clean_reads(dt)
})
names(lib_sizes) <- sample_table$sample
lib_sizes
stopifnot(!any(is.na(lib_sizes)), all(lib_sizes > 0))



#PART B Grouping: use class strings from SPORTS
assign_group <- function(class_str){
  x <- tolower(class_str)
  if (str_detect(x, "mirbase.*mirna")) return("miRNA")
  if (str_detect(x, "trna")) return("tsRNA")
  if (str_detect(x, "rrna|yrna|rny|\\b12s\\b|\\b16s\\b|\\b18s\\b|\\b28s\\b|5\\.8s|\\b5s\\b|4\\.5s|45s")) return("rsRNA")
  return("other")
}


#PART C Build all_dt using counts := round(reads)
all_dt <- rbindlist(lapply(seq_len(nrow(sample_table)), function(i){
  s <- sample_table$sample[i]
  dt <- read_summary(sample_table$file[i])

  dt[, sample := s]
  dt[, group := vapply(class, assign_group, character(1))]
  dt[, feature := sub_class]

  # IMPORTANT: use reads directly
  dt[, counts := as.integer(round(reads))]

  dt
}), fill=TRUE)



#PART D Make matrices (filter out the “-” subtotal rows)
make_matrix <- function(group_name){
  dtg <- all_dt[group == group_name]
  dtg <- dtg[feature != "-" & !is.na(feature)]

  dtg <- dtg[, .(counts = sum(counts, na.rm=TRUE)), by=.(feature, sample)]

  mat <- dcast(dtg, feature ~ sample, value.var="counts", fill=0)
  rn <- mat$feature
  mat$feature <- NULL

  m <- as.matrix(mat)
  rownames(m) <- rn
  m <- m[, sample_table$sample, drop=FALSE]
  storage.mode(m) <- "integer"
  m
}

counts_miRNA <- make_matrix("miRNA")
counts_tsRNA <- make_matrix("tsRNA")
counts_rsRNA <- make_matrix("rsRNA")

dim(counts_miRNA); dim(counts_tsRNA); dim(counts_rsRNA)

#Save countmatrix
out_dir <- "your_output_path"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

write.table(
  counts_miRNA,
  file = file.path(out_dir, "counts_miRNA.txt"),
  sep = "\t",
  quote = FALSE,
  col.names = NA
)

write.table(
  counts_tsRNA,
  file = file.path(out_dir, "counts_tsRNA.txt"),
  sep = "\t",
  quote = FALSE,
  col.names = NA
)

write.table(
  counts_rsRNA,
  file = file.path(out_dir, "counts_rsRNA.txt"),
  sep = "\t",
  quote = FALSE,
  col.names = NA
)


#PartE  DE
library(DESeq2)
library(data.table)
library(stringr)
#miRNA
count_data <- read.csv(
  "your_output_path/counts_miRNA.txt",
  sep = "\t",
  header = TRUE,
  row.names = 1,
  check.names = FALSE
)


col_data <- data.frame(
  Treatment = factor(
    ifelse(grepl("^L\\d+_", colnames(count_data)), "your_TREATMENT", "CONTROL"),
    levels = c("CONTROL", "your_TREATMENT")
  ),
  row.names = colnames(count_data)
)

table(col_data$Treatment)

dds <- DESeqDataSetFromMatrix(
  countData = round(count_data),
  colData   = col_data,
  design    = ~ Treatment
)

dds <- dds[rowSums(counts(dds) >= 10) >= 2, ]
dds <- DESeq(dds)

res <- results(dds, contrast = c("Treatment", "your_TREATMENT", "CONTROL"))

# Step 10: Filter significant results (adjusted p-value < 0.05)
sig_res <- res[which(res$padj < 0.05 & !is.na(res$padj)), ]
#sig_res <- res[which(res$padj < 0.10 & !is.na(res$padj)), ]

# Summary of significant results
total_clusters <- nrow(res)
significant_clusters <- nrow(sig_res)
upregulated_clusters <- sum(sig_res$log2FoldChange > 0)
downregulated_clusters <- sum(sig_res$log2FoldChange < 0)

cat("Total piRNA clusters analyzed:", total_clusters, "\n")
cat("Significant clusters (adjusted p-value < 0.05):", significant_clusters, "\n")
cat("Upregulated clusters:", upregulated_clusters, "\n")
cat("Downregulated clusters:", downregulated_clusters, "\n")

# *CHANGE padj to adj.P.Val in the res
colnames(res)[colnames(res) == "padj"] <- "adj.P.Val"
colnames(sig_res)[colnames(sig_res) == "padj"] <- "adj.P.Val"
# *Add a column with gene names
res$Gene <- rownames(res)
sig_res$Gene <- rownames(sig_res)

# Step 11: Save results to CSV
write.csv(as.data.frame(sig_res), "differentially_expressed_miRNA.csv", row.names = TRUE)
write.csv(as.data.frame(res), "Total_expressed_miRNA.csv", row.names = TRUE)





#tsRNA
count_data <- read.csv(
  "your_output_path/counts_tsRNA.txt",
  sep = "\t",
  header = TRUE,
  row.names = 1,
  check.names = FALSE
)


col_data <- data.frame(
  Treatment = factor(
    ifelse(grepl("^L\\d+_", colnames(count_data)), "your_TREATMENT", "CONTROL"),
    levels = c("CONTROL", "your_TREATMENT")
  ),
  row.names = colnames(count_data)
)

dds <- DESeqDataSetFromMatrix(
  countData = round(count_data),
  colData   = col_data,
  design    = ~ Treatment
)

dds <- dds[rowSums(counts(dds) >= 10) >= 2, ]
dds <- DESeq(dds)

res <- results(dds, contrast = c("Treatment", "your_TREATMENT", "CONTROL"))

# Step 10: Filter significant results (adjusted p-value < 0.05)
sig_res <- res[which(res$padj < 0.05 & !is.na(res$padj)), ]
#sig_res <- res[which(res$padj < 0.10 & !is.na(res$padj)), ]

# Summary of significant results
total_clusters <- nrow(res)
significant_clusters <- nrow(sig_res)
upregulated_clusters <- sum(sig_res$log2FoldChange > 0)
downregulated_clusters <- sum(sig_res$log2FoldChange < 0)

cat("Total piRNA clusters analyzed:", total_clusters, "\n")
cat("Significant clusters (adjusted p-value < 0.05):", significant_clusters, "\n")
cat("Upregulated clusters:", upregulated_clusters, "\n")
cat("Downregulated clusters:", downregulated_clusters, "\n")

# *CHANGE padj to adj.P.Val in the res
colnames(res)[colnames(res) == "padj"] <- "adj.P.Val"
colnames(sig_res)[colnames(sig_res) == "padj"] <- "adj.P.Val"
# *Add a column with gene names
res$Gene <- rownames(res)
sig_res$Gene <- rownames(sig_res)

# Step 11: Save results to CSV
write.csv(as.data.frame(sig_res), "differentially_expressed_tsRNA.csv", row.names = TRUE)
write.csv(as.data.frame(res), "Total_expressed_tsRNA.csv", row.names = TRUE)




#rsRNA
count_data <- read.csv(
  "your_output_path/counts_rsRNA.txt",
  sep = "\t",
  header = TRUE,
  row.names = 1,
  check.names = FALSE
)


col_data <- data.frame(
  Treatment = factor(
    ifelse(grepl("^L\\d+_", colnames(count_data)), "your_TREATMENT", "CONTROL"),
    levels = c("CONTROL", "your_TREATMENT")
  ),
  row.names = colnames(count_data)
)

dds <- DESeqDataSetFromMatrix(
  countData = round(count_data),
  colData   = col_data,
  design    = ~ Treatment
)

dds <- dds[rowSums(counts(dds) >= 10) >= 2, ]
dds <- DESeq(dds)

res <- results(dds, contrast = c("Treatment", "your_TREATMENT", "CONTROL"))

# Step 10: Filter significant results (adjusted p-value < 0.05)
sig_res <- res[which(res$padj < 0.05 & !is.na(res$padj)), ]
#sig_res <- res[which(res$padj < 0.10 & !is.na(res$padj)), ]

# Summary of significant results
total_clusters <- nrow(res)
significant_clusters <- nrow(sig_res)
upregulated_clusters <- sum(sig_res$log2FoldChange > 0)
downregulated_clusters <- sum(sig_res$log2FoldChange < 0)

cat("Total piRNA clusters analyzed:", total_clusters, "\n")
cat("Significant clusters (adjusted p-value < 0.05):", significant_clusters, "\n")
cat("Upregulated clusters:", upregulated_clusters, "\n")
cat("Downregulated clusters:", downregulated_clusters, "\n")

# *CHANGE padj to adj.P.Val in the res
colnames(res)[colnames(res) == "padj"] <- "adj.P.Val"
colnames(sig_res)[colnames(sig_res) == "padj"] <- "adj.P.Val"
# *Add a column with gene names
res$Gene <- rownames(res)
sig_res$Gene <- rownames(sig_res)

# Step 11: Save results to CSV
write.csv(as.data.frame(sig_res), "differentially_expressed_rsRNA.csv", row.names = TRUE)
write.csv(as.data.frame(res), "Total_expressed_rsRNA.csv", row.names = TRUE)







#Visualization
-------------------------------------------------------------------------------------
#MA
library(ggplot2)
library(ggrepel)

## 1.  Flag each gene
res$regulation <- with(res,
                       ifelse(adj.P.Val < 0.05 & log2FoldChange  > 0, "Upregulated",
                              ifelse(adj.P.Val < 0.05 & log2FoldChange  < 0, "Downregulated",
                                     "Not Significant"))
)

## keep a fixed order in the legend
res$regulation <- factor(res$regulation,
                         levels = c("Upregulated", "Downregulated", "Not Significant"))

## 2.  Calculate symmetric y‑axis
max_log2fc <- max(abs(res$log2FoldChange), na.rm = TRUE)

## 3.  MA plot
ma_plot <- ggplot(res, aes(x = log2(baseMean + 1), y = log2FoldChange)) +
  geom_point(aes(colour = regulation), alpha = 0.4) +
  scale_colour_manual(values = c("Upregulated"   = "firebrick",
                                 "Downregulated" = "royalblue",
                                 "Not Significant"= "grey70")) +
  geom_text_repel(
    data = subset(res, regulation != "Not Significant"),
    aes(label = Gene, colour = regulation),
    size = 3, box.padding = 0.35, point.padding = 0.5, max.overlaps = 10, show.legend = FALSE
  ) +
  labs(title = "MA plot", x = "Log2(Average Expression + 1)", y = "Log2 FoldChange") +
  theme_minimal() +
  theme(panel.background = element_rect(fill = "white", colour = "white"),
        panel.border     = element_rect(colour = "black", fill = NA, size = 1)) +
  xlim(0, max(log2(res$baseMean + 1))) +
  ylim(-max_log2fc, max_log2fc)

ma_plot




--------------------------------------------------------------------------------------------------
#** Volcano plot
#** Create the volcano plot WC -log10FDR as y axis
library(DESeq2)
library(ggplot2)
library(ggrepel)
# Calculate log10-transformed FDR values
res$log10FDR <- -log10(res$adj.P.Val)
# Convert DESeq Results to a data frame
res_df <- as.data.frame(res)
# Categorize genes based on log2FoldChange
res_df$Regulation <- ifelse(res_df$adj.P.Val < 0.05 & res_df$log2FoldChange > 0, "Upregulated",
                            ifelse(res_df$adj.P.Val < 0.05 & res_df$log2FoldChange < 0, "Downregulated", "Not Significant"))
# Determine the -log10 value for the significance level
significance_level <- -log10(0.05)
# Thresholds for log2FoldChange
logFC_threshold_pos <- 1
logFC_threshold_neg <- -1

#change title font and size
# Calculate the maximum absolute log2 fold change
max_abs_logFC <- max(abs(res_df$log2FoldChange))

nudge_y_amount <- 0.2

annotation_x_pos <- min(res_df$log2FoldChange)+0.5
annotation_y_pos <- significance_level -0.2

# Modify the annotation y-position to make the label visible below the line
volcano <- ggplot(res_df, aes(x=log2FoldChange, y=log10FDR)) +
  geom_point(aes(color=Regulation), alpha=0.8) +
  geom_text_repel(
    aes(label=Gene),
    data = res_df[res_df$adj.P.Val < 0.05 & (res_df$log2FoldChange > logFC_threshold_pos | res_df$log2FoldChange < logFC_threshold_neg), ],
    size = 4.5,
    box.padding = 0.35,
    point.padding = 0.3,
    max.overlaps = 10,
    nudge_y = nudge_y_amount
  ) +
  geom_hline(yintercept = significance_level, linetype="dashed", color="black") +
  geom_vline(xintercept = logFC_threshold_pos, linetype="dashed", color="black") +
  geom_vline(xintercept = logFC_threshold_neg, linetype="dashed", color="black") +
  geom_vline(xintercept = 0, color="black") +
  xlim(c(-max_abs_logFC, max_abs_logFC)) +
  ylim(0, 5) +  #max y value
  scale_color_manual(values=c("Upregulated"="darkseagreen", "Downregulated"="palevioletred2", "Not Significant"="grey50"),
                     breaks=c("Upregulated", "Downregulated", "Not Significant"),
                     labels=c("Upregulated", "Downregulated", "Not Significant")) +
  theme_minimal() +
  labs(title = "your_Treatment vs. Control", subtitle = "",
       x="Log2 Fold Change",
       y="-Log10 FDR",
       color="Regulation") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold"),
    panel.background = element_blank(),  # Removes background
    panel.grid.major = element_blank(),  # Removes major grid lines
    panel.grid.minor = element_blank(),  # Removes minor grid lines
    axis.line = element_line(color="black"),  # Adds black line to enclose the plot area
    legend.justification = c("center", "top"),  # Centers the legend horizontally at the top
    legend.position = "top",  # Moves the legend to the top
    legend.title.align = 0.5,  # Centers the legend title
    legend.text = element_text(size = 12, face = "bold"),  # Make legend text bigger and bold
    legend.title = element_text(size = 12, face = "bold")  # Make legend title bigger and bold
  )




# Display the plot
print(volcano)



--------------------------------------------------------------------------------
##heatmap


#*
selected_res <-head(sig_res,18)

# Extract the normalized counts
normalized_counts <- counts(dds, normalized=TRUE)
selected_counts <- normalized_counts[rownames(selected_res), ]
selected_counts <- selected_counts[complete.cases(selected_counts), ]

# 4. Plot Heatmap with desired layout
#install.packages("pheatmap")
library(pheatmap)


# Correct the column names in selected_counts to match desired_sample_order
Control_sample_order <- c(1:5)
your_treatment_sample_order <- c(1:5)



desired_sample_order <- c(
  paste("your_TREATMENT", your_treatment_sample_order, sep=" "),
  paste("CONTROL", Control_sample_order, sep=" ")
)
colnames(selected_counts) <- desired_sample_order

# Define the annotations for the heatmap columns
annotation_col = data.frame(Group = factor(c(rep("your_TREATMENT", 5), rep("CONTROL", 5))))

# 3. Plot Heatmap with desired layout
library(pheatmap)

# Log transform for better visualization
log_selected_counts <- log1p(selected_counts)

# Define custom color palette
my_palette <- colorRampPalette(c("palevioletred2", "white", "darkseagreen"))(100)

# Use pheatmap to plot
heatmap_plot <- pheatmap(log_selected_counts,
                         scale = "row",
                         clustering_distance_rows = "euclidean",
                         clustering_distance_cols = "euclidean",
                         clustering_method = "complete",
                         cluster_cols = FALSE, # Disable column clustering
                         show_rownames = TRUE,  # Display gene names on the left side
                         show_colnames = TRUE,  # Display column names at the top
                         annotation_col = NULL,
                         annotation_legend = FALSE,
                         fontsize_row = 12, # To reduce the size of row text
                         fontsize_col = 10, # To reduce the size of column text
                         angle_col = 45,  # Rotate column names for better readability
                         color = my_palette,
                         legend = TRUE,  #COLOR BAR
                         width = 10,
                         height = 8)

