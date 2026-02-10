#tDRmapper
setwd("your_tDRMapper_directory_here")

library(DESeq2)
library(dplyr)
library(ggplot2)

#Load the count data
count_data <- read.csv("tRNA_count_matrix.txt", header = TRUE, sep = "\t", row.names = 1)


# Save the modified data frame to a new CSV file
write.csv(count_data, "tRNA_count_matrix.csv", sep = ",", row.names = TRUE)

#Verify the data structure
head(count_data)

#Set up metadata for samples
col_data <- data.frame(
  Treatment = factor(c(rep("CONTROL", 5), rep("your_TREATMENT", 5))),
  row.names = c("C1", "C2", "C3", "C4", "C5", "your_TREATMENT_1", "your_TREATMENT_2", "your_TREATMENT_3", "your_TREATMENT_4", "your_TREATMENT_5")
)

#Ensure column names in count_data match row names in col_data
stopifnot(all(colnames(count_data) == rownames(col_data)))

#Create DESeq2 dataset
dds <- DESeqDataSetFromMatrix(countData = count_data, colData = col_data, design = ~ Treatment)

#Filter out low-count clusters (e.g., at least 5 counts in 4 samples)
dds <- dds[rowSums(counts(dds) >= 100) >= 2, ]



#Perform differential expression analysis
dds <- DESeq(dds)

#Extract results
res <- results(dds, contrast = c("Treatment", "your_TREATMENT", "CONTROL"))

#Filter significant results (adjusted p-value < 0.05)
sig_res <- res[which(res$padj < 0.05 & !is.na(res$padj)), ]

# Summary of significant results
total_clusters <- nrow(res)
significant_clusters <- nrow(sig_res)
upregulated_clusters <- sum(sig_res$log2FoldChange > 0)
downregulated_clusters <- sum(sig_res$log2FoldChange < 0)

cat("Total tRNA clusters analyzed:", total_clusters, "\n")
cat("Significant clusters (adjusted p-value < 0.05):", significant_clusters, "\n")
cat("Upregulated clusters:", upregulated_clusters, "\n")
cat("Downregulated clusters:", downregulated_clusters, "\n")

# *CHANGE padj to adj.P.Val in the res
colnames(res)[colnames(res) == "padj"] <- "adj.P.Val"
colnames(sig_res)[colnames(sig_res) == "padj"] <- "adj.P.Val"
# *Add a column with gene names
res$Gene <- rownames(res)
sig_res$Gene <- rownames(sig_res)

# Save results to CSV
write.csv(as.data.frame(sig_res), "differentially_expressed_tRNAs.csv", row.names = TRUE)
write.csv(as.data.frame(res), "Total_expressed_tRNAs.csv", row.names = TRUE)


#** Create the volcano plot with -log10FDR as y axis
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
  scale_color_manual(values=c("Upregulated"="red", "Downregulated"="green", "Not Significant"="grey50"),
                     breaks=c("Upregulated", "Downregulated", "Not Significant"),
                     labels=c("Upregulated", "Downregulated", "Not Significant")) +
  theme_minimal() +
  labs(title = "your_TREATMENT vs. Control", subtitle = "",
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


# Save the plot as a TIFF image with 300 DPI
ggsave("volcano_plot.tif", plot = volcano, dpi = 300, width = 7.09, height = 4.7, units = "in")

ggsave("volcano_plot_new.png", plot = volcano, dpi = 600, width = 8, height = 6, bg = "white")


--------------------------------------------------------------------------------
##heatmap for single DE:
# 1. Subset the results to your top DE genes
selected_res <- head(sig_res)

# 2. Extract the normalized counts
normalized_counts <- counts(dds, normalized=TRUE)

# 3. Subset using rownames from 'selected_res', preventing dimension drop
selected_counts <- normalized_counts[rownames(selected_res), , drop = FALSE]
selected_counts <- selected_counts[complete.cases(selected_counts), , drop = FALSE]

# 4. Define your sample order (CONTROL first, then TREATMENT)
your_TREATMENT_sample_order <- c(1:5)
Control_sample_order <- c(1:5)

desired_sample_order <- c(
  paste("your_TREATMENT", your_TREATMENT_sample_order, sep=" ") ,
  paste("CONTROL", Control_sample_order, sep=" ")
)

# 5. Rename columns of 'selected_counts' to match the desired sample order
#    Make sure 'selected_counts' has the same number of columns as 'desired_sample_order' length
colnames(selected_counts) <- desired_sample_order

# 6. Define column annotations so they match the order of 'desired_sample_order'
annotation_col <- data.frame(
  Group = factor(c(rep("your_TREATMENT", 5), rep("CONTROL", 5)))
)
rownames(annotation_col) <- desired_sample_order

# 7. Load pheatmap (install if needed: install.packages("pheatmap"))
library(pheatmap)

# 8.Log-transform the counts for better heatmap contrast
log_selected_counts <- log1p(selected_counts)

# 9. Define a custom color palette (choose one)
my_palette <- colorRampPalette(c( "moccasin", "white", "lightblue1"))(100)

# 10. Plot the heatmap-skip
pheatmap(
  log_selected_counts,
  scale = "row",
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method = "complete",
  cluster_rows = FALSE, #Need to disable row cluster for the single DE
  cluster_cols = FALSE,       # Disable column clustering
  show_rownames = TRUE,       # Show gene names
  show_colnames = TRUE,       # Show sample names
  annotation_col = NULL,      #Disable column bar
  annotation_legend = FALSE,  # Hide the little annotation legend
  fontsize_row = 10,          # Font size for row labels
  fontsize_col = 10,          # Font size for column labels
  angle_col = 45,             # Rotate column labels
  color = my_palette,
  legend = TRUE,              # Show color legend
  width = 10,
  height = 8
)

#use this if want to adjust cell size
pheatmap(
  log_selected_counts,
  scale = "row",
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method = "complete",
  cluster_rows = FALSE,  # Single DE gene -> no row clustering
  cluster_cols = FALSE,  # Disable column clustering
  show_rownames = TRUE,  # Show gene names
  show_colnames = TRUE,  # Show sample names
  annotation_col = NULL, # Disable column annotation bar
  annotation_legend = FALSE,
  # Slightly smaller fonts for a subtle zoom-out effect
  fontsize_row = 9,
  fontsize_col = 9,
  angle_col = 45,
  # Control cell size to shrink it a bit
  cellwidth = 35,
  cellheight = 55,
  color = my_palette,
  legend = TRUE,   # Show color legend
  width = 10,
  height = 8
)

#save high dpi
hm <- pheatmap(
  log_selected_counts,
  scale = "row",
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method = "complete",
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  show_rownames = TRUE,
  show_colnames = TRUE,
  annotation_col = NULL,
  annotation_legend = FALSE,
  fontsize_row = 9,
  fontsize_col = 9,
  angle_col = 45,
  cellwidth = 35,
  cellheight = 55,
  color = my_palette,
  legend = TRUE,
  width = 10,
  height = 8,
)


ggsave("heatmap600dpi.png", plot = hm$gtable, dpi = 600, width = 10, height = 8, bg = "white")


-------------------------------------------------------------------------------------
# Plot MA plot for visualization
library(ggplot2)
library(ggrepel)
# Check for NA or Inf in log2 fold change and filter them out
res_clean <- res[!is.na(res$log2FoldChange) & !is.infinite(res$log2FoldChange), ]


# MA Plot without the smoothing line and with minimal x-axis transformation
# Identify significant points
res$significant <- res$adj.P.Val < 0.05

# Calculate maximum log2 fold change for symmetric y-axis
max_log2fc <- max(abs(res$log2FoldChange), na.rm = TRUE)

# MA Plot with non-overlapping labels for significant points
#with panel border--USE THIS, next one same but define the ma plot to save, overlap Inf for annotation
ggplot(as.data.frame(res), aes(x = log2(baseMean + 1), y = log2FoldChange)) +
  geom_point(aes(color = significant), alpha = 0.4) +
  scale_color_manual(values = c("grey", "orangered")) +
  geom_text_repel(
    data = subset(res, significant),
    aes(label = Gene),
    size = 3,
    box.padding = 0.35,
    point.padding = 0.5,
    max.overlaps = Inf
  ) +
  labs(
    title = "MA Plot Highlighting Significant Points",
    x = "Log2(Average Expression + 1)",
    y = "Log2 Fold Change"
  ) +
  theme_minimal() +
  theme(
    # White panel background
    panel.background = element_rect(fill = "white", colour = "white"),
    # *Add* a black border (fill = NA to keep background white)
    panel.border = element_rect(colour = "black", fill = NA, size = 1)
  ) +
  xlim(0, max(log2(res$baseMean + 1))) +
  ylim(-max_log2fc, max_log2fc)

