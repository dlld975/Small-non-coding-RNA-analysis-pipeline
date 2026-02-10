setwd("your_protrac_results_directory")

library(DESeq2)
library(dplyr)
library(ggplot2)



# Step 1: Load the count data
count_data <- read.csv("piRNA_cluster_count_matrix_renamed.txt", header = TRUE, sep = "\t", row.names = 1)
#count_data <- read.csv("piRNA_cluster_count_matrix_merged_renamed.txt", header = TRUE, sep = "\t", row.names = 1)  #Merged cluster file


# Step 3: Verify the data structure
head(count_data)

# Step 4: Set up metadata for samples
col_data <- data.frame(
  Treatment = factor(c(rep("your_TREATMENT", 5), rep("CONTROL", 5))),
  row.names = c("your_TREATMENT1", "your_TREATMENT2", "your_TREATMENT3", "your_TREATMENT4", "your_TREATMENT5", "C1", "C2", "C3", "C4", "C5")  # Match updated names
)

# Step 5: Ensure column names in count_data match row names in col_data
stopifnot(all(colnames(count_data) == rownames(col_data)))

# Step 6: Create DESeq2 dataset
dds <- DESeqDataSetFromMatrix(countData = count_data, colData = col_data, design = ~ Treatment)

# Step 7: Filter out low-count clusters (e.g., at least 5 counts in 4 samples)
dds <- dds[rowSums(counts(dds) >= 10) >= 2, ]



# Step 8: Perform differential expression analysis
dds <- DESeq(dds)

# Step 9: Extract results
res <- results(dds, contrast = c("Treatment", "your_TREATMENT", "CONTROL"))

# Step 10: Filter significant results (adjusted p-value < 0.05)
sig_res <- res[which(res$padj < 0.05 & !is.na(res$padj)), ]
sig_res <- res[which(res$padj < 0.10 & !is.na(res$padj)), ]

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
write.csv(as.data.frame(sig_res), "differentially_expressed_clusters.csv", row.names = TRUE)
write.csv(as.data.frame(res), "Total_expressed_clusters.csv", row.names = TRUE)

#For merged clusters
write.csv(as.data.frame(sig_res), "differentially_expressed_clusters.csv2", row.names = TRUE)
write.csv(as.data.frame(res), "Total_expressed_clusters2.csv", row.names = TRUE)
#* save a copy of non-filtered res for circulargheatmap
write.csv(as.data.frame(res), "Total_expressed_clusters_for_circularheatmap.csv", row.names = TRUE)

-------------------------------------------------------------------------------------
#MA   FDR0.1
library(ggplot2)
library(ggrepel)

## 1.  Flag each gene
res$regulation <- with(res,
                       ifelse(adj.P.Val < 0.1 & log2FoldChange  > 0, "Upregulated",
                              ifelse(adj.P.Val < 0.1 & log2FoldChange  < 0, "Downregulated",
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

## 4.  Save--modify to include frame
ggsave("your name piRNA cluster FDR 0.1 MA_plot_custom_white1.png", plot = ma_plot,
       dpi = 600, width = 8, height = 6, units = "in")


--------------------------------------------------------------------------------------------------
#** Volcano plot--Assuming you've already run DESeq2 and have a DESeqDataSet object named 'dds_1 and dds_2'
#** Create the volcano plot WC -log10FDR as y axis FDR 0.1
library(DESeq2)
library(ggplot2)
library(ggrepel)

# Calculate log10-transformed FDR values
res$log10FDR <- -log10(res$adj.P.Val)
# Convert DESeq Results to a data frame
res_df <- as.data.frame(res)
# Categorize genes based on log2FoldChange
res_df$Regulation <- ifelse(res_df$adj.P.Val < 0.1 & res_df$log2FoldChange > 0, "Upregulated",
                            ifelse(res_df$adj.P.Val < 0.1 & res_df$log2FoldChange < 0, "Downregulated", "Not Significant"))

# Determine the -log10 value for the significance level
significance_level <- -log10(0.1)
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
    data = res_df[res_df$adj.P.Val < 0.1 & (res_df$log2FoldChange > logFC_threshold_pos | res_df$log2FoldChange < logFC_threshold_neg), ],
    size = 4.5,
    box.padding = 0.35,
    point.padding = 0.3,
    max.overlaps = 20,
    nudge_y = nudge_y_amount
  ) +
  geom_hline(yintercept = significance_level, linetype="dashed", color="black") +
  geom_vline(xintercept = logFC_threshold_pos, linetype="dashed", color="black") +
  geom_vline(xintercept = logFC_threshold_neg, linetype="dashed", color="black") +
  geom_vline(xintercept = 0, color="black") +
  xlim(c(-max_abs_logFC, max_abs_logFC)) +
  ylim(0, 5) +  #max y value
  scale_color_manual(values=c("Upregulated"="#5D6AB4", "Downregulated"="#CFA6C5", "Not Significant"="grey50"),
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
#ggsave("volcano_plotpiRNAcluster.tif", plot = volcano, dpi = 300, width = 7.09, height = 4.7, units = "in")

ggsave("volcano_plot_piRNAclusterFDR0.1.png", plot = volcano, dpi = 600, width = 8, height = 6, bg = "white")


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
your_TREATMENT_sample_order <- c(1:5)



desired_sample_order <- c(
  paste("your_TREATMENT", your_TREATMENT_sample_order, sep=" "),
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
my_palette <- colorRampPalette(c("#cfa6c5", "white", "#a2a3d1"))(100)

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

ggsave("HeatmappiRNAclusterFDR0.1.png", plot = heatmap_plot, dpi = 600, width = 8, height = 6, bg = "white")




