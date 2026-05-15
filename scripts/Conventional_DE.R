setwd("your_directory_here")


#with miRNA, tsRNA, piRNAs, rsRNAs transcripts
library(Rsubread)
buildindex(basename = "good_index", reference = "your_fasta.fa")

#treatment samples
subjunc(index = "good_index", readfile1 = "your_treatment_file.fastqsanger.gz", output_file = "subjuncts_result1.bam")
subjunc(index = "good_index", readfile1 = "your_treatment_file.fastqsanger.gz", output_file = "subjuncts_result2.bam")
subjunc(index = "good_index", readfile1 = "your_treatment_file.fastqsanger.gz", output_file = "subjuncts_result3.bam")
subjunc(index = "good_index", readfile1 = "your_treatment_file.fastqsanger.gz", output_file = "subjuncts_result4.bam")
subjunc(index = "good_index", readfile1 = "your_treatment_file.fastqsanger.gz", output_file = "subjuncts_result5.bam")

# Control samples
subjunc(index = "good_index", readfile1 = "your_control_file.fastqsanger.gz", output_file = "subjuncts_result6.bam")
subjunc(index = "good_index", readfile1 = "your_control_file.fastqsanger.gz", output_file = "subjuncts_result7.bam")
subjunc(index = "good_index", readfile1 = "your_control_file.fastqsanger.gz", output_file = "subjuncts_result8.bam")
subjunc(index = "good_index", readfile1 = "your_control_file.fastqsanger.gz", output_file = "subjuncts_result9.bam")
subjunc(index = "good_index", readfile1 = "your_control_file.fastqsanger.gz", output_file = "subjuncts_result10.bam")




#install and load seqinr package
#install.packages("seqinr")
library(seqinr)
combineanno <- read.fasta("miRNA reference.fa")
head(combineanno)

combine_anno2 <- data.frame(
  names = names(combineanno),
  length = sapply(combineanno, length)
)

#FASTA manipulation
#SETif (!requireNamespace("BiocManager", quietly = TRUE))
#install.packages("BiocManager")

#BiocManager::install("Biostrings")
library(Biostrings)

combineanno1<- readDNAStringSet("your_mirRNA.fa")
combineanno_df <- data.frame(width=width(combineanno1), seq=as.character(combineanno1), names=names(combineanno))
combineanno_df2 <- subset(combineanno_df, select = c("width", "names"))
head(combineanno_df2)
colnames(combineanno_df2)[1] <- "end"
colnames(combineanno_df2)[2] <- "GeneID"
combineanno_df2$Added_Column <- "0"
colnames(combineanno_df2)[3] <- "start"
combineanno_df2$Added_Column <- "*"
colnames(combineanno_df2)[4] <- "strand"
combineanno_df2$Added_Column <- combineanno_df2$GeneID
colnames(combineanno_df2)[5] <- "Chr"
rownames(combineanno_df2) <- NULL
combineanno_df3 <- combineanno_df2[, c(2, 5, 3, 1, 4)]



#*import targets
targets <- read.delim("target.txt")
targets



#*Feature Counts
library("Rsubread")
FC <- featureCounts(files = targets$File_name, annot.ext = combineanno_df3, countMultiMappingReads = FALSE, useMetaFeatures=TRUE)
save(FC, file="yourtreatmentnew_directlymapped.RData")
load("yourtreatmentnew_directlymapped.RData")

##*DESeq2 use raw counts instead of normalized counts so skip normalization step, so use after feature counts created FC
library(DESeq2)

#create the DESeqDataSet object using the modified count data
dds <- DESeqDataSetFromMatrix(countData = FC$counts, colData = targets, design = ~ Treatment)


##Filter Low Count Genes (similar to your rowSums filter)

dds <- dds[rowSums(counts(dds) > 10) > 2, ]
##*Run the DESeq Pipeline:
dds <- DESeq(dds)


##*Retrieve Differential Expression Results for treatment vs. CTRL:
res <- results(dds, contrast=c("Treatment", "your_treatment", "CONTROL"))
##  summary(res)


##if want 0.05 significance level, skip above res
#sig_res <- res[res$padj < 0.05, ]
##*or (if NA occurs)
sig_res <- res[which(res$padj < 0.05 & !is.na(res$padj)), ]
## summary(sig_res)
#*summary
total_genes <- nrow(res)
significant_genes <- nrow(sig_res)
upregulated_genes <- sum(sig_res$log2FoldChange > 0)
downregulated_genes <- sum(sig_res$log2FoldChange < 0)

cat("Total genes analyzed:", total_genes, "\n")
cat("Significant genes (adjusted p-value < 0.05):", significant_genes, "\n")
cat("Upregulated genes:", upregulated_genes, "\n")
cat("Downregulated genes:", downregulated_genes, "\n")



# View the top genes
head(res)
head(sig_res, 10)

# *CHANGE padj to adj.P.Val in the res
colnames(res)[colnames(res) == "padj"] <- "adj.P.Val"
colnames(sig_res)[colnames(sig_res) == "padj"] <- "adj.P.Val"
# *Add a column with gene names
res$Gene <- rownames(res)
sig_res$Gene <- rownames(sig_res)



# *if remove repeates follow both __ and -
clean_gene_name_3 <- function(gene) {
  # Split the gene name by underscores (single and double)
  parts <- unlist(strsplit(gene, "_+"))
  # Keep only unique parts
  unique_parts <- unique(parts)
  # Combine them back together using single underscores
  cleaned_name_3 <- paste(unique_parts, collapse = "_")

  return(cleaned_name_3)
}

# *Apply the function to the Gene column
res$Gene <- sapply(res$Gene, clean_gene_name_3)
sig_res$Gene <- sapply(sig_res$Gene, clean_gene_name_3)

# *Sort the results by adjusted p-value
res_sorted <- res[order(res$adj.P.Val), ]
sig_res_sorted <- sig_res[order(sig_res$adj.P.Val, decreasing = FALSE), ]


# *Save the results to an Excel file
data_miRNA <- as.data.frame(res_sorted)
library(writexl)
write_xlsx(data_miRNA, path = "data_miRNA.xlsx")

data_miRNA_deseq_top <- as.data.frame(sig_res)
library(writexl)
write_xlsx(data_miRNA_deseq_top, path = "data_miRNA_deseq_top.xlsx")

save(dds, file = "dds.RData")
save(res, sig_res, res_sorted, sig_res_sorted, file = "results.RData")
load("dds.RData")
load("results.RData")



##MA
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

## 4.  Save--modify to include frame
ggsave("conventional MA_plot_custom_white1.png", plot = ma_plot,
       dpi = 600, width = 8, height = 6, units = "in")




--------------------------------------------------------------------------------------------------
#** Volcano plot--Assuming you've already run DESeq2 and have a DESeqDataSet object named 'dds_1 and dds_2'
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
  scale_color_manual(values=c("Upregulated"="lightblue", "Downregulated"="yellow", "Not Significant"="grey50"),
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



# Save the plot as a TIFF image with 300 DPI
ggsave("conventional_volcano_plot.tif", plot = volcano, dpi = 300, width = 7.09, height = 4.7, units = "in")

ggsave("conventional volcano_plot_new.png", plot = volcano, dpi = 600, width = 8, height = 6, bg = "white")


--------------------------------------------------------------------------------
##heatmap


  # 1.Define a function to clean up the gene name
  clean_gene_4 <- function(gene) {
    # Check if the gene name starts with "Mus musculus" and remove it
    gene <- sub("^Mus musculus ", "", gene)

    # Split the gene name by underscores (single and double)
    parts <- unlist(strsplit(gene, "_+"))

    # Keep only unique parts to avoid redundancy
    unique_parts <- unique(parts)

    # Combine them back together using single underscores
    cleaned_name_4 <- paste(unique_parts, collapse = "_")

    return(cleaned_name_4)
  }

# Clean up the gene names in the dds object
rownames(dds) <- sapply(rownames(dds), clean_gene_4)

#*
selected_res <-head(sig_res_sorted, 14)


# Clean up the gene names in the selected_res object
rownames(selected_res) <- sapply(rownames(selected_res), clean_gene_4)

# Extract the normalized counts
normalized_counts <- counts(dds, normalized=TRUE)
selected_counts <- normalized_counts[rownames(selected_res), ]
selected_counts <- selected_counts[complete.cases(selected_counts), ]

# 4. Plot Heatmap with desired layout
#install.packages("pheatmap")
library(pheatmap)


# Correct the column names in selected_counts to match desired_sample_order
YourTreatment_sample_order <- c(1:5)
Control_sample_order <- c(1:5)

desired_sample_order <- c(
  paste("YourTreatment", LPS_sample_order, sep=" "),
  paste("CONTROL", Control_sample_order, sep=" ")
)
colnames(selected_counts) <- desired_sample_order

# Define the annotations for the heatmap columns
annotation_col = data.frame(Group = factor(rep(c("CONTROL", "your_Treatment"), each = 5)))

# 3. Plot Heatmap with desired layout
library(pheatmap)

# Log transform for better visualization
log_selected_counts <- log1p(selected_counts)

# Define custom color palette
my_palette <- colorRampPalette(c("yellow", "black", "lightblue"))(100)

# Use pheatmap to plot
hm <- pheatmap(log_selected_counts,
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

ggsave("conventional heatmap600dpi.png", plot = hm$gtable, dpi = 600, width = 10, height = 8, bg = "white")

