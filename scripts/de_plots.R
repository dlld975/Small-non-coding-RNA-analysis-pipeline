# Shared DE visualization functions: volcano plot, MA plot, and heatmap.

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggrepel)
  library(pheatmap)
})

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
  max_abs_logfc <- max(abs(res_df$log2FoldChange), na.rm = TRUE)

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
    xlim(c(-max_abs_logfc, max_abs_logfc)) +
    scale_color_manual(values = c("Upregulated" = "firebrick", "Downregulated" = "royalblue", "Not Significant" = "grey70")) +
    theme_minimal() +
    labs(title = title, x = "Log2 Fold Change", y = "-Log10 FDR", color = "Regulation")

  ggsave(output_file, plot = volcano, dpi = 600, width = 8, height = 6, bg = "white")
  volcano
}

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
