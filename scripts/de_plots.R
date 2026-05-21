# Shared DE visualization functions: volcano plot, MA plot, and heatmap.

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggrepel)
  library(pheatmap)
})

plot_volcano <- function(res,
                         output_file,
                         title = "your_TREATMENT vs CONTROL",
                         adjusted_p_value_cutoff = your_adjusted_p_value_cutoff,
                         label_log2fc_threshold = your_label_log2fc_threshold,
                         label_text_size = your_label_text_size,
                         max_label_overlaps = your_max_label_overlaps,
                         output_dpi = your_output_dpi,
                         output_width = your_output_width,
                         output_height = your_output_height) {
  res_df <- as.data.frame(res)
  res_df$Gene <- rownames(res_df)
  res_df$log10FDR <- -log10(res_df$adj.P.Val)
  res_df$Regulation <- ifelse(
    res_df$adj.P.Val < adjusted_p_value_cutoff & res_df$log2FoldChange > 0,
    "Upregulated",
    ifelse(res_df$adj.P.Val < adjusted_p_value_cutoff & res_df$log2FoldChange < 0, "Downregulated", "Not Significant")
  )

  significance_level <- -log10(adjusted_p_value_cutoff)
  max_abs_logfc <- max(abs(res_df$log2FoldChange), na.rm = TRUE)

  volcano <- ggplot(res_df, aes(x = log2FoldChange, y = log10FDR)) +
    geom_point(aes(color = Regulation), alpha = 0.8) +
    geom_text_repel(
      data = subset(res_df, adj.P.Val < adjusted_p_value_cutoff & abs(log2FoldChange) > label_log2fc_threshold),
      aes(label = Gene),
      size = label_text_size,
      max.overlaps = max_label_overlaps
    ) +
    geom_hline(yintercept = significance_level, linetype = "dashed", color = "black") +
    geom_vline(xintercept = c(-label_log2fc_threshold, label_log2fc_threshold), linetype = "dashed", color = "black") +
    geom_vline(xintercept = 0, color = "black") +
    xlim(c(-max_abs_logfc, max_abs_logfc)) +
    scale_color_manual(values = c("Upregulated" = "firebrick", "Downregulated" = "royalblue", "Not Significant" = "grey70")) +
    theme_minimal() +
    labs(title = title, x = "Log2 Fold Change", y = "-Log10 FDR", color = "Regulation")

  ggsave(output_file, plot = volcano, dpi = output_dpi, width = output_width, height = output_height, bg = "white")
  volcano
}

plot_ma <- function(res,
                    output_file,
                    adjusted_p_value_cutoff = your_adjusted_p_value_cutoff,
                    label_text_size = your_label_text_size,
                    max_label_overlaps = your_max_label_overlaps,
                    output_dpi = your_output_dpi,
                    output_width = your_output_width,
                    output_height = your_output_height) {
  res_df <- as.data.frame(res)
  res_df$Gene <- rownames(res_df)
  res_df$regulation <- ifelse(
    res_df$adj.P.Val < adjusted_p_value_cutoff & res_df$log2FoldChange > 0,
    "Upregulated",
    ifelse(res_df$adj.P.Val < adjusted_p_value_cutoff & res_df$log2FoldChange < 0, "Downregulated", "Not Significant")
  )

  max_log2fc <- max(abs(res_df$log2FoldChange), na.rm = TRUE)

  ma_plot <- ggplot(res_df, aes(x = log2(baseMean + 1), y = log2FoldChange)) +
    geom_point(aes(color = regulation), alpha = 0.4) +
    geom_text_repel(
      data = subset(res_df, regulation != "Not Significant"),
      aes(label = Gene, color = regulation),
      size = label_text_size,
      max.overlaps = max_label_overlaps,
      show.legend = FALSE
    ) +
    scale_color_manual(values = c("Upregulated" = "firebrick", "Downregulated" = "royalblue", "Not Significant" = "grey70")) +
    theme_minimal() +
    theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1)) +
    labs(title = "MA plot", x = "Log2(Average Expression + 1)", y = "Log2 Fold Change") +
    ylim(-max_log2fc, max_log2fc)

  ggsave(output_file, plot = ma_plot, dpi = output_dpi, width = output_width, height = output_height, bg = "white")
  ma_plot
}

plot_de_heatmap <- function(dds,
                            sig_res,
                            output_file,
                            top_n = your_number_of_top_features,
                            row_font_size = your_row_font_size,
                            column_font_size = your_column_font_size,
                            output_dpi = your_output_dpi,
                            output_width = your_output_width,
                            output_height = your_output_height) {
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
    fontsize_row = row_font_size,
    fontsize_col = column_font_size,
    angle_col = 45,
    color = my_palette,
    legend = TRUE
  )

  ggsave(output_file, plot = hm$gtable, dpi = output_dpi, width = output_width, height = output_height, bg = "white")
  hm
}
