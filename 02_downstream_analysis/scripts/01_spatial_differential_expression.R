rm(list = ls()); gc()
ORIGINAL_DIR <- "."
output <- file.path(ORIGINAL_DIR, "01_DEGs_Spatial")
if (!dir.exists(output)) {dir.create(output, recursive = TRUE)}
setwd(ORIGINAL_DIR)

library(DESeq2)
library(ggplot2)
library(RColorBrewer)
library(scales)
library(ggrepel)
library(ComplexHeatmap)
library(dplyr)

load('00_prepared_data/01.Spatial_dat_count.RData')
count<-round(dat1)
group<-read.csv('./00_prepared_data/01.Spatial_group.csv')
load('00_prepared_data/01.Spatial_dat_norm_ComBatSeq.RData')
vst_expr <- dat_norm
table(group$group)


# ===================== Fit one model using all samples =====================
all_groups <- c("Breast_normal", "Breast_normal_adjacent", "Breast_tumor",
                "Lung_metastasis", "Lung_normal", "Lung_normal_adjacent")

condition <- factor(group$group)
colData <- data.frame(row.names = colnames(count), condition, sequencing_batch = factor(group$sequencing_batch))

dds <- DESeqDataSetFromMatrix(count, colData, design = ~ sequencing_batch + condition)
keep <- rowSums(counts(dds) >= 1) >= 3
dds <- dds[keep, ]
dds <- DESeq(dds)  


# ===================== Define pairwise comparisons =====================
comparisons <- list(
  c("Breast_tumor", "Breast_normal"),
  c("Breast_tumor","Breast_normal_adjacent"),
  c("Breast_normal_adjacent", "Breast_normal"),
  c("Lung_normal_adjacent","Lung_normal"),
  c("Lung_metastasis","Lung_normal"),
  c("Lung_metastasis","Breast_tumor")
)


root_output <- "01_DEGs_Spatial"
dir.create(root_output, recursive = T, showWarnings = F)

# ===================== 4. Extract contrasts without refitting =====================
for (cmp in comparisons) {
  g1 <- cmp[1]
  g2 <- cmp[2]
  cmp_name <- paste0(g1, "_vs_", g2)
  

  output <- file.path(root_output, cmp_name)
  dir.create(output, recursive = T, showWarnings = F)
  message("\n===== Analyzing: ", cmp_name, " =====")
  
  # ===================== Extract pairwise contrast from the full model =====================
  DEG_all <- results(dds, contrast = c("condition", g1, g2))
  DEG_all <- as.data.frame(DEG_all[order(DEG_all$pvalue), ])
  

  logFC_cutoff <- 1.5
  DEG_all$change <- as.factor(
    ifelse(!is.na(DEG_all$padj) & DEG_all$padj < 0.01 & abs(DEG_all$log2FoldChange) > logFC_cutoff,
           ifelse(DEG_all$log2FoldChange > logFC_cutoff, "UP", "DOWN"), "NOT")
  )
  
  DEGs <- as.data.frame(DEG_all[, 1:7])
  DEGs_diff <- subset(DEGs, change %in% c("UP", "DOWN"))
  

  DEGs <- cbind(GeneSymbol = rownames(DEGs), DEGs)
  write.csv(DEGs, file.path(output, paste0("01.All_", cmp_name, ".csv")), row.names=F, quote=F)
  
  if(nrow(DEGs_diff) > 0){
    DEGs_sig <- cbind(GeneSymbol = rownames(DEGs_diff), DEGs_diff)
    write.csv(DEGs_sig, file.path(output, paste0("02.Sig_", cmp_name, ".csv")), row.names=F, quote=F)
  }
  

    if(nrow(DEGs_diff) > 0){
    # The volcano plot uses DESeq2 padj on the y axis and a fixed x range [-10, 10];
    # padj < 0.01 and |log2FC| > 1.5 define UP and DOWN genes.
    volcano_data <- DEGs
    volcano_data$plot_padj <- ifelse(
      is.na(volcano_data$padj), 1,
      pmax(volcano_data$padj, .Machine$double.xmin)
    )
    volcano_data$change <- factor(
      as.character(volcano_data$change),
      levels = c("DOWN", "NOT", "UP")
    )

    # Only label significant genes that are inside the displayed x-range.
    # Points outside [-10, 10] are intentionally not displayed, and their
    # labels are restricted to the displayed plotting range.
    volcano_xlim <- c(-10, 10)
    label_candidates <- DEGs_diff[
      is.finite(DEGs_diff$log2FoldChange) &
        DEGs_diff$log2FoldChange >= volcano_xlim[1] &
        DEGs_diff$log2FoldChange <= volcano_xlim[2], , drop = FALSE
    ]
    label_data <- rbind(
      head(label_candidates[order(label_candidates$log2FoldChange, decreasing = TRUE), , drop = FALSE], 10),
      head(label_candidates[order(label_candidates$log2FoldChange, decreasing = FALSE), , drop = FALSE], 10)
    )
    label_data <- label_data[!duplicated(rownames(label_data)), , drop = FALSE]
    label_data$plot_padj <- ifelse(
      is.na(label_data$padj), 1,
      pmax(label_data$padj, .Machine$double.xmin)
    )
    label_data$GeneSymbol <- rownames(label_data)

    cmp_name_show <- paste0(gsub("_", " ", g1), " vs ", gsub("_", " ", g2))
    volcano_plot <- ggplot(
      volcano_data,
      aes(x = log2FoldChange, y = -log10(plot_padj), color = change)
    ) +
      geom_point(size = 2.4, alpha = 0.45, na.rm = TRUE) +
      geom_vline(
        xintercept = c(-1.5, 1.5), linetype = 4,
        color = "gray40", linewidth = 0.6
      ) +
      geom_hline(
        yintercept = -log10(0.01), linetype = 4,
        color = "gray40", linewidth = 0.6
      ) +
      geom_label_repel(
        data = label_data,
        aes(label = GeneSymbol),
        max.overlaps = 20, size = 4,
        box.padding = unit(0.5, "lines"),
        min.segment.length = 0,
        point.padding = unit(0.8, "lines"),
        segment.color = "black",
        show.legend = FALSE
      ) +
      scale_color_manual(
        values = c(DOWN = "#1f77b5", NOT = "#7b7c7d", UP = "#db534c"),
        labels = c(DOWN = "DOWN", NOT = "NS", UP = "UP"),
        drop = FALSE
      ) +
      scale_x_continuous(
        breaks = seq(volcano_xlim[1], volcano_xlim[2], by = 5),
        expand = c(0, 0)
      ) +
      coord_cartesian(xlim = volcano_xlim, clip = "on") +
      theme_bw(base_size = 12) +
      theme(
        legend.position = "right",
        panel.grid = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(face = "bold", color = "black", size = 13),
        plot.title = element_text(hjust = 0.5),
        axis.text = element_text(face = "bold", color = "black", size = 15),
        axis.title = element_text(face = "bold", color = "black", size = 15)
      ) +
      labs(
        x = "log2(Fold Change)",
        y = expression(-log[10](adjusted~italic(P)~value)),
        title = cmp_name_show
      )

    base_volcano <- file.path(output, paste0("03.DEG_volcano_", cmp_name))
    ggsave(paste0(base_volcano, ".pdf"), volcano_plot,
           width = 8, height = 6, bg = "white")
    ggsave(paste0(base_volcano, ".png"), volcano_plot,
           width = 8, height = 6, dpi = 600, bg = "white")
    write.csv(
      volcano_data,
      paste0(base_volcano, "_source_data.csv"),
      row.names = FALSE, quote = FALSE
    )
    write.csv(
      label_data,
      paste0(base_volcano, "_labels_within_x10.csv"),
      row.names = FALSE, quote = FALSE
    )
  }

  # ===================== Heatmap =====================

    if(nrow(DEGs_diff) > 0){
    group1 <- group[group$group %in% c(g1,g2), ]
    count1 <- count[, group1$id, drop=F]
    
    # 04 density heatmap uses the complete candidate DEG set under the
    # padj/log2FC thresholds.
    dat_rep2 <- DEGs_diff
    mat <- t(scale(t(count1)))
    mat <- pmin(pmax(mat, -2), 2)
    mat <- mat[rownames(dat_rep2), group1$id, drop=F]
    
    
    group1_anno <- factor(gsub("_", " ", group1$group))
    
   
    group_col <- setNames(
      c("#db534c", "#1f77b5"),
      gsub("_", " ", c(g1, g2))  
    )
    
    pdf(file.path(output, paste0('04.DEG_DensityHeatmap_', cmp_name, '.pdf')), w = 8, h = 10)
    print(
      (if (nrow(mat) >= 2) densityHeatmap(mat, title = "Distribution as heatmap", ylab = " ", height = unit(2.5, "cm")) else Heatmap(mat, name = "expression", show_row_names = FALSE, show_column_names = FALSE, cluster_rows = FALSE, height = unit(12, "cm"), col = colorRampPalette(c('#1f77b5', "white", '#db534c'))(100))) %v%
        HeatmapAnnotation(
          Group = group1_anno,
          col = list(Group = group_col)  
        ) %v%
        Heatmap(
          mat, 
          row_names_gp = gpar(fontsize = 9), 
          show_column_names = F,
          show_row_names = F, 
          name = "expression", 
          height = unit(12, "cm"),
          cluster_rows = T, 
          col = colorRampPalette(c('#1f77b5', "white", '#db534c'))(100)
        )
    )
    dev.off()
    
    png(file.path(output, paste0('04.DEG_DensityHeatmap_', cmp_name, '.png')), w = 8, h = 10, units = 'in', res=600,)
    print(
      (if (nrow(mat) >= 2) densityHeatmap(mat, title = "Distribution as heatmap", ylab = " ", height = unit(2.5, "cm")) else Heatmap(mat, name = "expression", show_row_names = FALSE, show_column_names = FALSE, cluster_rows = FALSE, height = unit(12, "cm"), col = colorRampPalette(c('#1f77b5', "white", '#db534c'))(100))) %v%
        HeatmapAnnotation(
          Group = group1_anno,
          col = list(Group = group_col)  
        ) %v%
        Heatmap(
          mat, 
          row_names_gp = gpar(fontsize = 9), 
          show_column_names = F,
          show_row_names = F, 
          name = "expression", 
          height = unit(12, "cm"),
          cluster_rows = T, 
          col = colorRampPalette(c('#1f77b5', "white", '#db534c'))(100)
        )
    )
    dev.off()
    
    
    
  }
}

message("\n========== All comparisons completed ==========\n")

# ===================== Lung PCA =====================
# This PCA uses the spatial ComBat_seq -> VST matrix and selects the top 500
# genes by across-sample variance within the lung subset. Principal component
# analysis is run with center=TRUE and scale.=FALSE.
pca_output <- file.path(root_output, "05.PCA_Lung_ComBatSeq_VST")
dir.create(pca_output, recursive = TRUE, showWarnings = FALSE)
lung_meta <- group[group$group %in% c("Lung_metastasis", "Lung_normal_adjacent", "Lung_normal"), , drop = FALSE]
lung_meta <- lung_meta[match(colnames(vst_expr), lung_meta$id), , drop = FALSE]
lung_meta <- lung_meta[!is.na(lung_meta$id), , drop = FALSE]
lung_ids <- lung_meta$id
lung_expr <- vst_expr[, lung_ids, drop = FALSE]
lung_expr <- lung_expr[apply(lung_expr, 1, function(z) all(is.finite(z))), , drop = FALSE]
lung_expr <- lung_expr[apply(lung_expr, 1, sd) > 0, , drop = FALSE]
lung_variance <- apply(lung_expr, 1, var)
pca_feature_count <- 500L
lung_top <- order(lung_variance, decreasing = TRUE)[seq_len(min(pca_feature_count, nrow(lung_expr)))]
lung_expr_top <- lung_expr[lung_top, , drop = FALSE]
lung_pca <- prcomp(t(lung_expr_top), center = TRUE, scale. = FALSE)
lung_var_pct <- 100 * lung_pca$sdev^2 / sum(lung_pca$sdev^2)
lung_class <- c(
  Lung_metastasis = "Metastasis",
  Lung_normal = "Normal",
  Lung_normal_adjacent = "Normal adjacent"
)[lung_meta$group]
lung_score <- data.frame(
  sample_id = rownames(lung_pca$x),
  PC1 = lung_pca$x[, 1], PC2 = lung_pca$x[, 2],
  group = factor(unname(lung_class), levels = c("Metastasis", "Normal", "Normal adjacent")),
  sequencing_batch = lung_meta$sequencing_batch,
  stringsAsFactors = FALSE, row.names = NULL
)
lung_palette <- c(Metastasis = "#EF3B3B", Normal = "#2CB34A", `Normal adjacent` = "#1FA9D6")
lung_plot <- ggplot(lung_score, aes(PC1, PC2, colour = group)) +
  geom_point(size = 4.2, alpha = 0.9) +
  scale_colour_manual(values = lung_palette, drop = FALSE, name = NULL) +
  labs(title = "Lung PCA", x = paste0("PC1: ", round(lung_var_pct[1], 1), "% variance"),
       y = paste0("PC2: ", round(lung_var_pct[2], 1), "% variance")) +
  theme_bw(base_size = 11, base_family = "Arial") +
  theme(plot.title = element_text(size = 15, face = "bold", hjust = 0.5),
        axis.title = element_text(size = 12), axis.text = element_text(size = 10, colour = "black"),
        panel.grid.minor = element_blank(), legend.position = "right",
        legend.text = element_text(size = 10), plot.margin = margin(8, 10, 8, 8))
lung_stem <- file.path(pca_output, "PCA_Lung_ComBatSeq_VST")
ggsave(paste0(lung_stem, ".pdf"), lung_plot, width = 183, height = 125,
       units = "mm", device = grDevices::cairo_pdf, family = "Arial")
ggsave(paste0(lung_stem, ".png"), lung_plot, width = 183, height = 125,
       units = "mm", dpi = 600, bg = "white")
write.csv(lung_score, paste0(lung_stem, "_scores.csv"), row.names = FALSE)
write.csv(data.frame(PC = paste0("PC", seq_along(lung_var_pct)), variance_percent = lung_var_pct),
          paste0(lung_stem, "_variance.csv"), row.names = FALSE)
write.csv(data.frame(gene = rownames(lung_expr_top), variance = lung_variance[rownames(lung_expr_top)]),
          paste0(lung_stem, "_top500_genes.csv"), row.names = FALSE)
write.csv(data.frame(figure = "PCA_Lung_ComBatSeq_VST", input = "01.Spatial_dat_norm_ComBatSeq.RData",
                     correction = "existing ComBat_seq -> VST; no second correction",
                     gene_selection = "top 500 genes by VST row variance", center = TRUE,
                     scale = FALSE, feature_count_selection = "feature-count assessment; 500-gene setting",
                     n_samples = ncol(lung_expr_top), n_genes = nrow(lung_expr_top)),
          file.path(pca_output, "PCA_parameters.csv"), row.names = FALSE)
message("Lung PCA written to: ", pca_output)
