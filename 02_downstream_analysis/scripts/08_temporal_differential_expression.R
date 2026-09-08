rm(list = ls()); gc()
ORIGINAL_DIR <- "."
output <- file.path(ORIGINAL_DIR, "08_DEGs_Temporal")
if (!dir.exists(output)) {dir.create(output, recursive = TRUE)}
setwd(ORIGINAL_DIR)


library(DESeq2)
library(ggplot2)
library(RColorBrewer)
library(scales)
library(ggrepel)
library(ComplexHeatmap)
library(dplyr)


load('00_prepared_data/02.Temporal_dat_count.RData')
count<-round(dat1)
group<-read.csv('./00_prepared_data/02.Temporal_group.csv')
load('00_prepared_data/02.Temporal_dat_norm_ComBatSeq.RData')
vst_expr <- dat_norm
table(group$group)

# ===================== 1. Fit one model using all samples =====================
all_groups <- c("Breast_early", "Breast_late", "Breast_mid", "Lung_metastasis")


condition <- factor(group$group)
colData <- data.frame(row.names = colnames(count), condition, sequencing_batch = factor(group$sequencing_batch))

dds <- DESeqDataSetFromMatrix(count, colData, design = ~ sequencing_batch + condition)
keep <- rowSums(counts(dds) >= 1) >= 3
dds <- dds[keep, ]
dds <- DESeq(dds)  


# ===================== 2. Define pairwise comparisons =====================
comparisons <- list(
  c("Breast_mid", "Breast_early"),
  c("Breast_late", "Breast_mid"),
  c("Lung_metastasis","Breast_late")
)

# ===================== 3. Output directory =====================
root_output <- "08_DEGs_Temporal"
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

# Derived temporal figures use the temporal ComBat_seq -> VST matrix and the
# batch-aware differential-expression model.
derived_root <- file.path(root_output, "05.Derived")
dir.create(derived_root, recursive = TRUE, showWarnings = FALSE)
if (!requireNamespace("uwot", quietly = TRUE)) stop("uwot is required for the temporal gene UMAP")

# a) Gene UMAP coloured by the batch-aware log2FC for Lung metastasis versus
# Breast late. The sign is retained explicitly in the source table. The
# labelled genes are the fixed marker set used in the figure; they
# are not claimed to be the global top genes by log2FC or padj.
# The 7,000 most variable genes are retained.
umap_dir <- file.path(derived_root, "a.gene_UMAP_Breast_late_vs_Lung_metastasis")
dir.create(umap_dir, recursive = TRUE, showWarnings = FALSE)
expr <- vst_expr[, group$id, drop = FALSE]
expr <- expr[apply(expr, 1, function(z) all(is.finite(z))), , drop = FALSE]
gene_var <- apply(expr, 1, var)
umap_feature_count <- 7000L
top_genes <- names(sort(gene_var, decreasing = TRUE))[seq_len(min(umap_feature_count, length(gene_var)))]
umap_input <- expr[top_genes, , drop = FALSE]
set.seed(11)
umap_xy <- uwot::umap(umap_input, n_neighbors = min(85, nrow(umap_input) - 1), min_dist = .82,
                       spread = 1.35, metric = "cosine", n_threads = 1, verbose = FALSE)
umap_de <- results(dds, contrast = c("condition", "Lung_metastasis", "Breast_late"))
umap_dat <- data.frame(UMAP1 = umap_xy[, 1], UMAP2 = umap_xy[, 2], gene = rownames(umap_input),
                       log2FC_BreastLate_vs_LungMet = -as.numeric(umap_de[top_genes, "log2FoldChange"]),
                       padj = as.numeric(umap_de[top_genes, "padj"]), stringsAsFactors = FALSE)
label_genes <- c("Epcam", "Ccl11", "Pecam1", "Acvrl1", "Scgb1a1", "Sftpb", "Egln3", "Ddit4", "Krt8")
label_dat <- umap_dat[umap_dat$gene %in% label_genes, , drop = FALSE]
if (nrow(label_dat) < 5) label_dat <- umap_dat[order(abs(umap_dat$log2FC_BreastLate_vs_LungMet), decreasing = TRUE), , drop = FALSE][seq_len(min(9, nrow(umap_dat))), , drop = FALSE]
umap_plot <- ggplot(umap_dat, aes(UMAP1, UMAP2, colour = log2FC_BreastLate_vs_LungMet)) +
  geom_point(size = 1.05, alpha = .72, na.rm = TRUE) +
  geom_text_repel(data = label_dat, aes(label = gene), colour = "black", size = 4.4, family = "Arial",
                  min.segment.length = 0, box.padding = .65, point.padding = .25, segment.colour = "#333333",
                  segment.linewidth = .35, max.overlaps = Inf, show.legend = FALSE) +
  scale_colour_gradient2(low = "#4C89B8", mid = "#F1E6E0", high = "#B34B5F", midpoint = 0,
                         limits = c(-3.5, 3.5), oob = scales::squish, name = "log2FC") +
  labs(title = "Breast late vs Lung metastasis", x = "UMAP1", y = "UMAP2") +
  theme_classic(base_size = 11, base_family = "Arial") +
  theme(axis.line = element_line(linewidth = .7, colour = "black"), axis.ticks = element_line(linewidth = .65),
        axis.title = element_text(size = 13), axis.text = element_text(size = 10, colour = "black"),
        plot.title = element_text(size = 18, hjust = .5, face = "plain"), panel.grid = element_blank(),
        legend.title = element_text(size = 13), legend.text = element_text(size = 10),
        legend.key.height = grid::unit(1.8, "cm"), plot.margin = margin(8, 12, 8, 8))
ggsave(file.path(umap_dir, "gene_UMAP_log2FC.pdf"), umap_plot, width = 180, height = 145,
       units = "mm", device = grDevices::cairo_pdf, family = "Arial")
ggsave(file.path(umap_dir, "gene_UMAP_log2FC.png"), umap_plot, width = 180, height = 145,
       units = "mm", dpi = 600, bg = "white")
write.csv(umap_dat, file.path(umap_dir, "gene_UMAP_log2FC_source_data.csv"), row.names = FALSE)
write.csv(data.frame(n_genes = nrow(umap_dat), n_top_variable = length(top_genes),
                     n_neighbors = min(85, nrow(umap_input) - 1), min_dist = .82,
                     spread = 1.35, metric = "cosine", contrast = "Lung_metastasis_vs_Breast_late",
                     label_selection = "fixed marker set; fallback to top absolute log2FC only if fewer than 5 markers are present",
                     feature_count_selection = "feature-count assessment; 7,000-gene setting",
                     n_label_genes = nrow(label_dat)),
           file.path(umap_dir, "parameters.csv"), row.names = FALSE)
write.csv(data.frame(label_gene = label_dat$gene,
                     label_selection = "fixed reference marker set",
                     stringsAsFactors = FALSE),
          file.path(umap_dir, "labeled_genes.csv"), row.names = FALSE)

# b) Quadrant scatter: Lung metastasis versus Lung normal and Breast
# late versus Breast early. The latter is an explicit additional DESeq2
# contrast from the same model.
quad_dir <- file.path(derived_root, "b.DEG_quadrant_Breast_late_vs_Lung_metastasis")
dir.create(quad_dir, recursive = TRUE, showWarnings = FALSE)
late_early_de <- results(dds, contrast = c("condition", "Breast_late", "Breast_early"))
lung_normal_file <- file.path("01_DEGs_Spatial", "Lung_metastasis_vs_Lung_normal",
                              "01.All_Lung_metastasis_vs_Lung_normal.csv")
lung_normal <- read.csv(lung_normal_file, stringsAsFactors = FALSE, check.names = FALSE)
rownames(lung_normal) <- lung_normal$GeneSymbol
quad <- merge(
  data.frame(GeneSymbol = rownames(lung_normal), x = lung_normal$log2FoldChange,
             x_padj = lung_normal$padj),
  data.frame(GeneSymbol = rownames(late_early_de), y = late_early_de$log2FoldChange,
             y_padj = late_early_de$padj), by = "GeneSymbol"
)
quad <- quad[is.finite(quad$x) & is.finite(quad$y), , drop = FALSE]
# Match the reference quadrant definition: keep genes that exceed the effect
# threshold on both axes. This deliberately leaves the central cross empty.
quadrant_log2fc_threshold <- 1
quad <- quad[abs(quad$x) >= quadrant_log2fc_threshold &
               abs(quad$y) >= quadrant_log2fc_threshold, , drop = FALSE]
quad$quadrant <- ifelse(quad$x > 0 & quad$y > 0, "Q1: Lung up / Breast late up",
                 ifelse(quad$x < 0 & quad$y > 0, "Q2: Lung down / Breast late up",
                 ifelse(quad$x < 0 & quad$y < 0, "Q3: Lung down / Breast late down", "Q4: Lung up / Breast late down")))
quad_plot <- ggplot(quad, aes(x, y, colour = quadrant)) +
  geom_point(size = 1.8, alpha = .64) +
  geom_vline(xintercept = c(-1, 1), linetype = 2, linewidth = .65, colour = "#9A9A9A") +
  geom_hline(yintercept = c(-1, 1), linetype = 2, linewidth = .65, colour = "#9A9A9A") +
  geom_vline(xintercept = 0, linewidth = .5, colour = "#D0D0D0") + geom_hline(yintercept = 0, linewidth = .5, colour = "#D0D0D0") +
  scale_colour_manual(values = c("Q1: Lung up / Breast late up" = "#D84A3A", "Q2: Lung down / Breast late up" = "#A763B7", "Q3: Lung down / Breast late down" = "#3E88C5", "Q4: Lung up / Breast late down" = "#F28E2B"), guide = "none") +
  annotate("text", x = 6.75, y = 7.7, label = paste0("Q1: Lung up / Breast late up\nn = ", sum(quad$quadrant == "Q1: Lung up / Breast late up")), hjust = 1, vjust = 1, size = 3.8, fontface = "bold", family = "Arial") +
  annotate("text", x = -6.75, y = 7.7, label = paste0("Q2: Lung down / Breast late up\nn = ", sum(quad$quadrant == "Q2: Lung down / Breast late up")), hjust = 0, vjust = 1, size = 3.8, fontface = "bold", family = "Arial") +
  annotate("text", x = -6.75, y = -6.7, label = paste0("Q3: Lung down / Breast late down\nn = ", sum(quad$quadrant == "Q3: Lung down / Breast late down")), hjust = 0, vjust = 0, size = 3.8, fontface = "bold", family = "Arial") +
  annotate("text", x = 6.75, y = -6.7, label = paste0("Q4: Lung up / Breast late down\nn = ", sum(quad$quadrant == "Q4: Lung up / Breast late down")), hjust = 1, vjust = 0, size = 3.8, fontface = "bold", family = "Arial") +
  labs(x = "log2FC (Lung metastasis vs Lung normal)", y = "log2FC (Breast late vs Breast early)") +
  coord_cartesian(xlim = c(-7, 7), ylim = c(-7, 8)) +
  theme_classic(base_size = 11, base_family = "Arial") +
  theme(axis.line = element_line(linewidth = .7), axis.ticks = element_line(linewidth = .65),
        axis.title = element_text(size = 13, face = "bold"), axis.text = element_text(size = 10, colour = "black"),
        panel.grid.major = element_line(colour = "#E8E8E8", linewidth = .45), panel.grid.minor = element_blank(),
        plot.margin = margin(8, 8, 8, 8))
ggsave(file.path(quad_dir, "DEG_quadrant.pdf"), quad_plot, width = 185, height = 165,
       units = "mm", device = grDevices::cairo_pdf, family = "Arial")
ggsave(file.path(quad_dir, "DEG_quadrant.png"), quad_plot, width = 185, height = 165,
       units = "mm", dpi = 600, bg = "white")
write.csv(quad, file.path(quad_dir, "DEG_quadrant_source_data.csv"), row.names = FALSE)
write.csv(data.frame(threshold_log2FC = quadrant_log2fc_threshold,
                      x_contrast = "Lung_metastasis_vs_Lung_normal",
                      y_contrast = "Breast_late_vs_Breast_early", n_genes = nrow(quad)),
          file.path(quad_dir, "parameters.csv"), row.names = FALSE)

message("Temporal derived figures written to: ", derived_root)
