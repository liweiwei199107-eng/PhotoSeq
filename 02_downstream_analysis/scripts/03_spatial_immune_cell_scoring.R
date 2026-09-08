rm(list = ls()); gc()
ORIGINAL_DIR <- "."
output <- file.path(ORIGINAL_DIR, "03_Immune_Spatial")
if (!dir.exists(output)) {dir.create(output, recursive = TRUE)}
setwd(ORIGINAL_DIR)


library(tidyverse)
library(GSVA)
library(pheatmap)
library(dplyr)
library(tidyr)
library(ggplot2)

stat_compare_means_local <- function(data) {
  pdat <- data %>%
    group_by(Immune_Cell) %>%
    summarise(
      p = tryCatch(wilcox.test(score ~ Group)$p.value, error = function(e) NA_real_),
      y = max(score, na.rm = TRUE),
      .groups = "drop"
    )
  y_range <- diff(range(data$score, na.rm = TRUE))
  if (!is.finite(y_range) || y_range == 0) y_range <- 1
  pdat$y <- pdat$y + 0.08 * y_range
  pdat$padj <- p.adjust(pdat$p, method = "BH")
  pdat$label <- dplyr::case_when(
    is.na(pdat$padj) ~ "ns",
    pdat$padj < 0.0001 ~ "****",
    pdat$padj < 0.001 ~ "***",
    pdat$padj < 0.01 ~ "**",
    pdat$padj < 0.05 ~ "*",
    TRUE ~ "ns"
  )
  geom_text(data = pdat, aes(x = Immune_Cell, y = y, label = label),
            inherit.aes = FALSE, size = 3)
}

# ======================
# Read inputs
# ======================
load('00_prepared_data/01.Spatial_dat_norm_ComBatSeq.RData')
group <- read.csv('./00_prepared_data/01.Spatial_group.csv')

dat <- as.data.frame(dat_norm, check.names = FALSE)
dat <- dat[!duplicated(rownames(dat)), , drop = FALSE]

# ======================
# ssGSEA
# ======================
gene_set <- read.table("reference/immune/mmc3_ssGSEA_mouse.txt", header = T, sep ="\t", check.names = FALSE)
dat.final <- as.matrix(dat)
gene_list <- split(as.matrix(gene_set)[,1], gene_set[,2])

gsvapar <- gsvaParam(as.matrix(dat.final), gene_list, maxDiff = TRUE, kcdf = "Gaussian")
ssgsea_score <- gsva(gsvapar, verbose = TRUE)
ssgsea_score1 <- as.data.frame(ssgsea_score)
write.csv(ssgsea_score1, file.path(output, "01.ssgsea_result_cell.csv"), quote = F)

# ======================
# Normalize sample IDs
# ======================
normalize_sample_id <- function(x) {
  x %>%
    str_replace_all(c("\\s" = "\\.", "\\(" = "\\.", "\\)" = "\\.", "-" = "\\.", "/" = "\\.")) %>%
    str_replace_all("\\.+", "\\.") %>%
    str_remove("^\\.|\\.$")
}
group$id_normalized <- normalize_sample_id(group$id)
if (!all(group$id_normalized %in% colnames(ssgsea_score))) {
  group$id_normalized <- group$id
}
ssgsea_score <- ssgsea_score[, group$id_normalized, drop = FALSE]

# ======================
# Plot labels without underscores
# ======================
group$Group_plot <- gsub("_", " ", group$group)  # display labels

annotation_col <- data.frame(
  Group = group$Group_plot,  # labels without underscores
  row.names = group$id_normalized
)

# ======================
# Heatmap annotation colors
# ======================
ann_colors <- list(
  Group = c(
    "Breast normal"          = "#4575B4",
    "Breast normal adjacent" = "#F9804B",
    "Breast tumor"           = "#69bccE",
    "Lung metastasis"        = "#D73027",
    "Lung normal"            = "#91BFDB",
    "Lung normal adjacent"   = "#FEE090"
  )
)

# ======================
# Draw heatmap
# ======================
colors <- colorRampPalette(c("#176cdb", "white","#F9804B"))(100)
p <- pheatmap(
  ssgsea_score,
  border_color = NA,
  annotation_col = annotation_col,
  annotation_colors = ann_colors,
  labels_row = NULL,
  scale = "row",
  clustering_method = 'ward.D2',
  show_rownames = T,
  show_colnames = F,
  fontsize_col = 5,
  cluster_cols = F,
  cluster_rows = F
)
png(file.path(output, "01.heatmap.png"), height=7, width=8, units='in', res=600)
print(p)
dev.off()
pdf(file.path(output, "01.heatmap.pdf"), height=7, width=8)
print(p)
dev.off()

# ======================
# Generate boxplots
# ======================
tiics_result <- read.csv(file.path(output, "01.ssgsea_result_cell.csv"), check.names=F, row.names=1)
tiics_result <- as.matrix(tiics_result)
mode(tiics_result) <- "numeric"

comparisons <- list(
  c("Breast_tumor", "Breast_normal_adjacent"),
  c("Lung_metastasis", "Breast_tumor")
)

for (comp in comparisons) {
  group1 <- comp[1]
  group2 <- comp[2]
  
  # Plot labels without underscores
  group1_plot <- gsub("_", " ", group1)
  group2_plot <- gsub("_", " ", group2)
  
  sample1 <- group$id_normalized[group$group == group1]
  sample2 <- group$id_normalized[group$group == group2]
  samples_use <- c(sample1, sample2)
  dat_use <- tiics_result[, samples_use, drop=F]
  
  case_id <- sample1
  ctrl_id <- sample2
  
  pvalue <- numeric(nrow(dat_use))
  for (i in 1:nrow(dat_use)) {
    vec_case <- as.numeric(dat_use[i, case_id])
    vec_ctrl <- as.numeric(dat_use[i, ctrl_id])
    pvalue[i] <- wilcox.test(vec_case, vec_ctrl)$p.value
  }
  
  padj <- p.adjust(pvalue, "fdr")
  rTable <- data.frame(immune_cell = rownames(dat_use), pvalue, padj, row.names=rownames(dat_use))
  rTable$sig <- cut(rTable$padj,
                    breaks=c(-Inf,0.0001,0.001,0.01,0.05,Inf),
                    labels=c("****","***","**","*","ns"))
  
  write.csv(rTable, file.path(output, paste0("02.wilcox_reslut_", group1, "_vs_", group2, ".csv")), quote=F, row.names=F)
  
  diff_Table <- rTable
  plot.cell <- dat_use[rownames(diff_Table), , drop=F]
  
  box_dat <- data.frame(Immune_Cell = rownames(plot.cell), plot.cell) %>%
    gather(key="sample", value="score", -Immune_Cell) %>%
    mutate(Group = case_when(
      sample %in% ctrl_id ~ group2_plot,
      sample %in% case_id ~ group1_plot
    )) %>%
    mutate(Group = factor(Group, levels = c(group2_plot, group1_plot)))
  

  p <- ggplot(box_dat, aes(x=Immune_Cell, y=score, fill=Group)) +
    geom_boxplot(width=0.7, position=position_dodge(0.9), outlier.shape=NA, color="black") +
    scale_fill_manual(values=c("#5BBCD6", "#D51F26")) +
    labs(title = "Cell Infiltration", x = "", y = "Score") +
    stat_compare_means_local(box_dat) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust=0.5, face="bold", size=16),
      axis.text.x = element_text(angle=50, hjust=1, face="bold", size=12),
      axis.text.y = element_text(size=12),
      legend.text = element_text(face="bold", size=11),
      legend.title = element_text(face="bold", size=11),
      panel.grid = element_blank()
    )
  
  ggsave(file.path(output, paste0("02.boxplot_", group1, "_vs_", group2, ".pdf")), p, height=6, width=9)
  png(file.path(output, paste0("02.boxplot_", group1, "_vs_", group2, ".png")), height=6, width=9, units="in", res=300)
  print(p)
  dev.off()
  message(paste("Completed:", group1_plot, " vs ", group2_plot))
}

res3 <- read.csv(file.path(output, "01.ssgsea_result_cell.csv"),check.names = F, row.names = 1)

res3 <- t(res3) %>% as.data.frame()

cor_data <- cor(res3,method="spearman")
cor_pmat_local <- function(x) {
  p <- matrix(NA_real_, nrow = ncol(x), ncol = ncol(x),
              dimnames = list(colnames(x), colnames(x)))
  for (i in seq_len(ncol(x))) {
    for (j in seq_len(ncol(x))) {
      p[i, j] <- tryCatch(
        suppressWarnings(cor.test(x[[i]], x[[j]], method = "spearman", exact = FALSE)$p.value),
        error = function(e) NA_real_
      )
    }
  }
  p
}
bh_adjust_symmetric <- function(pmat) {
  out <- matrix(NA_real_, nrow = nrow(pmat), ncol = ncol(pmat),
                dimnames = dimnames(pmat))
  idx <- upper.tri(pmat, diag = FALSE) & is.finite(pmat)
  out[idx] <- p.adjust(pmat[idx], method = "BH")
  out[lower.tri(out)] <- t(out)[lower.tri(out)]
  diag(out) <- 0
  out
}
corp_raw <- cor_pmat_local(res3)
corp_bh <- bh_adjust_symmetric(corp_raw)

dat.r <- data.frame(cor_data)
dat.r$cell <- rownames(dat.r)
dat.r <- gather(dat.r,'Gene','r',-c(cell))

dat.p <- data.frame(corp_raw)
dat.p$cell <- rownames(dat.p)
dat.p <- gather(dat.p,'Gene','p_raw',-c(cell))
dat.p_bh <- data.frame(corp_bh)
dat.p_bh$cell <- rownames(dat.p_bh)
dat.p_bh <- gather(dat.p_bh,'Gene','p_bh',-c(cell))
corrp <- (cbind(dat.r, dat.p, dat.p_bh))[,c("cell","Gene",'r',"p_raw", "p_bh")]
colnames(corrp) <- c("cell1","cell2",'Correlation',"Pvalue", "Pvalue_BH")
corrp$cell2 <- gsub("\\.", " ",corrp$cell2)
write.csv(corrp,file.path(output,"03.cell_cor.csv"),quote=F)

loadNamespace("corrplot")
pdf(file.path(output,"03.cor_heatmap.pdf"),width = 14,height =14)
col1 <- colorRampPalette(c("#176cdb", "white","#FC8D62"))
corrplot::corrplot(corr = cor_data, 
         p.mat = corp_bh,
         method = "circle",
         type = "upper", 
         tl.pos = "lt", tl.cex = 1.2, tl.col = "black", tl.srt = 50, tl.offset=0.5,
         insig = "label_sig", sig.level = c(.001, .01, .05), pch.cex = 1.2, 
         col = col1(20)) 
corrplot::corrplot(corr = cor_data, 
         type="lower", 
         add=TRUE,
         method="number", 
         tl.pos = "n", 
         cl.pos = "n",
         diag=FALSE,
         number.cex = 1.0,
         col = col1(20))
dev.off()

png(file.path(output,"03.cor_heatmap.png"),height=14,width=14,units='in',res=600)           
col1 <- colorRampPalette(c("#176cdb", "white","#FC8D62"))
corrplot::corrplot(corr = cor_data, 
         p.mat = corp_bh,
         method = "circle",
         type = "upper", 
         tl.pos = "lt", tl.cex = 1.2, tl.col = "black", tl.srt = 50, tl.offset=0.5,
         insig = "label_sig", sig.level = c(.001, .01, .05), pch.cex = 1.2, 
         col = col1(20)) 
corrplot::corrplot(corr = cor_data, 
         type="lower", 
         add=TRUE,
         method="number", 
         tl.pos = "n", 
         cl.pos = "n",
         diag=FALSE,
         number.cex = 1.0,
         col = col1(20))
dev.off()
