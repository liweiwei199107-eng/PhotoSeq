rm(list = ls()); gc()
ORIGINAL_DIR <- "."
output <- file.path(ORIGINAL_DIR, "14_GSVA_All")
if (!dir.exists(output)) {dir.create(output, recursive = TRUE)}
setwd(ORIGINAL_DIR)


library(msigdbr)
library(GSVA)
library(RColorBrewer)
library(tidyverse)

load('00_prepared_data/03.All_dat_norm.RData')
group<-read.csv('./00_prepared_data/03.All_group.csv')
group_levels <- c(
  "Breast_early", "Breast_late", "Breast_tumor",
  "Breast_normal_adjacent", "Breast_normal",
  "Lung_metastasis", "Lung_normal_adjacent", "Lung_normal"
)
group$group <- factor(group$group, levels = group_levels)
if (any(is.na(group$group))) stop("Unknown group in 03.All_group.csv.")
group_sizes <- table(group$group)
if (length(group_sizes) != 8 || any(group_sizes != 3)) {
  stop("GSVA-All requires exactly 8 groups with N=3 each; observed: ",
       paste(names(group_sizes), as.integer(group_sizes), collapse = "; "))
}
group <- group[order(group$group, group$sample_id), , drop = FALSE]
write.csv(
  data.frame(group = names(group_sizes), N = as.integer(group_sizes)),
  file.path(output, "00.group_summary.csv"),
  quote = FALSE,
  row.names = FALSE
)
dat_norm <- dat_norm[,group$sample_id]

table(group$group)


msigdbr_species()
msigdbr_collections()
target_species <- "Mus musculus"

options(timeout = 9999999)

target_gene_sets = msigdbr(species = target_species,
                           category = "C2",
                           subcategory = "CP:REACTOME")
table(target_gene_sets$gs_cat)

gene_sets <- target_gene_sets[,c("gs_name","gene_symbol")]
#gene_sets <- target_gene_sets[,c("gs_description","gene_symbol")]
kegg_list <- split(gene_sets$gene_symbol,gene_sets$gs_name)
# The human MSigDB ortholog table does not contain the
# REACTOME_GLYCOLYSIS entry. Add the mouse-native Reactome set so
# the predefined 15-pathway panel remains complete for mouse data.
mouse_reactome <- msigdbr(db_species = "MM", species = "mouse",
                          collection = "M2", subcollection = "CP:REACTOME")
kegg_list[["REACTOME_GLYCOLYSIS"]] <- unique(
  mouse_reactome$gene_symbol[mouse_reactome$gs_name == "REACTOME_GLYCOLYSIS"]
)

params <- gsvaParam(
  exprData = as.matrix(dat_norm),  
  geneSets = kegg_list,         
  kcdf = "Gaussian",        
  maxDiff = TRUE)

gsva_mat <- gsva(
  params,
  verbose = TRUE)

gsva_df <- as.data.frame(gsva_mat) 
write.csv(gsva_df, file.path(output, "01.GSVA_result.csv"), row.names = FALSE)

es <- gsva_df[, group$sample_id]


library(dplyr)

all_groups <- levels(as.factor(group$group))
n_groups <- length(all_groups)

result <- apply(es, 1, function(x) {
  stats <- lapply(all_groups, function(g) {
    vals <- x[group$group == g]
    c(sd = sd(vals), mean = mean(vals))
  })
  
  unlist(stats)
})

result_df <- as.data.frame(t(result))

colnames(result_df) <- paste0(
  rep(c("SD_", "Mean_"), n_groups),
  rep(all_groups, each = 2)
)

result_df$Pathway <- rownames(es)

result_df <- result_df %>% dplyr::select(Pathway, everything())

mean_cols <- grep("Mean_", colnames(result_df), value = TRUE)
result_df$Max_Mean_Diff <- apply(result_df[, mean_cols], 1, function(x) max(x) - min(x))


write.csv(result_df, file.path(output,"02.GSVA_pathway_all_groups_statistics.csv"), row.names = FALSE)

pathway_list <- c(
  "REACTOME_APOPTOTIC_CLEAVAGE_OF_CELL_ADHESION_PROTEINS",
  "REACTOME_CELL_EXTRACELLULAR_MATRIX_INTERACTIONS",
  "REACTOME_G2_M_DNA_REPLICATION_CHECKPOINT",
  "REACTOME_CONDENSATION_OF_PROMETAPHASE_CHROMOSOMES",
  "REACTOME_E2F_ENABLED_INHIBITION_OF_PRE_REPLICATION_COMPLEX_FORMATION",
  "REACTOME_FASL_CD95L_SIGNALING",
  "REACTOME_CLEC7A_INFLAMMASOME_PATHWAY",
  "REACTOME_CD163_MEDIATING_AN_ANTI_INFLAMMATORY_RESPONSE",
  "REACTOME_GLYCOLYSIS",
  "REACTOME_CITRIC_ACID_CYCLE_TCA_CYCLE",
  "REACTOME_FATTY_ACID_METABOLISM",
  "REACTOME_FGFR1_LIGAND_BINDING_AND_ACTIVATION",
  "REACTOME_ERBB2_ACTIVATES_PTK6_SIGNALING",
  "REACTOME_REGULATION_OF_GENE_EXPRESSION_BY_HYPOXIA_INDUCIBLE_FACTOR",
  "REACTOME_VEGF_LIGAND_RECEPTOR_INTERACTIONS"
)

group$group1 <- gsub("_", " ", group$group)
annotation_col <- data.frame(Group = group$group1)
rownames(annotation_col) <- group$sample_id

ann_colors <- list(
  Group = c(
    "Breast early" = "#984EA3",
    "Breast late" = "#1B9E77",
    "Breast tumor" = "#33A02C",
    "Breast normal adjacent" = "#FF7F00",
    "Breast normal" = "#1F78B4",
    "Lung metastasis" = "#E31A1C",
    "Lung normal adjacent" = "#FDBF6F",
    "Lung normal" = "#A6CEE3"
  )
)
missing_pathways <- setdiff(pathway_list, rownames(es))
if (length(missing_pathways) > 0) {
  stop("Missing predefined GSVA pathways: ", paste(missing_pathways, collapse = ", "))
}
data <- as.matrix(es[pathway_list, , drop = FALSE])


rownames(data) <- c(
  "Apoptotic cleavage of cell adhesion proteins",
  "Cell extracellular matrix interactions",
  "G2/M DNA replication checkpoint",
  "Condensation of prometaphase chromosomes",
  "E2F-mediated inhibition of pre-replication complex formation",
  "FASL/CD95L signaling pathway",
  "CLEC7A inflammasome pathway",
  "CD163-mediated anti-inflammatory response",
  "Glycolysis",
  "Citric acid cycle (TCA cycle)",
  "Fatty acid metabolism",
  "FGFR1 ligand binding and activation",
  "ERBB2 activates PTK6 signaling",
  "HIF-dependent regulation of gene expression",
  "VEGF ligand-receptor interactions"
)
data_heatmap <- data
if (any(!is.finite(data_heatmap))) {
  stop("The predefined GSVA pathway panel contains non-finite values.")
}
# Inspect the result
head(rownames(data))


suppressPackageStartupMessages({
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
})

# Use the left-legends layout. The GSVA scores above are not
# batch-corrected here; they use the integrated normalized matrix generated by
# 00_prepare_analysis_data.R, consistent with the GSVA design.
stats_df <- read.csv(
  file.path(output, "02.GSVA_pathway_all_groups_statistics.csv"),
  check.names = FALSE,
  stringsAsFactors = FALSE
)
rownames(gsva_df) <- stats_df$Pathway
if (any(!pathway_list %in% rownames(gsva_df))) {
  stop("The GSVA-All result is missing one or more predefined pathways.")
}
data <- as.matrix(gsva_df[pathway_list, group$sample_id, drop = FALSE])
rownames(data) <- c(
  "Apoptotic cleavage of cell adhesion proteins",
  "Cell extracellular matrix interactions",
  "G2/M DNA replication checkpoint",
  "Condensation of prometaphase chromosomes",
  "E2F-mediated inhibition of pre-replication complex formation",
  "FASL/CD95L signaling pathway",
  "CLEC7A inflammasome pathway",
  "CD163-mediated anti-inflammatory response",
  "Glycolysis",
  "Citric acid cycle (TCA cycle)",
  "Fatty acid metabolism",
  "FGFR1 ligand binding and activation",
  "ERBB2 activates PTK6 signaling",
  "HIF-dependent regulation of gene expression",
  "VEGF ligand-receptor interactions"
)
z <- t(scale(t(data)))
z[!is.finite(z)] <- 0
z <- pmax(pmin(z, 3), -3)

ann_colors <- c(
  "Breast early" = "#984EA3",
  "Breast late" = "#1B9E77",
  "Breast tumor" = "#33A02C",
  "Breast normal adjacent" = "#FF7F00",
  "Breast normal" = "#1F78B4",
  "Lung metastasis" = "#E31A1C",
  "Lung normal adjacent" = "#FDBF6F",
  "Lung normal" = "#A6CEE3"
)
top_ha <- HeatmapAnnotation(
  Group = group$group1,
  col = list(Group = ann_colors),
  show_annotation_name = TRUE,
  annotation_name_gp = gpar(fontsize = 8, fontface = "bold"),
  simple_anno_size = unit(4.5, "mm")
)
col_fun <- colorRamp2(c(-3, 0, 3), c("#2166AC", "white", "#B2182B"))
ht <- Heatmap(
  z,
  name = "Z-score",
  col = col_fun,
  top_annotation = top_ha,
  cluster_columns = FALSE,
  cluster_rows = TRUE,
  clustering_method_rows = "ward.D2",
  show_column_names = FALSE,
  show_row_names = TRUE,
  row_names_side = "right",
  row_names_gp = gpar(fontsize = 9),
  column_title = "GSVA Pathway Enrichment Heatmap",
  column_title_gp = gpar(fontsize = 10, fontface = "bold"),
  show_heatmap_legend = FALSE,
  border = FALSE,
  rect_gp = gpar(col = NA)
)
lgd_heat <- Legend(
  col_fun = col_fun,
  at = c(-3, -2, -1, 0, 1, 2, 3),
  title = "",
  legend_height = unit(32, "mm"),
  labels_gp = gpar(fontsize = 7)
)
lgd_group <- Legend(
  title = "Group",
  legend_gp = gpar(fill = ann_colors),
  labels = names(ann_colors),
  labels_gp = gpar(fontsize = 8),
  title_gp = gpar(fontsize = 8, fontface = "bold")
)
packed_lgd <- packLegend(
  lgd_heat, lgd_group, direction = "horizontal", gap = unit(5, "mm")
)
draw_heatmap <- function() {
  grid.newpage()
  pushViewport(viewport(x = 0.61, y = 0.50, width = 0.74, height = 0.92))
  draw(
    ht,
    show_heatmap_legend = FALSE,
    show_annotation_legend = FALSE,
    newpage = FALSE,
    padding = unit(c(2, 2, 2, 2), "mm")
  )
  popViewport()
  draw(
    packed_lgd,
    x = unit(0.02, "npc"),
    y = unit(0.87, "npc"),
    just = c("left", "top")
  )
}
base <- file.path(output, "03.GSVA_heatmap_left_legends")
cairo_pdf(paste0(base, ".pdf"), width = 13, height = 6,
          family = "sans", bg = "white")
draw_heatmap()
dev.off()
png(paste0(base, ".png"), width = 13, height = 6, units = "in",
    res = 450, type = "cairo", bg = "white")
draw_heatmap()
dev.off()
