rm(list = ls()); gc()
ORIGINAL_DIR <- "."
output <- file.path(ORIGINAL_DIR, "11_GSVA_Temporal")
if (!dir.exists(output)) {dir.create(output, recursive = TRUE)}
setwd(ORIGINAL_DIR)


library(msigdbr)
library(GSVA)
library(pheatmap)
library(RColorBrewer)
library(tidyverse)

load('00_prepared_data/02.Temporal_dat_norm_ComBatSeq.RData')
group<-read.csv('./00_prepared_data/02.Temporal_group.csv')
group$group <- factor(group$group, 
                      levels = c("Breast_early","Breast_mid","Breast_late","Lung_metastasis"))
#genes <- read_csv("09_Venn_Temporal/02.Venn_DEGs.csv")
# table(group$group)
# dat <- dat_norm[rownames(dat_norm) %in% genes$x,]

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

es <- gsva_df[, group$id]


library(dplyr)

all_groups <- levels(group$group)
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

#sel <- order(result_df$Max_Mean_Diff, decreasing = TRUE)[1:20]
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
rownames(annotation_col) <- group$id

annotation_col <- data.frame(Group = group$group1)
rownames(annotation_col) <- group$id

ann_colors <- list(
  Group = c(
    "Breast early" = "#4575B4",  
    "Breast mid"   = "#F9804B",
    "Breast late" = "#69bcce",  
    "Lung metastasis"   = "#D73027" 
  )
)
missing_pathways <- setdiff(pathway_list, rownames(es))
if (length(missing_pathways) > 0) {
  stop("Missing predefined GSVA pathways: ", paste(missing_pathways, collapse = ", "))
}
data <- as.matrix(es[pathway_list, , drop = FALSE])
# rownames(data) <- rownames(data) %>%
#   # 1. Remove the REACTOME_ prefix
#   str_remove("^REACTOME_") %>%
#   # 2. Replace underscores with spaces
#   str_replace_all("_", " ") %>%
#   # 3. Convert to title case
#   str_to_lower() %>%
#   str_to_title()
# The data frame is named data
# Assign the 15 pathway labels

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
pdf(file.path(output, "03.GSVA_heatmap.pdf"),w = 12,h = 6)
pheatmap(
  data_heatmap,             
  scale = "row",         
  annotation_col = annotation_col,
  annotation_colors = ann_colors, 
  show_colnames = FALSE, 
  cluster_cols = FALSE,
  clustering_method = "ward.D2", 
  color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),  
  fontsize_row = 14,  
  row_names_gp = grid::gpar(fontface = "bold", fontsize = 14),
  border_color = NA,    
  main = "GSVA Pathway Enrichment Heatmap"  
)

dev.off()

png(file.path(output, "03.GSVA_heatmap.png"),w = 12,h = 6,units = "in",res = 300)
pheatmap(
  data_heatmap,             
  scale = "row",         
  annotation_col = annotation_col,
  annotation_colors = ann_colors, 
  show_colnames = FALSE, 
  cluster_cols = FALSE,
  clustering_method = "ward.D2", 
  color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),  
  fontsize_row = 14,  
  row_names_gp = grid::gpar(fontface = "bold"),
  border_color = NA,    
  main = "GSVA Pathway Enrichment Heatmap"  
)

dev.off()
