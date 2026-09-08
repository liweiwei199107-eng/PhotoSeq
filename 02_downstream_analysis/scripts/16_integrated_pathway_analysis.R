rm(list = ls()); gc()
ORIGINAL_DIR <- "."
output <- file.path(ORIGINAL_DIR, "16_Integrated_Pathway_Analysis")
if (!dir.exists(output)) {dir.create(output, recursive = TRUE)}
setwd(ORIGINAL_DIR)

library(tidyverse)
library(matrixStats)



load(file.path(ORIGINAL_DIR, "00_prepared_data/03.All_dat_norm.RData"))
group <- read.csv("00_prepared_data/03.All_group.csv", stringsAsFactors = FALSE)
if (!"id" %in% names(group) && "sample_id" %in% names(group)) {
  group$id <- group$sample_id
}
group$group <- gsub("_", " ", group$group)
group$group <- factor(group$group,levels = c("Breast early", "Breast late","Breast tumor",
                                             "Breast normal adjacent","Breast normal",
                                             "Lung metastasis","Lung normal adjacent","Lung normal"
))

# 01.Boxplot intentionally uses raw integrated counts -> DESeq2 VST only.
# It is a normalized-expression distribution plot, not a batch-corrected
# expression plot. Reconstruct the 24-sample deduplicated raw overlap from
# the spatial and temporal raw matrices rather than the ComBat_seq matrix.
load(file.path(ORIGINAL_DIR, "00_prepared_data/01.Spatial_dat_count.RData"))
spatial_raw_count <- round(dat1)
load(file.path(ORIGINAL_DIR, "00_prepared_data/02.Temporal_dat_count.RData"))
temporal_raw_count <- round(dat1)
spatial_ids <- intersect(group$sample_id, colnames(spatial_raw_count))
temporal_ids <- setdiff(intersect(group$sample_id, colnames(temporal_raw_count)), spatial_ids)
raw_common_genes <- intersect(rownames(spatial_raw_count), rownames(temporal_raw_count))
raw_integrated_count <- cbind(
  spatial_raw_count[raw_common_genes, spatial_ids, drop = FALSE],
  temporal_raw_count[raw_common_genes, temporal_ids, drop = FALSE]
)
raw_integrated_count <- raw_integrated_count[, group$sample_id, drop = FALSE]
box_col_data <- data.frame(
  row.names = colnames(raw_integrated_count),
  condition = group$group,
  sequencing_batch = factor(group$sequencing_batch)
)
box_dds <- DESeq2::DESeqDataSetFromMatrix(
  countData = raw_integrated_count,
  colData = box_col_data,
  design = ~ 1
)
box_dds <- box_dds[rowSums(DESeq2::counts(box_dds) >= 1) >= 3, ]
box_dds <- DESeq2::DESeq(box_dds, quiet = TRUE)
box_vsd <- DESeq2::vst(box_dds, blind = TRUE)
box_expr <- SummarizedExperiment::assay(box_vsd)

group_means <- sapply(levels(group$group), function(g) {
  sample_idx <- group$sample_id[group$group == g]
  rowMeans(box_expr[, sample_idx, drop = FALSE])
})

group_mat <- as.matrix(group_means)
colnames(group_mat) <- c(
  "Breast early", "Breast late", "Breast mid(tumor)",
  "Breast mid(normal_adjacent)", "Breast mid(normal)",
  "Lung metastasis", "Lung normal_adjacent", "Lung normal"
)
colors <- c(
 "#984EA3","#1B9E77",  "#4DAF4A",  "#FF7F00", "#1F78B4","#A6CEE3","#FDBF6F","#E31A1C"
)

pdf(file.path(output,"01.Boxplot.pdf"), width =8, height = 6)
par(mar = c(10, 8, 6, 8)) 
boxplot(
  group_mat,
  col = colors,
  outcol =colors,  
  outpch = 19,
  outcex = 0.6,
  las = 2,              
  cex.axis =0.9,
  outline = TRUE,
  ylab = "Gene-level mean VST-normalized expression",
  xaxt = "n"        
)

text(
  x = 1:ncol(group_mat),
  y = par("usr")[3] - 0.08 * diff(par("usr")[3:4]),
  labels = colnames(group_mat),
  srt = 45,   
  adj = 1,    
  xpd = TRUE,
  cex = 0.9   
)
dev.off()

png(file.path(output,"01.Boxplot.png"), width =8, height = 6, units = 'in',res = 600)
par(mar = c(10, 8, 6, 8)) 
boxplot(
  group_mat,
  col = colors,
  outcol =colors,  
  outpch = 19,
  outcex = 0.6,
  las = 2,              
  cex.axis =0.9,
  outline = TRUE,
  ylab = "Gene-level mean VST-normalized expression",
  xaxt = "n"        
)

text(
  x = 1:ncol(group_mat),
  y = par("usr")[3] - 0.08 * diff(par("usr")[3:4]),
  labels = colnames(group_mat),
  srt = 45,   
  adj = 1,    
  xpd = TRUE,
  cex = 0.9   
)
dev.off()





KEGG_Spatial <- read.csv("05_Enrich_Spatial/04.KEGG_res.csv", stringsAsFactors = FALSE)
KEGG_Temporal <- read.csv("10_Enrich_Temporal/04.KEGG_res.csv", stringsAsFactors = FALSE)
KEGG_Spatial <- KEGG_Spatial[!is.na(KEGG_Spatial$p.adjust) & KEGG_Spatial$p.adjust < 0.05, , drop = FALSE]
KEGG_Temporal <- KEGG_Temporal[!is.na(KEGG_Temporal$p.adjust) & KEGG_Temporal$p.adjust < 0.05, , drop = FALSE]

both <- intersect(KEGG_Spatial$Description,KEGG_Temporal$Description)
KEGG_both <-KEGG_Spatial[KEGG_Spatial$Description %in% both,]


Spatial <- setdiff(KEGG_Spatial$Description, both)
KEGG_Spatial1 <-KEGG_Spatial[KEGG_Spatial$Description %in% Spatial,]

Temporal <- setdiff(KEGG_Temporal$Description, both)
KEGG_Temporal1 <-KEGG_Temporal[KEGG_Temporal$Description %in% Temporal,]


KEGG_data1 <- KEGG_Spatial1 %>%
  dplyr::slice_head(n = 6) %>%
  dplyr::mutate(
    Category = "Spatial heterogeneity",
    neg_log10_pvalue = -log10(p.adjust),
    geneID = geneID
    )
KEGG_data2 <- KEGG_Temporal1 %>%
  dplyr::slice_head(n = 6) %>%
  dplyr::mutate(
    Category = "Temporal heterogeneity",
    neg_log10_pvalue = -log10(p.adjust),
    geneID = geneID
  )

KEGG_data3 <- KEGG_both %>%
  dplyr::slice_head(n = 6) %>%
  dplyr::mutate(
    Category = "Both",
    neg_log10_pvalue = -log10(p.adjust),
    geneID = geneID
  )

library(dplyr)
library(ggplot2)
library(patchwork)

all_result <- rbind(KEGG_data1, KEGG_data2, KEGG_data3) %>% 
  mutate(
    Description = sub(" - Mus musculus \\(house mouse\\)$", "", Description))
all_result <- all_result %>%
  mutate(
    Category = factor(Category, levels = c("Temporal heterogeneity", "Spatial heterogeneity", "Both"))
  ) %>%
  arrange(Category, p.adjust)

all_result$Description <- factor(all_result$Description, levels = unique(all_result$Description))

total_pathways <- nrow(all_result)
counts_per_category <- all_result %>% dplyr::count(Category) %>% pull(n)
p1 <- ggplot(all_result, aes(x = neg_log10_pvalue, y = Description, fill = Category)) +
  geom_bar(stat = "identity", width = 0.65) +
  geom_text(aes(x = 0.1, label = Description),
            fontface = "bold", size = 4, hjust = 0) +
  theme_classic() +
  theme(
    plot.title = element_text(size = 16, hjust = 0.5),
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.y = element_line(colour = '#7b7c7d', linewidth = 0.5),
    axis.line.x = element_line(colour = '#7b7c7d', linewidth = 0.5),
    axis.text.x = element_text(colour = '#7b7c7d', size = 10),
    axis.ticks.x = element_line(colour = '#7b7c7d'),
    axis.title.x = element_text(colour = 'black', size = 12),
    legend.position = "none"
  ) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_fill_manual(values = c('#037498', "#FDBF6F", '#db534c')) +
  scale_y_discrete(expand = c(0.05, 0)) +
  ggtitle("KEGG Enrichment")

all_result$Category <- factor(
  all_result$Category,
  levels = c("Temporal heterogeneity", "Spatial heterogeneity", "Both")
)

counts_per_category <- all_result %>% dplyr::count(Category) %>% pull(n)
df_anno <- all_result %>% 
  dplyr::count(Category) %>%
  dplyr::mutate(
    ymax = cumsum(n),
    ymin = lag(ymax, default = 0),
    position = (ymin + ymax) / 2,
    Category = factor(Category, levels = levels(all_result$Category))
  )

p2 <- ggplot() +
  geom_rect(data = df_anno,
            aes(xmin = 0, xmax = 1,
                ymin = ymin,
                ymax = ymax,
                fill = Category),
            alpha = 1) +
  geom_text(data = df_anno,
            aes(x = 0.5, y = position, label = Category),
            size = 4,
            fontface = 'bold.italic',
            angle = 90) +
  theme_classic() +
  theme(
    axis.text = element_blank(),
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    axis.line = element_blank(),
    legend.position = "none"
  ) +
  scale_fill_manual(values = c('#037498', "#FDBF6F", '#db534c')) +
  scale_y_continuous(limits = c(0, nrow(all_result)), expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0))

combined_plot <- p2 + p1 + plot_layout(widths = c(0.05, 1))
print(combined_plot)

ggsave(file.path(output, "02.KEGG_plot.pdf"), combined_plot, width = 10, height =8)
ggsave(file.path(output, "02.KEGG_plot.png"), combined_plot, width = 10, height = 8)




library(dplyr)
library(org.Mm.eg.db)  

# Use a local snapshot of the official KEGG REST tables so this step is
# reproducible and does not fail when clusterProfiler's HTTPS request fails.
kegg_link_file <- file.path(ORIGINAL_DIR, "reference/kegg/kegg_link_mmu_pathway.txt")
kegg_name_file <- file.path(ORIGINAL_DIR, "reference/kegg/kegg_list_pathway_mmu.txt")
if (!file.exists(kegg_link_file) || !file.exists(kegg_name_file)) {
  stop("Local KEGG annotation files are missing.")
}

path2gene <- read.delim(
  kegg_link_file,
  header = FALSE,
  col.names = c("from", "to"),
  stringsAsFactors = FALSE
)
path2gene$from <- sub("^path:", "", path2gene$from)
path2gene$to <- sub("^mmu:", "", path2gene$to)

path2name <- read.delim(
  kegg_name_file,
  header = FALSE,
  col.names = c("from", "to"),
  stringsAsFactors = FALSE
)
path2name$from <- sub("^path:", "", path2name$from)
# Keep the original KEGG description text so it matches the enrichment tables.

kegg_anno <- path2gene %>%
  dplyr::rename(pathway_id = from, gene_id = to) %>%
  inner_join(
    path2name %>% dplyr::rename(pathway_id = from, pathway_name = to),
    by = "pathway_id"
  ) %>%
  dplyr::select(pathway_id, pathway_name, gene_id)


symbol_map <- mapIds(
  x = org.Mm.eg.db,
  keys = as.character(kegg_anno$gene_id),  # use character identifiers
  keytype = "ENTREZID",
  column = "SYMBOL",
  multiVals = "first"  # use the first symbol when multiple mappings exist
)

kegg_anno$gene_symbol <- symbol_map[as.character(kegg_anno$gene_id)]
kegg_anno <- kegg_anno[!is.na(kegg_anno$gene_symbol) & kegg_anno$gene_symbol != "", , drop = FALSE]
head(kegg_anno)



library(GSVA)
library(dplyr)
library(tidyr)

# Select the spatial profiles by group name from the integrated matrix,
# avoiding duplicate Breast mid and Lung metastasis profiles.
group1 <- group[!group$group %in% c("Breast early", "Breast late"), , drop = FALSE]
dat <- dat_norm[,group1$sample_id]

all <- rbind(KEGG_Spatial, KEGG_Temporal) 
all_KEGG <- kegg_anno[kegg_anno$pathway_name %in% all$Description,]



pathway_list <- split(
  all_KEGG$gene_symbol, 
  all_KEGG$pathway_name
)

pathway_list <- pathway_list[lengths(pathway_list) >= 10]



message("Running ssGSEA...")
ssgsea_params <- ssgseaParam(
  exprData = as.matrix(dat),
  geneSets = pathway_list,
  minSize = 10
)

gsva_mat <- gsva(ssgsea_params, verbose = FALSE)

rownames(gsva_mat) <- sub(" - Mus musculus.*$", "", rownames(gsva_mat))
score_df <- as.data.frame(t(gsva_mat))

write.csv(score_df, file.path(output, '03.ssgsea_score.csv'), row.names = T)

score_df$id <- rownames(score_df)
score_df <- merge(score_df, group, by = "id")
score_df <- score_df %>%
  mutate(
    Tissue = case_when(
      grepl("^Lung", group) ~ "Lung",
      grepl("^Breast", group) ~ "Breast"
    ),
    Region = case_when(
      grepl("tumor", group, ignore.case = TRUE) ~ "Tumor",
      grepl("normal adjacent", group, ignore.case = TRUE) ~ "Normal_Adjacent",
      grepl("normal", group, ignore.case = TRUE) & !grepl("adjacent", group) ~ "Normal",
      grepl("metastasis", group, ignore.case = TRUE) ~ "Tumor"
    ),
    Tissue = factor(Tissue, levels = c("Breast", "Lung")),
    Region = factor(Region, levels = c("Normal", "Normal_Adjacent", "Tumor")),
    # N1/N2/N3 are biological samples within each tissue. The three regions
    # from one sample are repeated observations, not independent replicates.
    subject_id = paste(Tissue, replicate, sep = "__")
  )



pathway_names <- rownames(gsva_mat)
region_levels <- c("Normal", "Normal_Adjacent", "Tumor")
tissue_levels <- c("Breast", "Lung")

if (!requireNamespace("nlme", quietly = TRUE)) {
  stop("The repeated-measures 16_Integrated_Pathway_Analysis/05.boxplot analysis requires the nlme package.")
}

model_matrix <- function(tissue, region) {
  data.frame(
    Tissue = factor(tissue, levels = tissue_levels),
    Region = factor(region, levels = region_levels)
  )
}

# Estimate the Lung - Breast contrast from a repeated-measures model. When
# region is NULL, the contrast is averaged equally over the three regions.
model_contrast <- function(fit, region = NULL) {
  if (is.null(region)) {
    breast_x <- colMeans(model.matrix(
      ~ Tissue * Region,
      model_matrix(rep("Breast", length(region_levels)), region_levels)
    ))
    lung_x <- colMeans(model.matrix(
      ~ Tissue * Region,
      model_matrix(rep("Lung", length(region_levels)), region_levels)
    ))
  } else {
    breast_matrix <- model.matrix(~ Tissue * Region,
                                  model_matrix("Breast", region))
    lung_matrix <- model.matrix(~ Tissue * Region,
                                model_matrix("Lung", region))
    breast_x <- as.numeric(breast_matrix)
    lung_x <- as.numeric(lung_matrix)
    names(breast_x) <- colnames(breast_matrix)
    names(lung_x) <- colnames(lung_matrix)
  }

  beta <- nlme::fixef(fit)
  vbeta <- as.matrix(stats::vcov(fit))
  contrast <- numeric(length(beta))
  names(contrast) <- names(beta)
  common_names <- intersect(names(beta), colnames(vbeta))
  contrast[common_names] <- lung_x[common_names] - breast_x[common_names]
  estimate <- sum(contrast * beta)
  std_error <- sqrt(as.numeric(t(contrast) %*% vbeta %*% contrast))
  df <- suppressWarnings(min(as.numeric(fit$fixDF$X), na.rm = TRUE))
  if (!is.finite(df) || df <= 0) df <- max(1, nrow(fit$data) - length(beta))
  statistic <- estimate / std_error
  pvalue <- if (is.finite(statistic)) {
    2 * stats::pt(-abs(statistic), df = df)
  } else {
    NA_real_
  }
  data.frame(estimate = estimate, std_error = std_error,
             statistic = statistic, df = df, pvalue = pvalue,
             stringsAsFactors = FALSE)
}

fit_repeated_pathway <- function(temp) {
  fit <- tryCatch(
    nlme::lme(
      fixed = score ~ Tissue * Region,
      random = ~ 1 | subject_id,
      data = temp,
      method = "REML",
      na.action = na.omit,
      control = nlme::lmeControl(returnObject = TRUE)
    ),
    error = function(e) NULL
  )
  if (is.null(fit)) {
    return(list(fit = NULL, interaction_pvalue = NA_real_,
                contrasts = data.frame()))
  }

  an <- tryCatch(anova(fit), error = function(e) NULL)
  interaction_pvalue <- NA_real_
  if (!is.null(an) && "p-value" %in% colnames(an)) {
    interaction_row <- which(rownames(an) == "Tissue:Region")
    if (length(interaction_row) == 1L) {
      interaction_pvalue <- as.numeric(an[interaction_row, "p-value"])
    }
  }

  contrast_list <- list(
    Combined = model_contrast(fit, NULL),
    Normal = model_contrast(fit, "Normal"),
    `Normal adjacent` = model_contrast(fit, "Normal_Adjacent"),
    Tumor = model_contrast(fit, "Tumor")
  )
  contrast_df <- bind_rows(contrast_list, .id = "Compare")
  list(fit = fit, interaction_pvalue = interaction_pvalue,
       contrasts = contrast_df)
}

diff_result_list <- list()
model_test_list <- list()
interaction_result_list <- list()

for (p in pathway_names) {
  message("Calculating pathway: ", p)

  temp <- score_df[, c("id", "Tissue", "Region", "subject_id", p)]
  colnames(temp)[5] <- "score"
  temp <- temp[!is.na(temp$Tissue) & !is.na(temp$Region) &
                 !is.na(temp$score), , drop = FALSE]

  fit_result <- fit_repeated_pathway(temp)
  contrast_df <- fit_result$contrasts
  if (nrow(contrast_df) > 0) {
    contrast_df$pathway <- p
    model_test_list[[length(model_test_list) + 1L]] <- contrast_df
  }
  interaction_result_list[[length(interaction_result_list) + 1L]] <- data.frame(
    pathway = p,
    interaction_pvalue = fit_result$interaction_pvalue,
    stringsAsFactors = FALSE
  )

  region_means <- temp %>%
    group_by(Tissue, Region) %>%
    summarise(mean_score = mean(score, na.rm = TRUE), .groups = "drop")
  subject_combined <- temp %>%
    group_by(Tissue, subject_id) %>%
    summarise(score = mean(score, na.rm = TRUE), .groups = "drop") %>%
    group_by(Tissue) %>%
    summarise(combined_score = mean(score, na.rm = TRUE), .groups = "drop")
  get_region_mean <- function(tissue, region) {
    region_means$mean_score[region_means$Tissue == tissue &
                              region_means$Region == region][1]
  }
  get_combined_mean <- function(tissue) {
    subject_combined$combined_score[subject_combined$Tissue == tissue][1]
  }
  breast_normal <- get_region_mean("Breast", "Normal")
  breast_adjacent <- get_region_mean("Breast", "Normal_Adjacent")
  breast_tumor <- get_region_mean("Breast", "Tumor")
  lung_normal <- get_region_mean("Lung", "Normal")
  lung_adjacent <- get_region_mean("Lung", "Normal_Adjacent")
  lung_tumor <- get_region_mean("Lung", "Tumor")
  breast_all <- get_combined_mean("Breast")
  lung_all <- get_combined_mean("Lung")

  diff_result_list[[p]] <- data.frame(
    Breast_Normal = breast_normal,
    Breast_Normal_Adjacent = breast_adjacent,
    Breast_Tumor = breast_tumor,
    Lung_Normal = lung_normal,
    Lung_Normal_Adjacent = lung_adjacent,
    Lung_Tumor = lung_tumor,
    Lung_All = lung_all,
    Breast_All = breast_all,
    diff_combined = lung_all - breast_all,
    diff_normal = lung_normal - breast_normal,
    diff_adjacent = lung_adjacent - breast_adjacent,
    diff_tumor = lung_tumor - breast_tumor,
    abs_combined = abs(lung_all - breast_all),
    max_abs_region = max(abs(lung_normal - breast_normal),
                          abs(lung_adjacent - breast_adjacent),
                          abs(lung_tumor - breast_tumor)),
    pathway = p,
    interaction_pvalue = fit_result$interaction_pvalue,
    stringsAsFactors = FALSE
  )
}

diff_all <- bind_rows(diff_result_list)
model_tests <- bind_rows(model_test_list)
model_tests$padj <- p.adjust(model_tests$pvalue, method = "BH")
model_tests$sig <- dplyr::case_when(
  is.na(model_tests$padj) ~ "ns",
  model_tests$padj < 0.0001 ~ "****",
  model_tests$padj < 0.001 ~ "***",
  model_tests$padj < 0.01 ~ "**",
  model_tests$padj < 0.05 ~ "*",
  TRUE ~ "ns"
)
model_tests$method <- "nlme::lme; random intercept for subject_id"

interaction_tests <- bind_rows(interaction_result_list)
interaction_tests$interaction_padj <- p.adjust(
  interaction_tests$interaction_pvalue, method = "BH"
)
interaction_tests$interaction_sig <- ifelse(
  !is.na(interaction_tests$interaction_padj) &
    interaction_tests$interaction_padj < 0.05, "*", "ns"
)
diff_all <- diff_all %>%
  left_join(interaction_tests, by = c("pathway", "interaction_pvalue"))

# Effect-size screening is exploratory only. Formal selection requires a
# BH-significant Tissue x Region interaction and a larger regional effect than
# the subject-level Combined effect.
effect_screen <- diff_all %>%
  mutate(rank_combined = percent_rank(abs_combined),
         rank_region = percent_rank(-max_abs_region)) %>%
  filter(rank_combined <= 0.3, rank_region <= 0.3) %>%
  pull(pathway)
selected_pathways <- diff_all %>%
  filter(!is.na(interaction_padj), interaction_padj < 0.05,
         max_abs_region > abs_combined) %>%
  pull(pathway)
diff_all$effect_screen <- diff_all$pathway %in% effect_screen
diff_all$formal_selected <- diff_all$pathway %in% selected_pathways

message("Pathways passing the effect-size screen: ", length(effect_screen))
message("Pathways passing the interaction FDR threshold: ", length(selected_pathways))
print(selected_pathways)

write.csv(diff_all, file.path(output, "04.diff_all.csv"), row.names = FALSE)
write.csv(data.frame(pathway = selected_pathways),
          file.path(output, "04.selected_pathways_interaction_FDR.csv"),
          row.names = FALSE)
write.csv(model_tests, file.path(output, "05.boxplot_model_BH.csv"),
          row.names = FALSE, quote = TRUE)
# Historical compatibility alias; contents are now repeated-measures model
# contrasts, not independent-sample Wilcoxon tests.
write.csv(model_tests, file.path(output, "05.boxplot_wilcox_BH.csv"),
          row.names = FALSE, quote = TRUE)
write.csv(interaction_tests, file.path(output, "05.boxplot_interaction_BH.csv"),
          row.names = FALSE, quote = TRUE)



output_dir <- "./16_Integrated_Pathway_Analysis/05.boxplot"
dir.create(output_dir, showWarnings = F, recursive = T)

# Plot stars use the repeated-measures model contrasts calculated above. The
# BH correction covers the complete pathway-by-comparison contrast family.
boxplot_tests <- model_tests

for (i in 1:nrow(gsva_mat)) {
  
  pathway_name <- rownames(gsva_mat)[i]
  score <- as.data.frame(t(gsva_mat[i, , drop=F]))
  colnames(score) <- "score"
  
  score$id <- rownames(score)
  plot_df <- merge(score, group1, by = "id")
  
  plot_df <- plot_df %>%
    mutate(
      Tissue = case_when(
        grepl("^Lung", group) ~ "Lung",
        grepl("^Breast", group) ~ "Breast"
      ),
      Region = case_when(
        grepl("tumor", group, ignore.case=T) ~ "Tumor",
        grepl("normal adjacent", group, ignore.case=T) ~ "Normal adjacent",
        grepl("normal", group, ignore.case=T) ~ "Normal",
        grepl("metastasis", group, ignore.case=T) ~ "Tumor"
      )
    )
  
  # Combined is a subject-level mean across the three regions (n=3 per
  # tissue), not nine independent observations. Region panels retain the
  # three repeated regional observations per tissue. No paired lines are
  # drawn; the repeated-measures model is recorded in the statistical tables.
  combined_plot_df <- plot_df %>%
    mutate(subject_id = paste(Tissue, replicate, sep = "__")) %>%
    group_by(Tissue, subject_id) %>%
    summarise(score = mean(score, na.rm = TRUE), .groups = "drop") %>%
    mutate(Compare = "Combined")

  region_plot_df <- plot_df %>%
    mutate(
      Compare = case_when(
        Region == "Normal" ~ "Normal",
        Region == "Normal adjacent" ~ "Normal adjacent",
        Region == "Tumor" ~ "Tumor",
        TRUE ~ NA_character_
      )
    ) %>%
    filter(!is.na(Compare))

  plot_data <- bind_rows(combined_plot_df, region_plot_df) %>%
    mutate(Compare = factor(Compare, levels = c("Combined","Normal","Normal adjacent","Tumor")))
  
  
  p <- ggplot(plot_data, aes(x = Compare, y = score, fill = Tissue)) +
    geom_boxplot(width = 0.5, outlier.shape = NA, position = position_dodge(0.8)) +
    geom_jitter(alpha = 0.4, size = 1, position = position_dodge(0.8)) +
    geom_text(
      data = boxplot_tests[boxplot_tests$pathway == pathway_name, , drop = FALSE],
      aes(x = Compare, y = max(plot_data$score, na.rm = TRUE) * 1.1, label = sig),
      inherit.aes = FALSE, size = 5, fontface = "bold"
    ) +
    
    scale_fill_manual(values = c("Lung" = "#1f77b4", "Breast" = "#db534c")) +
    labs(title = pathway_name, x = "", y = "ssGSEA Score") +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, size=14, face="bold"),
      axis.text.x = element_text(angle=45, hjust=1, size=12, face="bold"),
      legend.position = "top"
    )

  fname <- gsub("[^a-zA-Z0-9_-]", "_", pathway_name)
  ggsave(file.path(output_dir, paste0(fname, ".pdf")), p, width=9, height=6)
  ggsave(file.path(output_dir, paste0(fname, ".png")), p, width=9, height=6)
  
  message("Completed: ", pathway_name)
}



message("16_integrated_pathway_analysis completed.")



