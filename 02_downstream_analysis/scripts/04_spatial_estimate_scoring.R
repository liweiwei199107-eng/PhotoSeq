rm(list = ls()); gc()
ORIGINAL_DIR <- "."
output <- file.path(ORIGINAL_DIR, "04_Estimate_Spatial")
if (!dir.exists(output)) {dir.create(output, recursive = TRUE)}
setwd(ORIGINAL_DIR)

library(utils)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(IOBR)
library(homologene)

load('00_prepared_data/01.Spatial_dat_norm_ComBatSeq.RData')
group<-read.csv('./00_prepared_data/01.Spatial_group.csv')

data <- data.frame(dat_norm)
data$ID <- rownames(data)
conv <- homologene(data$ID,
                   inTax = 10116, #rat
                   outTax = 9606)  # human

probe2symbol <- conv[,c(1,2)]
colnames(probe2symbol)<-c('ID','symbol')

dat <-data
dat<-dat %>%
  inner_join(probe2symbol,by='ID')%>% 
  dplyr::select(-ID)%>%     
  dplyr::select(symbol,everything())%>%     
  mutate(rowMean=rowMeans(.[grep('GSM.',names(.))]))%>%    
  dplyr::arrange(desc(rowMean))%>%       
  distinct(symbol,.keep_all = T)%>%      
  dplyr::select(-rowMean)%>%     
  tibble::column_to_rownames(colnames(.)[1])  


ESTIMATE_result <- deconvo_tme(eset = dat,
                               method = "estimate")

rownames(ESTIMATE_result) <- group$id

write.csv(ESTIMATE_result, file.path(output, "01.ESTIMATE_result.csv"), quote = FALSE, row.names = FALSE)

identical(rownames(group),rownames(ESTIMATE_result))
ESTIMATE_result$group <- factor(group$group)
colnames(ESTIMATE_result) <- c("ID","StromalScore","ImmuneScore","ESTIMATEScore",
                               "TumorPurity","group")

ESTIMATE_result$group_plot <- gsub("_", " ", ESTIMATE_result$group)

comparisons_list <- list(
  c("Breast_tumor", "Breast_normal_adjacent"),
  c("Lung_metastasis", "Breast_tumor")
)

metrics <- c("StromalScore","ImmuneScore","ESTIMATEScore","TumorPurity")
metric_names <- c("Stromal Score","Immune Score","ESTIMATE Score","Tumor Purity")
names(metric_names) <- metrics

fill_cols <- c('#db534c','#1f77b5')

sig_label <- function(p) {
  ifelse(is.na(p), "ns", ifelse(p < 0.0001, "****", ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", "ns")))))
}

estimate_stats <- do.call(rbind, lapply(comparisons_list, function(comp) {
  g1 <- comp[1]; g2 <- comp[2]
  sub <- ESTIMATE_result[ESTIMATE_result$group %in% c(g1, g2), , drop = FALSE]
  do.call(rbind, lapply(metrics, function(m) {
    data.frame(group1 = g1, group2 = g2, metric = m,
               pvalue = tryCatch(wilcox.test(sub[[m]] ~ sub$group)$p.value, error = function(e) NA_real_))
  }))
}))
estimate_stats$padj <- p.adjust(estimate_stats$pvalue, method = "BH")
estimate_stats$sig <- sig_label(estimate_stats$padj)
write.csv(estimate_stats, file.path(output, "02.ESTIMATE_wilcox_FDR.csv"), quote = FALSE, row.names = FALSE)

for (comp in comparisons_list) {
  g1 <- comp[1]
  g2 <- comp[2]
  
  g1_plot <- gsub("_", " ", g1)
  g2_plot <- gsub("_", " ", g2)
  
  fold_name <- paste0(g1,"_vs_",g2)
  fold_path <- file.path(output, fold_name)
  dir.create(fold_path, recursive = TRUE, showWarnings = FALSE)
  
  dat_sub <- ESTIMATE_result %>% 
    filter(group %in% c(g1,g2)) %>% droplevels()
  
  dat_sub$group_plot <- factor(dat_sub$group_plot, levels = c(g2_plot, g1_plot))
  
  plot_list <- list()
  
  for (m in metrics) {
    p <- ggviolin(
      dat_sub,
      x = "group_plot",       
      y = m,
      fill = "group_plot",     
      palette = fill_cols,
      add = "boxplot",
      add.params = list(fill = "white")
    ) +
      labs(y = metric_names[m]) +
      theme_bw() +
      theme(
        panel.grid = element_blank(),
        axis.title.x = element_blank(),
        axis.text = element_text(size = 14, face = "bold"),
        axis.title.y = element_text(size = 14, face = "bold"),
        legend.position = "none"
      )
    stat_row <- estimate_stats[estimate_stats$group1 == g1 & estimate_stats$group2 == g2 & estimate_stats$metric == m, , drop = FALSE]
    y_pos <- max(dat_sub[[m]], na.rm = TRUE)
    if (!is.finite(y_pos)) y_pos <- 0
    p <- p + annotate("text", x = 1.5, y = y_pos + max(abs(y_pos) * 0.08, 0.1), label = stat_row$sig, size = 5)
    
    ggsave(file.path(fold_path, paste0(m,"_violin.pdf")), p, width = 6, height = 5)
    ggsave(file.path(fold_path, paste0(m,"_violin.png")), p, width = 6, height = 5)
    plot_list[[m]] <- p
  }
  
  combine_p <- ggarrange(
    plot_list[[1]], plot_list[[2]],
    plot_list[[3]], plot_list[[4]],
    ncol = 2, nrow = 2, align = "hv"
  )
  
  ggsave(file.path(fold_path, paste0("ESTIMATE_4index_combined.pdf")),
         combine_p, width = 12, height = 10)
  ggsave(file.path(fold_path, paste0("ESTIMATE_4index_combined.png")),
         combine_p, width = 12, height = 10)
}
