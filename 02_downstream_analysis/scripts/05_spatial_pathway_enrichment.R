rm(list = ls()); gc()
ORIGINAL_DIR <- "."
output <- file.path(ORIGINAL_DIR, "05_Enrich_Spatial")
if (!dir.exists(output)) {dir.create(output, recursive = TRUE)}
setwd(ORIGINAL_DIR)

library(clusterProfiler)
library(org.Mm.eg.db)
library(ggplot2)
library(stringr)
library(ReactomePA)
library(readr)

###### Convert DEG identifiers to ENTREZID
genes <- read.csv("02_Venn_Spatial/02.Venn_DEGs.csv", stringsAsFactors = FALSE)

DEG.entrez_id = mapIds(x = org.Mm.eg.db,
                       keys = genes$x,
                       keytype = "SYMBOL",
                       column = "ENTREZID")

DEG.entrez_id = na.omit(DEG.entrez_id)

set.seed(1)
####### GO
erich.go.BP  = enrichGO(gene = DEG.entrez_id,
                        OrgDb = org.Mm.eg.db,
                        keyType = "ENTREZID",
                        ont = "BP",
                        pAdjustMethod = "BH",
                        pvalueCutoff = 1,
                        qvalueCutoff = 1,
                        readable = TRUE)
#View(erich.go.all@result)
erich.go.BP@result <- erich.go.BP@result %>% arrange(p.adjust)
df <- erich.go.BP@result
write.csv(df, file.path(output, '01.GO_BP_res.csv'), quote = T, row.names = TRUE)




BP <- c('cell-substrate adhesion',
        'cell-matrix adhesion',
        'proteasome-mediated ubiquitin-dependent protein catabolic process',
        'endothelial cell migration',
        'immune response-activating signaling pathway',
        'epidermal growth factor receptor signaling pathway',
        'immune response-regulating signaling pathway',
        'negative regulation of phosphorus metabolic process',
        'regulation of protein catabolic process',
        'dephosphorylation'
)
data <- df[df$Description %in% BP & !is.na(df$p.adjust) & df$p.adjust < 0.05,]%>%  arrange(p.adjust)


data$Description=factor(data$Description, levels = rev(data$Description))
p1 = ggplot(data,aes(x=Count,y=Description))+
  geom_point()+
  geom_point(aes(size=Count,color=-log10(p.adjust)))+
  scale_color_gradient(low='#db534c',high="#fbbcba", guide = guide_colorbar(display = "rectangles"))+
  scale_x_continuous(limits = c(0, max(data$Count)+2),  
                     expand = c(0, 0))+
  labs(title="Biological Process",x="Count",y="",color=expression(-log[10](adjusted~italic(P)~value)),size="Gene number")+
  scale_y_discrete(labels = function(x) stringr::str_wrap(x, width = 40))+
  theme_bw()+
  theme(aspect.ratio=2/1,
        plot.title = element_text(hjust = 0.5, size = 15, face = "bold"),
        #legend.position = "none",
        legend.title = element_text(size = 10, face = "bold", color = "black"),
        legend.text = element_text(size = 10, face = "bold", color = "black"),
        axis.title = element_text(size = 10, face = "bold", color = "black"),
        axis.text = element_text(size = 10, face = "bold", color = "black"),
        axis.text.y = element_text(size = 7, face = "plain", color = "black"),
        axis.text.x = element_text(angle = 0, hjust = 0.5, size = 12, face = "bold", color = "black")
  )+
  guides(size = 'none')
p1
ggsave(file.path(output, '01.GO_BP.pdf'), p1, width = 6, height = 5)
ggsave(file.path(output, '01.GO_BP.png'), p1, width = 6, height = 5)

erich.go.CC  = enrichGO(gene = DEG.entrez_id,
                        OrgDb = org.Mm.eg.db,
                        keyType = "ENTREZID",
                        ont = "CC",
                        pAdjustMethod = "BH",
                        pvalueCutoff = 1,
                        qvalueCutoff = 1,
                        readable = TRUE)

erich.go.CC@result <- erich.go.CC@result %>% arrange(p.adjust)
df <- erich.go.CC@result
write.csv(df, file.path(output, '02.GO_CC_res.csv'), quote = T, row.names = TRUE)


data <- head(df[!is.na(df$p.adjust) & df$p.adjust < 0.05, , drop = FALSE], 10)%>%  arrange(p.adjust)
data$Description=factor(data$Description, levels = rev(data$Description))
p2 = ggplot(data,aes(x=Count,y=Description))+
  geom_point()+
  geom_point(aes(size=Count,color=-log10(p.adjust)))+
  scale_color_gradient(low='#db534c',high="#fbbcba", guide = guide_colorbar(raster = FALSE))+
  scale_x_continuous(limits = c(0, max(data$Count)+5),  
                     expand = c(0, 0))+
  labs(title="Cellular Components",x="Count",y="",color=expression(-log[10](adjusted~italic(P)~value)),size="Gene number")+
  scale_y_discrete(labels = function(x) stringr::str_wrap(x, width = 40))+
  theme_bw()+
  theme(aspect.ratio=2/1,
        plot.title = element_text(hjust = 0.5, size = 15, face = "bold"),
        #legend.position = "none",
        legend.title = element_text(size = 10, face = "bold", color = "black"),
        legend.text = element_text(size = 10, face = "bold", color = "black"),
        axis.title = element_text(size = 10, face = "bold", color = "black"),
        axis.text = element_text(size = 10, face = "bold", color = "black"),
        axis.text.x = element_text(angle = 0, hjust = 0.5, size = 12, face = "bold", color = "black")
  )+
  guides(size = 'none')
p2
ggsave(file.path(output, '02.GO_CC.pdf'), p2, width = 6, height = 5)
ggsave(file.path(output, '02.GO_CC.png'), p2, width = 6, height = 5)

erich.go.MF  = enrichGO(gene = DEG.entrez_id,
                        OrgDb = org.Mm.eg.db,
                        keyType = "ENTREZID",
                        ont = "MF",
                        pAdjustMethod = "BH",
                        pvalueCutoff = 1,
                        qvalueCutoff = 1,
                        readable = TRUE)

erich.go.MF@result <- erich.go.MF@result %>% arrange(p.adjust)
df <- erich.go.MF@result
write.csv(df, file.path(output, '03.GO_MF_res.csv'), quote = T, row.names = TRUE)


data <- head(df[!is.na(df$p.adjust) & df$p.adjust < 0.05, , drop = FALSE], 10)%>%  arrange(p.adjust)
data$Description=factor(data$Description, levels = rev(data$Description))
p3 = ggplot(data,aes(x=Count,y=Description))+
  geom_point()+
  geom_point(aes(size=Count,color=-log10(p.adjust)))+
  scale_color_gradient(low='#db534c',high="#fbbcba")+
  scale_x_continuous(limits = c(0, max(data$Count)+5),  
                     expand = c(0, 0))+
  labs(title="Molecular Functions",x="Count",y="",color=expression(-log[10](adjusted~italic(P)~value)),size="Gene number")+
  scale_y_discrete(labels = function(x) stringr::str_wrap(x, width = 40))+
  theme_bw()+
  theme(aspect.ratio=2/1,
        plot.title = element_text(hjust = 0.5, size = 15, face = "bold"),
        #legend.position = "none",
        legend.title = element_text(size = 10, face = "bold", color = "black"),
        legend.text = element_text(size = 10, face = "bold", color = "black"),
        axis.title = element_text(size = 10, face = "bold", color = "black"),
        axis.text = element_text(size = 10, face = "bold", color = "black"),
        axis.text.x = element_text(angle = 0, hjust = 0.5, size = 12, face = "bold", color = "black")
  )+
  guides(size = 'none')
p3
ggsave(file.path(output, '03.GO_MF.pdf'), p3, width = 6, height = 5)
ggsave(file.path(output, '03.GO_MF.png'), p3, width = 6, height = 5)


set.seed(1)
######### KEGG
# Use the project-supplied KEGG snapshot so the rerun is reproducible and does
# not depend on a live KEGG REST request.
kegg_link <- read.delim("reference/kegg/kegg_link_mmu_pathway.txt", header = FALSE,
                        sep = "\t", stringsAsFactors = FALSE,
                        col.names = c("pathway", "gene"))
kegg_link$pathway <- sub("^path:", "", kegg_link$pathway)
kegg_link$gene <- sub("^mmu:", "", kegg_link$gene)
kegg_name <- read.delim("reference/kegg/kegg_list_pathway_mmu.txt", header = FALSE,
                        sep = "\t", stringsAsFactors = FALSE,
                        col.names = c("pathway", "Description"))
kegg_name$pathway <- sub("^path:", "", kegg_name$pathway)
enrich.KEGG <- clusterProfiler::enricher(
  gene = as.character(DEG.entrez_id),
  universe = unique(kegg_link$gene),
  TERM2GENE = kegg_link[, c("pathway", "gene")],
  TERM2NAME = kegg_name[, c("pathway", "Description")],
  pAdjustMethod = "BH", pvalueCutoff = 1, qvalueCutoff = 1,
  minGSSize = 1, maxGSSize = Inf)

enrich.KEGG = clusterProfiler::setReadable(enrich.KEGG,OrgDb = 'org.Mm.eg.db',keyType = 'ENTREZID')
#View(enrich.KEGG@result)
enrich.KEGG@result <- enrich.KEGG@result %>% arrange(p.adjust)
write.csv(enrich.KEGG@result, file.path(output, '04.KEGG_res.csv'), quote =T, row.names = TRUE)
df <- data.frame(enrich.KEGG@result)
#57

KEGG <- c("ECM-receptor interaction - Mus musculus (house mouse)",
          'Focal adhesion - Mus musculus (house mouse)',
          'Tight junction - Mus musculus (house mouse)',
          'PI3K-Akt signaling pathway - Mus musculus (house mouse)',
          'MAPK signaling pathway - Mus musculus (house mouse)',
          'Wnt signaling pathway - Mus musculus (house mouse)',
          'Notch signaling pathway - Mus musculus (house mouse)',
          'Glycerophospholipid metabolism - Mus musculus (house mouse)',
          'Selenocompound metabolism - Mus musculus (house mouse)',
          'Propanoate metabolism - Mus musculus (house mouse)')

data <- df[df$Description %in% KEGG & !is.na(df$p.adjust) & df$p.adjust < 0.05,]%>%  arrange(p.adjust)
data$Description <- sub(" - Mus musculus \\(house mouse\\)$", "", data$Description)
data$Description=factor(data$Description, levels = rev(data$Description))

p5 = ggplot(data,aes(x=Count,y=Description))+
  geom_point()+
  geom_point(aes(size=Count,color=-log10(p.adjust)))+
  scale_color_gradient(low='#db534c',high="#fbbcba")+
  scale_x_continuous(limits = c(0, max(data$Count)+5),  
                     expand = c(0, 0))+
  labs(title="KEGG",x="Count",y="",color=expression(-log[10](adjusted~italic(P)~value)),size="Gene number")+
  scale_y_discrete(labels = function(x) stringr::str_wrap(x, width = 40))+
  theme_bw()+
  theme(aspect.ratio=2/1,
        plot.title = element_text(hjust = 0.5, size = 15, face = "bold"),
        #legend.position = "none",
        legend.title = element_text(size = 10, face = "bold", color = "black"),
        legend.text = element_text(size = 10, face = "bold", color = "black"),
        axis.title = element_text(size = 10, face = "bold", color = "black"),
        axis.text = element_text(size = 10, face = "bold", color = "black"),
        axis.text.x = element_text(angle = 0, hjust = 0.5, size = 12, face = "bold", color = "black")
  )+
  guides(size = 'none')
p5
ggsave(file.path(output, '04.KEGG.pdf'), p5, width = 6, height = 5)
ggsave(file.path(output, '04.KEGG.png'), p5, width = 6, height = 5)


eReac <- enrichPathway(gene = DEG.entrez_id,
                       organism = 'mouse',
                       pAdjustMethod = "BH",
                       pvalueCutoff = 1,
                       qvalueCutoff = 1,
                       readable = TRUE)
eReac@result <- eReac@result %>% arrange(p.adjust)
df <- eReac@result
write.csv(df, file.path(output, '05.Reactome_res.csv'), quote = T, row.names = TRUE)


data <- head(df[!is.na(df$p.adjust) & df$p.adjust < 0.05, , drop = FALSE], 10)%>%  arrange(p.adjust)
data$Description=factor(data$Description, levels = rev(data$Description))
p6 = ggplot(data,aes(x=Count,y=Description))+
  geom_point()+
  geom_point(aes(size=Count,color=-log10(p.adjust)))+
  scale_color_gradient(low='#db534c',high="#fbbcba")+
  scale_x_continuous(limits = c(0, max(data$Count)+5),  
                     expand = c(0, 0))+
  labs(title="Reactome pathway enrichment",x="Count",y="",color=expression(-log[10](adjusted~italic(P)~value)),size="Gene number")+
  scale_y_discrete(labels = function(x) stringr::str_wrap(x, width = 40))+
  theme_bw()+
  theme(aspect.ratio=2/1,
        plot.title = element_text(hjust = 0.5, size = 15, face = "bold"),
        #legend.position = "none",
        legend.title = element_text(size = 10, face = "bold", color = "black"),
        legend.text = element_text(size = 10, face = "bold", color = "black"),
        axis.title = element_text(size = 10, face = "bold", color = "black"),
        axis.text = element_text(size = 10, face = "bold", color = "black"),
        axis.text.y = element_text(size = 7, face = "plain", color = "black"),
        axis.text.x = element_text(angle = 0, hjust = 0.5, size = 12, face = "bold", color = "black")
  )+
  guides(size = 'none')
p6
ggsave(file.path(output, '05.Reactome.pdf'), p6, width = 6, height = 5)
ggsave(file.path(output, '05.Reactome.png'), p6, width = 6, height = 5)

