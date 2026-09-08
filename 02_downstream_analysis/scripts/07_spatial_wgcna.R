rm(list = ls()); gc()
ORIGINAL_DIR <- "."
output <- file.path(ORIGINAL_DIR, "07_WGCNA_Spatial")
if (!dir.exists(output)) {dir.create(output, recursive = TRUE)}
setwd(ORIGINAL_DIR)

library(readr)
library(tibble)
library(dbplyr)
library(WGCNA)



load('00_prepared_data/01.Spatial_dat_norm_ComBatSeq.RData')
dat_expr<-dat_norm
head(dat_expr)[, 1:6]
range(dat_expr)
#[1]  0.00000 17.45999

group<-read.csv('./00_prepared_data/01.Spatial_group.csv')%>% setNames(c('sample', 'type', 'sequencing_batch'))
stopifnot(nrow(group) == 18, all(c('Batch1','Batch2') %in% group$sequencing_batch))
group$type <- gsub("_", " ", group$type)

mat <- as.data.frame(t(dat_expr))

temp1 <- goodSamplesGenes(mat, verbose = 3) 
temp1$allOK
if (!temp1$allOK){
  if (sum(!temp1$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(mat)[!temp1$goodGenes], collapse = ", ")));
  if (sum(!temp1$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(mat)[!temp1$goodSamples], collapse = ", ")));
  mat = mat[temp1$goodSamples, temp1$goodGenes]
}
temp1 <- goodSamplesGenes(mat, verbose = 3) 
temp1$allOK


sampleTree = hclust(dist(mat, method = 'euclidean'), 
                    method = "average")


pdf(file.path(output,file = "01.SampleTree.pdf"), width = 8, height = 6) 
par(cex = 0.6) 
par(mar = c(0,6,6,0)) 
plot(sampleTree, 
     xlab = "", 
     sub = "", 
     main = "Sample Clustering",
     # labels = FALSE,  # hide sample names
     cex = 2.0,  
     font = 2,
     cex.axis = 1.6,  
     cex.lab = 1.8,   
     cex.main = 1.8,   
     font.axis = 2,
     font.lab = 2,
     font.main = 2,
     font.sub = 2)
dev.off()

png(file.path(output,"01.SampleTree.png"),width = 10,height = 6,units = "in",res = 600)
# plotDendroAndColors(sampleTree, traitColors,
#                     groupLabels =names(datTraits),
#                     main ="Sample dendrogramand trait heatmap",
#                     dendroLabels = FALSE)
par(cex = 0.6) 
par(mar = c(0,6,6,0)) 
plot(sampleTree, 
     xlab = "", 
     sub = "", 
     main = "Sample Clustering",
     # labels = FALSE,  # hide sample names
     cex = 2.0,  
     font = 2,
     cex.axis = 1.6,  
     cex.lab = 1.8,   
     cex.main = 1.8,   
     font.axis = 2,
     font.lab = 2,
     font.main = 2,
     font.sub = 2)
dev.off()



clust = cutreeStatic(sampleTree, cutHeight = 250, minSize = 10)
table(clust) 
# The three Breast Mid(Tumor) samples form a biological group, not technical
# outliers. Keep all good samples so all six spatial groups are represented.
mat1 <- mat
ncol(mat1); nrow(mat1)
# [1] 13881
# [1] 18


# Enable limited parallelism for the same WGCNA calculation; this changes
# execution speed only, not the network or module parameters.
allowWGCNAThreads(nThreads = 4)
powers <- c(seq(1, 10, by = 1), seq(12, 20, by = 2)) 


num <- 0.85 # R-squared threshold; typical choices are 0.8, 0.85 or 0.9
sft <- pickSoftThreshold(mat1, RsquaredCut = num, powerVector = powers, verbose = 5)

sft$powerEstimate
#save(sft,file = "./03_WGCNA/powerEstimate.RData")
#load('./03_WGCNA/powerEstimate.RData')



pdf(file.path(output,'02.SoftThreshold.pdf'), w = 11, h = 7)
par(mfrow = c(1, 2), pin = c(4, 4), mar = c(6, 6, 6, 1))
cex1 = 1.5
plot(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3])*sft$fitIndices[, 2], 
     xlab = "Soft Threshold (power)", ylab = "Scale Free Topology Model Fit,signed R^2", type = "n",
     cex.axis = 1.4,
     cex.lab = 1.8,
     cex.main = 1.8,
     font.axis = 2,
     font.lab = 2,
     font.main = 2,
     font.sub = 2,
     main = paste("Scale independence"));
text(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3])*sft$fitIndices[, 2], labels = powers, cex = cex1, col = "red");
abline(h = num, col = "red")
plot(sft$fitIndices[, 1], sft$fitIndices[, 5], 
     xlab = "Soft Threshold (power)", ylab = "Mean Connectivity", type = "n", main = paste("Mean connectivity"),
     cex.axis = 1.6,
     cex.lab = 1.8,
     cex.main = 1.8,
     font.axis = 2,
     font.lab = 2,
     font.main = 2,
     font.sub = 2)
text(sft$fitIndices[, 1], sft$fitIndices[, 5], labels = powers, cex = cex1, col = "red")
dev.off()

png(file.path(output,'02.SoftThreshold.png'),w=10,h=6,units='in',res=600,bg='white')
par(mfrow = c(1, 2), pin = c(4, 4), mar = c(6, 6, 6, 1))
cex1 = 1.5
plot(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3])*sft$fitIndices[, 2], 
     xlab = "Soft Threshold (power)", ylab = "Scale Free Topology Model Fit,signed R^2", type = "n",
     cex.axis = 1.6,
     cex.lab = 1.8, 
     cex.main = 1.8, 
     font.axis = 2,
     font.lab = 2,
     font.main = 2,
     font.sub = 2,
     main = paste("Scale independence"));
text(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3])*sft$fitIndices[, 2], labels = powers, cex = cex1, col = "red");
abline(h = num, col = "red")
plot(sft$fitIndices[, 1], sft$fitIndices[, 5], 
     xlab = "Soft Threshold (power)", ylab = "Mean Connectivity", type = "n", main = paste("Mean connectivity"),
     cex.axis = 1.6,
     cex.lab = 1.8,
     cex.main = 1.8,
     font.axis = 2,
     font.lab = 2,
     font.main = 2,
     font.sub = 2)
text(sft$fitIndices[, 1], sft$fitIndices[, 5], labels = powers, cex = cex1, col = "red")
dev.off()

datExpr <- mat1
powers <- 8
cor <- WGCNA::cor
net <- blockwiseModules(
  datExpr,
  power = powers,
  minModuleSize = 250, # minimum genes per module
  deepSplit = 4,                       
  mergeCutHeight = 0.25, # module merge parameter; larger values produce fewer modules
  numericLabels = TRUE,                 
  networkType  = "signed",              
  maxBlockSize = ncol(datExpr),
  pamRespectsDendro = FALSE,
  # TOM files are generated separately for each input matrix.
  saveTOMs = FALSE,
  loadTOMs = FALSE,
  nThreads = 4,
  verbose = 3
)
#save(net,file = "./03_WGCNA/net_200.RData")
#load('./02_WGCNA/net_200.RData')


cor <- stats::cor
table(net$colors)

moduleColors <- labels2colors(net$colors)
MEs0 <- moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs <- orderMEs(MEs0)
useMEs <- subset(MEs, select = -c(MEgrey))


# Plot modules
mergedColors <-  labels2colors(net$colors)
mergedColors[net$blockGenes[[1]]]
pdf(file.path(output,"03.Wgcna_dendroColors.pdf"), height = 7, width = 9)
par(mfrow = c(1, 2), pin = c(4, 4), mar = c(6, 6, 6, 1))
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]], "Module colors",
                    dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05,
                    cex.axis = 1.6,
                    cex.lab = 1.8,
                    cex.main = 1.8,
                    font.axis = 2,
                    font.lab = 2,
                    font.main = 2,
                    cex.colorLabels = 1,
                    font.sub = 2)
dev.off()

png(file.path(output,"03.Wgcna_dendroColors.png"), height = 7, width = 9, units='in', res = 600)
par(mfrow = c(1, 2), pin = c(4, 4), mar = c(6, 6, 6, 1))
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]], "Module colors",
                    dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05,
                    cex.axis = 1.6,
                    cex.lab = 1.8,
                    cex.main = 1.8,
                    font.axis = 2,
                    font.lab = 2,
                    font.main = 2,
                    cex.colorLabels = 1,
                    font.sub = 2)
dev.off()


group2 <- read.csv('./00_prepared_data/01.Spatial_group.csv') %>% 
  setNames(c('ID', 'type', 'sequencing_batch'))
group2 <- group2[group2$ID %in% rownames(datExpr), ]
group2$type <- gsub("_", " ", group2$type)

trait_matrix <- model.matrix(~ 0 + type, data = group2)
colnames(trait_matrix) <- gsub("^type", "", colnames(trait_matrix)) # remove the type prefix

module <- useMEs
moduleTraitCor <- cor(module, trait_matrix, use = 'p')
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nrow(module))
moduleTraitFDR <- matrix(NA_real_, nrow = nrow(moduleTraitPvalue),
                         ncol = ncol(moduleTraitPvalue),
                         dimnames = dimnames(moduleTraitPvalue))
valid_p <- is.finite(moduleTraitPvalue)
moduleTraitFDR[valid_p] <- p.adjust(moduleTraitPvalue[valid_p], method = "BH")
write.csv(as.data.frame(as.table(moduleTraitCor)), file.path(output, "06.ModuleTrait_cor_long.csv"), row.names = FALSE)
module_trait_stats <- merge(as.data.frame(as.table(moduleTraitCor)), as.data.frame(as.table(moduleTraitPvalue)), by = c("Var1", "Var2"), suffixes = c("_r", "_p"))
fdr_long <- as.data.frame(as.table(moduleTraitFDR))
names(fdr_long) <- c("Var1", "Var2", "FDR")
module_trait_stats <- merge(module_trait_stats, fdr_long, by = c("Var1", "Var2"))
module_trait_stats$significant <- abs(module_trait_stats$Freq_r) > 0.3 & module_trait_stats$FDR < 0.05
write.csv(module_trait_stats, file.path(output, "06.ModuleTrait_statistics.csv"), row.names = FALSE)
write.csv(subset(module_trait_stats, significant), file.path(output, "07.ModuleTrait_significant_FDR.csv"), row.names = FALSE)

moduleTraitCor <- moduleTraitCor[!(row.names(moduleTraitCor) %in% c("MEgrey")), ]
moduleTraitPvalue <- moduleTraitPvalue[!(row.names(moduleTraitPvalue) %in% c("MEgrey")), ]

textMatrix <- matrix(
  paste(signif(moduleTraitCor, 2), '\n(', signif(moduleTraitPvalue, 1), ')', sep = ''),
  nrow = nrow(moduleTraitCor),
  dimnames = dimnames(moduleTraitCor)
)

pdf(file.path(output,"04.Module-trait_heatmap.pdf"), width = 12, height = 10)
par(mar = c(8, 10, 3, 3)) 
labeledHeatmap(
  Matrix = as.matrix(moduleTraitCor),
  main = paste('Module-trait relationships'),
  xLabels = colnames(moduleTraitCor),
  yLabels = rownames(moduleTraitCor),
  ySymbols = rownames(moduleTraitCor),
  colorLabels = FALSE,
  colors = blueWhiteRed(50),
  cex.text = 1, 
  zlim = c(-1,1),
  #textMatrix = textMatrix,
  setStdMargins = FALSE
)
dev.off()


png(file.path(output,"04.Module-trait_heatmap.png"),width = 10,height = 10,units = "in",res = 600)
par(mar = c(8, 10, 3, 3)) 
labeledHeatmap(
  Matrix = as.matrix(moduleTraitCor),
  main = paste('Module-trait relationships'),
  xLabels = colnames(moduleTraitCor),
  yLabels = rownames(moduleTraitCor),
  ySymbols = rownames(moduleTraitCor),
  colorLabels = FALSE,
  colors = blueWhiteRed(50),
  cex.text = 1, 
  zlim = c(-1,1),
  #textMatrix = textMatrix,
  setStdMargins = FALSE
)
dev.off()




gene_module <- data.frame(gene_name = colnames(datExpr),
                          module = mergedColors, stringsAsFactors = FALSE)
write.csv(gene_module, file.path(output,"05.WGCNA_gene_all.csv"),quote = F, row.names = F)
