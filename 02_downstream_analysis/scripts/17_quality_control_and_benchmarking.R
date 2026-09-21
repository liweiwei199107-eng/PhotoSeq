rm(list = ls()); gc()
ORIGINAL_DIR <- "."
out_root <- file.path(ORIGINAL_DIR, "17_Quality_Control")
dir.create(out_root, recursive = TRUE, showWarnings = FALSE)
suppressPackageStartupMessages({
  library(ggplot2); library(ComplexHeatmap); library(circlize)
})
if (!requireNamespace("readxl", quietly = TRUE)) stop("Package 'readxl' is required for the UV QC source workbook")

# QC contract: figures 1, 2, 3, 6 and 7 have no explicit batch correction;
# figure 4 uses the integrated ComBat_seq -> VST matrix; figure 5 uses
# batch-aware DE tables with padj < .01 and |log2FC| > 1.5.
save_plot <- function(p, stem, width = 180, height = 130) {
  ggsave(paste0(stem, ".pdf"), p, width = width, height = height, units = "mm",
         device = grDevices::cairo_pdf, family = "Arial", bg = "white")
  ggsave(paste0(stem, ".png"), p, width = width, height = height, units = "mm", dpi = 600, bg = "white")
}
load_object <- function(path, object) {
  e <- new.env(parent = emptyenv()); load(path, envir = e)
  if (!exists(object, envir = e, inherits = FALSE)) stop("Missing ", object, " in ", path)
  get(object, envir = e, inherits = FALSE)
}

spatial_counts <- round(load_object("00_prepared_data/01.Spatial_dat_count.RData", "dat1"))
temporal_counts <- round(load_object("00_prepared_data/02.Temporal_dat_count.RData", "dat1"))
integrated_vst <- as.matrix(load_object("00_prepared_data/03.All_dat_norm.RData", "dat_norm"))
all_group <- read.csv("00_prepared_data/03.All_group.csv", stringsAsFactors = FALSE, check.names = FALSE)

uv_source <- file.path("reference", "uv", "uv_pair_source_data.xlsx")
et_source <- file.path("input", "PhotoSeq_TechnicalReplicate_GeneCounts.csv")
ratio_source <- file.path(
  "reference", "benchmark",
  "Breast_Late_PhotoSeq_vs_Visium_ratio_by_transcript_length_gene_table.csv"
)
density_source <- file.path(
  "reference", "benchmark",
  "Breast_Late_PhotoSeq_vs_Visium_scatter_log2CPM_gene_table.csv"
)
overlap_source <- file.path(
  "reference", "benchmark", "detected_gene_overlap_gene_status.csv"
)
required_sources <- c(uv_source, et_source, ratio_source, density_source, overlap_source)
if (any(!file.exists(required_sources))) {
  stop("Missing required QC input(s): ", paste(required_sources[!file.exists(required_sources)], collapse = ", "))
}

write.csv(data.frame(
  figure = c("01_same_batch_technical_replicate_scatter", "02_PhotoSeq_vs_Visium_ratio",
             "03_PhotoSeq_vs_Visium_density", "04_integrated_Pearson",
             "05_DEG_zscore_heatmap", "06_detected_gene_overlap",
             "07_UV_power_condition_optimization", "07_UV_exposure_time_condition_optimization"),
  batch_handling = c("none; supplied same-batch ET technical replicates; raw count -> log2(count+1)", "none; supplied PhotoSeq CPM + external Visium reference",
                     "none; supplied PhotoSeq CPM + external Visium reference", "existing ComBat_seq -> VST",
                     "DESeq2 raw-count model includes sequencing_batch; expression display uses integrated ComBat_seq -> VST",
                     "none; supplied PhotoSeq CPM + external Visium reference",
                     "none; same-source condition optimization fluorescence intensity",
                     "none; same-source condition optimization fluorescence intensity"),
  padj = c("not applicable", "not applicable", "not applicable", "not applicable", "< 0.01", "not applicable", "not applicable", "not applicable"),
  log2fc = c("not applicable", "not applicable", "not applicable", "not applicable", "absolute value > 1.5", "not applicable", "not applicable", "not applicable"),
  stringsAsFactors = FALSE
), file.path(out_root, "QC_parameters.csv"), row.names = FALSE)

# External Visium source tables are used only on the Visium side of QC2/3/6;
# the PhotoSeq side is recalculated from the supplied PhotoSeq counts.
write.csv(data.frame(
  source = c("00_prepared_data/01.Spatial_dat_count.RData", "00_prepared_data/02.Temporal_dat_count.RData",
             "00_prepared_data/03.All_dat_norm.RData",
             ratio_source, density_source, overlap_source, et_source, uv_source),
  role = c("spatial counts", "temporal counts", "integrated ComBat_seq -> VST",
           "external Visium reference for QC2", "external Visium reference for QC3", "external Visium reference for QC6",
           "technical-replicate count source", "UV condition-optimization source workbook"),
  stringsAsFactors = FALSE
), file.path(out_root, "source_provenance.csv"), row.names = FALSE)

# 01: supplied same-batch technical-replicate scatter; no batch correction.
qc1_dir <- file.path(out_root, "01.same_batch_replicate_scatter")
dir.create(qc1_dir, recursive = TRUE, showWarnings = FALSE)
et_raw <- read.csv(et_source, check.names = FALSE, stringsAsFactors = FALSE)
if (ncol(et_raw) != 3L) stop("ET technical-replicate source must contain Gene plus exactly two replicate columns")
if (anyDuplicated(et_raw[[1]])) stop("Duplicate gene identifiers in ET technical-replicate source")
et_counts <- as.matrix(et_raw[, 2:3, drop = FALSE])
storage.mode(et_counts) <- "numeric"
rownames(et_counts) <- as.character(et_raw[[1]])
if (any(!is.finite(et_counts)) || any(et_counts < 0) || any(et_counts != round(et_counts))) {
  stop("ET technical-replicate source contains non-finite, negative, or non-integer counts")
}

make_scatter <- function(counts, ids, xlab, ylab, stem, point_colour = "#6FA8DC", point_alpha = .42,
                         width = 185, height = 185) {
  x <- counts[, ids, drop = FALSE]
  d <- data.frame(
    gene = rownames(x),
    count_1 = as.numeric(x[, 1]), count_2 = as.numeric(x[, 2]),
    x = log2(pmax(as.numeric(x[, 1]), 0) + 1),
    y = log2(pmax(as.numeric(x[, 2]), 0) + 1),
    stringsAsFactors = FALSE
  )
  r <- cor(d$x, d$y, method = "pearson", use = "complete.obs")
  lim <- max(c(d$x, d$y), na.rm = TRUE)
  lim <- ceiling(lim + .25)
  p <- ggplot(d, aes(x, y)) +
    geom_point(size = ifelse(nrow(d) > 20000, .42, .65), alpha = point_alpha, colour = point_colour) +
    geom_abline(slope = 1, intercept = 0, linetype = 2, linewidth = .55, colour = "#FF2B2B") +
    annotate("text", x = lim * .13, y = lim * .78,
             label = paste0("Pearson\nR = ", sprintf("%.3f", r)), hjust = 0, vjust = 1,
             size = 3.7, family = "Arial") +
    annotate("text", x = lim * .84, y = lim * .78, label = "y = x", hjust = .5, vjust = 0,
             size = 3.5, family = "Arial") +
    labs(x = xlab, y = ylab) +
    coord_equal(xlim = c(0, lim), ylim = c(0, lim), expand = FALSE) +
    theme_classic(base_size = 9.5, base_family = "Arial") +
    theme(axis.line = element_line(linewidth = .4, colour = "#111111"),
          axis.ticks = element_line(linewidth = .35, colour = "#111111"),
          axis.title = element_text(size = 9.5, colour = "#111111"),
          axis.text = element_text(size = 8.5, colour = "#111111"),
          panel.grid = element_blank(), plot.margin = margin(6, 6, 6, 6))
  save_plot(p, file.path(qc1_dir, stem), width, height)
  write.csv(d, file.path(qc1_dir, paste0(stem, "_source_data.csv")), row.names = FALSE)
  write.csv(data.frame(sample_1 = ids[1], sample_2 = ids[2], transform = "log2(count + 1)",
                       batch_correction = "none; same-batch technical replicates",
                       pearson_r = r, n_genes = nrow(d), zero_zero = sum(d$count_1 == 0 & d$count_2 == 0)),
            file.path(qc1_dir, paste0(stem, "_summary.csv")), row.names = FALSE)
}
make_scatter(et_counts, c("Breast_early_N1", "Breast_early_N1_techrep2"),
             "log2(count + 1) Breast_early_N1", "log2(count + 1) Breast_early_N1_techrep2",
             "ET_technical_replicate_scatter", point_colour = "#6FA8DC", point_alpha = .42,
             width = 185, height = 185)

# PhotoSeq late-stage CPM, reused by QC2/3/6.
late <- temporal_counts[, grep("^Breast_late_", colnames(temporal_counts)), drop = FALSE]
late_cpm_matrix <- sweep(late, 2, colSums(late), "/") * 1e6
late_cpm <- rowMeans(late_cpm_matrix)
photo_late <- data.frame(gene_name = names(late_cpm), PhotoSeq_mean_CPM = as.numeric(late_cpm), stringsAsFactors = FALSE)

# 02: transcript-length ratio boxplot; no batch correction.
qc2_dir <- file.path(out_root, "02.PhotoSeq_vs_Visium_ratio")
dir.create(qc2_dir, recursive = TRUE, showWarnings = FALSE)
ratio_ref <- read.csv(ratio_source, stringsAsFactors = FALSE)
ratio_dat <- merge(photo_late, ratio_ref[, c("gene_name", "Visium_mean_CPM", "transcript_length")], by = "gene_name")
ratio_dat <- ratio_dat[is.finite(ratio_dat$Visium_mean_CPM) & ratio_dat$Visium_mean_CPM > 0, , drop = FALSE]
ratio_dat$ratio <- ratio_dat$PhotoSeq_mean_CPM / ratio_dat$Visium_mean_CPM
ratio_dat <- ratio_dat[is.finite(ratio_dat$ratio) & ratio_dat$ratio > 0, , drop = FALSE]
ratio_dat$length_bin <- cut(ratio_dat$transcript_length, c(-Inf, 2000, 3500, 5000, 7000, Inf),
                            labels = c("0-2000", "2000-3500", "3500-5000", "5000-7000", ">7000"), right = FALSE)
ratio_cols <- c("0-2000" = "#F8E4AE", "2000-3500" = "#F2C278", "3500-5000" = "#E99A5B", "5000-7000" = "#D5774F", ">7000" = "#985345")
ratio_dat$length_bin <- factor(ratio_dat$length_bin, levels = names(ratio_cols))
# The ratio is plotted on its positive raw scale with a log10 transformation.
# Thus, ratio = 1 is shown at 10^0; a raw y value of 0 is undefined here.
ratio_y_breaks <- 10^(-2:4)
ratio_y_labels <- parse(text = paste0("10^", -2:4))
ratio_plot <- ggplot(ratio_dat, aes(length_bin, ratio, fill = length_bin)) +
  geom_boxplot(outlier.shape = 18, outlier.size = 1.25, outlier.colour = "#3E3E3E",
               colour = "#555555", linewidth = .55, alpha = .9,
               whisker.linewidth = .55, staple.linewidth = .55, staplewidth = .5) +
  geom_hline(aes(yintercept = 1, linetype = "Ratio = 1"), colour = "#D6272D", linewidth = 1.0) +
  scale_fill_manual(values = ratio_cols, guide = "none") +
  scale_linetype_manual(values = c("Ratio = 1" = "dashed"), name = NULL) +
  scale_y_log10(limits = c(1e-2, 2e4), breaks = ratio_y_breaks, labels = ratio_y_labels, expand = c(0, 0)) +
  labs(title = "Transcript length-dependent PhotoSeq-to-Visium expression ratio",
       x = "Transcript Length (nt)", y = "Mean CPM ratio (PhotoSeq / Visium)") +
  theme_classic(base_size = 12, base_family = "Arial") +
  theme(axis.line = element_line(linewidth = 1.05, colour = "black"), axis.ticks = element_line(linewidth = 1.0, colour = "black"),
        axis.title = element_text(size = 14), axis.text = element_text(size = 12, colour = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        plot.title = element_text(size = 14, hjust = .5, face = "plain"),
        panel.grid.major.y = element_line(colour = "#E5E5E5", linewidth = .5), panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.05),
        legend.position = c(.91, .95), legend.justification = c(1, 1),
        legend.background = element_rect(fill = "white", colour = "#BDBDBD", linewidth = .8),
        legend.key.width = unit(1.2, "cm"), legend.text = element_text(size = 11),
        plot.margin = margin(8, 12, 8, 10))
save_plot(ratio_plot, file.path(qc2_dir, "Breast_Late_PhotoSeq_vs_Visium_ratio_by_transcript_length"), 190, 140)
write.csv(ratio_dat, file.path(qc2_dir, "Breast_Late_ratio_source_data.csv"), row.names = FALSE)

# 03: all-gene expression density; no batch correction.
qc3_dir <- file.path(out_root, "03.PhotoSeq_vs_Visium_density")
dir.create(qc3_dir, recursive = TRUE, showWarnings = FALSE)
density_ref <- read.csv(density_source, stringsAsFactors = FALSE)
density_dat <- merge(photo_late, density_ref[, c("gene_name", "Visium_mean_CPM")], by = "gene_name")
density_dat$PhotoSeq_log2_CPM_plus1 <- log2(density_dat$PhotoSeq_mean_CPM + 1)
density_dat$Visium_log2_CPM_plus1 <- log2(density_dat$Visium_mean_CPM + 1)
density_dat <- density_dat[is.finite(density_dat$PhotoSeq_log2_CPM_plus1) & is.finite(density_dat$Visium_log2_CPM_plus1), , drop = FALSE]
density_r <- cor(density_dat$PhotoSeq_log2_CPM_plus1, density_dat$Visium_log2_CPM_plus1, method = "pearson")
density_s <- cor(density_dat$PhotoSeq_log2_CPM_plus1, density_dat$Visium_log2_CPM_plus1, method = "spearman")
density_n <- nrow(density_dat)
density_hex <- function(d, nx = 55) {
  xr <- range(d$PhotoSeq_log2_CPM_plus1); yr <- range(d$Visium_log2_CPM_plus1)
  dx <- diff(xr) / nx; dy <- dx * sqrt(3) / 2
  ix <- floor((d$PhotoSeq_log2_CPM_plus1 - xr[1]) / dx) + 1
  iy <- floor((d$Visium_log2_CPM_plus1 - yr[1]) / dy) + 1
  tab <- as.data.frame(table(ix = ix, iy = iy), stringsAsFactors = FALSE)
  tab$ix <- as.integer(tab$ix); tab$iy <- as.integer(tab$iy); tab$count <- as.numeric(tab$Freq)
  tab <- tab[tab$count > 0, , drop = FALSE]
  tab$cx <- xr[1] + (tab$ix - .5) * dx
  tab$cy <- yr[1] + (tab$iy - .5) * dy + ifelse(tab$ix %% 2 == 0, dy / 2, 0)
  a <- seq(0, 2 * pi, length.out = 7)
  do.call(rbind, lapply(seq_len(nrow(tab)), function(i) data.frame(hex_id = i, count = tab$count[i], x = tab$cx[i] + .58 * dx * cos(a), y = tab$cy[i] + .58 * dx * sin(a))))
}
density_hex_dat <- density_hex(density_dat)
density_plot <- ggplot(density_dat, aes(PhotoSeq_log2_CPM_plus1, Visium_log2_CPM_plus1)) +
  geom_polygon(data = density_hex_dat, aes(x, y, group = hex_id, fill = count), inherit.aes = FALSE, colour = NA) +
  # Darker teal ramp improves visibility while preserving the log-scaled density range.
  scale_fill_gradientn(colours = c("#D4ECEB", "#9BCDC8", "#5C9E98", "#28777A", "#084656"), trans = "log10", name = "Genes per hexagon\n(log scale)") +
  annotate("label", x = 1.0, y = 13.8, hjust = 0, vjust = 1,
           label = paste0("Pearson R = ", sprintf("%.3f", density_r), "\nSpearman R = ", sprintf("%.3f", density_s), "\nn = ", format(density_n, big.mark = ","), " matched gene pairs"),
           size = 4.4, label.size = .45, fill = "white", colour = "black", label.r = unit(.15, "lines")) +
  labs(title = "Gene-level expression concordance between PhotoSeq and Visium",
       x = "PhotoSeq log2(CPM + 1)", y = "10x Genomics Visium log2(CPM + 1)") +
  coord_cartesian(xlim = c(0, 14.5), ylim = c(0, 14.5), expand = FALSE) +
  theme_classic(base_size = 11, base_family = "Arial") +
  theme(axis.line = element_line(linewidth = .8), axis.ticks = element_line(linewidth = .7),
        axis.title = element_text(size = 13), axis.text = element_text(size = 10, colour = "black"),
        plot.title = element_text(size = 14, hjust = .5, face = "plain"), legend.title = element_text(size = 11),
        legend.text = element_text(size = 9), panel.grid = element_blank(), plot.margin = margin(12, 22, 12, 22))
save_plot(density_plot, file.path(qc3_dir, "Breast_Late_PhotoSeq_vs_Visium_all_gene_expression_density"), 190, 170)
write.csv(density_dat, file.path(qc3_dir, "Breast_Late_density_source_data.csv"), row.names = FALSE)
write.csv(data.frame(
  metric = c("PhotoSeq late gene rows in input matrix", "PhotoSeq late genes with mean CPM > 1", "Visium reference rows", "matched gene pairs used for density", "batch_correction"),
  value = c(nrow(photo_late), sum(late_cpm > 1), nrow(density_ref), density_n, "none"),
  stringsAsFactors = FALSE
), file.path(qc3_dir, "Breast_Late_density_summary.csv"), row.names = FALSE)

# 04: complete numeric Pearson R-squared matrix after integrated ComBat_seq.
qc4_dir <- file.path(out_root, "04.integrated_Pearson_correlation")
dir.create(qc4_dir, recursive = TRUE, showWarnings = FALSE)
group_order <- c("Breast_early", "Breast_tumor", "Breast_late", "Breast_normal_adjacent", "Breast_normal", "Lung_metastasis", "Lung_normal_adjacent", "Lung_normal")
group_labels <- c(Breast_early = "Breast\nearly", Breast_tumor = "Breast mid\n(tumor)", Breast_late = "Breast\nlate", Breast_normal_adjacent = "Breast mid\n(adjacent)", Breast_normal = "Breast mid\n(normal)", Lung_metastasis = "Lung\nmetastasis", Lung_normal_adjacent = "Lung\nnormal adjacent", Lung_normal = "Lung\nnormal")
group_mean <- sapply(group_order, function(g) rowMeans(integrated_vst[, all_group$sample_id[all_group$group == g], drop = FALSE]))
rownames(group_mean) <- rownames(integrated_vst)
pearson <- cor(group_mean, method = "pearson", use = "pairwise.complete.obs")
pearson_r2 <- pearson^2
write.csv(pearson_r2, file.path(qc4_dir, "PhotoSeq_8group_Pearson_R2_matrix.csv"))
cor_dat <- expand.grid(row_group = group_order, col_group = group_order, stringsAsFactors = FALSE)
cor_dat$R2 <- as.vector(pearson_r2[cbind(match(cor_dat$row_group, group_order), match(cor_dat$col_group, group_order))])
cor_dat$row_group <- factor(cor_dat$row_group, levels = group_order); cor_dat$col_group <- factor(cor_dat$col_group, levels = group_order)
cor_dat$row_i <- as.integer(cor_dat$row_group); cor_dat$col_i <- as.integer(cor_dat$col_group)
cor_dat$number_colour <- ifelse(cor_dat$R2 >= .72, "white", "black")
cor_plot <- ggplot(cor_dat, aes(col_group, row_group, fill = R2)) + geom_tile(colour = "white", linewidth = .55) +
  geom_text(aes(label = ifelse(R2 == 1, "1", sprintf("%.3f", R2)), colour = number_colour), size = 3.1, show.legend = FALSE) +
  scale_colour_identity() + scale_fill_gradientn(colours = c("#F2F0FF", "#B7A9EE", "#6154E5", "#0A00E8"), limits = c(0, 1), name = expression(Pearson~R^2)) +
  scale_x_discrete(labels = function(x) unname(group_labels[x])) + scale_y_discrete(labels = function(x) unname(group_labels[x])) +
  labs(title = "Pearson correlation", x = NULL, y = NULL) + coord_fixed() + theme_minimal(base_size = 10) +
  theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, size = 10), axis.text.y = element_text(size = 10),
        plot.title = element_text(hjust = .5, face = "plain", size = 20),
        legend.position = "left", legend.justification = "center", legend.margin = margin(r = 6),
        legend.title = element_text(size = 12), legend.text = element_text(size = 10),
        plot.margin = margin(8, 8, 8, 12))
save_plot(cor_plot, file.path(qc4_dir, "PhotoSeq_8group_Pearson_R2_full_numeric"), 205, 180)
write.csv(cor_dat, file.path(qc4_dir, "PhotoSeq_8group_Pearson_R2_plot_data.csv"), row.names = FALSE)

# 05: DEG union, requiring padj < .01 and absolute log2FC > 1.5.
qc5_dir <- file.path(out_root, "05.DEG_zscore_heatmap"); dir.create(qc5_dir, recursive = TRUE, showWarnings = FALSE)
sig_files <- c(list.files("01_DEGs_Spatial", pattern = "^02\\.Sig_.*\\.csv$", recursive = TRUE, full.names = TRUE), list.files("08_DEGs_Temporal", pattern = "^02\\.Sig_.*\\.csv$", recursive = TRUE, full.names = TRUE))
sig_list <- lapply(sig_files, function(f) read.csv(f, stringsAsFactors = FALSE, check.names = FALSE)); sig_list <- sig_list[vapply(sig_list, nrow, integer(1)) > 0]
# QC5 uses the batch-aware DESeq2 result set (padj < .01, |log2FC| > 1.5)
# but displays a balanced representative panel rather than letting the largest
# Breast-mid comparisons dominate a variance-ranked top-350 list. For each
# comparison, retain the 35 smallest-padj up genes and 35 smallest-padj down
# genes, then deduplicate across comparisons.
sig_genes_by_comparison <- lapply(sig_list, function(z) {
  z <- z[!is.na(z$padj) & z$padj < .01 & is.finite(z$log2FoldChange) & abs(z$log2FoldChange) > 1.5, , drop = FALSE]
  z <- z[!is.na(z$GeneSymbol) & nzchar(z$GeneSymbol), , drop = FALSE]
  up <- z[z$log2FoldChange > 1.5, , drop = FALSE]
  down <- z[z$log2FoldChange < -1.5, , drop = FALSE]
  up <- up[order(up$padj, -abs(up$log2FoldChange)), , drop = FALSE]
  down <- down[order(down$padj, -abs(down$log2FoldChange)), , drop = FALSE]
  unique(c(head(up$GeneSymbol, 35), head(down$GeneSymbol, 35)))
})
sig_genes <- unique(unlist(sig_genes_by_comparison, use.names = FALSE))
sig_genes <- intersect(sig_genes, rownames(integrated_vst))
z <- t(scale(t(group_mean[, group_order, drop = FALSE][sig_genes, , drop = FALSE])))
z <- z[apply(z, 1, function(v) all(is.finite(v))), , drop = FALSE]
z <- pmax(pmin(z, 2), -2)
heatmap_order <- c("Lung_normal", "Lung_metastasis", "Lung_normal_adjacent", "Breast_tumor", "Breast_late", "Breast_normal_adjacent", "Breast_early", "Breast_normal")
z <- z[, heatmap_order, drop = FALSE]
colnames(z) <- unname(group_labels[heatmap_order])
write.csv(z, file.path(qc5_dir, "DEG_zscore_heatmap_matrix.csv"))
write.csv(data.frame(
  metric = c("complete significant DEG union", "per-comparison up genes", "per-comparison down genes", "deduplicated display genes", "display_transform", "expression_matrix"),
  value = c(length(unique(unlist(lapply(sig_list, function(z) z$GeneSymbol[z$padj < .01 & abs(z$log2FoldChange) > 1.5])))), 35, 35, nrow(z), "row z-score, clipped to [-2, 2]", "integrated ComBat_seq -> VST group means"),
  stringsAsFactors = FALSE
), file.path(qc5_dir, "DEG_zscore_heatmap_selection_parameters.csv"), row.names = FALSE)
heatmap_cols <- c(Lung_normal = "#9DB8D1", Lung_metastasis = "#8E2E5A", Lung_normal_adjacent = "#C6E7AF", Breast_tumor = "#9A87C8", Breast_late = "#5E3999", Breast_normal_adjacent = "#E5C8A0", Breast_early = "#D7D2E8", Breast_normal = "#9DB8D1")
hm <- function() Heatmap(z, name = "Z score", col = circlize::colorRamp2(c(-2, 0, 2), c("#2166AC", "white", "#D73027")),
                          cluster_rows = TRUE, cluster_columns = FALSE, show_row_names = FALSE, column_names_rot = 45,
                          column_names_gp = grid::gpar(fontsize = 10), column_title = "Representative DEG patterns across PhotoSeq groups",
                          column_title_gp = grid::gpar(fontsize = 16, fontface = "plain"),
                          top_annotation = HeatmapAnnotation(Group = heatmap_order, col = list(Group = heatmap_cols), show_legend = FALSE,
                                                             annotation_height = grid::unit(4, "mm")),
                          heatmap_legend_param = list(title = "Z score", title_gp = grid::gpar(fontsize = 11), labels_gp = grid::gpar(fontsize = 9)))
pdf(file.path(qc5_dir, "Differential_gene_expression_zscore_heatmap_reference_balanced.pdf"), width = 9, height = 12, useDingbats = FALSE)
if (nrow(z) > 0) print(hm()) else { plot.new(); text(.5, .5, "No genes passed padj < 0.01 and |log2FC| > 1.5") }; dev.off()
png(file.path(qc5_dir, "Differential_gene_expression_zscore_heatmap_reference_balanced.png"), width = 1800, height = 2400, res = 250)
if (nrow(z) > 0) print(hm()) else { plot.new(); text(.5, .5, "No genes passed thresholds") }; dev.off()

# 06: area-proportional detected-gene overlap at mean CPM > 1; no batch correction.
qc6_dir <- file.path(out_root, "06.detected_gene_overlap"); dir.create(qc6_dir, recursive = TRUE, showWarnings = FALSE)
status_ref <- read.csv(overlap_source, stringsAsFactors = FALSE)
gene_col <- intersect(c("gene_name", "Gene", "gene", "gene_id"), names(status_ref))[1]; if (is.na(gene_col)) gene_col <- names(status_ref)[1]
visium_mean_col <- intersect(c("Visium_mean_CPM", "visium_mean_CPM"), names(status_ref))[1]
if (is.na(visium_mean_col)) stop("Visium mean-CPM column is required for QC6")
visium_mean <- suppressWarnings(as.numeric(status_ref[[visium_mean_col]]))
visium_genes <- unique(as.character(status_ref[[gene_col]][is.finite(visium_mean) & visium_mean > 1]))
detected_photo_genes <- names(late_cpm)[late_cpm > 1]
shared_genes <- intersect(detected_photo_genes, visium_genes); union_genes <- union(detected_photo_genes, visium_genes)
photo_n <- length(detected_photo_genes); visium_n <- length(visium_genes); shared_n <- length(shared_genes)
shared_visium_pct <- 100 * shared_n / max(1, visium_n)
jaccard_pct <- 100 * shared_n / max(1, length(union_genes))
overlap_coefficient_pct <- 100 * shared_n / max(1, min(photo_n, visium_n))
write.csv(data.frame(PhotoSeq = photo_n, Visium = visium_n, Shared = shared_n, Union = length(union_genes), Shared_percent_Visium = shared_visium_pct, Jaccard_percent_union = jaccard_pct, Overlap_coefficient = overlap_coefficient_pct), file.path(qc6_dir, "detected_gene_overlap_summary.csv"), row.names = FALSE)
write.csv(data.frame(
  metric = c("PhotoSeq genes with mean CPM > 1", "Detected genes with mean CPM exactly 0", "Detected genes with at least one sample CPM = 0", "Detected genes with all three samples CPM > 0", "Detected genes with all three samples CPM > 1", "Minimum detected mean CPM", "Detection rule"),
  value = c(photo_n, sum(late_cpm[detected_photo_genes] == 0),
            sum(rowSums(late_cpm_matrix[detected_photo_genes, , drop = FALSE] == 0) > 0),
            sum(rowSums(late_cpm_matrix[detected_photo_genes, , drop = FALSE] > 0) == ncol(late_cpm_matrix)),
            sum(rowSums(late_cpm_matrix[detected_photo_genes, , drop = FALSE] > 1) == ncol(late_cpm_matrix)),
            sprintf("%.8f", min(late_cpm[detected_photo_genes])),
            "mean CPM across Breast_late_N1/N2/N3 > 1"),
  stringsAsFactors = FALSE
), file.path(qc6_dir, "photo_seq_late_detection_summary.csv"), row.names = FALSE)

# Solve the circle separation so the geometric intersection is proportional to shared_n,
# while each circle area remains proportional to its detected-gene count.
circle_intersection <- function(r1, r2, d) {
  if (d >= r1 + r2) return(0)
  if (d <= abs(r1 - r2)) return(pi * min(r1, r2)^2)
  a1 <- acos((d^2 + r1^2 - r2^2) / (2 * d * r1))
  a2 <- acos((d^2 + r2^2 - r1^2) / (2 * d * r2))
  .5 * r1^2 * (2 * a1 - sin(2 * a1)) + .5 * r2^2 * (2 * a2 - sin(2 * a2))
}
plot_scale <- .35 / sqrt(pi)
r_photo <- sqrt(photo_n) * plot_scale; r_visium <- sqrt(visium_n) * plot_scale
target_intersection <- shared_n * pi * plot_scale^2
d <- uniroot(function(x) circle_intersection(r_photo, r_visium, x) - target_intersection,
             lower = abs(r_photo - r_visium) + 1e-8, upper = r_photo + r_visium - 1e-8)$root
circle_data <- function(cx, cy, r, name) { a <- seq(0, 2 * pi, length.out = 721); data.frame(x = cx + r * cos(a), y = cy + r * sin(a), set = name) }
venn_shape <- rbind(circle_data(-d / 2, 0, r_photo, "PhotoSeq"), circle_data(d / 2, 0, r_visium, "10x Genomics Visium"))
xmax <- max(abs(venn_shape$x)) + 18; ymax <- max(r_photo, r_visium) + 20
venn_plot <- ggplot() +
  geom_polygon(data = venn_shape, aes(x, y, group = set, fill = set), colour = "#777777", linewidth = .75, alpha = .62) +
  scale_fill_manual(values = c("PhotoSeq" = "#9BBBD6", "10x Genomics Visium" = "#F3CE91"), guide = "none") +
  annotate("text", x = -d / 2, y = r_photo + 8, label = paste0("PhotoSeq\nn = ", format(photo_n, big.mark = ","), " detected"), colour = "#35668F", size = 4.5, fontface = "bold", family = "Arial", lineheight = .95) +
  annotate("text", x = d / 2, y = r_visium + 8, label = paste0("10x Genomics Visium\nn = ", format(visium_n, big.mark = ","), " detected"), colour = "#A6661F", size = 4.5, fontface = "bold", family = "Arial", lineheight = .95) +
  annotate("text", x = 0, y = .10 * min(r_photo, r_visium), label = paste0(format(shared_n, big.mark = ","), "\n(", sprintf("%.2f", shared_visium_pct), "% of Visium)"), size = 5.4, fontface = "bold", family = "Arial", lineheight = .95) +
  annotate("text", x = 0, y = -ymax + 8, label = paste0("Detected genes: mean CPM > 1; Jaccard = ", sprintf("%.2f", jaccard_pct), "%"), size = 4.2, colour = "#4B4B4B", family = "Arial") +
  labs(title = "Overlap of genes detected by PhotoSeq and Visium") +
  coord_fixed(xlim = c(-xmax, xmax), ylim = c(-ymax, ymax), expand = FALSE, clip = "off") +
  theme_void(base_family = "Arial") + theme(plot.title = element_text(size = 18, hjust = .5, face = "bold", margin = margin(b = 8)), plot.margin = margin(10, 20, 10, 20))
save_plot(venn_plot, file.path(qc6_dir, "Supplementary_detected_gene_overlap_CPM_gt1"), 205, 150)

# 07: UV condition optimization using the supplied fluorescence table;
# no batch correction or significance testing is applied. The figure design
# uses mean bars, mean +/- SD, three black replicate points, and a fixed
# blue-grey -> pale-red condition palette. Only the Y-axis quantity changes
# from detected gene count to mean fluorescence intensity (a.u.).
qc7_dir <- file.path(out_root, "07.UV_condition_optimization")
dir.create(qc7_dir, recursive = TRUE, showWarnings = FALSE)
uv_raw <- as.data.frame(readxl::read_excel(uv_source, sheet = 1, .name_repair = "minimal"), stringsAsFactors = FALSE)
required_uv_cols <- c("panel", "condition", "replicate", "Mean fluorescence intensity (a.u.)")
if (!all(required_uv_cols %in% names(uv_raw))) {
  stop("UV workbook is missing required columns: ", paste(setdiff(required_uv_cols, names(uv_raw)), collapse = ", "))
}
uv_dat <- uv_raw[, required_uv_cols, drop = FALSE]
names(uv_dat) <- c("panel", "condition", "replicate", "fluorescence_intensity")
uv_dat$panel <- as.character(uv_dat$panel)
uv_dat$condition <- as.character(uv_dat$condition)
uv_dat$replicate <- suppressWarnings(as.integer(uv_dat$replicate))
uv_dat$fluorescence_intensity <- suppressWarnings(as.numeric(uv_dat$fluorescence_intensity))
if (any(!is.finite(uv_dat$fluorescence_intensity)) || any(is.na(uv_dat$replicate))) stop("UV workbook contains non-numeric fluorescence or replicate values")
if (anyDuplicated(uv_dat[, c("panel", "condition", "replicate")])) stop("UV workbook contains duplicate panel-condition-replicate rows")

uv_specs <- list(
  uv_power = list(
    order = c("100% UV", "80% UV", "60% UV", "40% UV", "20% UV"),
    xlabel = "Relative UV power (%, 100% = 20 mW)",
    colors = c("#91A6C7", "#AFCDE8", "#C9D6E8", "#F5E6DF", "#F1978F"),
    stem = "UV_power_condition_optimization_fluorescence_intensity"
  ),
  exposure_time = list(
    order = c("5 s", "10 s", "20 s", "30 s", "40 s", "50 s"),
    xlabel = "UV Exposure Time at 40% UV power (s)",
    colors = c("#91A6C7", "#AFCDE8", "#C9D6E8", "#D7D8EC", "#F5E6DF", "#F1978F"),
    stem = "UV_exposure_time_condition_optimization_fluorescence_intensity"
  )
)
if (!all(names(uv_specs) %in% unique(uv_dat$panel))) stop("UV workbook is missing one or more expected panels")

summarise_uv <- function(d, condition_order) {
  out <- do.call(rbind, lapply(condition_order, function(label) {
    v <- d$fluorescence_intensity[d$condition == label]
    if (!length(v)) stop("UV workbook is missing condition: ", label)
    data.frame(condition = label, mean = mean(v), sd = if (length(v) > 1) sd(v) else 0, n = length(v), stringsAsFactors = FALSE)
  }))
  out$condition <- factor(out$condition, levels = condition_order)
  out
}

make_uv_optimization_plot <- function(panel_name, spec) {
  d <- uv_dat[uv_dat$panel == panel_name, , drop = FALSE]
  if (!nrow(d)) stop("UV workbook contains no rows for panel: ", panel_name)
  d$condition <- factor(d$condition, levels = spec$order)
  sm <- summarise_uv(d, spec$order)
  upper <- max(sm$mean + sm$sd, na.rm = TRUE)
  y_limit <- ceiling(upper / 10) * 10 + 10
  y_ticks <- seq(0, y_limit - 10, by = 20)
  p <- ggplot(sm, aes(condition, mean, fill = condition)) +
    geom_col(width = .56, colour = "#4A4A4A", linewidth = .55) +
    geom_errorbar(aes(ymin = pmax(0, mean - sd), ymax = mean + sd), width = .22, linewidth = .65, colour = "black") +
    geom_point(data = d, aes(condition, fluorescence_intensity), inherit.aes = FALSE,
               position = position_jitter(width = .055, height = 0, seed = 11), shape = 16, size = 2.25, colour = "black") +
    scale_fill_manual(values = setNames(spec$colors, spec$order), guide = "none") +
    scale_y_continuous(limits = c(0, y_limit), breaks = y_ticks, expand = c(0, 0)) +
    labs(x = spec$xlabel, y = "Mean fluorescence intensity (a.u.)") +
    theme_classic(base_size = 10, base_family = "Arial") +
    theme(axis.line = element_line(linewidth = .75, colour = "#4A4A4A"),
          axis.ticks = element_line(linewidth = .75, colour = "black"),
          axis.title = element_text(size = 10.5, colour = "black"),
          axis.text = element_text(size = 9.2, colour = "black"),
          panel.grid = element_blank(),
          panel.border = element_rect(colour = "#4A4A4A", fill = NA, linewidth = .75),
          plot.margin = margin(5, 6, 5, 5))
  save_plot(p, file.path(qc7_dir, spec$stem), 144, 108)
  write.csv(d[order(d$condition, d$replicate), c("panel", "condition", "replicate", "fluorescence_intensity")],
            file.path(qc7_dir, paste0(spec$stem, "_source_data.csv")), row.names = FALSE)
  write.csv(data.frame(panel = panel_name, condition = as.character(sm$condition), mean_fluorescence_intensity = sm$mean,
                       sd_fluorescence_intensity = sm$sd, n = sm$n, y_axis_limit = y_limit,
                       y_axis = "Mean fluorescence intensity (a.u.)", batch_correction = "none",
                       stringsAsFactors = FALSE),
            file.path(qc7_dir, paste0(spec$stem, "_summary.csv")), row.names = FALSE)
}
make_uv_optimization_plot("uv_power", uv_specs$uv_power)
make_uv_optimization_plot("exposure_time", uv_specs$exposure_time)
write.csv(data.frame(
  parameter = c("source_workbook", "source_sheet", "Y_axis", "summary", "replicates_per_condition", "batch_correction", "padj", "image_formats"),
  value = c(uv_source, "uv_pair_source_data", "Mean fluorescence intensity (a.u.)", "mean with SD and three replicate points", "3", "none", "not applicable", "PNG and PDF only"),
  stringsAsFactors = FALSE
), file.path(qc7_dir, "UV_condition_optimization_parameters.csv"), row.names = FALSE)

message("17_quality_control_and_benchmarking completed: ", out_root)
