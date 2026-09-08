rm(list = ls()); gc()
ORIGINAL_DIR <- "."
output <- file.path(ORIGINAL_DIR, "12_Venn_All")
if (!dir.exists(output)) {dir.create(output, recursive = TRUE)}
setwd(ORIGINAL_DIR)

source(file.path(ORIGINAL_DIR, "scripts", "area_proportional_venn_helpers.R"))

genes1 <- read.csv("02_Venn_Spatial/02.Venn_DEGs.csv", stringsAsFactors = FALSE)
genes2 <- read.csv("09_Venn_Temporal/02.Venn_DEGs.csv", stringsAsFactors = FALSE)

sets <- list(
  DEG1 = unique(genes1$x),
  DEG2 = unique(genes2$x)
)

p <- plot_area_venn2(
  sets,
  colors = c("#94be98", "#fcbe6e"),
  title = NULL
)

save_venn_plot(p, file.path(output, "01.Venn_DEGs"), width = 6, height = 6)

genes <- intersect(sets[[1]], sets[[2]])
write.csv(data.frame(x = genes), file.path(output, "02.Venn_DEGs.csv"), quote = FALSE, row.names = FALSE)

summary_df <- data.frame(
  set = c("DEG1", "DEG2", "DEG1_and_DEG2"),
  count = c(length(sets[[1]]), length(sets[[2]]), length(genes))
)
write.csv(summary_df, file.path(output, "01.Venn_DEGs_area_summary.csv"), quote = FALSE, row.names = FALSE)


# -------------------- 03 shared-gene horizontal heatmap --------------------
# Landscape layout for the shared-gene heatmap.
# The analysis uses the shared-gene set, integrated normalized matrix,
# row-wise z-score, clipping to [-2, 2], fixed sample order, gene clustering
# and group annotation.
# Backend: R / ComplexHeatmap only.

suppressPackageStartupMessages({
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
})

root <- "."
out_dir <- file.path(root, "12_Venn_All", "shared_gene_heatmap_horizontal")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

genes <- read.csv(file.path(root, "12_Venn_All", "02.Venn_DEGs.csv"),
                  stringsAsFactors = FALSE)[["x"]]
genes <- unique(genes)
stopifnot(length(genes) > 0L)

load(file.path(root, "00_prepared_data", "03.All_dat_norm.RData"))
group <- read.csv(file.path(root, "00_prepared_data", "03.All_group.csv"),
                  stringsAsFactors = FALSE, check.names = FALSE)
group_levels <- c(
  "Breast_early", "Breast_late", "Breast_tumor",
  "Breast_normal_adjacent", "Breast_normal",
  "Lung_metastasis", "Lung_normal_adjacent", "Lung_normal"
)
group$group <- factor(group$group, levels = group_levels)
group <- group[order(group$group, group$sample_id), , drop = FALSE]

missing_genes <- setdiff(genes, rownames(dat_norm))
if (length(missing_genes) > 0) {
  stop("Shared genes missing from normalized matrix: ",
       paste(missing_genes, collapse = ", "))
}
expr <- as.matrix(dat_norm[genes, group$sample_id, drop = FALSE])
row_zscore <- t(apply(expr, 1, function(x) {
  s <- sd(x, na.rm = TRUE)
  if (!is.finite(s) || s == 0) rep(0, length(x)) else
    (x - mean(x, na.rm = TRUE)) / s
}))
rownames(row_zscore) <- rownames(expr)
colnames(row_zscore) <- colnames(expr)
row_zscore[!is.finite(row_zscore)] <- 0
row_zscore <- pmax(pmin(row_zscore, 2), -2)

source_data <- data.frame(GeneSymbol = rownames(row_zscore),
                          row_zscore, check.names = FALSE)
write.csv(source_data,
          file.path(out_dir, "shared_gene_heatmap_horizontal_source_data.csv"),
          row.names = FALSE, quote = FALSE)

group_labels <- gsub("_", " ", as.character(group$group))
group_colors <- c(
  "Breast early" = "#984EA3",
  "Breast late" = "#1B9E77",
  "Breast tumor" = "#33A02C",
  "Breast normal adjacent" = "#FF7F00",
  "Breast normal" = "#1F78B4",
  "Lung metastasis" = "#E31A1C",
  "Lung normal adjacent" = "#FDBF6F",
  "Lung normal" = "#A6CEE3"
)
group_annotation <- HeatmapAnnotation(
  Group = group_labels,
  col = list(Group = group_colors),
  annotation_legend_param = list(
    Group = list(
      title_gp = gpar(fontsize = 17, fontface = "bold"),
      labels_gp = gpar(fontsize = 16),
      grid_height = unit(8, "mm"),
      grid_width = unit(8, "mm")
    )
  ),
  show_annotation_name = TRUE,
  annotation_name_gp = gpar(fontsize = 12, fontface = "bold"),
  simple_anno_size = unit(6, "mm")
)
expression_colors <- colorRamp2(c(-2, 0, 2), c("#2166AC", "white", "#B2182B"))

heatmap_density <- densityHeatmap(
  row_zscore,
  title = "Distribution as heatmap",
  ylab = " ",
  height = unit(3.2, "cm"),
  show_column_names = FALSE,
  title_gp = gpar(fontsize = 16),
  ylab_gp = gpar(fontsize = 10),
  tick_label_gp = gpar(fontsize = 9),
  quantile_gp = gpar(fontsize = 9),
  heatmap_legend_param = list(
    title_gp = gpar(fontsize = 17, fontface = "bold"),
    labels_gp = gpar(fontsize = 16),
    grid_height = unit(8, "mm"),
    grid_width = unit(8, "mm")
  )
)

heatmap_expression <- Heatmap(
  row_zscore,
  name = "expression",
  col = expression_colors,
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  show_column_names = FALSE,
  show_row_names = TRUE,
  row_names_side = "right",
  row_names_gp = gpar(fontsize = 14, fontface = "italic"),
  column_title = NULL,
  heatmap_legend_param = list(
    at = c(-2, -1, 0, 1, 2),
    title = "expression",
    labels_gp = gpar(fontsize = 16),
    title_gp = gpar(fontsize = 17, fontface = "bold"),
    grid_height = unit(8, "mm"),
    grid_width = unit(8, "mm")
  ),
  border = FALSE,
  rect_gp = gpar(col = NA)
)

heatmap_all <- heatmap_density %v% group_annotation %v% heatmap_expression

draw_landscape <- function() {
  grid.newpage()
  draw(
    heatmap_all,
    column_title = paste0(length(genes), " shared spatiotemporal DEGs"),
    column_title_gp = gpar(fontsize = 28, fontface = "bold"),
    merge_legends = TRUE,
    heatmap_legend_side = "right",
    annotation_legend_side = "right",
    padding = unit(c(12, 12, 12, 12), "mm")
  )
}

base <- file.path(out_dir, "shared_gene_heatmap_horizontal")
cairo_pdf(paste0(base, ".pdf"), width = 18, height = 13.5,
          family = "Arial", bg = "white")
draw_landscape()
dev.off()

png(paste0(base, ".png"), width = 18, height = 13.5,
    units = "in", res = 450, type = "cairo", bg = "white")
draw_landscape()
dev.off()

message("Wrote landscape ", length(genes), "-shared-gene heatmap to ", out_dir)



# -------------------- quadrant scatter --------------------
suppressPackageStartupMessages({
  library(ggplot2)
  library(ggrepel)
})

root <- "."
out_dir <- file.path(root, "12_Venn_All", "quadrant_scatter")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

shared_genes <- read.csv(
  file.path(root, "12_Venn_All", "02.Venn_DEGs.csv"),
  stringsAsFactors = FALSE
)[["x"]]

read_log2fc <- function(path) {
  df <- read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  gene_col <- if ("GeneSymbol" %in% names(df)) "GeneSymbol" else names(df)[1]
  lfc_col <- if ("log2FoldChange" %in% names(df)) {
    "log2FoldChange"
  } else if ("log2FC" %in% names(df)) {
    "log2FC"
  } else {
    names(df)[grepl("log2", names(df), ignore.case = TRUE)][1]
  }
  p_col <- if ("pvalue" %in% names(df)) "pvalue" else names(df)[grepl("^p", names(df), ignore.case = TRUE)][1]
  padj_col <- if ("padj" %in% names(df)) "padj" else p_col
  out <- df[, c(gene_col, lfc_col, p_col, padj_col)]
  names(out) <- c("GeneSymbol", "log2FC", "pvalue", "padj")
  out
}

spatial <- read_log2fc(file.path(
  root, "01_DEGs_Spatial", "Lung_metastasis_vs_Breast_tumor",
  "01.All_Lung_metastasis_vs_Breast_tumor.csv"
))
spatial <- spatial[spatial$GeneSymbol %in% shared_genes, ]
names(spatial)[names(spatial) == "log2FC"] <- "spatial_log2FC"
names(spatial)[names(spatial) == "pvalue"] <- "spatial_pvalue"
names(spatial)[names(spatial) == "padj"] <- "spatial_padj"

comparisons <- list(
  list(
    id = "05.Scatter_Spatial_Breast_mid_vs_Breast_early",
    title = "Breast mid vs Breast early",
    xlabel = "Breast mid vs Breast early (log2FC)",
    path = file.path(root, "08_DEGs_Temporal", "Breast_mid_vs_Breast_early", "01.All_Breast_mid_vs_Breast_early.csv")
  ),
  list(
    id = "06.Scatter_Spatial_Breast_late_vs_Breast_mid",
    title = "Breast late vs Breast mid",
    xlabel = "Breast late vs Breast mid (log2FC)",
    path = file.path(root, "08_DEGs_Temporal", "Breast_late_vs_Breast_mid", "01.All_Breast_late_vs_Breast_mid.csv")
  ),
  list(
    id = "07.Scatter_Spatial_Lung_metastasis_vs_Breast_late",
    title = "Lung metastasis vs Breast late",
    xlabel = "Lung metastasis vs Breast late (log2FC)",
    path = file.path(root, "08_DEGs_Temporal", "Lung_metastasis_vs_Breast_late", "01.All_Lung_metastasis_vs_Breast_late.csv")
  )
)

make_plot <- function(comp) {
  temporal <- read_log2fc(comp$path)
  temporal <- temporal[temporal$GeneSymbol %in% shared_genes, ]
  names(temporal)[names(temporal) == "log2FC"] <- "temporal_log2FC"
  names(temporal)[names(temporal) == "pvalue"] <- "temporal_pvalue"
  names(temporal)[names(temporal) == "padj"] <- "temporal_padj"
  dat <- merge(spatial, temporal, by = "GeneSymbol", all = FALSE)
  dat <- dat[match(shared_genes, dat$GeneSymbol), ]
  if (nrow(dat) != length(unique(shared_genes))) {
    stop(comp$title, ": the quadrant plot lost shared genes during merging (",
         nrow(dat), " of ", length(unique(shared_genes)), ").")
  }
  dat$regulation_class <- ifelse(
    dat$spatial_log2FC > 0 & dat$temporal_log2FC > 0,
    "Both up",
    ifelse(
      dat$spatial_log2FC < 0 & dat$temporal_log2FC < 0,
      "Both down",
      "Discordant"
    )
  )
  dat$x_p_score <- -log10(pmax(ifelse(is.na(dat$temporal_padj), 1, dat$temporal_padj), .Machine$double.xmin))
  dat$y_p_score <- -log10(pmax(ifelse(is.na(dat$spatial_padj), 1, dat$spatial_padj), .Machine$double.xmin))
  dat$point_size <- scales::rescale(dat$x_p_score, to = c(1.7, 4.4))
  dat$point_alpha <- scales::rescale(dat$y_p_score, to = c(0.25, 0.95))
  dat$point_color <- ifelse(
    dat$regulation_class == "Both up", "#d62728",
    ifelse(dat$regulation_class == "Both down", "#1f4e79", "#9c6b43")
  )

  # Data-derived axis limits keep every shared gene visible.
  axis_limits <- function(x) {
    r <- range(x, finite = TRUE)
    pad <- max(0.25, diff(r) * 0.05)
    c(floor(r[1] - pad), ceiling(r[2] + pad))
  }
  x_limits <- axis_limits(dat$temporal_log2FC)
  y_limits <- axis_limits(dat$spatial_log2FC)

  p <- ggplot(dat, aes(x = temporal_log2FC, y = spatial_log2FC)) +
    geom_hline(yintercept = 0, linewidth = 0.75, color = "#8f8f8f") +
    geom_vline(xintercept = 0, linewidth = 0.75, color = "#8f8f8f") +
    geom_point(
      aes(size = point_size, color = point_color, alpha = point_alpha),
      stroke = 0.55
    ) +
    geom_text_repel(
      aes(label = GeneSymbol),
      family = "sans",
      size = 2.75,
      fontface = "italic",
      color = "black",
      box.padding = 0.20,
      point.padding = 0.14,
      min.segment.length = 0,
      segment.size = 0.22,
      segment.alpha = 0.75,
      max.overlaps = Inf,
      max.iter = 20000,
      force = 1.8,
      force_pull = 0.08,
      seed = 7
    ) +
    scale_color_identity() +
    scale_alpha_identity() +
    scale_size_identity() +
    labs(
      title = comp$title,
      subtitle = paste0("n = ", nrow(dat), " shared genes"),
      x = comp$xlabel,
      y = "Lung metastasis vs Breast tumor (log2FC)"
    ) +
    coord_cartesian(clip = "off") +
    theme_bw(base_family = "sans", base_size = 10) +
    theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5, color = "black"),
      plot.subtitle = element_text(size = 10.5, hjust = 0.5, color = "#444444"),
      axis.title = element_text(size = 11.5, face = "bold", color = "black"),
      axis.text = element_text(size = 10.5, face = "bold", color = "black"),
      panel.grid.major = element_line(color = "#e5e5e5", linewidth = 0.45),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.75),
      plot.margin = margin(10, 14, 10, 12)
    )

  write.csv(dat, file.path(out_dir, paste0(comp$id, "_source_data.csv")), quote = FALSE, row.names = FALSE)
  p <- p +
    scale_x_continuous(
      limits = x_limits,
      breaks = pretty(x_limits, n = 7),
      expand = expansion(mult = c(0, 0))
    ) +
    scale_y_continuous(
      limits = y_limits,
      breaks = pretty(y_limits, n = 7),
      expand = expansion(mult = c(0, 0))
    )
  ggsave(file.path(out_dir, paste0(comp$id, ".pdf")), p, width = 6.8, height = 6.2, bg = "white", device = cairo_pdf)
  ggsave(file.path(out_dir, paste0(comp$id, ".png")), p, width = 6.8, height = 6.2, dpi = 600, bg = "white")
  p
}

plots <- lapply(comparisons, make_plot)

summary_df <- do.call(rbind, lapply(comparisons, function(comp) {
  dat <- read.csv(file.path(out_dir, paste0(comp$id, "_source_data.csv")), stringsAsFactors = FALSE)
  axis_limits <- function(x) {
    r <- range(x, finite = TRUE)
    pad <- max(0.25, diff(r) * 0.05)
    c(floor(r[1] - pad), ceiling(r[2] + pad))
  }
  x_limits <- axis_limits(dat$temporal_log2FC)
  y_limits <- axis_limits(dat$spatial_log2FC)
  data.frame(
      comparison = comp$title,
      n_shared_genes = nrow(dat),
      x_min = min(dat$temporal_log2FC, na.rm = TRUE),
      x_max = max(dat$temporal_log2FC, na.rm = TRUE),
      y_min = min(dat$spatial_log2FC, na.rm = TRUE),
      y_max = max(dat$spatial_log2FC, na.rm = TRUE),
      x_axis_min = x_limits[1],
      x_axis_max = x_limits[2],
      y_axis_min = y_limits[1],
      y_axis_max = y_limits[2]
  )
}))
write.csv(summary_df, file.path(out_dir, "quadrant_scatter_axis_summary.csv"), quote = FALSE, row.names = FALSE)



# -------------------- state-transition summary --------------------
# The shared genes from the spatial and temporal DEG sets
# occupy distinct quadrant states that change across three comparisons.
# The three quadrant source-data tables contain each gene
# retained and assigned to a mathematical quadrant from the two log2FC signs.
# The alluvial summary reports the corresponding state counts.

suppressPackageStartupMessages(library(grid))

root <- "."
out_dir <- file.path(root, "12_Venn_All", "state_transition_50_shared_genes")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

shared_genes <- read.csv(
  file.path(root, "12_Venn_All", "02.Venn_DEGs.csv"),
  stringsAsFactors = FALSE
)[["x"]]
shared_genes <- unique(shared_genes)
stopifnot(length(shared_genes) > 0L)
shared_count <- length(shared_genes)
state_tag <- paste0("state_transition_", shared_count, "_shared_genes")

source_files <- c(
  file.path(root, "12_Venn_All", "quadrant_scatter",
            "05.Scatter_Spatial_Breast_mid_vs_Breast_early_source_data.csv"),
  file.path(root, "12_Venn_All", "quadrant_scatter",
            "06.Scatter_Spatial_Breast_late_vs_Breast_mid_source_data.csv"),
  file.path(root, "12_Venn_All", "quadrant_scatter",
            "07.Scatter_Spatial_Lung_metastasis_vs_Breast_late_source_data.csv")
)
comparison_names <- c(
  "Breast mid vs Breast early",
  "Breast late vs Breast mid",
  "Lung metastasis vs Breast late"
)

quadrant <- function(x, y) {
  ifelse(x > 0 & y > 0, "Q1",
    ifelse(x < 0 & y > 0, "Q2",
      ifelse(x < 0 & y < 0, "Q3", "Q4")))
}

source_data <- lapply(source_files, function(path) {
  dat <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  dat <- dat[match(shared_genes, dat$GeneSymbol), , drop = FALSE]
  stopifnot(nrow(dat) == length(shared_genes), !anyNA(dat$GeneSymbol))
  dat
})

states <- do.call(cbind, lapply(source_data, function(dat) {
  quadrant(as.numeric(dat$temporal_log2FC), as.numeric(dat$spatial_log2FC))
}))
colnames(states) <- c("state_1", "state_2", "state_3")
gene_states <- data.frame(
  GeneSymbol = shared_genes,
  states,
  stringsAsFactors = FALSE,
  check.names = FALSE
)
write.csv(gene_states,
          file.path(out_dir, paste0(state_tag, "_gene_states.csv")),
          row.names = FALSE, quote = FALSE)

state_defs <- c(
  Q1 = "concordant upregulation (x>0, y>0)",
  Q2 = "discordant: temporal down, spatial up (x<0, y>0)",
  Q3 = "concordant downregulation (x<0, y<0)",
  Q4 = "discordant: temporal up, spatial down (x>0, y<0)"
)
state_order <- names(state_defs)

state_counts <- do.call(rbind, lapply(seq_along(comparison_names), function(i) {
  tab <- table(factor(gene_states[[paste0("state_", i)]], levels = state_order))
  data.frame(comparison = comparison_names[i], state = state_order,
             definition = unname(state_defs), n = as.integer(tab),
             stringsAsFactors = FALSE)
}))
write.csv(state_counts,
          file.path(out_dir, paste0(state_tag, "_state_counts.csv")),
          row.names = FALSE, quote = FALSE)

transition_counts <- do.call(rbind, lapply(1:2, function(i) {
  tab <- table(
    factor(gene_states[[paste0("state_", i)]], levels = state_order),
    factor(gene_states[[paste0("state_", i + 1)]], levels = state_order)
  )
  out <- as.data.frame(tab, stringsAsFactors = FALSE)
  names(out) <- c("from_state", "to_state", "n")
  out$from_comparison <- comparison_names[i]
  out$to_comparison <- comparison_names[i + 1]
  out[, c("from_comparison", "to_comparison", "from_state", "to_state", "n")]
}))
write.csv(transition_counts,
          file.path(out_dir, paste0(state_tag, "_transition_counts.csv")),
          row.names = FALSE, quote = FALSE)

trajectory_key <- apply(gene_states[, c("state_1", "state_2", "state_3")],
                        1, paste, collapse = " ")
trajectory_counts <- as.data.frame(table(trajectory_key), stringsAsFactors = FALSE)
names(trajectory_counts) <- c("trajectory", "n")
trajectory_counts$genes <- vapply(trajectory_counts$trajectory, function(path) {
  ids <- gene_states$GeneSymbol[trajectory_key == path]
  paste(ids, collapse = ";")
}, character(1))
write.csv(trajectory_counts,
          file.path(out_dir, paste0(state_tag, "_trajectory_counts.csv")),
          row.names = FALSE, quote = FALSE)

green <- "#2c8b1f"
green_fill <- "#eef8e9"
orange <- "#f26400"
orange_fill <- "#fff2df"
blue <- "#1f4e79"
blue_fill <- "#e9f1fb"
purple <- "#7b3f98"
purple_fill <- "#f3eafa"
black <- "#111111"
state_col <- c(Q1 = green, Q2 = orange, Q3 = blue, Q4 = purple)
state_fill <- c(Q1 = green_fill, Q2 = orange_fill, Q3 = blue_fill, Q4 = purple_fill)
font_main <- "Arial"

push_fig_viewport <- function() {
  grid.newpage()
  pushViewport(viewport(width = unit(1, "npc"), height = unit(1, "npc")))
}

txt <- function(label, x, y, size = 12, col = black, face = "plain",
                just = "centre", lineheight = 1.05, rot = 0) {
  grid.text(label, x = unit(x, "npc"), y = unit(y, "npc"), just = just, rot = rot,
            gp = gpar(fontfamily = font_main, fontsize = size, col = col,
                      fontface = face, lineheight = lineheight))
}

box <- function(x, y, w, h, state, n) {
  grid.roundrect(x = unit(x, "npc"), y = unit(y, "npc"),
                 width = unit(w, "npc"), height = unit(h, "npc"),
                 r = unit(0.020, "snpc"),
                 gp = gpar(col = state_col[[state]], fill = state_fill[[state]], lwd = 1.8))
  txt(paste0(state, "\n(", n, ")"), x, y, size = 13.5,
      col = state_col[[state]], face = "bold")
}

layout_boxes <- function(counts, y_top = 0.755, y_bottom = 0.245,
                         gap = 0.014, min_h = 0.065) {
  avail <- y_top - y_bottom - gap * (length(counts) - 1)
  extra <- max(0, avail - min_h * length(counts))
  heights <- min_h + extra * as.numeric(counts) / sum(counts)
  centers <- numeric(length(counts))
  cur <- y_top
  for (i in seq_along(counts)) {
    centers[i] <- cur - heights[i] / 2
    cur <- cur - heights[i] - gap
  }
  data.frame(state = state_order, n = as.integer(counts), y = centers,
             h = heights, stringsAsFactors = FALSE)
}

draw_flow <- function(x1, y1, x2, y2, col, n) {
  if (n <= 0) return(invisible(NULL))
  mid <- (x1 + x2) / 2
  grid.bezier(x = unit(c(x1, x1 + (x2 - x1) * 0.35, mid + (x2 - x1) * 0.15, x2), "npc"),
              y = unit(c(y1, y1, y2, y2), "npc"),
              gp = gpar(col = adjustcolor(col, alpha.f = 0.28),
                        lwd = max(1.8, 3.2 * n), lineend = "round"))
}

draw_panel <- function() {
  push_fig_viewport()
  txt("a", 0.030, 0.965, size = 22, face = "bold")
  txt(paste0("State-transition summary of ", length(shared_genes), " shared genes"), 0.525, 0.962,
      size = 25, face = "bold")
  txt("Quadrant trajectories across three pairwise comparisons", 0.525, 0.918,
      size = 16, face = "italic")

  xs <- c(0.145, 0.500, 0.855)
  txt("Breast mid vs\nBreast early", xs[1], 0.842, size = 14.5, face = "bold")
  txt("Breast late vs\nBreast mid", xs[2], 0.842, size = 14.5, face = "bold")
  txt("Lung metastasis vs\nBreast late", xs[3], 0.842, size = 14.5, face = "bold")

  layouts <- lapply(seq_along(comparison_names), function(i) {
    counts <- state_counts$n[state_counts$comparison == comparison_names[i]]
    layout_boxes(counts)
  })

  for (i in 1:2) {
    mat <- matrix(0, nrow = 4, ncol = 4,
                  dimnames = list(state_order, state_order))
    tr <- transition_counts[transition_counts$from_comparison == comparison_names[i], ]
    for (j in seq_len(nrow(tr))) mat[tr$from_state[j], tr$to_state[j]] <- tr$n[j]
    src_pos <- layouts[[i]]$y
    dst_pos <- layouts[[i + 1]]$y
    src_cursor <- src_pos - layouts[[i]]$h / 2
    dst_cursor <- dst_pos - layouts[[i + 1]]$h / 2
    src_slots <- matrix(NA_real_, 4, 4)
    dst_slots <- matrix(NA_real_, 4, 4)
    for (s in seq_len(4)) {
      if (sum(mat[s, ]) > 0) {
        src_slots[s, ] <- src_cursor[s] + (cumsum(mat[s, ]) - mat[s, ] / 2) /
          sum(mat[s, ]) * layouts[[i]]$h[s]
      }
    }
    for (t in seq_len(4)) {
      if (sum(mat[, t]) > 0) {
        dst_slots[, t] <- dst_cursor[t] + (cumsum(mat[, t]) - mat[, t] / 2) /
          sum(mat[, t]) * layouts[[i + 1]]$h[t]
      }
    }
    for (s in seq_len(4)) for (t in seq_len(4)) if (mat[s, t] > 0) {
      draw_flow(xs[i] + 0.055, src_slots[s, t],
                xs[i + 1] - 0.055, dst_slots[s, t],
                state_col[state_order[s]], mat[s, t])
    }
  }

  for (i in seq_along(layouts)) {
    for (j in seq_len(nrow(layouts[[i]]))) {
      box(xs[i], layouts[[i]]$y[j], 0.110, layouts[[i]]$h[j],
          layouts[[i]]$state[j], layouts[[i]]$n[j])
    }
  }

  top <- trajectory_counts[order(-trajectory_counts$n, trajectory_counts$trajectory), ]
  top <- head(top, 5)
  top_label <- paste(paste0(top$trajectory, " (", top$n, ")"), collapse = "   ")
  txt("Most frequent trajectories (n):", 0.055, 0.155, size = 10.5, face = "bold", just = "left")
  txt(top_label, 0.055, 0.122, size = 9.2, just = "left")

  legend_y <- c(0.180, 0.145, 0.110, 0.075)
  legend_x <- 0.705
  legend_labels <- c(
    "Q1 concordant up", "Q2 temporal down / spatial up",
    "Q3 concordant down", "Q4 temporal up / spatial down"
  )
  for (i in seq_along(state_order)) {
    grid.roundrect(x = unit(legend_x, "npc"), y = unit(legend_y[i], "npc"),
                   width = unit(0.022, "npc"), height = unit(0.018, "npc"),
                   r = unit(0.006, "snpc"),
                   gp = gpar(col = state_col[state_order[i]], fill = state_fill[state_order[i]], lwd = 1.2))
    txt(legend_labels[i], legend_x + 0.020, legend_y[i], size = 8.8, just = "left")
  }
  popViewport()
}

base <- file.path(out_dir, state_tag)
cairo_pdf(paste0(base, ".pdf"), width = 13.2, height = 10.2,
          family = font_main, bg = "white")
draw_panel(); dev.off()
png(paste0(base, ".png"), width = 13.2, height = 10.2,
    units = "in", res = 450, type = "cairo", bg = "white")
draw_panel(); dev.off()

message("Wrote: ", base, ".{pdf,png}")


# -------------------- 04 PPI network scoring and four-panel figure --------------------
# The PPI analysis is part of the five-way shared-DEG workflow. The Python
# helper reads 12_Venn_All/02.Venn_DEGs.csv, reconstructs the frozen STRING
# network, and calculates weighted PPI strength, Degree, MCC and MCODE scores.
# The R plotting helper then generates the four-panel PDF, SVG and PNG outputs.

project_root <- normalizePath(".", winslash = "/", mustWork = TRUE)
python_candidates <- Sys.which(c("python", "python3"))
python_candidates <- unname(python_candidates[nzchar(python_candidates)])
if (!length(python_candidates)) {
  stop(
    "Python was not found. Install Python 3.12 and the packages in ",
    "environment/python_requirements.txt."
  )
}
python <- python_candidates[1]
score_script <- file.path(project_root, "scripts", "ppi_network_scoring.py")
score_arguments <- shQuote(score_script)
if (tolower(Sys.getenv("PHOTOSEQ_REFRESH_STRING", "false")) %in%
    c("1", "true", "yes")) {
  score_arguments <- c(score_arguments, "--refresh-string")
}
status <- system2(python, args = score_arguments)
if (!identical(status, 0L)) {
  stop("PPI network scoring failed (exit status ", status, ")")
}

source(
  file.path(project_root, "scripts", "ppi_network_plot.R"),
  local = new.env(parent = globalenv())
)
message("Integrated shared-gene and PPI analyses completed.")
