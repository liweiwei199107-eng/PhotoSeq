rm(list = ls()); gc()

suppressPackageStartupMessages(library(grid))

project_root <- normalizePath(".", winslash = "/", mustWork = TRUE)
ppi_dir <- file.path(project_root, "12_Venn_All", "04.PPI")
network_dir <- file.path(ppi_dir, "network")
output_dir <- file.path(ppi_dir, "figure")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

edges <- read.csv(
  file.path(network_dir, "ppi_edges.csv"),
  stringsAsFactors = FALSE
)
scores <- read.csv(
  file.path(network_dir, "ppi_node_scores.csv"),
  stringsAsFactors = FALSE
)
required_edge_columns <- c("gene_a", "gene_b", "combined_score")
required_score_columns <- c(
  "gene", "ppi_strength", "degree", "mcc",
  "mcode_module", "mcode_module_score"
)
if (!all(required_edge_columns %in% names(edges))) stop("PPI edge columns are incomplete")
if (!all(required_score_columns %in% names(scores))) stop("PPI score columns are incomplete")

edges <- edges[edges$combined_score > 0.4, , drop = FALSE]
network_nodes <- sort(unique(c(edges$gene_a, edges$gene_b)))
scores <- scores[match(network_nodes, scores$gene), , drop = FALSE]
if (anyNA(scores$gene)) stop("One or more network nodes are missing from the score table")

rank_top5 <- function(metric) {
  values <- setNames(scores[[metric]], scores$gene)
  names(sort(values, decreasing = TRUE))[seq_len(min(5L, length(values)))]
}
rank_top10 <- function(metric) {
  order_index <- order(-scores[[metric]], -scores$ppi_strength, scores$gene)
  head(scores$gene[order_index], 10L)
}

top_ppi <- rank_top5("ppi_strength")
top_mcc <- rank_top5("mcc")
top_degree <- rank_top5("degree")
top_mcode <- scores$gene[scores$mcode_module == "MCODE_1"]
top_mcode <- top_mcode[
  order(-scores$degree[match(top_mcode, scores$gene)], top_mcode)
]
top10_mcc <- rank_top10("mcc")
top10_degree <- rank_top10("degree")
consensus_candidates <- Reduce(intersect, list(top_mcc, top_degree, top_mcode))

# Fixed coordinates preserve the reported four-panel layout.
coordinates <- data.frame(
  gene = c(
    "Atp8a1", "B2m", "Pygm", "Cmya5", "Actg1", "Itga3", "Lama4",
    "Itga8", "Itga9", "Ctdspl", "Rpl10", "Mrpl50", "Dpysl3", "Cfl1"
  ),
  x = c(-0.42, -0.23, 0.22, 0.42, 0.00, 0.23, -0.19, -0.02, 0.21, 0.43, -0.39, -0.24, 0.05, 0.41),
  y = c(0.48, 0.51, 0.51, 0.46, 0.20, 0.18, -0.01, -0.17, -0.08, 0.03, -0.29, -0.45, -0.46, -0.31),
  stringsAsFactors = FALSE
)

fallback_ring <- function(missing, radius = 0.46) {
  if (!length(missing)) {
    return(data.frame(gene = character(), x = numeric(), y = numeric()))
  }
  angle <- seq(0, 2 * pi, length.out = length(missing) + 1)[-1]
  data.frame(
    gene = missing,
    x = radius * cos(angle),
    y = radius * sin(angle),
    stringsAsFactors = FALSE
  )
}
coordinates <- rbind(
  coordinates[coordinates$gene %in% network_nodes, ],
  fallback_ring(setdiff(network_nodes, coordinates$gene))
)

palette <- c(
  peach = "#F6C899", gray = "#BDBDBD", edge = "#9B9B9B",
  red1 = "#F23824", red2 = "#F04D23", orange1 = "#F56A1F",
  orange2 = "#F88A2A", orange3 = "#F6A55C", cyan = "#10AEE0"
)
score_palette <- colorRampPalette(
  c(palette["orange3"], palette["orange2"], palette["orange1"], palette["red1"])
)(101)

named_scores <- function(metric) setNames(scores[[metric]], scores$gene)
scale_radius <- function(values, low = 0.038, high = 0.066) {
  values[!is.finite(values)] <- 0
  if (!length(values) || max(values) == min(values)) {
    return(rep((low + high) / 2, length(values)))
  }
  low + (values - min(values)) / (max(values) - min(values)) * (high - low)
}

draw_network <- function(title, nodes, xy, highlight, metric, ppi_panel = FALSE) {
  xy <- xy[match(nodes, xy$gene), , drop = FALSE]
  x_values <- setNames(xy$x, xy$gene)
  y_values <- setNames(xy$y, xy$gene)
  panel_edges <- edges[
    edges$gene_a %in% nodes & edges$gene_b %in% nodes,
    ,
    drop = FALSE
  ]
  x_at <- function(gene) unit(0.5 + x_values[gene] * 0.88, "npc")
  y_at <- function(gene) unit(0.58 + y_values[gene] * 0.72, "npc")

  for (row_index in seq_len(nrow(panel_edges))) {
    gene_a <- panel_edges$gene_a[row_index]
    gene_b <- panel_edges$gene_b[row_index]
    grid.lines(
      x = unit(c(0.5 + x_values[gene_a] * 0.88, 0.5 + x_values[gene_b] * 0.88), "npc"),
      y = unit(c(0.58 + y_values[gene_a] * 0.72, 0.58 + y_values[gene_b] * 0.72), "npc"),
      gp = gpar(
        col = palette["edge"],
        lwd = 0.45 + pmax(0, pmin(1, (panel_edges$combined_score[row_index] - 0.4) / 0.6)) * 1.15,
        alpha = 0.72,
        lineend = "round"
      )
    )
  }

  values <- named_scores(metric)[nodes]
  radii <- setNames(
    scale_radius(
      values,
      if (ppi_panel) 0.035 else 0.038,
      if (ppi_panel) 0.060 else 0.066
    ),
    nodes
  )
  highlight_order <- highlight[highlight %in% nodes]
  highlight_values <- values[highlight_order]
  if (
    length(highlight_values) &&
      max(highlight_values, na.rm = TRUE) > min(highlight_values, na.rm = TRUE)
  ) {
    highlight_z <- (
      highlight_values - min(highlight_values, na.rm = TRUE)
    ) / (
      max(highlight_values, na.rm = TRUE) - min(highlight_values, na.rm = TRUE)
    )
  } else {
    highlight_z <- rep(0.55, length(highlight_values))
  }
  highlight_colors <- setNames(
    score_palette[1L + round(highlight_z * 100)],
    highlight_order
  )

  for (gene in nodes) {
    is_highlighted <- gene %in% highlight_order
    fill <- if (is_highlighted) {
      highlight_colors[gene]
    } else if (ppi_panel) {
      palette["peach"]
    } else {
      palette["gray"]
    }
    border <- if (is_highlighted) "#D84A19" else NA
    radius <- radii[gene]
    if (identical(gene, "Itga8")) {
      grid.circle(
        x = x_at(gene), y = y_at(gene), r = unit(radius + 0.013, "npc"),
        gp = gpar(fill = NA, col = palette["cyan"], lwd = 2.8)
      )
    }
    grid.circle(
      x = x_at(gene), y = y_at(gene), r = unit(radius, "npc"),
      gp = gpar(fill = fill, col = border, lwd = 0.65)
    )
    grid.text(
      gene, x = x_at(gene), y = y_at(gene),
      gp = gpar(fontfamily = "Arial", fontsize = 6.6, col = "#222222")
    )
  }
  grid.text(
    title, x = unit(0.5, "npc"), y = unit(0.055, "npc"),
    gp = gpar(fontfamily = "Arial", fontsize = 11.5, col = "#111111")
  )
}

draw_figure <- function() {
  grid.newpage()
  panel_specifications <- list(
    list(x = 0.27, y = 0.73, title = "PPI network", nodes = network_nodes, highlight = top_ppi, metric = "ppi_strength", ppi = TRUE),
    list(x = 0.73, y = 0.73, title = "MCC", nodes = top10_mcc, highlight = top_mcc, metric = "mcc", ppi = FALSE),
    list(x = 0.27, y = 0.31, title = "Degree", nodes = top10_degree, highlight = top_degree, metric = "degree", ppi = FALSE),
    # The complete network is retained here; MCODE_1 members are highlighted.
    list(x = 0.73, y = 0.31, title = "MCODE", nodes = network_nodes, highlight = top_mcode, metric = "mcode_module_score", ppi = FALSE)
  )
  for (panel in panel_specifications) {
    pushViewport(viewport(
      x = unit(panel$x, "npc"), y = unit(panel$y, "npc"),
      width = unit(0.43, "npc"), height = unit(0.38, "npc")
    ))
    draw_network(
      panel$title,
      panel$nodes,
      coordinates,
      panel$highlight,
      panel$metric,
      panel$ppi
    )
    popViewport()
  }

  if ("Itga8" %in% consensus_candidates) {
    grid.circle(
      x = unit(0.075, "npc"), y = unit(0.055, "npc"), r = unit(0.015, "npc"),
      gp = gpar(fill = NA, col = palette["cyan"], lwd = 2.8)
    )
    grid.circle(
      x = unit(0.075, "npc"), y = unit(0.055, "npc"), r = unit(0.0105, "npc"),
      gp = gpar(fill = palette["orange2"], col = "#D84A19", lwd = 0.7)
    )
    grid.text(
      "Itga8", x = unit(0.095, "npc"), y = unit(0.055, "npc"),
      just = c("left", "center"),
      gp = gpar(fontfamily = "Arial", fontsize = 8.2, fontface = "bold", col = "#111111")
    )
    grid.lines(
      x = unit(c(0.135, 0.175), "npc"), y = unit(c(0.055, 0.055), "npc"),
      arrow = arrow(type = "closed", length = unit(0.06, "in")),
      gp = gpar(col = "#111111", lwd = 1.0)
    )
    grid.text(
      "Highlighted candidate gene",
      x = unit(0.188, "npc"), y = unit(0.055, "npc"),
      just = c("left", "center"),
      gp = gpar(fontfamily = "Arial", fontsize = 7.2, fontface = "bold", col = "#1F4E79")
    )
  }
}

panel_rows <- list(
  list(panel = "PPI network", genes = network_nodes, highlight = top_ppi, metric = "ppi_strength"),
  list(panel = "MCC", genes = top10_mcc, highlight = top_mcc, metric = "mcc"),
  list(panel = "Degree", genes = top10_degree, highlight = top_degree, metric = "degree"),
  list(panel = "MCODE", genes = network_nodes, highlight = top_mcode, metric = "mcode_module_score")
)
source_nodes <- do.call(rbind, lapply(panel_rows, function(panel) {
  xy <- coordinates[coordinates$gene %in% panel$genes, , drop = FALSE]
  values <- named_scores(panel$metric)[xy$gene]
  data.frame(
    panel = panel$panel,
    gene = xy$gene,
    x = xy$x,
    y = xy$y,
    highlighted = xy$gene %in% panel$highlight,
    metric = panel$metric,
    metric_value = as.numeric(values),
    stringsAsFactors = FALSE
  )
}))
write.csv(
  source_nodes,
  file.path(output_dir, "ppi_four_panel_nodes.csv"),
  row.names = FALSE,
  quote = FALSE
)
write.csv(
  edges,
  file.path(output_dir, "ppi_four_panel_edges.csv"),
  row.names = FALSE,
  quote = FALSE
)
write.csv(
  data.frame(
    panel = c("PPI network", "MCC", "Degree", "MCODE"),
    highlighted = c(
      paste(top_ppi, collapse = ";"),
      paste(top_mcc, collapse = ";"),
      paste(top_degree, collapse = ";"),
      paste(top_mcode, collapse = ";")
    ),
    gray_nodes = c(
      paste(setdiff(network_nodes, top_ppi), collapse = ";"),
      paste(setdiff(top10_mcc, top_mcc), collapse = ";"),
      paste(setdiff(top10_degree, top_degree), collapse = ";"),
      paste(setdiff(network_nodes, top_mcode), collapse = ";")
    )
  ),
  file.path(output_dir, "ppi_four_panel_highlights.csv"),
  row.names = FALSE,
  quote = FALSE
)

base_path <- file.path(output_dir, "ppi_four_panel")
cairo_pdf(paste0(base_path, ".pdf"), width = 7.2, height = 6.1, family = "Arial", bg = "white")
draw_figure()
dev.off()
svg(paste0(base_path, ".svg"), width = 7.2, height = 6.1, family = "Arial", bg = "white")
draw_figure()
dev.off()
png(paste0(base_path, ".png"), width = 7.2, height = 6.1, units = "in", res = 600, type = "cairo", bg = "white")
draw_figure()
dev.off()

message(
  "PPI network: ", length(network_nodes), " non-isolated nodes, ",
  nrow(edges), " edges"
)
message("MCODE_1 is highlighted within the complete PPI network.")
message("PPI four-panel figure completed.")
