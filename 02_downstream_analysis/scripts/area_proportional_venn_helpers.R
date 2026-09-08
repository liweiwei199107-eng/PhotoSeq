suppressPackageStartupMessages(library(ggplot2))

circle_polygon <- function(cx, cy, r, set_name, n = 720) {
  theta <- seq(0, 2 * pi, length.out = n)
  data.frame(
    x = cx + r * cos(theta),
    y = cy + r * sin(theta),
    set = set_name
  )
}

circle_overlap_area <- function(r1, r2, d) {
  if (d >= r1 + r2) return(0)
  if (d <= abs(r1 - r2)) return(pi * min(r1, r2)^2)
  part1 <- r1^2 * acos((d^2 + r1^2 - r2^2) / (2 * d * r1))
  part2 <- r2^2 * acos((d^2 + r2^2 - r1^2) / (2 * d * r2))
  part3 <- 0.5 * sqrt((-d + r1 + r2) * (d + r1 - r2) * (d - r1 + r2) * (d + r1 + r2))
  part1 + part2 - part3
}

distance_for_overlap <- function(r1, r2, target_area) {
  max_overlap <- pi * min(r1, r2)^2
  target_area <- max(0, min(target_area, max_overlap))
  if (target_area <= 0) return(r1 + r2)
  if (target_area >= max_overlap) return(abs(r1 - r2))
  f <- function(d) circle_overlap_area(r1, r2, d) - target_area
  uniroot(f, interval = c(abs(r1 - r2), r1 + r2), tol = 1e-7)$root
}

save_venn_plot <- function(plot, base, width = 6, height = 6, dpi = 600) {
  ggsave(paste0(base, ".pdf"), plot, width = width, height = height, bg = "white")
  png(paste0(base, ".png"), width = width, height = height, units = "in", res = dpi, bg = "white")
  print(plot)
  dev.off()
  svg(paste0(base, ".svg"), width = width, height = height, family = "sans")
  print(plot)
  dev.off()
}

adaptive_text_size <- function(values, min_size = 4.5, max_size = 9) {
  values <- as.numeric(values)
  if (length(values) == 0 || max(values, na.rm = TRUE) <= 0) {
    return(rep(min_size, length(values)))
  }
  min_size + (max_size - min_size) * sqrt(pmax(values, 0) / max(values, na.rm = TRUE))
}

format_region_labels <- function(values, percent_threshold = 0.07) {
  values <- as.numeric(values)
  total <- sum(values, na.rm = TRUE)
  pct <- if (total > 0) values / total else rep(0, length(values))
  ifelse(
    pct >= percent_threshold,
    sprintf("%s\n(%.2f%%)", values, pct * 100),
    as.character(values)
  )
}

split_label_layers <- function(labels, max_radius, pct_scale = 0.54) {
  labels$count_label <- sub("\\n.*$", "", labels$label)
  labels$pct_label <- ifelse(grepl("\\n", labels$label), sub("^.*\\n", "", labels$label), "")
  labels$has_pct <- labels$pct_label != ""
  labels$count_y <- ifelse(labels$has_pct, labels$y + 0.055 * max_radius, labels$y)
  labels$pct_y <- labels$y - 0.105 * max_radius
  labels$pct_size <- pmax(labels$size * pct_scale, 2.2)
  labels
}

find_region_label_positions <- function(centers, patterns, labels, sizes, grid_n = 360) {
  xmin <- min(centers$x - centers$r)
  xmax <- max(centers$x + centers$r)
  ymin <- min(centers$y - centers$r)
  ymax <- max(centers$y + centers$r)
  grid <- expand.grid(
    x = seq(xmin, xmax, length.out = grid_n),
    y = seq(ymin, ymax, length.out = grid_n)
  )
  dist_mat <- sapply(seq_len(nrow(centers)), function(i) {
    sqrt((grid$x - centers$x[i])^2 + (grid$y - centers$y[i])^2)
  })
  if (is.null(dim(dist_mat))) {
    dist_mat <- matrix(dist_mat, ncol = 1)
  }
  inside_mat <- sweep(dist_mat, 2, centers$r, "<=")

  out <- lapply(seq_len(nrow(patterns)), function(i) {
    pattern <- as.logical(patterns[i, ])
    keep <- rep(TRUE, nrow(grid))
    margins <- matrix(NA_real_, nrow = nrow(grid), ncol = length(pattern))
    for (j in seq_along(pattern)) {
      if (pattern[j]) {
        keep <- keep & inside_mat[, j]
        margins[, j] <- centers$r[j] - dist_mat[, j]
      } else {
        keep <- keep & !inside_mat[, j]
        margins[, j] <- dist_mat[, j] - centers$r[j]
      }
    }
    candidates <- which(keep)
    if (length(candidates) == 0) {
      # Fallback to a weighted centre when a mathematically exact region is
      # too small for the sampling grid.
      positive <- which(pattern)
      if (length(positive) == 0) positive <- seq_len(nrow(centers))
      return(data.frame(
        x = mean(centers$x[positive]),
        y = mean(centers$y[positive]),
        label = labels[i],
        size = sizes[i]
      ))
    }
    score <- apply(margins[candidates, , drop = FALSE], 1, min)
    best <- candidates[which.max(score)]
    data.frame(
      x = grid$x[best],
      y = grid$y[best],
      label = labels[i],
      size = sizes[i]
    )
  })
  do.call(rbind, out)
}

plot_area_venn2 <- function(sets, colors = c("#94be98", "#fcbe6e"), title = NULL) {
  stopifnot(length(sets) == 2)
  sets <- lapply(sets, unique)
  names(sets) <- names(sets) %||% c("Set 1", "Set 2")
  a <- length(sets[[1]])
  b <- length(sets[[2]])
  ab <- length(intersect(sets[[1]], sets[[2]]))
  scale_factor <- 1 / sqrt(pi)
  r1 <- sqrt(a) * scale_factor
  r2 <- sqrt(b) * scale_factor
  d <- distance_for_overlap(r1, r2, ab * scale_factor^2 * pi)
  c1 <- c(0, 0)
  c2 <- c(d, 0)
  circles <- rbind(
    circle_polygon(c1[1], c1[2], r1, names(sets)[1]),
    circle_polygon(c2[1], c2[2], r2, names(sets)[2])
  )
  centers <- data.frame(
    set = names(sets),
    x = c(c1[1], c2[1]),
    y = c(c1[2], c2[2]),
    r = c(r1, r2)
  )
  region_counts <- c(length(setdiff(sets[[1]], sets[[2]])), ab, length(setdiff(sets[[2]], sets[[1]])))
  labels <- find_region_label_positions(
    centers = centers,
    patterns = data.frame(A = c(TRUE, TRUE, FALSE), B = c(FALSE, TRUE, TRUE)),
    labels = format_region_labels(region_counts, percent_threshold = 0.07),
    sizes = adaptive_text_size(region_counts, min_size = 4.8, max_size = 9.2)
  )
  labels <- split_label_layers(labels, max_radius = max(r1, r2), pct_scale = 0.52)
  set_labels <- data.frame(
    x = c(c1[1], c2[1]),
    y = c(-max(r1, r2) * 1.18, -max(r1, r2) * 1.18),
    label = names(sets)
  )
  p <- ggplot(circles, aes(x, y, group = set, fill = set)) +
    geom_polygon(alpha = 0.72, color = NA) +
    geom_path(color = "white", linewidth = 0.45, alpha = 0.7) +
    geom_text(
      data = labels,
      aes(x, count_y, label = count_label, size = size),
      inherit.aes = FALSE,
      fontface = "bold"
    ) +
    geom_text(
      data = labels[labels$has_pct, , drop = FALSE],
      aes(x, pct_y, label = pct_label, size = pct_size),
      inherit.aes = FALSE,
      fontface = "bold"
    ) +
    geom_text(data = set_labels, aes(x, y, label = label), inherit.aes = FALSE, size = 6.2) +
    scale_fill_manual(values = colors) +
    coord_equal(clip = "off") +
    scale_size_identity() +
    theme_void(base_size = 12) +
    theme(
      legend.position = "none",
      plot.title = element_text(hjust = 0.5, size = 15),
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA)
    ) +
    labs(title = title)
  p
}

plot_area_venn3 <- function(sets, colors = c("#94be98", "#fcbe6e", "#9ecae1"), title = NULL) {
  stopifnot(length(sets) == 3)
  sets <- lapply(sets, unique)
  names(sets) <- names(sets) %||% c("Set 1", "Set 2", "Set 3")
  n <- vapply(sets, length, numeric(1))
  scale_factor <- 1 / sqrt(pi)
  r <- sqrt(n) * scale_factor
  ab <- length(intersect(sets[[1]], sets[[2]]))
  ac <- length(intersect(sets[[1]], sets[[3]]))
  bc <- length(intersect(sets[[2]], sets[[3]]))
  d12 <- distance_for_overlap(r[1], r[2], ab * scale_factor^2 * pi)
  d13 <- distance_for_overlap(r[1], r[3], ac * scale_factor^2 * pi)
  d23 <- distance_for_overlap(r[2], r[3], bc * scale_factor^2 * pi)
  x3 <- (d13^2 - d23^2 + d12^2) / (2 * d12)
  y3_sq <- d13^2 - x3^2
  if (!is.finite(y3_sq) || y3_sq <= 0) {
    x3 <- d12 / 2
    y3 <- max(r) * 0.85
  } else {
    y3 <- sqrt(y3_sq)
  }
  centers <- data.frame(
    set = names(sets),
    x = c(0, d12, x3),
    y = c(0, 0, y3),
    r = r
  )
  circles <- do.call(rbind, lapply(seq_len(3), function(i) {
    circle_polygon(centers$x[i], centers$y[i], centers$r[i], centers$set[i])
  }))
  A <- sets[[1]]
  B <- sets[[2]]
  C <- sets[[3]]
  abc <- Reduce(intersect, sets)
  only_a <- setdiff(A, union(B, C))
  only_b <- setdiff(B, union(A, C))
  only_c <- setdiff(C, union(A, B))
  only_ab <- setdiff(intersect(A, B), C)
  only_ac <- setdiff(intersect(A, C), B)
  only_bc <- setdiff(intersect(B, C), A)
  region_counts <- c(length(only_a), length(only_b), length(only_c), length(only_ab), length(only_ac), length(only_bc), length(abc))
  labels <- find_region_label_positions(
    centers = centers,
    patterns = data.frame(
      A = c(TRUE, FALSE, FALSE, TRUE, TRUE, FALSE, TRUE),
      B = c(FALSE, TRUE, FALSE, TRUE, FALSE, TRUE, TRUE),
      C = c(FALSE, FALSE, TRUE, FALSE, TRUE, TRUE, TRUE)
    ),
    labels = format_region_labels(region_counts, percent_threshold = 0.07),
    sizes = adaptive_text_size(region_counts, min_size = 3.8, max_size = 6.9),
    grid_n = 420
  )
  labels <- split_label_layers(labels, max_radius = max(r), pct_scale = 0.38)
  set_labels <- data.frame(
    x = centers$x,
    y = c(
      centers$y[1] - centers$r[1] - 0.14 * max(r),
      centers$y[2] - centers$r[2] - 0.14 * max(r),
      centers$y[3] + centers$r[3] + 0.14 * max(r)
    ),
    label = names(sets)
  )
  ggplot(circles, aes(x, y, group = set, fill = set)) +
    geom_polygon(alpha = 0.68, color = NA) +
    geom_path(color = "white", linewidth = 0.45, alpha = 0.7) +
    geom_text(
      data = labels,
      aes(x, count_y, label = count_label, size = size),
      inherit.aes = FALSE,
      fontface = "bold"
    ) +
    geom_text(
      data = labels[labels$has_pct, , drop = FALSE],
      aes(x, pct_y, label = pct_label, size = pct_size),
      inherit.aes = FALSE,
      fontface = "bold"
    ) +
    geom_text(data = set_labels, aes(x, y, label = label), inherit.aes = FALSE, size = 5.2) +
    scale_fill_manual(values = colors) +
    coord_equal(clip = "off") +
    scale_size_identity() +
    theme_void(base_size = 12) +
    theme(
      legend.position = "none",
      plot.title = element_text(hjust = 0.5, size = 15),
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA)
    ) +
    labs(title = title)
}

`%||%` <- function(x, y) {
  if (is.null(x) || any(is.na(x)) || any(x == "")) y else x
}
