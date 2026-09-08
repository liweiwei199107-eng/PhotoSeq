rm(list = ls()); gc()

ORIGINAL_DIR <- normalizePath(".", winslash = "/", mustWork = TRUE)
output <- file.path(ORIGINAL_DIR, "00_prepared_data")
dir.create(output, recursive = TRUE, showWarnings = FALSE)
validation_dir <- file.path(ORIGINAL_DIR, "validation")
dir.create(validation_dir, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(DESeq2)
  library(sva)
})

metadata <- read.csv(
  file.path("input", "PhotoSeq_sample_metadata.csv"),
  stringsAsFactors = FALSE, check.names = FALSE
)
required_meta <- c("sample_id", "condition", "replicate", "sequencing_batch")
stopifnot(all(required_meta %in% names(metadata)))
stopifnot(!anyDuplicated(metadata$sample_id), nrow(metadata) == 24)
stopifnot(sum(metadata$sequencing_batch == "Batch1") == 10,
          sum(metadata$sequencing_batch == "Batch2") == 14)

sample_id_from_header <- function(x, kind) {
  rep <- sub(".*_([0-9]+)$", "\\1", x)
  label <- sub("_[0-9]+$", "", x)
  if (kind == "temporal") {
    label <- c(
      BreastEarly = "Breast_early",
      BreastMid = "Breast_mid",
      BreastLate = "Breast_late",
      LungMetastasis = "Lung_metastasis",
      `Breast Early` = "Breast_early",
      `Breast Mid` = "Breast_mid",
      `Breast Late` = "Breast_late",
      `Lung Metastasis` = "Lung_metastasis"
    )[label]
  } else {
    label <- c(
      BreastMidTumor = "Breast_mid",
      BreastMidNormalAdjacent = "Breast_normal_adjacent",
      BreastMidNormal = "Breast_normal",
      LungMetastasis = "Lung_metastasis",
      LungNormalAdjacent = "Lung_normal_adjacent",
      LungNormal = "Lung_normal",
      `Breast Mid(Tumor)` = "Breast_mid",
      `Breast Mid(Normal adjacent)` = "Breast_normal_adjacent",
      `Breast Mid(Normal)` = "Breast_normal",
      `Lung Metastasis` = "Lung_metastasis",
      `Lung Normal adjacent` = "Lung_normal_adjacent",
      `Lung Normal` = "Lung_normal"
    )[label]
  }
  if (anyNA(label) || anyNA(rep) || any(!grepl("^[0-9]+$", rep))) {
    stop("Unrecognized matrix sample header: ", paste(x[is.na(label) | is.na(rep)], collapse = ", "))
  }
  paste0(unname(label), "_N", rep)
}

read_matrix <- function(path, kind) {
  raw <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE,
                  na.strings = c("", "NA"))
  if (ncol(raw) < 2) stop("Count matrix has no sample columns: ", path)
  headers <- colnames(raw)[-1]
  ids <- sample_id_from_header(headers, kind)
  if (anyDuplicated(ids)) stop("Duplicate sample IDs after explicit matrix mapping: ", path)
  if (!all(ids %in% metadata$sample_id)) {
    stop("Matrix samples missing from metadata: ", paste(setdiff(ids, metadata$sample_id), collapse = ", "))
  }
  genes <- as.character(raw[[1]])
  keep_gene <- !is.na(genes) & nzchar(genes)
  m <- as.matrix(raw[keep_gene, -1, drop = FALSE])
  suppressWarnings(storage.mode(m) <- "numeric")
  if (anyNA(m) || any(!is.finite(m)) || any(m < 0) || any(m != floor(m))) {
    stop("Invalid non-negative integer count values in: ", path)
  }
  genes <- genes[keep_gene]
  m <- rowsum(m, group = genes, reorder = FALSE)
  high_load_loci <- c("ENSMUSG00000136525.1", "ENSMUSG00000119584.1")
  m <- m[!rownames(m) %in% high_load_loci, , drop = FALSE]
  colnames(m) <- ids
  meta <- metadata[match(ids, metadata$sample_id), , drop = FALSE]
  stopifnot(all(meta$sample_id == colnames(m)))
  list(counts = m, meta = meta)
}

spatial_input <- file.path("input", "PhotoSeq_Spatial_GeneCounts.csv")
temporal_input <- file.path("input", "PhotoSeq_Temporal_GeneCounts.csv")
spatial <- read_matrix(spatial_input, "spatial")
temporal <- read_matrix(temporal_input, "temporal")

spatial$meta$group <- c(
  Breast_mid = "Breast_tumor",
  Breast_normal_adjacent = "Breast_normal_adjacent",
  Breast_normal = "Breast_normal",
  Lung_metastasis = "Lung_metastasis",
  Lung_normal_adjacent = "Lung_normal_adjacent",
  Lung_normal = "Lung_normal"
)[spatial$meta$condition]
temporal$meta$group <- temporal$meta$condition

write_subset <- function(x, prefix) {
  x$meta <- x$meta[order(x$meta$condition, x$meta$replicate), , drop = FALSE]
  x$counts <- x$counts[, x$meta$sample_id, drop = FALSE]
  dat1 <- x$counts
  group <- x$meta[, c("sample_id", "group", "sequencing_batch")]
  names(group)[1] <- "id"
  save(dat1, file = file.path(output, paste0(prefix, "_dat_count.RData")))
  write.csv(dat1, file.path(output, paste0(prefix, "_dat_count.csv")), quote = FALSE, row.names = TRUE)
  write.csv(group, file.path(output, paste0(prefix, "_group.csv")), quote = FALSE, row.names = FALSE)
  x
}

spatial <- write_subset(spatial, "01.Spatial")
temporal <- write_subset(temporal, "02.Temporal")

temporal_unique <- temporal$meta[!temporal$meta$condition %in% c("Breast_mid", "Lung_metastasis"), , drop = FALSE]
common_rows <- intersect(rownames(spatial$counts), rownames(temporal$counts))
exp <- cbind(
  spatial$counts[common_rows, , drop = FALSE],
  temporal$counts[common_rows, temporal_unique$sample_id, drop = FALSE]
)
integrated_meta <- rbind(
  spatial$meta[, c("sample_id", "group", "sequencing_batch", "replicate")],
  temporal_unique[, c("sample_id", "group", "sequencing_batch", "replicate")]
)
stopifnot(nrow(integrated_meta) == 24, !anyDuplicated(integrated_meta$sample_id))
integrated_meta <- integrated_meta[match(colnames(exp), integrated_meta$sample_id), , drop = FALSE]
stopifnot(all(integrated_meta$sample_id == colnames(exp)))
integrated_levels <- c(
  "Breast_early", "Breast_late", "Breast_tumor", "Breast_normal_adjacent",
  "Breast_normal", "Lung_metastasis", "Lung_normal_adjacent", "Lung_normal"
)
integrated_meta$group <- factor(integrated_meta$group, levels = integrated_levels)
integrated_meta <- integrated_meta[order(integrated_meta$group, integrated_meta$replicate), , drop = FALSE]
exp <- exp[, integrated_meta$sample_id, drop = FALSE]

exp2 <- sva::ComBat_seq(
  counts = as.matrix(round(exp)),
  batch = factor(integrated_meta$sequencing_batch, levels = c("Batch1", "Batch2")),
  group = integrated_meta$group
)
colnames(exp2) <- integrated_meta$sample_id
write.csv(integrated_meta, file.path(output, "03.All_group.csv"), quote = FALSE, row.names = FALSE)
write.csv(integrated_meta, file.path(ORIGINAL_DIR, "validation", "integrated_metadata_used.csv"), quote = FALSE, row.names = FALSE)
save(exp2, file = file.path(output, "03.All_dat_count.RData"))
write.csv(exp2, file.path(output, "03.All_dat_count.csv"), quote = FALSE, row.names = TRUE)

make_vst <- function(counts, meta, apply_combat_seq = FALSE) {
  counts <- as.matrix(round(counts[, meta$sample_id, drop = FALSE]))
  if (apply_combat_seq) {
    if (length(unique(meta$sequencing_batch)) < 2) stop("ComBat_seq requires at least two sequencing batches.")
    counts <- sva::ComBat_seq(
      counts = counts,
      batch = factor(meta$sequencing_batch, levels = c("Batch1", "Batch2")),
      group = factor(meta$group)
    )
    colnames(counts) <- meta$sample_id
  }
  cd <- S4Vectors::DataFrame(
    row.names = meta$sample_id,
    condition = factor(meta$group),
    sequencing_batch = factor(meta$sequencing_batch, levels = c("Batch1", "Batch2"))
  )
  dds <- DESeq2::DESeqDataSetFromMatrix(counts, cd, design = ~ sequencing_batch + condition)
  dds <- dds[rowSums(DESeq2::counts(dds) >= 1) >= 3, ]
  dds <- DESeq2::DESeq(dds, quiet = TRUE)
  list(counts = counts, vst = as.matrix(SummarizedExperiment::assay(DESeq2::vst(dds, blind = FALSE))))
}

spatial_norm <- make_vst(spatial$counts, spatial$meta, apply_combat_seq = TRUE)
combat_counts <- spatial_norm$counts
dat_norm <- spatial_norm$vst
save(combat_counts, file = file.path(output, "01.Spatial_dat_count_ComBatSeq.RData"))
write.csv(combat_counts, file.path(output, "01.Spatial_dat_count_ComBatSeq.csv"), quote = FALSE, row.names = TRUE)
save(dat_norm, file = file.path(output, "01.Spatial_dat_norm_ComBatSeq.RData"))
write.csv(dat_norm, file.path(output, "01.Spatial_dat_norm_ComBatSeq.csv"), quote = FALSE, row.names = TRUE)

temporal_norm <- make_vst(temporal$counts, temporal$meta, apply_combat_seq = TRUE)
combat_counts <- temporal_norm$counts
dat_norm <- temporal_norm$vst
save(combat_counts, file = file.path(output, "02.Temporal_dat_count_ComBatSeq.RData"))
write.csv(combat_counts, file.path(output, "02.Temporal_dat_count_ComBatSeq.csv"), quote = FALSE, row.names = TRUE)
save(dat_norm, file = file.path(output, "02.Temporal_dat_norm_ComBatSeq.RData"))
write.csv(dat_norm, file.path(output, "02.Temporal_dat_norm_ComBatSeq.csv"), quote = FALSE, row.names = TRUE)

integrated_norm <- make_vst(exp2, integrated_meta, apply_combat_seq = FALSE)
dat_norm <- integrated_norm$vst
save(dat_norm, file = file.path(output, "03.All_dat_norm.RData"))
write.csv(dat_norm, file.path(output, "03.All_dat_norm.csv"), quote = FALSE, row.names = TRUE)

write.csv(data.frame(
  matrix = c("spatial", "temporal", "integrated"),
  input = c(spatial_input, temporal_input, "deduplicated overlap of the two supplied matrices"),
  combat_seq = TRUE, combat_seq_stage = "before VST", preserve_condition = TRUE,
  vst_design = rep("~ sequencing_batch + condition", 3), stringsAsFactors = FALSE
), file.path(ORIGINAL_DIR, "validation", "batch_correction_summary.csv"), quote = FALSE, row.names = FALSE)
write.csv(data.frame(
  check = c("unique_profiles", "Batch1", "Batch2", "spatial_profiles", "spatial_Batch1", "spatial_Batch2", "temporal_profiles", "temporal_Batch1", "temporal_Batch2"),
  value = c(24, sum(metadata$sequencing_batch == "Batch1"), sum(metadata$sequencing_batch == "Batch2"), nrow(spatial$meta), sum(spatial$meta$sequencing_batch == "Batch1"), sum(spatial$meta$sequencing_batch == "Batch2"), nrow(temporal$meta), sum(temporal$meta$sequencing_batch == "Batch1"), sum(temporal$meta$sequencing_batch == "Batch2"))
), file.path(ORIGINAL_DIR, "validation", "metadata_integrity_summary.csv"), quote = FALSE, row.names = FALSE)
