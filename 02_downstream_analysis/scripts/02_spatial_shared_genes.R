rm(list = ls()); gc()
ORIGINAL_DIR <- "."
output <- file.path(ORIGINAL_DIR, "02_Venn_Spatial")
if (!dir.exists(output)) {dir.create(output, recursive = TRUE)}
setwd(ORIGINAL_DIR)

source(file.path(ORIGINAL_DIR, "scripts", "area_proportional_venn_helpers.R"))

root_dir <- file.path(ORIGINAL_DIR, "01_DEGs_Spatial")

if (!dir.exists(root_dir)) {
  stop("Path does not exist: ", root_dir)
}

folders <- list.dirs(root_dir, recursive = FALSE, full.names = TRUE)
folders <- folders[folders != root_dir]

if (length(folders) == 0) {
  stop("No comparison subfolders found under: ", root_dir)
}

deg_list <- list()

for (folder in folders) {
  folder_name <- basename(folder)
  sig_file <- list.files(folder, pattern = "^02\\.Sig.*\\.csv$", full.names = TRUE)
  if (length(sig_file) == 0) {
    message("No 02.Sig file found; skipping: ", folder_name)
    next
  }
  if (length(sig_file) > 1) {
    sig_file <- sig_file[1]
  }
  deg <- read.csv(sig_file, stringsAsFactors = FALSE)
  if (!"GeneSymbol" %in% colnames(deg)) {
    message("GeneSymbol column missing; skipping: ", folder_name)
    next
  }
  deg_list[[folder_name]] <- unique(deg$GeneSymbol)
  message("Read: ", folder_name, " | genes: ", length(deg_list[[folder_name]]))
}

set_a <- deg_list[["Breast_tumor_vs_Breast_normal_adjacent"]]
set_b <- deg_list[["Lung_metastasis_vs_Breast_tumor"]]

if (is.null(set_a) || is.null(set_b)) {
  stop("Required DEG sets are missing.")
}

sets <- list(
  DEG1 = set_a,
  DEG2 = set_b
)

p <- plot_area_venn2(
  sets,
  colors = c("#94be98", "#fcbe6e"),
  title = NULL
)

save_venn_plot(p, file.path(output, "01.Venn_DEGs"), width = 6, height = 6)

genes <- intersect(set_a, set_b)
write.csv(data.frame(x = genes), file.path(output, "02.Venn_DEGs.csv"), quote = FALSE, row.names = FALSE)

summary_df <- data.frame(
  set = c("DEG1", "DEG2", "DEG1_and_DEG2"),
  count = c(length(set_a), length(set_b), length(genes))
)
write.csv(summary_df, file.path(output, "01.Venn_DEGs_area_summary.csv"), quote = FALSE, row.names = FALSE)
