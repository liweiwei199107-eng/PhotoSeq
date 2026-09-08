project_root <- normalizePath(".", winslash = "/", mustWork = TRUE)
script_dir <- file.path(project_root, "scripts")

if (file.exists(file.path(project_root, "renv.lock"))) {
  if (!requireNamespace("renv", quietly = TRUE)) {
    stop("Install the renv package and run renv::restore() before starting the analysis.")
  }
  renv::load(project = project_root)
  Sys.setenv(R_LIBS = paste(.libPaths(), collapse = .Platform$path.sep))
}

rscript <- file.path(
  R.home("bin"),
  if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript"
)

python_candidates <- Sys.which(c("python", "python3"))
python_candidates <- unname(python_candidates[nzchar(python_candidates)])
if (!length(python_candidates)) {
  stop("Python was not found. Install Python 3.12 and the packages in environment/python_requirements.txt.")
}
python <- python_candidates[1]

steps <- list(
  c("R", "00_prepare_analysis_data.R"),
  c("R", "01_spatial_differential_expression.R"),
  c("R", "02_spatial_shared_genes.R"),
  c("R", "03_spatial_immune_cell_scoring.R"),
  c("R", "04_spatial_estimate_scoring.R"),
  c("R", "05_spatial_pathway_enrichment.R"),
  c("R", "06_spatial_gsva.R"),
  c("R", "07_spatial_wgcna.R"),
  c("R", "08_temporal_differential_expression.R"),
  c("R", "09_temporal_shared_genes.R"),
  c("R", "10_temporal_pathway_enrichment.R"),
  c("R", "11_temporal_gsva.R"),
  c("R", "12_integrated_shared_gene_analysis.R"),
  c("R", "13_integrated_wgcna.R"),
  c("R", "14_integrated_gsva.R"),
  c("R", "15_temporal_immune_cell_scoring.R"),
  c("R", "16_integrated_pathway_analysis.R"),
  c("R", "17_quality_control_and_benchmarking.R")
)

for (step in steps) {
  interpreter <- step[1]
  script_name <- step[2]
  script_path <- file.path(script_dir, script_name)
  if (!file.exists(script_path)) stop("Missing analysis script: ", script_path)
  message("Running ", script_name)
  executable <- if (identical(interpreter, "R")) rscript else python
  status <- system2(executable, args = shQuote(script_path))
  if (!identical(status, 0L)) stop("Analysis stopped after ", script_name, " (exit status ", status, ")")
}

message("All PhotoSeq downstream scripts completed.")
