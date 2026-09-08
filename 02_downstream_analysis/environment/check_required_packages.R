manifest_path <- file.path("environment", "required_R_packages.csv")
if (!file.exists(manifest_path)) {
  stop("Run this script from the 02_downstream_analysis directory.")
}

manifest <- read.csv(manifest_path, stringsAsFactors = FALSE, check.names = FALSE)
installed <- vapply(manifest$Package, requireNamespace, logical(1), quietly = TRUE)
installed_version <- vapply(
  manifest$Package,
  function(package_name) {
    if (!requireNamespace(package_name, quietly = TRUE)) return(NA_character_)
    as.character(utils::packageVersion(package_name))
  },
  character(1)
)

report <- data.frame(
  Package = manifest$Package,
  Installed = installed,
  InstalledVersion = installed_version,
  RecordedVersion = manifest$RecordedVersion,
  stringsAsFactors = FALSE
)

print(report, row.names = FALSE)
if (any(!installed)) {
  stop("Missing required packages: ", paste(report$Package[!installed], collapse = ", "))
}

message("All packages called directly by the analysis scripts are installed.")
