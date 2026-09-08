reference_dir <- file.path("reference", "immune")
gmt_path <- file.path(reference_dir, "CellReports_TISIDB_GMT.txt")
homology_path <- file.path(reference_dir, "HOM_MouseHumanSequence.rpt")
human_output <- file.path(reference_dir, "mmc3_ssGSEA.txt")
mouse_output <- file.path(reference_dir, "mmc3_ssGSEA_mouse.txt")

required_files <- c(gmt_path, homology_path)
if (any(!file.exists(required_files))) {
  stop("Missing immune-signature source file(s): ", paste(required_files[!file.exists(required_files)], collapse = ", "))
}

gmt_fields <- strsplit(readLines(gmt_path, warn = FALSE), "\t", fixed = TRUE)
human_signature <- do.call(rbind, lapply(gmt_fields, function(fields) {
  if (length(fields) < 3L) return(NULL)
  data.frame(
    Metagene = fields[-c(1, 2)],
    `Cell type` = fields[1],
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}))
human_signature <- unique(human_signature)

homology <- read.delim(
  homology_path,
  stringsAsFactors = FALSE,
  check.names = FALSE,
  quote = ""
)
required_columns <- c("DB Class Key", "Common Organism Name", "Symbol")
if (!all(required_columns %in% names(homology))) {
  stop("Homology table is missing required columns: ", paste(setdiff(required_columns, names(homology)), collapse = ", "))
}

mouse_rows <- lapply(seq_len(nrow(human_signature)), function(i) {
  human_gene <- human_signature$Metagene[i]
  cell_type <- human_signature[["Cell type"]][i]
  class_keys <- unique(homology[["DB Class Key"]][
    homology[["Common Organism Name"]] == "human" & homology$Symbol == human_gene
  ])
  mouse_genes <- unlist(lapply(class_keys, function(class_key) {
    homology$Symbol[
      homology[["Common Organism Name"]] == "mouse, laboratory" &
        homology[["DB Class Key"]] == class_key
    ]
  }), use.names = FALSE)
  mouse_genes <- unique(mouse_genes[nzchar(mouse_genes)])
  if (!length(mouse_genes)) return(NULL)
  data.frame(
    Metagene = mouse_genes,
    `Cell type` = cell_type,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
})
mouse_signature <- unique(do.call(rbind, mouse_rows))

write.table(human_signature, human_output, sep = "\t", quote = FALSE, row.names = FALSE)
write.table(mouse_signature, mouse_output, sep = "\t", quote = FALSE, row.names = FALSE)

message("Wrote ", nrow(human_signature), " human and ", nrow(mouse_signature), " mouse signature rows.")
