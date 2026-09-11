#' Save a CONTENTO-style annotation data.frame to CSV and RDS with validation
#' 
#' This helper validates an annotation data.frame and writes both CSV and .rds
#' representations. It enforces the required `Geneid` column, ensures `Geneid`
#' is character (coercing if necessary with a warning), replaces commas in common
#' GO/function columns with the `!!!` separator, removes duplicated Geneid rows
#' (keeping the first), and returns the cleaned data.frame invisibly.
#' 
#' @param annot A data.frame containing gene annotations
#' @param prefix File prefix (without extension) to write to
#' @param dir Directory to write files into (default: "inst/extdata/genomes")
#' @param write_csv Logical; whether to write a CSV file
#' @param write_rds Logical; whether to write a .rds file
#' @param compress Compression method passed to saveRDS when using base::saveRDS ("gzip", "bzip2", "xz" or NULL). Default: "xz"
#' @return Invisibly returns the cleaned annotation data.frame on success, or stops with an informative error on failure.
#' @export
save_annotation <- function(annot, prefix, dir = "inst/extdata/genomes", write_csv = TRUE, write_rds = TRUE, compress = "xz") {
  if (!is.data.frame(annot)) stop("annot must be a data.frame")
  if (!is.character(prefix) || length(prefix) != 1) stop("prefix must be a single string")
  if (!dir.exists(dir)) dir.create(dir, recursive = TRUE, showWarnings = FALSE)

  # Required column
  if (!("Geneid" %in% colnames(annot))) {
    stop("Annotation file must contain a 'Geneid' column. Found columns: ", paste(colnames(annot), collapse = ", "))
  }

  # Ensure Geneid is character (coerce with a warning if necessary)
  if (!is.character(annot$Geneid)) {
    warning("Coercing 'Geneid' column to character from ", class(annot$Geneid)[1])
    annot$Geneid <- as.character(annot$Geneid)
  }

  # Replace commas in common GO/function-like columns with '!!!' so multi-value fields
  # are represented consistently with CONTENTO expectations
  go_cols <- intersect(c("Ontology_term", "go_function", "go_process", "go_component", "Regulator"), colnames(annot))
  for (col in go_cols) {
    # Protect NAs
    annot[[col]] <- vapply(annot[[col]], function(x) {
      if (is.na(x)) return(NA_character_)
      gsub(",", "!!!", as.character(x), fixed = TRUE)
    }, FUN.VALUE = character(1))
  }

  # Remove duplicated Geneid rows, keeping the first occurrence
  if (any(duplicated(annot$Geneid))) {
    annot <- annot[!duplicated(annot$Geneid), , drop = FALSE]
    message("Removed duplicated Geneid rows; kept the first occurrence for each Geneid.")
  }

  csv_path <- file.path(dir, paste0(prefix, "_annot.csv"))
  rds_path <- file.path(dir, paste0(prefix, "_annot.rds"))

  if (write_csv) {
    # Use readr if available, otherwise fallback to utils::write.table
    if (requireNamespace("readr", quietly = TRUE)) {
      readr::write_csv(annot, csv_path)
    } else {
      utils::write.table(annot, file = csv_path, sep = ",", row.names = FALSE, col.names = TRUE, quote = TRUE)
    }
    message("Wrote CSV: ", csv_path)
  }

  if (write_rds) {
    # Prefer readr::write_rds for compatibility with tidyverse, fallback to saveRDS
    if (requireNamespace("readr", quietly = TRUE)) {
      readr::write_rds(annot, rds_path, compress = compress)
    } else {
      saveRDS(annot, file = rds_path, compress = compress)
    }
    message("Wrote RDS: ", rds_path)
  }

  invisible(annot)
}
