# Script to convert the test CSV annotation files into .rds files and to intentionally set types
# Run from the package root (project) directory: Rscript inst/extdata/test_annotations/generate_annotations.R

# Source local helper if available (allows using save_annotation above)
if (file.exists("R/save_annotation.R")) {
  source("R/save_annotation.R")
}

library(readr)
library(dplyr)

files <- list(
  good = "inst/extdata/test_annotations/annotation_good.csv",
  missing_gid = "inst/extdata/test_annotations/annotation_missing_geneid.csv",
  numeric_gid = "inst/extdata/test_annotations/annotation_geneid_numeric.csv",
  multiline = "inst/extdata/test_annotations/annotation_multiline_names.csv",
  zerosets = "inst/extdata/test_annotations/annotation_zerosets.csv"
)

outdir <- "inst/extdata/test_annotations/rds"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# 1. Good annotation: ensure Geneid is character, prefer save_annotation when available
ann <- read_csv(files$good, show_col_types = FALSE)
ann$Geneid <- as.character(ann$Geneid)
if (exists("save_annotation")) {
  # write into the outdir with prefix 'annotation_good' (save_annotation writes prefix_annot.rds)
  save_annotation(ann, prefix = "annotation_good", dir = outdir, write_csv = FALSE, write_rds = TRUE)
} else {
  readr::write_rds(ann, file.path(outdir, "annotation_good.rds"))
}

# 2. Missing Geneid: write as-is (no Geneid column) — deliberately malformed for testing
ann2 <- read_csv(files$missing_gid, show_col_types = FALSE)
readr::write_rds(ann2, file.path(outdir, "annotation_missing_geneid.rds"))

# 3. Numeric Geneid: coerce Geneid to numeric type to trigger validation failure
ann3 <- read_csv(files$numeric_gid, show_col_types = FALSE)
ann3$Geneid <- as.numeric(ann3$Geneid)
readr::write_rds(ann3, file.path(outdir, "annotation_geneid_numeric.rds"))

# 4. Multiline names: ensure Regulator contains literal newlines (already in CSV)
ann4 <- read_csv(files$multiline, show_col_types = FALSE)
readr::write_rds(ann4, file.path(outdir, "annotation_multiline_names.rds"))

# 5. Zero-sets: unique tiny_group values so that filtering by minimum gene set size yields zero gene sets
ann5 <- read_csv(files$zerosets, show_col_types = FALSE)
readr::write_rds(ann5, file.path(outdir, "annotation_zerosets.rds"))

cat("Wrote .rds files to:", normalizePath(outdir), "\\n")
