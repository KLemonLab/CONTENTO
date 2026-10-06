# ==============================================================================
# build_bacterial_annotations.R
#
# Batch version of create-genome-annotations.qmd (Steps 2.2 + 3.2) for
# bacterial genomes. No custom/RegPrecise merge.
#
# Expects, per genome, in INPUT_DIR:
#   <NAME>.gff          NCBI/Prokka GFF3 (required)
#   <NAME>_annot.txt    anvi-export-functions output (required)
#   <NAME>_calls.txt    anvi-script-process-genbank gene calls (optional, adds aa_sequence)
#
# Writes, per genome, to OUTPUT_DIR:
#   <NAME><OUT_SUFFIX>.rds / .csv
#
# Usage:
#   Rscript build_bacterial_annotations.R [INPUT_DIR] [OUTPUT_DIR]
#   or source() it from R after editing the CONFIG block.
# ==============================================================================

suppressPackageStartupMessages(library(tidyverse))

# ------------------------------------------------------------------------------
# CONFIG
# ------------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)

INPUT_DIR  <- if (length(args) >= 1) args[1] else here::here("inst/extdata/genomes_temp")
OUTPUT_DIR <- if (length(args) >= 2) args[2] else INPUT_DIR

GENOMES    <- NULL              # NULL = auto-detect every <NAME>.gff with a matching <NAME>_annot.txt
FEATURES   <- "CDS"             # GFF3 feature types kept, e.g. c("CDS", "tRNA", "rRNA", "tmRNA", "ncRNA")
OUT_SUFFIX <- "_annot_Anvio"    # output filename suffix
WRITE_CSV  <- TRUE

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------------------------
# HELPERS
# ------------------------------------------------------------------------------

# Parse a GFF3 attributes string into a named list
parse_gff3_attributes <- function(attributes_string) {
  pairs <- str_split(attributes_string, ";")[[1]]
  attr_list <- list()
  for (pair in pairs) {
    if (pair == "") next
    parts <- str_split(pair, "=", n = 2)[[1]]
    if (length(parts) == 2) attr_list[[parts[1]]] <- parts[2]
  }
  attr_list
}

# Step 2.2: GFF3 -> one row per gene, Geneid = locus_tag
read_gff3_annotations <- function(gff_path, features = FEATURES) {
  read_delim(
    gff_path,
    delim = "\t",
    col_names = c("seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attributes"),
    comment = "#",
    show_col_types = FALSE
  ) |>
    filter(feature %in% features) |>
    rowwise() |>
    mutate(attr_list = list(parse_gff3_attributes(attributes))) |>
    ungroup() |>
    unnest_wider(attr_list, names_repair = "universal") |>
    select(-attributes) |>
    dplyr::rename(Geneid = locus_tag) |>
    mutate(across(any_of(c("Ontology_term", "go_function", "go_process", "go_component")),
                  ~ gsub(",", "!!!", .x))) |>
    group_by(Geneid) |>
    summarise(                     # collapse CDS fragments (programmed frameshifts)
      start  = min(start),
      end    = max(end),
      strand = first(strand),
      across(-c(start, end, strand, frame), first),
      .groups = "drop"
    )
}

# Step 3.2: anvi'o functions -> wide, one row per gene, KEGG hierarchy split
read_anvio_annotations <- function(annot_path, calls_path = NULL) {
  anvio_annot <- read_delim(annot_path, delim = "\t", show_col_types = FALSE) |>
    mutate(source = as.factor(source)) |>
    mutate(across(where(is.character), trimws)) |>
    select(-e_value) |>
    dplyr::rename("acc" = "accession", "func" = "function") |>
    pivot_wider(
      names_from = source,
      values_from = c(acc, func),
      values_fn = list(
        acc  = ~ paste(unique(.), collapse = "!!!"),
        func = ~ paste(unique(.), collapse = "!!!"))
    )
  
  if (!"func_GENBANK_LOCUS_TAG" %in% colnames(anvio_annot)) {
    stop("No GENBANK_LOCUS_TAG source in ", basename(annot_path),
         " (was anvi-script-process-genbank run with --include-locus-tags-as-functions?)")
  }
  
  if ("func_KEGG_Class" %in% colnames(anvio_annot)) {
    anvio_annot <- anvio_annot |>
      separate_rows(func_KEGG_Class, sep = ",") |>
      separate(func_KEGG_Class,
               into = c("func_KEGG_Class", "func_KEGG_Subclass", "func_KEGG_Pathway"),
               sep = ";", extra = "drop", fill = "right") |>
      group_by(across(-starts_with("func_KEGG_"))) |>
      summarise(
        func_KEGG_Class    = paste(unique(na.omit(func_KEGG_Class)),    collapse = "!!!"),
        func_KEGG_Subclass = paste(unique(na.omit(func_KEGG_Subclass)), collapse = "!!!"),
        func_KEGG_Pathway  = paste(unique(na.omit(func_KEGG_Pathway)),  collapse = "!!!"),
        .groups = "drop"
      ) |>
      mutate(across(starts_with("func_KEGG_"), ~ na_if(., "")))
  }
  
  if (!is.null(calls_path) && file.exists(calls_path)) {
    anvio_calls <- read_delim(calls_path, delim = "\t", show_col_types = FALSE) |>
      select(gene_callers_id, aa_sequence)
    anvio_annot <- anvio_annot |>
      left_join(anvio_calls, by = "gene_callers_id")
  }
  
  anvio_annot
}

# Full pipeline for one genome
process_bacterial_genome <- function(GENOME_NAME) {
  gff_path   <- file.path(INPUT_DIR, paste0(GENOME_NAME, ".gff"))
  annot_path <- file.path(INPUT_DIR, paste0(GENOME_NAME, "_annot.txt"))
  calls_path <- file.path(INPUT_DIR, paste0(GENOME_NAME, "_calls.txt"))
  
  annotations <- read_gff3_annotations(gff_path)
  if (any(duplicated(annotations$Geneid))) stop("Duplicated Geneid after aggregation")
  
  anvio_annot <- read_anvio_annotations(annot_path, calls_path)
  
  annotations_merged <- annotations |>
    left_join(anvio_annot, by = c("Geneid" = "func_GENBANK_LOCUS_TAG")) |>
    rename_with(~ paste0("anvio_", .), -all_of(colnames(annotations)))
  
  out_base <- file.path(OUTPUT_DIR, paste0(GENOME_NAME, OUT_SUFFIX))
  saveRDS(annotations_merged, paste0(out_base, ".rds"))
  if (WRITE_CSV) write_csv(annotations_merged, paste0(out_base, ".csv"))
  
  tibble(
    genome       = GENOME_NAME,
    n_genes      = nrow(annotations_merged),
    n_anvio      = sum(annotations$Geneid %in% anvio_annot$func_GENBANK_LOCUS_TAG),
    has_aa       = "anvio_aa_sequence" %in% colnames(annotations_merged),
    n_anvio_cols = sum(startsWith(colnames(annotations_merged), "anvio_")),
    status       = "ok"
  )
}

# ------------------------------------------------------------------------------
# RUN
# ------------------------------------------------------------------------------

# Detect genomes
if (is.null(GENOMES)) {
  GENOMES <- list.files(INPUT_DIR, pattern = "\\.gff$") |> str_remove("\\.gff$")
}

inputs <- tibble(genome = GENOMES) |>
  mutate(
    gff   = file.exists(file.path(INPUT_DIR, paste0(genome, ".gff"))),
    annot = file.exists(file.path(INPUT_DIR, paste0(genome, "_annot.txt"))),
    calls = file.exists(file.path(INPUT_DIR, paste0(genome, "_calls.txt")))
  )

cat("Input:", INPUT_DIR, "\nOutput:", OUTPUT_DIR, "\n\n")
print(inputs, n = Inf)

skipped <- inputs |> filter(!gff | !annot)
if (nrow(skipped) > 0) {
  cat("\nSkipping (missing .gff or _annot.txt):", paste(skipped$genome, collapse = ", "), "\n")
}
to_run <- inputs |> filter(gff, annot) |> pull(genome)
if (length(to_run) == 0) stop("No genomes with both <NAME>.gff and <NAME>_annot.txt in ", INPUT_DIR)

# Process; one failing genome does not stop the batch
summary_tbl <- map_dfr(to_run, function(g) {
  cat("Processing", g, "...\n")
  tryCatch(
    process_bacterial_genome(g),
    error = function(e) tibble(genome = g, status = paste("ERROR:", conditionMessage(e)))
  )
})

cat("\n=== Summary ===\n")
print(summary_tbl, n = Inf, width = Inf)

low_match <- summary_tbl |> filter(status == "ok", n_anvio < 0.5 * n_genes)
if (nrow(low_match) > 0) {
  warning("Low anvi'o match (<50% of genes) for: ", paste(low_match$genome, collapse = ", "),
          ". GFF3 locus_tag and anvi'o GENBANK_LOCUS_TAG probably differ.")
}
