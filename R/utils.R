#' Get base filename for downloads
#' 
#' Extracts a meaningful base name for downloaded files from the SE object
#' or uploaded filename.
#' 
#' @param input Shiny input object
#' @param state App state list containing reactive values
#' @return Character string to use as base filename
get_download_filename <- function(input, state) {
  # Try SE metadata first
  se_name <- tryCatch(
    metadata(state$se_obj())$dataset_name,
    error = function(e) NULL
  )
  if (!is.null(se_name) && nzchar(trimws(se_name))) {
    return(se_name)
  }
  
  # Fall back to uploaded filename
  if (!is.null(input$seFile$name)) {
    return(tools::file_path_sans_ext(basename(input$seFile$name)))
  }
  
  # Last resort
  "RNASeq_results"
}


#' Find built-in annotation file path for a given annotation name
#'
#' Searches first in the installed package, then falls back to the local
#' inst/annotations directory (useful during development).
#'
#' @param annotation Character string with the annotation name (without extension)
#'
#' @return Character path to the .rds file, or NULL if not found
find_annotation_path <- function(annotation) {
  path <- system.file("annotations", paste0(annotation, ".rds"), package = "RNASeqApp")
  if (nzchar(path)) {
    path
  } else {
    local_path <- file.path("inst", "annotations", paste0(annotation, ".rds"))
    if (file.exists(local_path)) local_path else NULL
  }
}


#' Extract all contrasts from a SummarizedExperiment into long format
#'
#' Reads rowData column prefixes (baseMean_, log2FoldChange_, etc.) and
#' metadata(se)$contrasts to build a tidy data frame with one row per
#' gene-contrast combination.
#'
#' @param se A SummarizedExperiment object with contrast results stored in
#'   rowData and contrast names in metadata(se)$contrasts
#'
#' @return A long-format data frame with columns: Geneid, contrast, baseMean,
#'   log2FC, log2FC_shrunk, lfcSE, stat, pvalue, padj
extract_de_results <- function(se) {
  rd       <- as.data.frame(rowData(se))
  gene_ids <- rownames(se)
  
  meta <- tryCatch(
    metadata(se),
    error = function(e) list()
  )
  
  contrast_names <- meta$contrasts
  
  if (is.null(contrast_names) || length(contrast_names) == 0) {
    stop("No contrasts found in metadata(se)$contrasts. ",
         "Ensure the SE object includes metadata(se)$contrasts as a character vector of contrast names.")
  }
  
  de_list <- lapply(contrast_names, function(cname) {
    get_col <- function(prefix) {
      col <- paste0(prefix, cname)
      if (col %in% colnames(rd)) rd[[col]] else rep(NA_real_, nrow(rd))
    }
    data.frame(
      Geneid        = gene_ids,
      contrast      = cname,
      baseMean      = get_col("baseMean_"),
      log2FC        = get_col("log2FoldChange_"),
      log2FC_shrunk = get_col("log2FoldChange_shrunk_"),
      lfcSE         = get_col("lfcSE_"),
      stat          = get_col("stat_"),
      pvalue        = get_col("pvalue_"),
      padj          = get_col("padj_"),
      stringsAsFactors = FALSE
    )
  })
  
  do.call(rbind, de_list)
}


#' Merge annotation data frame into a contrast data frame
#'
#' Performs a left join on Geneid with input validation. Returns a named list
#' so callers can distinguish success from failure without try/catch.
#'
#' @param de_df Data frame of DE results containing a Geneid column
#' @param annot Data frame of gene annotations containing a Geneid column
#'
#' @return Named list with either:
#'   \item{result}{Merged data frame on success}
#'   \item{message}{Diagnostic string on failure (result will be NULL)}
merge_annotation <- function(de_df, annot) {
  if (!is.data.frame(de_df) || !is.data.frame(annot)) {
    return(list(result = NULL, message = "Invalid data frame structure"))
  }
  
  if (nrow(annot) == 0) {
    return(list(result = NULL, message = "Annotation file is empty"))
  }
  
  if (!("Geneid" %in% colnames(de_df))) {
    return(list(result = NULL, message = "DE data missing 'Geneid' column (internal error)"))
  }
  
  if (!("Geneid" %in% colnames(annot))) {
    return(list(
      result = NULL,
      message = paste0("Annotation file must contain a 'Geneid' column. ",
                       "Found columns: ", paste(colnames(annot), collapse = ", "))
    ))
  }
  
  if (!is.character(annot[["Geneid"]])) {
    return(list(
      result = NULL,
      message = paste0("Column 'Geneid' must be character type, got ",
                       class(annot[["Geneid"]])[1], ". Cannot join.")
    ))
  }
  
  merged <- tryCatch(
    dplyr::left_join(de_df, annot, by = "Geneid"),
    error = function(e) NULL
  )
  
  if (!is.null(merged)) {
    list(result = merged)
  } else {
    list(result = NULL, message = "Join operation failed (internal error)")
  }
}


#' Get all potential gsea columns available in dataframe
#'
#' Returns all column names except common metadata columns (Geneid, symbol, locus_tag, etc.)
#' These are candidate columns for functional enrichment analysis.
#'
#' @param annot_df Data frame containing annotations
#'
#' @return Named list where names are column names and values are also column names
#'         (suitable for use in selectInput choices)
get_gsea_columns <- function(annot_df) {
  if (!is.data.frame(annot_df) || nrow(annot_df) == 0) {
    return(list())
  }
  
  # Exclude metadata/ID/structural columns not useful for enrichment
  exclude_cols <- c(
    # Structural columns
    "Geneid", "symbol", "seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attributes",
    # GFF3 metadata attributes
    "ID", "Parent", "Dbxref", "Name", "Ontology_term", "gbkey", "gene", "inference",
    "locus_tag", "product", "protein_id", "transl_table", "Note", "partial", "pseudo",
    "start_range", "end_range", "exception",
    # anvi'o internal columns (not useful for enrichment)
    "anvio_gene_callers_id", "anvio_acc_prodigal", "anvio_acc_GENBANK_LOCUS_TAG",
    "anvio_func_prodigal", "anvio_aa_sequence"
  )
  
  gsea_cols <- setdiff(colnames(annot_df), exclude_cols)
  
  # Create named list: display name -> column name
  as.list(gsea_cols) |>
    stats::setNames(gsea_cols)
}

#' Build gene sets from any annotation column
#'
#' Converts ANY functional annotation column into a named list of gene vectors.
#' Handles multiple values per gene separated by commas (GFF3 format).
#'
#' @param annot_df Data frame containing gene annotations with a Geneid column
#' @param column_name Character string specifying which annotation column to use
#'
#' @return Named list where each element is a character vector of Geneids
build_bacterial_genesets <- function(annot_df, column_name) {
  if (!column_name %in% colnames(annot_df)) {
    return(list())
  }
  
  annot_df |>
    dplyr::filter(!is.na(.data[[column_name]]) & .data[[column_name]] != "") |>
    dplyr::select(Geneid, pathway = dplyr::all_of(column_name)) |>
    dplyr::mutate(pathway = strsplit(as.character(pathway), "!!!")) |>
    tidyr::unnest(pathway) |>
    dplyr::mutate(pathway = trimws(pathway)) |>
    dplyr::filter(pathway != "") |>
    dplyr::group_by(pathway) |>
    dplyr::summarise(genes = list(Geneid), .groups = "drop") |>
    tibble::deframe()
}

#' Determine the default symbol column from available annotation columns
#'
#' Checks a prioritised list of common gene name column names and returns
#' the first match found, falling back to the first available column.
#'
#' @param annot_cols Character vector of column names from the annotation data frame
#'
#' @return Character string with the name of the best candidate symbol column
default_symbol_col <- function(annot_cols) {
  found <- intersect(c("gene", "Gene", "product", "symbol", "hgnc_symbol", "gene_name"), annot_cols)
  if (length(found) > 0) found[1] else annot_cols[1]
}


#' Apply a symbol column to a data frame with optional fallback
#'
#' Creates or overwrites a 'symbol' column using the specified primary column,
#' falling back to a secondary column when the primary value is NA or empty.
#'
#' @param df Data frame to modify
#' @param primary_col Character string naming the column to use as symbol
#' @param secondary_col Character string naming the fallback column, or "none"
#'   to disable fallback (default: "none")
#'
#' @return The input data frame with a 'symbol' column added or updated.
#'   Returns df unchanged if primary_col is not present.
apply_symbol <- function(df, primary_col, secondary_col = "none") {
  if (!primary_col %in% colnames(df)) return(df)
  if (!is.null(secondary_col) && secondary_col != "none" &&
      secondary_col %in% colnames(df)) {
    df |>
      dplyr::mutate(symbol = ifelse(
        is.na(.data[[primary_col]]) | as.character(.data[[primary_col]]) == "",
        as.character(.data[[secondary_col]]),
        as.character(.data[[primary_col]])
      ))
  } else {
    df |>
      dplyr::mutate(symbol = as.character(.data[[primary_col]]))
  }
}