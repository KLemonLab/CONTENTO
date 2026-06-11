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
    # GFF3 structural columns
    "seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attributes",
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
  as.list(gsea_cols) %>%
    setNames(gsea_cols)
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
  
  annot_df %>%
    dplyr::filter(!is.na(.data[[column_name]]) & .data[[column_name]] != "") %>%
    dplyr::select(Geneid, pathway = dplyr::all_of(column_name)) %>%
    dplyr::mutate(pathway = strsplit(as.character(pathway), ",")) %>%
    tidyr::unnest(pathway) %>%
    dplyr::mutate(pathway = trimws(pathway)) %>%
    dplyr::filter(pathway != "") %>%
    dplyr::group_by(pathway) %>%
    dplyr::summarise(genes = list(Geneid), .groups = "drop") %>%
    tibble::deframe()
}


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
    return(file_path_sans_ext(basename(input$seFile$name)))
  }
  
  # Last resort
  "RNASeq_results"
}