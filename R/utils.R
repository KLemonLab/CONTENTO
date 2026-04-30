#' Build gene sets from bacterial annotation columns
#'
#' Converts functional annotation columns (COG24, KEGG, etc.) into a named list
#' of gene vectors, matching the format returned by msigdbr for compatibility
#' with fgsea functions.
#'
#' @param annot_df Data frame containing gene annotations with a Geneid column
#' @param column_name Character string specifying which annotation column to use
#'
#' @return Named list where each element is a character vector of Geneids
build_bacterial_genesets <- function(annot_df, column_name) {
  # Special handling for KEGG_Class hierarchical levels
  if (column_name == "func_KEGG_Class_L2") {
    return(build_kegg_class_genesets(annot_df, level = 2))
  }
  
  if (column_name == "func_KEGG_Class_L3") {
    return(build_kegg_class_genesets(annot_df, level = 3))
  }
  
  # Check if column exists for standard annotations
  if (!column_name %in% colnames(annot_df)) {
    return(list())
  }
  
  # Standard handling for COG24 and KEGG_Module (separated by !!!)
  annot_df %>%
    dplyr::filter(!is.na(.data[[column_name]]) & .data[[column_name]] != "") %>%
    dplyr::select(Geneid, pathway = dplyr::all_of(column_name)) %>%
    dplyr::mutate(pathway = strsplit(as.character(pathway), "!!!")) %>%
    tidyr::unnest(pathway) %>%
    dplyr::mutate(pathway = trimws(pathway)) %>%
    dplyr::filter(pathway != "") %>%
    dplyr::group_by(pathway) %>%
    dplyr::summarise(genes = list(Geneid), .groups = "drop") %>%
    tibble::deframe()
}

#' Build KEGG Class gene sets at specific hierarchical level
#'
#' Parses hierarchical KEGG Class annotations (format: "Level1; Level2; Level3")
#' Extracts a specific level for gene set analysis
#'
#' @param annot_df Data frame with Geneid and func_KEGG_Class columns
#' @param level Integer (2 or 3) specifying which hierarchical level to extract
#' @return Named list of gene vectors
build_kegg_class_genesets <- function(annot_df, level = 2) {
  if (!"func_KEGG_Class" %in% colnames(annot_df)) {
    return(list())
  }
  
  if (!level %in% c(2, 3)) {
    stop("level must be 2 or 3")
  }
  
  annot_df %>%
    dplyr::filter(!is.na(func_KEGG_Class) & func_KEGG_Class != "") %>%
    dplyr::select(Geneid, kegg_class = func_KEGG_Class) %>%
    # Split multiple KEGG classes separated by !!!
    dplyr::mutate(kegg_class = strsplit(as.character(kegg_class), "!!!")) %>%
    tidyr::unnest(kegg_class) %>%
    dplyr::mutate(kegg_class = trimws(kegg_class)) %>%
    dplyr::filter(kegg_class != "") %>%
    # Split hierarchical levels by ;
    dplyr::mutate(levels = strsplit(kegg_class, ";")) %>%
    dplyr::mutate(
      pathway = purrr::map_chr(levels, ~ if(length(.x) >= level) trimws(.x[level]) else NA_character_)
    ) %>%
    dplyr::select(Geneid, pathway) %>%
    dplyr::filter(!is.na(pathway) & pathway != "") %>%
    dplyr::group_by(pathway) %>%
    dplyr::summarise(genes = list(Geneid), .groups = "drop") %>%
    tibble::deframe()
}