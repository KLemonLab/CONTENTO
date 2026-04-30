#' Build gene sets from bacterial annotation columns
#'
#' Converts functional annotation columns (COG24, KEGG, etc.) into a named list
#' of gene vectors, matching the format returned by msigdbr for compatibility
#' with fgsea functions.
#'
#' @param annot_df Data frame containing gene annotations with a Geneid column
#' @param column_name Character string specifying which annotation column to use
#'   (e.g., "func_COG24_CATEGORY", "func_COG24_PATHWAY", "func_KEGG_Module", "func_KEGG_Class")
#'
#' @return Named list where each element is a character vector of Geneids
#'   belonging to that pathway/category. Names are pathway/category labels.
#'   Returns empty list if column not found.
#'
#' @details
#' The function handles multiple annotations per gene separated by "!!!"
#' (e.g., "Pathway A!!!Pathway B"), splitting them and assigning the gene
#' to all relevant pathways.
#'
#' @examples
#' \dontrun{
#' annot <- data.frame(
#'   Geneid = c("gene1", "gene2", "gene3"),
#'   func_COG24_CATEGORY = c("Energy production", 
#'                           "Energy production!!!Translation", 
#'                           "Translation")
#' )
#' genesets <- build_bacterial_genesets(annot, "func_COG24_CATEGORY")
#' # Returns list with two elements:
#' # $`Energy production` = c("gene1", "gene2")
#' # $Translation = c("gene2", "gene3")
#' }
#'
#' @export
build_bacterial_genesets <- function(annot_df, column_name) {
  if (!column_name %in% colnames(annot_df)) {
    return(list())
  }
  
  # Split multiple annotations by "!!!" and create gene sets
  annot_df %>%
    dplyr::filter(!is.na(.data[[column_name]]) & .data[[column_name]] != "") %>%
    dplyr::select(Geneid, pathway = dplyr::all_of(column_name)) %>%
    # Split multiple pathways separated by !!!
    dplyr::mutate(pathway = strsplit(as.character(pathway), "!!!")) %>%
    tidyr::unnest(pathway) %>%
    dplyr::mutate(pathway = trimws(pathway)) %>%
    dplyr::filter(pathway != "") %>%
    dplyr::group_by(pathway) %>%
    dplyr::summarise(genes = list(Geneid), .groups = "drop") %>%
    tibble::deframe()  # Convert to named list
}