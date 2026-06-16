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


#' Add regulated/DE flag columns to a DE data frame
#'
#' @param df   Data frame with \code{log2FC} and \code{padj} columns
#' @param lfc_cut  Numeric log2FC threshold
#' @param padj_cut Numeric adjusted p-value threshold
#' @return \code{df} with two new columns: \code{regulated} ("up"/"down"/NA)
#'   and \code{DE} (logical)
add_de_flags <- function(df, lfc_cut, padj_cut) {
  df |>
    mutate(
      regulated = case_when(
        !is.na(padj) & !is.na(log2FC) & padj < padj_cut & log2FC >  lfc_cut ~ "up",
        !is.na(padj) & !is.na(log2FC) & padj < padj_cut & log2FC < -lfc_cut ~ "down",
        TRUE ~ NA_character_
      ),
      DE = !is.na(regulated)
    )
}


#' Get all potential gsea columns available in dataframe
#'
#' Returns all column names except common metadata columns (Geneid, symbol, locus_tag, etc.)
#' These are candidate columns for functional enrichment analysis.
#'
#' @param annot_df Data frame containing annotations
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


#' Get top DE genes as a character vector of Geneids
#' @param selected_data Data frame from selected_data() reactive with DE flags
#' @param top_n Integer number of top genes to return
#' @return Character vector of Geneids ordered by abs(log2FC)
get_top_de_genes <- function(selected_data, top_n) {
  selected_data |>
    filter(DE) |>
    arrange(desc(abs(log2FC))) |>
    distinct(Geneid) |>
    slice_head(n = top_n) |>
    pull(Geneid)
}


#' Extract and scale VST matrix for a set of genes
#' @param se_obj SummarizedExperiment with a vst assay
#' @param geneids Character vector of Geneids
#' @param de_df Optional data frame with Geneid and symbol columns for row labelling.
#'   If NULL or no symbol column, Geneids are used as rownames.
#' @return Scaled matrix with Geneids as rownames
get_vst_matrix <- function(se_obj, geneids, de_df = NULL) {
  mat <- assay(se_obj, "vst")[geneids[geneids %in% rownames(se_obj)], ]
  if (!is.null(de_df) && "symbol" %in% colnames(de_df)) {
    sym_lookup <- de_df |> distinct(Geneid, symbol)
    rownames(mat) <- sym_lookup$symbol[match(rownames(mat), sym_lookup$Geneid)]
  }
  scale(mat)
}


#' Build a wide comparison table of log2FC and DE status per contrast
#'
#' Reshapes a long-format differential expression data frame into a wide table
#' with one row per gene and paired columns for each contrast:
#' log2 fold-change and DE flag. Column order follows the user-provided
#' contrast vector, ensuring consistent display in tables.
#'
#' This function assumes DE flags (column `DE`) have already been added,
#' e.g. via `add_de_flags()`.
#'
#' @param df Data frame containing differential expression results.
#'   Must include columns: `Geneid`, `contrast`, `log2FC`, and `DE`.
#'   Optionally includes `symbol`.
#' @param contrasts Character vector of selected contrast names to include in the output.
#' @return A data frame in wide format with pairs of `log2FC_<contrast>` and `DE_<contrast>` for each contrast.
build_compare_table <- function(df, contrasts) {
  if (is.null(df)) return(NULL)
  
  has_symbol <- "symbol" %in% colnames(df)
  id_cols <- if (has_symbol) c("Geneid", "symbol") else "Geneid"
  
  wide <- df |>
    select(all_of(c(id_cols, "contrast", "log2FC", "DE"))) |>
    pivot_wider(
      id_cols = all_of(id_cols),
      names_from = contrast,
      values_from = c(log2FC, DE),
      names_sep = "_"
    ) |>
    mutate(
      across(starts_with("log2FC_"), ~ round(., 2)),
      across(starts_with("DE_"), ~ tidyr::replace_na(., FALSE))
    ) |> 
    filter(rowSums(across(starts_with("DE_"))) > 0)
  
  ordered_cols <- id_cols
  for (ct in contrasts) {
    ordered_cols <- c(
      ordered_cols,
      paste0("log2FC_", ct),
      paste0("DE_", ct)
    )
  }
  
  wide |>
    select(any_of(ordered_cols)) |>
    arrange(Geneid)
}


#' Load gene sets as a named list of gene ID vectors
#'
#' Single entry point for geneset loading across the app. Returns a named list
#' suitable for direct use in fgseaMultilevel() and geseca(). Handles both
#' Bacteria (annotation column) and Human (MSigDB) organisms.
#'
#' For Bacteria, delegates to build_bacterial_genesets().
#' For Human, fetches from msigdbr and splits into a named list by gs_name,
#' using ensembl_gene as the gene identifier.
#'
#' @param organism       Character; "Bacteria", "Human", or NULL/other
#' @param annotation_df  Data frame; required when organism == "Bacteria"
#' @param bacterial_source Character; annotation column name for gene set membership
#' @param gs_collection  Character; MSigDB collection code (e.g. "H", "C2")
#' @param gs_subcollection Character; MSigDB subcollection (e.g. "CP:REACTOME"),
#'   or "" / NULL to omit
#' @return Named list of character vectors (Geneids per pathway).
#'   Returns an empty list if inputs are invalid or nothing is found.
load_genesets <- function(organism,
                          annotation_df    = NULL,
                          bacterial_source = NULL,
                          gs_collection    = NULL,
                          gs_subcollection = NULL) {
  if (!is.null(organism) && organism == "Bacteria") {
    if (is.null(annotation_df) || is.null(bacterial_source)) return(list())
    build_bacterial_genesets(annotation_df, bacterial_source)
    
  } else {
    if (is.null(gs_collection)) return(list())
    use_sub <- !is.null(gs_subcollection) && nzchar(gs_subcollection)
    genesets <- if (use_sub) {
      msigdbr(species = "Homo sapiens",
              collection    = gs_collection,
              subcollection = gs_subcollection)
    } else {
      msigdbr(species = "Homo sapiens",
              collection = gs_collection)
    }
    genesets |>
      split(genesets$gs_name) |>
      lapply(function(x) x$ensembl_gene)
  }
}




