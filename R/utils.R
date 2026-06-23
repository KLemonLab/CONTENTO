#' Get base filename for downloads
#' 
#' Extracts a meaningful base name for downloaded files from the SE object
#' or uploaded filename.
#' 
#' @param input Shiny input object
#' @param state App state list containing reactive values
#' @return Character string to use as base filename
#' @export
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
#' @export
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
#' @export
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
#' @export
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


#' Add DE and regulation flags
#'
#' Adds two columns to a data frame of differential expression results:
#' \itemize{
#'   \item \code{regulated}: "up", "down", or NA
#'   \item \code{DE}: TRUE/FALSE
#' }
#' A row is considered DE if \code{padj < padj_cut} and
#' \code{|log2FC| > lfc_cut}.
#' 
#' @param df Data frame with \code{padj} and log2 fold-change column
#' @param lfc_cut Numeric log2FC threshold
#' @param padj_cut Numeric adjusted p-value threshold
#' @param lfc_col Column name for log2FC (default: "log2FC")
#' @return Data frame with added \code{regulated} and \code{DE} columns
#' @export
add_de_flags <- function(df, lfc_cut, padj_cut, lfc_col = "log2FC") {
  # Fall back to log2FC if chosen column doesn't exist
  if (!lfc_col %in% colnames(df)) lfc_col <- "log2FC"
  df |>
    mutate(
      regulated = case_when(
        !is.na(padj) & !is.na(.data[[lfc_col]]) & padj < padj_cut & .data[[lfc_col]] >  lfc_cut ~ "up",
        !is.na(padj) & !is.na(.data[[lfc_col]]) & padj < padj_cut & .data[[lfc_col]] < -lfc_cut ~ "down",
        TRUE ~ NA_character_
      ),
      DE = !is.na(regulated)
    )
}

#' Determine the default symbol column from available annotation columns
#'
#' Checks a prioritised list of common gene name column names and returns
#' the first match found, falling back to the first available column.
#'
#' @param annot_cols Character vector of column names from the annotation data frame
#' @return Character string with the name of the best candidate symbol column
#' @export
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
#' @export
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
#' @param lfc_col Character string naming the log2 fold-change column to rank by (default: "log2FC")
#' @return Character vector of Geneids ordered by abs(log2FC)
#' @export
get_top_de_genes <- function(selected_data, top_n, lfc_col = "log2FC") {
  if (!lfc_col %in% colnames(selected_data)) lfc_col <- "log2FC"
  selected_data |>
    filter(DE) |>
    arrange(desc(abs(.data[[lfc_col]]))) |>
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
#' @export
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
#' @param lfc_col Character string naming the log2 fold-change column in `df` (default: "log2FC").
#' @return A data frame in wide format with pairs of `log2FC_<contrast>` and `DE_<contrast>` for each contrast.
#' @export
build_compare_table <- function(df, contrasts, lfc_col = "log2FC") {
  if (is.null(df)) return(NULL)
  
  has_symbol <- "symbol" %in% colnames(df)
  id_cols <- if (has_symbol) c("Geneid", "symbol") else "Geneid"
  
  # Fall back if column doesn't exist
  if (!lfc_col %in% colnames(df)) lfc_col <- "log2FC"
  
  wide <- df |>
    select(all_of(c(id_cols, "contrast", lfc_col, "DE"))) |>
    (\(x) { names(x)[names(x) == lfc_col] <- "FC"; x })() |>
    pivot_wider(
      id_cols = all_of(id_cols),
      names_from = contrast,
      values_from = c("FC", "DE"),
      names_sep = "_"
    ) |>
    mutate(
      across(starts_with(paste0("FC", "_")), ~ round(., 2)),
      across(starts_with("DE_"), ~ tidyr::replace_na(., FALSE))
    ) |>
    filter(rowSums(across(starts_with("DE_"))) > 0)
  
  ordered_cols <- id_cols
  for (ct in contrasts) {
    ordered_cols <- c(ordered_cols, paste0("FC_", ct), paste0("DE_", ct))
  }
  
  wide |>
    select(any_of(ordered_cols)) |>
    arrange(Geneid)
}

#' Get the msigdbr db_species code for a given organism name
#'
#' Returns "MM" for Mus musculus and "HS" for all other species supported by
#' msigdbr (non-human/mouse species use human gene sets with ortholog mapping).
#' Returns NULL if the organism is not in the msigdbr species list.
#'
#' An optional `override` argument allows you to manually specify the db_species
#' code. If provided as "HS" or "MM", it takes precedence over automatic detection
#' based on the organism name.
#'
#' @param organism Character; organism name as stored in metadata(se)$organism,
#'   e.g. "Homo sapiens", "Mus musculus", "Rattus norvegicus"
#' @param override Optional character; manually force the db_species value.
#'   Must be either "HS" or "MM". If supplied and valid, this value is returned
#'   directly, bypassing organism-based logic.
#' @return Character "HS", "MM", or NULL
#' @export
get_msigdbr_db_species <- function(organism, override = NULL) {
  if (!is.null(override) && override %in% c("HS", "MM")) return(override)
  if (is.null(organism) || !nzchar(trimws(organism))) return(NULL)
  valid <- tryCatch(msigdbr::msigdbr_species()$species_name, error = function(e) character(0))
  if (!organism %in% valid) return(NULL)
  if (organism == "Mus musculus") "MM" else "HS"
}


#' Get available MSigDB collections for a given db_species
#'
#' Calls msigdbr_collections() and returns a tidy data frame of available
#' collections and subcollections for building the UI picker.
#'
#' @param db_species Character; "HS" or "MM"
#' @return Data frame with columns gs_collection, gs_subcollection,
#'   gs_collection_name, num_genesets; or NULL on error
#' @export
get_msigdbr_collections <- function(db_species) {
  tryCatch(
    msigdbr::msigdbr_collections(db_species = db_species),
    error = function(e) NULL
  )
}


#' Detect the Ensembl gene ID column in an annotation data frame
#'
#' Scans column values looking for Ensembl-style IDs (ENSG..., ENSMUSG...,
#' ENSRNOG..., etc.). Returns the first matching column name, or NULL.
#'
#' @param annot_df Data frame of gene annotations
#' @return Character column name, or NULL if none detected
#' @export
detect_ensembl_col <- function(annot_df) {
  if (!is.data.frame(annot_df) || nrow(annot_df) == 0) return(NULL)
  for (col in colnames(annot_df)) {
    vals <- as.character(annot_df[[col]])
    # Sample up to 20 non-NA values for speed
    sample_vals <- vals[!is.na(vals) & nzchar(vals)]
    if (length(sample_vals) == 0) next
    sample_vals <- head(sample_vals, 20)
    if (any(grepl("^ENS[A-Z]*G[0-9]{11}", sample_vals))) return(col)
  }
  NULL
}

#' Get all potential gsea columns available in dataframe
#'
#' Returns all column names except common metadata columns (Geneid, symbol, locus_tag, etc.)
#' These are candidate columns for functional enrichment analysis.
#'
#' @param annot_df Data frame containing annotations
#' @return Named list where names are column names and values are also column names
#'         (suitable for use in selectInput choices)
#' @export
get_gsea_columns <- function(annot_df) {
  if (!is.data.frame(annot_df) || nrow(annot_df) == 0) {
    return(list())
  }
  
  # Exclude metadata/ID/structural columns not useful for enrichment
  exclude_cols <- c(
    # Structural columns
    "Geneid", "symbol", "seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attributes", "biotype",
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
#' Handles multiple values per gene separated by \code{!!!}.
#' This works with any organism's annotation data frame.
#'
#' @param annot_df Data frame containing gene annotations with a Geneid column
#' @param column_name Character string specifying which annotation column to use
#' @return Named list where each element is a character vector of Geneids
#' @export
build_annotation_genesets <- function(annot_df, column_name) {
  if (!column_name %in% colnames(annot_df)) return(list())
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

#' Load gene sets as a named list of gene ID vectors
#'
#' Single entry point for geneset loading across the app. Returns a named list
#' suitable for direct use in fgseaMultilevel() and geseca().
#'
#' Two source types are supported:
#' \itemize{
#'   \item \strong{MSigDB}: pass \code{msigdbr_species}, \code{db_species},
#'     \code{gs_collection}, and optionally \code{gs_subcollection}.
#'     Gene sets are keyed by Ensembl gene ID (the \code{ensembl_gene} column
#'     from msigdbr), so SE rownames must be Ensembl IDs.
#'   \item \strong{Annotation column}: pass \code{annotation_df} and
#'     \code{annotation_source}. The column may contain multiple values per gene
#'     separated by \code{!!!}.
#' }
#'
#' @param source_type    Character; "msigdb" or "annotation"
#' @param msigdbr_species Character; full species name for msigdbr(), e.g.
#'   "Homo sapiens". Only used when source_type == "msigdb".
#' @param db_species     Character; "HS" or "MM" for msigdbr db_species.
#'   Only used when source_type == "msigdb".
#' @param gs_collection  Character; MSigDB collection code (e.g. "H", "C2").
#'   Only used when source_type == "msigdb".
#' @param gs_subcollection Character; MSigDB subcollection (e.g. "CP:REACTOME"),
#'   or "" / NULL to omit. Only used when source_type == "msigdb".
#' @param annotation_df  Data frame; annotation with Geneid column.
#'   Only used when source_type == "annotation".
#' @param annotation_source Character; column name for gene set membership.
#'   Only used when source_type == "annotation".
#' @return Named list of character vectors (Geneids per pathway).
#'   Returns an empty list if inputs are invalid or nothing is found.
#' @export
load_genesets <- function(source_type,
                          msigdbr_species  = NULL,
                          db_species       = NULL,
                          gs_collection    = NULL,
                          gs_subcollection = NULL,
                          annotation_df    = NULL,
                          annotation_source = NULL) {
  
  if (source_type == "annotation") {
    if (is.null(annotation_df) || is.null(annotation_source)) return(list())
    return(build_annotation_genesets(annotation_df, annotation_source))
  }
  
  if (source_type == "msigdb") {
    if (is.null(msigdbr_species) || is.null(db_species) || is.null(gs_collection)) return(list())
    use_sub <- !is.null(gs_subcollection) && nzchar(gs_subcollection)
    genesets <- if (use_sub) {
      msigdbr::msigdbr(db_species = db_species, species = msigdbr_species,
                       collection = gs_collection, subcollection = gs_subcollection)
    } else {
      msigdbr::msigdbr(db_species = db_species, species = msigdbr_species,
                       collection = gs_collection)
    }
    if (is.null(genesets) || nrow(genesets) == 0) return(list())
    genesets |>
      split(genesets$gs_name) |>
      lapply(function(x) x$ensembl_gene)
  }
}

#' Resolve the active gene-set source into load_genesets() arguments
#'
#' Shared logic for both Explore-by-Contrast (single source set) and
#' Compare-Contrast (separate GSEA/GESECA source sets). Inspects available
#' MSigDB/annotation sources from state, determines which one is active
#' (via the radio input when both are available), and returns a ready-to-use
#' argument list for load_genesets(), plus a display label.
#'
#' @param input            Shiny input object
#' @param state            App state list (reactiveVals)
#' @param organism         Character or NULL; full organism name (for msigdbr_species)
#' @param source_input_id  Shiny input ID for the radio source selector
#' @param collection_id    Shiny input ID for the MSigDB collection selectInput
#' @param subcollection_id Shiny input ID for the MSigDB subcollection textInput
#' @param annot_source_id  Shiny input ID for the annotation column selectInput
#' @return Named list with all load_genesets() arguments plus `label`,
#'   or NULL if no source is available / required inputs aren't ready yet.
#'   Use req() on the result, then strip `label` before do.call(load_genesets, .).
#' @export
resolve_gsea_source <- function(input, state, organism,
                                source_input_id, collection_id,
                                subcollection_id, annot_source_id) {
  db_sp      <- state$db_species()
  ensembl_c  <- state$ensembl_col()
  avail_cols <- state$available_gsea_columns()
  
  has_msigdb <- !is.null(db_sp) && !is.null(ensembl_c)
  has_annot  <- length(avail_cols) > 0
  
  use_source <- if (has_msigdb && has_annot) {
    input[[source_input_id]] %||% "msigdb"
  } else if (has_msigdb) {
    "msigdb"
  } else if (has_annot) {
    "annotation"
  } else {
    return(NULL)
  }
  
  if (use_source == "msigdb") {
    if (is.null(input[[collection_id]])) return(NULL)
    list(
      source_type      = "msigdb",
      msigdbr_species  = organism,
      db_species       = db_sp,
      gs_collection    = input[[collection_id]],
      gs_subcollection = input[[subcollection_id]],
      label = paste0(input[[collection_id]],
                     if (nzchar(input[[subcollection_id]] %||% ""))
                       paste0("_", input[[subcollection_id]]))
    )
  } else {
    if (is.null(state$annotation_df()) || is.null(input[[annot_source_id]])) return(NULL)
    list(
      source_type       = "annotation",
      annotation_df     = state$annotation_df(),
      annotation_source = input[[annot_source_id]],
      label             = input[[annot_source_id]]
    )
  }
}

#' Build a standardised download filename for Shiny downloadHandler
#'
#' Derives the base name from the SE object metadata or uploaded filename,
#' then appends contrast, type, filters, and suffix segments.
#'
#' @param input    Shiny input object
#' @param state    App state list containing reactive values
#' @param type     High-level label, e.g. "all_genes", "filtered", "GSEA"; NULL to omit
#' @param ext      File extension without dot (default "csv")
#' @param contrast Character scalar or vector of contrast names; NULL to omit
#' @param filters  Named list of filter values, e.g. list(FC = 1.5, FDR = 0.05)
#' @param suffix   Extra free-form suffix inserted before the extension; NULL to omit
#' @return A single filename string
#' @export
build_download_filename <- function(input, state,
                                    type     = NULL,
                                    ext      = "csv",
                                    contrast = NULL,
                                    filters  = NULL,
                                    suffix   = NULL) {
  base <- tryCatch(
    metadata(state$se_obj())$dataset_name,
    error = function(e) NULL
  )
  
  parts <- list(
    base = if (!is.null(base) && nzchar(trimws(base))) {
      base
    } else if (!is.null(input$seFile$name)) {
      tools::file_path_sans_ext(basename(input$seFile$name))
    } else {
      "RNASeq_results"
    },
    
    contrast = if (!is.null(contrast) && any(nzchar(contrast)))
      gsub("[^A-Za-z0-9._-]+", "__", paste(contrast, collapse = "-")),
    
    type = if (!is.null(type) && nzchar(type)) type,
    
    filters = if (!is.null(filters) && length(filters) > 0)
      paste(mapply(function(k, v) paste0(k, gsub("\\.", "p", as.character(v))),
                   names(filters), filters), collapse = "__"),
    
    suffix = if (!is.null(suffix) && nzchar(suffix))
      gsub("[^A-Za-z0-9._-]+", "-", as.character(suffix))
  )
  
  paste0(paste(Filter(Negate(is.null), parts), collapse = "__"), ".", ext)
}

