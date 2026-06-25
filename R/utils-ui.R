#' SE metadata summary panel
#'
#' Renders a compact bullet list of key SE properties for display in the
#' sidebar after a file is loaded.
#'
#' @param se            SummarizedExperiment object
#' @param organism      Character; organism string from metadata, or NULL
#' @param annotation    Character; annotation name from metadata, or NULL
#' @return A tagList containing an unordered list
#' @export
se_info_ui <- function(se, organism = NULL, annotation = NULL) {
  meta        <- tryCatch(metadata(se), error = function(e) list())
  n_contrasts <- length(meta$contrasts)
  
  info_rows <- tagList(
    tags$li(icon("dna"),         strong("Genes: "),     nrow(se)),
    tags$li(icon("vials"),       strong("Samples: "),   ncol(se)),
    tags$li(icon("layer-group"), strong("Contrasts: "),  n_contrasts),
    if (!is.null(organism))
      tags$li(icon("bug"),  strong("Organism: "),    organism),
    if (!is.null(annotation))
      tags$li(icon("book"), strong("Annotation: "),  annotation)
  )
  
  tagList(
    tags$ul(style = "list-style: none; padding-left: 15px; margin: 5px 0;",
            info_rows)
  )
}


empty_state_msg <- function() {
  div(class = "alert alert-warning",
      icon("exclamation-triangle"),
      " Please upload a SE file to begin.")
}

#' Gene symbol / label column picker
#'
#' Builds a selectInput pre-populated with annotation columns, excluding
#' structural / coordinate columns, with the detected primary column selected.
#'
#' @param annot_cols  Character vector of all annotation column names
#' @param primary_sel Character; the column to pre-select
#' @return A tagList containing a selectInput
#' @export
build_symbol_select_ui <- function(annot_cols, primary_sel) {
  exclude <- c("seqname", "source", "feature", "start", "end", "score", "strand", 
               "frame", "attributes", "Geneid")
  filtered_cols <- setdiff(annot_cols, exclude)
  ordered_cols  <- c(primary_sel, setdiff(filtered_cols, primary_sel))
  tagList(
    selectInput("symbolPrimaryCol",
                "Select Gene symbol/label column:",
                choices  = ordered_cols,
                selected = primary_sel)
  )
}

#' MSigDB gene set picker
#'
#' Generates a single selectInput listing every queryable MSigDB gene set
#' for \code{db_species}, grouped by top-level collection (optgroups).
#' Collections without subcollections (e.g. C1, H) appear as a single choice;
#' collections with subcollections (e.g. C2, C5) list each subcollection by
#' its descriptive name. Defaults to Hallmark (H) when available.
#'
#' @param collection_input_id Shiny input ID for the gene set selectInput
#' @param db_species          Character; "HS" or "MM"
#' @param organism_label      Character; display name shown above the picker,
#'   e.g. "Homo sapiens" or "Rattus norvegicus". NULL to omit.
#' @return A tagList of Shiny UI elements, or an error alert div on failure
#' @export
msigdb_controls_ui <- function(collection_input_id,
                               db_species,
                               organism_label = NULL) {
  cols <- get_msigdbr_collections(db_species)
  
  if (is.null(cols)) {
    return(div(
      class = "alert alert-warning",
      icon("exclamation-triangle"),
      "Could not load MSigDB collections. Check msigdbr installation."
    ))
  }
  
  cols <- cols |>
    dplyr::mutate(
      has_sub = !is.na(gs_subcollection) & gs_subcollection != "",
      # value encodes both pieces needed later to query msigdbr()
      value   = ifelse(has_sub, paste0(gs_collection, "|", gs_subcollection), gs_collection),
      # label is always the descriptive name (collection name for top-level,
      # subcollection's own name for nested entries)
      label   = gs_collection_name
    )
  
  # Build optgroup-style nested choices: named list of named vectors
  grouped_choices <- cols |>
    dplyr::group_by(gs_collection) |>
    dplyr::group_map(~ setNames(.x$value, .x$label)) |>
    setNames(unique(cols$gs_collection))
  
  default_val <- cols$value[cols$gs_collection == "H"][1] %||% cols$value[1]
  
  # Show the queried db_species only when it differs from a native match —
  # i.e. only when ortholog mapping is actually happening.
  is_native <- (db_species == "HS" && identical(organism_label, "Homo sapiens")) ||
    (db_species == "MM" && identical(organism_label, "Mus musculus"))
  
  db_note <- if (is_native) {
    NULL
  } else {
    db_label <- if (db_species == "MM") "Mus musculus (MM)" else "Homo sapiens (HS)"
    tags$span(
      paste0(" \u2014 ortholog-mapped from ", db_label),
      style = "font-style: italic; color: #555;"
    )
  }
  
  organism_note <- if (!is.null(organism_label)) {
    tags$p(
      icon("globe"),
      strong("Select from MSigDB for: "),
      organism_label,
      db_note,
      style = "margin-bottom:4px; font-size: 0.95em;"
    )
  }
  
  tagList(
    organism_note,
    
    div(
      style = "margin-bottom:6px;",
      selectInput(
        collection_input_id,
        label = NULL,
        choices  = grouped_choices,
        selected = default_val,
        width = "100%"
      )
    )
  )
}


#' Source-aware GSEA controls panel
#'
#' Builds the full gene set source picker used in both the single-contrast
#' (Explore by Contrast) and multi-contrast (Compare Contrast) tabs.
#'
#' Behaviour:
#' \itemize{
#'   \item \strong{Both MSigDB and annotation available:} radio button to choose
#'     source, then the relevant sub-controls.
#'   \item \strong{Only MSigDB:} MSigDB controls shown directly.
#'   \item \strong{Only annotation:} annotation column picker shown directly.
#'   \item \strong{Neither:} an informative warning alert.
#' }
#'
#' All reactive values must be resolved by the caller before passing them in.
#'
#' @param db_species        Character or NULL; from \code{state$db_species()}
#' @param ensembl_col       Character or NULL; from \code{state$ensembl_col()}
#' @param avail_cols        Named list; from \code{state$available_gsea_columns()}
#' @param organism_name     Character or NULL; from \code{organism()}
#' @param source_input_id   Shiny input ID for the radio source selector
#' @param collection_id     Shiny input ID for the MSigDB collection selectInput.
#'   The selected value encodes both collection and optional subcollection as
#'   \code{"COLLECTION|SUBCOLLECTION"} (parsed by \code{resolve_gsea_source}).
#'   No separate subcollection input is used.
#' @param annot_source_id   Shiny input ID for the annotation column selectInput
#' @param current_source    Character; current value of \code{input[[source_input_id]]},
#'   used to preserve selection across re-renders. Pass \code{NULL} to default to "msigdb".
#' @return A tagList of Shiny UI elements
#' @export
gsea_source_ui <- function(db_species,
                           ensembl_col,
                           avail_cols,
                           organism_name,
                           source_input_id,
                           collection_id,
                           annot_source_id,
                           current_source = NULL) {
  
  has_msigdb <- !is.null(db_species) && !is.null(ensembl_col)
  has_annot  <- length(avail_cols) > 0
  
  if (!has_msigdb && !has_annot) {
    msg <- if (is.null(organism_name)) {
      "Upload an SE file to enable functional enrichment."
    } else if (!is.null(db_species) && is.null(ensembl_col)) {
      paste0(
        "MSigDB requires an Ensembl ID column in your annotation file. ",
        "Upload an annotation with Ensembl IDs, or add functional columns ",
        "for annotation-based GSEA."
      )
    } else {
      paste0(
        "No gene sets available for organism: \u2018", organism_name, "\u2019. ",
        "Upload an annotation file with functional columns to enable GSEA."
      )
    }
    return(div(class = "alert alert-warning", icon("exclamation-triangle"), msg))
  }
  
  tagList(
    
    if (has_msigdb && has_annot)
      div(
        style = "display:flex; align-items:center; gap:12px; margin-bottom:6px;",
        
        tags$label(
          "Gene sets:",
          style = "margin:0; font-weight:600;"
        ),
        
        radioButtons(
          source_input_id,
          label = NULL,
          choices = c(
            "MSigDB" = "msigdb",
            "Annotation" = "annotation"
          ),
          selected = current_source %||% "msigdb",
          inline = TRUE
        )
      ),
    
    # MSigDB
    if (has_msigdb)
      conditionalPanel(
        condition = if (has_annot)
          paste0("input.", source_input_id, " == 'msigdb'")
        else "true",
        
        msigdb_controls_ui(
          collection_input_id = collection_id,
          db_species = db_species,
          organism_label = organism_name
        )
      ),
    
    # Annotation (match MSigDB style)
    if (has_annot)
      conditionalPanel(
        condition = if (has_msigdb)
          paste0("input.", source_input_id, " == 'annotation'")
        else "true",
        
        tags$div(
          icon("tags"),   
          strong("Select from Annotation"),
          style = "margin-bottom:4px; font-size: 0.95em;"
        ),
        
        div(
          style = "margin-bottom:6px;",
          selectInput(
            annot_source_id,
            label = NULL,
            choices  = avail_cols,
            selected = avail_cols[[1]],
            width = "100%"
          )
        )
      )
  )
}