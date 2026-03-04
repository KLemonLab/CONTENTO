load_files_server <- function(input, output, session, state) {

  # ---- Helpers ----------------------------------------------------------------

  # Locate a built-in annotation file for the given annotation name.
  find_annotation_path <- function(annotation) {
    # Try installed-package path first.
    path <- system.file("annotations",
                        paste0(annotation, ".rds"),
                        package = "RNASeqApp")
    if (nchar(path) > 0) return(path)

    # Fall back to relative path when running the app directly from source.
    local_path <- file.path("inst", "annotations",
                            paste0(annotation, ".rds"))
    if (file.exists(local_path)) return(local_path)

    return(NULL)
  }

  # Extract all contrasts from SE rowData into a long-format data frame.
  extract_contrasts <- function(se) {
    rd       <- as.data.frame(SummarizedExperiment::rowData(se))
    gene_ids <- rownames(se)

    # Get contrast names from metadata (SOURCE OF TRUTH!)
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

  # Merge an annotation data frame into the contrast data frame.
  merge_annotation <- function(de_df, annot) {
    if (!"Geneid" %in% colnames(annot)) {
      showModal(modalDialog(
        title = "Annotation error",
        "Annotation file is missing required 'Geneid' column.",
        easyClose = TRUE,
        footer = NULL
      ))
      return(de_df)
    }

    # Store full annotation for export functionality.
    state$annotation_df(annot)

    dplyr::left_join(de_df, annot, by = "Geneid")
  }

  # Build an info panel from SE metadata and basic dimension info.
  # If organism and/or annotation are provided they are included in the list
  # with the same style as genes/samples/contrasts.
  se_info_ui <- function(se, organism = NULL, annotation = NULL) {
    meta  <- tryCatch(metadata(se), error = function(e) list())
    n_contrasts <- length(meta$contrasts)

    # Core stats always shown
    info_rows <- tagList(
      tags$li(icon("dna"),         strong("Genes: "),    nrow(se)),
      tags$li(icon("vials"),       strong("Samples: "),  ncol(se)),
      tags$li(icon("layer-group"), strong("Contrasts: "), n_contrasts),
      if (!is.null(organism))
        tags$li(icon("bug"), strong("Organism: "), organism),
      if (!is.null(annotation))
        tags$li(icon("book"), strong("Annotation: "), annotation)
    )

    tagList(
      tags$ul(style = "list-style: none; padding-left: 15px; margin: 5px 0;",
              info_rows)
    )
  }

  observeEvent(input$seFile, {
    req(input$seFile)

    se <- tryCatch(
      readRDS(input$seFile$datapath),
      error = function(e) {
        showModal(modalDialog(
          title = "File error",
          paste("Error reading file:", e$message),
          easyClose = TRUE,
          footer = NULL
        ))
        NULL
      }
    )
    req(se)

    if (!inherits(se, "SummarizedExperiment")) {
      showModal(modalDialog(
        title = "File error",
        "Uploaded file is not a SummarizedExperiment object.",
        easyClose = TRUE,
        footer = NULL
      ))
      return(NULL)
    }

    # Store SE (compatible with DESeq2/SummarizedExperiment accessor functions).
    state$dds_obj(se)

    # Extract contrasts from rowData.
    de_df <- tryCatch(
      extract_contrasts(se),
      error = function(e) {
        showModal(modalDialog(
          title = "Contrast extraction error",
          e$message,
          easyClose = TRUE,
          footer = NULL
        ))
        NULL
      }
    )
    req(de_df)

    # Extract variance partition if varpart_* columns are present.
    rd <- as.data.frame(SummarizedExperiment::rowData(se))
    varpart_cols <- grep("^varpart_", colnames(rd), value = TRUE)
    if (length(varpart_cols) > 0) {
      vp_mat <- as.data.frame(rd[, varpart_cols, drop = FALSE])
      colnames(vp_mat) <- sub("^varpart_", "", colnames(vp_mat))
      rownames(vp_mat) <- rownames(se)
      state$varpart_obj(list(varPart = vp_mat))
    }

    # Determine organism and annotation from SE metadata only.
    se_organism <- tryCatch(
      metadata(se)$organism,
      error = function(e) NULL
    )
    organism <- if (!is.null(se_organism) && nzchar(trimws(se_organism))) {
      se_organism
    } else {
      NULL
    }

    se_annotation <- tryCatch(
      metadata(se)$annotation,
      error = function(e) NULL
    )
    annotation <- if (!is.null(se_annotation) && nzchar(trimws(se_annotation))) {
      se_annotation
    } else {
      NULL
    }

    # Track annotation status for UI rendering after de_df is set.
    annot_status  <- "none"
    annot_message <- NULL

    if (!is.null(annotation)) {
      # Try to load built-in annotation.
      annot_path <- find_annotation_path(annotation)

      if (!is.null(annot_path)) {
        annot <- tryCatch(
          readRDS(annot_path),
          error = function(e) {
            showNotification(paste("Annotation load error:", e$message),
                             type = "error")
            NULL
          }
        )
        if (!is.null(annot)) {
          de_df         <- merge_annotation(de_df, annot)
          annot_status  <- "loaded"
          annot_message <- paste("Annotation loaded:", basename(annot_path))
        }
      } else {
        annot_status  <- "not_found"
        annot_message <- paste("No built-in annotation found for:", annotation)
      }
    }

    # Set de_df exactly once, after all annotation logic is complete.
    state$de_df(de_df)

    info_panel <- tags$div(
      style = "padding-left: 15px; margin-bottom: 8px;",
      tags$p(style = "color: steelblue; margin: 0;",
             icon("info-circle"), strong("SE loaded")),
      se_info_ui(se, organism = organism, annotation = annotation)
    )

    output$annotationStatus <- renderUI({
      if (annot_status == "loaded") {
        tagList(
          info_panel,
          tags$p(icon("check-circle"), annot_message,
                 style = "color: green; padding-left: 15px; margin: 2px 0;")
        )
      } else if (annot_status == "not_found") {
        tagList(
          info_panel,
          tags$p(icon("exclamation-triangle"),
                 HTML(paste0(annot_message, ". Upload an annotation file or fix the annotation name in your SE object metadata.")),
                 style = "color: orange; padding-left: 15px;"),
          fileInput("annotFile", "Upload Annotation (.rds)", accept = ".rds")
        )
      } else {
        # annot_status == "none": no annotation in SE metadata
        tagList(
          info_panel,
          tags$p(icon("exclamation-triangle"),
                 HTML("Annotation not found in SE metadata. Upload an annotation file or add the annotation to your SE object metadata."),
                 style = "color: orange; padding-left: 15px;"),
          fileInput("annotFile", "Upload Annotation (.rds)", accept = ".rds")
        )
      }
    })
  })

  # ---- Custom annotation upload -----------------------------------------------

  observeEvent(input$annotFile, {
    req(input$annotFile, state$de_df())

    annot <- tryCatch(
      readRDS(input$annotFile$datapath),
      error = function(e) {
        showModal(modalDialog(
          title = "Annotation error",
          paste("Error reading annotation file:", e$message),
          easyClose = TRUE,
          footer = NULL
        ))
        NULL
      }
    )
    req(annot)

    de_df <- merge_annotation(state$de_df(), annot)
    state$de_df(de_df)

    # Rebuild the info panel using the stored SE so metadata stays visible.
    se_current <- state$dds_obj()
    output$annotationStatus <- renderUI({
      tagList(
        tags$div(style = "padding-left: 15px; margin-bottom: 8px;",
                 tags$p(style = "color: steelblue; margin: 0;",
                        icon("info-circle"), strong("SE loaded")),
                 se_info_ui(se_current)),
        tags$p(icon("check-circle"), "Custom annotation loaded successfully.",
               style = "color: green; padding-left: 15px;")
      )
    })
    showNotification("Custom annotation loaded successfully.", type = "message")
  })
}