load_files_server <- function(input, output, session, state) {

  # ---- Helpers ----------------------------------------------------------------

  # Locate a built-in annotation file for the given organism name.
  find_annotation_path <- function(organism) {
    # Try installed-package path first.
    path <- system.file("annotations",
                        paste0(organism, "_annot.rds"),
                        package = "RNASeqApp")
    if (nchar(path) > 0) return(path)

    # Fall back to relative path when running the app directly from source.
    local_path <- file.path("inst", "annotations",
                            paste0(organism, "_annot.rds"))
    if (file.exists(local_path)) return(local_path)

    return(NULL)
  }

  # Extract all contrasts from SE rowData into a long-format data frame.
  extract_contrasts <- function(se) {
    rd       <- as.data.frame(SummarizedExperiment::rowData(se))
    gene_ids <- rownames(se)

    lfc_cols <- grep("^log2FoldChange_", colnames(rd), value = TRUE)
    if (length(lfc_cols) == 0) {
      stop("No contrast columns found in rowData. ",
           "Expected columns like 'log2FoldChange_ContrastName'.")
    }

    contrast_names <- sub("^log2FoldChange_", "", lfc_cols)

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
  # If organism is provided it is included in the list with the same style as
  # genes/samples/contrasts.
  se_info_ui <- function(se, organism = NULL) {
    meta  <- tryCatch(SummarizedExperiment::metadata(se), error = function(e) list())
    rd    <- as.data.frame(SummarizedExperiment::rowData(se))
    n_contrasts <- length(grep("^log2FoldChange_", colnames(rd), value = TRUE))

    # Core stats always shown
    info_rows <- tagList(
      tags$li(icon("dna"),         strong("Genes: "),    nrow(se)),
      tags$li(icon("vials"),       strong("Samples: "),  ncol(se)),
      tags$li(icon("layer-group"), strong("Contrasts: "), n_contrasts),
      if (!is.null(organism))
        tags$li(icon("bug"), strong("Organism: "), organism)
    )

    # Add any named metadata fields (skip 'organism', already shown above)
    extra <- meta[setdiff(names(meta), "organism")]
    extra_items <- Filter(Negate(is.null), lapply(names(extra), function(k) {
      val <- extra[[k]]
      if (is.character(val) || is.numeric(val)) {
        tags$li(strong(paste0(k, ": ")), as.character(val))
      }
    }))

    tagList(
      tags$ul(style = "list-style: none; padding-left: 15px; margin: 5px 0;",
              info_rows,
              if (length(extra_items) > 0) extra_items)
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

    # Determine organism from SE metadata only.
    se_organism <- tryCatch(
      metadata(se)$organism,
      error = function(e) NULL
    )
    organism <- if (!is.null(se_organism) && nzchar(trimws(se_organism))) {
      se_organism
    } else {
      NULL
    }

    if (is.null(organism)) {
      state$de_df(de_df)
      output$annotationStatus <- renderUI({
        tagList(
          tags$div(style = "padding-left: 15px; margin-bottom: 8px;",
                   tags$p(style = "color: steelblue; margin: 0;",
                          icon("info-circle"), strong("SE loaded")),
                   se_info_ui(se)),
          tags$p(icon("exclamation-triangle"),
                 HTML("Organism not found in SE metadata. Upload an annotation file or add the organism to your SE object metadata."),
                 style = "color: orange; padding-left: 15px;"),
          fileInput("annotFile", "Upload Annotation (.rds)", accept = ".rds")
        )
      })
      return(NULL)
    }

    # Try to load built-in annotation.
    annot_path <- find_annotation_path(organism)

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
        de_df <- merge_annotation(de_df, annot)
      }
      state$de_df(de_df)
      output$annotationStatus <- renderUI({
        tagList(
          tags$div(style = "padding-left: 15px; margin-bottom: 8px;",
                   tags$p(style = "color: steelblue; margin: 0;",
                          icon("info-circle"), strong("SE loaded")),
                   se_info_ui(se, organism = organism)),
          tags$p(icon("check-circle"),
                 paste("Annotation loaded:", basename(annot_path)),
                 style = "color: green; padding-left: 15px; margin: 2px 0;")
        )
      })
    } else {
      # No built-in annotation found -- show upload UI.
      state$de_df(de_df)
      output$annotationStatus <- renderUI({
        tagList(
          tags$div(style = "padding-left: 15px; margin-bottom: 8px;",
                   tags$p(style = "color: steelblue; margin: 0;",
                          icon("info-circle"), strong("SE loaded")),
                   se_info_ui(se, organism = organism)),
          tags$p(icon("exclamation-triangle"),
                 paste("No built-in annotation found for:", organism, ". Upload an annotation file or fix the organism name in your SE object metadata"),
                 style = "color: orange; padding-left: 15px;"),
          fileInput("annotFile", "Upload Annotation (.rds)", accept = ".rds")
        )
      })
    }
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