load_files_server <- function(input, output, session, state) {

  #== == == == == == == == == == == == == == == == ==
  #===== HELPERS ====================================
  #== == == == == == == == == == == == == == == == ==

  # Locate a built-in annotation file for the given annotation name.
  find_annotation_path <- function(annotation) {
    path <- system.file("annotations", paste0(annotation, ".rds"), package = "RNASeqApp")
    if (nzchar(path)) {
      path
    } else {
      local_path <- file.path("inst", "annotations", paste0(annotation, ".rds"))
      if (file.exists(local_path)) local_path else NULL
    }
  }

  # Extract all contrasts from SE rowData into a long-format data frame.
  extract_contrasts <- function(se) {
    rd       <- as.data.frame(SummarizedExperiment::rowData(se))
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

  # Merge annotation data frame into the contrast data frame.
  # Detects the join column automatically (ensembl_gene_id or Geneid).
  # Does NOT create a symbol column — that is handled via apply_symbol().
  # Returns the merged data frame, or NULL if no usable join column is found.
  merge_annotation <- function(de_df, annot) {
    if ("ensembl_gene_id" %in% colnames(annot)) {
      dplyr::left_join(de_df, annot, by = c("Geneid" = "ensembl_gene_id"))
    } else if ("Geneid" %in% colnames(annot)) {
      dplyr::left_join(de_df, annot, by = "Geneid")
    } else {
      NULL
    }
  }

  # Create (or update) the 'symbol' column in a data frame using primary and
  # optional fallback columns.  When the primary value is NA or empty the
  # fallback is used.  Returns the data frame unchanged if primary_col is not
  # present.
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

  # Determine default primary and secondary symbol columns from a set of
  # available column names.
  default_symbol_cols <- function(annot_cols) {
    primary <- {
      found <- intersect(c("gene", "hgnc_symbol"), annot_cols)
      if (length(found) > 0) found[1] else annot_cols[1]
    }
    secondary <- {
      found <- intersect(c("Geneid", "ensembl_gene_id"), annot_cols)
      if (length(found) > 0) found[1] else "none"
    }
    list(primary = primary, secondary = secondary)
  }

  # Build the two-dropdown symbol-column selection UI from annotation columns.
  build_symbol_select_ui <- function(annot_cols, primary_sel, secondary_sel) {
    secondary_choices <- c("none" = "none", setNames(annot_cols, annot_cols))

    tagList(
      tags$hr(style = "margin: 6px 0;"),
      tags$p(icon("tag"), strong("Select symbol columns:"),
             style = "margin: 4px 0 2px 0; font-size: 0.9em;"),
      selectInput("symbolPrimaryCol",
                  "Primary symbol column:",
                  choices  = annot_cols,
                  selected = primary_sel),
      selectInput("symbolSecondaryCol",
                  "Fallback column (used when primary is NA):",
                  choices  = secondary_choices,
                  selected = secondary_sel),
      actionButton("applySymbolColumns", "Apply", class = "btn-primary btn-sm")
    )
  }

  # Build the SE metadata info panel.
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

  #== == == == == == == == == == == == == == == == ==
  #===== EVENT HANDLERS =============================
  #== == == == == == == == == == == == == == == == ==

  ## ---- SE file upload: validation, contrast extraction, and initial annotation -----
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

    state$se_obj(se)

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
    rd           <- as.data.frame(SummarizedExperiment::rowData(se))
    varpart_cols <- grep("^varpart_", colnames(rd), value = TRUE)
    if (length(varpart_cols) > 0) {
      vp_mat <- as.data.frame(rd[, varpart_cols, drop = FALSE])
      colnames(vp_mat) <- sub("^varpart_", "", colnames(vp_mat))
      rownames(vp_mat) <- rownames(se)
      state$varpart_obj(list(varPart = vp_mat))
    }

    # Determine organism and annotation from SE metadata.
    se_organism <- tryCatch(metadata(se)$organism, error = function(e) NULL)
    organism    <- if (!is.null(se_organism) && nzchar(trimws(se_organism))) se_organism else NULL
    
    state$se_organism(organism)

    se_annotation <- tryCatch(metadata(se)$annotation, error = function(e) NULL)
    annotation    <- if (!is.null(se_annotation) && nzchar(trimws(se_annotation))) se_annotation else NULL

    annot_status  <- "none"
    annot_message <- NULL
    annot_cols    <- NULL
    sym_primary   <- NULL
    sym_secondary <- NULL

    if (!is.null(annotation)) {
      annot_path <- find_annotation_path(annotation)

      if (!is.null(annot_path)) {
        annot <- tryCatch(
          readRDS(annot_path),
          error = function(e) {
            showNotification(paste("Annotation load error:", e$message), type = "error")
            NULL
          }
        )
        if (!is.null(annot)) {
          merged <- merge_annotation(de_df, annot)
          if (!is.null(merged)) {
            # Determine and auto-apply default symbol columns.
            annot_cols    <- colnames(annot)
            defaults      <- default_symbol_cols(annot_cols)
            sym_primary   <- defaults$primary
            sym_secondary <- defaults$secondary

            de_df         <- apply_symbol(merged, sym_primary, sym_secondary)
            annot         <- apply_symbol(annot,  sym_primary, sym_secondary)
            state$annotation_df(annot)

            annot_status  <- "loaded"
            annot_message <- paste("Annotation loaded:", basename(annot_path))
          } else {
            annot_status  <- "no_join_col"
            annot_message <- paste0(
              "Annotation '", basename(annot_path), "' has no recognized join column ",
              "(expected 'Geneid' or 'ensembl_gene_id'). Upload a different annotation file."
            )
          }
        }
      } else {
        annot_status  <- "not_found"
        annot_message <- paste("No built-in annotation found for:", annotation)
      }
    }

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
                 style = "color: green; padding-left: 15px; margin: 2px 0;"),
          tags$div(style = "padding-left: 15px;",
                   build_symbol_select_ui(annot_cols, sym_primary, sym_secondary))
        )
      } else {
        tagList(
          info_panel,
          tags$p(icon("exclamation-triangle"),
                 HTML(
                   if (annot_status == "no_join_col") annot_message
                   else if (annot_status == "not_found")
                     paste0(annot_message,
                            ". Upload an annotation file or fix the annotation name in your SE object metadata.")
                   else
                     "Annotation not found in SE metadata. Upload an annotation file or add the annotation to your SE object metadata."
                 ),
                 style = "color: orange; padding-left: 15px;"),
          fileInput("annotFile", "Upload Annotation (.rds)", accept = ".rds")
        )
      }
    })
  })

  ## ---- Custom annotation upload -----
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

    se_current <- state$se_obj()
    organism   <- tryCatch(metadata(se_current)$organism, error = function(e) NULL)
    organism   <- if (!is.null(organism) && nzchar(trimws(organism))) organism else NULL

    merged <- merge_annotation(state$de_df(), annot)

    if (!is.null(merged)) {
      annot_cols    <- colnames(annot)
      defaults      <- default_symbol_cols(annot_cols)
      sym_primary   <- defaults$primary
      sym_secondary <- defaults$secondary

      merged <- apply_symbol(merged, sym_primary, sym_secondary)
      annot  <- apply_symbol(annot,  sym_primary, sym_secondary)
      state$annotation_df(annot)
      state$de_df(merged)

      output$annotationStatus <- renderUI({
        tagList(
          tags$div(style = "padding-left: 15px; margin-bottom: 8px;",
                   tags$p(style = "color: steelblue; margin: 0;",
                          icon("info-circle"), strong("SE loaded")),
                   se_info_ui(se_current, organism = organism)),
          tags$p(icon("check-circle"), "Custom annotation loaded successfully.",
                 style = "color: green; padding-left: 15px;"),
          tags$div(style = "padding-left: 15px;",
                   build_symbol_select_ui(annot_cols, sym_primary, sym_secondary))
        )
      })
      showNotification("Custom annotation loaded successfully.", type = "message")
    } else {
      output$annotationStatus <- renderUI({
        tagList(
          tags$div(style = "padding-left: 15px; margin-bottom: 8px;",
                   tags$p(style = "color: steelblue; margin: 0;",
                          icon("info-circle"), strong("SE loaded")),
                   se_info_ui(se_current, organism = organism)),
          tags$p(
            icon("exclamation-triangle"),
            HTML(paste0(
              "Annotation has no recognized join column ",
              "(<strong>Geneid</strong> or <strong>ensembl_gene_id</strong>). ",
              "Please upload an annotation file that contains one of these columns."
            )),
            style = "color: red; padding-left: 15px;"
          ),
          fileInput("annotFile", "Upload Annotation (.rds)", accept = ".rds")
        )
      })
      showNotification(
        "Annotation could not be joined: no 'Geneid' or 'ensembl_gene_id' column found.",
        type = "error"
      )
    }
  })

  ## ---- Symbol-column selection (apply button) -----
  observeEvent(input$applySymbolColumns, {
    req(input$symbolPrimaryCol, state$de_df())

    primary_col   <- input$symbolPrimaryCol
    secondary_col <- input$symbolSecondaryCol  # may be "none"

    if (!primary_col %in% colnames(state$de_df())) {
      showNotification(paste("Column not found in data:", primary_col), type = "error")
      return()
    }

    state$de_df(apply_symbol(state$de_df(), primary_col, secondary_col))

    # Mirror the symbol column in annotation_df for full-annotation downloads.
    annot <- state$annotation_df()
    if (!is.null(annot) && primary_col %in% colnames(annot)) {
      state$annotation_df(apply_symbol(annot, primary_col, secondary_col))
    }

    showNotification(
      paste0("Symbol column set to '", primary_col, "'",
             if (!is.null(secondary_col) && secondary_col != "none")
               paste0(" (fallback: '", secondary_col, "')") else ""),
      type = "message"
    )
  })

  #== == == == == == == == == == == == == == == == ==
  #===== REACTIVES ==================================
  #== == == == == == == == == == == == == == == == ==

  ## ---- Filtered DE table with global cutoffs -----
  filtered_de_df <- reactive({
    req(state$de_df())
    df <- state$de_df()
    if (!"log2FC" %in% colnames(df)) df$log2FC <- NA_real_
    if (!"padj"   %in% colnames(df)) df$padj   <- NA_real_
    padj_cut <- if (!is.null(input$global_padj_cutoff)   && !is.na(input$global_padj_cutoff))   input$global_padj_cutoff   else 0.05
    lfc_cut  <- if (!is.null(input$global_log2FC_cutoff) && !is.na(input$global_log2FC_cutoff)) input$global_log2FC_cutoff else 2
    df |>
      dplyr::filter(!is.na(padj) & !is.na(log2FC) & padj <= padj_cut & abs(log2FC) >= lfc_cut)
  })
  state$filtered_de_df <- filtered_de_df
}
