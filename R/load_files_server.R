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

    # Get contrast names from metadata
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
  merge_annotation <- function(de_df, annot, organism = NULL, symbol_column = NULL) {

    # Step 1: Check for required columns based on organism.
    id_column            <- NULL
    symbol_source_column <- NULL
    missing_columns      <- c()

    if (!is.null(organism) && organism == "Human") {
      # Human: need ensembl_gene_id AND hgnc_symbol
      if (!("ensembl_gene_id" %in% colnames(annot))) {
        missing_columns <- c(missing_columns, "ensembl_gene_id")
      } else {
        id_column <- "ensembl_gene_id"
      }

      if (!("hgnc_symbol" %in% colnames(annot))) {
        missing_columns <- c(missing_columns, "hgnc_symbol")
      } else {
        symbol_source_column <- "hgnc_symbol"
      }

    } else if (!is.null(organism) && organism == "Bacteria") {
      # Bacterial: need Geneid AND gene
      if (!("Geneid" %in% colnames(annot))) {
        missing_columns <- c(missing_columns, "Geneid")
      } else {
        id_column <- "Geneid"
      }

      if (!("gene" %in% colnames(annot))) {
        missing_columns <- c(missing_columns, "gene")
      } else {
        symbol_source_column <- "gene"
      }

    } else {
      # Unknown organism - try to auto-detect an ID column.
      if ("Geneid" %in% colnames(annot)) {
        id_column <- "Geneid"
      } else if ("ensembl_gene_id" %in% colnames(annot)) {
        id_column <- "ensembl_gene_id"
      } else {
        missing_columns <- c(missing_columns, "Geneid", "ensembl_gene_id")
      }
    }

    # If required columns are missing and user hasn't supplied a symbol column,
    # return NULL to trigger interactive column selection.
    if (length(missing_columns) > 0 && is.null(symbol_column)) {
      return(NULL)
    }

    # Step 2: Create standardised 'symbol' column if not already present.
    if (!"symbol" %in% colnames(annot)) {
      if (!is.null(symbol_column)) {
        if (symbol_column %in% colnames(annot)) {
          # User manually selected a column.
          annot <- annot |>
            mutate(symbol = .data[[symbol_column]])
        } else {
          # User-specified column not found - signal failure.
          return(NULL)
        }
      } else if (!is.null(organism) && organism == "Human") {
        # Human: use hgnc_symbol.
        annot <- annot |>
          mutate(symbol = hgnc_symbol)
      } else if (!is.null(organism) && organism == "Bacteria") {
        # Bacterial: use gene, fall back to Geneid when gene is NA/empty.
        annot <- annot |>
          mutate(symbol = ifelse(is.na(gene) | gene == "", Geneid, gene))
      } else {
        # Unknown organism with no symbol source - need user input.
        return(NULL)
      }
    }

    # Store full annotation
    state$annotation_df(annot)

    # Step 3: Join using the detected ID column.
    if (!is.null(id_column) && id_column == "ensembl_gene_id") {
      # Human: join ensembl_gene_id from annotation to Geneid from de_df.
      dplyr::left_join(de_df, annot, by = c("Geneid" = "ensembl_gene_id"))
    } else {
      # Bacterial / unknown: join on Geneid.
      dplyr::left_join(de_df, annot, by = "Geneid")
    }
  }

  # Build metadata panel from SE.
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

    # Store SE (compatible with DESeq2/SummarizedExperiment accessor functions).
    state$se_obj(se)

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
    missing_cols  <- NULL

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
          merged <- merge_annotation(de_df, annot, organism = organism)
          if (!is.null(merged)) {
            de_df         <- merged
            annot_status  <- "loaded"
            annot_message <- paste("Annotation loaded:", basename(annot_path))
          } else {
            # Required columns missing - store pending data for user selection.
            state$pending_annotation(annot)
            state$pending_de_df(de_df)
            state$se_organism(organism)
            annot_status <- "missing_cols"
            missing_cols <- setdiff(
              if (!is.null(organism) && organism == "Human")
                c("ensembl_gene_id", "hgnc_symbol")
              else if (!is.null(organism) && organism == "Bacteria")
                c("Geneid", "gene")
              else
                c("Geneid", "ensembl_gene_id"),
              colnames(annot)
            )
          }
        }
      } else {
        annot_status  <- "not_found"
        annot_message <- paste("No built-in annotation found for:", annotation)
      }
    }

    # Set de_df after all annotation logic is complete.
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
      } else if (annot_status == "missing_cols") {
        annot_cols <- colnames(isolate(state$pending_annotation()))
        tagList(
          info_panel,
          tags$p(
            icon("exclamation-triangle"),
            HTML(paste0(
              "Annotation is missing required columns: <strong>",
              paste(missing_cols, collapse = ", "),
              "</strong>. Please select a column to use as the gene symbol."
            )),
            style = "color: orange; padding-left: 15px;"
          ),
          tags$div(
            style = "padding-left: 15px;",
            selectInput("symbolColumnSelect", "Use column as symbol:",
                        choices = annot_cols),
            actionButton("applySymbolColumn", "Apply", class = "btn-primary btn-sm")
          )
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

    # Get organism from stored SE.
    se_current <- state$se_obj()
    organism   <- tryCatch(metadata(se_current)$organism, error = function(e) NULL)
    organism   <- if (!is.null(organism) && nzchar(trimws(organism))) organism else NULL

    merged <- merge_annotation(state$de_df(), annot, organism = organism)

    if (!is.null(merged)) {
      state$de_df(merged)

      output$annotationStatus <- renderUI({
        tagList(
          tags$div(style = "padding-left: 15px; margin-bottom: 8px;",
                   tags$p(style = "color: steelblue; margin: 0;",
                          icon("info-circle"), strong("SE loaded")),
                   se_info_ui(se_current, organism = organism)),
          tags$p(icon("check-circle"), "Custom annotation loaded successfully.",
                 style = "color: green; padding-left: 15px;")
        )
      })
      showNotification("Custom annotation loaded successfully.", type = "message")
    } else {
      # Required columns missing - store pending data for interactive selection.
      state$pending_annotation(annot)
      state$pending_de_df(state$de_df())
      state$se_organism(organism)

      annot_cols <- colnames(annot)
      output$annotationStatus <- renderUI({
        tagList(
          tags$div(style = "padding-left: 15px; margin-bottom: 8px;",
                   tags$p(style = "color: steelblue; margin: 0;",
                          icon("info-circle"), strong("SE loaded")),
                   se_info_ui(se_current, organism = organism)),
          tags$p(
            icon("exclamation-triangle"),
            HTML("Annotation is missing required columns. Please select a column to use as the gene symbol."),
            style = "color: orange; padding-left: 15px;"
          ),
          tags$div(
            style = "padding-left: 15px;",
            selectInput("symbolColumnSelect", "Use column as symbol:",
                        choices = annot_cols),
            actionButton("applySymbolColumn", "Apply", class = "btn-primary btn-sm")
          )
        )
      })
    }
  })

  ## ---- Custom symbol-column selection -----
  observeEvent(input$applySymbolColumn, {
    req(input$symbolColumnSelect, state$pending_annotation(), state$pending_de_df())

    annot    <- state$pending_annotation()
    de_df    <- state$pending_de_df()
    organism <- state$se_organism()

    de_df_merged <- merge_annotation(
      de_df,
      annot,
      organism      = organism,
      symbol_column = input$symbolColumnSelect
    )

    if (!is.null(de_df_merged)) {
      state$de_df(de_df_merged)
      # Clear pending state.
      state$pending_annotation(NULL)
      state$pending_de_df(NULL)

      se_current <- state$se_obj()
      output$annotationStatus <- renderUI({
        tagList(
          tags$div(
            style = "padding-left: 15px; margin-bottom: 8px;",
            tags$p(style = "color: steelblue; margin: 0;",
                   icon("info-circle"), strong("SE loaded")),
            se_info_ui(se_current, organism = organism)
          ),
          tags$p(icon("check-circle"),
                 paste0("Annotation loaded. Using '", input$symbolColumnSelect, "' as symbol column."),
                 style = "color: green; padding-left: 15px; margin: 2px 0;")
        )
      })

      showNotification(
        paste("Annotation applied successfully using", input$symbolColumnSelect),
        type = "message"
      )
    }
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
