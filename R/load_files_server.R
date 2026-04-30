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
  # Returns a list with:
  #   $result: merged data frame or NULL
  #   $message: diagnostic message (only present if result is NULL)
  merge_annotation <- function(de_df, annot) {
    # Validate inputs
    if (!is.data.frame(de_df) || !is.data.frame(annot)) {
      return(list(result = NULL, message = "Invalid data frame structure"))
    }
    
    if (nrow(annot) == 0) {
      return(list(result = NULL, message = "Annotation file is empty"))
    }
    
    # Check for exact "Geneid" match first
    if ("Geneid" %in% colnames(de_df) && "Geneid" %in% colnames(annot)) {
      # Check type compatibility
      if (!is.character(annot[["Geneid"]])) {
        return(list(
          result = NULL,
          message = paste0("Column 'Geneid' is type ", class(annot[["Geneid"]])[1], 
                          " (expected character). Cannot join.")
        ))
      }
      
      merged <- tryCatch(
        dplyr::left_join(de_df, annot, by = "Geneid"),
        error = function(e) NULL
      )
      if (!is.null(merged)) {
        return(list(result = merged))
      } else {
        return(list(result = NULL, message = "Join by 'Geneid' failed"))
      }
    }
    
    # Search for gene ID column using strict pattern
    # Matches: geneid, gene_id, ensembl_gene_id, etc.
    # Does NOT match: gene_callers_id, other_id, etc.
    gene_id_pattern <- "^(.*_)?gene[_\\s]?id$"
    annot_gene_col <- grep(gene_id_pattern, colnames(annot), 
                           ignore.case = TRUE, value = TRUE)
    
    if (length(annot_gene_col) > 0) {
      annot_gene_col <- annot_gene_col[1]
      
      # Check type compatibility
      if (!is.character(annot[[annot_gene_col]])) {
        return(list(
          result = NULL, 
          message = paste0("Column '", annot_gene_col, "' is type ", 
                           class(annot[[annot_gene_col]])[1], 
                           " (expected character). Cannot join with Geneid.")
        ))
      }
      
      # Attempt join
      merged <- tryCatch(
        dplyr::left_join(de_df, annot, by = c("Geneid" = annot_gene_col)),
        error = function(e) NULL
      )
      
      if (!is.null(merged)) {
        return(list(result = merged))
      } else {
        return(list(result = NULL, message = paste0("Join by '", annot_gene_col, "' failed")))
      }
    }
    
    # No usable join column found
    return(list(result = NULL, message = "No compatible gene ID column found (expected 'Geneid' or similar pattern)"))
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
  
  # Determine default primary symbol column from available annotation columns.
  default_symbol_col <- function(annot_cols) {
    found <- intersect(c("Gene", "gene", "hgnc_symbol", "gene_name", "symbol"), annot_cols)
    if (length(found) > 0) found[1] else annot_cols[1]
  }
  
  # Build the single-dropdown symbol-column selection UI from annotation columns.
  build_symbol_select_ui <- function(annot_cols, primary_sel) {
    tagList(
      selectInput("symbolPrimaryCol",
                  "Select symbol column (Geneid used as fallback):",
                  choices  = annot_cols,
                  selected = primary_sel)
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
          merged_result <- merge_annotation(de_df, annot)
          
          if (!is.null(merged_result$result)) {
            merged <- merged_result$result
            # Determine and auto-apply default symbol column (Geneid always fallback).
            annot_cols  <- colnames(annot)
            sym_primary <- default_symbol_col(annot_cols)

            de_df       <- apply_symbol(merged, sym_primary, "Geneid")
            annot       <- apply_symbol(annot,  sym_primary, "Geneid")
            state$annotation_df(annot)

            annot_status  <- "loaded"
            annot_message <- "Annotation loaded"
          } else {
            annot_status  <- "no_join_col"
            annot_message <- merged_result$message
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
                   build_symbol_select_ui(annot_cols, sym_primary))
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

    merged_result <- merge_annotation(state$de_df(), annot)

    if (!is.null(merged_result$result)) {
      merged <- merged_result$result
      annot_cols  <- colnames(annot)
      sym_primary <- default_symbol_col(annot_cols)

      merged <- apply_symbol(merged, sym_primary, "Geneid")
      annot  <- apply_symbol(annot,  sym_primary, "Geneid")
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
                   build_symbol_select_ui(annot_cols, sym_primary))
        )
      })
      showNotification("Custom annotation loaded", type = "message", duration = 3)
    } else {
      # Join failed - show specific error
      output$annotationStatus <- renderUI({
        tagList(
          tags$div(style = "padding-left: 15px; margin-bottom: 8px;",
                   tags$p(style = "color: steelblue; margin: 0;",
                          icon("info-circle"), strong("SE loaded")),
                   se_info_ui(se_current, organism = organism)),
          tags$p(
            icon("exclamation-triangle"),
            merged_result$message,
            style = "color: red; padding-left: 15px;"
          ),
          fileInput("annotFile", "Upload Annotation (.rds)", accept = ".rds")
        )
      })
    }
  })

  ## ---- Symbol-column selection -----
  observeEvent(input$symbolPrimaryCol, {
    req(input$symbolPrimaryCol, state$de_df())
    
    primary_col   <- input$symbolPrimaryCol
    secondary_col <- "Geneid"
    
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
      paste0("Symbol column: '", primary_col, "'"),
      type = "message",
      duration = 2
    )
  })
  
  #== == == == == == == == == == == == == == == == ==
  #===== OUTPUTS FOR UI CONDITIONALS =================
  #== == == == == == == == == == == == == == == == ==
  
  # Output organism for conditional UI in app.R
  output$se_organism <- reactive({
    state$se_organism()
  })
  outputOptions(output, "se_organism", suspendWhenHidden = FALSE)
  

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