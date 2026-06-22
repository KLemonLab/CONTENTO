load_files_server <- function(input, output, session, state) {
  
  # == == == == == == == == == == == == == == == == == == == == == == == == ==
  #### SECTION 1: UI RENDERING ####
  # == == == == == == == == == == == == == == == == == == == == == == == == ==
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ##### UI: Symbol Column Selection #####
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ##### UI: SE Information Display #####
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ##### Output: Annotation Status Panel #####
  # output$annotationStatus — It is rendered inside observeEvents
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ##### Output: Organism (conditional UI bridge) #####
  # Exposes se_organism to conditionalPanel() in the UI layer
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  output$se_organism <- reactive({ state$se_organism() })
  outputOptions(output, "se_organism", suspendWhenHidden = FALSE)
  
  # == == == == == == == == == == == == == == == == == == == == == == == == ==
  #### SECTION 2: ERROR HANDLERS ####
  # == == == == == == == == == == == == == == == == == == == == == == == == ==
  
  handle_se_read_error <- \(e) {
    showModal(modalDialog(
      title = "SE File error",
      paste("Error reading SE file:", e$message),
      easyClose = TRUE,
      footer = NULL
    ))
    NULL
  }
  
  handle_contrast_extraction_error <- \(e) {
    showModal(modalDialog(
      title = "Contrast extraction error",
      paste("Error extracting contrasts:", e$message),
      easyClose = TRUE,
      footer = NULL
    ))
    NULL
  }
  
  handle_annotation_read_error <- \(e) {
    showModal(modalDialog(
      title = "Annotation error",
      paste("Error reading annotation file:", e$message),
      easyClose = TRUE,
      footer = NULL
    ))
    NULL
  }
  
  handle_annotation_load_error <- \(e) {
    showNotification(paste("Annotation load error:", e$message), type = "error")
    NULL
  }
  
  # == == == == == == == == == == == == == == == == == == == == == == == == ==
  #### SECTION 3: FILE LOADING & DATA PROCESSING ####
  # == == == == == == == == == == == == == == == == == == == == == == == == ==
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ##### Event: SE File Upload #####
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  observeEvent(input$seFile, {
    req(input$seFile)
    
    se <- tryCatch(readRDS(input$seFile$datapath), error = handle_se_read_error)
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
    
    de_df <- tryCatch(extract_de_results(se), error = handle_contrast_extraction_error)
    
    req(de_df)
    
    # Extract variance partition if varpart_* columns are present.
    rd           <- as.data.frame(rowData(se))
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
    annot         <- NULL
    
    if (!is.null(annotation)) {
      annot_path <- find_annotation_path(annotation)
      
      if (!is.null(annot_path)) {
        annot <- tryCatch(readRDS(annot_path), error = handle_annotation_load_error)
        
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
            annot_status  <- "failed"
            annot_message <- merged_result$message
          }
        }
      } else {
        annot_status  <- "not_found"
        annot_message <- paste("No built-in annotation found for:", annotation)
      }
    }
    
    state$de_df(de_df)
    
    # Detect available GSEA columns for gene set selection
    if (!is.null(annot)) {
      available_gsea_cols <- get_gsea_columns(annot)
      state$available_gsea_columns(available_gsea_cols)
    } else {
      state$available_gsea_columns(list())
    }
    
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
                 annot_message,
                 style = "color: orange; padding-left: 15px;"),
          fileInput("annotFile", "Upload Annotation (.rds)", accept = ".rds")
        )
      }
    })
  })
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ##### Event: Custom Annotation Upload #####
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  observeEvent(input$annotFile, {
    req(input$annotFile, state$de_df())
    
    annot <- tryCatch(readRDS(input$annotFile$datapath), error = handle_annotation_read_error)
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
      
      # Detect available GSEA columns for gene set selection
      available_gsea_cols <- get_gsea_columns(annot)
      state$available_gsea_columns(available_gsea_cols)
      
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
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ##### Event: Symbol Column Selection #####
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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
  
}
