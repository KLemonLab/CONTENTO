load_files_server <- function(input, output, session, state) {
  
  # == == == == == == == == == == == == == == == == == == == == == == == == ==
  #### SECTION 1: UI RENDERING ####
  # == == == == == == == == == == == == == == == == == == == == == == == == ==
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ##### Output: Organism (conditional UI bridge) #####
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  output$se_organism <- reactive({ state$se_organism() })
  outputOptions(output, "se_organism", suspendWhenHidden = FALSE)
  
  # == == == == == == == == == == == == == == == == == == == == == == == == ==
  #### SECTION 2: ERROR HANDLERS ####
  # == == == == == == == == == == == == == == == == == == == == == == == == ==
  
  handle_se_read_error <- \(e) { showModal(modalDialog(title = "SE File error", paste("Error reading SE file:", e$message), 
                                                       easyClose = TRUE, footer = NULL)); NULL }
  handle_contrast_extraction_error <- \(e) { showModal(modalDialog(title = "Contrast extraction error", paste("Error extracting contrasts:", e$message), 
                                                                   easyClose = TRUE, footer = NULL)); NULL }
  handle_annotation_read_error <- \(e) { showModal(modalDialog(title = "Annotation error", paste("Error reading annotation file:", e$message), 
                                                               easyClose = TRUE, footer = NULL)); NULL }
  handle_annotation_load_error <- \(e) { showNotification(paste("Annotation load error:", e$message), type = "error"); NULL }

  
  # == == == == == == == == == == == == == == == == == == == == == == == == ==
  #### SECTION 3: INTERNAL HELPERS ####
  # == == == == == == == == == == == == == == == == == == == == == == == == ==

  # After annotation is loaded / updated, refresh all derived state:
  # available_gsea_columns, ensembl_col, db_species (which depends on organism).
  refresh_annotation_state <- function(annot, organism) {
    if (!is.null(annot)) {
      state$available_gsea_columns(get_gsea_columns(annot))
    } else {
      state$available_gsea_columns(list())
    }
    
    # Determine Ensembl column (if any)
    state$ensembl_col(
      if (se_rownames_are_ensembl(state$se_obj())) {
        "rownames"                          # SE rownames already Ensembl
      } else if (!is.null(annot)) {
        detect_ensembl_col(annot)           # Look for an Ensembl column in annotation 
      } else {
        NULL                                # No annotation & no Ensembl rownames -> MSigDB unavailable
      }
    )
    
    # Determine msigdbr db_species
    db_override <- tryCatch(metadata(state$se_obj())$db_species, error = function(e) NULL)
    db_sp <- get_msigdbr_db_species(organism, db_override)
    state$db_species(db_sp)
  }
  
  
  # == == == == == == == == == == == == == == == == == == == == == == == == ==
  #### SECTION 4: FILE LOADING & DATA PROCESSING ####
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
        easyClose = TRUE, footer = NULL
      ))
      return(NULL)
    }
    
    # Reset state for new SE upload
    state$annotation_df(NULL)
    state$varpart_obj(NULL)
    
    state$se_obj(se)
    
    de_df <- tryCatch(extract_de_results(se), error = handle_contrast_extraction_error)
    req(de_df)
    
    # Extract variance partition if varpart_* columns are present
    rd           <- as.data.frame(rowData(se))
    varpart_cols <- grep("^varpart_", colnames(rd), value = TRUE)
    if (length(varpart_cols) > 0) {
      vp_mat <- as.data.frame(rd[, varpart_cols, drop = FALSE])
      colnames(vp_mat) <- sub("^varpart_", "", colnames(vp_mat))
      rownames(vp_mat) <- rownames(se)
      state$varpart_obj(list(varPart = vp_mat))
    }
    
    # Organism from SE metadata
    se_organism <- tryCatch(metadata(se)$organism, error = function(e) NULL)
    organism    <- if (!is.null(se_organism) && nzchar(trimws(se_organism))) se_organism else NULL
    state$se_organism(organism)
    
    # Annotation name from SE metadata
    se_annotation <- tryCatch(metadata(se)$annotation, error = function(e) NULL)
    annotation    <- if (!is.null(se_annotation) && nzchar(trimws(se_annotation))) se_annotation else NULL
    
    annot_status  <- "none"
    annot_message <- "No annotation specified in SE metadata."
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
            merged      <- merged_result$result
            annot_cols  <- colnames(annot)
            sym_primary <- default_symbol_col(annot_cols)
            
            de_df <- apply_symbol(merged, sym_primary, "Geneid")
            annot <- apply_symbol(annot,  sym_primary, "Geneid")
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
    refresh_annotation_state(state$annotation_df(), state$se_organism())
    
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
          tags$p(icon("exclamation-triangle"), annot_message,
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
      merged      <- merged_result$result
      annot_cols  <- colnames(annot)
      sym_primary <- default_symbol_col(annot_cols)
      
      merged <- apply_symbol(merged, sym_primary, "Geneid")
      annot  <- apply_symbol(annot,  sym_primary, "Geneid")
      state$annotation_df(annot)
      state$de_df(merged)
      
      refresh_annotation_state(state$annotation_df(), state$se_organism())
      
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
      output$annotationStatus <- renderUI({
        tagList(
          tags$div(style = "padding-left: 15px; margin-bottom: 8px;",
                   tags$p(style = "color: steelblue; margin: 0;",
                          icon("info-circle"), strong("SE loaded")),
                   se_info_ui(se_current, organism = organism)),
          tags$p(icon("exclamation-triangle"), merged_result$message,
                 style = "color: red; padding-left: 15px;"),
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
    
    annot <- state$annotation_df()
    if (!is.null(annot) && primary_col %in% colnames(annot)) {
      state$annotation_df(apply_symbol(annot, primary_col, secondary_col))
    }
    
    showNotification(
      paste0("Symbol column: '", primary_col, "'"),
      type = "message", duration = 2
    )
  })
  
}
