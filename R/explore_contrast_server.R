explore_contrast_server <- function(input, output, session, state, organism) {
  
  # == == == == == == == == == == == == == == == == == == == == == == == == ==
  #### SECTION 1: UI RENDERING ####
  # Dynamic UI controls that respond to user inputs and state changes
  # == == == == == == == == == == == == == == == == == == == == == == == == ==
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ##### UI: Contrast Dropdown #####
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  output$contrastSelect <- renderUI({
    req(state$de_df())
    contrast_choices <- unique(state$de_df()$contrast)
    selectInput("contrast", "Select Contrast from Summarized Experiment:", choices = contrast_choices)
  })
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ##### UI: Conditional Tab Layout #####
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  output$contrastSubTabs <- renderUI({
    tabs <- list(
      tabPanel("Selected Genes", 
               withSpinner(DTOutput("DETable"), type = 5),
               tags$div(
                 style = "margin-top: 15px; display: flex; gap: 10px;",
                 downloadButton("downloadAllContrastData", "Download All Genes for Selected Contrast", class = "btn btn-warning"),
                 downloadButton("downloadFilteredContrastData", "Download Filtered Genes in Table", class = "btn btn-warning")
               )
      ),
      tabPanel("Volcano Plot", 
               withSpinner(plotlyOutput("volcanoPlot", height = "500px"), type = 5)
      )
    )
    
    if (!is.null(organism()) && organism() %in% c("Human", "Bacteria")) {
      tabs <- append(tabs, 
                     list(
                       tabPanel("Functional Enrichment ",
                                fluidRow(
                                  # LEFT: explanation
                                  column(
                                    width = 8,
                                    div(
                                      span("Choose an annotation column to perform  GSEA."),
                                      tags$ul(
                                        style = "margin: 5px 0 0 15px; padding:0;",
                                        tags$li("Table includes Gene Sets ranked by NES (Normalized Enrichment Score)"),
                                        tags$li("Plot represents Top 20 significant pathways (FDR < 0.05) with leading edge genes highlighted")
                                      )
                                    )
                                  ),
                                  
                                  # RIGHT: controls
                                  column(
                                    width = 4,
                                    div(
                                      style = "padding-left:10px;",
                                      uiOutput("gseaControlsUI")
                                    )
                                  )
                                ),
                                
                                hr(),
                                withSpinner(DTOutput("gseaResultsTable"), type = 5),
                                hr(),
                                h4("Top 20 Enriched Pathways (FDR < 0.05)"),
                                withSpinner(plotOutput("gseaTablePlot", height = "600px"), type = 5)
                       )
                     )
      )
      
    }
    
    do.call(tabsetPanel, c(list(type = "pills"), tabs))
  }) 
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ##### UI: Dynamic GSEA Controls #####
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  output$gseaControlsUI <- renderUI({
    if (!is.null(organism()) && organism() == "Human") {
      # MSigDB for Human
      return(tagList(
        selectInput("gs_collection", "MSigDB Collection",
                    choices = c("H", "C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8"),
                    selected = "H"),
        conditionalPanel(
          condition = "!(input.gs_collection == 'H' || input.gs_collection == 'C1' || input.gs_collection == 'C6' || input.gs_collection == 'C8')",
          textInput("gs_subcollection", "Subcollection (optional)", placeholder = "e.g., CP:REACTOME")
        ),
        tags$div(style = "margin-top: -10px; margin-bottom: 15px;",
                 tags$a(href = "https://www.gsea-msigdb.org/gsea/msigdb/human/annotate.jsp",
                        target = "_blank",
                        icon("external-link-alt"),
                        "MSigDB"))
      ))
    } else if (!is.null(organism()) && organism() == "Bacteria") {
      # All annotation columns for Bacteria
      avail_cols <- state$available_gsea_columns()
      
      if (length(avail_cols) == 0) {
        return(div(class = "alert alert-warning",
                   icon("exclamation-triangle"),
                   "No annotation columns available. Upload annotation file with functional information."))
      }
      
      return(selectInput("bacterial_geneset_source", "Select Gene Sets from Annotation",
                         choices = avail_cols,
                         selected = avail_cols[[1]]))
    } else {
      return(div(class = "alert alert-info",
                 icon("info-circle"),
                 "Awaiting SE file upload..."))
    }
  })
  
  # == == == == == == == == == == == == == == == == == == == == == == == == ==
  #### SECTION 2: DATA PROCESSING & FILTERING ####
  # Reactive expressions that filter, transform, and prepare data for display
  # == == == == == == == == == == == == == == == == == == == == == == == == ==
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ##### Reactive: All Contrast Data Merged with Annotations #####
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  all_contrast_data_export <- reactive({
    req(state$de_df(), input$contrast, state$annotation_df())
    
    # Get ALL genes for this contrast (no filtering)
    all_data <- state$de_df() %>%
      filter(contrast == input$contrast) %>%
      select(Geneid, log2FC, log2FC_shrunk, stat, padj, regulated) %>%
      arrange(desc(abs(log2FC)))
    
    # Merge with full annotations
    merged <- all_data %>%
      left_join(state$annotation_df(), by = "Geneid")
    
    # Put Geneid first, then symbol if it exists, then everything else
    first_cols <- intersect(c("Geneid", "symbol"), colnames(merged))
    merged %>%
      select(all_of(first_cols), everything())
  })

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ##### Reactive: DE Filtered by User-Defined Cutoffs #####
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  selected_data <- reactive({
    req(state$de_df(), input$contrast)
    
    lfc_cut  <- if (!is.null(input$global_log2FC_cutoff) && !is.na(input$global_log2FC_cutoff)) input$global_log2FC_cutoff else 2
    padj_cut <- if (!is.null(input$global_padj_cutoff)   && !is.na(input$global_padj_cutoff))   input$global_padj_cutoff   else 0.05
    
    tryCatch({
      df <- state$de_df() %>%
        filter(contrast == input$contrast)
      
      # Validate required columns exist
      required_cols <- c("log2FC", "padj")
      missing_cols <- setdiff(required_cols, colnames(df))
      
      if (length(missing_cols) > 0) {
        stop("SE object is missing required DE statistics: ", paste(missing_cols, collapse = ", "))
      }
      
      df %>%
        mutate(
          tooltip = paste0(
            if ("symbol" %in% names(.)) .data$symbol else .data$Geneid, " (", Geneid, ")",
            "\nlog2FC: ", round(log2FC, 2),
            "\nFDR: ", signif(padj, 3)
          ),
          regulated = case_when(
            !is.na(padj) & !is.na(log2FC) & padj < padj_cut & log2FC >  lfc_cut ~ "up",
            !is.na(padj) & !is.na(log2FC) & padj < padj_cut & log2FC < -lfc_cut ~ "down",
            TRUE ~ NA_character_
          ),
          DE = !is.na(regulated)
        )
    }, error = handle_filter_data_error)
  })
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ##### Reactive: DE Filtered Merged with Annotations #####
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  selected_data_export <- reactive({
    req(selected_data(), state$annotation_df())
    
    # Get DE genes with their stats (select only columns that exist)
    stat_cols <- intersect(c("Geneid", "log2FC", "log2FC_shrunk", "stat", "padj", "regulated"),
                           colnames(selected_data()))
    de_data <- selected_data() %>%
      filter(DE) %>%
      select(all_of(stat_cols)) %>%
      arrange(desc(abs(log2FC)))
    
    # Merge with full annotations
    merged <- de_data %>%
      left_join(state$annotation_df(), by = "Geneid")
    
    # Put Geneid first, then symbol if it exists, then everything else
    first_cols <- intersect(c("Geneid", "symbol"), colnames(merged))
    merged %>%
      select(all_of(first_cols), everything())
  })
  
  
  # == == == == == == == == == == == == == == == == == == == == == == == == ==
  #### SECTION 3: GENE SET ENRICHMENT ANALYSIS ####
  # Load gene sets and perform GSEA
  # == == == == == == == == == == == == == == == == == == == == == == == == ==
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ##### Reactive: Load Gene Sets (MSigDB or Bacterial) #####
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  genesets <- reactive({
    # For Bacteria: use user-provided annotation columns
    if (!is.null(organism()) && organism() == "Bacteria") {
      req(state$annotation_df(), input$bacterial_geneset_source)
      tryCatch(
        build_bacterial_genesets(state$annotation_df(), input$bacterial_geneset_source),
        error = handle_geneset_error
      )
      
    } else {
      # For Human: use MSigDB
      req(input$gs_collection)
      
      tryCatch({
        if (nzchar(input$gs_subcollection)) {
          msigdbr(species = "Homo sapiens", collection = input$gs_collection, subcollection = input$gs_subcollection)
        } else {
          msigdbr(species = "Homo sapiens", collection = input$gs_collection)
        }
      }, error = handle_msigdb_error)
    }
  })
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ##### Reactive: GSEA Analysis #####
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  gsea_result <- reactive({
    req(state$de_df(), input$contrast, genesets())
    
    tryCatch({
      # Prepare ranked gene list
      ranks <- state$de_df() %>%
        filter(contrast == input$contrast) %>%
        filter(!is.na(stat))
      
      if (nrow(ranks) == 0) {
        showNotification("No valid statistics found for GSEA", type = "warning")
        return(NULL)
      }
      
      ranks_vec <- setNames(ranks$stat, ranks$Geneid)
      
      # Convert gene sets to named list of gene vectors
      # For bacteria: already in correct format
      # For human: need to convert from msigdbr format
      if (!is.null(organism()) && organism() == "Bacteria") {
        pathways_list <- genesets()
      } else {
        pathways_list <- genesets() %>%
          split(.$gs_name) %>%
          lapply(function(x) x$ensembl_gene)
      }
      
      if (length(pathways_list) == 0) {
        showNotification("No pathways found in selected gene set", type = "warning")
        return(NULL)
      }
      
      # Run GSEA
      fgseaRes <- fgseaMultilevel(
        pathways = pathways_list,
        stats = ranks_vec,
        minSize = 15,
        maxSize = 500
      ) %>%
        arrange(padj, pval)
      
      # Select top pathways by absolute NES
      topPathways <- fgseaRes %>%
        filter(padj < 0.05) %>%
        arrange(desc(abs(NES))) %>%
        slice_head(n = 20) %>%
        pull(pathway)
      
      # Generate table plot (only if there are significant pathways)
      tableplot <- if (length(topPathways) > 0) {
        plotGseaTable(
          pathways = pathways_list[topPathways], 
          stats = ranks_vec, 
          fgseaRes = fgseaRes, 
          gseaParam = 0.5
        )
      } else {
        NULL
      }
      
      list(tableplot = tableplot, fgseaRes = fgseaRes, ranks_vec = ranks_vec, pathways_list = pathways_list)
    },     
    error = handle_gsea_error)
  })
  
  
  # == == == == == == == == == == == == == == == == == == == == == == == == ==
  #### SECTION 4: OUTPUT RENDERING ####
  # Display tables, plots, and interactive visualizations
  # == == == == == == == == == == == == == == == == == == == == == == == == ==
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ##### Output: DE Results Table #####
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  output$DETable <- renderDT(
    {
      req(selected_data())
      
      file_base <- get_download_filename(input, state)
      contrast_str <- if (!is.null(input$contrast) && nzchar(input$contrast)) input$contrast else "contrast"
      contrast_str <- gsub("[^A-Za-z0-9._-]+", "__", contrast_str)
      lfc_cut <- if (!is.null(input$global_log2FC_cutoff) && !is.na(input$global_log2FC_cutoff)) input$global_log2FC_cutoff else 2
      fc_str <- paste0("FC", gsub("\\.", "p", as.character(lfc_cut)))
      file_name <- paste(file_base, contrast_str, fc_str, sep = "__")
      
      display_cols <- intersect(c("Geneid", "symbol", "log2FC", "log2FC_shrunk", "stat", "padj", "regulated"),
                                colnames(selected_data()))
      
      df <- selected_data() %>%
        filter(DE) %>%
        mutate(
          across(any_of(c("log2FC", "log2FC_shrunk", "stat")), ~ round(.x, 2)),
          padj = formatC(padj, format = "e", digits = 2)
        ) %>%
        arrange(desc(abs(log2FC))) %>%
        select(all_of(display_cols))
      
      if (nrow(df) == 0) {
        showNotification("No differentially expressed genes found with current cutoffs", type = "warning")
      }
      
      dt <- datatable(
        df,
        extensions = 'Buttons',
        rownames = FALSE,
        filter = 'top',
        options = list(
          pageLength = 15,
          scrollX = TRUE,
          dom = 'frtip'  
        )
      )
      
      if ("regulated" %in% colnames(df)) {
        dt <- dt %>%
          formatStyle(
            'regulated',
            target = 'row',
            backgroundColor = DT::styleEqual(c("up", "down"), c("lightblue", "pink"))
          )
      }
      dt
    }, server = FALSE)
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ##### Output: Volcano Plot Top DE Genes #####
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  output$volcanoPlot <- renderPlotly({
    req(selected_data())
    
    tryCatch({
      lfc_cut  <- if (!is.null(input$global_log2FC_cutoff) && !is.na(input$global_log2FC_cutoff)) input$global_log2FC_cutoff else 2
      padj_cut <- if (!is.null(input$global_padj_cutoff)   && !is.na(input$global_padj_cutoff))   input$global_padj_cutoff   else 0.05
      
      x_col <- if ("log2FC_shrunk" %in% colnames(selected_data())) "log2FC_shrunk" else "log2FC"
      
      gg <- ggplot(selected_data(), aes(x = .data[[x_col]], y = -log10(padj), text = tooltip)) +
        geom_point(aes(color = DE), alpha = 0.6) +
        scale_color_manual(values = c("TRUE" = "red", "FALSE" = "gray"), guide = "none") +
        geom_vline(xintercept = c(-lfc_cut, lfc_cut), linetype = "dashed") +
        geom_hline(yintercept = -log10(padj_cut), linetype = "dashed") +
        labs(x = paste0("log2 Fold Change", if (x_col == "log2FC_shrunk") " (shrunken)" else ""), y = "-log10(FDR)") +
        theme_minimal()
      
      ggplotly(gg, tooltip = "text")
    }, error = handle_volcano_plot_error)
  })
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ##### Output: GSEA Table Plot #####
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  output$gseaTablePlot <- renderPlot({
    req(gsea_result())
    
    if (is.null(gsea_result()$tableplot)) {
      plot.new()
      text(0.5, 0.5, "No significant pathways found (FDR < 0.05)", cex = 1.5)
    } else {
      gsea_result()$tableplot
    }
  })
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ##### Output: GSEA Results Table #####
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  output$gseaResultsTable <- renderDT(
    {
      req(gsea_result())
      
      # Prepare file name
      file_base <- get_download_filename(input, state)
      contrast_str <- if (!is.null(input$contrast) && nzchar(input$contrast)) input$contrast else "contrast"
      contrast_str <- gsub("[^A-Za-z0-9._-]+", "__", contrast_str)
      
      # Get gene set source name (handles both Human and Bacteria)
      if (!is.null(organism()) && organism() == "Bacteria") {
        gs_str <- if (!is.null(input$bacterial_geneset_source)) {
          gsub("[^A-Za-z0-9._-]+", "__", input$bacterial_geneset_source)
        } else {
          "unknown"
        }
      } else {
        gs_str <- paste0(input$gs_collection, if (nzchar(input$gs_subcollection)) paste0("_", input$gs_subcollection) else "")
      }
      file_name <- paste(file_base, contrast_str, "GSEA", gs_str, sep = "__")
      
      # Format the dataframe
      df <- gsea_result()$fgseaRes %>%
        mutate(
          across(c(pval, padj), ~ formatC(.x, format = "e", digits = 2)),
          across(c(ES, NES), ~ round(.x, 3)),
          leadingEdge = sapply(leadingEdge, function(x) paste(head(x, 10), collapse = "; "))
        ) %>%
        select(pathway, pval, padj, ES, NES, size, leadingEdge)
      
      datatable(
        df,
        extensions = 'Buttons',
        rownames = FALSE,
        filter = 'top',
        options = list(
          pageLength = 15,
          scrollX = TRUE,
          dom = 'rtip',
          columnDefs = list(
            list(
              targets = 6,  
              width = '300px',
              render = JS(
                "function(data, type, row, meta) {",
                "  return '<div style=\"max-width:300px; overflow-x:auto; white-space:nowrap;\">' + data + '</div>';",
                "}"
              )
            )
          )
        )
      )
    }, server = FALSE)

  
  # == == == == == == == == == == == == == == == == == == == == == == == == ==
  #### SECTION 5: DOWNLOAD HANDLERS ####
  # Manage file downloads with proper naming and filtering
  # == == == == == == == == == == == == == == == == == == == == == == == == ==
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ##### Download Handler: Download All Genes for Selected Contrast #####
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  output$downloadAllContrastData <- downloadHandler(
    filename = function() {
      file_base <- get_download_filename(input, state)
      contrast_str <- if (!is.null(input$contrast) && nzchar(input$contrast)) input$contrast else "contrast"
      contrast_str <- gsub("[^A-Za-z0-9._-]+", "__", contrast_str)
      
      paste(file_base, contrast_str, "all_genes.csv", sep = "__")
    },
    content = function(file) {
      # Get all genes for this contrast with annotations
      full_data <- all_contrast_data_export() %>%
        mutate(
          across(any_of(c("log2FC", "log2FC_shrunk")), ~ round(.x, 2)),
          padj = formatC(padj, format = "e", digits = 2)
        )
      
      write.csv(full_data, file, row.names = FALSE)
    }
  )
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ##### Download Handler: Download Filtered Genes in Table #####
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  output$downloadFilteredContrastData <- downloadHandler(
    filename = function() {
      file_base <- get_download_filename(input, state)
      contrast_str <- if (!is.null(input$contrast) && nzchar(input$contrast)) input$contrast else "contrast"
      contrast_str <- gsub("[^A-Za-z0-9._-]+", "__", contrast_str)
      lfc_cut <- if (!is.null(input$global_log2FC_cutoff) && !is.na(input$global_log2FC_cutoff)) input$global_log2FC_cutoff else 2
      fc_str <- paste0("FC", gsub("\\.", "p", as.character(lfc_cut)))
      
      paste(file_base, contrast_str, fc_str, "full.csv", sep = "__")
    },
    content = function(file) {
      req(input$DETable_rows_all)
      
      # Get full export data
      full_data <- selected_data_export() %>%
        mutate(
          across(any_of(c("log2FC", "log2FC_shrunk")), ~ round(.x, 2)),
          padj = formatC(padj, format = "e", digits = 2)
        )
      
      # Get filtered row indices from DT
      filtered_indices <- input$DETable_rows_all
      
      # Get the display data to extract Geneids
      display_data <- selected_data() %>%
        filter(DE) %>%
        arrange(desc(abs(log2FC)))
      
      # Extract Geneids from filtered rows
      filtered_geneids <- display_data[filtered_indices, "Geneid", drop = TRUE]
      
      # Filter full data to match user's selection
      filtered_full <- full_data %>%
        filter(Geneid %in% filtered_geneids)
      
      write.csv(filtered_full, file, row.names = FALSE)
    }
  )
  
  
  # == == == == == == == == == == == == == == == == == == == == == == == == ==
  #### SECTION 6: ERROR HANDLERS ####
  # Helper functions for error management
  # == == == == == == == == == == == == == == == == == == == == == == == == ==
  
  handle_filter_data_error <- \(e) {
    showNotification(paste("Error filtering contrast data:", e$message), type = "error")
    NULL
  }
  
  handle_volcano_plot_error <- \(e) {
    showNotification(paste("Error creating volcano plot:", e$message), type = "error")
    NULL
  }
  
  handle_geneset_error <- \(e) {
    showNotification(paste("Error loading gene sets:", e$message), type = "error")
    NULL
  }
  
  handle_msigdb_error <- \(e) {
    showNotification(paste("Error loading MSigDB:", e$message), type = "error")
    NULL
  }
  
  handle_gsea_error <- \(e) {
    showNotification(paste("Error running GSEA:", e$message), type = "error")
    NULL
  }
  
  handle_enrichment_plot_error <- \(e) {
    plot.new()
    text(0.5, 0.5, paste("Error creating plot:", e$message), cex = 1)
  }
}