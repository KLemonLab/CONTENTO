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
                 downloadButton("downloadFilteredContrastData", "Download Differentially Expresssed Genes", class = "btn btn-warning")
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
  #### SECTION 2: DATA PROCESSING & ANALYSIS ####
  # Reactive expressions that transform and analyze data based on user selections and inputs
  # == == == == == == == == == == == == == == == == == == == == == == == == ==

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ##### Reactive: DEG for the Selected Contrast with DE based on Global Cutoffs #####
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  selected_data <- reactive({
    req(state$de_df(), input$contrast)
    
    tryCatch({
      state$de_df() |>
        filter(contrast == input$contrast) |>
        add_de_flags(input$global_log2FC_cutoff, input$global_padj_cutoff) |>
        mutate(
          tooltip = paste0(
            dplyr::coalesce(.data$symbol, Geneid),
            " (", Geneid, ")",
            "\nlog2FC: ", round(log2FC, 2),
            "\nFDR: ", signif(padj, 3)
          )
        )
    }, error = handle_filter_data_error)
  })
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ##### Reactive: Load Gene Sets (MSigDB or Bacterial) #####
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  genesets <- reactive({
    if (!is.null(organism()) && organism() == "Bacteria") {
      req(state$annotation_df(), input$bacterial_geneset_source)
      tryCatch(
        load_genesets(organism(), annotation_df = state$annotation_df(),
                      bacterial_source = input$bacterial_geneset_source),
        error = handle_geneset_error
      )
    } else {
      req(input$gs_collection)
      tryCatch(
        load_genesets(organism(), gs_collection = input$gs_collection,
                      gs_subcollection = input$gs_subcollection),
        error = handle_msigdb_error
      )
    }
  })
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ##### Reactive: GSEA Analysis #####
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  gsea_result <- reactive({
    req(state$de_df(), input$contrast, genesets())
    
    tryCatch({
      # Prepare ranked gene list for GSEA
      ranks <- state$de_df() |>
        filter(contrast == input$contrast) |>
        filter(!is.na(stat))
      
      if (nrow(ranks) == 0) {
        showNotification("No valid statistics found for GSEA", type = "warning")
        return(NULL)
      }
      ranks_vec <- setNames(ranks$stat, ranks$Geneid)
      
      # Get gene sets for GSEA
      pathways_list <- genesets()
      
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
      ) |>
        arrange(padj, pval)
      
      # Select top pathways by absolute NES
      topPathways <- fgseaRes |>
        filter(padj < 0.05) |>
        arrange(desc(abs(NES))) |>
        slice_head(n = 20) |>
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
  #### SECTION 3: OUTPUT RENDERING ####
  # Display tables, plots, and interactive visualizations
  # == == == == == == == == == == == == == == == == == == == == == == == == ==
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ##### Output: DE Results Table #####
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  output$DETable <- renderDT(
    {
      req(selected_data())
    
      display_cols <- intersect(c("Geneid", "symbol", "log2FC", "log2FC_shrunk", "stat", "padj", "regulated"),
                                colnames(selected_data()))
      
      df <- selected_data() |>
        filter(DE) |>
        mutate(
          across(any_of(c("log2FC", "log2FC_shrunk", "stat")), ~ round(.x, 2)),
          padj = formatC(padj, format = "e", digits = 2)
        ) |>
        arrange(desc(abs(log2FC))) |>
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
          dom = 'rtip'  
        )
      )
      
      if ("regulated" %in% colnames(df)) {
        dt <- dt |>
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
      x_col <- if ("log2FC_shrunk" %in% colnames(selected_data())) "log2FC_shrunk" else "log2FC"
      
      gg <- ggplot(selected_data(), aes(x = .data[[x_col]], y = -log10(padj), text = tooltip)) +
        geom_point(aes(color = DE), alpha = 0.6) +
        scale_color_manual(values = c("TRUE" = "red", "FALSE" = "gray"), guide = "none") +
        geom_vline(xintercept = c(-input$global_log2FC_cutoff, input$global_log2FC_cutoff), linetype = "dashed") +
        geom_hline(yintercept = -log10(input$global_padj_cutoff), linetype = "dashed") +
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
      df <- gsea_result()$fgseaRes |>
        mutate(
          across(c(pval, padj), ~ formatC(.x, format = "e", digits = 2)),
          across(c(ES, NES), ~ round(.x, 3)),
          leadingEdge = sapply(leadingEdge, function(x) paste(head(x, 10), collapse = "; "))
        ) |>
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
  #### SECTION 4: DOWNLOAD HANDLERS ####
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
      req(selected_data())
      
      selected_data() |>
        arrange(desc(abs(log2FC))) |>
        select(-any_of(c("tooltip", "DE"))) |>
        write.csv(file, row.names = FALSE)
    }
  )
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ##### Download Handler: Download DE Genes based on Filters #####
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  output$downloadFilteredContrastData <- downloadHandler(
    filename = function() {
      file_base <- get_download_filename(input, state)
      contrast_str <- if (!is.null(input$contrast) && nzchar(input$contrast)) input$contrast else "contrast"
      contrast_str <- gsub("[^A-Za-z0-9._-]+", "__", contrast_str)
      fc_str  <- paste0("FC",  gsub("\\.", "p", as.character(input$global_log2FC_cutoff)))
      pvl_str <- paste0("FDR", gsub("\\.", "p", as.character(input$global_padj_cutoff)))
      paste(file_base, contrast_str, fc_str, pvl_str, "filtered.csv", sep = "__")
    },
    content = function(file) {
      req(selected_data(), nrow(selected_data()) > 0)
      selected_data() |>
        filter(DE) |>
        arrange(desc(abs(log2FC))) |>
        select(-any_of(c("tooltip", "DE"))) |>
        write.csv(file, row.names = FALSE)
    }
  )
  
  
  # == == == == == == == == == == == == == == == == == == == == == == == == ==
  #### SECTION 5: ERROR HANDLERS ####
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