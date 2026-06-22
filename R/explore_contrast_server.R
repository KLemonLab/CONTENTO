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
    
    selectInput(
      "contrast",
      NULL,
      choices = contrast_choices,
      width = "90%"
    )
  })
  
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ##### UI: Tab Layout #####
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  output$contrastSubTabs <- renderUI({
    tabs <- list(
      tabPanel("Differentially Expressed Genes",
               fluidRow(
                 column(
                   width = 7,
                   div(
                     style = "padding: 15px; background-color: #f8f9fa; border-radius: 10px; border: 1px solid #e3e6ea;",
                     h4(strong("Explore Differentially Expressed Genes (DEGs)"),
                        style = "margin-top: 0; margin-bottom: 8px;"),
                     tags$ul(
                       style = "margin: 5px 0 0 15px; padding:0;",
                       tags$li("Adjust log2FC column, fold-change and p-value cutoffs in the sidebar; download the full gene list or just the DEGs"),
                       tags$li("The volcano plot visualizes all genes, highlighting significant DEGs based on your threshold"),
                       tags$li("Hover over points to view gene details (based on the selected symbol column), and download the plot in the interactive viewer")
                     )
                   )
                 ),
                 column(
                   width = 5,
                   div(
                     style = "padding: 15px; background-color: #f8f9fa; border-radius: 10px; border: 1px solid #e3e6ea;",
                     strong("Download Options"),
                     tags$div(
                       style = "margin-top: 15px; display: flex; gap: 10px; justify-content: center; flex-wrap: wrap;",
                       downloadButton("downloadAllContrastData", "Download All Genes", class = "btn btn-info"),
                       downloadButton("downloadFilteredContrastData", "Download DEGs",  class = "btn btn-info")
                     )
                   )
                 )
               ),
               hr(),
               withSpinner(DTOutput("DETable"), type = 5),
               hr(),
               h4("Volcano Plot of All Genes in the Selected Contrast"),
               withSpinner(plotlyOutput("volcanoPlot", height = "500px",  width = "100%"), type = 5),
      ),
      tabPanel("Functional Enrichment ",
               fluidRow(
                 column(
                   width = 7,
                   div(
                     style = "padding: 15px; background-color: #f8f9fa; border-radius: 10px; border: 1px solid #e3e6ea;",
                     h4(strong("Perform Gene Set Enrichment Analysis (GSEA)"),
                        style = "margin-top: 0; margin-bottom: 8px;"),
                     tags$ul(
                       style = "margin: 5px 0 0 15px; padding:0;",
                       tags$li(tags$b("Human:"), " select curated biological gene sets from MSigDB collections"),
                       tags$li(tags$b("Bacteria:"), " use gene sets defined in annotation columns"),
                       tags$li("Gene sets are ranked by NES (Normalized Enrichment Score)"),
                       tags$li("Leading-edge genes indicate core contributors to enrichment")
                       
                     )
                   )
                 ),
                 column(
                   width = 5,
                   div(
                     style = "padding: 15px; background-color: #f8f9fa; border-radius: 10px; border: 1px solid #e3e6ea;",
                     uiOutput("gseaControlsUI"),
                     tags$div(
                       style = "margin-top: 15px; display: flex; gap: 10px; justify-content: center; flex-wrap: wrap;",
                       downloadButton("downloadGseaResults", "Download GSEA Table", class = "btn btn-info"),
                       downloadButton("downloadGseaPlot",    "Download GSEA Plot",  class = "btn btn-info")
                     )
                   )
                 )
               ),
               hr(),
               withSpinner(DTOutput("gseaResultsTable"), type = 5),
               hr(),
               h4("Top 20 Enriched Gene Sets (FDR < 0.05)"),
               withSpinner(plotOutput("gseaTablePlot", height = "600px"), type = 5)
      )
    )
    do.call(tabsetPanel, c(list(type = "pills"), tabs))
  }) 
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ##### UI: Dynamic GSEA Controls #####
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  output$gseaControlsUI <- renderUI({
    if (!is.null(organism()) && organism() == "Human") {
      # MSigDB for Human
      return(tagList(
        selectInput("gs_collection", "Select Gene Sets from MSigDB",
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
      has_symbol <- "symbol" %in% colnames(state$de_df())
      
      state$de_df() |>
        filter(contrast == input$contrast) |>
        add_de_flags(input$global_log2FC_cutoff, input$global_padj_cutoff, input$global_lfc_col) |>
        mutate(
          tooltip = paste0(
            if (has_symbol) dplyr::coalesce(.data$symbol, Geneid) else Geneid,
            " (", Geneid, ")",
            "\n", if (input$global_lfc_col == "log2FC_shrunk") "log2FC (shrunken)" else "log2FC", 
            ": ", round(.data[[input$global_lfc_col]], 2),
            "\nFDR: ", signif(padj, 3)
          )
        )
    }, error = handle_filter_data_error)
  })
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ##### Reactive: Load Gene Sets (MSigDB or Bacterial) #####
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  genesets <- reactive({
    if (organism() == "Bacteria") {
      req(state$annotation_df(), input$bacterial_geneset_source)
      tryCatch({
        list(
          data  = load_genesets(organism(), annotation_df = state$annotation_df(),
                                bacterial_source = input$bacterial_geneset_source),
          label = input$bacterial_geneset_source
        )
      }, error = handle_geneset_error)
      
    } else if (organism() == "Human") {
      req(input$gs_collection)
      tryCatch({
        list(
          data  = load_genesets(organism(), gs_collection = input$gs_collection,
                                gs_subcollection = input$gs_subcollection),
          label = paste0(input$gs_collection,
                         if (nzchar(input$gs_subcollection %||% "")) paste0("_", input$gs_subcollection))
        )
      }, error = handle_msigdb_error)
      
    } else {
      NULL
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
      pathways_list <- genesets()$data
      
      if (length(pathways_list) == 0) {
        showNotification("No gene sets found in selected column", type = "warning")
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
        arrange(desc(abs(.data[[input$global_lfc_col]]))) |>
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
            backgroundColor = DT::styleEqual(c("up", "down"), c("#d9f2f9", "#f8d7da"))
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
      x_col <- if (input$global_lfc_col %in% colnames(selected_data())) input$global_lfc_col else "log2FC"
      
      gg <- ggplot(selected_data(), aes(x = .data[[x_col]], y = -log10(padj), text = tooltip)) +
        
        geom_point(aes(color = regulated), alpha = 0.7) +
        
        scale_color_manual(
          values = c(
            "up"   = "#00a9cf",  
            "down" = "#d9534f"  
          ),
          na.value = "#b0b0b0", 
          name = NULL
        ) +
        
        geom_vline(
          xintercept = c(-input$global_log2FC_cutoff, input$global_log2FC_cutoff), linetype = "dashed", color = "#888888") +
        geom_hline(
          yintercept = -log10(input$global_padj_cutoff), linetype = "dashed", color = "#888888") +
        
        labs(
          x = case_when(
            x_col == "log2FC_shrunk" ~ "log2 Fold Change (shrunken)",
            x_col == "log2FC"        ~ "log2 Fold Change",
            TRUE                     ~ x_col  
          ),
          y = "-log10(FDR)"
        ) +
        
        theme_minimal()
      
      ggplotly(gg, tooltip = "text") |>
        config(toImageButtonOptions = list(
          format   = "png",
          filename = build_download_filename(input, state,
                                             type = "volcano", contrast = input$contrast,
                                             filters  = list(FC = input$global_log2FC_cutoff, FDR = input$global_padj_cutoff)),
          width    = 1800,
          height   = 900,
          scale    = 2
        ))
    }, error = handle_volcano_plot_error)
  })
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ##### Output: GSEA Table Plot #####
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  output$gseaTablePlot <- renderPlot({
    req(gsea_result())
    
    if (is.null(gsea_result()$tableplot)) {
      plot.new()
      text(0.5, 0.5, "No significant gene sets found (FDR < 0.05)", cex = 1.5)
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
    function() build_download_filename(input, state, 
                                       type = "DE", contrast = input$contrast),
    function(file) {
      req(selected_data())
      selected_data() |>
        arrange(desc(abs(.data[[input$global_lfc_col]]))) |>
        select(-any_of(c("tooltip", "DE"))) |>
        write.csv(file, row.names = FALSE)
    }
  )
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ##### Download Handler: Download DE Genes based on Filters #####
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  output$downloadFilteredContrastData <- downloadHandler(
    function() build_download_filename(input, state, 
                                       type = "DEfiltered", contrast = input$contrast,
                                       filters = list(FC = input$global_log2FC_cutoff, FDR = input$global_padj_cutoff)),
    function(file) {
      req(selected_data(), nrow(selected_data()) > 0)
      selected_data() |>
        filter(DE) |>
        arrange(desc(abs(.data[[input$global_lfc_col]]))) |>
        select(-any_of(c("tooltip", "DE"))) |>
        write.csv(file, row.names = FALSE)
    }
  )
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ##### Download Handler: Download GSEA Table Results #####
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  output$downloadGseaResults <- downloadHandler(
    function() build_download_filename(input, state,
                                       type = "GSEA", contrast = input$contrast,
                                       suffix = genesets()$label),
    function(file) {
      req(gsea_result())
      gsea_result()$fgseaRes |>
        mutate(leadingEdge = sapply(leadingEdge, paste, collapse = "; ")) |>
        write.csv(file, row.names = FALSE)
    }
  )
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ##### Download Handler: Download GSEA Plot Results #####
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  output$downloadGseaPlot <- downloadHandler(
    function() build_download_filename(input, state, ext = "png",
                                       type = "GSEA", contrast = input$contrast,
                                       suffix = genesets()$label),
    function(file) {
      req(gsea_result())
      png(file, width = 1800, height = 900, res = 150)
      print(gsea_result()$tableplot)
      dev.off()
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