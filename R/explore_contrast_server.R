explore_contrast_server <- function(input, output, session, state, organism) {
  
  #==============================
  # UI: Contrast dropdown 
  #==============================
  output$contrastSelect <- renderUI({
    req(state$de_df())
    contrast_choices <- unique(state$de_df()$contrast)
    selectInput("contrast", "Select Contrast", choices = contrast_choices)
  })
  
  #==============================
  # UI: Conditional Sub-tabs
  #==============================
  output$contrastSubTabs <- renderUI({
    tabs <- list(
      tabPanel("Selected Genes", withSpinner(DTOutput("DETable"), type = 5)),
      tabPanel("Volcano Plot", withSpinner(plotlyOutput("volcanoPlot", height = "500px"), type = 5))
    )
    
    if (!is.null(organism()) && organism() %in% c("Human", "Bacteria")) {
      tabs <- append(tabs, 
                     list(
                       tabPanel("GSEA",
                                tabsetPanel(
                                  tabPanel("Overview",
                                           h4("Top 20 GSEA Results"),
                                           withSpinner(plotOutput("gseaTablePlot", height = "600px"), type = 5),
                                           hr(),
                                           h4("All GSEA Results"),
                                           withSpinner(DTOutput("gseaResultsTable"), type = 5)
                                  ),
                                  tabPanel("Pathway Detail",
                                           fluidRow(
                                             column(12,
                                                    uiOutput("pathwaySelectUI_gsea"),
                                                    hr(),
                                                    withSpinner(plotOutput("enrichmentPlot", height = "400px"), type = 5)
                                             )
                                           )
                                  )
                                )
                       )
                     )
      )
    }
    
    do.call(tabsetPanel, tabs)
  }) 
  
  #==============================
  # UI: Selectors for GSEA
  #==============================
  output$pathwaySelectUI_gsea <- renderUI({
    req(gsea_result())
    
    pathways <- gsea_result()$fgseaRes %>%
      filter(padj < 0.05) %>%
      arrange(padj) %>%
      pull(pathway)
    
    if (length(pathways) > 0) {
      selectInput("selected_pathway_gsea", "Select Pathway:", 
                  choices = pathways, 
                  selected = pathways[1],
                  width = "100%")
    } else {
      div(
        class = "alert alert-warning",
        icon("exclamation-triangle"),
        "No significant pathways found (FDR < 0.05)"
      )
    }
  })
  
  #==============================
  # Reactive: Filtered DE data
  #==============================
  selected_data <- reactive({
    req(state$de_df(), input$contrast)
    
    lfc_cut  <- if (!is.null(input$global_log2FC_cutoff) && !is.na(input$global_log2FC_cutoff)) input$global_log2FC_cutoff else 2
    padj_cut <- if (!is.null(input$global_padj_cutoff)   && !is.na(input$global_padj_cutoff))   input$global_padj_cutoff   else 0.05
    
    tryCatch({
      df <- state$de_df() %>%
        filter(contrast == input$contrast)
      
      if (!"log2FC" %in% colnames(df)) df$log2FC <- NA_real_
      if (!"padj"   %in% colnames(df)) df$padj   <- NA_real_
      
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
    }, error = function(e) {
      showNotification(paste("Error filtering data:", e$message), type = "error")
      NULL
    })
  })
  
  #==============================
  # Reactive: Export table with full annotations for single contrast
  #==============================
  selected_data_export <- reactive({
    req(selected_data(), state$annotation_df())
    
    # Get DE genes with their stats (select only columns that exist)
    stat_cols <- intersect(c("Geneid", "log2FC", "log2FC_shrunk", "padj", "regulated"),
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
  
  #==============================
  # Reactive: Gene sets from MSigDB or Bacterial Annotations
  #==============================
  genesets <- reactive({
    # For Bacteria: use functional annotations
    if (!is.null(organism()) && organism() == "Bacteria") {
      req(state$annotation_df(), input$bacterial_geneset_source)
      
      tryCatch({
        build_bacterial_genesets(state$annotation_df(), input$bacterial_geneset_source)
      }, error = function(e) {
        showNotification(paste("Error loading bacterial gene sets:", e$message), type = "error")
        NULL
      })
      
    } else {
      # For Human: use MSigDB
      req(input$gs_collection)
      
      tryCatch({
        if (nzchar(input$gs_subcollection)) {
          msigdbr(species = "Homo sapiens", collection = input$gs_collection, subcollection = input$gs_subcollection)
        } else {
          msigdbr(species = "Homo sapiens", collection = input$gs_collection)
        }
      }, error = function(e) {
        showNotification(paste("Error loading gene sets:", e$message), type = "error")
        NULL
      })
    }
  })
  
  #==============================
  # Reactive: GSEA Analysis
  #==============================
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
      
    }, error = function(e) {
      showNotification(paste("Error running GSEA:", e$message), type = "error")
      NULL
    })
  })
  
  #==============================
  # Output: DE Table
  #==============================
  output$DETable <- renderDT(
    {
      req(selected_data())
      
      file_base <- get_download_filename(input, state)
      contrast_str <- if (!is.null(input$contrast) && nzchar(input$contrast)) input$contrast else "contrast"
      contrast_str <- gsub("[^A-Za-z0-9._-]+", "__", contrast_str)
      lfc_cut <- if (!is.null(input$global_log2FC_cutoff) && !is.na(input$global_log2FC_cutoff)) input$global_log2FC_cutoff else 2
      fc_str <- paste0("FC", gsub("\\.", "p", as.character(lfc_cut)))
      file_name <- paste(file_base, contrast_str, fc_str, sep = "__")
      
      display_cols <- intersect(c("Geneid", "symbol", "log2FC", "log2FC_shrunk", "padj", "regulated"),
                                colnames(selected_data()))
      
      df <- selected_data() %>%
        filter(DE) %>%
        mutate(
          across(any_of(c("log2FC", "log2FC_shrunk")), ~ round(.x, 2)),
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
          dom = 'Bfrtip',
          buttons = list(
            list(
              extend = 'csv',
              text = 'Download DE Genes',
              filename = file_name,
              exportOptions = list(modifier = list(page = "all"))
            )
          )
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
    },
    server = FALSE   
  )
  
  #==============================
  # Download Handler: Full annotations with user's filter
  #==============================
  output$downloadDETableFull <- downloadHandler(
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
  
  #==============================
  # Output: Volcano Plot of Top DE Genes
  #==============================
  output$volcanoPlot <- renderPlotly({
    req(selected_data())
    
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
  })
  
  #==============================
  # Output: GSEA Table Plot
  #==============================
  output$gseaTablePlot <- renderPlot({
    req(gsea_result())
    
    if (is.null(gsea_result()$tableplot)) {
      plot.new()
      text(0.5, 0.5, "No significant pathways found (FDR < 0.05)", cex = 1.5)
    } else {
      gsea_result()$tableplot
    }
  })
  
  #==============================
  # Output: GSEA Results Table
  #==============================
  output$gseaResultsTable <- renderDT(
    {
      req(gsea_result())
      
      # Prepare file name
      file_base <- get_download_filename(input, state)
      contrast_str <- if (!is.null(input$contrast) && nzchar(input$contrast)) input$contrast else "contrast"
      contrast_str <- gsub("[^A-Za-z0-9._-]+", "__", contrast_str)
      gs_str <- paste0(input$gs_collection, if (nzchar(input$gs_subcollection)) paste0("_", input$gs_subcollection) else "")
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
          dom = 'Bfrtip',
          columnDefs = list(
            list(
              targets = 6,  # leadingEdge column (0-indexed)
              width = '300px',
              render = JS(
                "function(data, type, row, meta) {",
                "  return '<div style=\"max-width:300px; overflow-x:auto; white-space:nowrap;\">' + data + '</div>';",
                "}"
              )
            )
          ),
          buttons = list(
            list(
              extend = 'csv',
              text = 'Download Full GSEA Results',
              filename = file_name,
              exportOptions = list(modifier = list(page = "all"))
            )
          )
        )
      )
    },
    server = FALSE
  )

  
  #==============================
  # Output: GSEA Enrichment Plot
  #==============================
  output$enrichmentPlot <- renderPlot({
    req(gsea_result(), input$selected_pathway_gsea)
    
    tryCatch({
      # Get pathway genes
      pathway_genes <- gsea_result()$pathways_list[[input$selected_pathway_gsea]]
      
      if (is.null(pathway_genes) || length(pathway_genes) == 0) {
        plot.new()
        text(0.5, 0.5, "Pathway not found", cex = 1.2)
        return()
      }
      
      # Create enrichment plot
      plotEnrichment(pathway_genes, gsea_result()$ranks_vec) +
        labs(title = input$selected_pathway_gsea) +
        theme_minimal() +
        theme(plot.title = element_text(size = 10))
      
    }, error = function(e) {
      plot.new()
      text(0.5, 0.5, paste("Error creating plot:", e$message), cex = 1)
    })
  })
}