explore_contrast_server <- function(input, output, session, state, organism) {
  
  # == == == == == == == == == == == == == == == == == == == == == == == == ==
  #### SECTION 1: UI RENDERING ####
  # == == == == == == == == == == == == == == == == == == == == == == == == ==
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ##### UI: Tab Layout #####
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  output$contrastSubTabs <- renderUI({
    req(state$de_df())
    tabs <- list(
      
      # ── Tab 1: DEGs ────────────────────────────────
      tabPanel("Differentially Expressed Genes",
               fluidRow(
                 column(width = 7,
                        div(style = "padding: 15px; background-color: #f8f9fa; border-radius: 10px; border: 1px solid #e3e6ea;",
                            h4(strong("Explore Differentially Expressed Genes (DEGs)"),
                               style = "margin-top: 0; margin-bottom: 8px;"),
                            tags$ul(
                              style = "margin: 5px 0 0 15px; padding:0;",
                              tags$li("Adjust log2FC column, fold-change and p-value cutoffs in the sidebar; download the full gene list or just the DEGs"),
                              tags$li("Genes are color-coded by regulation direction (up/down-regulated)"),
                              tags$li("Volcano plot highlights DEGs based on your cutoffs; hover over points to view gene details, and download the plot in the interactive viewer")
                            )
                        )
                 ),
                 column(width = 5,
                        div(style = "padding: 15px; background-color: #f8f9fa; border-radius: 10px; border: 1px solid #e3e6ea;",
                            strong("Download Options"),
                            tags$div(
                              style = "margin-top: 15px; display: flex; gap: 10px; justify-content: center; flex-wrap: wrap;",
                              downloadButton("downloadAllContrastData",      "Download All Genes", class = "btn btn-info"),
                              downloadButton("downloadFilteredContrastData", "Download DEGs",      class = "btn btn-info")
                            )
                        )
                 )
               ),
               hr(),
               withSpinner(DTOutput("DETable"), type = 5),
               hr(),
               h4(strong("Volcano Plot of All Genes in the Selected Contrast")),
               withSpinner(plotlyOutput("volcanoPlot", height = "500px", width = "100%"), type = 5)
      ),
      
      # ── Tab 2: Functional Enrichment ────────────────────────────────
      tabPanel("Functional Enrichment",
               fluidRow(
                 column(width = 7,
                        div(style = "padding: 15px; background-color: #f8f9fa; border-radius: 10px; border: 1px solid #e3e6ea;",
                            h4(strong("Perform Gene Set Enrichment Analysis (GSEA)"),
                               style = "margin-top: 0; margin-bottom: 8px;"),
                            tags$ul(
                              style = "margin: 5px 0 0 15px; padding:0;",
                              tags$li("Choose gene set collections from:"),
                              tags$li(style = "margin-left:10px;",
                                      tags$b("MSigDB: "), 
                                      "Curated pathways (e.g. Hallmark, GO, Reactome). ",
                                      "Available for human and mouse through ",
                                      tags$a("MSigDB",
                                             href = "https://www.gsea-msigdb.org/gsea/msigdb",
                                             target = "_blank"),
                                      ", and for additional species via orthology mapping (see ",
                                      tags$a("supported species",
                                             href = "https://igordot.github.io/msigdbr/reference/msigdbr_species.html",
                                             target = "_blank"),
                                      ")."
                              ),
                              tags$li(style = "margin-left:10px;",
                                      tags$b("Annotation: "), 
                                      "Gene sets from your annotation columns (e.g. KEGG, COG)."),
                              tags$li("Gene sets are ranked by NES (Normalized Enrichment Score) and Leading-edge genes indicate core contributors to enrichment"),
                              tags$li("Plot shows top 20 enriched gene sets (FDR < 0.05) based on selected gene set source and collection"),
                              tags$li("Volcano plot highlights DEGs colored by gene set (annotation-based only); hover over points to view gene details, and download the plot in the interactive viewer")
                            )
                        )
                 ),
                 column(width = 5,
                        div(style = "padding: 15px; background-color: #f8f9fa; border-radius: 10px; border: 1px solid #e3e6ea;",
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
               h4(strong("Top 20 Enriched Gene Sets (FDR < 0.05)")),
               withSpinner(plotOutput("gseaTablePlot", height = "600px"), type = 5),
               uiOutput("volcanoGeneSetPlotUI")

      )
    )
    do.call(tabsetPanel, c(list(type = "pills"), tabs))
  })
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ##### UI: Contrast Selector #####
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  output$contrastSelect <- renderUI({
    if (is.null(state$de_df())) return(empty_state_msg())
    div(
      style = "margin-top: -10px;",  
      selectInput(
        "contrast",
        label   = NULL,
        choices = unique(state$de_df()$contrast),
        width   = "90%"
      )
    )
  })
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ##### UI: GSEA Selectors #####
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  output$gseaControlsUI <- renderUI({
    req(state$se_obj())
    gsea_source_ui(
      db_species        = state$db_species(),
      ensembl_col       = state$ensembl_col(),
      avail_cols        = state$available_gsea_columns(),
      organism_name     = organism(),
      source_input_id   = "gsea_source",
      collection_id     = "gs_collection",
      annot_source_id   = "annot_source",
      current_source    = input$gsea_source
    )
  })
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ##### UI: Volcano Plot Colored by Gene Set #####
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  output$volcanoGeneSetPlotUI <- renderUI({
    req(state$se_obj())
    has_annot <- length(state$available_gsea_columns()) > 0
    if (!has_annot) return(NULL)
    
    # Only show once the user has actually selected the annotation source
    # (covers the "both available" case where a radio does exist)
    src <- resolve_gsea_source(input, state, organism(), "gsea_source", "gs_collection", "annot_source")
    req_ok <- !is.null(src) && src$source_type == "annotation"
    if (!req_ok) return(NULL)
    
    tagList(
      hr(),
      h4(strong("Volcano Plot Colored by Gene Set")),
      withSpinner(plotlyOutput("volcanoGeneSetPlot", height = "500px"), type = 5)
    )
  })
  
  # == == == == == == == == == == == == == == == == == == == == == == == == ==
  #### SECTION 2: DATA PROCESSING & ANALYSIS ####
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
  ##### Reactive: DEG Colored by Gene Set #####
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  volcano_geneset_data <- reactive({
    req(selected_data(), genesets())
    req(genesets()$source_type == "annotation")
    
    col <- genesets()$label
    annot <- state$annotation_df()
    req(col %in% colnames(annot))
    
    lookup <- annot |>
      dplyr::filter(!is.na(.data[[col]]) & .data[[col]] != "") |>
      dplyr::select(Geneid, geneset_val = dplyr::all_of(col)) |>
      dplyr::mutate(geneset_val = strsplit(as.character(geneset_val), "!!!")) |>
      tidyr::unnest(geneset_val) |>
      dplyr::mutate(geneset_val = trimws(geneset_val)) |>
      dplyr::filter(geneset_val != "") |>
      dplyr::distinct()
    
    de_only <- selected_data() |> dplyr::filter(DE)
    
    de_colored <- de_only |>
      dplyr::inner_join(lookup, by = "Geneid") |>
      dplyr::mutate(tooltip2 = paste0(tooltip, "\n", col, ": ", geneset_val))
    
    de_na <- de_only |>
      dplyr::filter(!Geneid %in% de_colored$Geneid) |>
      dplyr::mutate(tooltip2 = paste0(tooltip, "\n", col, ": NA"))
    
    list(
      not_de  = selected_data() |> dplyr::filter(!DE),
      de_na   = de_na,
      colored = de_colored
    )
  })
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ##### Reactive: Gene Sets for GSEA #####
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  genesets <- reactive({
    src <- req(resolve_gsea_source(input, state, organism(),
                                   "gsea_source", "gs_collection", "annot_source"))
    tryCatch({
      list(
        data  = do.call(load_genesets, src[names(src) != "label"]),
        label = src$label,
        source_type = src$source_type
      )
    }, error = handle_geneset_error)
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
      
      # Remap gene IDs to Ensembl if rownames are not Ensembl and MSigDB is used
      if (genesets()$source_type == "msigdb" &&
          !is.null(state$ensembl_col()) && !identical(state$ensembl_col(), "rownames")) {
        ranks <- remap_to_ensembl(ranks, state$annotation_df(), state$ensembl_col())
      }
      
      ranks_vec <- setNames(ranks$stat, ranks$Geneid)
      
      # Check for sufficient unique statistics and mixed signs
      if (length(unique(ranks_vec)) < 10) {
        showNotification("Not enough unique statistics for GSEA in this contrast", type = "warning")
        return(NULL)
      }
      if (all(ranks_vec > 0) || all(ranks_vec < 0)) {
        showNotification("All statistics have the same sign; GSEA requires both directions", type = "warning")
        return(NULL)
      }
      
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
        gridExtra::arrangeGrob(
          plotGseaTable(pathways = pathways_list[topPathways], stats = ranks_vec,
                        fgseaRes = fgseaRes, gseaParam = 0.5),
          top = grid::textGrob(
            paste0("Contrast: ", input$contrast, " | Gene Set Collection: ", genesets()$label),
            gp = grid::gpar(fontsize = 14, fontface = "bold")
          )
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
        extensions = "Buttons",
        rownames = FALSE,
        filter = "top",
        options = list(
          pageLength = 15,
          scrollX = TRUE,
          dom = "rtip"  
        )
      )
      
      if ("regulated" %in% colnames(df)) {
        dt <- dt |>
          formatStyle(
            "regulated",
            target = "row",
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
      fc_label <- case_when(
        x_col == "log2FC_shrunk" ~ "log2FC (shrunken)",
        x_col == "log2FC"        ~ "log2FC",
        TRUE                     ~ x_col
      )

      gg <- ggplot(selected_data(), aes(x = .data[[x_col]], y = -log10(padj), text = tooltip)) +
        geom_point(aes(color = regulated), alpha = 0.7) +
        scale_color_manual(values = c("up" = "#00a9cf", "down" = "#d9534f"),
                           na.value = "#b0b0b0", name = NULL) +
        geom_vline(
          xintercept = c(-input$global_log2FC_cutoff, input$global_log2FC_cutoff), linetype = "dashed", color = "#888888") +
        geom_hline(
          yintercept = -log10(input$global_padj_cutoff), linetype = "dashed", color = "#888888") +
        labs(
          title = paste0(
            "Contrast: ", input$contrast, " | ", fc_label, " > ", input$global_log2FC_cutoff, " | FDR cutoff = ", input$global_padj_cutoff),
          x = fc_label,
          y = "-log10(FDR)"
        ) + theme_minimal()
      
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
      grid::grid.draw(gsea_result()$tableplot) 
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
        extensions = "Buttons",
        rownames = FALSE,
        filter = "top",
        options = list(
          pageLength = 15,
          scrollX = TRUE,
          dom = "rtip",
          columnDefs = list(
            list(
              targets = 6,  
              width = "300px",
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
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ##### Output: Volcano Plot Colored by Gene Set ###
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  output$volcanoGeneSetPlot <- renderPlotly({
    req(volcano_geneset_data())
    tryCatch({
      x_col <- if (input$global_lfc_col %in% colnames(selected_data())) input$global_lfc_col else "log2FC"
      d <- volcano_geneset_data()
      has_colored <- nrow(d$colored) > 0
      
      gg <- ggplot() +
        geom_point(data = d$not_de, aes(x = .data[[x_col]], y = -log10(padj), text = tooltip),
                   color = "#d9d9d9", alpha = 0.5) +
        geom_point(data = d$de_na, aes(x = .data[[x_col]], y = -log10(padj), text = tooltip2),
                   color = "#6e6e6e", alpha = 0.7)
      
      if (has_colored) {
        n_cat <- length(unique(d$colored$geneset_val))
        gg <- gg +
          geom_point(data = d$colored, aes(x = .data[[x_col]], y = -log10(padj),
                                           color = geneset_val, text = tooltip2),
                     alpha = 0.8, show.legend = FALSE) +
          scale_color_manual(values = colorRampPalette(brewer.pal(8, "Set2"))(n_cat))
      }
      
      gg <- gg +
        geom_vline(xintercept = c(-input$global_log2FC_cutoff, input$global_log2FC_cutoff), linetype = "dashed", color = "#888888") +
        geom_hline(yintercept = -log10(input$global_padj_cutoff), linetype = "dashed", color = "#888888") +
        labs(title = paste0("DEGs colored by: ", genesets()$label), x = x_col, y = "-log10(FDR)") +
        theme_minimal()
      
      if (!has_colored) {
        showNotification("No DE genes matched any category in this column", type = "warning")
      }
      
      ggplotly(gg, tooltip = "text") |> layout(showlegend = FALSE)
    }, error = handle_volcano_plot_error)
  })

  
  # == == == == == == == == == == == == == == == == == == == == == == == == ==
  #### SECTION 4: DOWNLOAD HANDLERS ####
  # == == == == == == == == == == == == == == == == == == == == == == == == ==
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ##### Download: All Genes for Selected Contrast #####
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
    })
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ##### Download: DE Genes based on Filters #####
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
    })
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ##### Download: GSEA Table Results #####
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
    })
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ##### Download: GSEA Plot Results #####
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  output$downloadGseaPlot <- downloadHandler(
    function() build_download_filename(input, state, ext = "png",
                                       type = "GSEA", contrast = input$contrast,
                                       suffix = genesets()$label),
    
    function(file) {
      req(gsea_result())
      png(file, width = 1800, height = 900, res = 150)
      if (is.null(gsea_result()$tableplot)) {
        plot.new(); text(0.5, 0.5, "No significant gene sets found (FDR < 0.05)", cex = 1.5)
      } else {
        grid::grid.draw(gsea_result()$tableplot)
      }
      dev.off()
    })
  
  # == == == == == == == == == == == == == == == == == == == == == == == == ==
  #### SECTION 5: ERROR HANDLERS ####
  # == == == == == == == == == == == == == == == == == == == == == == == == ==
  
  handle_filter_data_error   <- \(e) { showNotification(paste("Error filtering contrast data:", e$message), type = "error"); NULL }
  handle_volcano_plot_error  <- \(e) { showNotification(paste("Error creating volcano plot:", e$message), type = "error"); NULL }
  handle_geneset_error       <- \(e) { showNotification(paste("Error loading gene sets:", e$message), type = "error"); NULL }
  handle_gsea_error          <- \(e) { showNotification(paste("Error running GSEA:", e$message), type = "error"); NULL }
}