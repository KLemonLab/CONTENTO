explore_contrast_server <- function(input, output, session, state) {
  
  #==============================
  # UI: Contrast dropdown 
  #==============================
  output$contrastSelect <- renderUI({
    req(state$de_df())
    contrast_choices <- unique(state$de_df()$contrast)
    selectInput("contrast", "Select Contrast", choices = contrast_choices)
  })
  
  #==============================
  # Reactive: Filtered and annotated DE data
  #==============================
  selected_data <- reactive({
    req(state$de_df(), input$contrast, input$fc_cutoff)
    
    state$de_df() %>%
      filter(contrast == input$contrast) %>%
      mutate(
        tooltip = paste0(
          symbol, " (", Geneid, ")",
          "\nlog2FC: ", round(log2FC, 2),
          "\nFDR: ", signif(padj, 3)
        ),
        regulated = case_when(
          padj < 0.05 & log2FC >  input$fc_cutoff  ~ "up",
          padj < 0.05 & log2FC < -input$fc_cutoff ~ "down",
          TRUE ~ NA_character_
        ),
        DE = !is.na(regulated)
      )
  })
  
  #==============================
  # Output: DE Table
  #==============================
  output$DETable <- renderDT(
    {
      req(selected_data())
      
      file_base <- if (!is.null(input$deFile) && !is.null(input$deFile$name)) {
        file_path_sans_ext(basename(input$deFile$name))
      } else {
        "contrasts"
      }
      contrast_str <- if (!is.null(input$contrast) && nzchar(input$contrast)) input$contrast else "contrast"
      contrast_str <- gsub("[^A-Za-z0-9._-]+", "__", contrast_str)
      fc_str <- paste0("FC", gsub("\\.", "p", as.character(input$fc_cutoff)))
      file_name <- paste(file_base, contrast_str, fc_str, sep = "__")
      
      df <- selected_data() %>%
        filter(DE) %>%
        mutate(
          across(c(log2FC, log2FC_shrunk), ~ round(.x, 2)),
          padj = formatC(padj, format = "e", digits = 2)
        ) %>%
        arrange(desc(abs(log2FC))) %>%
        select(Geneid, symbol, biotype, log2FC, log2FC_shrunk, padj, regulated)
      
      datatable(
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
              text = 'Download CSV',
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
  # Output: Volcano Plot
  #==============================
  output$volcanoPlot <- renderPlotly({
    req(selected_data())
    gg <- ggplot(selected_data(), aes(x = log2FC_shrunk, y = -log10(padj), text = tooltip)) +
      geom_point(aes(color = DE), alpha = 0.6) +
      scale_color_manual(values = c("TRUE" = "red", "FALSE" = "gray"), guide = "none") +
      geom_vline(xintercept = c(-input$fc_cutoff, input$fc_cutoff), linetype = "dashed") +
      geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
      theme_minimal()
    ggplotly(gg, tooltip = "text")
  })
  
  #==============================
  # Output: Heatmap of Top DE Genes
  #==============================
  output$heatmapPlot <- renderPlot({
    req(selected_data(), state$dds_obj())
    
    top_genes <- selected_data() %>%
      filter(DE) %>%
      slice_max(order_by = abs(log2FC), n = input$top_n)
    
    vsd_mat <- assay(state$dds_obj(), "vst")[rownames(state$dds_obj()) %in% top_genes$Geneid, ]
    vsd_mat <- vsd_mat[match(top_genes$Geneid, rownames(vsd_mat)), ]
    
    rownames(vsd_mat) <- top_genes$symbol
    
    Heatmap(
      scale(vsd_mat),
      col = viridis(100, option = input$viridis_palette),
      column_names_gp = grid::gpar(fontsize = 12, rot = 45),
      row_names_gp = grid::gpar(fontsize = 12),
      heatmap_legend_param = list(title = "Z-scores")
    )
  })
}