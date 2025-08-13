explore_gene_server <- function(input, output, session, state) {
  
  #==============================
  # UI: Update gene plot selectInputs when DDS loaded
  #==============================
  observe({
    req(state$dds_obj())
    col_vars <- colnames(colData(state$dds_obj()))
    
    updateSelectInput(session, "x_col", choices = col_vars)
    updateSelectInput(session, "color_col", choices = col_vars)
    updateSelectInput(session, "shape_col", choices = col_vars)
  })
  
  #==============================
  # UI: Display gene symbol
  #==============================
  output$geneSymbol <- renderText({
    req(input$gene_select, state$de_df())
    gene_symbol <- state$de_df() %>%
      filter(Geneid == input$gene_select) %>%
      pull(symbol) %>%
      unique()
    if (length(gene_symbol) > 0 && !is.na(gene_symbol)) {
      paste("Gene:", gene_symbol)
    } else {
      "Gene name: not found"
    }
  })
  
  #==============================
  # UI: Conditional Sub-tab for varPartPlot
  #==============================
  output$geneSubTabs <- renderUI({
    tabs <- list(
      tabPanel("Expression Plot", plotOutput("genePlot", height = "600px")),
      tabPanel("Gene Table", DTOutput("geneDetails"))
    )
    if (!is.null(state$varpart_obj())) {
      tabs <- append(
        tabs,
        list(tabPanel("Variance Decomposition", plotOutput("varPartPlot", height = "600px")))
      )
    }
    do.call(tabsetPanel, tabs)
  })
  
  #==============================
  # Output: Gene expression plot
  #==============================
  output$genePlot <- renderPlot({
    req(input$gene_select, input$x_col, state$dds_obj())
    vst_mat <- assay(state$dds_obj(), "vst")
    
    validate(
      need(input$gene_select %in% rownames(vst_mat), "Gene not found in dataset")
    )
    
    expr_values <- vst_mat[input$gene_select, ]
    meta <- as.data.frame(colData(state$dds_obj()))
    meta$expression <- expr_values
    
    n_colors <- length(unique(meta[[input$color_col]]))
    palette_colors <- colorRampPalette(RColorBrewer::brewer.pal(8, "Dark2"))(n_colors)
    
    ggplot(meta, aes_string(x = input$x_col, y = "expression")) +
      geom_boxplot(aes_string(color = input$color_col), outliers = FALSE, show.legend = FALSE) +
      geom_jitter(aes_string(color = input$color_col, shape = input$shape_col),
                  width = 0.2, size = 3, alpha = 0.9) +
      scale_color_manual(values = palette_colors) +
      labs(y = "VST expression", x = input$x_col) +
      theme_bw(base_size = 20) +
      theme(
        axis.text = element_text(angle = 45, hjust = 1),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()
      )
  })
  
  #==============================
  # Output: Gene table for all contrasts
  #==============================
  output$geneDetails <- renderDT({
    req(input$gene_select, state$de_df())
    gene_table <- state$de_df() %>%
      filter(Geneid == input$gene_select) %>%
      mutate(across(c(log2FC, log2FC_shrunk), ~ round(.x, 2))) %>%
      mutate(padj = formatC(padj, format = "e", digits = 2)) %>%
      arrange(desc(abs(log2FC))) %>%
      select(contrast, log2FC, log2FC_shrunk, padj, sign, DE, regulated)
    
    datatable(
      gene_table,
      options = list(dom = 't', ordering = TRUE, pageLength = nrow(gene_table)),
      rownames = FALSE
    ) %>%
      formatStyle(
        'regulated',
        target = 'row',
        backgroundColor = DT::styleEqual(c("up", "down"), c("lightblue", "pink"))
      )
  })
  
  #==============================
  # Output: Variance decomposition plot
  #==============================
  output$varPartPlot <- renderPlot({
    req(input$gene_select, state$varpart_obj())
    gene_id <- input$gene_select
    vp_gene <- state$varpart_obj()$varPart[gene_id, ]
    vp_df <- data.frame(
      Factor = names(vp_gene),
      Variance = as.numeric(vp_gene)
    )
    
    ggplot(vp_df, aes(x = reorder(Factor, -Variance), y = Variance)) +
      geom_col(fill = "steelblue") +
      scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
      labs(x = NULL, y = "Fraction of Variance") +
      theme_bw(base_size = 20) +
      theme(
        axis.text = element_text(angle = 45, hjust = 1),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()
      )
  })
}