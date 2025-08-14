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
  # Reactive: Filter gene data and calculate DE based on FC cutoff
  #==============================
  selected_gene_data <- reactive({
    req(input$gene_select, state$de_df(), input$gene_fc_cutoff)
    
    state$de_df() %>%
      subset(Geneid == input$gene_select) %>%
      transform(
        regulated = ifelse(padj < 0.05 & log2FC_shrunk >  input$gene_fc_cutoff, "up",
                           ifelse(padj < 0.05 & log2FC_shrunk < -input$gene_fc_cutoff, "down", NA)),
        DE = !is.na(regulated)
      )
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
  # UI: Conditional Sub-tabs
  #==============================
  output$geneSubTabs <- renderUI({
    tabs <- list(
      tabPanel("Expression Plot", plotOutput("genePlot", height = "600px")),
      tabPanel("Gene Table", DTOutput("geneDetails"))
    )
    
    if (!is.null(state$varpart_obj())) {
      tabs <- append(tabs, list(tabPanel("Variance Decomposition", plotOutput("varPartPlot", height = "600px"))))
    }
    
    if (input$organism == "Bacteria") {
      tabs <- append(tabs, list(tabPanel("Neighbourhood Analysis", plotOutput("neighbourhoodPlot", height = "600px"))))
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
    palette_colors <- colorRampPalette(brewer.pal(8, "Dark2"))(n_colors)
    
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
    req(selected_gene_data())
    
    gene_table <- selected_gene_data() %>%
      mutate(across(c(log2FC, log2FC_shrunk), ~ round(.x, 2))) %>%
      mutate(padj = formatC(padj, format = "e", digits = 2)) %>%
      arrange(desc(abs(log2FC_shrunk))) %>%
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
  
  #==============================
  # Output: Neighbourhood Analysis
  #==============================
  output$neighbourhoodPlot <- renderPlot({
    req(selected_gene_data(), state$de_df(), input$gene_select, input$neigh_window)
    
    # Define window around central gene
    central_gene <- selected_gene_data()
    window_start <- central_gene$start[1] - input$neigh_window
    window_end   <- central_gene$end[1] + input$neigh_window
    
    # Filter all DE data in window
    de_df <- state$de_df()
    neigh_genes <- de_df %>%
      subset(start <= window_end & end >= window_start) %>%
      # Apply FC cutoff dynamically
      transform(
        regulated = ifelse(padj < 0.05 & log2FC_shrunk > input$gene_fc_cutoff, "up",
                           ifelse(padj < 0.05 & log2FC_shrunk < -input$gene_fc_cutoff, "down", NA))
      )
    
    req(nrow(neigh_genes) > 0)  # Ensure there is data
    
    # Order genes by start
    neigh_genes <- neigh_genes[order(neigh_genes$start), ]
    
    # Rectangle aesthetics
    neigh_genes$xmin <- neigh_genes$start
    neigh_genes$xmax <- neigh_genes$end
    neigh_genes$ymin <- as.numeric(factor(neigh_genes$contrast)) - 0.4
    neigh_genes$ymax <- as.numeric(factor(neigh_genes$contrast)) + 0.4
    
    # Arrow coordinates based on strand
    neigh_genes$arrow_x <- ifelse(neigh_genes$strand == "+", neigh_genes$xmax, neigh_genes$xmin)
    neigh_genes$arrow_dir <- ifelse(neigh_genes$strand == "+", 1, -1)
    
    ggplot(neigh_genes) +
      geom_rect(aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = log2FC_shrunk), color = "black") +
      geom_segment(aes(x = arrow_x, xend = arrow_x + arrow_dir*50, y = (ymin + ymax)/2, yend = (ymin + ymax)/2),
                   arrow = arrow(length = unit(0.1, "inches")), color = "black") +
      scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, name = "log2FC") +
      scale_y_continuous(breaks = 1:length(unique(neigh_genes$contrast)),
                         labels = unique(neigh_genes$contrast),
                         expand = expansion(add = 0.5)) +
      labs(x = "Genomic Position", y = "Contrast", title = paste("Neighbourhood around", input$gene_select)) +
      theme_bw(base_size = 16)
  })
  
 
  
}