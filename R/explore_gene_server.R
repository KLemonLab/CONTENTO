explore_gene_server <- function(input, output, session, state, organism) {
  
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
        regulated = ifelse(padj < 0.05 & log2FC >  input$gene_fc_cutoff, "up",
                           ifelse(padj < 0.05 & log2FC < -input$gene_fc_cutoff, "down", NA)),
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
      tabPanel("Gene Info", DTOutput("geneDetails")),
      tabPanel("Expression Plot", 
               tags$h4(textOutput("geneSymbol"), style = "margin-top: 10px; margin-bottom: 20px;"),
               plotOutput("genePlot", height = "600px")),
      tabPanel("Gene Table", DTOutput("geneContrasts"))
    )
    
    if (!is.null(state$varpart_obj())) {
      tabs <- append(tabs, list(tabPanel("Variance Decomposition", plotOutput("varPartPlot", height = "600px"))))
    }
    
    if (!is.null(organism()) && organism() == "Bacteria") {
      tabs <- append(tabs, list(tabPanel("Neighbourhood Analysis", girafeOutput("neighbourhoodPlot", height = "600px"))))
    }
    
    do.call(tabsetPanel, tabs)
  })
  
  #==============================
  # Output: Gene table Info
  #==============================
  output$geneDetails <- renderDT({
    req(selected_gene_data())
    
    gene_table <- selected_gene_data() %>%
      select(match("biotype", names(.)):match("symbol", names(.))) %>%
      select(symbol, everything()) %>%
      distinct()
    
    # Transpose and convert to data frame
    transposed <- as.data.frame(t(gene_table))
    colnames(transposed) <- "Value"
    transposed$Field <- rownames(transposed)
    transposed <- transposed[, c("Field", "Value")]
    
    datatable(
      transposed,
      options = list(dom = 't', ordering = FALSE, pageLength = nrow(transposed)),
      rownames = FALSE
    )
  })
  
  
  #==============================
  # UI: Contrast dropdown 
  #==============================
  output$contrastSelectGene <- renderUI({
    req(state$de_df())
    contrast_choices <- unique(state$de_df()$contrast)
    selectInput("contrast", "Select Contrast", choices = contrast_choices)
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
  output$geneContrasts <- renderDT({
    req(selected_gene_data())
    
    gene_contrasts <- selected_gene_data() %>%
      mutate(across(c(log2FC, log2FC_shrunk), ~ round(.x, 2))) %>%
      mutate(padj = formatC(padj, format = "e", digits = 2)) %>%
      arrange(desc(abs(log2FC))) %>%
      select(contrast, log2FC, log2FC_shrunk, padj, sign, DE, regulated)
    
    datatable(
      gene_contrasts,
      options = list(dom = 't', ordering = TRUE, pageLength = nrow(gene_contrasts)),
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
  output$neighbourhoodPlot <- renderGirafe({
    req(state$de_df(), input$contrast, input$gene_select, input$neigh_window)
    
    de <- state$de_df()
    
    # Find central gene
    central_gene <- de %>% 
      filter(contrast == input$contrast, Geneid == input$gene_select)
    req(nrow(central_gene) >= 1)
    
    # Define genomic window
    window_start <- central_gene$start[1] - input$neigh_window
    window_end   <- central_gene$end[1] + input$neigh_window
    
    # Filter genes in window
    plot_df <- de %>%
      filter(contrast == input$contrast) %>%
      filter(start <= window_end & end >= window_start) %>%
      mutate(
        is_central = Geneid == input$gene_select,
        tooltip = paste0(
          "GeneID: ", Geneid, "<br>",
          "Symbol: ", symbol, "<br>",
          "log2FC: ", round(log2FC, 2)
        )
      )
    
    # Assign tracks separately per strand
    plot_df <- plot_df %>%
      arrange(strand, start) %>%
      group_by(strand) %>%
      mutate(track = NA_integer_)
    
    for (s in c("+", "-")) {
      strand_rows <- which(plot_df$strand == s)
      tracks <- list()
      for (i in strand_rows) {
        placed <- FALSE
        for (t in seq_along(tracks)) {
          if (plot_df$start[i] > tracks[[t]]) {
            plot_df$track[i] <- t
            tracks[[t]] <- plot_df$end[i]
            placed <- TRUE
            break
          }
        }
        if (!placed) {
          tracks[[length(tracks) + 1]] <- plot_df$end[i]
          plot_df$track[i] <- length(tracks)
        }
      }
    }
    
    plot_df <- ungroup(plot_df) %>%
      mutate(
        strand = factor(strand, levels = c("+", "-"),
                        labels = c("Forward", "Reverse")),
        # Make track a factor with consistent levels
        track = factor(track)
      )
    
    # Create ggplot with interactive rectangles
    gg <- ggplot(plot_df) +
      geom_rect_interactive(aes(
        xmin = start, xmax = end,
        ymin = as.numeric(track) - 0.4,
        ymax = as.numeric(track) + 0.4,
        fill = log2FC,
        tooltip = tooltip
      ), color = "black") +
      scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
      labs(y = NULL, x = "Genomic Position", fill = "log2FC") +
      facet_grid(strand ~ ., scales = "free_y", space = "free_y", drop = TRUE) +
      theme_minimal() +
      theme(
        panel.grid = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        strip.background = element_rect(fill = "grey90", color = "black", size = 1),
        strip.text = element_text(face = "bold", size = 12),
        #panel.border = element_rect(color = "black", fill = NA, size = 1),
        panel.spacing = unit(0.5, "lines")  
      )
    
    # Render interactive plot with tooltips
    girafe(
      ggobj = gg,
      options = list(
        opts_tooltip(opacity = 0.9, offx = 10, offy = -10),
        opts_sizing(rescale = TRUE)
      )
    )
  })
  
}