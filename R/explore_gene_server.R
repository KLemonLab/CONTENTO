explore_gene_server <- function(input, output, session, state, organism) {
  
  # == == == == == == == == == == == == == == == == == == == == == == == == ==
  #### SECTION 1: UI RENDERING ####
  # Dynamic UI controls that respond to user inputs and state changes
  # == == == == == == == == == == == == == == == == == == == == == == == == ==
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ##### UI: Gene Select #####
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  output$geneSelect <- renderUI({
    if (is.null(state$de_df())) return(empty_state_msg())
    
    div(
      style = "margin-top: -10px;",  
      selectizeInput(
        "geneSelect",
        label = NULL,
        choices = NULL, 
        options = list(
          placeholder = 'Start typing gene name or ID...',
          maxOptions = 20,
          loadThrottle = 200
        ),
        width = "90%"
      )
    )
  })
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ##### UI: Update Gene Select Choices #####
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  observe({
    req(state$de_df())
    gene_df <- state$de_df() |>
      select(any_of(c("Geneid", "symbol"))) |>
      distinct()
    
    if ("symbol" %in% colnames(gene_df)) {
      choices_df <- gene_df |>
        mutate(display = ifelse(
          is.na(symbol) | symbol == "" | symbol == Geneid,
          Geneid,
          paste0(symbol, " [", Geneid, "]")
        )) |>
        arrange(display)
      choice_vec <- setNames(choices_df$Geneid, choices_df$display)
    } else {
      choice_vec <- setNames(gene_df$Geneid, gene_df$Geneid)
    }
    
    updateSelectizeInput(session, "geneSelect",
                         choices  = choice_vec,
                         server   = TRUE,
                         selected = character(0))
  })
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ##### UI: Conditional Tab Layout #####
  # Neighbourhood tab: shown for Bacteria only when annotation with
  # genomic coordinates (start/end/strand) is available.
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  output$geneSubTabs <- renderUI({
    req(state$de_df())
    tabs <- list(
      tabPanel("Gene Info",
               DTOutput("geneDetails")),
      tabPanel("Expression Plot",
               selectInput("x_col",     "X-axis", choices = NULL),
               selectInput("color_col", "Color",  choices = NULL),
               selectInput("shape_col", "Shape",  choices = NULL),
               tags$h4(textOutput("geneSymbol"), style = "margin-top: 10px; margin-bottom: 20px;"),
               plotOutput("genePlot", height = "600px")),
      tabPanel("Gene Table",
               DTOutput("geneContrasts"))
    )
    
    if (!is.null(state$varpart_obj())) {
      tabs <- append(tabs,
                     list(tabPanel("Variance Decomposition",
                                   plotOutput("varPartPlot", height = "600px"))))
    }
    
    # Neighbourhood tab intentionally restricted to Bacteria.
    # Show even without annotation — but display a clear message if missing coords
    is_bacteria <- !is.null(organism()) && organism() == "Bacteria"
    if (is_bacteria) {
      annot    <- state$annotation_df()
      has_coords <- !is.null(annot) &&
        all(c("start", "end", "strand") %in% colnames(annot))
      
      neigh_content <- if (has_coords) {
        tagList(
          uiOutput("contrastSelectGene"),
          numericInput("neigh_window", "Neighbourhood window (nt)",
                       value = 10000, step = 100, min = 0),
          girafeOutput("neighbourhoodPlot", height = "600px")
        )
      } else {
        div(class = "alert alert-warning", style = "margin-top: 20px;",
            icon("exclamation-triangle"),
            strong("Annotation required for Neighbourhood Analysis."),
            tags$p("Upload an annotation file that includes genomic coordinate columns ",
                   code("start"), ", ", code("end"), ", and ", code("strand"),
                   " (e.g. derived from a GFF3 file) to enable this feature.")
        )
      }
      
      tabs <- append(tabs,
                     list(tabPanel("Neighbourhood Analysis", neigh_content)))
    }
    
    do.call(tabsetPanel, c(list(type = "pills"), tabs))
  })

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ##### UI: Update Plot Variable Choices #####
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  observe({
    req(state$se_obj())
    col_vars <- colnames(colData(state$se_obj()))
    
    selected_x <- if (!is.null(input$x_col) && input$x_col %in% col_vars) {
      input$x_col
    } else if ("condition" %in% col_vars) {
      "condition"
    } else {
      col_vars[1]
    }
    
    selected_color <- if (!is.null(input$color_col) && input$color_col %in% col_vars) {
      input$color_col
    } else if ("condition" %in% col_vars) {
      "condition"
    } else {
      col_vars[1]
    }
    
    selected_shape <- if (!is.null(input$shape_col) && input$shape_col %in% col_vars) {
      input$shape_col
    } else if ("exp" %in% col_vars) {
      "exp"
    } else {
      col_vars[1]
    }
    
    updateSelectInput(session, "x_col",
                      choices = col_vars,
                      selected = selected_x)
    
    updateSelectInput(session, "color_col",
                      choices = col_vars,
                      selected = selected_color)
    
    updateSelectInput(session, "shape_col",
                      choices = col_vars,
                      selected = selected_shape)
  })
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ##### UI: Contrast Dropdown for Neighbourhood #####
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  output$contrastSelectGene <- renderUI({
    req(state$de_df())
    selectInput("contrast_gene", "Select Contrast", choices = unique(state$de_df()$contrast))
  })
  
  
  # == == == == == == == == == == == == == == == == == == == == == == == == ==
  #### SECTION 2: DATA PROCESSING & ANALYSIS ####
  # Reactive expressions that transform and analyze data based on user selections and inputs
  # == == == == == == == == == == == == == == == == == == == == == == == == ==
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ##### Reactive: DE Filtered by User-Defined Cutoffs #####
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  selected_gene_data <- reactive({
    req(input$geneSelect, state$de_df())
    tryCatch({
      state$de_df() |>
        filter(Geneid == input$geneSelect) |>
        add_de_flags(input$global_log2FC_cutoff, input$global_padj_cutoff, input$global_lfc_col)
    }, error = handle_gene_filter_error)
  })
  
  # == == == == == == == == == == == == == == == == == == == == == == == == ==
  #### SECTION 3: OUTPUT RENDERING ####
  # Display tables, plots, and interactive visualizations
  # == == == == == == == == == == == == == == == == == == == == == == == == ==
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ##### Output: Gene Symbol Text #####
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  output$geneSymbol <- renderText({
    req(input$geneSelect, state$de_df())
    if (!"symbol" %in% colnames(state$de_df())) return(paste("ID:", input$geneSelect))
    sym <- state$de_df() |>
      filter(Geneid == input$geneSelect) |>
      pull(symbol) |> unique()
    if (length(sym) > 0 && !all(is.na(sym))) paste("Gene:", sym[!is.na(sym)][1])
    else "Gene name: not found"
  })
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ##### Output: Gene Info Table #####
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  output$geneDetails <- renderDT({
    req(input$geneSelect, state$de_df()) 
    gene_row <- state$de_df() |>
      filter(Geneid == input$geneSelect) |>
      select(-any_of(c("contrast", "baseMean", "log2FC", "log2FC_shrunk",
                       "lfcSE", "stat", "pvalue", "padj", "regulated", "DE", "tooltip"))) |>
      distinct()
    if (nrow(gene_row) == 0) gene_row <- data.frame(Geneid = input$geneSelect)
    transposed        <- as.data.frame(t(gene_row))
    colnames(transposed) <- "Value"
    transposed$Field  <- rownames(transposed)
    transposed        <- transposed[, c("Field", "Value")]
    datatable(transposed,
              options = list(dom = 't', ordering = FALSE, pageLength = nrow(transposed)),
              rownames = FALSE)
  })
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ##### Output: Gene Expression Plot #####
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  output$genePlot <- renderPlot({
    req(input$geneSelect, input$x_col, state$se_obj())
    tryCatch({
      vst_mat <- assay(state$se_obj(), "vst")
      if (!input$geneSelect %in% rownames(vst_mat)) stop("Gene not found in dataset")
      meta             <- as.data.frame(colData(state$se_obj()))
      meta$expression  <- vst_mat[input$geneSelect, ]
      n_colors         <- length(unique(meta[[input$color_col]]))
      palette_colors   <- colorRampPalette(brewer.pal(8, "Dark2"))(n_colors)
      ggplot(meta, aes(.data[[input$x_col]], expression)) +
        geom_boxplot(aes(color = .data[[input$color_col]]),
                     outliers = FALSE, show.legend = FALSE) +
        geom_jitter(aes(color = .data[[input$color_col]], shape = .data[[input$shape_col]]),
                        width = 0.2, size = 3, alpha = 0.9) +
        scale_color_manual(values = palette_colors) +
        labs(y = "VST expression", x = input$x_col) +
        theme_bw(base_size = 20) +
        theme(axis.text = element_text(angle = 45, hjust = 1),
              panel.grid.major.x = element_blank(),
              panel.grid.minor.x = element_blank())
    }, error = handle_gene_plot_error)
  })
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ##### Output: Gene Contrasts Table #####
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  output$geneContrasts <- renderDT({
    req(selected_gene_data())
    display_cols <- intersect(c("contrast", "log2FC", "log2FC_shrunk", "padj", "DE", "regulated"),
                              colnames(selected_gene_data()))
    gene_contrasts <- selected_gene_data() |>
      mutate(across(any_of(c("log2FC", "log2FC_shrunk")), ~ round(.x, 2)),
             padj = formatC(padj, format = "e", digits = 2)) |>
      arrange(desc(abs(.data[[input$global_lfc_col]]))) |>
      select(all_of(display_cols))
    dt <- datatable(gene_contrasts,
                    options = list(dom = 't', ordering = TRUE, pageLength = nrow(gene_contrasts)),
                    rownames = FALSE)
    if ("regulated" %in% colnames(gene_contrasts)) {
      dt <- dt |>
        formatStyle('regulated', target = 'row',
                    backgroundColor = DT::styleEqual(c("up", "down"), c("#d9f2f9", "#f8d7da")))
    }
    dt
  })
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ##### Output: Variance Decomposition Plot #####
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  output$varPartPlot <- renderPlot({
    req(input$geneSelect, state$varpart_obj())
    vp_gene <- state$varpart_obj()$varPart[input$geneSelect, ]
    vp_df   <- data.frame(Factor = names(vp_gene), Variance = as.numeric(vp_gene))
    ggplot(vp_df, aes(x = reorder(Factor, -Variance), y = Variance)) +
      geom_col(fill = "steelblue") +
      scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
      labs(x = NULL, y = "Fraction of Variance") +
      theme_bw(base_size = 20) +
      theme(axis.text = element_text(angle = 45, hjust = 1),
            panel.grid.major.x = element_blank(),
            panel.grid.minor.x = element_blank())
  })
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ##### Output: Neighbourhood Analysis #####
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  output$neighbourhoodPlot <- renderGirafe({
    req(state$de_df(), input$contrast_gene, input$geneSelect, input$neigh_window)
    tryCatch({
      de           <- state$de_df()
      central_gene <- de |> filter(contrast == input$contrast_gene, Geneid == input$geneSelect)
      req(nrow(central_gene) >= 1)
      
      window_start <- central_gene$start[1] - input$neigh_window
      window_end   <- central_gene$end[1]   + input$neigh_window
      
      plot_df <- de |>
        filter(contrast == input$contrast_gene,
               start <= window_end & end >= window_start) |>
        mutate(
          is_central = Geneid == input$geneSelect,
          tooltip    = paste0(
            "GeneID: ", Geneid, "<br>",
            if ("symbol" %in% colnames(de))
              paste0("Symbol: ", symbol, "<br>")
            else "",
            if (input$global_lfc_col == "log2FC_shrunk") "log2FC (shrunken): " else "log2FC: ",
            round(.data[[input$global_lfc_col]], 2)
          )
        )
      
      # Assign tracks per strand (greedy interval packing)
      plot_df <- plot_df |>
        arrange(strand, start) |>
        group_by(strand) |>
        mutate(track = NA_integer_)
      
      for (s in c("+", "-")) {
        strand_rows <- which(plot_df$strand == s)
        tracks <- list()
        for (i in strand_rows) {
          placed <- FALSE
          for (t in seq_along(tracks)) {
            if (plot_df$start[i] > tracks[[t]]) {
              plot_df$track[i] <- t; tracks[[t]] <- plot_df$end[i]; placed <- TRUE; break
            }
          }
          if (!placed) {
            tracks[[length(tracks) + 1]] <- plot_df$end[i]
            plot_df$track[i] <- length(tracks)
          }
        }
      }
      
      plot_df <- ungroup(plot_df) |>
        mutate(strand = factor(strand, levels = c("+", "-"),
                               labels = c("Forward", "Reverse")),
               track = factor(track))
      
      gg <- ggplot(plot_df) +
        geom_rect_interactive(aes(
          xmin = start, xmax = end,
          ymin = as.numeric(track) - 0.4,
          ymax = as.numeric(track) + 0.4,
          fill = .data[[input$global_lfc_col]], tooltip = tooltip
        ), color = "black") +
        scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
        labs(y = NULL, x = "Genomic Position",
             fill = if (input$global_lfc_col == "log2FC_shrunk") "log2FC (shrunken)" else "log2FC") +
        facet_grid(strand ~ ., scales = "free_y", space = "free_y", drop = TRUE) +
        theme_minimal() +
        theme(panel.grid     = element_blank(),
              axis.title.y   = element_blank(),
              axis.text.y    = element_blank(),
              axis.ticks.y   = element_blank(),
              strip.background = element_rect(fill = "grey90", color = "black", linewidth = 1),
              strip.text     = element_text(face = "bold", size = 12),
              panel.spacing  = unit(0.5, "lines"))
      
      girafe(ggobj = gg,
             options = list(opts_tooltip(opacity = 0.9, offx = 10, offy = -10),
                            opts_sizing(rescale = TRUE)))
    }, error = handle_neighbourhood_plot_error)
  })
  
  
  # == == == == == == == == == == == == == == == == == == == == == == == == ==
  #### SECTION 4: ERROR HANDLERS ####
  # Helper functions for error management
  # == == == == == == == == == == == == == == == == == == == == == == == == ==
  
  handle_gene_filter_error        <- \(e) { showNotification(paste("Error filtering gene data:", e$message), type = "error"); NULL }
  handle_gene_plot_error          <- \(e) { showNotification(paste("Error creating expression plot:", e$message), type = "error"); NULL }
  handle_neighbourhood_plot_error <- \(e) { showNotification(paste("Error creating neighbourhood plot:", e$message), type = "error"); NULL }
  
}