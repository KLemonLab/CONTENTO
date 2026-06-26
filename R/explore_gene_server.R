explore_gene_server <- function(input, output, session, state, organism) {
  
  # == == == == == == == == == == == == == == == == == == == == == == == == ==
  #### SECTION 1: UI RENDERING ####
  # == == == == == == == == == == == == == == == == == == == == == == == == ==
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ##### UI: Tab Layout #####
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  output$geneSubTabs <- renderUI({
    req(state$de_df())
    tabs <- list(
      
      # ── Tab 1: Gene Info ───────────────────────────────────────────
      tabPanel("Gene Info",
               fluidRow(
                 column(width = 12,
                        div(style = "padding: 15px; background-color: #f8f9fa; border-radius: 10px; border: 1px solid #e3e6ea;",
                            h4(strong("Gene Annotation Information"),
                               style = "margin-top: 0; margin-bottom: 8px;"),
                            tags$ul(
                              style = "margin: 5px 0 0 15px; padding:0;",
                              tags$li("Displays all annotation fields available for the selected gene from the loaded annotation file"),
                              tags$li("Fields include genomic coordinates, functional annotations, and any custom columns from your annotation")
                            )
                        )
                 )
               ),
               hr(),
               DTOutput("geneDetails")
      ),
      
      # ── Tab 2: Expression Plot (merged with Gene Table) ─────────────
      tabPanel("Expression Plot",
               fluidRow(
                 column(width = 7,
                        div(style = "padding: 15px; background-color: #f8f9fa; border-radius: 10px; border: 1px solid #e3e6ea;",
                            h4(strong("VST-Normalized Expression Across Conditions"),
                               style = "margin-top: 0; margin-bottom: 8px;"),
                            tags$ul(
                              style = "margin: 5px 0 0 15px; padding:0;",
                              tags$li("Table shows DE results for the selected gene across all contrasts, colored by regulation direction (up/down)"),
                              tags$li("Adjust fold-change and p-value cutoffs in the sidebar to update DE calls in the table"),
                              tags$li("Plot shows VST-normalized expression for the selected gene across all samples"),
                              tags$li("Use X-axis, Color, and Shape controls to group samples by any experimental variable from the SE metadata")

                            )
                        )
                 ),
                 column(width = 5,
                        div(style = "padding: 15px; background-color: #f8f9fa; border-radius: 10px; border: 1px solid #e3e6ea;",
                            strong("Download Options"),
                            tags$div(
                              style = "margin-top: 15px; display: flex; gap: 10px; justify-content: center; flex-wrap: wrap;",
                              downloadButton("downloadGenePlot", "Download Expression Plot for Selected Gene", class = "btn btn-info")
                            )
                        )
                 )
               ),
               hr(),
               h4(strong("DE Results Across All Contrasts")),
               DTOutput("geneContrasts"),
               hr(),
               fluidRow(
                 column(width = 6,
                        h4(strong("Expression Plot"), style = "margin-top: 10px;"),
                        tags$h5(textOutput("geneSymbol"), style = "margin-top: 5px; color: #555;")
                 ),
                 column(width = 6,
                        div(
                          style = "display: flex; justify-content: flex-end; gap: 10px; margin-bottom: 10px; flex-wrap: wrap;",
                          div(style = "min-width: 140px;",
                              selectInput("x_col",     "X-axis", choices = NULL, width = "100%")),
                          div(style = "min-width: 140px;",
                              selectInput("color_col", "Color",  choices = NULL, width = "100%")),
                          div(style = "min-width: 140px;",
                              selectInput("shape_col", "Shape",  choices = NULL, width = "100%"))
                        )
                 )
               ),
               plotOutput("genePlot", height = "600px")
      )
    )
    
    # ── Tab 3 (conditional): Variance Decomposition ──────────────────
    if (!is.null(state$varpart_obj())) {
      tabs <- append(tabs,
        list(tabPanel("Variance Decomposition",
               fluidRow(
                 column(width = 7,
                        div(style = "padding: 15px; background-color: #f8f9fa; border-radius: 10px; border: 1px solid #e3e6ea;",
                            h4(strong("Variance Decomposition for Selected Gene"),
                               style = "margin-top: 0; margin-bottom: 8px;"),
                            tags$ul(
                              style = "margin: 5px 0 0 15px; padding:0;",
                              tags$li("Bar chart shows the fraction of total expression variance attributed to each experimental factor"),
                              tags$li("Variance partition is computed from the full VST expression matrix using variancePartition"),
                              tags$li("Factors with higher bars contribute more to explaining expression differences across samples")
                            )
                        )
                 ),
                 column(width = 5,
                        div(style = "padding: 15px; background-color: #f8f9fa; border-radius: 10px; border: 1px solid #e3e6ea;",
                            strong("Download Options"),
                            tags$div(
                              style = "margin-top: 15px; display: flex; gap: 10px; justify-content: center; flex-wrap: wrap;",
                              downloadButton("downloadVarPartPlot", "Download Variance Plot", class = "btn btn-info")
                            )
                        )
                 )
               ),
               hr(),
               plotOutput("varPartPlot", height = "600px")
        ))
      )
    }
    
    # ── Tab 4 (conditional): Neighbourhood Analysis (Bacteria only) ──
    is_bacteria <- !is.null(organism()) && organism() == "Bacteria"
    if (is_bacteria) {
      annot      <- state$annotation_df()
      has_coords <- !is.null(annot) &&
        all(c("start", "end", "strand") %in% colnames(annot))
      
      neigh_content <- if (has_coords) {
        tagList(
          fluidRow(
            column(width = 7,
                   div(style = "padding: 15px; background-color: #f8f9fa; border-radius: 10px; border: 1px solid #e3e6ea;",
                       h4(tagList(
                         strong("Genomic Neighbourhood Analysis"),
                         tags$span("BETA",
                                   style = "font-size: 0.6em; vertical-align: middle; margin-left: 8px;
                                background-color: #f0ad4e; color: white; padding: 2px 7px;
                                border-radius: 10px; font-weight: bold; letter-spacing: 0.05em;")
                       ), style = "margin-top: 0; margin-bottom: 8px;"),
                       tags$ul(
                         style = "margin: 5px 0 0 15px; padding:0;",
                         tags$li("Displays genes upstream and downstream of the selected gene within a configurable window"),
                         tags$li("Genes are colored by fold-change for the selected contrast; hover to see gene details"),
                         tags$li("Forward and reverse strand genes are shown in separate panels with greedy interval packing to avoid overlap"),
                         tags$li("Requires genomic coordinate columns (start, end, strand) in the annotation file and Organism set to 'Bacteria' in the SE")
                       )
                   )
            ),
            column(width = 5,
                   div(style = "padding: 15px; background-color: #f8f9fa; border-radius: 10px; border: 1px solid #e3e6ea;",
                       div(style = "margin-top: 0;",
                           uiOutput("contrastSelectGene"),
                           numericInput("neigh_window", "Neighbourhood window (nt)",
                                        value = 10000, step = 100, min = 0)
                       )
                   )
            )
          ),
          hr(),
          girafeOutput("neighbourhoodPlot", height = "600px")
        )
      } else {
        tagList(
          fluidRow(
            column(width = 12,
                   div(style = "padding: 15px; background-color: #f8f9fa; border-radius: 10px; border: 1px solid #e3e6ea;",
                       h4(strong("Genomic Neighbourhood Analysis"),
                          style = "margin-top: 0; margin-bottom: 8px;"),
                       tags$ul(
                         style = "margin: 5px 0 0 15px; padding:0;",
                         tags$li("Displays genes upstream and downstream of the selected gene within a configurable window"),
                         tags$li("Genes are colored by log2 fold-change; hover to see gene details"),
                         tags$li("Forward and reverse strand genes are shown in separate panels")
                       )
                   )
            )
          ),
          hr(),
          div(class = "alert alert-warning", style = "margin-top: 20px;",
              icon("exclamation-triangle"),
              strong("Annotation required for Neighbourhood Analysis."),
              tags$p("Upload an annotation file that includes genomic coordinate columns ",
                     code("start"), ", ", code("end"), ", and ", code("strand"),
                     " (e.g. derived from a GFF3 file) to enable this feature.")
          )
        )
      }
      
      tabs <- append(tabs,
                     list(tabPanel("Neighbourhood Analysis", neigh_content)))
    }
    
    do.call(tabsetPanel, c(list(type = "pills"), tabs))
  })
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ##### UI: Gene Select #####
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  output$geneSelect <- renderUI({
    if (is.null(state$de_df())) return(empty_state_msg())
    
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
    
    div(
      style = "margin-top: -10px;",
      selectizeInput(
        "geneSelect",
        label   = NULL,
        choices = c("", choice_vec),
        selected = character(0),
        options = list(
          placeholder  = "Start typing gene name or ID...",
          maxOptions   = 20,
          loadThrottle = 200
        ),
        width = "90%"
      )
    )
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
              options = list(dom = "t", ordering = FALSE, pageLength = nrow(transposed)),
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
                    options = list(dom = "t", ordering = TRUE, pageLength = nrow(gene_contrasts)),
                    rownames = FALSE)
    if ("regulated" %in% colnames(gene_contrasts)) {
      dt <- dt |>
        formatStyle("regulated", target = "row",
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
  #### SECTION 4: DOWNLOAD HANDLERS ####
  # == == == == == == == == == == == == == == == == == == == == == == == == ==
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ##### Download: Expression Plot for Gene #####
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  output$downloadGenePlot <- downloadHandler(
    function() build_download_filename(input, state, ext = "png",
                                       type = "GenePlot", suffix = input$geneSelect),
    function(file) {
      req(input$geneSelect, input$x_col, state$se_obj())
      tryCatch({
        vst_mat         <- assay(state$se_obj(), "vst")
        meta            <- as.data.frame(colData(state$se_obj()))
        meta$expression <- vst_mat[input$geneSelect, ]
        n_colors        <- length(unique(meta[[input$color_col]]))
        palette_colors  <- colorRampPalette(brewer.pal(8, "Dark2"))(n_colors)
        p <- ggplot(meta, aes(.data[[input$x_col]], expression)) +
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
        png(file, width = 1800, height = 1200, res = 150)
        print(p)
        dev.off()
      }, error = function(e) message("Download plot error: ", e$message))
    }
  )

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ##### Download: Variance Decomposition Plot #####
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  output$downloadVarPartPlot <- downloadHandler(
    function() build_download_filename(input, state, ext = "png",
                                       type = "VarPart", suffix = input$geneSelect),
    function(file) {
      req(input$geneSelect, state$varpart_obj())
      vp_gene <- state$varpart_obj()$varPart[input$geneSelect, ]
      vp_df   <- data.frame(Factor = names(vp_gene), Variance = as.numeric(vp_gene))
      p <- ggplot(vp_df, aes(x = reorder(Factor, -Variance), y = Variance)) +
        geom_col(fill = "steelblue") +
        scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
        labs(x = NULL, y = "Fraction of Variance") +
        theme_bw(base_size = 20) +
        theme(axis.text = element_text(angle = 45, hjust = 1),
              panel.grid.major.x = element_blank(),
              panel.grid.minor.x = element_blank())
      png(file, width = 1600, height = 900, res = 150)
      print(p)
      dev.off()
    }
  )
  
  # == == == == == == == == == == == == == == == == == == == == == == == == ==
  #### SECTION 5: ERROR HANDLERS ####
  # == == == == == == == == == == == == == == == == == == == == == == == == ==
  
  handle_gene_filter_error        <- \(e) { showNotification(paste("Error filtering gene data:", e$message), type = "error"); NULL }
  handle_gene_plot_error          <- \(e) { showNotification(paste("Error creating expression plot:", e$message), type = "error"); NULL }
  handle_neighbourhood_plot_error <- \(e) { showNotification(paste("Error creating neighbourhood plot:", e$message), type = "error"); NULL }
  
}