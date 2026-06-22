compare_contrast_server <- function(input, output, session, state, organism) {
  
  # == == == == == == == == == == == == == == == == == == == == == == == == ==
  #### SECTION 1: UI RENDERING ####
  # Dynamic UI controls that respond to user inputs and state changes
  # == == == == == == == == == == == == == == == == == == == == == == == == ==
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ##### UI: Contrast Checkbox Selection #####
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  output$multiContrastSelectHeader <- renderUI({
    req(state$de_df())
    div(
      style = "display: flex; align-items: center; gap: 10px;",
      strong("Select Contrasts:"),
      actionLink("select_all_contrasts", "Select All"),
      actionLink("clear_all_contrasts", "Clear All")
    )
  })
  
  output$multiContrastCheckboxes <- renderUI({
    req(state$de_df())
    contrast_choices <- unique(state$de_df()$contrast)
    checkboxGroupInput(
      inputId = "compare_contrasts",
      label = NULL,
      choices = contrast_choices,
      selected = isolate(input$compare_contrasts)  
    )
  })
  
  observeEvent(input$select_all_contrasts, {
    req(state$de_df())
    contrast_choices <- unique(state$de_df()$contrast)
    updateCheckboxGroupInput(session, "compare_contrasts", selected = contrast_choices)
  })
  
  observeEvent(input$clear_all_contrasts, {
    updateCheckboxGroupInput(session, "compare_contrasts", selected = character(0))
  })
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ##### UI: Dynamic GSEA Controls #####
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Shared helper — called by both tab-specific outputs below
  gsea_controls_ui <- function() {
    if (!is.null(organism()) && organism() == "Human") {
      tagList(
        selectInput("compare_gs_collection", "MSigDB Collection",
                    choices = c("H", "C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8"),
                    selected = "H"),
        conditionalPanel(
          condition = "!(input.compare_gs_collection == 'H' || input.compare_gs_collection == 'C1' || input.compare_gs_collection == 'C6' || input.compare_gs_collection == 'C8')",
          textInput("compare_gs_subcollection", "Subcollection (optional)", placeholder = "e.g., CP:REACTOME")
        ),
        tags$div(style = "margin-top: -10px; margin-bottom: 15px;",
                 tags$a(href = "https://www.gsea-msigdb.org/gsea/msigdb/human/annotate.jsp",
                        target = "_blank",
                        icon("external-link-alt"),
                        "MSigDB"))
      )
    } else if (!is.null(organism()) && organism() == "Bacteria") {
      avail_cols <- state$available_gsea_columns()
      
      if (length(avail_cols) == 0) {
        return(div(class = "alert alert-warning",
                   icon("exclamation-triangle"),
                   "No annotation columns available. Upload annotation file."))
      }
      
      selectInput("bacterial_geneset_source_compare", "",
                  choices  = avail_cols,
                  selected = avail_cols[[1]])
    } else {
      div(class = "alert alert-info",
          icon("info-circle"),
          "Awaiting SE file upload...")
    }
  }
  
  output$compareGseaControlsUI_gsea   <- renderUI({ gsea_controls_ui() })
  output$compareGseaControlsUI_geseca <- renderUI({ gsea_controls_ui() })
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ##### UI: Dynamic Plot Heights Based on Data Size #####
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  output$heatmapPlotUI <- renderUI({
    req(input$top_n)
    h <- max(400, input$top_n * 14)
    plotOutput("heatmapPlot", height = paste0(h, "px"))
  })
  
  output$compareGSEAPlotUI <- renderUI({
    req(gsea_results())
    h <- max(400, nrow(gsea_results()$nes_matrix) * 18)
    plotOutput("compareGSEAPlot", height = paste0(h, "px"))
  })
  
  output$gesecaTablePlotUI <- renderUI({
    req(geseca_result())
    n <- sum(geseca_result()$gesecaRes$padj < 0.05, na.rm = TRUE)
    h <- max(400, min(n, 20) * 40)
    plotOutput("gesecaTablePlot", height = paste0(h, "px"))
  })
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ##### UI: Conditional Sub-tabs #####
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  output$compareSubTabs <- renderUI({
    tabs <- list(
      tabPanel("Expression Heatmap",
               uiOutput("heatmapControlsUI"),
               hr(),
               withSpinner(uiOutput("heatmapPlotUI"), type = 5)),
      tabPanel("DEG Overlap", 
               withSpinner(plotOutput("compareUpsetPlot", height = "500px"), type = 5),  
               h4("Table of DEGs in All Selected Contrasts"),
               withSpinner(DTOutput("compareTable"), type = 5),
               tags$div(
                 style = "margin-top: 15px; display: flex; gap: 10px;",
                 downloadButton("downloadAllContrastsData", "Download All Genes for All Selected Contrasts", class = "btn btn-warning")
               )
      )
    )
    
    if (!is.null(organism()) && organism() %in% c("Human", "Bacteria")) {
      tabs <- append(tabs, 
                     list(
                       tabPanel("GSEA Overlap",
                                hr(style = "margin: 8px 0;"),
                                h5("Gene Sets (GSEA / GESECA)", style = "margin-top: 0; margin-bottom: 8px; font-weight: bold;"),
                                uiOutput("compareGseaControlsUI_gsea"),
                                tabsetPanel(
                                  type = "pills",
                                  tabPanel("Overview",
                                           h4("GSEA Heatmap: Pathways Significant in At Least One Contrast"),
                                           fluidRow(
                                             column(6,
                                                    numericInput("max_pathways", "Max pathways to show", 
                                                                 value = 100, min = 10, max = 500, step = 10)
                                             ),
                                             column(6,
                                                    numericInput("pathway_name_length", "Max pathway name length", 
                                                                 value = 50, min = 25, max = 500, step = 10)
                                             )
                                           ),
                                           hr(),
                                           withSpinner(uiOutput("compareGSEAPlotUI"), type = 5),
                                           hr(),
                                           h4("NES Values for Significant Pathways"),
                                           withSpinner(DTOutput("compareGSEATable"), type = 5)
                                  ),
                                  tabPanel("Pathway Detail",
                                           fluidRow(
                                             column(12,
                                                    uiOutput("pathwaySelectUI_compare"),
                                                    hr(),
                                                    h4("Leading Edge Genes Overlap"),
                                                    withSpinner(plotOutput("leadingEdgeUpsetPlot", height = "500px"), type = 5),
                                                    hr(),
                                                    h4("Leading Edge Summary Statistics"),
                                                    verbatimTextOutput("leadingEdgeSummary"),
                                                    hr(),
                                                    h4("All Leading Edge Genes"),
                                                    withSpinner(DTOutput("leadingEdgeTable"), type = 5)
                                             )
                                           )
                                  )
                                )
                       ),
                       tabPanel("GESECA",
                                hr(style = "margin: 8px 0;"),
                                h5("Gene Sets (GSEA / GESECA)", style = "margin-top: 0; margin-bottom: 8px; font-weight: bold;"),
                                uiOutput("compareGseaControlsUI_geseca"),
                                tabsetPanel(
                                  type = "pills",
                                  tabPanel("Overview",
                                           h4("Top 20 GESECA Results"),
                                           withSpinner(uiOutput("gesecaTablePlotUI"), type = 5),
                                           hr(),
                                           h4("All GESECA Results"),
                                           withSpinner(DTOutput("gesecaResultsTable"), type = 5)
                                  ),
                                  tabPanel("Pathway Detail",
                                           fluidRow(
                                             column(12,
                                                    uiOutput("pathwaySelectUI_geseca"),
                                                    uiOutput("conditionSelectUI_geseca"),
                                                    hr(),
                                                    withSpinner(plotOutput("CoregulationPlot", height = "400px"), type = 5)
                                             )
                                           )
                                  )
                                )
                       )
                     )
      )
    }
    
    do.call(tabsetPanel, c(list(type = "pills"), tabs))
  })

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ##### UI: Heatmap Controls #####
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  output$heatmapControlsUI <- renderUI({
    fluidRow(
      column(6,
             numericInput("top_n", "Top N DE genes", value = 50, min = 10, max = 500, step = 10)
      ),
      column(6,
             selectInput("viridis_palette", "Viridis palette",
                         choices = c("viridis", "magma", "plasma", "inferno", "cividis", "mako", "rocket", "turbo"),
                         selected = "viridis")
      )
    )
  })
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ##### UI: Pathway Selector for GSEA Detail #####
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  output$pathwaySelectUI_compare <- renderUI({
    req(gsea_results())
    
    pathways <- rownames(gsea_results()$nes_matrix)
    
    if (length(pathways) > 0) {
      selectInput("selected_pathway_compare", "Select Pathway:", 
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
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ##### UI: Selectors for GESECA Detail #####
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  output$pathwaySelectUI_geseca <- renderUI({
    req(geseca_result())
    
    pathways <- geseca_result()$gesecaRes |>
      arrange(padj) |>
      pull(pathway)
    
    if (length(pathways) > 0) {
      selectInput("selected_pathway_geseca", "Select Pathway:", 
                  choices = pathways, 
                  selected = pathways[1],
                  width = "100%")
    } else {
      div(
        class = "alert alert-warning",
        icon("exclamation-triangle"),
        "No pathways found"
      )
    }
  })
  
  output$conditionSelectUI_geseca <- renderUI({
    req(state$se_obj())
    
    available_vars <- colnames(colData(state$se_obj()))
    
    default_var <- tryCatch({
      meta <- metadata(state$se_obj())
      if (!is.null(meta$design_formula)) {
        all.vars(as.formula(meta$design_formula))[1]
      } else {
        available_vars[1]
      }
    }, error = function(e) available_vars[1])
    
    fluidRow(
      column(6,
             selectInput("geseca_color_var", "Color by:", 
                         choices = available_vars,
                         selected = default_var)
      ),
      column(6,
             selectInput("geseca_sort_var", "Sort by:", 
                         choices = available_vars,
                         selected = default_var)
      )
    )
  })
  
  
  # == == == == == == == == == == == == == == == == == == == == == == == == ==
  #### SECTION 2: DATA PROCESSING & ANALYSIS ####
  # Reactive expressions that transform and analyze data based on user selections and inputs
  # == == == == == == == == == == == == == == == == == == == == == == == == ==
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ##### Reactive: DEG for the Selected Contrasts with DE based on Global Cutoffs #####
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  selected_contrast_data <- reactive({
    req(state$de_df(), input$compare_contrasts)
    
    tryCatch({
      state$de_df() |>
        filter(contrast %in% input$compare_contrasts) |>
        add_de_flags(input$global_log2FC_cutoff, input$global_padj_cutoff, input$global_lfc_col)
    }, error = handle_filter_contrasts_error)
  })
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ##### Reactive: DEG Comparison Table #####
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  compare_data <- reactive({
    req(selected_contrast_data(), input$compare_contrasts)
    
    tryCatch({
      build_compare_table(selected_contrast_data(), input$compare_contrasts, input$global_lfc_col)
    }, error = handle_build_compare_table_error)
  })
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ##### Reactive: Load Gene Sets (MSigDB or Bacterial) #####
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  genesets_compare <- reactive({
    gs <- if (!is.null(organism()) && organism() == "Bacteria") {
      req(state$annotation_df(), input$bacterial_geneset_source_compare)
      
      tryCatch(
        load_genesets(
          organism(),
          annotation_df    = state$annotation_df(),
          bacterial_source = input$bacterial_geneset_source_compare
        ),
        error = handle_compare_geneset_error
      )
      
    } else {
      req(input$compare_gs_collection)
      
      tryCatch(
        load_genesets(
          organism(),
          gs_collection    = input$compare_gs_collection,
          gs_subcollection = input$compare_gs_subcollection
        ),
        error = handle_compare_msigdb_error
      )
    }
    
    if (length(gs) == 0) {
      showNotification("No pathways found in selected gene set", type = "warning")
      return(NULL)
    }
    
    gs
  })
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ##### Reactive: Multi-Contrast GSEA Analysis #####
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  gsea_results <- reactive({
    req(state$de_df(), input$compare_contrasts, genesets_compare())
    pathways_list <- genesets_compare()
    
    # Run fgsea per contrast (robust handling)
    fgsea_results <- purrr::map(input$compare_contrasts, function(ct) {
      
      ranks <- state$de_df() |>
        dplyr::filter(contrast == ct, !is.na(stat))
      
      if (nrow(ranks) == 0) return(NULL)
      
      ranks_vec <- setNames(ranks$stat, ranks$Geneid)
      
      # Optional safety checks
      if (length(unique(ranks_vec)) < 10) return(NULL)
      if (all(ranks_vec > 0) || all(ranks_vec < 0)) return(NULL)
      
      tryCatch({
        suppressWarnings(
          fgseaMultilevel(
            pathways = pathways_list,
            stats    = ranks_vec,
            minSize  = 15,
            maxSize  = 500
          )
        )
      }, error = function(e) handle_fgsea_error(e, ct))
    })
    
    names(fgsea_results) <- input$compare_contrasts
    
    # Remove failed contrasts cleanly
    fgsea_results <- purrr::compact(fgsea_results)
    
    if (length(fgsea_results) == 0) {
      showNotification("No valid GSEA results for selected contrasts", type = "warning")
      return(NULL)
    }
    
    # Combine results
    nes_df <- purrr::map_dfr(names(fgsea_results), function(ct) {
      fgsea_results[[ct]] |>
        tibble::as_tibble() |>
        dplyr::select(pathway, NES, padj) |>
        dplyr::mutate(contrast = ct)
    })
    
    # Keep only significant pathways
    sig_pathways <- nes_df |>
      dplyr::filter(padj < 0.05) |>
      dplyr::pull(pathway) |>
      unique()
    
    if (length(sig_pathways) == 0) {
      showNotification("No significant pathways found across contrasts", type = "warning")
      return(NULL)
    }
    
    # Build matrices
    nes_matrix <- nes_df |>
      dplyr::filter(pathway %in% sig_pathways) |>
      dplyr::select(pathway, contrast, NES) |>
      tidyr::pivot_wider(names_from = contrast, values_from = NES, values_fill = 0) |>
      tibble::column_to_rownames("pathway") |>
      as.matrix()
    
    padj_matrix <- nes_df |>
      dplyr::filter(pathway %in% sig_pathways) |>
      dplyr::select(pathway, contrast, padj) |>
      tidyr::pivot_wider(names_from = contrast, values_from = padj, values_fill = 1) |>
      tibble::column_to_rownames("pathway") |>
      as.matrix()
    
    # Order by variability
    nes_range <- apply(nes_matrix, 1, function(x) max(x) - min(x))
    sort_order <- order(nes_range, decreasing = TRUE)
    
    nes_matrix  <- nes_matrix[sort_order, , drop = FALSE]
    padj_matrix <- padj_matrix[sort_order, , drop = FALSE]
    
    # Limit pathways if requested
    if (!is.null(input$max_pathways) && nrow(nes_matrix) > input$max_pathways) {
      nes_matrix  <- nes_matrix[1:input$max_pathways, , drop = FALSE]
      padj_matrix <- padj_matrix[1:input$max_pathways, , drop = FALSE]
    }
    
    # Final output
    list(
      nes_matrix    = nes_matrix,
      padj_matrix   = padj_matrix,
      nes_df        = nes_df |> dplyr::filter(pathway %in% sig_pathways),
      fgsea_results = fgsea_results
    )
  })
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ##### Reactive: Leading Edge Analysis for Selected Pathway #####
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  leading_edge_data <- reactive({
    req(gsea_results(), input$selected_pathway_compare)
    
    tryCatch({
      fgsea_results <- gsea_results()$fgsea_results
      pathway_name  <- input$selected_pathway_compare
      
      leading_edge_list <- map(
        names(fgsea_results),
        ~ {
          res <- fgsea_results[[.x]] |>
            as_tibble() |>
            filter(pathway == pathway_name)
          
          if (nrow(res) == 0) {
            return(tibble(contrast = .x, gene = character(), in_leading_edge = logical()))
          }
          
          le_genes <- res$leadingEdge[[1]]
          
          tibble(
            contrast        = .x,
            gene            = le_genes,
            in_leading_edge = TRUE
          )
        }
      )
      
      le_df <- bind_rows(leading_edge_list)
      
      if (nrow(le_df) == 0) {
        showNotification("No leading edge genes found for this pathway", type = "warning")
        return(NULL)
      }
      
      le_matrix <- le_df |>
        mutate(present = 1) |>
        pivot_wider(names_from = contrast, values_from = present, values_fill = 0) |>
        column_to_rownames("gene") |>
        select(-in_leading_edge)
      
      all_genes  <- unique(le_df$gene)
      contrasts  <- names(fgsea_results)
      genes_all  <- rownames(le_matrix)[rowSums(le_matrix) == length(contrasts)]
      
      genes_unique <- map(contrasts, ~ {
        genes_in_contrast <- le_df |> filter(contrast == .x) |> pull(gene)
        genes_in_others   <- le_df |> filter(contrast != .x) |> pull(gene) |> unique()
        setdiff(genes_in_contrast, genes_in_others)
      }) |> set_names(contrasts)
      
      overlap_stats <- list(
        total_genes    = length(all_genes),
        genes_in_all   = length(genes_all),
        genes_in_any   = length(all_genes),
        genes_in_all_list = genes_all,
        genes_unique   = genes_unique
      )
      
      upset_df <- le_matrix |>
        as.data.frame() |>
        rownames_to_column("gene") |>
        mutate(across(-gene, ~ as.logical(.)))
      
      list(
        leading_edge_df     = le_df,
        leading_edge_matrix = le_matrix,
        upset_df            = upset_df,
        overlap_stats       = overlap_stats
      )
      
    }, error = handle_leading_edge_error)
  })
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ##### Reactive: GESECA Analysis #####
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  geseca_result <- reactive({
    req(state$se_obj(), genesets_compare())
    pathways_list <- genesets_compare()
    
    tryCatch({
      # Get VST matrix from SE object
      vst_matrix <- assays(state$se_obj())[["vst"]]
      
      if (is.null(vst_matrix)) {
        showNotification("VST matrix not found in SummarizedExperiment object", type = "error")
        return(NULL)
      }

      # Run GESECA 
      gesecaRes <- geseca(
        pathways = pathways_list,
        E        = vst_matrix,
        minSize  = 15,
        maxSize  = 500
      ) |>
        arrange(padj, pval)
      
      # Select top pathways by variance explained
      topPathways <- gesecaRes |>
        filter(padj < 0.05) |>
        arrange(desc(abs(pctVar))) |>
        slice_head(n = 20) |>
        pull(pathway)
      
      # Generate table plot (only if there are significant pathways)
      tableplot <- if (length(topPathways) > 0) {
        plotGesecaTable(
          gesecaRes = gesecaRes,
          pathways  = pathways_list[topPathways], 
          E         = vst_matrix
        )
      } else {
        NULL
      }
      
      list(tableplot = tableplot, gesecaRes = gesecaRes, vst_matrix = vst_matrix, pathways_list = pathways_list)
      
    }, error = handle_geseca_error)
  })
  
  
  # == == == == == == == == == == == == == == == == == == == == == == == == ==
  #### SECTION 3: OUTPUT RENDERING ####
  # Display tables, plots, and interactive visualizations
  # == == == == == == == == == == == == == == == == == == == == == == == == ==
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ##### Output: Expression Heatmap #####
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  output$heatmapPlot <- renderPlot({
    req(selected_contrast_data(), state$se_obj())
    req(input$top_n, input$viridis_palette)
    
    geneids <- get_top_de_genes(selected_contrast_data(), input$top_n, input$global_lfc_col)
    
    if (length(geneids) == 0) {
      plot.new()
      text(0.5, 0.5, "No DE genes found for heatmap\nAdjust cutoffs or select contrasts", cex = 1.5)
      return()
    }
    
    mat <- get_vst_matrix(state$se_obj(), geneids, selected_contrast_data())
    
    Heatmap(
      scale(mat),
      col = viridis(100, option = input$viridis_palette),
      column_names_gp    = grid::gpar(fontsize = 12),
      row_names_gp       = grid::gpar(fontsize = 10),
      heatmap_legend_param = list(title = "Z-scores"),
      cluster_rows       = TRUE,
      cluster_columns    = TRUE,
      show_row_dend      = TRUE,
      show_column_dend   = TRUE
    )
  })
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ##### Output: DEG Overlap UpSet Plot #####
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  output$compareUpsetPlot <- renderPlot({
    req(compare_data())
    
    df <- compare_data()    
    compare_data_boolean <- df[, sapply(df, is.logical), drop = FALSE]
    
    if (nrow(compare_data_boolean) == 0) {
      plot.new()
      text(0.5, 0.5, "No DE genes found for selected contrasts\nAdjust cutoffs or selected contrasts", cex = 1.5)
      return()
    }
    
    upset(
      compare_data_boolean,
      colnames(compare_data_boolean),
      name = "DEGs",
      min_size = 1,
      base_annotations = list(
        'Intersection size' = intersection_size(
          text = list(size = 5)  
        )
      ),
      themes = upset_default_themes(
        text       = element_text(size = 16),
        axis.title = element_text(size = 16),
        axis.text  = element_text(size = 14)
      )
    )
  })
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ##### Output: DEG Comparison Table #####
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  output$compareTable <- renderDT(
    {
      req(compare_data())
      
      df_display <- compare_data()
      id_cols    <- intersect(c("Geneid", "symbol"), colnames(df_display))
      lfc_cols   <- setdiff(colnames(df_display), id_cols)
      
      datatable(
        df_display,
        extensions = 'Buttons',
        filter = 'top',
        options = list(
          pageLength = 20,
          scrollX = TRUE,
          dom = 'Bfrtip',
          buttons = list(
            list(
              extend = 'csv',
              text = 'Download Filtered (Simple)',
              exportOptions = list(modifier = list(page = "all"))
            )
          )
        ),
        rownames = FALSE
      ) |>
        formatStyle(
          columns = lfc_cols,
          backgroundColor = styleInterval(
            0,
            c("#d9f2f9", "#f8d7da")
          )
        )
    },
    server = FALSE
  )
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ##### Output: GSEA Heatmap #####
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  output$compareGSEAPlot <- renderPlot({
    req(gsea_results())
    
    nes_mat  <- gsea_results()$nes_matrix
    padj_mat <- gsea_results()$padj_matrix
    
    original_names <- rownames(nes_mat)
    if (!is.null(input$pathway_name_length) && input$pathway_name_length > 0) {
      display_names <- ifelse(
        nchar(original_names) > input$pathway_name_length,
        paste0(substr(original_names, 1, input$pathway_name_length), "..."),
        original_names
      )
      rownames(nes_mat)  <- display_names
      rownames(padj_mat) <- display_names
    }
    
    sig_text <- matrix("", nrow = nrow(padj_mat), ncol = ncol(padj_mat))
    sig_text[padj_mat < 0.05] <- "*"
    
    col_fun <- circlize::colorRamp2(
      c(-3, 0, 3),
      c("#009ad1", "#fefbea", "#AD1457")
    )
    
    n_pathways   <- nrow(nes_mat)
    row_fontsize <- max(8, min(12, 400 / n_pathways))
    n_contrasts  <- ncol(nes_mat)
    
    ht <- Heatmap(
      nes_mat,
      name = "NES",
      col  = col_fun,
      
      width        = unit(15, "cm"),
      column_gap   = unit(2, "mm"),
      
      cluster_rows    = FALSE,
      cluster_columns = FALSE,
      show_row_dend   = FALSE,
      show_column_dend = FALSE,
      
      row_names_side      = "left",
      row_names_gp        = grid::gpar(fontsize = row_fontsize),
      row_names_max_width = unit(12, "cm"),
      column_names_gp     = grid::gpar(fontsize = 11),
      column_names_rot    = 45,
      column_names_centered = FALSE,
      
      cell_fun = function(j, i, x, y, width, height, fill) {
        if (sig_text[i, j] == "*") {
          grid::grid.text("*", x, y, gp = grid::gpar(fontsize = 14, col = "black"))
        }
      },
      
      heatmap_legend_param = list(
        title          = "NES",
        direction      = "vertical",
        title_position = "topcenter",
        legend_height  = unit(4, "cm")
      ),
      
      border  = TRUE,
      rect_gp = grid::gpar(col = "grey60", lwd = 0.5),
      
      column_title = if (!is.null(organism()) && organism() == "Bacteria") {
        paste0("GSEA: ", input$bacterial_geneset_source_compare, " (FDR < 0.05)")
      } else {
        paste0("GSEA: ", input$compare_gs_collection, 
               if (!is.null(input$compare_gs_subcollection) && nzchar(input$compare_gs_subcollection)) 
                 paste0(" - ", input$compare_gs_subcollection) else "",
               " (FDR < 0.05)")
      },
      column_title_gp = grid::gpar(fontsize = 14, fontface = "bold")
    )
    
    draw(ht, heatmap_legend_side = "right")
  })
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ##### Output: GSEA Comparison Table #####
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  output$compareGSEATable <- renderDT(
    {
      req(gsea_results())
      
      file_base <- get_download_filename(input, state)
      
      contrasts_str <- if (!is.null(input$compare_contrasts) && length(input$compare_contrasts) > 0) {
        paste(input$compare_contrasts, collapse = "-")
      } else {
        "contrasts"
      }
      contrasts_str <- gsub("[^A-Za-z0-9._-]+", "__", contrasts_str)
      
      if (!is.null(organism()) && organism() == "Bacteria") {
        gs_str <- if (!is.null(input$bacterial_geneset_source_compare)) {
          gsub("[^A-Za-z0-9._-]+", "__", input$bacterial_geneset_source_compare)
        } else {
          "unknown"
        }
      } else {
        gs_str <- paste0(input$compare_gs_collection, if (nzchar(input$compare_gs_subcollection)) paste0("_", input$compare_gs_subcollection) else "")
      }
      file_name <- paste(file_base, contrasts_str, "GSEA", gs_str, sep = "__")
      
      df <- gsea_results()$nes_df |>
        mutate(
          padj_fmt    = formatC(padj, format = "e", digits = 2),
          NES_display = paste0(round(NES, 2), " (", padj_fmt, ")")
        ) |>
        select(pathway, contrast, NES_display) |>
        pivot_wider(
          names_from  = contrast,
          values_from = NES_display,
          values_fill = "NS"
        )
      
      datatable(
        df,
        extensions = 'Buttons',
        rownames   = FALSE,
        filter     = 'top',
        options = list(
          pageLength = 20,
          scrollX = TRUE,
          dom = 'Bfrtip',
          buttons = list(
            list(
              extend = 'csv',
              text = 'Download GSEA Comparison',
              filename = file_name,
              exportOptions = list(modifier = list(page = "all"))
            )
          )
        )
      )
    },
    server = FALSE
  )
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ##### Output: Leading Edge UpSet Plot #####
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  output$leadingEdgeUpsetPlot <- renderPlot({
    req(leading_edge_data())
    
    upset(
      leading_edge_data()$upset_df,
      colnames(leading_edge_data()$leading_edge_matrix),
      name     = "Leading Edge Genes",
      min_size = 1,
      base_annotations = list(
        'Intersection size' = intersection_size(
          text = list(size = 5)
        )
      ),
      themes = upset_default_themes(
        text       = element_text(size = 16),
        axis.title = element_text(size = 16),
        axis.text  = element_text(size = 14)
      )
    )
  })
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ##### Output: Leading Edge Summary #####
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  output$leadingEdgeSummary <- renderText({
    req(leading_edge_data())
    
    stats <- leading_edge_data()$overlap_stats
    
    summary_text <- paste0(
      "Total Leading Edge Genes: ", stats$total_genes, "\n",
      "Genes in ALL Contrasts: ",   stats$genes_in_all, "\n\n",
      "Genes Unique to Each Contrast:\n"
    )
    
    unique_counts <- map_chr(names(stats$genes_unique), ~ {
      paste0("  ", .x, ": ", length(stats$genes_unique[[.x]]))
    })
    
    paste0(summary_text, paste(unique_counts, collapse = "\n"))
  })
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ##### Output: Leading Edge Table #####
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  output$leadingEdgeTable <- renderDT(
    {
      req(leading_edge_data(), state$de_df())
      
      file_base    <- get_download_filename(input, state)
      pathway_str  <- gsub("[^A-Za-z0-9._-]+", "__", input$selected_pathway_compare)
      file_name    <- paste(file_base, pathway_str, "LeadingEdge", sep = "__")
      
      gene_symbols <- state$de_df() |>
        select(any_of(c("Geneid", "symbol"))) |>
        distinct()
      
      df <- leading_edge_data()$leading_edge_matrix |>
        as.data.frame() |>
        rownames_to_column("gene") |>
        left_join(gene_symbols, by = c("gene" = "Geneid")) |>
        select(any_of(c("gene", "symbol")), everything()) |>
        mutate(across(-any_of(c("gene", "symbol")), ~ ifelse(. == 1, "TRUE", "FALSE")))
      
      datatable(
        df,
        extensions = 'Buttons',
        rownames   = FALSE,
        filter     = 'top',
        options = list(
          pageLength = 20,
          scrollX = TRUE,
          dom = 'Bfrtip',
          buttons = list(
            list(
              extend = 'csv',
              text = 'Download Leading Edge Genes',
              filename = file_name,
              exportOptions = list(modifier = list(page = "all"))
            )
          )
        )
      )
    },
    server = FALSE
  )
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ##### Output: GESECA Table Plot #####
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  output$gesecaTablePlot <- renderPlot({
    req(geseca_result())
    
    if (is.null(geseca_result()$tableplot)) {
      plot.new()
      text(0.5, 0.5, "No significant pathways found (FDR < 0.05)", cex = 1.5)
    } else {
      geseca_result()$tableplot
    }
  })
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ##### Output: GESECA Results Table #####
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  output$gesecaResultsTable <- renderDT(
    {
      req(geseca_result())
      
      file_base <- get_download_filename(input, state)
      
      gs_str <- if (!is.null(organism()) && organism() == "Bacteria") {
        gsub("func_", "", input$bacterial_geneset_source_compare)
      } else {
        paste0(input$compare_gs_collection, 
               if (!is.null(input$compare_gs_subcollection) && nzchar(input$compare_gs_subcollection)) 
                 paste0("_", input$compare_gs_subcollection) else "")
      }
      
      file_name <- paste(file_base, "GESECA", gs_str, sep = "__")
      
      df <- geseca_result()$gesecaRes |>
        mutate(
          across(c(pval, padj),       ~ formatC(.x, format = "e", digits = 2)),
          across(c(pctVar, log2err),  ~ round(.x, 3))
        ) 
      
      datatable(
        df,
        extensions = 'Buttons',
        rownames   = FALSE,
        filter     = 'top',
        options = list(
          pageLength = 15,
          scrollX = TRUE,
          dom = 'Bfrtip',
          buttons = list(
            list(
              extend = 'csv',
              text = 'Download Full GESECA Results',
              filename = file_name,
              exportOptions = list(modifier = list(page = "all"))
            )
          )
        )
      )
    },
    server = FALSE
  )
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ##### Output: GESECA Co-regulation Plot #####
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  output$CoregulationPlot <- renderPlot({
    req(geseca_result(), input$selected_pathway_geseca)
    
    tryCatch({
      pathway_genes <- geseca_result()$pathways_list[[input$selected_pathway_geseca]]
      
      if (is.null(pathway_genes) || length(pathway_genes) == 0) {
        plot.new()
        text(0.5, 0.5, "Pathway not found", cex = 1.2)
        return()
      }
      
      color_variable <- if (!is.null(input$geseca_color_var) && nzchar(input$geseca_color_var)) {
        input$geseca_color_var
      } else {
        colnames(colData(state$se_obj()))[1]
      }
      
      sort_variable <- if (!is.null(input$geseca_sort_var) && nzchar(input$geseca_sort_var)) {
        input$geseca_sort_var
      } else {
        colnames(colData(state$se_obj()))[1]
      }
      
      color_conditions <- colData(state$se_obj())[[color_variable]]
      sort_conditions  <- colData(state$se_obj())[[sort_variable]]
      
      sample_order            <- order(sort_conditions)
      vst_sorted              <- geseca_result()$vst_matrix[, sample_order]
      color_conditions_sorted <- color_conditions[sample_order]
      
      plotCoregulationProfile(
        pathway_genes, 
        vst_sorted, 
        conditions = color_conditions_sorted,
        scale      = TRUE
      ) +
        labs(title = input$selected_pathway_geseca) +
        theme_minimal() +
        theme(
          plot.title  = element_text(size = 10),
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
        )
      
    }, error = handle_coregulation_plot_error)
  })
  
  
  # == == == == == == == == == == == == == == == == == == == == == == == == ==
  #### SECTION 4: DOWNLOAD HANDLERS ####
  # Manage file downloads with proper naming and filtering
  # == == == == == == == == == == == == == == == == == == == == == == == == ==
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ##### Download Handler: Download All Genes for Selected Contrast #####
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  output$downloadAllContrastsData <- downloadHandler(
    function() build_download_filename(input, state, 
                                       type = "DE", contrast = input$contrasts),
    
    function(file) {
      req(selected_contrast_data())
      selected_contrast_data() |>
        arrange(desc(abs(.data[[input$global_lfc_col]]))) |>
        select(-any_of(c("tooltip", "DE"))) |>
        write.csv(file, row.names = FALSE)
    }
  )
  
  
  # == == == == == == == == == == == == == == == == == == == == == == == == ==
  #### SECTION 5: ERROR HANDLERS ####
  # Helper functions for error management
  # == == == == == == == == == == == == == == == == == == == == == == == == ==
  
  handle_filter_contrasts_error <- \(e) {
    showNotification(paste("Error filtering contrasts data:", e$message), type = "error")
    NULL
  }
  
  handle_build_compare_table_error <- \(e) {
    showNotification(paste("Error building comparison table:", e$message), type = "error")
    NULL
  }
  
  handle_compare_geneset_error <- \(e) {
    showNotification(paste("Error loading gene sets:", e$message), type = "error")
    NULL
  }
  
  handle_compare_msigdb_error <- \(e) {
    showNotification(paste("Error loading MSigDB:", e$message), type = "error")
    NULL
  }
  
  handle_fgsea_error <- function(e, ct) {
    message("fgsea failed for contrast: ", ct, " | ", e$message)
    NULL
  }
  
  handle_leading_edge_error <- \(e) {
    showNotification(paste("Error extracting leading edge:", e$message), type = "error")
    NULL
  }
  
  handle_geseca_error <- \(e) {
    showNotification(paste("Error running GESECA:", e$message), type = "error")
    NULL
  }
  
  handle_coregulation_plot_error <- \(e) {
    plot.new()
    text(0.5, 0.5, paste("Error creating plot:", e$message), cex = 1)
  }
  
}
