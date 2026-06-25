compare_contrast_server <- function(input, output, session, state, organism) {
  
  # == == == == == == == == == == == == == == == == == == == == == == == == ==
  #### SECTION 1: UI RENDERING ####
  # == == == == == == == == == == == == == == == == == == == == == == == == ==
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ##### UI: Tab Layout #####
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  output$compareSubTabs <- renderUI({
    req(state$de_df())
    db_sp      <- state$db_species()
    ensembl_c  <- state$ensembl_col()
    avail_cols <- state$available_gsea_columns()
    
    has_msigdb   <- !is.null(db_sp) && !is.null(ensembl_c)
    has_annot    <- length(avail_cols) > 0
    has_genesets <- has_msigdb || has_annot
    
    tabs <- list(
      
      # ── Tab 1: DEGs ──────────────────────────────────────────────
      tabPanel("Differentially Expressed Genes",
               fluidRow(
                 column(width = 7,
                        div(style = "padding: 15px; background-color: #f8f9fa; border-radius: 10px; border: 1px solid #e3e6ea;",
                            h4(strong("Compare Differentially Expressed Genes Across Contrasts"),
                               style = "margin-top: 0; margin-bottom: 8px;"),
                            tags$ul(
                              style = "margin: 5px 0 0 15px; padding:0;",
                              tags$li("Adjust log2FC column, fold-change and p-value cutoffs in the sidebar"),
                              tags$li("UpSet plot shows overlap of significant DEGs between selected contrasts"),
                              tags$li("Table includes genes DE in at least one contrast; fold-change (FC) columns are color-coded (up/down-regulated)"),
                              tags$li("Heatmap shows top DE genes across all samples using VST Z-scores")
                            )
                        )
                 ),
                 column(width = 5,
                        div(style = "padding: 15px; background-color: #f8f9fa; border-radius: 10px; border: 1px solid #e3e6ea;",
                            strong("Download Options"),
                            div(style = "margin-top: 15px; display: flex; gap: 10px; justify-content: center; flex-wrap: wrap;",
                                downloadButton("downloadAllContrastsData",
                                               "Download All Genes for Selected Contrasts", class = "btn btn-info"),
                                downloadButton("downloadCompareUpsetPlot",
                                               "Download DEGs UpSet Plot", class = "btn btn-info"),
                                downloadButton("downloadHeatmap",
                                               "Download Expression Heatmap Data", class = "btn btn-info")
                            )
                        )
                 )
               ),
               hr(),
               h4(strong("DEG Overlap")),
               withSpinner(plotOutput("compareUpsetPlot", height = "500px"), type = 5),
               hr(),
               h4(strong("Table of DEGs in All Selected Contrasts")),
               withSpinner(DTOutput("compareTable"), type = 5),
               hr(),
               fluidRow(
                 column(
                   width = 6,
                   h4(strong("Expression Heatmap"), style = "margin-top: 10px;")
                 ),
                 column(
                   width = 6,
                   div(
                     style = "display: flex; justify-content: flex-end; gap: 10px; margin-bottom: 20px;",
                     numericInput(
                       "top_n", "Top N DE genes",
                       value = 50, min = 10, max = 500, step = 10,
                       width = "180px"
                     ),
                     selectInput(
                       "viridis_palette", "Palette",
                       choices = c("viridis", "magma", "plasma", "inferno", "cividis", "mako", "rocket", "turbo"),
                       selected = "viridis",
                       width = "180px"
                     )
                   )
                 )
               ),
               withSpinner(uiOutput("heatmapPlotUI"), type = 5)
      ),
      
      # ── Tab 2: Functional Enrichment ─────────────────────────   
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
                              tags$li("Heatmap shows NES for pathways significant (FDR < 0.05) in at least one contrast"),
                              tags$li("Select a specific gene set to inspect leading-edge gene overlap across contrasts")
                            )
                        )
                 ),
                 column(width = 5,
                        div(style = "padding: 15px; background-color: #f8f9fa; border-radius: 10px; border: 1px solid #e3e6ea;",
                            uiOutput("compareGseaControlsUI_gsea"),
                            div(style = "margin-top: 10px; display: flex; gap: 10px; justify-content: center; flex-wrap: wrap;",
                                downloadButton("downloadCompareGSEAPlot", "Download GSEA Heatmap", class = "btn btn-info"),
                                downloadButton("downloadCompareGseaResults", "Download GSEA Table", class = "btn btn-info")
                            )
                        )
                 )
               ),
               hr(),
               fluidRow(
                 column(
                   width = 6,
                   h4(strong("Heatmap for Gene Sets Significant in At Least One Contrast"),
                      style = "margin-top: 10px;")
                 ),
                 column(
                   width = 6,
                   div(
                     style = "display: flex; justify-content: flex-end; gap: 10px; margin-bottom: 20px;",
                     numericInput(
                       "max_pathways", "Max pathways",
                       value = 100, min = 10, max = 500, step = 10,
                       width = "180px"
                     ),
                     numericInput(
                       "pathway_name_length", "Name length",
                       value = 50, min = 25, max = 500, step = 10,
                       width = "180px"
                     )
                   )
                 )
               ),
               withSpinner(uiOutput("compareGSEAPlotUI"), type = 5),
               hr(),
               h4(strong("Normalized Enrichment Score (NES) for Gene Sets Significant in At Least One Contrast")),
               withSpinner(DTOutput("compareGSEATable"), type = 5),
               hr(),
               fluidRow(
                 column(
                   width = 12,
                   div(
                     style = "display: flex; justify-content: space-between; align-items: center; gap: 15px; flex-wrap: wrap;",
                     div(
                       style = "display: flex; align-items: center; gap: 10px;",
                       h4(strong("Leading Edge Gene Overlap for a Specific Gene Set"),
                          style = "margin: 0;"),
                       downloadButton("downloadLeadingEdgeUpsetPlot", "Download Leading Edge Genes UpSet Plot", class = "btn btn-info")
                     ),
                     div(
                       style = "min-width: 250px;", uiOutput("pathwaySelectUI_compare")
                     )
                   )
                 )
               ),
               withSpinner(plotOutput("leadingEdgeUpsetPlot", height = "500px"), type = 5),
               hr(),
               withSpinner(DTOutput("leadingEdgeTable"), type = 5)
      ),
      
      # ── Tab 3: Co-expression Analysis ─────────────────────────
      tabPanel("Co-Expression Analysis",
               fluidRow(
                 column(width = 7,
                        div(style = "padding: 15px; background-color: #f8f9fa; border-radius: 10px; border: 1px solid #e3e6ea;",
                            h4(strong("Gene Set Co-expression Analysis (GESECA)"),
                               style = "margin-top: 0; margin-bottom: 8px;"),
                            tags$ul(
                              style = "margin: 5px 0 0 15px; padding:0;",
                              tags$li("Select a gene set collection to view its co-regulation profile across samples"),
                              tags$li("GESECA uses the VST expression matrix directly (no contrast statistics required)"),
                              tags$li("Plot inlcudes top 20 significant gene sets (FDR < 0.05)"),
                              tags$li("Select a specific gene set to view its co-regulation profile across samples")

                            )
                        )
                 ),
                 column(width = 5,
                        div(style = "padding: 15px; background-color: #f8f9fa; border-radius: 10px; border: 1px solid #e3e6ea;",
                            uiOutput("compareGseaControlsUI_geseca"),
                            div(style = "margin-top: 10px; display: flex; gap: 10px; justify-content: center; flex-wrap: wrap;",
                                downloadButton("downloadGesecaTablePlot", "Download GESECA Plot", class = "btn btn-info"),
                                downloadButton("downloadGesecaResults", "Download GESECA Table", class = "btn btn-info")
                            )
                        )
                 )
               ),
               hr(),
               h4(strong("Top 20 GESECA Results")),
               withSpinner(uiOutput("gesecaTablePlotUI"), type = 5),
               hr(),
               h4(strong("All GESECA Results")),
               withSpinner(DTOutput("gesecaResultsTable"), type = 5),
               hr(),
               fluidRow(
                 column(
                   width = 12,
                   div(
                     style = "display: flex; justify-content: space-between; align-items: center; gap: 15px; flex-wrap: wrap;",
                     div(
                       style = "display: flex; align-items: center; gap: 10px;",
                       h4(strong("Co-regulation Profile for a Specific Gene Set"),
                          style = "margin: 0;"),
                       downloadButton("downloadCoregulationPlot", "Download Co-regulation Plot", class = "btn btn-info")
                     ),
                     div(
                       style = "display: flex; gap: 10px; flex-wrap: wrap;",
                       div(style = "min-width: 220px;", uiOutput("pathwaySelectUI_geseca")),
                       div(style = "min-width: 220px;", uiOutput("conditionSelectUI_geseca"))
                     )
                   )
                 )
               ),
               withSpinner(plotOutput("CoregulationPlot", height = "400px"), type = 5)
      )
    )
    
    do.call(tabsetPanel, c(list(type = "pills"), tabs))
  })
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ##### UI: Dynamic Plot Heights Based on Data Size #####
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  output$heatmapPlotUI <- renderUI({
    req(input$top_n)
    h <- max(400, input$top_n * 15)
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
  ##### UI: Contrast Checkbox Selectors #####
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  output$multiContrastSelectHeader <- renderUI({
    req(state$de_df())
    div(style = "display: flex; align-items: center; gap: 10px;",
        actionLink("select_all_contrasts", "Select All"),
        actionLink("clear_all_contrasts",  "Clear All"))
  })
  
  output$multiContrastCheckboxes <- renderUI({
    if (is.null(state$de_df())) return(empty_state_msg())
    checkboxGroupInput(
      inputId  = "compare_contrasts",
      label    = NULL,
      choices  = unique(state$de_df()$contrast),
      selected = isolate(input$compare_contrasts)
    )
  })
  
  observeEvent(input$select_all_contrasts, {
    req(state$de_df())
    updateCheckboxGroupInput(session, "compare_contrasts",
                             selected = unique(state$de_df()$contrast))
  })
  
  observeEvent(input$clear_all_contrasts, {
    updateCheckboxGroupInput(session, "compare_contrasts", selected = character(0))
  })
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ##### UI: GSEA/GESECA Selectors #####
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  output$compareGseaControlsUI_gsea <- renderUI({
    req(state$se_obj())
    gsea_source_ui(
      db_species        = state$db_species(),
      ensembl_col       = state$ensembl_col(),
      avail_cols        = state$available_gsea_columns(),
      organism_name     = organism(),
      source_input_id   = "compare_gsea_source",
      collection_id     = "compare_gs_collection",
      annot_source_id   = "compare_annot_source_gsea",
      current_source    = input$compare_gsea_source
    )
  })
  
  output$compareGseaControlsUI_geseca <- renderUI({
    req(state$se_obj())
    gsea_source_ui(
      db_species        = state$db_species(),
      ensembl_col       = state$ensembl_col(),
      avail_cols        = state$available_gsea_columns(),
      organism_name     = organism(),
      source_input_id   = "geseca_gsea_source",
      collection_id     = "geseca_gs_collection",
      annot_source_id   = "geseca_annot_source",
      current_source    = input$geseca_gsea_source
    )
  })
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ##### UI: Individual Pathway Selectors #####
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  output$pathwaySelectUI_compare <- renderUI({
    req(gsea_results())
    pathways <- rownames(gsea_results()$nes_matrix)
    if (length(pathways) > 0) {
      selectInput("selected_pathway_compare", "Select Gene Set:",
                  choices = pathways, selected = pathways[1], width = "100%")
    } else {
      div(class = "alert alert-warning", icon("exclamation-triangle"),
          "No significant gene sets found (FDR < 0.05)")
    }
  })
  
  output$pathwaySelectUI_geseca <- renderUI({
    req(geseca_result())
    pathways <- geseca_result()$gesecaRes |> arrange(padj) |> pull(pathway)
    if (length(pathways) > 0) {
      selectInput("selected_pathway_geseca", "Select Gene Set:",
                  choices = pathways, selected = pathways[1], width = "100%")
    } else {
      div(class = "alert alert-warning", icon("exclamation-triangle"), "No pathways found")
    }
  })
  
  output$conditionSelectUI_geseca <- renderUI({
    req(state$se_obj())
    available_vars <- colnames(colData(state$se_obj()))
    default_var <- tryCatch({
      meta <- metadata(state$se_obj())
      if (!is.null(meta$design_formula)) all.vars(as.formula(meta$design_formula))[1]
      else available_vars[1]
    }, error = function(e) available_vars[1])
    fluidRow(
      column(6, selectInput("geseca_color_var", "Color by:", choices = available_vars, selected = default_var)),
      column(6, selectInput("geseca_sort_var",  "Sort by:",  choices = available_vars, selected = default_var))
    )
  })
  
  # == == == == == == == == == == == == == == == == == == == == == == == == ==
  #### SECTION 2: DATA PROCESSING & ANALYSIS ####
  # == == == == == == == == == == == == == == == == == == == == == == == == ==
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ##### Reactive: DEG for Selected Contrasts #####
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
  ##### Reactive: DEG UpSet Plot #####
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  deg_upset_plot <- reactive({
    req(compare_data())
    df <- compare_data()
    set_cols <- names(df)[sapply(df, is.logical)]
    validate(need(length(set_cols) >= 2, "Need ≥2 contrasts"),
             need(nrow(df) > 0, "No DEGs found"))
    upset(
      df, set_cols,
      name = "DEGs",
      min_size = 1,
      base_annotations = list(
        "Intersection size" = intersection_size(text = list(size = 5))
      ),
      themes = upset_default_themes(
        text = element_text(size = 16),
        axis.title = element_text(size = 16),
        axis.text  = element_text(size = 14)
      )
    )
  })
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ##### Reactive: Expression Heatmap #####
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  heatmap_plot <- reactive({
    req(selected_contrast_data(), state$se_obj(), input$top_n, input$viridis_palette)
    
    geneids <- get_top_de_genes(selected_contrast_data(), input$top_n, input$global_lfc_col)
    
    if (length(geneids) == 0) {
      return(function() {
        plot.new()
        text(0.5, 0.5,
             "No DE genes found for heatmap\nAdjust cutoffs or select contrasts",
             cex = 1.5)
      })
    }
    mat <- subset_scale_vst_matrix(state$se_obj(), geneids, selected_contrast_data())
    hm <- Heatmap(
      mat,
      col = viridis(100, option = input$viridis_palette),
      column_names_gp = grid::gpar(fontsize = 12),
      row_names_gp    = grid::gpar(fontsize = 10),
      heatmap_legend_param = list(title = "Z-scores"),
      cluster_rows = TRUE,
      cluster_columns = TRUE
    )
    
    function() draw(hm)
  })  
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ##### Reactive: Gene Sets for GSEA #####
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  genesets_compare <- reactive({
    src <- req(resolve_gsea_source(input, state, organism(),
                                   "compare_gsea_source", "compare_gs_collection", "compare_annot_source_gsea"))
    tryCatch({
      gs <- do.call(load_genesets, src[names(src) != "label"])
      if (length(gs) == 0) {
        showNotification("No gene sets found in selected collection", type = "warning")
        return(NULL)
      }
      list(data = gs, label = src$label)
    }, error = handle_compare_geneset_error)
  })
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ##### Reactive: Gene Sets for GESECA #####
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  genesets_geseca <- reactive({
    src <- req(resolve_gsea_source(input, state, organism(),
                                   "geseca_gsea_source", "geseca_gs_collection", "geseca_annot_source"))
    tryCatch({
      gs <- do.call(load_genesets, src[names(src) != "label"])
      if (length(gs) == 0) {
        showNotification("No gene sets found in selected collection", type = "warning")
        return(NULL)
      }
      list(data = gs, label = src$label)
    }, error = handle_compare_geneset_error)
  })
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ##### Reactive: Multi-Contrast GSEA Results #####
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  gsea_results <- reactive({
    req(state$de_df(), input$compare_contrasts, genesets_compare())
    pathways_list <- genesets_compare()$data
    
    fgsea_results <- purrr::map(input$compare_contrasts, function(ct) {
      ranks <- state$de_df() |>
        dplyr::filter(contrast == ct, !is.na(stat))
      if (nrow(ranks) == 0) return(NULL)
      ranks_vec <- setNames(ranks$stat, ranks$Geneid)
      if (length(unique(ranks_vec)) < 10) return(NULL)
      if (all(ranks_vec > 0) || all(ranks_vec < 0)) return(NULL)
      tryCatch(
        suppressWarnings(fgseaMultilevel(pathways = pathways_list, stats = ranks_vec,
                                         minSize = 15, maxSize = 500)),
        error = function(e) handle_fgsea_error(e, ct)
      )
    })
    names(fgsea_results) <- input$compare_contrasts
    fgsea_results <- purrr::compact(fgsea_results)
    
    if (length(fgsea_results) == 0) {
      showNotification("No valid GSEA results for selected contrasts", type = "warning")
      return(NULL)
    }
    
    nes_df <- purrr::map_dfr(names(fgsea_results), function(ct) {
      fgsea_results[[ct]] |>
        tibble::as_tibble() |>
        dplyr::select(pathway, NES, padj) |>
        dplyr::mutate(contrast = ct)
    })
    
    sig_pathways <- nes_df |> dplyr::filter(padj < 0.05) |> dplyr::pull(pathway) |> unique()
    
    if (length(sig_pathways) == 0) {
      showNotification("No significant gene sets found across contrasts", type = "warning")
      return(NULL)
    }
    
    make_wide_matrix <- function(value_col, fill_val) {
      nes_df |>
        dplyr::filter(pathway %in% sig_pathways) |>
        dplyr::select(pathway, contrast, !!value_col) |>
        tidyr::pivot_wider(names_from = contrast, values_from = !!value_col,
                           values_fill = fill_val) |>
        tibble::column_to_rownames("pathway") |>
        as.matrix()
    }
    
    nes_mat  <- make_wide_matrix("NES", 0)
    padj_mat <- make_wide_matrix("padj", 1)
    
    sort_order <- order(apply(nes_mat, 1, function(x) max(x) - min(x)), decreasing = TRUE)
    nes_mat  <- nes_mat[sort_order, , drop = FALSE]
    padj_mat <- padj_mat[sort_order, , drop = FALSE]
    
    if (!is.null(input$max_pathways) && nrow(nes_mat) > input$max_pathways) {
      nes_mat  <- nes_mat[seq_len(input$max_pathways), , drop = FALSE]
      padj_mat <- padj_mat[seq_len(input$max_pathways), , drop = FALSE]
    }
    
    list(nes_matrix = nes_mat, padj_matrix = padj_mat,
         nes_df = nes_df |> dplyr::filter(pathway %in% sig_pathways),
         fgsea_results = fgsea_results)
  })
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ##### Reactive: GSEA Heatmap #####
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  gsea_heatmap <- reactive({
    req(gsea_results())
    
    nes_mat  <- gsea_results()$nes_matrix
    padj_mat <- gsea_results()$padj_matrix
    
    if (!is.null(input$pathway_name_length) && input$pathway_name_length > 0) {
      orig <- rownames(nes_mat)
      disp <- ifelse(nchar(orig) > input$pathway_name_length,
                     paste0(substr(orig, 1, input$pathway_name_length), "..."),
                     orig)
      rownames(nes_mat)  <- disp
      rownames(padj_mat) <- disp
    }
    
    sig_text <- matrix("", nrow = nrow(padj_mat), ncol = ncol(padj_mat))
    sig_text[padj_mat < 0.05] <- "*"
    
    col_fun <- circlize::colorRamp2(c(-3, 0, 3),
                                    c("#009ad1", "#fefbea", "#AD1457"))
    
    n_pathways   <- nrow(nes_mat)
    row_fontsize <- max(8, min(12, 400 / n_pathways))
    
    hm <- Heatmap(
      nes_mat, name = "NES", col = col_fun,
      cluster_rows = FALSE, cluster_columns = FALSE,
      show_row_dend = FALSE, show_column_dend = FALSE,
      row_names_gp = grid::gpar(fontsize = row_fontsize),
      column_names_gp = grid::gpar(fontsize = 11),
      column_names_rot = 45,
      cell_fun = function(j, i, x, y, width, height, fill) {
        if (sig_text[i, j] == "*") {
          grid::grid.text("*", x, y,
                          gp = grid::gpar(fontsize = 14, col = "black"))
        }
      },
      heatmap_legend_param = list(title = "NES")
    )
    function() draw(hm, heatmap_legend_side = "right")
  })
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ##### Reactive: Leading Edge Results #####
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  leading_edge_data <- reactive({
    req(gsea_results(), input$selected_pathway_compare)
    
    # Guard: silently wait if the selected pathway is stale (e.g. from a previous
    # gene set collection). input$selected_pathway_compare will update shortly
    # after pathwaySelectUI_compare re-renders with the new gsea_results().
    available_pathways <- unique(gsea_results()$nes_df$pathway)
    req(input$selected_pathway_compare %in% available_pathways)
    
    tryCatch({
      fgsea_results <- gsea_results()$fgsea_results
      pathway_name  <- input$selected_pathway_compare
      contrasts     <- names(fgsea_results)
      
      le_list <- map(contrasts, ~ {
        res <- fgsea_results[[.x]] |> as_tibble() |> filter(pathway == pathway_name)
        # leadingEdge may be missing entirely OR an empty list — both yield character(0)
        genes <- if (nrow(res) == 0) character(0) else res$leadingEdge[[1]]
        tibble(contrast = .x, gene = genes)
      })
      le_df <- bind_rows(le_list)
      
      if (nrow(le_df) == 0) {
        showNotification("No leading edge genes found for this gene set", type = "warning")
        return(NULL)
      }
      
      # Build one column per contrast explicitly so contrasts with zero
      # leading-edge genes still get an all-zero column in the matrix
      all_genes <- unique(le_df$gene)
      le_matrix <- map_dfc(contrasts, function(ct) {
        genes_in_ct <- le_df |> filter(contrast == ct) |> pull(gene)
        tibble(!!ct := as.integer(all_genes %in% genes_in_ct))
      }) |>
        as.data.frame() |>
        `rownames<-`(all_genes)
      
      # ncol(le_matrix) == length(contrasts) by construction, but using ncol
      # makes the intent explicit and is robust if contrasts ever changes shape
      genes_all    <- rownames(le_matrix)[rowSums(le_matrix) == ncol(le_matrix)]
      genes_unique <- map(contrasts, ~ {
        setdiff(le_df |> filter(contrast == .x) |> pull(gene),
                le_df |> filter(contrast != .x) |> pull(gene) |> unique())
      }) |> set_names(contrasts)
      
      upset_df <- le_matrix |>
        rownames_to_column("gene") |>
        mutate(across(-gene, as.logical))
      
      list(leading_edge_df = le_df, leading_edge_matrix = le_matrix, upset_df = upset_df,
           overlap_stats = list(total_genes = length(all_genes), genes_in_all = length(genes_all),
                                genes_in_any = length(all_genes),
                                genes_in_all_list = genes_all, genes_unique = genes_unique))
    }, error = handle_leading_edge_error)
  })
 
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ##### Reactive: Leading Edge Plot #####
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  leading_edge_upset_plot <- reactive({
    req(leading_edge_data())
    df <- leading_edge_data()$upset_df
    set_cols <- names(df)[sapply(df, is.logical)]
    validate(need(length(set_cols) >= 2, "Need ≥2 contrasts"),
             need(nrow(df) > 0, "No leading edge genes found"))
    upset(
      df, set_cols,
      name = "Leading Edge Genes",
      min_size = 1,
      base_annotations = list(
        "Intersection size" = intersection_size(text = list(size = 5))
      ),
      themes = upset_default_themes(
        text = element_text(size = 16),
        axis.title = element_text(size = 16),
        axis.text  = element_text(size = 14)
      )
    )
  })
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ##### Reactive: GESECA Results #####
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  geseca_result <- reactive({
    req(state$se_obj(), genesets_geseca())
    pathways_list <- genesets_geseca()$data 
    
    tryCatch({
      vst_matrix <- assays(state$se_obj())[["vst"]]
      if (is.null(vst_matrix)) {
        showNotification("VST matrix not found in SummarizedExperiment object", type = "error")
        return(NULL)
      }
      
      gesecaRes <- geseca(pathways = pathways_list, E = vst_matrix,
                          minSize = 15, maxSize = 500) |>
        arrange(padj, pval)
      
      topPathways <- gesecaRes |> filter(padj < 0.05) |>
        arrange(desc(abs(pctVar))) |> slice_head(n = 20) |> pull(pathway)
      
      tableplot <- if (length(topPathways) > 0) {
        plotGesecaTable(gesecaRes = gesecaRes,
                        pathways  = pathways_list[topPathways],
                        E         = vst_matrix)
      } else { NULL }
      
      list(tableplot = tableplot, gesecaRes = gesecaRes,
           vst_matrix = vst_matrix, pathways_list = pathways_list)
    }, error = handle_geseca_error)
  })
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ##### Reactive: GESECA Table Plot #####
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  geseca_table_plot <- reactive({
    req(geseca_result())
    
    tp <- geseca_result()$tableplot
    
    if (is.null(tp)) {
      return(function() {
        plot.new()
        text(0.5, 0.5,
             "No significant gene sets found (FDR < 0.05)",
             cex = 1.5)
      })
    }
    
    function() print(tp)
  })
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ##### Reactive: GESECA Co-regulation Plot #####
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  coregulation_plot <- reactive({
    req(geseca_result(), input$selected_pathway_geseca)
    
    pathway_genes <- geseca_result()$pathways_list[[input$selected_pathway_geseca]]
    req(!is.null(pathway_genes), length(pathway_genes) > 0)
    
    color_var <- input$geseca_color_var %||% colnames(colData(state$se_obj()))[1]
    sort_var  <- input$geseca_sort_var  %||% colnames(colData(state$se_obj()))[1]
    
    color_cond <- colData(state$se_obj())[[color_var]]
    sort_cond  <- colData(state$se_obj())[[sort_var]]
    
    sample_order <- order(sort_cond)
    vst_sorted   <- geseca_result()$vst_matrix[, sample_order]
    
    plotCoregulationProfile(
      pathway_genes, vst_sorted,
      conditions = color_cond[sample_order],
      scale = TRUE
    ) +
      labs(title = input$selected_pathway_geseca) +
      theme_minimal()
  })
  
  
  # == == == == == == == == == == == == == == == == == == == == == == == == ==
  #### SECTION 3: OUTPUT RENDERING ####
  # == == == == == == == == == == == == == == == == == == == == == == == == ==
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ##### Output: DEG Comparison Table #####
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  output$compareTable <- renderDT({
    req(compare_data())
    df_display <- compare_data()
    id_cols  <- intersect(c("Geneid", "symbol"), colnames(df_display))
    lfc_cols <- grep("^FC_", colnames(df_display), value = TRUE)
    datatable(df_display, extensions = "Buttons", filter = "top",
              options = list(pageLength = 20, scrollX = TRUE, dom = "Bfrtip",
                             buttons = list(list(extend = "csv", text = "Download DEGs (genes can be filtered based on DE condition (true/false) across different contrasts)",
                                                 exportOptions = list(modifier = list(page = "all"))))),
              rownames = FALSE) |>
      formatStyle(columns = lfc_cols,
                  backgroundColor = styleInterval(0, c("#d9f2f9", "#f8d7da")))
  }, server = FALSE)
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ##### Output: GSEA Comparison Table #####
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  output$compareGSEATable <- renderDT({
    req(gsea_results())
    df <- gsea_results()$nes_df |>
      mutate(padj_fmt    = formatC(padj, format = "e", digits = 2),
             NES_display = paste0(round(NES, 2), " (", padj_fmt, ")")) |>
      select(pathway, contrast, NES_display) |>
      pivot_wider(names_from = contrast, values_from = NES_display, values_fill = "NS")
    datatable(df, rownames = FALSE,
              options = list(pageLength = 20, scrollX = TRUE, dom = "rtip"))
  }, server = FALSE)

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ##### Output: Leading Edge Comparison Table #####
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  output$leadingEdgeTable <- renderDT({
    req(leading_edge_data(), state$de_df())
    pathway_str <- gsub("[^A-Za-z0-9._-]+", "_", input$selected_pathway_compare)
    filename <- build_download_filename(input, state, type = "LeadingEdge", ext = NULL,
                            contrast = input$compare_contrasts,
                            suffix   = paste(genesets_compare()$label, pathway_str, sep = "__"))
    gene_symbols <- state$de_df() |> select(any_of(c("Geneid", "symbol"))) |> distinct()
    df <- leading_edge_data()$leading_edge_matrix |>
      as.data.frame() |> rownames_to_column("gene") |>
      left_join(gene_symbols, by = c("gene" = "Geneid")) |>
      select(any_of(c("gene", "symbol")), everything()) |>
      mutate(across(-any_of(c("gene", "symbol")), ~ ifelse(. == 1, "TRUE", "FALSE")))
    datatable(df, extensions = "Buttons", rownames = FALSE, filter = "top",
              options = list(pageLength = 20, scrollX = TRUE, dom = "Bfrtip",
                             buttons = list(list(extend = "csv", text = "Download Leading Edge Genes (genes can be filtered based on presence (true/false) across different contrasts)",
                                                 filename = filename,
                                                 exportOptions = list(modifier = list(page = "all"))))))
  }, server = FALSE)
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ##### Output: GESECA Results Table #####
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  output$gesecaResultsTable <- renderDT({
    req(geseca_result())
    df <- geseca_result()$gesecaRes |>
      mutate(across(c(pval, padj),      ~ formatC(.x, format = "e", digits = 2)),
             across(c(pctVar, log2err), ~ round(.x, 3)))
    datatable(df, rownames = FALSE,
              options = list(pageLength = 20, scrollX = TRUE, dom = "rtip"))
  }, server = FALSE)
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ##### Output: Plots from Reactives #####
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  output$heatmapPlot <- renderPlot({heatmap_plot()()})
  output$compareUpsetPlot <- renderPlot({deg_upset_plot()})
  output$compareGSEAPlot <- renderPlot({gsea_heatmap()()})
  output$leadingEdgeUpsetPlot <- renderPlot({leading_edge_upset_plot()})
  output$gesecaTablePlot <- renderPlot({geseca_table_plot()()})
  output$CoregulationPlot <- renderPlot({coregulation_plot()})
  
  # == == == == == == == == == == == == == == == == == == == == == == == == ==
  #### SECTION 4: DOWNLOAD HANDLERS ####
  # == == == == == == == == == == == == == == == == == == == == == == == == ==
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ##### Download: All Genes for Selected Contrast #####
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  output$downloadAllContrastsData <- downloadHandler(
    function() build_download_filename(input, state,
                                       type = "DE", contrast = input$compare_contrasts),
    function(file) {
      req(selected_contrast_data())
      selected_contrast_data() |>
        arrange(desc(abs(.data[[input$global_lfc_col]]))) |>
        select(-any_of(c("tooltip", "DE"))) |>
        write.csv(file, row.names = FALSE)
    }
  )
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ##### Download: DEG UpSet Plot #####
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  output$downloadCompareUpsetPlot <- downloadHandler(
    function() build_download_filename(input, state, ext = "png",
                                       type = "UpSet", contrast = input$compare_contrasts,
                                       filters = list(FC = input$global_log2FC_cutoff, FDR = input$global_padj_cutoff)),
    function(file) {
      png(file, width = 1800, height = 900, res = 150)
      print(deg_upset_plot())
      dev.off()
    }
  )
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ##### Download: Expression Heatmap #####
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  output$downloadHeatmap <- downloadHandler(
    function() build_download_filename(input, state, ext = "png",
                                       type = "Heatmap", contrast = input$compare_contrasts,
                                       filters = list(FC = input$global_log2FC_cutoff, FDR = input$global_padj_cutoff)),
    function(file) {
      req(input$top_n)
      h_px <- max(600, input$top_n * 15 + 200)
      png(file, width = 1800, height = h_px, res = 150)
      heatmap_plot()()
      dev.off()
    }
  )
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ##### Download: GSEA Comparison Table #####
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  output$downloadCompareGseaResults <- downloadHandler(
    function() {
      build_download_filename(input, state,
                              type = "GSEA", contrast = input$compare_contrasts,
                              suffix = genesets_compare()$label)
    },
    function(file) {
      req(gsea_results())
      gsea_results()$nes_df |>
        mutate(padj = formatC(padj, format = "e", digits = 2)) |>
        write.csv(file, row.names = FALSE)
    }
  )
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ##### Download: GSEA Heatmap #####
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  output$downloadCompareGSEAPlot <- downloadHandler(
    function() {
      build_download_filename(input, state, ext = "png",
                              type = "GSEA", contrast = input$compare_contrasts,
                              suffix = genesets_compare()$label)
    },
    function(file) {
      png(file, width = 1800, height = 900, res = 150)
      gsea_heatmap()()
      dev.off()
    }
  )
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ##### Download: Leading Edge UpSet Plot #####
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  output$downloadLeadingEdgeUpsetPlot <- downloadHandler(
    function() {
      pathway_str <- gsub("[^A-Za-z0-9._-]+", "_", input$selected_pathway_compare)
      build_download_filename(
        input, state,
        ext      = "png",
        type     = "LeadingEdge",
        contrast = input$compare_contrasts,
        suffix   = paste(genesets_compare()$label, pathway_str, sep = "__")
      )
    },
    function(file) {
      png(file, width = 1800, height = 900, res = 150)
      print(leading_edge_upset_plot())
      dev.off()
    }
  )
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ##### Download: GESECA Table #####
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  output$downloadGesecaResults <- downloadHandler(
    function() {
      build_download_filename(input, state,
                              type = "GESECA",
                              suffix = genesets_geseca()$label)
    },
    function(file) {
      req(geseca_result())
      geseca_result()$gesecaRes |>
        mutate(across(c(pval, padj),      ~ formatC(.x, format = "e", digits = 2)),
               across(c(pctVar, log2err), ~ round(.x, 3))) |>
        write.csv(file, row.names = FALSE)
    }
  )
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ##### Download: GESECA Plot #####
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  output$downloadGesecaTablePlot <- downloadHandler(
    function() {
      build_download_filename(input, state, ext = "png",
                              type = "GESECA",
                              suffix = genesets_geseca()$label)
    },
    function(file) {
      req(geseca_result())
      n   <- sum(geseca_result()$gesecaRes$padj < 0.05, na.rm = TRUE)
      h_px <- max(600, min(n, 20) * 40 + 200)
      png(file, width = 1800, height = h_px, res = 150)
      geseca_table_plot()()
      dev.off()
    }
  )
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ##### Download: GESECA Co-regulation Plot #####
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  output$downloadCoregulationPlot <- downloadHandler(
    function() {
      pathway_str <- gsub("[^A-Za-z0-9._-]+", "_", input$selected_pathway_geseca)
      build_download_filename(
        input, state,
        ext    = "png",
        type   = "GESECA",
        suffix = paste(genesets_geseca()$label, pathway_str, sep = "__")
      )
    },
    function(file) {
      png(file, width = 1800, height = 600, res = 150)
      print(coregulation_plot())
      dev.off()
    }
  )
  
  # == == == == == == == == == == == == == == == == == == == == == == == == ==
  #### SECTION 5: ERROR HANDLERS ####
  # == == == == == == == == == == == == == == == == == == == == == == == == ==
  
  handle_filter_contrasts_error    <- \(e) { showNotification(paste("Error filtering contrast data:", e$message), type = "error"); NULL }
  handle_build_compare_table_error <- \(e) { showNotification(paste("Error building comparison table:", e$message), type = "error"); NULL }
  handle_compare_geneset_error     <- \(e) { showNotification(paste("Error loading gene sets:", e$message), type = "error"); NULL }
  handle_fgsea_error <- function(e, ct) { message("fgsea failed for contrast: ", ct, " | ", e$message); NULL }
  handle_leading_edge_error        <- \(e) { showNotification(paste("Error extracting leading edge:", e$message), type = "error"); NULL }
  handle_geseca_error              <- \(e) { showNotification(paste("Error running GESECA:", e$message), type = "error"); NULL }
  handle_coregulation_plot_error   <- \(e) { plot.new(); text(0.5, 0.5, paste("Error creating plot:", e$message), cex = 1) }
  
}
