compare_contrast_server <- function(input, output, session, state, organism) {
  
  # == == == == == == == == == == == == == == == == == == == == == == == == ==
  #### SECTION 1: UI RENDERING ####
  # Dynamic UI controls that respond to user inputs and state changes
  # == == == == == == == == == == == == == == == == == == == == == == == == ==
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ##### UI: Contrast Checkbox Selection #####
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  output$multiContrastSelect <- renderUI({
    req(state$de_df())
    contrast_choices <- unique(state$de_df()$contrast)
    
    tagList(
      div(
        style = "display: flex; align-items: center; gap: 10px;",
        strong("Select Contrasts:"),
        actionLink("select_all_contrasts", "Select All"),
        actionLink("clear_all_contrasts", "Clear All")
      ),
      checkboxGroupInput(
        inputId = "compare_contrasts",
        label = NULL,
        choices = contrast_choices,
        selected = NULL 
      )
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
  output$compareGseaControlsUI <- renderUI({
    if (!is.null(organism()) && organism() == "Human") {
      return(tagList(
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
      ))
    } else if (!is.null(organism()) && organism() == "Bacteria") {
      avail_cols <- state$available_gsea_columns()
      
      if (length(avail_cols) == 0) {
        return(div(class = "alert alert-warning",
                   icon("exclamation-triangle"),
                   "No annotation columns available. Upload annotation file."))
      }
      
      return(selectInput("bacterial_geneset_source_compare", "Functional Annotation",
                         choices = avail_cols,
                         selected = avail_cols[[1]]))
    } else {
      return(div(class = "alert alert-info",
                 icon("info-circle"),
                 "Awaiting SE file upload..."))
    }
  })
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ##### UI: Conditional Sub-tabs #####
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  output$compareSubTabs <- renderUI({
    tabs <- list(
      tabPanel("Expression Heatmap",
               uiOutput("heatmapControlsUI"),
               hr(),
               withSpinner(plotOutput("heatmapPlot", height = "700px"), type = 5)),
      tabPanel("DEG Overlap", 
               withSpinner(plotOutput("compareUpsetPlot", height = "500px"), type = 5),  
               h4("Table of DEGs in All Selected Contrasts"),
               withSpinner(DTOutput("compareTable"), type = 5))
    )
    
    if (!is.null(organism()) && organism() %in% c("Human", "Bacteria")) {
      tabs <- append(tabs, 
                     list(
                       tabPanel("GSEA Overlap",
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
                                           withSpinner(plotOutput("compareGSEAPlot", height = "1000px"), type = 5),
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
                                tabsetPanel(
                                  type = "pills",
                                  tabPanel("Overview",
                                           h4("Top 20 GESECA Results"),
                                           withSpinner(plotOutput("gesecaTablePlot", height = "600px"), type = 5),
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
             numericInput("top_n", "Top N genes", value = 50, min = 10, max = 500, step = 10)
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
    req(compare_gsea_data())
    
    pathways <- rownames(compare_gsea_data()$nes_matrix)
    
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
    
    pathways <- geseca_result()$gesecaRes %>%
      arrange(padj) %>%
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
  #### SECTION 2: DATA PROCESSING & FILTERING ####
  # Reactive expressions that filter, transform, and prepare data for display
  # == == == == == == == == == == == == == == == == == == == == == == == == ==
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ##### Reactive: DEG Data for Upset Plot & Heatmap #####
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  compare_data <- reactive({
    req(state$de_df(), input$compare_contrasts, input$global_log2FC_cutoff)
    
    padj_cut <- if (!is.null(input$global_padj_cutoff) && !is.na(input$global_padj_cutoff)) {
      input$global_padj_cutoff
    } else {
      0.05
    }
    lfc_cut <- if (!is.null(input$global_log2FC_cutoff) && !is.na(input$global_log2FC_cutoff)) {
      input$global_log2FC_cutoff
    } else {
      2
    }
    
    deg_list <- map(input$compare_contrasts, function(ct) {
      state$de_df() %>%
        filter(
          contrast == ct,
          !is.na(padj) & !is.na(log2FC),
          padj <= padj_cut,
          abs(log2FC) >= lfc_cut
        ) %>%
        pull(Geneid)
    }) %>% set_names(input$compare_contrasts)
    
    if (all(lengths(deg_list) == 0)) {
      return(NULL)
    }
    
    enframe(deg_list, name = "contrast", value = "Geneid") %>%
      unnest(Geneid) %>%
      mutate(value = TRUE) %>%
      pivot_wider(names_from = contrast, values_from = value, values_fill = FALSE)
  })
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ##### Reactive: DEG Comparison Table Data #####
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  compare_table_data <- reactive({
    req(state$de_df(), input$compare_contrasts, input$global_log2FC_cutoff)
    
    padj_cut <- if (!is.null(input$global_padj_cutoff) && !is.na(input$global_padj_cutoff)) {
      input$global_padj_cutoff
    } else {
      0.05
    }
    lfc_cut <- if (!is.null(input$global_log2FC_cutoff) && !is.na(input$global_log2FC_cutoff)) {
      input$global_log2FC_cutoff
    } else {
      2
    }
    
    has_symbol <- "symbol" %in% colnames(state$de_df())
    id_cols <- if (has_symbol) c("Geneid", "symbol") else "Geneid"
    
    passing_genes <- state$de_df() %>%
      filter(
        contrast %in% input$compare_contrasts,
        !is.na(padj) & !is.na(log2FC),
        padj <= padj_cut,
        abs(log2FC) >= lfc_cut
      ) %>%
      pull(Geneid) %>%
      unique()
    
    if (length(passing_genes) == 0) {
      return(NULL)
    }
    
    filtered <- state$de_df() %>%
      filter(
        contrast %in% input$compare_contrasts,
        Geneid %in% passing_genes
      ) %>%
      select(all_of(c(id_cols, "contrast", "log2FC")))
    
    wide_data <- filtered %>%
      pivot_wider(
        id_cols = all_of(id_cols),
        names_from = contrast,
        values_from = log2FC
      ) %>%
      mutate(across(-any_of(id_cols), ~ round(., 2)))
    
    contrast_cols <- setdiff(names(wide_data), id_cols)
    
    for(col in contrast_cols) {
      wide_data[[paste0(col, "_DE")]] <- !is.na(wide_data[[col]])
    }
    
    col_order <- id_cols
    for(col in contrast_cols) {
      col_order <- c(col_order, col, paste0(col, "_DE"))
    }
    
    wide_data %>%
      select(all_of(col_order)) %>%
      arrange(Geneid)
  })
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ##### Reactive: Full Annotation Export Data #####
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  compare_table_export <- reactive({
    req(compare_table_data(), state$annotation_df())
    
    display_data <- compare_table_data()
    
    display_data %>%
      select(Geneid) %>%
      distinct() %>%
      left_join(
        display_data %>% select(-any_of("symbol")),
        by = "Geneid"
      ) %>%
      left_join(state$annotation_df(), by = "Geneid")
  })
  
  
  # == == == == == == == == == == == == == == == == == == == == == == == == ==
  #### SECTION 3: GENE SET ENRICHMENT ANALYSIS (GSEA) ####
  # Load gene sets and perform multi-contrast GSEA and GESECA
  # == == == == == == == == == == == == == == == == == == == == == == == == ==
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ##### Reactive: Multi-Contrast GSEA Results #####
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  compare_gsea_data <- reactive({
    req(state$de_df(), input$compare_contrasts)
    
    tryCatch({
      if (!is.null(organism()) && organism() == "Bacteria") {
        req(input$bacterial_geneset_source_compare, state$annotation_df())
        pathways_list <- build_bacterial_genesets(state$annotation_df(), input$bacterial_geneset_source_compare)
      } else {
        req(input$compare_gs_collection)
        genesets <- if (!is.null(input$compare_gs_subcollection) && nzchar(input$compare_gs_subcollection)) {
          msigdbr(species = "Homo sapiens", collection = input$compare_gs_collection, subcollection = input$compare_gs_subcollection)
        } else {
          msigdbr(species = "Homo sapiens", collection = input$compare_gs_collection)
        }
        
        pathways_list <- genesets %>%
          split(.$gs_name) %>%
          lapply(function(x) x$ensembl_gene)
      }
      
      fgsea_results <- map(input$compare_contrasts, function(ct) {
        ranks <- state$de_df() %>%
          filter(contrast == ct, !is.na(stat))
        
        if (nrow(ranks) == 0) return(NULL)
        
        ranks_vec <- setNames(ranks$stat, ranks$Geneid)
        
        fgseaMultilevel(
          pathways = pathways_list,
          stats = ranks_vec,
          minSize = 15,
          maxSize = 500
        )
      }) %>% set_names(input$compare_contrasts)
      
      fgsea_results <- fgsea_results[!sapply(fgsea_results, is.null)]
      
      if (length(fgsea_results) == 0) {
        showNotification("No valid GSEA results for selected contrasts", type = "warning")
        return(NULL)
      }
      
      nes_df <- map_dfr(names(fgsea_results), ~ {
        fgsea_results[[.x]] %>%
          as_tibble() %>%
          select(pathway, NES, padj) %>%
          mutate(contrast = .x)
      })
      
      sig_pathways <- nes_df %>%
        filter(padj < 0.05) %>%
        pull(pathway) %>%
        unique()
      
      if (length(sig_pathways) == 0) {
        showNotification("No significant pathways found across contrasts", type = "warning")
        return(NULL)
      }
      
      nes_matrix <- nes_df %>%
        filter(pathway %in% sig_pathways) %>%
        select(pathway, contrast, NES) %>%
        pivot_wider(names_from = contrast, values_from = NES, values_fill = 0) %>%
        column_to_rownames("pathway") %>%
        as.matrix()
      
      padj_matrix <- nes_df %>%
        filter(pathway %in% sig_pathways) %>%
        select(pathway, contrast, padj) %>%
        pivot_wider(names_from = contrast, values_from = padj, values_fill = 1) %>%
        column_to_rownames("pathway") %>%
        as.matrix()
      
      nes_range <- apply(nes_matrix, 1, function(x) max(x) - min(x))
      sort_order <- order(nes_range, decreasing = TRUE)
      
      nes_matrix  <- nes_matrix[sort_order,  , drop = FALSE]
      padj_matrix <- padj_matrix[sort_order, , drop = FALSE]
      
      if (!is.null(input$max_pathways) && nrow(nes_matrix) > input$max_pathways) {
        nes_matrix  <- nes_matrix[1:input$max_pathways,  , drop = FALSE]
        padj_matrix <- padj_matrix[1:input$max_pathways, , drop = FALSE]
      }
      
      list(
        nes_matrix    = nes_matrix,
        padj_matrix   = padj_matrix,
        nes_df        = nes_df %>% filter(pathway %in% sig_pathways),
        fgsea_results = fgsea_results
      )
      
    }, error = handle_compare_gsea_error)
  })
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ##### Reactive: Leading Edge Analysis for Selected Pathway #####
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  leading_edge_data <- reactive({
    req(compare_gsea_data(), input$selected_pathway_compare)
    
    tryCatch({
      fgsea_results <- compare_gsea_data()$fgsea_results
      pathway_name  <- input$selected_pathway_compare
      
      leading_edge_list <- map(
        names(fgsea_results),
        ~ {
          res <- fgsea_results[[.x]] %>%
            as_tibble() %>%
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
      
      le_matrix <- le_df %>%
        mutate(present = 1) %>%
        pivot_wider(names_from = contrast, values_from = present, values_fill = 0) %>%
        column_to_rownames("gene") %>%
        select(-in_leading_edge)
      
      all_genes  <- unique(le_df$gene)
      contrasts  <- names(fgsea_results)
      genes_all  <- rownames(le_matrix)[rowSums(le_matrix) == length(contrasts)]
      
      genes_unique <- map(contrasts, ~ {
        genes_in_contrast <- le_df %>% filter(contrast == .x) %>% pull(gene)
        genes_in_others   <- le_df %>% filter(contrast != .x) %>% pull(gene) %>% unique()
        setdiff(genes_in_contrast, genes_in_others)
      }) %>% set_names(contrasts)
      
      overlap_stats <- list(
        total_genes    = length(all_genes),
        genes_in_all   = length(genes_all),
        genes_in_any   = length(all_genes),
        genes_in_all_list = genes_all,
        genes_unique   = genes_unique
      )
      
      upset_df <- le_matrix %>%
        as.data.frame() %>%
        rownames_to_column("gene") %>%
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
    req(state$se_obj())
    
    tryCatch({
      if (!is.null(organism()) && organism() == "Bacteria") {
        req(input$bacterial_geneset_source_compare, state$annotation_df())
        pathways_list <- build_bacterial_genesets(state$annotation_df(), input$bacterial_geneset_source_compare)
      } else {
        req(input$compare_gs_collection)
        genesets <- if (!is.null(input$compare_gs_subcollection) && nzchar(input$compare_gs_subcollection)) {
          msigdbr(species = "Homo sapiens", collection = input$compare_gs_collection, subcollection = input$compare_gs_subcollection)
        } else {
          msigdbr(species = "Homo sapiens", collection = input$compare_gs_collection)
        }
        
        pathways_list <- genesets %>%
          split(.$gs_name) %>%
          lapply(function(x) x$ensembl_gene)
      }
      
      if (length(pathways_list) == 0) {
        showNotification("No pathways found in selected gene set", type = "warning")
        return(NULL)
      }
      
      vst_matrix <- assays(state$se_obj())[["vst"]]
      
      if (is.null(vst_matrix)) {
        showNotification("VST matrix not found in SummarizedExperiment object", type = "error")
        return(NULL)
      }
      
      gesecaRes <- geseca(
        pathways = pathways_list,
        E        = vst_matrix,
        minSize  = 15,
        maxSize  = 500
      ) %>%
        arrange(padj, pval)
      
      topPathways <- gesecaRes %>%
        filter(padj < 0.05) %>%
        arrange(desc(abs(pctVar))) %>%
        slice_head(n = 20) %>%
        pull(pathway)
      
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
  #### SECTION 4: OUTPUT RENDERING ####
  # Display tables, plots, and interactive visualizations
  # == == == == == == == == == == == == == == == == == == == == == == == == ==
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ##### Output: Expression Heatmap #####
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  output$heatmapPlot <- renderPlot({
    req(compare_data(), state$se_obj())
    req(input$top_n, input$viridis_palette)
    
    all_genes <- compare_data()$Geneid
    
    if (length(all_genes) == 0) {
      plot.new()
      text(0.5, 0.5, "No DE genes found for heatmap\nAdjust cutoffs or select contrasts", cex = 1.5)
      return()
    }
    
    top_genes <- state$de_df() %>%
      filter(
        Geneid %in% all_genes,
        contrast %in% input$compare_contrasts
      ) %>%
      arrange(desc(abs(log2FC))) %>%
      slice_head(n = input$top_n) %>%
      distinct(Geneid, .keep_all = TRUE)
    
    if (nrow(top_genes) == 0) {
      plot.new()
      text(0.5, 0.5, "No genes to display", cex = 1.5)
      return()
    }
    
    vsd_mat <- assay(state$se_obj(), "vst")[rownames(state$se_obj()) %in% top_genes$Geneid, ]
    vsd_mat <- vsd_mat[match(top_genes$Geneid, rownames(vsd_mat)), ]
    
    rownames(vsd_mat) <- if ("symbol" %in% colnames(top_genes)) top_genes$symbol else top_genes$Geneid
    
    Heatmap(
      scale(vsd_mat),
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
    
    upset(
      compare_data(),
      input$compare_contrasts,
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
      req(compare_table_data())
      
      df_display <- compare_table_data()
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
      ) %>%
        formatStyle(
          columns = lfc_cols,
          backgroundColor = styleInterval(
            0,
            c("pink", "lightblue")
          )
        )
    },
    server = FALSE
  )
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ##### Output: GSEA Heatmap #####
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  output$compareGSEAPlot <- renderPlot({
    req(compare_gsea_data())
    
    nes_mat  <- compare_gsea_data()$nes_matrix
    padj_mat <- compare_gsea_data()$padj_matrix
    
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
      req(compare_gsea_data())
      
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
      
      df <- compare_gsea_data()$nes_df %>%
        mutate(
          padj_fmt    = formatC(padj, format = "e", digits = 2),
          NES_display = paste0(round(NES, 2), " (", padj_fmt, ")")
        ) %>%
        select(pathway, contrast, NES_display) %>%
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
      
      gene_symbols <- state$de_df() %>%
        select(Geneid, symbol) %>%
        distinct()
      
      df <- leading_edge_data()$leading_edge_matrix %>%
        as.data.frame() %>%
        rownames_to_column("gene") %>%
        left_join(gene_symbols, by = c("gene" = "Geneid")) %>%
        select(gene, symbol, everything()) %>%
        mutate(across(-c(gene, symbol), ~ ifelse(. == 1, "TRUE", "FALSE")))
      
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
      
      df <- geseca_result()$gesecaRes %>%
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
  #### SECTION 5: DOWNLOAD HANDLERS ####
  # Manage file downloads with proper naming and filtering
  # == == == == == == == == == == == == == == == == == == == == == == == == ==
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ##### Download Handler: DEG Comparison Table (Full Annotations) #####
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  output$downloadCompareTableFull <- downloadHandler(
    filename = function() {
      file_base <- get_download_filename(input, state)
      
      contrasts_str <- if (!is.null(input$compare_contrasts) && length(input$compare_contrasts) > 0) {
        paste(input$compare_contrasts, collapse = "-")
      } else {
        "contrasts"
      }
      contrasts_str <- gsub("[^A-Za-z0-9._-]+", "__", contrasts_str)
      
      lfc_cut <- if (!is.null(input$global_log2FC_cutoff) && !is.na(input$global_log2FC_cutoff)) input$global_log2FC_cutoff else 2
      fc_str  <- paste0("FC", gsub("\\.", "p", as.character(lfc_cut)))
      
      paste(file_base, contrasts_str, fc_str, "full.csv", sep = "__")
    },
    content = function(file) {
      req(input$compareTable_rows_all)
      
      full_data       <- compare_table_export()
      filtered_indices <- input$compareTable_rows_all
      display_data    <- compare_table_data()
      
      filtered_geneids <- display_data[filtered_indices, "Geneid", drop = TRUE]
      
      filtered_full <- full_data %>%
        filter(Geneid %in% filtered_geneids)
      
      write.csv(filtered_full, file, row.names = FALSE)
    }
  )
  
  
  # == == == == == == == == == == == == == == == == == == == == == == == == ==
  #### SECTION 6: ERROR HANDLERS ####
  # Helper functions for error management
  # == == == == == == == == == == == == == == == == == == == == == == == == ==
  
  handle_compare_gsea_error <- \(e) {
    showNotification(paste("Error running multi-contrast GSEA:", e$message), type = "error")
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
