compare_contrast_server <- function(input, output, session, state) {
  
  #==============================
  # UI: Contrast selection
  #==============================
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
  
  #==============================
  # UI: Conditional Sub-tabs
  #==============================
  output$compareSubTabs <- renderUI({
    tabs <- list(
      tabPanel("DEG Overlap", 
               withSpinner(plotOutput("compareUpsetPlot", height = "500px"), type = 5),  
               h4("Table of DEGs in All Selected Contrasts"),
               withSpinner(DTOutput("compareTable"), type = 5))
    )
    
    if (input$organism == "Human") {
      tabs <- append(tabs, 
                     list(
                       tabPanel("Gene Sets Overlap",
                                tabsetPanel(
                                  tabPanel("Overview",
                                           h4("GSEA Heatmap: Pathways Significant in At Least One Contrast"),
                                           withSpinner(plotOutput("compareGSEAPlot", height = "700px"), type = 5),
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
                       )
                     )
      )
    }
    
    do.call(tabsetPanel, tabs)
  }) 
    
  #==============================
  # Select All / Clear All actions
  #==============================
  observeEvent(input$select_all_contrasts, {
    req(state$de_df())
    contrast_choices <- unique(state$de_df()$contrast)
    updateCheckboxGroupInput(session, "compare_contrasts", selected = contrast_choices)
  })
  
  observeEvent(input$clear_all_contrasts, {
    updateCheckboxGroupInput(session, "compare_contrasts", selected = character(0))
  })
  
  #==============================
  # UI: Pathway selector for comparison
  #==============================
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
  
  #==============================
  # Reactive: Data for DEGs upset plot 
  #==============================
  compare_data <- reactive({
    req(state$de_df(), input$compare_contrasts)
    
    deg_list <- map(input$compare_contrasts, function(ct) {
      state$de_df() %>%
        filter(
          contrast == ct,
          padj < 0.05,
          abs(log2FC) > input$compare_fc_cutoff
        ) %>%
        pull(Geneid)
    }) %>% set_names(input$compare_contrasts)
    
    enframe(deg_list, name = "contrast", value = "Geneid") %>%
      unnest(Geneid) %>%
      mutate(value = TRUE) %>%
      pivot_wider(names_from = contrast, values_from = value, values_fill = FALSE)
  })
  
  #==============================
  # Reactive: Table of DEGs in all selected contrasts
  #==============================
  compare_table_data <- reactive({
    req(state$de_df(), input$compare_contrasts)
    
    # Filter DE results for selected contrasts and cutoff
    filtered <- state$de_df() %>%
      filter(
        contrast %in% input$compare_contrasts,
        padj < 0.05,
        abs(log2FC) > input$compare_fc_cutoff
      ) %>%
      select(Geneid, symbol, contrast, log2FC)
    
    # Keep only genes that appear in all selected contrasts
    genes_in_all <- filtered %>%
      group_by(Geneid) %>%
      summarize(n_contrasts = n_distinct(contrast), .groups = "drop") %>%
      filter(n_contrasts == length(input$compare_contrasts)) %>%
      pull(Geneid)
    
    # Filter again and pivot to wide format
    filtered %>%
      filter(Geneid %in% genes_in_all) %>%
      pivot_wider(
        id_cols = c(Geneid, symbol),  
        names_from = contrast,
        values_from = log2FC
      ) %>%
      mutate(across(-c(Geneid, symbol), ~ round(., 2))) %>%
      arrange(Geneid)
  })
  
  #==============================
  # Reactive: GSEA Results for Multiple Contrasts
  #==============================
  compare_gsea_data <- reactive({
    req(state$de_df(), input$compare_contrasts, input$compare_gs_collection)
    
    tryCatch({
      # Load gene sets
      genesets <- if (nzchar(input$compare_gs_subcollection)) {
        msigdbr(species = "Homo sapiens", collection = input$compare_gs_collection, subcollection = input$compare_gs_subcollection)
      } else {
        msigdbr(species = "Homo sapiens", collection = input$compare_gs_collection)
      }
      
      # Convert to pathway list
      pathways_list <- genesets %>%
        split(.$gs_name) %>%
        lapply(function(x) x$ensembl_gene)
      
      # Run GSEA for each selected contrast
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
      
      # Remove NULL results
      fgsea_results <- fgsea_results[!sapply(fgsea_results, is.null)]
      
      if (length(fgsea_results) == 0) {
        showNotification("No valid GSEA results for selected contrasts", type = "warning")
        return(NULL)
      }
      
      # Extract NES and padj for all pathways
      nes_df <- map_dfr(names(fgsea_results), ~ {
        fgsea_results[[.x]] %>%
          as_tibble() %>%
          select(pathway, NES, padj) %>%
          mutate(contrast = .x)
      })
      
      # Find pathways significant in at least one contrast
      sig_pathways <- nes_df %>%
        filter(padj < 0.05) %>%
        pull(pathway) %>%
        unique()
      
      if (length(sig_pathways) == 0) {
        showNotification("No significant pathways found across contrasts", type = "warning")
        return(NULL)
      }
      
      # Create NES matrix
      nes_matrix <- nes_df %>%
        filter(pathway %in% sig_pathways) %>%
        select(pathway, contrast, NES) %>%
        pivot_wider(names_from = contrast, values_from = NES, values_fill = 0) %>%
        column_to_rownames("pathway") %>%
        as.matrix()
      
      # Create padj matrix
      padj_matrix <- nes_df %>%
        filter(pathway %in% sig_pathways) %>%
        select(pathway, contrast, padj) %>%
        pivot_wider(names_from = contrast, values_from = padj, values_fill = 1) %>%
        column_to_rownames("pathway") %>%
        as.matrix()
      
      list(
        nes_matrix = nes_matrix,
        padj_matrix = padj_matrix,
        nes_df = nes_df %>% filter(pathway %in% sig_pathways),
        fgsea_results = fgsea_results
      )
      
    }, error = function(e) {
      showNotification(paste("Error running multi-contrast GSEA:", e$message), type = "error")
      NULL
    })
  })
  
  #==============================
  # Reactive: Leading Edge Analysis for Selected Pathway
  #==============================
  leading_edge_data <- reactive({
    req(compare_gsea_data(), input$selected_pathway_compare)
    
    tryCatch({
      fgsea_results <- compare_gsea_data()$fgsea_results
      pathway_name <- input$selected_pathway_compare
      
      # Extract leading edge genes for each contrast
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
            contrast = .x,
            gene = le_genes,
            in_leading_edge = TRUE
          )
        }
      )
      
      # Combine all contrasts
      le_df <- bind_rows(leading_edge_list)
      
      if (nrow(le_df) == 0) {
        showNotification("No leading edge genes found for this pathway", type = "warning")
        return(NULL)
      }
      
      # Create presence/absence matrix
      le_matrix <- le_df %>%
        mutate(present = 1) %>%
        pivot_wider(names_from = contrast, values_from = present, values_fill = 0) %>%
        column_to_rownames("gene") %>%
        select(-in_leading_edge)
      
      # Calculate overlap statistics
      all_genes <- unique(le_df$gene)
      contrasts <- names(fgsea_results)
      
      genes_all <- rownames(le_matrix)[rowSums(le_matrix) == length(contrasts)]
      
      # Genes unique to each contrast
      genes_unique <- map(contrasts, ~ {
        genes_in_contrast <- le_df %>% filter(contrast == .x) %>% pull(gene)
        genes_in_others <- le_df %>% filter(contrast != .x) %>% pull(gene) %>% unique()
        setdiff(genes_in_contrast, genes_in_others)
      }) %>% set_names(contrasts)
      
      overlap_stats <- list(
        total_genes = length(all_genes),
        genes_in_all = length(genes_all),
        genes_in_any = length(all_genes),
        genes_in_all_list = genes_all,
        genes_unique = genes_unique
      )
      
      # Convert to upset format
      upset_df <- le_matrix %>%
        as.data.frame() %>%
        rownames_to_column("gene") %>%
        mutate(across(-gene, ~ as.logical(.)))
      
      list(
        leading_edge_df = le_df,
        leading_edge_matrix = le_matrix,
        upset_df = upset_df,
        overlap_stats = overlap_stats
      )
      
    }, error = function(e) {
      showNotification(paste("Error extracting leading edge:", e$message), type = "error")
      NULL
    })
  })
  
  #==============================
  # Output: DEGs Upset plot
  #==============================
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
        text = element_text(size = 16),       # general text
        axis.title = element_text(size = 16), # axis titles
        axis.text = element_text(size = 14)   # axis labels
      )
    )
  })
  
  #==============================
  # Output: DEGs Comparison Table
  #==============================
  output$compareTable <- renderDT(
    {
      req(compare_table_data())
      
      file_base <- if (!is.null(input$deFile) && !is.null(input$deFile$name)) {
        file_path_sans_ext(basename(input$deFile$name))
      } else {
        "contrasts"
      }
      
      contrasts_str <- if (!is.null(input$compare_contrasts) && length(input$compare_contrasts) > 0) {
        paste(input$compare_contrasts, collapse = "-")
      } else {
        "contrasts"
      }
      contrasts_str <- gsub("[^A-Za-z0-9._-]+", "__", contrasts_str)
      
      fc_str <- paste0("FC", gsub("\\.", "p", as.character(input$fc_cutoff)))
      
      file_name <- paste(file_base, contrasts_str, fc_str, sep = "__")
      
      df <- compare_table_data()
      lfc_cols <- setdiff(colnames(df), c("Geneid", "symbol"))
      
      datatable(
        df,
        extensions = 'Buttons',
        options = list(
          pageLength = 20,
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
  
  #==============================
  # Output: GSEA Heatmap
  #==============================
  output$compareGSEAPlot <- renderPlot({
    req(compare_gsea_data())
    
    nes_mat <- compare_gsea_data()$nes_matrix
    padj_mat <- compare_gsea_data()$padj_matrix
    
    # Create significance markers
    sig_markers <- ifelse(padj_mat < 0.05, "*", "")
    
    # Create custom color palette
    my_colors <- colorRampPalette(c("#009ad1", "#fefbea", "#AD1457"))(100)
    
    pheatmap::pheatmap(
      nes_mat,
      color = my_colors,
      breaks = seq(-3, 3, length.out = 101),
      cluster_cols = FALSE,
      cluster_rows = FALSE,
      fontsize_row = 10,
      fontsize_col = 12,
      main = paste0("GSEA: ", input$compare_gs_collection, 
                    if (nzchar(input$compare_gs_subcollection)) paste0(" - ", input$compare_gs_subcollection) else "",
                    " (FDR < 0.05)"),
      border_color = "grey60",
      display_numbers = sig_markers,
      number_color = "black",
      fontsize_number = 14
    )
  })
  
  #==============================
  # Output: GSEA Comparison Table
  #==============================
  output$compareGSEATable <- renderDT(
    {
      req(compare_gsea_data())
      
      # Prepare file name
      file_base <- if (!is.null(input$deFile) && !is.null(input$deFile$name)) {
        file_path_sans_ext(basename(input$deFile$name))
      } else {
        "contrasts"
      }
      
      contrasts_str <- if (!is.null(input$compare_contrasts) && length(input$compare_contrasts) > 0) {
        paste(input$compare_contrasts, collapse = "-")
      } else {
        "contrasts"
      }
      contrasts_str <- gsub("[^A-Za-z0-9._-]+", "__", contrasts_str)
      
      gs_str <- paste0(input$compare_gs_collection, if (nzchar(input$compare_gs_subcollection)) paste0("_", input$compare_gs_subcollection) else "")
      file_name <- paste(file_base, contrasts_str, "GSEA", gs_str, sep = "__")
      
      # Prepare wide format table with NES values
      df <- compare_gsea_data()$nes_df %>%
        mutate(
          padj_fmt = formatC(padj, format = "e", digits = 2),
          NES_display = paste0(round(NES, 2), " (", padj_fmt, ")")
        ) %>%
        select(pathway, contrast, NES_display) %>%
        pivot_wider(
          names_from = contrast,
          values_from = NES_display,
          values_fill = "NS"
        )
      
      datatable(
        df,
        extensions = 'Buttons',
        rownames = FALSE,
        filter = 'top',
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
  
  #==============================
  # Output: Leading Edge UpSet Plot
  #==============================
  output$leadingEdgeUpsetPlot <- renderPlot({
    req(leading_edge_data())
    
    upset(
      leading_edge_data()$upset_df,
      colnames(leading_edge_data()$leading_edge_matrix),
      name = "Leading Edge Genes",
      min_size = 1,
      base_annotations = list(
        'Intersection size' = intersection_size(
          text = list(size = 5)
        )
      ),
      themes = upset_default_themes(
        text = element_text(size = 16),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14)
      )
    )
  })
  
  #==============================
  # Output: Leading Edge Summary
  #==============================
  output$leadingEdgeSummary <- renderText({
    req(leading_edge_data())
    
    stats <- leading_edge_data()$overlap_stats
    
    summary_text <- paste0(
      "Total Leading Edge Genes: ", stats$total_genes, "\n",
      "Genes in ALL Contrasts: ", stats$genes_in_all, "\n\n",
      "Genes Unique to Each Contrast:\n"
    )
    
    unique_counts <- map_chr(names(stats$genes_unique), ~ {
      paste0("  ", .x, ": ", length(stats$genes_unique[[.x]]))
    })
    
    paste0(summary_text, paste(unique_counts, collapse = "\n"))
  })
  
  #==============================
  # Output: Leading Edge Table
  #==============================
  output$leadingEdgeTable <- renderDT(
    {
      req(leading_edge_data())
      
      # Prepare file name
      file_base <- if (!is.null(input$deFile) && !is.null(input$deFile$name)) {
        file_path_sans_ext(basename(input$deFile$name))
      } else {
        "contrasts"
      }
      
      pathway_str <- gsub("[^A-Za-z0-9._-]+", "__", input$selected_pathway_compare)
      file_name <- paste(file_base, pathway_str, "LeadingEdge", sep = "__")
      
      # Create display table with presence markers
      df <- leading_edge_data()$leading_edge_matrix %>%
        as.data.frame() %>%
        rownames_to_column("gene") %>%
        mutate(across(-gene, ~ ifelse(. == 1, "✓", "")))
      
      datatable(
        df,
        extensions = 'Buttons',
        rownames = FALSE,
        filter = 'top',
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
  
}