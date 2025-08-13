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
  # Reactive: Data for upset plot (binary wide format)
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
  # Output: Upset plot
  #==============================
  output$compareUpsetPlot <- renderPlot({
    req(compare_data())
    upset(
      compare_data(),
      input$compare_contrasts,
      name = "DEGs",
      min_size = 1,
      base_annotations = list(
        'Intersection size' = intersection_size(text = list(size = 3))
      )
    )
  })
  
  #==============================
  # Output: Table of DEGs
  #==============================
  library(circlize)
  
  output$compareTable <- renderDT({
    req(compare_table_data())
    df <- compare_table_data()
    lfc_cols <- setdiff(colnames(df), c("Geneid", "symbol"))
    
    # Determine min and max
    lfc_min <- min(df[lfc_cols], na.rm = TRUE)
    lfc_max <- max(df[lfc_cols], na.rm = TRUE)
    
    datatable(df, options = list(pageLength = 20, scrollX = TRUE), rownames = FALSE) %>%
      formatStyle(
        columns = lfc_cols,
        backgroundColor = styleInterval(
          0,
          c("lightblue", "pink")
        )
      )
  })
}