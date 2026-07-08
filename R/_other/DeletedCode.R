# == == == == == == == == == == == == == == == == == == == == == == == == ==
#### DELETED CODE ####
# == == == == == == == == == == == == == == == == == == == == == == == == ==

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##### UI: Pathway Selection for GSEA Detail Tab #####
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
output$pathwaySelectUI_gsea <- renderUI({
  req(gsea_result())
  
  pathways <- gsea_result()$fgseaRes %>%
    filter(padj < 0.05) %>%
    arrange(padj) %>%
    pull(pathway)
  
  if (length(pathways) > 0) {
    selectInput("selected_pathway_gsea", "Select Pathway:", 
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
##### Output: GSEA Enrichment Plot #####
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
output$enrichmentPlot <- renderPlot({
  req(gsea_result(), input$selected_pathway_gsea)
  
  tryCatch({
    # Get pathway genes
    pathway_genes <- gsea_result()$pathways_list[[input$selected_pathway_gsea]]
    
    if (is.null(pathway_genes) || length(pathway_genes) == 0) {
      plot.new()
      text(0.5, 0.5, "Pathway not found", cex = 1.2)
      return()
    }
    
    # Create enrichment plot
    plotEnrichment(pathway_genes, gsea_result()$ranks_vec) +
      labs(title = input$selected_pathway_gsea) +
      theme_minimal() +
      theme(plot.title = element_text(size = 10))
    
  }, error = handle_enrichment_plot_error)
})