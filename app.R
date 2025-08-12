library(shiny)
library(tidyverse)
library(viridis)
library(ComplexHeatmap)
library(DESeq2)
library(plotly)
library(DT)
library(RColorBrewer)

ui <- fluidPage(
  titlePanel("DESeq2 Interactive Viewer"),
  
  fluidRow(
    column(4, fileInput("ddsFile", "Upload dds (.rds)", accept = ".rds")),
    column(4, fileInput("deFile", "Upload Contrasts (.rds)", accept = ".rds")),
    column(4, fileInput("varPartFile", "Upload VarPart (.rds) — Optional", accept = ".rds"))
  ),
  
  sidebarLayout(
    sidebarPanel(
      width = 2,
      conditionalPanel(
        condition = "input.mainTab == 'Explore Contrast'",
        sliderInput("fc_cutoff", "Fold Change cutoff", min = 0, max = 10, value = 3, step = 0.5),
        numericInput("top_n", "Top DE genes (heatmap)", value = 25, min = 15, step = 5),
        uiOutput("contrastSelect")
      ),
      conditionalPanel(
        condition = "input.mainTab == 'Explore Gene'",
        textInput("gene_select", "Enter Gene ID", placeholder = "e.g., ENSG00000141510"),
        tags$h4(textOutput("geneSymbol"), style = "margin-top: 10px; margin-bottom: 20px;"),
        selectInput("x_col", "X-axis", choices = NULL),
        selectInput("color_col", "Color", choices = NULL),
        selectInput("shape_col", "Shape", choices = NULL)
      )
    ),
    
    mainPanel(
      tabsetPanel(id = "mainTab",
                  tabPanel(title = "Explore Contrast", value = "Explore Contrast",
                           tabsetPanel(
                             tabPanel("Table", DTOutput("DETable")),
                             tabPanel("Volcano", plotlyOutput("volcanoPlot", height = "600px")),
                             tabPanel("Heatmap", plotOutput("heatmapPlot", height = "700px"))
                           )
                  ),
                  tabPanel(title = "Explore Gene", value = "Explore Gene",
                           uiOutput("geneSubTabs")
                  )
      )
    )
  )
)

server <- function(input, output, session) {
  options(shiny.maxRequestSize = 200 * 1024^2)  
  
  dds_obj <- reactiveVal()
  de_df <- reactiveVal()
  varpart_obj <- reactiveVal()
  
  observeEvent(input$ddsFile, {
    req(input$ddsFile)
    dds_obj(readRDS(input$ddsFile$datapath))
  })
  
  observeEvent(input$deFile, {
    req(input$deFile)
    de_df(readRDS(input$deFile$datapath))
  })
  
  observeEvent(input$varPartFile, {
    req(input$varPartFile)
    varpart_obj(readRDS(input$varPartFile$datapath))
  })
  
  observe({
    req(dds_obj())
    col_vars <- colnames(colData(dds_obj()))
    updateSelectInput(session, "x_col", choices = col_vars)
    updateSelectInput(session, "color_col", choices = col_vars)
    updateSelectInput(session, "shape_col", choices = col_vars)
  })
  
  output$contrastSelect <- renderUI({
    req(de_df())
    contrast_choices <- unique(de_df()$contrast)
    selectInput("contrast", "Select Contrast", choices = contrast_choices)
  })
  
  selected_data <- reactive({
    req(de_df(), input$contrast)
    de_df() %>%
      filter(contrast == input$contrast) %>%
      mutate(
        tooltip = paste0(symbol, " (", Geneid, ")\nlog2FC: ", round(log2FoldChange_shrunk, 2), "\nFDR: ", signif(padj, 3)),
        DE = padj < 0.05 & abs(log2FoldChange_shrunk) > input$fc_cutoff,
        regulated = case_when(
          padj < 0.05 & log2FoldChange_shrunk > input$fc_cutoff  ~ "up",
          padj < 0.05 & log2FoldChange_shrunk < -input$fc_cutoff ~ "down",
          TRUE ~ NA_character_
        )
      )
  })
  
  output$DETable <- renderDT({
    req(selected_data())
    selected_data() %>%
      filter(DE) %>%
      mutate(across(c(log2FoldChange, log2FoldChange_shrunk), ~ round(.x, 2))) %>%
      mutate(padj = formatC(padj, format = "e", digits = 2)) %>%
      arrange(desc(abs(log2FoldChange_shrunk))) %>%
      select(Geneid, symbol, biotype, log2FoldChange, log2FoldChange_shrunk, padj, sign, DE, regulated)
  })
  
  output$volcanoPlot <- renderPlotly({
    req(selected_data())
    gg <- ggplot(selected_data(), aes(x = log2FoldChange_shrunk, y = -log10(padj), text = tooltip)) +
      geom_point(aes(color = DE), alpha = 0.6) +
      scale_color_manual(values = c("TRUE" = "red", "FALSE" = "gray"), guide = "none") +
      geom_vline(xintercept = c(-input$fc_cutoff, input$fc_cutoff), linetype = "dashed") +
      geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
      theme_minimal()
    
    ggplotly(gg, tooltip = "text")
  })
  
  output$heatmapPlot <- renderPlot({
    req(selected_data(), dds_obj())
    
    top_genes <- selected_data() %>%
      filter(DE) %>%
      slice_max(order_by = abs(log2FoldChange_shrunk), n = input$top_n)
    
    vsd_mat <- assay(dds_obj(), "vst")[rownames(dds_obj()) %in% top_genes$Geneid, ]
    vsd_mat <- vsd_mat[match(top_genes$Geneid, rownames(vsd_mat)), ]
    
    validate(
      need(nrow(vsd_mat) > 1, "No DE genes pass filter or match rownames in dds")
    )
    
    rownames(vsd_mat) <- top_genes$symbol
    
    Heatmap(
      scale(vsd_mat),
      col = viridis(100),
      column_names_gp = gpar(fontsize = 12, rot = 45),  
      row_names_gp = gpar(fontsize = 12),
      heatmap_legend_param = list(title = "Z-scores")
    )
  })
  
  output$geneSymbol <- renderText({
    req(input$gene_select, de_df())
    
    gene_symbol <- de_df() %>%
      filter(Geneid == input$gene_select) %>%
      pull(symbol) %>%
      unique()
    
    if (length(gene_symbol) > 0 && !is.na(gene_symbol)) {
      paste("Symbol:", gene_symbol)
    } else {
      "Symbol: not found"
    }
  })
  
  output$geneSubTabs <- renderUI({
    tabs <- list(
      tabPanel("Expression Plot", plotOutput("genePlot", height = "600px", width = "1000px")),
      tabPanel("Gene Table", DTOutput("geneDetails"))
    )
    
    if (!is.null(varpart_obj())) {
      tabs <- append(
        tabs,
        list(tabPanel("Variance Decomposition", plotOutput("varPartPlot", height = "600px", width = "1000px")))
      )
    }
    
    do.call(tabsetPanel, tabs)
  })
  
  output$genePlot <- renderPlot({
    req(input$gene_select, input$x_col, dds_obj())
    
    vst_mat <- assay(dds_obj(), "vst")
    
    validate(
      need(input$gene_select %in% rownames(vst_mat), "Gene not found in dataset")
    )
    
    expr_values <- vst_mat[input$gene_select, ]
    meta <- as.data.frame(colData(dds_obj()))
    meta$expression <- expr_values
    
    n_colors <- length(unique(meta[[input$color_col]]))
    palette_colors <- colorRampPalette(brewer.pal(8, "Dark2"))(n_colors)
    
    ggplot(meta, aes_string(x = input$x_col, y = "expression")) +
      geom_boxplot(aes_string(color = input$color_col), outliers = FALSE, show.legend = FALSE) +
      geom_jitter(aes_string(color = input$color_col, shape = input$shape_col), width = 0.2, size = 3, alpha = 0.9) +
      scale_color_manual(values = palette_colors) +
      labs(title = NULL, y = "VST expression", x = input$x_col) +
      theme_bw(base_size = 20) +
      theme(
        axis.text = element_text(angle = 45, hjust = 1),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()
      )
  })
  
  output$geneDetails <- renderDT({
    req(input$gene_select, de_df())
    
    gene_table <- de_df() %>%
      filter(Geneid == input$gene_select) %>%
      mutate(across(c(log2FoldChange, log2FoldChange_shrunk), ~ round(.x, 2))) %>%
      mutate(padj = formatC(padj, format = "e", digits = 2)) %>%
      arrange(desc(abs(log2FoldChange))) %>%
      select(contrast, log2FoldChange, log2FoldChange_shrunk, padj, sign, DE, regulated)
    
    datatable(
      gene_table,
      options = list(
        dom = 't',        
        ordering = TRUE,
        pageLength = nrow(gene_table)
      ),
      rownames = FALSE
    ) %>%
      formatStyle(
        'regulated',
        target = 'row',
        backgroundColor = styleEqual(
          c("up", "down"),
          c("#d0f0c0", "#f4cccc")
        )
      )
  })
  
  output$varPartPlot <- renderPlot({
    req(input$gene_select, varpart_obj())
    
    gene_id <- input$gene_select
    result <- varpart_obj()
    
    vp_gene <- result$varPart[gene_id, ]
    
    vp_df <- data.frame(
      Factor = names(vp_gene),
      Variance = as.numeric(vp_gene)
    )
    
    ggplot(vp_df, aes(x = reorder(Factor, -Variance), y = Variance)) +
      geom_col(fill = "steelblue") +
      scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
      labs(title = NULL, x = NULL, y = "Fraction of Variance") +
      theme_bw(base_size = 20) +
      theme(
        axis.text = element_text(angle = 45, hjust = 1),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()
      )
  })

}
shinyApp(ui, server)

