# ==========================
# DESeq2 Explorer Dashboard
# ==========================

library(shiny)
library(shinydashboard)
library(shinycssloaders)
library(DT)
library(plotly)
library(tidyverse)
library(purrr)
library(viridis)
library(ComplexHeatmap)
library(DESeq2)
library(plotly)
library(DT)
library(RColorBrewer)
library(ComplexUpset) 
library(scales)

# Source helper server modules
source("R/load_files_server.R")
source("R/explore_contrast_server.R")
source("R/compare_contrast_server.R")
source("R/explore_gene_server.R")

# ==========================
# UI
# ==========================

ui <- dashboardPage(
  dashboardHeader(title = "RNASeq Explorer"),
  
  dashboardSidebar(
    collapsed = FALSE,
    sidebarMenu(
      id = "tabs",
      menuItem("Explore by Contrast", tabName = "contrast", icon = icon("chart-line")),
      menuItem("Compare Contrast", tabName = "compare", icon = icon("exchange-alt")),
      menuItem("Explore by Gene", tabName = "gene", icon = icon("dna"))
    ),
    hr(),
    h4("Upload Files", style = "padding-left: 20px;"),
    fileInput("ddsFile", "DESeq2", accept = ".rds", placeholder = "dds object"),
    fileInput("deFile", "Contrasts", accept = ".rds", placeholder = "Contrast object"),
    fileInput("varPartFile", "VarPart", accept = ".rds", placeholder = "Optional VarPart")
  ),
  
  dashboardBody(
    
    tabItems(
      
      # --- Explore by Contrast ---
      tabItem(tabName = "contrast",
              fluidRow(
                box(title = "Controls", width = 3, status = "info", collapsible = TRUE,
                    sliderInput("fc_cutoff", "Fold Change cutoff", min = 0, max = 8, value = 3, step = 0.5),
                    numericInput("top_n", "Top DE genes (heatmap)", value = 25, min = 15, step = 5),
                    selectInput("viridis_palette", "Color palette (heatmap)",
                                choices = c("viridis", "magma", "plasma", "inferno", "cividis"), selected = "viridis"),
                    uiOutput("contrastSelect")
                ),
                
                box(title = "Explore Results by Contrast", width = 9, status = "primary", collapsible = TRUE,
                    tabsetPanel(
                      tabPanel("Selected Genes", withSpinner(DTOutput("DETable"), type = 5)),
                      tabPanel("Volcano Plot", withSpinner(plotlyOutput("volcanoPlot", height = "500px"), type = 5)),
                      tabPanel("Heatmap", withSpinner(plotOutput("heatmapPlot", height = "700px"), type = 5))
                    )
                )
              )
      ),
      
      # --- Compare Contrast ---
      tabItem(tabName = "compare",
              fluidRow(
                box(title = "Controls", width = 3, status = "info", collapsible = TRUE,
                    sliderInput("compare_fc_cutoff", "Fold Change cutoff", min = 0, max = 8, value = 3, step = 0.5),
                    uiOutput("multiContrastSelect")
                ),

                box(title = "Compare Contrast Results", width = 9, status = "primary", collapsible = TRUE,
                    tabsetPanel(
                      tabPanel("Upset Plot", withSpinner(plotOutput("compareUpsetPlot", height = "500px"), type = 5)),
                      tabPanel("Selected Contrasts", withSpinner(DTOutput("compareTable"), type = 5))
                    )
                )
              )
      ),
      
      # --- Explore by Gene ---
      tabItem(tabName = "gene",
              fluidRow(
                box(title = "Controls", width = 3, status = "info", collapsible = TRUE,
                    textInput("gene_select", "Enter Gene ID", placeholder = "e.g., ENSG00000141510"),
                    tags$h4(textOutput("geneSymbol"), style = "margin-top: 10px; margin-bottom: 20px;"),
                    selectInput("x_col", "X-axis", choices = NULL),
                    selectInput("color_col", "Color", choices = NULL),
                    selectInput("shape_col", "Shape", choices = NULL)
                ),
                
                box(title = "Explore Results by Gene", width = 9, status = "primary", collapsible = TRUE,
                    uiOutput("geneSubTabs")
                )
              )
      )
    )
  )
)

# ==========================
# Server
# ==========================

server <- function(input, output, session) {
  options(shiny.maxRequestSize = 200 * 1024^2)
  
  # Shared state across helpers
  state <- list(
    dds_obj = reactiveVal(NULL),
    de_df = reactiveVal(NULL),
    varpart_obj = reactiveVal(NULL)
  )
  
  # Wire helpers
  load_files_server(input, output, session, state)
  explore_contrast_server(input, output, session, state)
  compare_contrast_server(input, output, session, state)
  explore_gene_server(input, output, session, state)
}

# ==========================
# Run App
# ==========================

shinyApp(ui, server)
