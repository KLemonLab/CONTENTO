library(shiny)
library(tidyverse)
library(purrr)
library(viridis)
library(ComplexHeatmap)
library(DESeq2)
library(plotly)
library(DT)
library(RColorBrewer)
library(shinycssloaders)
library(ComplexUpset) 
library(scales)

# Source helpers
source("R/load_files_server.R")
source("R/explore_contrast_server.R")
source("R/compare_contrast_server.R")
source("R/explore_gene_server.R")

# ==========================
# ========   UI   ==========
# ==========================

ui <- fluidPage(
  titlePanel("DESeq2 Interactive Viewer"),
  
  fluidRow(
    column(4, fileInput("ddsFile", "Upload dds (.rds)", accept = ".rds")),
    column(4, fileInput("deFile", "Upload Contrasts (.rds)", accept = ".rds")),
    column(4, fileInput("varPartFile", "Upload VarPart (.rds) — Optional", accept = ".rds"))
  ),
  
  sidebarLayout(
    sidebarPanel(
      width = 3,
      
      # --- Explore by Contrast Sidebar ---
      conditionalPanel(
        condition = "input.mainTab == 'Explore by Contrast'",
        sliderInput("fc_cutoff", "Fold Change cutoff", min = 0, max = 8, value = 3, step = 0.5),
        numericInput("top_n", "Top DE genes (heatmap)", value = 25, min = 15, step = 5),
        selectInput("viridis_palette", "Color palette (heatmap)",
                    choices = c("viridis", "magma", "plasma", "inferno", "cividis"),
                    selected = "viridis"),
        uiOutput("contrastSelect")
      ),
      
      # --- Compare Contrast Sidebar ---
      conditionalPanel(
        condition = "input.mainTab == 'Compare Contrast'",
        sliderInput("compare_fc_cutoff", "Fold Change cutoff", min = 0, max = 8, value = 3, step = 0.5),
        uiOutput("multiContrastSelect")
      ),
      
      # --- Explore by Gene Sidebar ---
      conditionalPanel(
        condition = "input.mainTab == 'Explore by Gene'",
        textInput("gene_select", "Enter Gene ID", placeholder = "e.g., ENSG00000141510"),
        tags$h4(textOutput("geneSymbol"), style = "margin-top: 10px; margin-bottom: 20px;"),
        selectInput("x_col", "X-axis", choices = NULL),
        selectInput("color_col", "Color", choices = NULL),
        selectInput("shape_col", "Shape", choices = NULL)
      )
    ),
    
    mainPanel(
      tabsetPanel(id = "mainTab",
                  
                  # --- Explore Contrast Tab ---
                  tabPanel(title = "Explore by Contrast", value = "Explore by Contrast",
                           tabsetPanel(
                             tabPanel("Selected Genes", withSpinner(DTOutput("DETable"), type = 5)),
                             tabPanel("Volcano Plot", withSpinner(plotlyOutput("volcanoPlot", height = "500px"), type = 5)),
                             tabPanel("Heatmap", withSpinner(plotOutput("heatmapPlot", height = "700px"), type = 5))
                           )
                  ),
                  
                  # --- Compare Contrast Tab ---
                  tabPanel(title = "Compare Contrast", value = "Compare Contrast",
                           tabsetPanel(
                             tabPanel("Upset Plot", withSpinner(plotOutput("compareUpsetPlot", height = "500px"), type = 5)),
                             tabPanel("Selected Contrasts", withSpinner(DTOutput("compareTable"), type = 5))
                           )
                  ),
                  
                  # --- Explore Gene Tab ---
                  tabPanel(title = "Explore by Gene", value = "Explore by Gene",
                           uiOutput("geneSubTabs")
                  )
      )
    )
  )
)

# ==========================
# ======  SERVER  ==========
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

shinyApp(ui, server)
