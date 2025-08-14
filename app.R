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
library(tools)
library(grid)
library(ggiraph)

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
      menuItem("Getting Started", tabName = "intro", icon = icon("info-circle")),
      menuItem("Explore by Contrast", tabName = "contrast", icon = icon("chart-line")),
      menuItem("Compare Contrast", tabName = "compare", icon = icon("exchange-alt")),
      menuItem("Explore by Gene", tabName = "gene", icon = icon("dna"))
    ),
    hr(),
    h4("Select Organism", style = "padding-left: 20px;"),
    selectInput("organism", "Organism type:", choices = c("Human", "Bacteria"), selected = "Bacteria"),
    hr(),
    h4("Upload Files", style = "padding-left: 20px;"),
    fileInput("ddsFile", "DESeq2", accept = ".rds", placeholder = "dds object"),
    fileInput("deFile", "Contrasts", accept = ".rds", placeholder = "Contrast object"),
    fileInput("varPartFile", "VarPart", accept = ".rds", placeholder = "Optional VarPart")
  ),
  
  dashboardBody(
    
    tabItems(
      
      # --- Getting Started ---
      tabItem(tabName = "intro",
              
              # --- First box: Required Inputs ---
              fluidRow(
                box(title = "Welcome to the KLemon Lab RNASeq Explorer!", width = 12, status = "info", solidHeader = TRUE,
                    tags$h4("Required Inputs (Upload .rds files on the left sidebar):"),
                    tags$ul(
                      tags$li(strong("Organism type :"), " Select based on the RNASeq reads analyzed by DESeq2."),
                      tags$li(strong("DESeq2 file:"), " Contains the DESeqDataSet object with count data and metadata."),
                      tags$li(strong("Contrasts file:"), " Contains results for one or more contrasts from the DESeq2 analysis."),
                      tags$li(strong("VarPart file:"), " Optional — contains variance partition analysis results.")
                    )
                )
              ),
              
              # --- Second row: three main app sections ---
              fluidRow(
                box(title = tagList(icon("chart-line"), "Explore by Contrast"), width = 12, status = "success", solidHeader = TRUE,
                    p("Select a contrast to explore its results and adjust fold-change cutoff for both table and plots. Controls also allow adjusting number of top genes and colors for the heatmap."),
                    tags$ul(
                      tags$li(strong("Selected Genes:"), " Table of differential expression results. You can save a .csv file with the filtered results."),
                      tags$li(strong("Volcano Plot:"), " Interactive plot of significance vs. fold-change. The log2 fold-changes are shrinked to stabilize variance across genes."),
                      tags$li(strong("Heatmap:"), " Visualizes top DE genes across samples. Expression values are normalized using a variance-stabilizing transformation (VST) and scaled to Z-scores for visualization.")
                    )
                )
              ),
              fluidRow(
                box(title = tagList(icon("exchange-alt"), "Compare Contrast"), width = 12, status = "warning", solidHeader = TRUE,
                    p("Compare multiple contrasts simultaneously. Use the fold-change cutoff and contrast selector to control displayed results."),
                    tags$ul(
                      tags$li(strong("Upset Plot:"), " Shows overlap of significant genes between selected contrasts."),
                      tags$li(strong("Selected Contrasts:"), " Table of results filtered by your settings and colored by up/down regulation. You can save a .csv file.")
                    )
                )
              ),
              fluidRow(
                box(title = tagList(icon("dna"), "Explore by Gene"), width = 12, status = "danger", solidHeader = TRUE,
                    p("Visualize expression for a single gene. Use a Gene ID (from your DESeq2 contrasts analysis) to select the gene you want to explore."),
                    tags$ul(
                      tags$li(strong("Expression Plot:"), " Shows the VST-normalized expression for the selected gene. You can select the variables used for the X-axis, color and shapes from the metadata columns in your DESeq2 object."),
                      tags$li(strong("Gene Table:"), " Displays DE results for the selected gene across all contrasts, with rows colored by up/down regulation."),
                      tags$li(strong("Variance Decomposition:"), " If a variance partition object is loaded, this plot shows the fraction of expression variance explained by different factors for the selected gene."),
                      tags$li(strong("Neighbourhood Analysis:"), " Allows exploration of the genes up/downstreams the selected one in bacterial datasets.")
                    )
                )
              )
      ),
      
      # --- Explore by Contrast ---
      tabItem(tabName = "contrast",
              fluidRow(
                box(title = "Controls", width = 3, status = "info", collapsible = TRUE,
                    uiOutput("contrastSelect"),
                    sliderInput("fc_cutoff", "log2FC cutoff", min = 0, max = 8, value = 3, step = 0.5),
                    numericInput("top_n", "Top DE genes (heatmap)", value = 25, min = 15, step = 5),
                    selectInput("viridis_palette", "Color palette (heatmap)",
                                choices = c("viridis", "magma", "plasma", "inferno", "cividis"), selected = "viridis")
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
                    sliderInput("compare_fc_cutoff", "log2FC cutoff", min = 0, max = 8, value = 3, step = 0.5),
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
                    textInput("gene_select", "Enter Gene ID", placeholder = "e.g., ENSG00000141510", value = "HMPREF1290_RS08565"),
                    tags$h4(textOutput("geneSymbol"), style = "margin-top: 10px; margin-bottom: 20px;"),
                    selectInput("x_col", "X-axis", choices = NULL),
                    selectInput("color_col", "Color", choices = NULL),
                    selectInput("shape_col", "Shape", choices = NULL),
                    sliderInput("gene_fc_cutoff", "log2FC cutoff", min = 0, max = 8, value = 3, step = 0.5),
                    uiOutput("contrastSelectGene"),
                    numericInput("neigh_window", "Neighbourhood window (nt)", value = 1000, step = 100, min = 0)
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
