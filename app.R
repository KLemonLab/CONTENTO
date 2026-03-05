#== == == == == == == == == == == == == == == == ==
#===== DESeq2 Explorer App  =======================
#== == == == == == == == == == == == == == == == ==

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
library(fgsea)
library(msigdbr)
library(pheatmap)
library(circlize)
library(SummarizedExperiment)


# Source helper server modules
source("R/load_files_server.R")
source("R/explore_contrast_server.R")
source("R/compare_contrast_server.R")
source("R/explore_gene_server.R")

#== == == == == == == == == == == == == == == == ==
#===== UI =========================================
#== == == == == == == == == == == == == == == == ==

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
    h4("Upload Files", style = "padding-left: 20px;"),
    fileInput("seFile", "Summarized Experiment",
              accept = ".rds",
              placeholder = "SE object (.rds)"),
    uiOutput("annotationStatus"),
    hr(),
    h4("DE Cutoffs", style = "padding-left: 20px;"),
    sliderInput("global_log2FC_cutoff", "log2FC cutoff", min = 0, max = 8, value = 2, step = 0.5),
    numericInput("global_padj_cutoff", "p-value FDR cutoff", value = 0.05, min = 0, max = 1, step = 0.01)
  ),
  
  dashboardBody(
    
    tabItems(
      
      ## ---- Getting Started -----
      tabItem(tabName = "intro",
              
              # First box: Required Inputs
              fluidRow(
                box(title = "Welcome to the KLemon Lab RNASeq Explorer!", width = 12, status = "info", solidHeader = TRUE,
                    tags$h4("Required Inputs (Upload on the left sidebar):"),
                    tags$ul(
                      tags$li(strong("SE file:"), " SummarizedExperiment object containing counts, contrasts, and optional variance partition data."),
                      tags$li(strong("Annotation:"), " Automatically loaded from your annotations folder to match the metadata in your SE object. Upload if needed.")
                    )
                )
              ),
              
              # Second row: three main app sections
              fluidRow(
                box(title = tagList(icon("layer-group"), "Explore by Contrast"), width = 12, status = "success", solidHeader = TRUE,
                    p("Select a contrast to explore its results and adjust fold-change cutoff for both table and plots. Controls also allow adjusting number of top genes and colors for the heat as well as gene set collections for enrichment analyses."),
                    tags$ul(
                      tags$li(strong("Selected Genes:"), " Table of differential expression results. You can save a .csv file with the filtered results."),
                      tags$li(strong("Volcano Plot:"), " Interactive plot of significance vs. fold-change. The log2 fold-changes are shrinked to stabilize variance across genes."),
                      tags$li(strong("Heatmap:"), " Visualizes top DE genes across samples. Expression values are normalized using a variance-stabilizing transformation (VST) and scaled to Z-scores for visualization."),
                      tags$li(strong("GSEA:"), " Gene Set Enrichment Analysis for human datasets. Ranks genes by their DE statistics and identifies enriched pathways using data sets from MSigDB."),
                      tags$li(strong("GESECA:"), " Gene Set Co-expression Analysis for human datasets. Identifies co-expressed gene sets based on VST-normalized expression values.")
                    )
                )
              ),
              fluidRow(
                box(title = tagList(icon("exchange-alt"), "Compare Contrast"), width = 12, status = "warning", solidHeader = TRUE,
                    p("Compare multiple contrasts simultaneously. Use contrast selector to control displayed results."),
                    tags$ul(
                      tags$li(strong("DEG Overlap:"), " Shows overlap of significant genes between selected contrasts using the fold-change cutoff."),
                      tags$li(strong("Gene Sets Overlap:"), " Shows overlap of enriched gene sets between selected contrasts (only for human datasets).")
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
      
      ## ---- Explore by Contrast -----
      tabItem(tabName = "contrast",
              fluidRow(
                box(title = "Controls", width = 3, status = "info", collapsible = TRUE,
                    uiOutput("contrastSelect"),
                    numericInput("top_n", "Top DE genes (Heatmap)", value = 25, min = 15, step = 5),
                    selectInput("viridis_palette", "Color palette (Heatmap)",
                                choices = c("viridis", "magma", "plasma", "inferno", "cividis"), selected = "viridis"),
                    selectInput("gs_collection", "Gene Set (GSEA/GESECA)",
                                choices = c("H", "C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8"), selected = "H"),
                    conditionalPanel(condition = "!(input.gs_collection == 'H' || input.gs_collection == 'C1' || input.gs_collection == 'C6' || input.gs_collection == 'C8')",
                                     textInput("gs_subcollection", "Subcollection (optional)", placeholder = "e.g., CP:REACTOME")),
                    tags$div(style = "margin-top: -10px; margin-bottom: 15px;",
                             tags$a(href = "https://www.gsea-msigdb.org/gsea/msigdb/human/annotate.jsp",
                                    target = "_blank",
                                    icon("external-link-alt"),
                                    "MSigDB"))
                ),
                
                box(title = "Explore Results by Contrast", width = 9, status = "primary", collapsible = TRUE,
                    uiOutput("contrastSubTabs"),
                    hr(),
                    tags$div(style = "text-align: center;",
                             downloadButton("downloadDETableFull",
                                            "Download Filtered DE Genes (Full Annotations)",
                                            class = "btn-primary")
                    )
                )
              )
      ),
      
      ## ---- Compare Contrast -----
      tabItem(tabName = "compare",
              fluidRow(
                box(title = "Controls", width = 3, status = "info", collapsible = TRUE,
                    uiOutput("multiContrastSelect"),
                    selectInput("compare_gs_collection", "Gene Set (GSEA/GESECA)",
                                choices = c("H", "C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8"), selected = "H"),
                    conditionalPanel(condition = "!(input.compare_gs_collection == 'H' || input.compare_gs_collection == 'C1' || input.compare_gs_collection == 'C6' || input.compare_gs_collection == 'C8')",
                                     textInput("compare_gs_subcollection", "Subcollection (optional)", placeholder = "e.g., CP:REACTOME")),
                    tags$div(style = "margin-top: -10px; margin-bottom: 15px;",
                             tags$a(href = "https://www.gsea-msigdb.org/gsea/msigdb/human/annotate.jsp",
                                    target = "_blank",
                                    icon("external-link-alt"),
                                    "MSigDB")),
                    numericInput("max_pathways", "Max pathways to show", value = 100, min = 10, max = 500, step = 10),
                    numericInput("pathway_name_length", "Max pathway name length", value = 50, min = 25, max = 500, step = 10)
                ),
                
                box(title = "Compare Contrast Results", width = 9, status = "primary", collapsible = TRUE,
                    uiOutput("compareSubTabs"),
                    hr(),
                    tags$div(style = "text-align: center;",
                             downloadButton("downloadCompareTableFull",
                                            "Download Filtered DEGs (Full Annotations)",
                                            class = "btn-primary")
                    )
                )
              )
              
      ),
      
      ## ---- Explore by Gene -----
      tabItem(tabName = "gene",
              fluidRow(
                box(title = "Controls", width = 3, status = "info", collapsible = TRUE,
                    textInput("gene_select", "Enter Gene ID", placeholder = "e.g., ENSG00000141510"),
                    selectInput("x_col", "X-axis", choices = NULL),
                    selectInput("color_col", "Color", choices = NULL),
                    selectInput("shape_col", "Shape", choices = NULL),
                    uiOutput("contrastSelectGene"),
                    numericInput("neigh_window", "Neighbourhood window (nt)", value = 10000, step = 100, min = 0)
                ),
                
                box(title = "Explore Results by Gene", width = 9, status = "primary", collapsible = TRUE,
                    uiOutput("geneSubTabs")
                )
              )
      )
    )
  )
)

#== == == == == == == == == == == == == == == == ==
#===== Server =====================================
#== == == == == == == == == == == == == == == == ==

server <- function(input, output, session) {
  options(shiny.maxRequestSize = 200 * 1024^2)
  
  # Shared state across helpers
  state <- list(
    se_obj = reactiveVal(NULL),
    de_df = reactiveVal(NULL),
    varpart_obj = reactiveVal(NULL),
    annotation_df = reactiveVal(NULL),
    se_organism = reactiveVal(NULL),
    filtered_de_df = reactiveVal(NULL)
  )
  
  # Reactive organism type from SE metadata
  organism <- reactive({
    req(state$se_obj())
    meta_organism <- tryCatch(
      metadata(state$se_obj())$organism,
      error = function(e) NULL
    )
    if (!is.null(meta_organism) && nzchar(trimws(meta_organism))) {
      meta_organism
    } else {
      NULL
    }
  })
  
  # Wire helpers
  load_files_server(input, output, session, state)
  explore_contrast_server(input, output, session, state, organism)
  compare_contrast_server(input, output, session, state, organism)
  explore_gene_server(input, output, session, state, organism)
}

#== == == == == == == == == == == == == == == == ==
#===== Run App ====================================
#== == == == == == == == == == == == == == == == ==

shinyApp(ui, server)
