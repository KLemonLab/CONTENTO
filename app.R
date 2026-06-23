# == == == == == == == == == == == == == == == == == == == == == == == == ==
# ========================= CONTENTO App  ===========================
# == == == == == == == == == == == == == == == == == == == == == == == == ==

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


# Source utilities and server modules
source("R/utils.R")
source("R/utils-ui.R")
source("R/utils-global.R")
source("R/load_files_server.R")
source("R/explore_contrast_server.R")
source("R/compare_contrast_server.R")
source("R/explore_gene_server.R")

# == == == == == == == == == == == == == == == == == == == == == == == == ==
# ============================= UI =========================================
# == == == == == == == == == == == == == == == == == == == == == == == == ==

ui <- dashboardPage(
  dashboardHeader(title = "RNASeq Explorer"),
  
  dashboardSidebar(
    collapsed = FALSE,
    sidebarMenu(
      id = "tabs",
      menuItem("Getting Started",     tabName = "intro",    icon = icon("info-circle")),
      menuItem("Explore by Contrast", tabName = "contrast", icon = icon("chart-line")),
      menuItem("Compare Contrast",    tabName = "compare",  icon = icon("exchange-alt")),
      menuItem("Explore by Gene",     tabName = "gene",     icon = icon("dna"))
    ),
    hr(),
    h4("Upload Files", style = "padding-left: 20px;"),
    fileInput("seFile", "Summarized Experiment",
              accept = ".rds",
              placeholder = "SE object (.rds)"),
    uiOutput("annotationStatus"),
    hr(),
    h4("DE Cutoffs", style = "padding-left: 20px;"),
    selectInput("global_lfc_col", "fold-change column",
                choices = c("Shrunken (log2FC_shrunk)" = "log2FC_shrunk",
                            "Unshrunken (log2FC)"      = "log2FC"),
                selected = "log2FC_shrunk"),
    sliderInput("global_log2FC_cutoff", "fold-change cutoff",
                min = 0, max = 8, value = 2, step = 0.5),
    numericInput("global_padj_cutoff", "p-value FDR cutoff",
                 value = 0.05, min = 0, max = 1, step = 0.01)
  ),
  
  dashboardBody(
    
    includeCSS("www/styles.css"),
    
    tabItems(
      
      ## ---- Getting Started -----
      tabItem(tabName = "intro",
              
              fluidRow(
                box(title = "Welcome to the KLemon Lab RNASeq Explorer!", width = 12,
                    status = "info", solidHeader = TRUE,
                    tags$h4("Required Inputs (Upload on the left sidebar):"),
                    tags$ul(
                      tags$li(strong("SE file:"), " SummarizedExperiment object containing counts, contrasts, and optional variance partition data."),
                      tags$li(strong("Annotation:"), " Automatically loaded from your annotations folder to match the metadata in your SE object. Upload if needed."),
                      tags$li(strong("DE Cutoffs:"), " Select fold-change column and adjust cutoffs for fold-change and p-value to define significance.")
                    )
                )
              ),
              
              fluidRow(
                box(title = tagList(icon("chart-line"), "Explore by Contrast"), width = 12,
                    status = "success", solidHeader = TRUE,
                    p("Select a contrast to explore its results."),
                    tags$ul(
                      tags$li(strong("Selected Genes:"), " Table of differential expression results."),
                      tags$li(strong("Volcano Plot:"), " Interactive plot of significance vs. fold-change."),
                      tags$li(strong("GSEA:"), " Gene Set Enrichment Analysis ranked by DE statistics.")
                    )
                )
              ),
              fluidRow(
                box(title = tagList(icon("exchange-alt"), "Compare Contrast"), width = 12,
                    status = "warning", solidHeader = TRUE,
                    p("Compare multiple contrasts simultaneously."),
                    tags$ul(
                      tags$li(strong("Heatmap:"), " Visualizes top DE genes across all samples using VST Z-scores."),
                      tags$li(strong("DEG Overlap:"), " Shows overlap of significant genes between selected contrasts."),
                      tags$li(strong("Gene Sets Overlap:"), " Shows overlap of enriched gene sets between selected contrasts (GSEA)."),
                      tags$li(strong("GESECA:"), " Gene Set Co-expression Analysis across all samples.")
                    )
                )
              ),
              fluidRow(
                box(title = tagList(icon("dna"), "Explore by Gene"), width = 12,
                    status = "danger", solidHeader = TRUE,
                    p("Visualize expression for a single gene across conditions."),
                    tags$ul(
                      tags$li(strong("Gene Info:"), " Displays gene annotations from your annotation file."),
                      tags$li(strong("Expression Plot:"), " VST-normalized expression for the selected gene."),
                      tags$li(strong("Gene Table:"), " DE results across all contrasts, colored by regulation."),
                      tags$li(strong("Variance Decomposition:"), " Fraction of variance explained by experimental factors."),
                      tags$li(strong("Neighbourhood Analysis:"), " Genes upstream/downstream in bacterial datasets.")
                    )
                )
              )
      ),
      
      
      ## ---- Explore by Contrast -----
      tabItem(tabName = "contrast",
              fluidRow(
                box(width = 12, status = "info", solidHeader = TRUE,
                    collapsible = TRUE, collapsed = FALSE,
                    title = tagList(icon("chart-line"), "Explore by Contrast"),
                    fluidRow(
                      column(width = 7,
                             h5(strong("Choose a contrast to explore differential expression results"))),
                      column(width = 5,
                             div(style = "padding-left:10px;", uiOutput("contrastSelect")))
                    )
                )
              ),
              fluidRow(
                box(width = 12, status = "info", solidHeader = TRUE,
                    uiOutput("contrastSubTabs"))
              )
      ),
      
      ## ---- Compare Contrast -----
      tabItem(tabName = "compare",
              fluidRow(
                box(width = 12, status = "info", solidHeader = TRUE,
                    collapsible = TRUE, collapsed = FALSE,
                    title = tagList(icon("exchange-alt"), "Compare Contrast"),
                    fluidRow(
                      column(width = 4,
                             div(style = "padding-right: 20px;",
                                 span("Select multiple contrasts to compare results side by side."),
                                 tags$ul(
                                   style = "margin: 5px 0 12px 15px; padding:0;",
                                   tags$li("Adjust log2FC column, fold-change and p-value cutoffs in the sidebar"),
                                   tags$li("Explore heatmaps, DEG overlap, and cross-contrast enrichment (GSEA / GESECA)")
                                 )
                             )
                      ),
                      column(width = 8,
                             div(style = "padding-left: 10px;",
                                 uiOutput("multiContrastSelectHeader"),
                                 div(style = "max-height: 160px; overflow-y: auto; border: 1px solid #ddd; border-radius: 4px; padding: 10px 12px; margin-top: 4px;",
                                     uiOutput("multiContrastCheckboxes"))
                             )
                      )
                    )
                )
              ),
              fluidRow(
                box(width = 12, status = "info", solidHeader = TRUE,
                    uiOutput("compareSubTabs"), hr())
              )
      ),
      
      ## ---- Explore by Gene -----
      tabItem(tabName = "gene",
              fluidRow(
                box(width = 12, status = "info", solidHeader = TRUE,
                    collapsible = TRUE, collapsed = FALSE,
                    title = tagList(icon("dna"), "Explore by Gene"),
                    fluidRow(
                      column(width = 8,
                             div(
                               span("Choose a gene to explore expression and annotations."),
                               tags$ul(
                                 style = "margin: 5px 0 0 15px; padding:0;",
                                 tags$li("Adjust plotting variables to customize visualization of gene expression across samples"),
                                 tags$li("Review gene-level statistics across contrasts, including DE results and variance partitioning (if available)"),
                                 tags$li("Examine genomic neighbourhood and gene context in bacterial datasets")
                               )
                             )
                      ),
                      column(width = 4,
                             div(style = "padding-left:10px;",
                                 selectizeInput("gene_select",
                                                "Select ID (Search by GeneID or Symbol selected column):",
                                                choices = NULL,
                                                options = list(
                                                  placeholder = 'Start typing gene name or ID...',
                                                  maxOptions = 20,
                                                  loadThrottle = 200
                                                ))
                             )
                      )
                    )
                )
              ),
              fluidRow(
                box(width = 12, status = "info", solidHeader = TRUE,
                    uiOutput("geneSubTabs"))
              )
      )
    )
  )
)

# == == == == == == == == == == == == == == == == == == == == == == == == ==
# ============================== Server ====================================
# == == == == == == == == == == == == == == == == == == == == == == == == ==

server <- function(input, output, session) {
  options(shiny.maxRequestSize = 200 * 1024^2)
  
  # ------------------------------------------------------------------
  # Shared state across all server modules
  # ------------------------------------------------------------------
  state <- list(
    # Core data objects
    se_obj        = reactiveVal(NULL),
    de_df         = reactiveVal(NULL),
    varpart_obj   = reactiveVal(NULL),
    annotation_df = reactiveVal(NULL),
    
    # Organism & MSigDB routing
    se_organism           = reactiveVal(NULL),   # raw string from metadata(se)$organism
    db_species            = reactiveVal(NULL),   # "HS", "MM", or NULL (not in msigdbr)
    
    # Annotation-derived
    available_gsea_columns = reactiveVal(list()), # named list of functional columns
    ensembl_col            = reactiveVal(NULL)    # column name holding Ensembl IDs, or NULL
  )
  
  # Convenience reactive: organism string (same as state$se_organism but reactive)
  organism <- reactive({ state$se_organism() })
  
  # Wire helpers
  load_files_server(input, output, session, state)
  explore_contrast_server(input, output, session, state, organism)
  compare_contrast_server(input, output, session, state, organism)
  explore_gene_server(input, output, session, state, organism)
}

# == == == == == == == == == == == == == == == == == == == == == == == == ==
# ============================= Run App ====================================
# == == == == == == == == == == == == == == == == == == == == == == == == ==

shinyApp(ui, server)