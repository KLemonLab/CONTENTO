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
  dashboardHeader(
    title = strong("CONTENTO"),
    tags$li(class = "dropdown",
            tags$a(icon("github"), " GitHub",
                   href = "https://github.com/KLemonLab/CONTENTO",
                   target = "_blank")),
    tags$li(class = "dropdown",
            tags$a(icon("flask"), " Lemon Lab",
                   href = "https://www.bcm.edu/research/faculty-labs/katherine-lemon-lab",
                   target = "_blank"))
  ),

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
    h4("Differential Expression Settings", style = "padding-left: 20px;"),
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
                column(width = 2,
                       tags$img(src = "logo_hex.svg", style = "width: 100%; max-width: 160px; margin-top: 10px;")),
                column(width = 10,
                       box(title = strong("CONTENTO: CONTrast ExploratioN Toolkit for Omics"), width = 12,
                           status = "info", solidHeader = TRUE,
                           p(strong("CONTENTO"), " is an interactive dashboard for exploring and comparing differential ",
                             "expression and functional enrichment results across any number of contrasts from ",
                             "DESeq2-processed RNAseq data. Supporting both microbial and host transcriptomes, ",
                             "it is a great resource for dual host-microbe experimental designs."),
                           p("Upload a", strong("SummarizedExperiment (SE) object"), "on the left sidebar to begin. ",
                             "An annotation file is optional but unlocks gene labels, functional enrichment, ",
                             "and (for bacteria) neighborhood analysis."),
                           p("Once an annotation is loaded, use ", strong("Select Gene symbol/label column"), " to choose ",
                             "which annotation column is displayed as the gene label throughout tables and plots."),
                           p("Use the ", strong("Differential Expression Settings"), " in the sidebar to choose your ",
                             "fold-change column and set fold-change and FDR cutoffs. These apply globally across ",
                             "all tabs and update DE calls, tables, and plots in real time."),
                           tags$div(style = "display: flex; gap: 10px; flex-wrap: wrap; margin-top: 10px;",
                                    tags$a(icon("file-alt"), "How to install CONTENTO",
                                           href = "articles/install.html", target = "_blank",
                                           class = "btn btn-outline-secondary"),
                                    tags$a(icon("file-alt"), " How to create your SE file",
                                           href = "articles/create-se-file.html", target = "_blank",
                                           class = "btn btn-outline-secondary"),
                                    tags$a(icon("file-alt"), " How to create annotation files",
                                           href = "articles/create-genome-annotations.html", target = "_blank",
                                           class = "btn btn-outline-secondary")
                           )
                       ))
              ),
              fluidRow(
                style = "margin-top: 15px;",
                column(width = 4,
                       actionButton("go_contrast", class = "info-card-btn",
                                    label = tagList(
                                      h4(icon("chart-line"), strong(" Explore by Contrast")),
                                      p("DE table, volcano plot, and single-contrast GSEA.")
                                    ))),
                column(width = 4,
                       actionButton("go_compare", class = "info-card-btn",
                                    label = tagList(
                                      h4(icon("exchange-alt"), strong(" Compare Contrast")),
                                      p("DEG overlap, heatmap, multi-contrast GSEA, GESECA.")
                                    ))),
                column(width = 4,
                       actionButton("go_gene", class = "info-card-btn",
                                    label = tagList(
                                      h4(icon("dna"), strong(" Explore by Gene")),
                                      p("Gene info, expression plot, variance decomposition.")
                                    )))
              )
      ),
      
      ## ---- Explore by Contrast -----
      tabItem(tabName = "contrast",
              fluidRow(
                box(width = 12, status = "info", solidHeader = TRUE,
                    collapsible = TRUE, collapsed = FALSE,
                    title = tagList(icon("chart-line"), strong(" Explore by Contrast")),
                    fluidRow(
                      column(width = 4,
                             h5(strong("Choose a contrast to explore results:"))),
                      column(width = 8,
                             div(style = "padding-left: 10px; padding-top: 10px;", uiOutput("contrastSelect")))
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
                    title = tagList(icon("exchange-alt"), strong(" Compare Contrast")),
                    fluidRow(
                      column(width = 4,
                             h5(strong("Select multiple contrasts to compare results:")),
                             uiOutput("multiContrastSelectHeader")),
                      column(width = 8,
                             div(style = "padding-left: 10px; padding-top: 10px;", uiOutput("multiContrastCheckboxes")))
                    )
                )
              ),
              fluidRow(
                box(width = 12, status = "info", solidHeader = TRUE,
                    uiOutput("compareSubTabs"))
              )
      ),
      
      ## ---- Explore by Gene -----
      tabItem(tabName = "gene",
              fluidRow(
                box(width = 12, status = "info", solidHeader = TRUE,
                    collapsible = TRUE, collapsed = FALSE,
                    title = tagList(icon("dna"), strong(" Explore by Gene")),
                    fluidRow(
                      column(width = 4,
                             h5(strong("Choose a gene to explore:"))),
                      column(width = 8,
                             div(style = "padding-left: 10px; padding-top: 10px;", uiOutput("geneSelect")))
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
  
  observeEvent(input$go_contrast, updateTabItems(session, "tabs", "contrast"))
  observeEvent(input$go_compare,  updateTabItems(session, "tabs", "compare"))
  observeEvent(input$go_gene,     updateTabItems(session, "tabs", "gene"))

  # Shared state across all server modules
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