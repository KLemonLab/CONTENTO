# Install BiocManager if not already present
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# CRAN packages
cran_packages <- c(
  "shiny", "shinydashboard", "shinycssloaders", "DT", "plotly", 
  "tidyverse", "purrr", "viridis", "RColorBrewer", 
  "scales", "tools", "grid", "ggiraph", "fgsea", "msigdbr", "pheatmap", "circlize"
)

# Install missing CRAN packages
install.packages(setdiff(cran_packages, installed.packages()[,"Package"]))

# Bioconductor packages
bioc_packages <- c("ComplexHeatmap", "DESeq2", "ComplexUpset")

# Install missing Bioconductor packages
BiocManager::install(setdiff(bioc_packages, installed.packages()[,"Package"]))
