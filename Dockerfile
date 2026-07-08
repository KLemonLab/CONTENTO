FROM bioconductor/bioconductor_docker:RELEASE_3_19

ENV DEBIAN_FRONTEND=noninteractive

# System libraries commonly needed by tidyverse, plotly, DT, ggiraph, etc.
RUN apt-get update && apt-get install -y --no-install-recommends \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libfontconfig1-dev \
    libfreetype6-dev \
    libpng-dev \
    libtiff5-dev \
    libjpeg-dev \
    libcairo2-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libgit2-dev \
    libglpk-dev \
    pandoc \
    pandoc-citeproc \
    && rm -rf /var/lib/apt/lists/*

# CRAN packages
RUN R -e "install.packages(c( \
    'shiny', \
    'shinydashboard', \
    'shinycssloaders', \
    'DT', \
    'plotly', \
    'tidyverse', \
    'purrr', \
    'viridis', \
    'RColorBrewer', \
    'ComplexUpset', \
    'scales', \
    'ggiraph' \
), repos = 'https://cloud.r-project.org')"

# Bioconductor packages
RUN R -e "BiocManager::install(c( \
    'ComplexHeatmap', \
    'DESeq2', \
    'fgsea', \
    'msigdbr', \
    'circlize', \
    'SummarizedExperiment' \
), ask = FALSE, update = FALSE)"

# App directory
WORKDIR /srv/shiny-server/app

# Copy app
COPY . .

# Avoid running as root
VOLUME ["inst"]

EXPOSE 3838

HEALTHCHECK CMD curl -f http://localhost:3838 || exit 1

CMD ["R","-e","shiny::runApp('/srv/shiny-server/app',host='0.0.0.0',port=3838)"]