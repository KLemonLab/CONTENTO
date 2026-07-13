FROM bioconductor/bioconductor_docker:RELEASE_3_19

ENV DEBIAN_FRONTEND=noninteractive

# System libraries commonly needed by tidyverse, plotly, DT, ggiraph, etc.
RUN apt-get update && apt-get install -y --no-install-recommends \
    zlib1g-dev \
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
RUN R -e " \
    pkgs <- c('shiny','shinydashboard','shinycssloaders','DT','plotly','tidyverse', \
              'purrr','viridis','RColorBrewer','ComplexUpset','scales','ggiraph', \
              'gridExtra','patchwork'); \
    install.packages(pkgs, repos = 'https://cloud.r-project.org'); \
    missing <- pkgs[!sapply(pkgs, requireNamespace, quietly = TRUE)]; \
    if (length(missing) > 0) stop('Failed to install: ', paste(missing, collapse = ', '))"


# Bioconductor packages
# fgsea declares SystemRequirements: C++11 in its DESCRIPTION, but its BH
# (Boost) dependency now requires C++14+ internals. Overriding CXX11STD
# forces anything that requests C++11 to compile with C++17 instead.
RUN mkdir -p /root/.R && echo "CXX11STD = -std=gnu++17" > /root/.R/Makevars

RUN R -e " \
    pkgs <- c('ComplexHeatmap','DESeq2','fgsea','msigdbr','circlize','SummarizedExperiment'); \
    BiocManager::install(pkgs, ask = FALSE, update = FALSE); \
    missing <- pkgs[!sapply(pkgs, requireNamespace, quietly = TRUE)]; \
    if (length(missing) > 0) stop('Failed to install: ', paste(missing, collapse = ', '))"

# App directory
WORKDIR /srv/shiny-server/app

# Copy app
COPY . .

# Avoid running as root
VOLUME ["inst"]

EXPOSE 3838

HEALTHCHECK CMD curl -f http://localhost:3838 || exit 1

CMD ["R","-e","shiny::runApp('/srv/shiny-server/app',host='0.0.0.0',port=3838)"]