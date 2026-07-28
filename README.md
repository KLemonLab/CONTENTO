# CONTENTO

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.21629493.svg)](https://doi.org/10.5281/zenodo.21629493)

> *CONTrast ExploratioN Toolkit for Omics* 

## Intro

Modern omics experiments in microbial research rarely conform to simple two-group comparisons. Experimental designs frequently involve multiple variables tested simultaneously (microbial species, growth conditions, treatments, and their interaction effects) generating multiple pairwise contrasts that are difficult to navigate with standard tools. Identifying which genes or pathways are uniquely regulated in one condition, shared across several, or driven by interaction effects requires an analytical framework that scales with experimental complexity. We present CONTENTO (CONTrast ExploratioN Toolkit for Omics), an interactive R Shiny dashboard designed to help researchers explore, visualize, and compare results across any number of contrasts, with the broader vision of supporting multiple omics data types. CONTENTO currently focuses on DESeq2-processed RNAseq data, supporting both microbial and host transcriptomes within a unified interface, making it well-suited for experimental systems ranging from in vitro microbial cultures to complex in vivo and coculture host-microbe models.

CONTENTO integrates per-contrast visualization of differential expression results, interactive volcano plots, expression heatmaps, and gene set enrichment analysis. Gene sets are drawn from MSigDB collections for host transcriptomes, while microbial transcriptomes are supported through COG and KEGG functional annotations, enabling pathway-level enrichment analysis across both organisms in dual-host-microbe experimental systems. CONTENTO’s multi-contrast comparison framework allows users to identify differentially expressed genes and pathway overlap across any combination of conditions, revealing shared and condition-specific transcriptional responses. Single-gene expression visualization across user-defined experimental metadata enables detailed interrogation of individual genes simultaneously across contrasts. CONTENTO provides an accessible and scalable platform for extracting biological insight from experimentally complex omics datasets.

## Documentation

Full documentation, including tutorials for creating SE and annotation files, is available at: <https://klemonlab.github.io/CONTENTO/>

## Running CONTENTO

CONTENTO can be executed in three different ways depending on your computational environment and familiarity with containerized software. We strongly recommend the Docker-based workflows, as they guarantee complete reproducibility across operating systems and avoid dependency conflicts.

### *Option 1:* Run using Docker (recommended)

If Docker is not yet installed on your system, we recommend you install Docker Desktop from the [official website](https://www.docker.com/products/docker-desktop/), since it provides a friendly graphic interface to manage your images and keep your Docker Engine up-to-date. If you already have Docker installed, running CONTENTO requires a single command from your terminal *(make sure your Docker Desktop is running)*:

``` bash
docker run \
    --rm \
    -p 3838:3838 \
    -v $(pwd)/inst:/srv/shiny-server/app/inst \
    ghcr.io/klemonlab/contento:latest
```

The terminal window will show the app startup and soon will start accepting connections. When you see the message that the app is listening, open your web browser (any browser) and navigate to:

```         
http://localhost:3838
```

This approach provides a fully reproducible computational environment with all R, Bioconductor, and system dependencies pre-installed.

> **Rebuilding image from scratch**
>
> A Dockerfile is included in the repository so that users familiar with building container images can self-build with `docker build -t contento .` if necessary.

### *Option 2:* Run using Docker Compose

For users who prefer a persistent local setup and the convenience of a preconfigured container mount, **CONTENTO** can also be launched using [Docker Compose](https://docs.docker.com/compose/). To do this, make sure you have the `compose` plugin installed (it is bundled with modern distributions of Docker Desktop by default) by running

``` bash
docker compose --help
```

The helper string for the `compose` plugin should come up.

Then, clone the CONTENTO repository:

``` bash
git clone https://github.com/klemonlab/contento.git
```

To start the application, navigate to the working directory and simply run:

``` bash
cd contento
docker compose up
```

Or run it in the background:

``` bash
docker compose up -d
```

Once the container starts, access the application again at:

```         
http://localhost:3838
```

### *Option 3:* Run locally within R/RStudio

Some users may prefer to execute **CONTENTO** directly from an existing R installation, usually used in conjunction with the popular IDE [RStudio](https://posit.co/download/rstudio-desktop).

Clone the CONTENTO repository:

``` bash
git clone https://github.com/klemonlab/contento.git
```

After opening the project (`CONTENTO.Rproj`) within Rstudio, install the dependencies by running:

``` r
install.packages(c(
    "shiny",
    "shinydashboard",
    "shinycssloaders",
    "DT",
    "plotly",
    "tidyverse",
    "purrr",
    "viridis",
    "RColorBrewer",
    "ComplexUpset",
    "scales",
    "ggiraph",
    "gridExtra",
    "patchwork"
))

if (!requireNamespace("BiocManager"))
    install.packages("BiocManager")

BiocManager::install(c(
    "ComplexHeatmap",
    "DESeq2",
    "fgsea",
    "msigdbr",
    "circlize",
    "SummarizedExperiment"
))
```

Open `app.R` in RStudio and click the interactive **:arrow_forward: Run App** button, or execute:

``` r
shiny::runApp("app.R")
```

The application will then become available in your web browser.

------------------------------------------------------------------------

### Which option should I choose?

| Method | Recommended for | Widely Reproducible | Requires installation |
|----|----|----|----|
| Docker | Most users | ✓ | Docker |
| Docker Compose | Persistent local deployments | ✓ | Docker + Compose |
| R/RStudio | RStudio users and developers | Depends on environment | R + package installation |

## Citation

If you use this application in scientific work, please cite:

Isabel FE, et al. "CONTENTO: CONTrast ExploratioN Toolkit for Omics". <https://github.com/KLemonLab/CONTENTO> 

DOI: https://doi.org/10.5281/zenodo.21629493

## License

This project is covered under the GPL-3 license.

## Contact

Issues, bug reports, and feature requests are welcome through the GitHub issue tracker.
