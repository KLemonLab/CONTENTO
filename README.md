# CONTENTO

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.21629492.svg)](https://doi.org/10.5281/zenodo.21629492)

> *CONTrast ExploratioN Toolkit for Omics*

## Overview

CONTENTO reads a SummarizedExperiment (SE) file with differential expression results for any number of contrasts, plus an optional annotation file with gene labels and gene sets (A), and explores them in three modules: Explore by Contrast, Compare Contrasts, and Explore by Gene (B).

![CONTENTO inputs and outputs](docs-src/images/contento_overview.png)

## Documentation

Full documentation, including tutorials for creating SE and annotation files, is available at: <https://klemonlab.github.io/CONTENTO/>

## Running CONTENTO

CONTENTO can be run in two ways. We recommend the container, since it ships with R, Bioconductor, and every dependency pre-installed, and behaves identically on any operating system.

CONTENTO is built with shiny and uses DESeq2, fgsea, msigdbr, ComplexHeatmap, circlize and ComplexUpset for its analyses and visualizations.

### *Option 1:* Run the container (recommended)

A container bundles the app together with everything it needs to run. You only need to install a container runtime once:

- **macOS and Windows:** [Podman Desktop](https://podman-desktop.io/), free and open source. Docker Desktop also works, but requires a paid subscription at larger organizations (more than 250 employees or more than \$10 million in annual revenue), which usually includes universities and hospitals.
- **Linux:** [Docker Engine](https://docs.docker.com/engine/install/) or Podman, both free and available from your distribution's package manager.

Then start CONTENTO from a terminal *(with Podman Desktop running)*:

``` bash
podman run --rm -p 3838:3838 ghcr.io/klemonlab/contento:latest
```

If you use Docker, replace `podman` with `docker`; the rest of the command is identical. The first run downloads the image, which takes a few minutes. When the terminal shows `Listening on http://0.0.0.0:3838`, open your web browser and navigate to:

```         
http://localhost:3838
```

Press `Ctrl+C` in the terminal to stop the app.

**Loading your data.** Upload your SE file in the app's sidebar, and an annotation file from the load panel if needed. Built-in annotations are included in the container and load automatically when `metadata(se)$annotation` names one of them. Uploaded files are only kept for the current session.

**Trying CONTENTO with example data.** The example SE files are not included in the container. Download them from [`inst/extdata/SE`](https://github.com/KLemonLab/CONTENTO/tree/main/inst/extdata/SE) and upload them in the sidebar. See [Create SE for CONTENTO](https://klemonlab.github.io/CONTENTO/create-se-file.html) for how they were built.

**Using your own annotation files.** To load your own annotation files automatically instead of uploading them each session, mount the folder that contains them and point `CONTENTO_ANNOTATION_DIR` to it:

``` bash
podman run --rm -p 3838:3838 \
    -v /path/to/annotations:/annotations:ro \
    -e CONTENTO_ANNOTATION_DIR=/annotations \
    ghcr.io/klemonlab/contento:latest
```

**Pinning a version.** `latest` always points to the newest release. For analyses you intend to publish, replace `latest` with a specific release version (see [Releases](https://github.com/KLemonLab/CONTENTO/releases)) and report that version in your methods:

``` bash
podman run --rm -p 3838:3838 ghcr.io/klemonlab/contento:1.0.0
```

**Building the image yourself.** A Dockerfile is included in the repository, so you can build and run the image locally from the cloned repository:

``` bash
podman build -t contento .
podman run --rm -p 3838:3838 contento
```

### *Option 2:* Run locally within R/RStudio

Some users may prefer to execute **CONTENTO** directly from an existing R installation, usually used in conjunction with the popular IDE [RStudio](https://posit.co/download/rstudio-desktop). This requires R ≥ 4.6 and Bioconductor ≥ 3.23, the versions used in the container.

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

## Citation

If you use this application in scientific work, please cite:

F. Escapa, I. & Lemon, K. P. "CONTENTO: CONTrast ExploratioN Toolkit for Omics". <https://github.com/KLemonLab/CONTENTO>

DOI: <https://doi.org/10.5281/zenodo.21629492>

This DOI always resolves to the latest release; the DOI for each specific version is listed on the Zenodo page.

## Funding

This work was supported by the National Institute of General Medical Sciences (NIGMS) of the National Institutes of Health under award number R35GM141806, and by the National Institute of Allergy and Infectious Diseases (NIAID) under award number U19AI157981. The content is solely the responsibility of the authors and does not necessarily represent the official views of the National Institutes of Health.

## License

This project is licensed under the GNU General Public License v3.0 or later (GPL-3.0-or-later).

## Contact

Issues, bug reports, and feature requests are welcome through the [GitHub issue tracker](https://github.com/KLemonLab/CONTENTO/issues).
