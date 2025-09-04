FROM r-base:4.5.0

LABEL maintainer = "Método añadido por Julieth López <julalopezcas@unal.edu.co>"

USER root

WORKDIR /

# Install system dependencies for R packages
# These packages are required for various Bioconductor and CRAN packages
RUN apt-get update && apt-get install -y \
    time \
    libgsl-dev \
    gsl-bin \
    libxml2-dev \
    libssl-dev \
    libcurl4-openssl-dev \
    libfontconfig1-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libfreetype6-dev \
    libpng-dev \
    libtiff5-dev \
    libjpeg-dev \
    libgit2-dev \
    libblas-dev \
    liblapack-dev \
    gfortran \
    build-essential \
    git \
    && rm -rf /var/lib/apt/lists/*

# Install BiocManager
RUN R -e "install.packages('BiocManager', repos='https://cloud.r-project.org/')"
# Set Bioconductor version to 3.21
RUN R -e "BiocManager::install(version='3.21')"

# Install Bioconductor packages
RUN R -e "BiocManager::install(c( \
    'DirichletMultinomial', \
    'TFBSTools', \
    'motifmatchr', \
    'GenomicRanges', \
    'EnsDb.Hsapiens.v86', \
    'BSgenome.Hsapiens.UCSC.hg38', \
    'JASPAR2020', \
    'JASPAR2024', \
    'devtools', \
    'UCell', \
    'GenomicRanges', \
    'GeneOverlap' \
), ask=FALSE, update=FALSE)"

# Install CRAN packages (xgboost, tidyverse, cowplot, patchwork, magrittr, igraph)
##RUN R -e "install.packages('https://cran.r-project.org/src/contrib/BTR_1.2.4.tar.gz', type = 'source', dependencies = TRUE)"
RUN R -e "install.packages(c( \
    'xgboost', \
    'tidyverse', \
    'cowplot', \
    'patchwork', \
    'magrittr', \
    'igraph' \
), repos='https://cloud.r-project.org/')"

# Install additional Bioconductor/CRAN packages for hdWGCNA workflow
RUN R -e "BiocManager::install(c('Seurat', 'WGCNA'), ask=FALSE, update=FALSE)"

# Install packages with devtools
RUN R -e "devtools::install_github('smorabit/hdWGCNA', ref='dev')"

COPY runTF_hdWGCNA.R /

RUN mkdir data/