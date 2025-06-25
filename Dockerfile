# Base image with R and Seurat preinstalled (e.g., Satijalab's Seurat image)
FROM satijalab/seurat:latest

# Install system libraries required for Signac, rtracklayer, Rsamtools, etc.
RUN apt-get update && apt-get install -y \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libhdf5-dev \
    libbz2-dev \
    liblzma-dev \
    libz-dev \
    libhts-dev \
    libncurses5-dev \
    libncursesw5-dev \
    libfontconfig1-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libfreetype6-dev \
    libpng-dev \
    libtiff5-dev \
    libjpeg-dev \
    && apt-get clean && rm -rf /var/lib/apt/lists/*

# Install R packages from CRAN and Bioconductor
RUN R -e "install.packages(c('remotes', 'ggplot2', 'dplyr', 'data.table', 'Matrix'), repos='https://cloud.r-project.org')" \
 && R -e "if (!requireNamespace('BiocManager', quietly=TRUE)) install.packages('BiocManager')" \
 && R -e "BiocManager::install(c('GenomicRanges', 'rtracklayer', 'Rsamtools', 'BiocGenerics', 'IRanges', 'S4Vectors', 'GenomeInfoDb', 'BiocParallel'))" \
 && R -e "remotes::install_github('timoast/signac')"

# Set working directory inside the container
WORKDIR /usr/src/app

# Copy your R script into the image
COPY Final_Project.R .

# Run the script 
CMD ["Rscript", "Final_Project.R"]


