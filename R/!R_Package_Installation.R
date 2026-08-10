#!/usr/bin/env Rscript

.libPaths(new = "/hpc/home/kailasamms/R/x86_64-pc-linux-gnu-library/4.1")
.libPaths()
chooseCRANmirror(ind=1)

# NOTE: Install all packages from submit node NOT computing node. 

# NOTE: Some packages need hdf5r package which in turn needs HDF5 library files
# for installation. Others like SSGSEA need gcc compiler. Rcpp which is needed 
# by most R packages doesnt get installed if you skip "module load hdf5/1.8.18"
# So, before installing do the following:
# conda activate R
# module load hdf5/1.8.18
# module load gcc/11.1.0
# module load zlib/1.3
# module load libtiff/4.2.0
# module load cmake/3.19.5        #3.2.3 doesnt work
# R

# NOTE: If you get "libxml not found", "Cannot find xml2-config",
# "ERROR: configuration failed for package ‘XML’", then, 
# Exit R, then "conda activate R", then "conda install -c conda-forge r-xml"
# IMPORTANT: DO NOT UPDATE XML even if outdated. 

#******************************************************************************#
#                          INSTALL NECESSARY PACKAGES                          #
#******************************************************************************#

### 1. Non-CRAN Repository Managers (CRAN)
# Installs BiocManager and remotes (needed for Bioconductor and GitHub packages)
install.packages(pkgs = c("BiocManager", "remotes"),
                 repos = "https://cloud.r-project.org",
                 INSTALL_opts = '--no-lock')

# Load httr library to set up proxy configuration before GitHub installs
library("httr")

# NOTE: The proxy settings below (8.8.8.8:8080) are placeholders. 
# You may need to remove this line or replace it with your organization's specific proxy.
# httr::set_config(use_proxy("8.8.8.8", port = 8080))


### 2. Annotation Packages (Bioconductor)
BiocManager::install(pkgs = c("AnnotationHub", "ensembldb", 
                              "org.Hs.eg.db", "org.Mm.eg.db"),
                     update = FALSE,
                     INSTALL_opts = '--no-lock')


### 3. Functional Data Analysis Packages (Bioconductor)
BiocManager::install(pkgs = c("fgsea", "clusterProfiler", "decoupleR", 
                              "DESeq2", "sva", "GSVA", "glmGamPoi", "tximport",
                              "ashr", "OmnipathR"),
                     update = FALSE,
                     INSTALL_opts = '--no-lock')


### 4. Single-Cell Data Analysis Packages (CRAN)
install.packages(pkgs = c("Seurat", "harmony", "hdf5r", "clustree"),
                 repos = 'https://cloud.r-project.org',
                 INSTALL_opts = '--no-lock')


### 5. Single-Cell Data Analysis Packages (Bioconductor)
BiocManager::install(pkgs = c("infercnv", "UCell", "scDblFinder", "RcisTarget",
                              "DropletUtils"),
                     INSTALL_opts = '--no-lock')


### 6. Single-Cell Data Analysis Packages (GitHub via remotes)
remotes::install_github(repo = c("mojaveazure/seurat-disk",
                                 "jinworks/CellChat", 
                                 "satijalab/seurat-wrappers", 
                                 "immunogenomics/presto", 
                                 "chris-mcginnis-ucsf/DoubletFinder",
                                 "satijalab/seurat-data"),
                        INSTALL_opts = '--no-lock')


### 7. Microarray & GeoMx Analysis Packages (Bioconductor)
BiocManager::install(pkgs = c("oligo", "oligoData", "illuminaHumanv4.db",
                              "hgu133plus2.db", "GEOquery", "affy", "lumi",
                              "NanoStringNCTools", "GeomxTools", 
                              "GeoMxWorkflows"),
                     update = FALSE,
                     INSTALL_opts = '--no-lock')


### 8. Data Wrangling Packages (CRAN)
install.packages(pkgs = c("openxlsx", "dplyr", "tibble", "stringr", "purrr"),
                 repos = 'https://cloud.r-project.org',
                 INSTALL_opts = '--no-lock')


### 9. Visualization & Plotting Packages (CRAN)
install.packages(pkgs = c("ggplot2", "ggplotify", "ggrepel", "ggpubr", 
                          "ggfortify", "ggridges", "ggbeeswarm", "pheatmap", 
                          "VennDiagram", "survival", "survminer", "UpSetR", 
                          "umap", "plot3D", "cowplot", "viridis", 
                          "RColorBrewer", "colorspace", "ragg"),
                 repos = 'https://cloud.r-project.org',
                 INSTALL_opts = '--no-lock')


### 10. Visualization Packages (Bioconductor)
BiocManager::install(pkgs = c("enrichplot", "ComplexHeatmap"),
                     update = FALSE,
                     INSTALL_opts = '--no-lock')

# List outdated packages
utils::old.packages(lib.loc = .libPaths(),
             repos = 'https://cloud.r-project.org')

# Update outdated packages
utils::update.packages(lib.loc = .libPaths(),
                    repos = 'https://cloud.r-project.org')

# List all packages
all_pkgs <- c(
  # === 1. Core Repository Managers (CRAN) ===
  "BiocManager", # For installing Bioconductor packages
  "remotes",     # For installing GitHub packages
  "httr",        # For handling proxy/http settings
  
  # === 2. Annotation and Database Packages (Bioconductor) ===
  "AnnotationHub",
  "ensembldb",
  "org.Hs.eg.db", # Human annotation
  "org.Mm.eg.db", # Mouse annotation
  
  # === 3. Functional & Bulk Data Analysis (Bioconductor) ===
  "fgsea",
  "clusterProfiler",
  "decoupleR",
  "DESeq2",
  "sva",
  "GSVA",
  "glmGamPoi",
  "tximport",
  "ashr",
  "OmnipathR",
  
  # === 4. Single-Cell Data Analysis (CRAN/Bioconductor/GitHub) ===
  # -- CRAN --
  "Seurat",
  "harmony",
  "hdf5r",
  "clustree",
  # -- Bioconductor --
  "infercnv",
  "UCell",
  "scDblFinder",
  "RcisTarget",
  "DropletUtils",
  # -- GitHub Repositories (install via remotes) --
  "SeuratDisk",
  "CellChat",
  "SeuratWrappers",
  "presto",
  "DoubletFinder",
  "SeuratData",
  
  # === 5. Microarray & Spatial (GeoMx) Analysis (Bioconductor) ===
  "oligo",
  "oligoData",
  "illuminaHumanv4.db",
  "hgu133plus2.db",
  "GEOquery",
  "affy",
  "lumi",
  "NanoStringNCTools",
  "GeomxTools",
  "GeoMxWorkflows",
  
  # === 6. Data Wrangling & Manipulation (CRAN) ===
  "openxlsx",
  "dplyr",
  "tibble",
  "stringr",
  "purrr",
  
  # === 7. Visualization & Plotting (CRAN/Bioconductor) ===
  # -- CRAN --
  "ggplot2",
  "ggplotify",
  "ggrepel",
  "ggpubr",
  "ggfortify",
  "ggridges",
  "ggbeeswarm",
  "pheatmap",
  "VennDiagram",
  "survival",
  "survminer",
  "UpSetR",
  "umap",
  "plot3D",
  "cowplot",
  "viridis",
  "RColorBrewer",
  "colorspace",
  "ragg",
  # -- Bioconductor --
  "enrichplot",
  "ComplexHeatmap"
)

# Display packages that couldn't be installed
cat("These pacakges have NOT yet been installed:", sort(all_pkgs[!(all_pkgs %in% installed.packages()[,1])]), sep="\n")


#BiocManager::install(pkgs = c("BiocNeighbors", ),
# "aertslab/SCENIC",
# "aertslab/SCopeLoomR",
# "neurorestore/Augur"
# "ChIPQC",  "metap", "ktplots", , "wordcloud",, 
# "pathview", "multtest",  
