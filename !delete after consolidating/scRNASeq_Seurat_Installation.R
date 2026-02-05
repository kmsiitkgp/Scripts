#!/usr/bin/env Rscript

.libPaths(new = "/hpc/home/kailasamms/R/x86_64-pc-linux-gnu-library/4.1")
.libPaths()
chooseCRANmirror(ind=1)

# NOTE: Install all packages from submit node NOT computing node. 

# NOTE: Some packages need hdf5r package which in turn needs HDF5 library files
# for installation as well as to work. Also, Rcpp which is needed by most basic
# R packages doesnt get installed if you skip "module load hdf5/1.8.18"
# (i) qrsh -> conda activate R -> module load hdf5/1.8.18 -> R -> install... 

# NOTE: If you get "libxml not found", "Cannot find xml2-config",
# "ERROR: configuration failed for package ‘XML’", then, 
# Exit R, then "conda activate R", then "conda install -c conda-forge r-xml"
# IMPORTANT: DO NOT UPDATE XML even if outdated. 

#******************************************************************************#
#                          INSTALL NECESSARY PACKAGES                          #
#******************************************************************************#

# Non-CRAN repository managing packages
install.packages(pkgs = c("BiocManager", "remotes"),
                 repos ='http://cran.us.r-project.org',
                 force = FALSE,
                 INSTALL_opts = '--no-lock')

# Data analysis packages
BiocManager::install(pkgs = c("AnnotationHub", "ensembldb", "fgsea", 
                              "clusterProfiler", "org.Hs.eg.db", "org.Mm.eg.db",
                              "DESeq2", "progeny", "dorothea", "viper", "sva",
                              "preprocessCore", "MAGeCKFlute"),
                     force = FALSE,
                     INSTALL_opts = '--no-lock')

install.packages(pkgs = c("msigdbr"),
                 repos ='http://cran.us.r-project.org',
                 force = FALSE,
                 INSTALL_opts = '--no-lock')

# Data wrangling packages
install.packages(pkgs = c("openxlsx", "dplyr", "tibble", "stringr", "purrr"),
                 repos ='http://cran.us.r-project.org',
                 force = FALSE,
                 INSTALL_opts = '--no-lock')

# Graph plotting packages
install.packages(pkgs = c("ggplot2", "cowplot", "viridis", "RColorBrewer"),
                 repos = 'http://cran.us.r-project.org',
                 force = FALSE,
                 INSTALL_opts = '--no-lock')

# Specialized Graph plotting packages
install.packages(pkgs = c("pheatmap", "VennDiagram", "survival", "survminer", 
                          "ggridges", "ggbeeswarm"),
                 repos = 'http://cran.us.r-project.org',
                 force = FALSE,
                 INSTALL_opts = '--no-lock')

# Single cell analysis packages
install.packages(pkgs = c("Seurat", "harmony"),
                 repos ='http://cran.us.r-project.org',
                 force = FALSE,
                 INSTALL_opts = '--no-lock')

# "BiocNeighbors" is key dependency for CellChat package
# "RcisTarget" is key dependency for SCENIC package
BiocManager::install(pkgs = c("BiocNeighbors", "RcisTarget", "UCell", 
                              "glmGamPoi", "scCustomize"),
                     force = FALSE,
                     INSTALL_opts = '--no-lock')

# Sometimes connection is timed out when connecting to GitHub. So, use a proxy.
library("httr")
httr::set_config(use_proxy("8.8.8.8", port = 8080))
remotes::install_github(repo = c("hhoeflin/hdf5r",
                                 "mojaveazure/seurat-disk",
                                 "aertslab/SCENIC",
                                 "aertslab/SCopeLoomR",
                                 "sqjin/CellChat",
                                 "stephens999/ashr",
                                 "chris-mcginnis-ucsf/DoubletFinder",
                                 "neurorestore/Augur"),
                        INSTALL_opts = '--no-lock')

old.packages(lib.loc = .libPaths(),
             repos = 'http://cran.us.r-project.org')

pkgs <- c("BiocManager", "remotes", "TCGAbiolinks", "ensembldb", "AnnotationHub",
          "affy", "lumi", "ChIPQC", "fgsea", "enrichplot", "clusterProfiler",
          "org.Hs.eg.db", "org.Mm.eg.db", "DESeq2", "progeny", "dorothea", 
          "viper", "msigdbr", "openxlsx", "dplyr", "tibble", "stringr", "purrr",
          "ggplot2", "cowplot", "viridis", "RColorBrewer", "pheatmap", 
          "ggridges", "VennDiagram", "survival", "survminer", "ggbeeswarm", "arrow", 
          "Seurat", "metap", "hdf5r", "SeuratDisk", "DoubletFinder", "Augur", 
          "SCENIC", "SCopeLoomR", "ktplots", "CellChat", "wordcloud", "UCell", 
          "infercnv", "illuminaHumanv4.db", "pathview", "colorspace", 
          "reticulate", "ashr", "multtest", "RcisTarget", "BiocNeighbors", 
          "harmony", "glmGamPoi")

# Display packages that couldn't be installed
cat("These pacakges have NOT yet been installed:", sort(pkgs[!(pkgs %in% installed.packages()[,1])]), sep="\n")

# arrow package is key dependency for SCENIC
# Install using compute node. It might fail using submit node.
#BiocNeighbors is key dependency for CellChat
Sys.setenv("ARROW_R_DEV"=TRUE, "LIBARROW_BINARY"=FALSE,
           "ARROW_WITH_ZSTD"="ON", "ARROW_DEPENDENCY_SOURCE"="BUNDLED")
install.packages(pkgs = c("arrow"),
                 repos = 'http://cran.us.r-project.org',
                 force = FALSE,
                 INSTALL_opts = '--no-lock')