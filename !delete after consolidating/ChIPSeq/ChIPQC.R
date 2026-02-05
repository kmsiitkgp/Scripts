#!/usr/bin/env Rscript

#********************STEP 1: SET UP THE WORKING ENVIRONMENT********************#

#*********************STEP 1A: INSTALL NECESSARY PACKAGES**********************#

.libPaths("/hpc/home/kailasamms/NGSTools/R_packages")
.libPaths()
chooseCRANmirror(ind=1)

# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("ChIPQC")

# install.packages("openxlsx")


#*********************STEP 1B: LOAD THE NECESSARY PACKAGES*********************#

library("ChIPQC")     # Useful for QC of MAC2 Peak call results
library("openxlsx")		# Useful for reading and writing xlsx files

#***************STEP 1C: STORE PATHS OF DIRECTORIES IN VARIABLES***************#

# Store path of parent directory containing individual folders for analysis
# Store path of results directory where seurat output will be stored
# Store path of scripts directory containing scripts, input files like cell
# cycle markers, cell markers, gene annotations etc..
# Store path of diagnostics directory where QC figures of Seurat and DESeq2 will
# be stored

parent_path <-  "/hpc/home/kailasamms/scratch/ChIPSeq/"
results_path <- "/hpc/home/kailasamms/scratch/ChIPSeq/ChIPQC_results/"
scripts_path <- "/hpc/home/kailasamms/scratch/ChIPSeq/scripts/"
diagnostics_path <- "/hpc/home/kailasamms/scratch/ChIPSeq/diagnostics/"

# parent_path <-  "C:/Users/KailasammS/Box/Saravana@cedars/Single Cell Scripts/test/"
# results_path <- "C:/Users/KailasammS/Box/Saravana@cedars/Single Cell Scripts/results/"
# scripts_path <- "C:/Users/KailasammS/Box/Saravana@cedars/[SCRIPTS] Single Cell/"
# diagnostics_path <- "C:/Users/KailasammS/Box/Saravana@cedars/Single Cell Scripts/results/"

#***************STEP 1D: READ THE METADATA***************#

# Read the meta data
meta_data <- read.xlsx(xlsxFile = paste0(scripts_path,"AR_ChIP_Metadata.xlsx"),
                       colNames = TRUE)

## Create ChIPQC object
chipObj <- ChIPQC(samples, annotation="hg19")