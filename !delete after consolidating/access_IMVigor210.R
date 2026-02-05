#!/usr/bin/env Rscript

#******************************************************************************#
#                          INSTALL NECESSARY PACKAGES                          #
#******************************************************************************#

# IMvigor210 packages
BiocManager::install("edgeR")
BiocManager::install("genefilter")
BiocManager::install("geneplotter")
install.packages(c("lsmeans", "spatstat"))
install.packages("https://www.bioconductor.org/packages//2.10/bioc/src/contrib/DESeq_1.8.3.tar.gz", repos = NULL)
install.packages("https://github.com/SiYangming/IMvigor210CoreBiologies/releases/download/v2.0.0/IMvigor210CoreBiologies_2.0.0.tar.gz", repos = NULL)

#******************************************************************************#
#                           LOAD NECESSARY PACKAGES                            #
#******************************************************************************#

library("openxlsx")
library("dplyr")
library("tibble")
library("edgeR")
library("lsmeans")
library("DESeq")
library("IMvigor210CoreBiologies")

#******************************************************************************#
#                                  IMPORT DATA                                 #
#******************************************************************************#

# Check all data available. Type q to exit interactive mode.
data()

# Load the data containing counts for ImVigor trial
utils::data(cds)

# Import count data
read_data <- data.frame(DESeq2::counts(cds)) %>%
  tibble::rownames_to_column("ENTREZID")

# Import meta data
meta_data <- data.frame(Biobase::pData(cds)) %>%
  tibble::rownames_to_column("Sample_id")
  
# Save the data in excel file
wb <- openxlsx::createWorkbook()
openxlsx::addWorksheet(wb, sheetName = "Raw_counts")
openxlsx::writeData(wb, sheet = "Raw_counts", x = read_data)
openxlsx::addWorksheet(wb, sheetName = "Clinical_data")
openxlsx::writeData(wb, sheet = "Clinical_data", x = meta_data)

openxlsx::saveWorkbook(wb, file = "/hpc/home/kailasamms/IMVigor210.xlsx", overwrite=TRUE)

# WE keep ENTREZID and map to GENE SYMBOLS at the end of DESeq2 analysis to 
# preserve most data. By mapping to ENSEMBL or GENE SYMBOL before DESeq2 analysis,
# we lose many genes. So, do this at last.

#******************************************************************************#