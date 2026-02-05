#!/usr/bin/env Rscript

#******************************************************************************#
#                       STEP 2: LOAD TCGA EXPRESSION DATA                      #
#******************************************************************************#

# # (ii) TCGAbiolinks (OK as it is updated frequently)
# # Check GDC server status using the api https://api.gdc.cancer.gov/status
# GDC_info <- TCGAbiolinks::getGDCInfo()
# GDC_info
# 
# # Check all the available projects at TCGA
# GDC_projects <- TCGAbiolinks::getGDCprojects()
# 
# # Check project summary
# GDC_project_summary <- TCGAbiolinks:::getProjectSummary(project="TCGA-BLCA",
#                                                         legacy=FALSE)
# 
# # Download clinical data using TCGAbiolinks
# TCGAbiolinks_clinical <- TCGAbiolinks::GDCquery_clinic(project="TCGA-BLCA",
#                                                        type="clinical",
#                                                        save.csv=TRUE)
# 
# # (iii) RTCGA (NOT RECOMMENDED as it is not updated frequently) 
# # Download clinical data using RTCGA
# RTCGA_clinical <- RTCGA::survivalTCGA(BLCA.clinical,
#                                       extract.cols="admin.disease_code",
#                                       extract.names=FALSE,
#                                       barcode.name="patient.bcr_patient_barcode",
#                                       event.name="patient.vital_status",
#                                       days.to.followup.name="patient.days_to_last_followup",
#                                       days.to.death.name="patient.days_to_death")
#
# # You will notice the clinical data is very different. Two main differences:
# # (i) "RTCGA_clinical$patient.vital_status" is same as
# # "TCGAbiolinks_clinical$year_of_death". Clearly, RTCGA is wrong because you can
# # notice some patients are dead based on "TCGAbiolonks_clinical$vital_status"
# # (ii) "RTCGA_clinical$times" is combination of
# # "TCGAbiolonks_clinical$days_to_last_follow_up" and "TCGAbiolonks_clinical$days_to_death".


#******************************************************************************#
#              STEP 4A: NORMALIZE EACH CANCER DATA INDIVIDUALLY                #
#******************************************************************************#

parent_path <- "/hpc/home/kailasamms/scratch/TCGA_GDC/"
results_path <- parent_path
species <- "Homo sapiens"
Variable <- "Treatment"
Comparisons <- list(Target=c("Radiation"),
                    Reference=c("No Treatment"))
padj.cutoff <- 0.1
lfc.cutoff <- 0  
annotations <- get_annotations(species)

meta_data_full <- openxlsx::read.xlsx(xlsxFile=paste0(parent_path, "Meta_data_TCGA.xlsx"))

for(proj in unique(meta_data_full$Project_ID)){
  
  read_data <- openxlsx::read.xlsx(xlsxFile=paste0(parent_path, "Read_data_", proj, ".xlsx")) %>% 
    dplyr::mutate(SYMBOL = make.names(SYMBOL, unique=TRUE))
  
  meta_data <- meta_data_full %>% 
    dplyr::filter(make.names(Sample_ID) %in% colnames(read_data)) %>%
    dplyr::mutate(Batch=Project_ID)
  
  # Save the cancer specific meta_data
  wb <- openxlsx::createWorkbook()
  openxlsx::addWorksheet(wb, sheetName="Metadata")
  openxlsx::writeData(wb, sheet="Metadata", x=meta_data)
  openxlsx::saveWorkbook(wb, file=paste0(parent_path, "Meta_data_", proj, ".xlsx"), overwrite=TRUE)
  
  meta_data <- prep_metadata(meta_data, Variable)
  read_data <- prep_readdata(read_data, meta_data)
  l <- check_data(read_data, meta_data)
  meta_data <- l[[2]]
  read_data <- l[[1]]
  sva_dds <- svaseq_batch(read_data, meta_data)
  n <- 1
  
  #Perform DESeq2() using in-built batch modelling
  approach <- "DESeq2"
  if (length(unique(meta_data$Batch)) > 1){
    dds <- DESeq2::DESeqDataSetFromMatrix(countData=read_data,
                                          colData=meta_data, 
                                          design=~ Batch+id)
  } else {
    dds <- DESeq2::DESeqDataSetFromMatrix(countData=read_data,
                                          colData=meta_data, 
                                          design=~ id)
  }
  dds <- run_deseq2(dds, meta_data, annotations, Comparisons, n, lfc.cutoff, padj.cutoff, approach)
  deseq2_norm_counts(dds, annotations, approach) # batch corrected if you more than 1 batch
  plot_qc(dds, meta_data, approach)
  
  # Perform DESeq2() using sva modelled surrogate variables SV1 and SV2
  approach <- "sva_modelled"
  sva_dds <- run_deseq2(sva_dds, meta_data, annotations, Comparisons, n, lfc.cutoff, padj.cutoff, approach)
  # calc_norm_counts(sva_dds, annotations, approach)   # uncorrected
  svaseq_norm_counts(sva_dds, annotations, approach)   # sva seq batch corrected
  plot_qc(sva_dds, meta_data, approach)
}


#*****************************************************************************#

# The below 2 posts recommend using all samples from a single experiment for 
# normalizing but avoiding using all samples from multiple experiments.
# https://support.bioconductor.org/p/92879/, https://www.biostars.org/p/9560478/
# https://support.bioconductor.org/p/59711/
# https://support.bioconductor.org/p/98765/
# https://www.biostars.org/p/336298/
# https://support.bioconductor.org/p/75260/

# LINK #1: In general, it is recommended to use all the available samples from 
# an experiment in the analysis, even if you are not interested in differential 
# expression for some of those samples. The main reason for this is that having
# more samples allows more accurate estimation of the gene dispersion values.

# LINK #2: If you are not going to test for differential expression between
# experiments, then there is no purpose in normalizing them together. The more
# worrying problem with analyzing all your experiments as a single data set is 
# that a single dispersion value will be estimated for each gene across all 
# experiments. This is only ok if you believe that every gene has equal 
# biological variability in all your experiments, which is unlikely to be case.

# Since, it is not appropriate to merge all samples from 33 different TCGA 
# projects and normalize them together due to issues explained above, we will
# perform normalization on each TCGA project individually ASSUMING there are
# no batch effects within each experiment.

# Get list of tsv files we have in folder
tsv_files <- list.files(data_path)
l <- c()
for (proj in unique(sample_sheet$Project.ID)){
  
  # Get list of tsv files for each Project_ID i.e. each cancer from Sample sheet
  files <- sample_sheet %>% 
    dplyr::filter(Project.ID == proj) %>% 
    dplyr::select(File.Name) %>%
    unlist(use.names=FALSE)
  
  # Get list of tsv files we actually have for each cancer
  files <- setdiff(files, tsv_files)
  
  l <- c(l,files)
}