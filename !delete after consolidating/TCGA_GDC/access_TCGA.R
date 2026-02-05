#******************************************************************************#
#     STEP 1 (INDIVIDUAL CANCER): DOWNLOAD TCGA CLINICAL DATA & RAW COUNTS     #
#******************************************************************************#

# There are multiple ways to download TCGA data. 
# (i) Directly download from GDC Portal (SAFEST, BEST & EASIEST) 
# (ii) TCGAbiolinks (OK as it is updated frequently)
# (iii) RTCGA (NOT RECOMMENDED as it is not updated frequently)

# This script covers ONLY (i)
# Go to https://portal.gdc.cancer.gov/ ; click "projects" and filter as below:
# "Program"               -> "TCGA" ; 
# "Data Category"         -> "sequencing reads"
# "Experimental Strategy" -> "RNA-Seq"
# "Primary Site"          -> "Bladder" (or) choose all 33 TCGA projects
# Click "Save New cohort", give a name and select your cohort.

# Now, click "Repository" and filter as below:
# "Experimental Strategy" -> "RNA-Seq"
# "Workflow Type"         -> "STAR-Counts"
# "Access"                -> "open"

# Verify that filenames are star_gene_counts.tsv
# (i) Click "Manifest" to download manifest info
# (ii) Click "Download Associated Data" and download "Sample sheet" to get info 
# linking sample name to manifest info.
# NOTE: "Sample Sheet" is NECESSARY to match the downloaded count files with 
# correct patient as the downloaded files have random names.

# (iii) Click on "11428 cases" (top right corner), "Clinical" and "tsv" to 
# download clinical data for all 33 TCGA projects. The clinical data will be in
# tar.gz format. Extract the tsv files.

# Next, upload the manifest txt file to "projects/TCGA_GDC" folder in HPC cluster. 
# "projects/TCGA_GDC" folder also has "01_TCGA_GDC.sh" file. 
# Adjust manifest file name in "01_TCGA_GDC.sh" and run it to download count data.
# Follow steps in "01_TCGA_GDC.sh" to move files to appropriate location.

# Once download from GDC portal is complete and all files have been moved to 
# appropriate location, run "02_TCGA_GDC.sh"

#******************************************************************************#
#       STEP 2 (INDIVIDUAL CANCER): COMBINE CLINICAL DATA WITH FILE INFO       #
#******************************************************************************#

source("/hpc/home/kailasamms/projects/RNASeq/RNASeq_DESeq2_Functions.R")

data_path <- "/hpc/home/kailasamms/scratch/ICB_TCGA_datasets/original/"
output_path <- "/hpc/home/kailasamms/scratch/ICB_TCGA_datasets/processed/"

# sample_sheet.tsv                          : 11499 gene count files with Case.ID
# clinical.tsv                              : clinical data of 11428  patients
# follow_up.tsv                             : follow_up data of 11283 patients
# TCGA-CDR-SupplementalTableS1.xlsx         : clinical data of 11160  patients
# clinical_PANCAN_patient_with_followup.tsv : follow_up data of 10956 patients

# NOTE: TCGA-CDR-SupplementalTableS1.xlsx & clinical_PANCAN_patient_with_followup.tsv
# are PanCancer data and have nicely summarized OS, DSS, DFI, PFI and other info 
# for every patient which are unavailable in clinical.tsv and follow_up.tsv

# NOTE: We ignore family_history.tsv, exposure.tsv, pathology_detail.tsv as they
# are mostly blank

# NOTE: Since clinical.tsv and follow_up.tsv have more patient info than 
# PanCancer data, we are trying to merge all clinical data to get as much 
# complete info as possible

sample_sheet <- utils::read.table(file= paste0(data_path, "gdc_sample_sheet.2025-03-06.tsv"), 
                                  header=TRUE, sep="\t", quote="", skip=0, fill=TRUE)

clinical <- utils::read.table(file= paste0(data_path, "clinical.tsv"), 
                              header=TRUE, sep="\t", quote="", skip=0, fill=TRUE) %>%
  dplyr::distinct_at("case_submitter_id", .keep_all = TRUE)

follow_up <- utils::read.table(file= paste0(data_path, "follow_up.tsv"), 
                               header=TRUE, sep="\t", quote="", skip=0, fill=TRUE) %>%
  dplyr::distinct_at("case_submitter_id", .keep_all = TRUE)

# Merge clinical.tsv and follow_up.tsv
df1 <- clinical %>% 
  dplyr::full_join(follow_up %>% dplyr::select(everything(), -c("case_id", "project_id")), 
                   by=c("case_submitter_id"="case_submitter_id")) %>%
  dplyr::select(case_id, case_submitter_id, project_id)

# Merge the two PanCancer data files
pan_tcga_meta <- openxlsx::read.xlsx(xlsxFile=paste0(data_path, "TCGA-CDR-SupplementalTableS1.xlsx"))
pan_tcga_follow <- utils::read.table(file = paste0(data_path, "clinical_PANCAN_patient_with_followup.tsv"),
                                     header=TRUE, sep="\t", quote="", skip=0, fill=TRUE)

df2 <- pan_tcga_meta %>%
  dplyr::full_join(pan_tcga_follow %>% dplyr::select(bcr_patient_barcode, c(setdiff(colnames(pan_tcga_follow), colnames(pan_tcga_meta)))), 
                   by=c("bcr_patient_barcode"="bcr_patient_barcode"))

# Merge df1 and df2
df3 <- df1 %>%
  dplyr::full_join(df2, by=c("case_submitter_id"="bcr_patient_barcode"))

# Merge sample_sheet and df3
metadata <- df3 %>%
  dplyr::full_join(sample_sheet, by=c("case_submitter_id"="Case.ID")) 

# Remove all patients who dont have expression data
metadata <- metadata %>% 
  dplyr::filter(nchar(File.Name) != 0)

# Identify blank entries
metadata[metadata == "[Not Applicable]"] <- ""
metadata[metadata == "[Not Available]"] <- ""
metadata[metadata == "[Not Evaluated]"] <- ""
metadata[metadata == "[Not evaluated]"] <- ""
metadata[metadata == "[Unknown]"] <- ""
metadata[is.na(metadata)] <- ""

# Format metadata
metadata <-metadata %>%
  dplyr::rename(Sample_ID     = case_submitter_id,
                Project_ID    = project_id,
                Case.ID       = case_id,
                Sex           = gender,
                Race          = race,
                Ethnicity     = ethnicity,
                pS            = pathologic_stage,
                pT            = pathologic_T,
                pN            = pathologic_N,
                pM            = pathologic_M) %>%
  dplyr::mutate(Sample_ID = make.names(Sample_ID), # do not make them unique here
                Sex = str_to_sentence(Sex),
                Race = str_to_sentence(Race),
                Ethnicity  = str_to_sentence(Ethnicity),
                vital_status  = str_to_sentence(vital_status),
                Time = round(as.numeric(OS.time)/30,2),
                Status = as.numeric(OS), 
                OS = as.numeric(OS),
                OS.time = as.numeric(OS.time),
                DSS = as.numeric(DSS),
                DSS.time = as.numeric(DSS.time),
                DFI = as.numeric(DFI),
                DFI.time = as.numeric(DFI.time),
                PFI = as.numeric(PFI),
                PFI.time = as.numeric(PFI.time)) %>%
  dplyr::select(Sample_ID, Project_ID, Sample.Type, Time, Status, OS, OS.time, 
                DSS, DSS.time, DFI, DFI.time, PFI, PFI.time, Sex, Race, 
                Ethnicity, pS, pT, pN, pM, tumor_tissue_site, vital_status, 
                last_contact_days_to,	death_days_to,	new_tumor_event_dx_days_to,
                new_tumor_event_after_initial_treatment,
                treatment_outcome_first_course,
                primary_therapy_outcome_success,
                history_of_neoadjuvant_treatment,	radiation_therapy,	prior_dx,
                File.ID, File.Name, Sample.ID, Case.ID, Project.ID) %>%
  dplyr::arrange(Project_ID, Sample_ID)

# Create 2 dataframe : one for normal samples and other for tumors
# NOTE: While downloading from GDC portal, you might have noticed there are ONLY
# 10511 cases but 11499 files. This is because multiple samples from certain
# patients have been sequenced. 
tumor_df <- metadata %>% 
  dplyr::filter(Sample.Type != "Solid Tissue Normal") %>%
  dplyr::mutate(Sample_ID = make.names(Sample_ID, unique=TRUE))
normal_df <- metadata %>% 
  dplyr::filter(Sample.Type == "Solid Tissue Normal") %>%
  dplyr::mutate(Sample_ID = make.names(Sample_ID, unique=TRUE))

# Save excel
wb <- openxlsx::createWorkbook()
openxlsx::addWorksheet(wb, sheetName="Tumor")
openxlsx::writeData(wb, sheet="Tumor", x=tumor_df)
openxlsx::addWorksheet(wb, sheetName="Normal")
openxlsx::writeData(wb, sheet="Normal", x=normal_df)
openxlsx::saveWorkbook(wb, file=paste0(output_path, "TCGA.PanAtlas.Metadata.xlsx"), overwrite=TRUE)

#******************************************************************************#
#           STEP 3 (INDIVIDUAL CANCER): COMBINE TCGA RAW COUNTS                #
#******************************************************************************#

source("/hpc/home/kailasamms/projects/RNASeq/RNASeq_DESeq2_Functions.R")

file_path <- "/hpc/home/kailasamms/scratch/ICB_TCGA_datasets/original/"
data_path <- "/hpc/home/kailasamms/scratch/ICB_TCGA_datasets/original/raw_counts/"
output_path <- "/hpc/home/kailasamms/scratch/ICB_TCGA_datasets/processed/"

# Read the sample sheet
sample_sheet <- utils::read.table(file= paste0(file_path, "gdc_sample_sheet.2025-03-06.tsv"), 
                                  header=TRUE, sep="\t", quote="", skip=0, fill=TRUE)
colnames(sample_sheet) <- make.names(colnames(sample_sheet))

# Confirm that all files have been downloaded
files <- list.files(data_path)
sum(sort(sample_sheet$File.Name) == sort(files)) == length(files)

# Read metadata for tumor samples ONLY
metadata <- read.xlsx(paste0(output_path, "TCGA.PanAtlas.Metadata.xlsx"),
                      sheet = "Tumor")

# Read count files for each Project_ID and save as xlsx file
# NOTE: We use the value in "unstranded" column i.e. column 4 as counts since
# the kits used for library preparation may or may not be strand specific.
# https://www.biostars.org/p/9536035/
for (proj in unique(metadata$Project_ID)){
  
  # Get list of tsv files for each Project_ID i.e. each cancer from Sample sheet
  files <- metadata %>% 
    dplyr::filter(Project_ID == proj) %>% 
    dplyr::select(File.Name) %>%
    unlist(use.names=FALSE)
  
  # Get list of tsv files we have in folder
  tsv_files <- list.files(data_path)
  
  # Get list of tsv files we actually have for each cancer
  files <- intersect(files, tsv_files)
  
  # Load one of the files to extract gene ids and names
  temp_file <- utils::read.table(file=paste0(data_path, files[1]),
                                 header=TRUE, sep="\t", quote="", skip=0, fill=TRUE)
  
  # Generate a dummy dataframe to store raw counts
  counts <- data.frame(ENSEMBL_ID=temp_file[,1],
                       SYMBOL=temp_file[,2])
  
  # Generate a dummy dataframe to store tpm values
  tpm <- data.frame(ENSEMBL_ID=temp_file[,1],
                    SYMBOL=temp_file[,2])
  
  # Read each tsv file one by one
  for (i in 1:length(files)){
    
    # Keep track of file being processed
    cat(proj, ":", i, ":", files[i], "\n")
    
    # Reformat the tsv file
    # Keep ONLY columns that are needed and rename columns 
    temp_data <- utils::read.table(file=paste0(data_path, files[i]),
                                   header=TRUE, sep="\t", quote="", skip=0, fill=TRUE) %>%
      dplyr::select(gene_id, unstranded, tpm_unstranded) %>%
      dplyr::rename(ENSEMBL_ID = gene_id)
    
    # Append raw counts to count dataframe
    counts <- counts %>% 
      dplyr::left_join(temp_data %>% dplyr::select(ENSEMBL_ID, unstranded), 
                       by=c("ENSEMBL_ID"="ENSEMBL_ID")) %>%
      dplyr::rename(!!rlang::sym(files[i]) := unstranded)
    
    # Append tpm values to tpm dataframe
    tpm <- tpm %>% 
      dplyr::left_join(temp_data %>% dplyr::select(ENSEMBL_ID, tpm_unstranded), 
                       by=c("ENSEMBL_ID"="ENSEMBL_ID")) %>%
      dplyr::rename(!!rlang::sym(files[i]) := tpm_unstranded)
  } 
  
  # Remove rows containing counts for N_unmapped, N_multimapping, N_noFeature, etc
  # Remove version number from Ensembl_IDs
  counts <- counts %>%
    dplyr::filter(!(ENSEMBL_ID %in% c("N_unmapped", "N_multimapping", "N_noFeature", "N_ambiguous"))) %>%
    dplyr::mutate(ENSEMBL_ID = base::gsub(pattern=".[0-9]+$", replacement="", x=ENSEMBL_ID))
  
  tpm <- tpm %>%
    dplyr::filter(!(ENSEMBL_ID %in% c("N_unmapped", "N_multimapping", "N_noFeature", "N_ambiguous"))) %>%
    dplyr::mutate(ENSEMBL_ID = base::gsub(pattern=".[0-9]+$", replacement="", x=ENSEMBL_ID))
  
  # Rename column names to match with Sample_ID from metadata
  cols <- colnames(counts)
  
  # If exactly one match is found, rename column name with appropriate Case.ID
  for (i in 1:ncol(counts)){
    if(sum(stringr::str_detect(string=make.names(colnames(counts)[i]),
                               pattern=make.names(metadata$File.Name))) == 1){
      cols[i] <- metadata$Sample_ID[stringr::str_detect(string=make.names(colnames(counts)[i]),
                                                        pattern=make.names(metadata$File.Name))]
    }
    #   colnames(counts)[i] <- dplyr::if_else(colnames(counts)[i] %in% metadata$File.Name,
    #                                     metadata$Sample_ID[which(metadata$File.Name == colnames(counts)[i], arr.ind=TRUE)[1]],
    #                                     colnames(counts)[i])
  }
  
  colnames(counts) <- cols
  colnames(tpm) <- cols
  
  # Save raw counts & tpm values
  wb <- openxlsx::createWorkbook()
  openxlsx::addWorksheet(wb, sheetName="raw_counts")
  openxlsx::writeData(wb, sheet="raw_counts", x=counts)
  openxlsx::addWorksheet(wb, sheetName="tpm values")
  openxlsx::writeData(wb, sheet="tpm values", x=tpm)
  openxlsx::saveWorkbook(wb, file=paste0(output_path, make.names(proj), ".raw.counts.xlsx"), 
                         overwrite=TRUE)
}

#******************************************************************************#
#        STEP 4 (INDIVIDUAL CANCER): CALCULATE DESEQ2 NORMALIZED COUNTS        #
#******************************************************************************#

source("/hpc/home/kailasamms/projects/RNASeq/RNASeq_DESeq2_Functions.R")

file_path <- "/hpc/home/kailasamms/scratch/ICB_TCGA_datasets/original/"
data_path <- "/hpc/home/kailasamms/scratch/ICB_TCGA_datasets/processed/"
output_path <- "/hpc/home/kailasamms/scratch/ICB_TCGA_datasets/processed/"

meta_data_full <- openxlsx::read.xlsx(xlsxFile=paste0(data_path, "TCGA.PanAtlas.Metadata.xlsx"))
species <- "Homo sapiens"
annotations <- get_annotations(species)

for(proj in unique(meta_data_full$Project_ID)){
  
  # Read the raw counts
  read_data <- openxlsx::read.xlsx(xlsxFile=paste0(data_path, make.names(proj), ".raw.counts.xlsx"),
                                   sheet = "raw_counts") %>% 
    dplyr::mutate(SYMBOL = ENSEMBL_ID) %>%
    dplyr::select(everything(), -c("ENSEMBL_ID"))
  
  # Extract cancer specific metadata
  meta_data <- meta_data_full %>% 
    dplyr::filter(make.names(Sample_ID) %in% make.names(colnames(read_data)))
  
  # Perform QC
  Comparisons <- list(Variable =c(NA),
                      Target   =c(NA),
                      Reference=c(NA))
  
  meta_data <- prep_metadata(meta_data, read_data)
  read_data <- prep_readdata(read_data, meta_data)
  l <- check_data(read_data, meta_data)
  meta_data <- l[[2]]
  read_data <- l[[1]]
  
  # Prepare DESeq2 object
  dds <- DESeq2::DESeqDataSetFromMatrix(countData=read_data,
                                        colData=meta_data, 
                                        design=~ 1)
  
  # Calculate normalized counts using DESeq2()
  norm_counts <- deseq2_norm_counts(dds, meta_data, annotations)
  
  # Save batch corrected normalized counts
  wb <- openxlsx::createWorkbook()
  openxlsx::addWorksheet(wb, sheetName="Normalized counts")
  openxlsx::writeData(wb, sheet="Normalized counts", x=norm_counts)
  openxlsx::saveWorkbook(wb, file=paste0(output_path, make.names(proj), ".Normalized.counts.xlsx"), 
                         overwrite=TRUE)
  
  # Save the cancer specific meta_data
  wb <- openxlsx::createWorkbook()
  openxlsx::addWorksheet(wb, sheetName="Metadata")
  openxlsx::writeData(wb, sheet="Metadata", x=meta_data %>% tibble::rownames_to_column("Sample_ID"))
  openxlsx::saveWorkbook(wb, file=paste0(output_path, make.names(proj), ".Metadata.xlsx"), 
                         overwrite=TRUE)
}

#******************************************************************************#
#     STEP 5 (PAN CANCER): DOWNLOAD TCGA CLINICAL DATA & NORMALIZED COUNTS     #
#******************************************************************************#

# Download batch corrected normalized pan cancer RNA expression data from below
# https://gdc.cancer.gov/about-data/publications/pancanatlas

# TCGA-Clinical Data Resource (CDR) Outcome* - TCGA-CDR-SupplementalTableS1.xlsx
# Clinical with Follow-up - clinical_PANCAN_patient_with_followup.tsv
# RNA (Final) - EBPlusPlusAdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.tsv

# Check below link to understand clinical data columns
# https://gdc.cancer.gov/about-data/gdc-data-processing/clinical-data-standardization

#******************************************************************************#
#               STEP 6 (PAN CANCER): REFORMAT TCGA CLINICAL DATA               #
#******************************************************************************#

# Already generated in Step 1

#******************************************************************************#
#               STEP 7 (PAN CANCER): REFORMAT TCGA NORMALIZED COUNTS           #
#******************************************************************************#

source("/hpc/home/kailasamms/projects/RNASeq/RNASeq_DESeq2_Functions.R")

file_path <- "/hpc/home/kailasamms/scratch/ICB_TCGA_datasets/original/"
data_path <- "/hpc/home/kailasamms/scratch/ICB_TCGA_datasets/processed/"
output_path <- "/hpc/home/kailasamms/scratch/ICB_TCGA_datasets/processed/"

species <- "Homo sapiens"
annotations <- get_annotations(species)

# Read count data
norm_counts <- utils::read.table(file = paste0(file_path, "EBPlusPlusAdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.tsv"),
                                 header=TRUE, sep="\t", quote="", skip=0, fill=TRUE)

# Read metadata we made in Step 1
metadata <- openxlsx::read.xlsx(xlsxFile=paste0(data_path, "TCGA.PanAtlas.Metadata.xlsx"),
                                sheet = "Tumor")

# Rename column names  of norm_counts to match with Sample_ID from metadata
cols <- colnames(norm_counts)

# NOTE: str_detect is better than grepl
# grepl(pattern, replacement) where both pattern and replacement MUST be string
# str_detect(string, pattern) either string or pattern can be a vector but not BOTH.
for (i in 1:ncol(norm_counts)){
  
  if(i %% 500 == 0){
    cat(i, "samples processed\n")
  }
  
  if(sum(stringr::str_detect(string=colnames(norm_counts)[i],
                             pattern=make.names(metadata$Sample_ID))) == 1){
    cols[i] <- metadata$Sample_ID[stringr::str_detect(string=colnames(norm_counts)[i],
                                                      pattern=make.names(metadata$Sample_ID))]
  }
}

# Replace with metadata$Sample_ID as column names
colnames(norm_counts) <- cols
colnames(norm_counts)[1] <- "ENTREZ_ID"

# Keep only samples that have info in metadata
norm_counts <- norm_counts %>% 
  dplyr::select(ENTREZ_ID, intersect(colnames(norm_counts), metadata$Sample_ID))

# Separate the first column which contains SYMBOL and ENTREZ_ID joined together
# using piping character | into 2 columns
norm_counts <- norm_counts %>% 
  tidyr::separate(col=ENTREZ_ID, into=c("SYMBOL", "ENTREZ_ID"), sep="\\s*\\|\\s*")

# Remove the \" from SYMBOL and ENTREZ_ID
# NOTE: cat(norm_counts$ENTREZ_ID[100]) will show 5826"
# NOTE: print(norm_counts$ENTREZ_ID[100]) will show "5826\""
norm_counts$SYMBOL <- gsub(pattern="[\"\"]", replacement="", x=norm_counts$SYMBOL)
norm_counts$ENTREZ_ID <- gsub(pattern="[\"\"]", replacement="", x=norm_counts$ENTREZ_ID)

# Replace all NA values with 0
norm_counts <- norm_counts %>%
  base::replace(is.na(.), 0)

# Remove genes with no expression in all samples
norm_counts <- norm_counts[rowSums(norm_counts[,c(-1,-2)]) !=0,]

# Replace ? in SYMBOL with ENTREZ_ID
norm_counts <- norm_counts %>%
  dplyr::mutate(SYMBOL = dplyr::case_when(SYMBOL == "?" ~ ENTREZ_ID,
                                          TRUE ~ SYMBOL))

# Get latest gene symbols for the Entrez ids
normalized_counts <- add_annotation(norm_counts %>% dplyr::select(everything(), -c("SYMBOL")) %>% dplyr::rename(ID = ENTREZ_ID), 
                                    annotations)

# Save the reformatted norm_counts(DO NOT SAVE xlsx as file is too large)
write.table(x=normalized_counts, file=paste0(output_path, "TCGA.PanAtlas.Normalized.counts.tsv"), 
            quote=FALSE, sep='\t', row.names = FALSE)

#**********************************************