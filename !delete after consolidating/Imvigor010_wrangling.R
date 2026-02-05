library(openxlsx)
library(dplyr)

sheets <- c("EGAD00001007574-signatera", "EGAD00001007576-20210423", 
            "EGAD00001007576-20210301", "EGAD00001011108-rnaseq", 
            "EGAD00001011108-signatera")


df1 <- openxlsx::read.xlsx(xlsxFile = "F:/Imvigor010_Metadata_merged.xlsx",
                           sheet = sheets[1])
df2 <- openxlsx::read.xlsx(xlsxFile = "F:/Imvigor010_Metadata_merged.xlsx",
                           sheet = sheets[2])
df3 <- openxlsx::read.xlsx(xlsxFile = "F:/Imvigor010_Metadata_merged.xlsx",
                           sheet = sheets[3])
df4 <- openxlsx::read.xlsx(xlsxFile = "F:/Imvigor010_Metadata_merged.xlsx",
                           sheet = sheets[4])
df5 <- openxlsx::read.xlsx(xlsxFile = "F:/Imvigor010_Metadata_merged.xlsx",
                           sheet = sheets[5])

# Generate ctDNA metadata
# Some patients have multiple samples: C1D1(Cycle 1,Day 1), C3D1(Cycle 3,Day 1)
all(sort(df1$anon_sampleID) == sort(df5$alias))
df_ctdna <- df5 %>% 
  dplyr::left_join(df1, by=c("alias" = "anon_sampleID")) %>%
  dplyr::mutate(gender = dplyr::case_when(gender == "male" ~ "Male",
                                          gender == "female" ~ "Female")) %>%
  dplyr::rename(Sample_id = subjectId,
                Sex = gender,
                Primary_disease = phenotype,
                PATIENT_ID = anon_patientID,
                FULL_PATIENT_ID = anon_full_patientID,
                CTDNA_SAMPLE_ID = alias,
                FULL_CTDNA_SAMPLE_ID = anon_full_sampleID) %>%
  dplyr::select(Sample_id, Sex, Primary_disease, FULL_PATIENT_ID, PATIENT_ID, 
                everything())

empty_cols <- c()
for (i in 1:ncol(df_ctdna)){
  if(all(is.na(df_ctdna[,i]))){
    empty_cols <- c(empty_cols, i)
  }
}
df_ctdna <- df_ctdna[, -empty_cols]

# Generate clinical metadata for patients (unique patient level data)
df3 <- df3 %>%
  dplyr::mutate(UNI_ID = paste0("PAT-", UNI_ID))

all(sort(df2$FULL_PATIENT_ID) == sort(df3$UNI_ID))
df_clinical <- df3 %>% 
  dplyr::left_join(df2, by=c("UNI_ID" = "FULL_PATIENT_ID")) %>%
  dplyr::mutate(Sample_id = PATIENT_ID,
                Time = OS_months*30,
                Status =dplyr::case_when(OS_event == "Alive" ~ 0,
                                         OS_event == "Death" ~ 1),
                SEX = dplyr::case_when(SEX == "M" ~ "Male",
                                       SEX == "F" ~ "Female")) %>%
  dplyr::rename(FULL_PATIENT_ID = UNI_ID,
                Stage = tumor_stage,
                Primary_disease = PRIDIS,
                Race = RACE,
                Sex = SEX) %>%
  dplyr::select(Sample_id, Sex, Time, Status, Stage, Race, Primary_disease,
                FULL_PATIENT_ID, PATIENT_ID, everything())

# Now, add RNAseq file names to clinical metadata so we can link each fastq file
# with a patient

empty_cols <- c()
for (i in 1:ncol(df4)){
  if(all(is.na(df4[,i]))){
    empty_cols <- c(empty_cols, i)
  }
}
df4 <- df4[, -empty_cols]

df_clinical <- df_clinical %>% 
  dplyr::left_join(df4, by=c("Sample_id" = "subjectId")) %>%
  dplyr::select(everything(), -c(gender, phenotype)) %>%
  dplyr::rename(RNA_filename = alias) %>%
  dplyr::select(Sample_id, Sex, Time, Status, Stage, Race, Primary_disease,
                RNA_filename, FULL_PATIENT_ID, PATIENT_ID, everything())

# Save batch corrected normalized counts for entire dataset
wb <- openxlsx::createWorkbook()
openxlsx::addWorksheet(wb, sheetName = "Clinical")
openxlsx::writeData(wb, sheet = "Clinical", x = df_clinical, rowNames = FALSE)
openxlsx::addWorksheet(wb, sheetName = "Ct_DNA")
openxlsx::writeData(wb, sheet = "Ct_DNA", x = df_ctdna, rowNames = FALSE)
openxlsx::saveWorkbook(wb, file = "F:/Imvigor010_Metadata.xlsx", overwrite = TRUE)

