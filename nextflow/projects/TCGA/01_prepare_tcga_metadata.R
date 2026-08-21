# =============================================================================
# 01_prepare_tcga_metadata.R
# Merges GDC manifest, sample sheet, CDR survival data, and clinical data
# into clean per-cancer and pan-cancer metadata files.
#
# Run order: 01 → 02 → 03
# =============================================================================

library(dplyr)
library(stringr)
library(readr)
library(openxlsx)

# =============================================================================
# INPUTS / OUTPUTS
# =============================================================================

BASE_DIR         <- "/home/kailasamms/scripts/nextflow/resources/"

# Input files
MANIFEST_TXT     <- file.path(BASE_DIR, "Datasets", "gdc_manifest.2026-05-01.txt")
SAMPLE_SHEET_TSV <- file.path(BASE_DIR, "Datasets", "gdc_sample_sheet.2026-05-01.tsv")
SURVIVAL_XLSX    <- file.path(BASE_DIR, "Datasets", "TCGA-CDR-SupplementalTableS1.xlsx")
CLINICAL_TSV     <- file.path(BASE_DIR, "Datasets", "clinical.tsv")

# Output directory
METADATA_DIR     <- file.path(BASE_DIR, "Datasets", "01.Metadata")

dir.create(METADATA_DIR, recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# STEP 1: READ INPUT FILES
# =============================================================================

# Values treated as NA across all GDC files — GDC uses many placeholder strings
NA_STRINGS <- c("'--", "--", "[Not Available]", "[Not Applicable]",
                "[Not Evaluated]", "[Unknown]", "#N/A",
                "Not Reported", "not reported", "Unknown", "")

cat("Reading manifest...\n")
manifest <- readr::read_tsv(MANIFEST_TXT, show_col_types = FALSE)
cat("  Rows:", nrow(manifest), "\n")

cat("Reading sample sheet...\n")
sample_sheet <- readr::read_tsv(SAMPLE_SHEET_TSV, show_col_types = FALSE) %>%
  dplyr::rename(
    File_ID          = `File ID`,
    File_Name        = `File Name`,
    Project_ID       = `Project ID`,
    Case_ID          = `Case ID`,
    Sample_ID        = `Sample ID`,
    Tissue_Type      = `Tissue Type`,
    Tumor_Descriptor = `Tumor Descriptor`,
    Specimen_Type    = `Specimen Type`,
    Preservation     = `Preservation Method`
  )
cat("  Rows:", nrow(sample_sheet), "\n")

cat("Reading CDR survival data...\n")
# CDR (Liu et al. 2018 Cell) is the preferred source for survival endpoints —
# it is peer-reviewed and cleaner than GDC's own clinical.tsv survival columns.
# clinical.tsv is used ONLY for pT/pN/pM/grade/treatment (not in CDR).
cdr <- openxlsx::read.xlsx(
  xlsxFile = SURVIVAL_XLSX,
  sheet    = "TCGA-CDR",
  startRow = 1,
  colNames = TRUE
) %>%
  dplyr::select(-1) %>%
  dplyr::mutate(dplyr::across(dplyr::everything(),
                              ~ifelse(. %in% NA_STRINGS, NA_character_,
                                      as.character(.)))) %>%
  dplyr::select(-Redaction) %>%
  dplyr::rename(
    Case_ID   = bcr_patient_barcode,
    Cancer    = type,
    Age       = age_at_initial_pathologic_diagnosis,
    Sex       = gender,
    Race      = race,
    Stage     = ajcc_pathologic_tumor_stage,
    Histology = histological_type,
    Grade_CDR = histological_grade
  ) %>%
  dplyr::mutate(
    Sex       = stringr::str_to_sentence(Sex),
    Race      = stringr::str_to_sentence(Race),
    Stage     = stringr::str_to_sentence(Stage),
    Histology = stringr::str_to_sentence(Histology),
    OS        = as.numeric(OS),
    OS.time   = as.numeric(OS.time),
    DSS       = as.numeric(DSS),
    DSS.time  = as.numeric(DSS.time),
    DFI       = as.numeric(DFI),
    DFI.time  = as.numeric(DFI.time),
    PFI       = as.numeric(PFI),
    PFI.time  = as.numeric(PFI.time),
    # Convert days → months for easier clinical interpretation
    OS.months  = round(as.numeric(OS.time)  / 30.44, 2),
    DSS.months = round(as.numeric(DSS.time) / 30.44, 2),
    DFI.months = round(as.numeric(DFI.time) / 30.44, 2),
    PFI.months = round(as.numeric(PFI.time) / 30.44, 2)
  )
cat("  Rows:", nrow(cdr), "\n")

cat("Reading clinical.tsv...\n")
clinical <- readr::read_tsv(CLINICAL_TSV, show_col_types = FALSE) %>%
  dplyr::mutate(dplyr::across(dplyr::everything(),
                              ~ifelse(. %in% NA_STRINGS, NA_character_,
                                      as.character(.)))) %>%
  # Drop internal GDC UUID columns — not useful for analysis
  dplyr::select(-dplyr::any_of(c("cases.case_id",
                                 "demographic.submitter_id",
                                 "demographic.demographic_id",
                                 "diagnoses.submitter_id",
                                 "diagnoses.diagnosis_id",
                                 "diagnoses.tumor_of_origin"))) %>%
  # Drop columns >90% missing — not informative at pan-cancer scale
  dplyr::select(dplyr::where(~mean(is.na(.)) < 0.90)) %>%
  # Drop columns already covered by CDR — avoid duplication
  dplyr::select(-dplyr::any_of(c("project.project_id",
                                 "cases.consent_type",
                                 "cases.days_to_consent",
                                 "cases.lost_to_followup",
                                 "cases.index_date",
                                 "demographic.gender",
                                 "demographic.race",
                                 "demographic.vital_status",
                                 "demographic.days_to_death",
                                 "demographic.days_to_birth",
                                 "demographic.age_at_index",
                                 "demographic.age_is_obfuscated",
                                 "diagnoses.ajcc_pathologic_stage",
                                 "diagnoses.age_at_diagnosis",
                                 "diagnoses.days_to_last_follow_up")))
cat("  Rows:", nrow(clinical), "| Cols:", ncol(clinical), "\n")

# =============================================================================
# STEP 2: DEDUPLICATE SAMPLE SHEET — KEEP LARGEST FILE PER SAMPLE
# =============================================================================

# Some samples have multiple entries in the sample sheet (e.g. reprocessed files).
# Keep the largest file per sample — larger file = more complete quantification.
cat("Deduplicating sample sheet...\n")
sample_sheet_dedup <- sample_sheet %>%
  dplyr::left_join(manifest %>% dplyr::select(id, size),
                   by = c("File_ID" = "id")) %>%
  dplyr::group_by(Sample_ID) %>%
  dplyr::slice_max(size, n = 1, with_ties = FALSE) %>%
  dplyr::ungroup() %>%
  dplyr::select(-size)

cat("  Before:", nrow(sample_sheet), "| After:", nrow(sample_sheet_dedup), "\n")

# =============================================================================
# STEP 3: MATCH CLINICAL ROWS TO RNA-SEQ SAMPLES
# =============================================================================

cat("Matching clinical rows to samples...\n")

# Map Tumor_Descriptor → classification_of_tumor for row-selection priority
tumor_map <- c(
  "Primary"        = "primary",
  "Metastatic"     = "metastasis",
  "Recurrence"     = "recurrence",
  "New Primary"    = "primary",
  "Not Applicable" = NA_character_   # Normal samples have no tumor classification
)

# Collapse treatments per patient — clinical.tsv has one row per treatment,
# which would multiply rows during join. Collapse to one row per patient first.
treatments_collapsed <- clinical %>%
  dplyr::group_by(cases.submitter_id) %>%
  dplyr::summarise(
    Therapeutic_Agents = paste(unique(na.omit(treatments.therapeutic_agents)),
                               collapse = " | "),
    Treatment_Type     = paste(unique(na.omit(treatments.treatment_type)),
                               collapse = " | "),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    Therapeutic_Agents = dplyr::na_if(Therapeutic_Agents, ""),
    Treatment_Type     = dplyr::na_if(Treatment_Type, "")
  )

clinical_no_treat <- clinical %>%
  dplyr::select(-dplyr::starts_with("treatments."))

# Join sample sheet → clinical. One patient can have multiple diagnosis rows
# (e.g. primary + recurrence). Pick the best-matching row per sample using
# a priority hierarchy: matching tumor descriptor > primary disease > most data.
clinical_matched <- sample_sheet_dedup %>%
  dplyr::left_join(clinical_no_treat,
                   by = c("Case_ID" = "cases.submitter_id")) %>%
  dplyr::mutate(
    expected_classification = dplyr::recode(Tumor_Descriptor, !!!tumor_map)
  ) %>%
  dplyr::group_by(File_ID) %>%
  dplyr::arrange(
    dplyr::desc(
      !is.na(expected_classification) &
        !is.na(diagnoses.classification_of_tumor) &
        stringr::str_to_lower(diagnoses.classification_of_tumor) ==
        stringr::str_to_lower(expected_classification)
    ),
    dplyr::desc(diagnoses.diagnosis_is_primary_disease == "true"),
    dplyr::desc(rowSums(!is.na(dplyr::across(dplyr::everything())))),
    .by_group = TRUE
  ) %>%
  dplyr::slice(1) %>%
  dplyr::ungroup() %>%
  dplyr::select(-expected_classification,
                -diagnoses.diagnosis_is_primary_disease) %>%
  dplyr::left_join(treatments_collapsed,
                   by = c("Case_ID" = "cases.submitter_id"))

cat("  Matched rows:", nrow(clinical_matched), "\n")

# =============================================================================
# STEP 4: JOIN CDR SURVIVAL DATA
# =============================================================================

cat("Joining CDR...\n")
metadata <- clinical_matched %>%
  dplyr::left_join(cdr, by = "Case_ID") %>%
  dplyr::filter(!is.na(Sample_ID))

cat("  Rows after CDR join:", nrow(metadata), "\n")

# =============================================================================
# STEP 5: RENAME, ORDER, AND CLEAN COLUMNS
# =============================================================================

metadata <- metadata %>%
  dplyr::rename(
    pT                   = diagnoses.ajcc_pathologic_t,
    pN                   = diagnoses.ajcc_pathologic_n,
    pM                   = diagnoses.ajcc_pathologic_m,
    Tumor_Grade          = diagnoses.tumor_grade,
    Prior_Malig          = diagnoses.prior_malignancy,
    Prior_Treat          = diagnoses.prior_treatment,
    FIGO_Stage           = diagnoses.figo_stage,
    Morphology           = diagnoses.morphology,
    Primary_Diagnosis    = diagnoses.primary_diagnosis,
    Classification_Tumor = diagnoses.classification_of_tumor,
    Primary_Site         = cases.primary_site,
    Ethnicity            = demographic.ethnicity,
    Sex_At_Birth         = demographic.sex_at_birth
  ) %>%
  dplyr::mutate(
    Classification_Tumor = stringr::str_to_sentence(Classification_Tumor),
    # Flag samples missing from CDR — 40 samples lack CDR entries.
    # These still get expression data but have no survival endpoints.
    CDR_available = dplyr::if_else(is.na(Cancer), FALSE, TRUE),
    # Fill Cancer from Project_ID for the ~40 samples not in CDR
    Cancer        = dplyr::coalesce(Cancer,
                                    stringr::str_remove(Project_ID, "TCGA-"))
  ) %>%
  dplyr::select(
    Case_ID, Cancer, Project_ID, CDR_available,
    File_ID, File_Name,
    Sample_ID, Tissue_Type, Tumor_Descriptor,
    Specimen_Type, Preservation,
    Age, Sex, Sex_At_Birth, Race, Ethnicity,
    Stage, FIGO_Stage, pT, pN, pM,
    Histology, Grade_CDR, Tumor_Grade,
    Morphology, Primary_Diagnosis,
    Classification_Tumor, Primary_Site,
    Prior_Malig, Prior_Treat,
    Therapeutic_Agents, Treatment_Type,
    vital_status, tumor_status,
    OS, OS.time, OS.months,
    DSS, DSS.time, DSS.months,
    DFI, DFI.time, DFI.months,
    PFI, PFI.time, PFI.months,
    last_contact_days_to, death_days_to,
    cause_of_death, new_tumor_event_type,
    new_tumor_event_dx_days_to,
    treatment_outcome_first_course,
    margin_status, residual_tumor
  ) %>%
  dplyr::arrange(Project_ID, Case_ID)

cat("CDR available:", sum(metadata$CDR_available), "\n")
cat("CDR missing  :", sum(!metadata$CDR_available), "\n")
cat("Final metadata:", nrow(metadata), "rows x", ncol(metadata), "cols\n")

# =============================================================================
# STEP 6: SAVE PER-CANCER METADATA FILES
# =============================================================================

cat("Saving per-cancer metadata files...\n")

for (proj in sort(unique(metadata$Project_ID))) {

  proj_meta <- metadata %>% dplyr::filter(Project_ID == proj)

  out_file <- file.path(METADATA_DIR, paste0(proj, "_metadata.xlsx"))

  wb <- openxlsx::createWorkbook()
  openxlsx::addWorksheet(wb, sheetName = proj)
  openxlsx::writeData(wb, sheet = proj, x = proj_meta)
  openxlsx::saveWorkbook(wb, file = out_file, overwrite = TRUE)

  cat("  Saved:", basename(out_file), "\n")
}

# =============================================================================
# STEP 7: SAVE PAN-CANCER METADATA FILE
# =============================================================================

cat("Saving pan-cancer metadata file...\n")

out_file <- file.path(METADATA_DIR, "TCGA-PanCancer_metadata.xlsx")

wb_pan <- openxlsx::createWorkbook()
openxlsx::addWorksheet(wb_pan, sheetName = "PanCancer")
openxlsx::writeData(wb_pan, sheet = "PanCancer", x = metadata)
openxlsx::saveWorkbook(wb_pan, file = out_file, overwrite = TRUE)

cat("Done. Outputs in:", METADATA_DIR, "\n")
cat("  Per-cancer files:", length(unique(metadata$Project_ID)), "\n")
cat("  Pan-cancer file : TCGA-PanCancer_metadata.xlsx\n")
