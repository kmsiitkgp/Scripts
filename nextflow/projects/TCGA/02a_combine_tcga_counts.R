# =============================================================================
# 02_combine_tcga_counts.R
# Reads per-sample GDC STAR count TSV files and combines them into
# per-cancer raw count and TPM matrices.
#
# Run order: 01 → 02 → 03
# =============================================================================

library(dplyr)
library(readr)
library(purrr)
library(openxlsx)

# =============================================================================
# INPUTS / OUTPUTS
# =============================================================================

BASE_DIR      <- "/home/kailasamms/scripts/nextflow/resources/"

# Input directories
DOWNLOADS_DIR <- file.path(BASE_DIR, "Datasets", "00.Downloads", "TCGA")  # GDC STAR TSV files
METADATA_DIR  <- file.path(BASE_DIR, "Datasets", "01.Metadata")   # output of Script 01

# Output directory
COUNTS_DIR    <- file.path(BASE_DIR, "Datasets", "02.Counts")

dir.create(COUNTS_DIR, recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# CONSTANTS
# =============================================================================

# STAR/GDC prepends 4 summary rows before gene-level counts — exclude these
SUMMARY_ROWS <- c("N_unmapped", "N_multimapping", "N_noFeature", "N_ambiguous")

# =============================================================================
# STEP 1: FIND PER-CANCER METADATA FILES
# =============================================================================

cat("Scanning metadata directory...\n")

metadata_files <- list.files(
  path       = METADATA_DIR,
  pattern    = "^TCGA-.*_metadata\\.xlsx$",
  full.names = TRUE
)

# Exclude pan-cancer file — process per cancer only
metadata_files <- metadata_files[!grepl("PanCancer", metadata_files)]

if (length(metadata_files) == 0) {
  stop("No per-cancer metadata files found in: ", METADATA_DIR)
}

cat("  Found", length(metadata_files), "cancer(s)\n")

# =============================================================================
# STEP 2: PROCESS EACH CANCER
# =============================================================================

for (meta_file in sort(metadata_files)) {

  proj <- sub("_metadata\\.xlsx$", "", basename(meta_file))
  cat("\n", proj, "\n", sep = "")

  # ---------------------------------------------------------------------------
  # 2a. Read metadata — extract File_Name → Sample_ID mapping
  # ---------------------------------------------------------------------------
  metadata <- openxlsx::read.xlsx(
    xlsxFile = meta_file,
    sheet    = proj,
    colNames = TRUE
  )

  if (!all(c("File_Name", "Sample_ID") %in% colnames(metadata))) {
    cat("  SKIPPING", proj, "— missing File_Name or Sample_ID\n")
    next
  }

  # One-to-one lookup: each GDC TSV filename maps to one Sample_ID
  file_to_sample <- metadata %>%
    dplyr::select(File_Name, Sample_ID) %>%
    dplyr::filter(!is.na(File_Name), !is.na(Sample_ID)) %>%
    dplyr::distinct(File_Name, .keep_all = TRUE)

  cat("  Samples in metadata:", nrow(file_to_sample), "\n")

  # ---------------------------------------------------------------------------
  # 2b. Find TSV files on disk for this cancer
  # ---------------------------------------------------------------------------
  tsv_files <- list.files(
    path       = DOWNLOADS_DIR,
    pattern    = "*_star_gene_counts\\.tsv$",
    full.names = TRUE
  )

  # Keep only files listed in this cancer's metadata
  tsv_files <- tsv_files[basename(tsv_files) %in% file_to_sample$File_Name]

  if (length(tsv_files) == 0) {
    cat("  SKIPPING", proj, "— no TSV files found on disk\n")
    next
  }

  cat("  TSV files on disk:", length(tsv_files), "\n")

  missing_files <- setdiff(file_to_sample$File_Name, basename(tsv_files))
  if (length(missing_files) > 0) {
    cat("  WARNING:", length(missing_files), "metadata files not found on disk\n")
  }

  # ---------------------------------------------------------------------------
  # 2c. Read each TSV and extract raw counts + TPM
  # ---------------------------------------------------------------------------
  # GDC STAR TSV format:
  #   Line 1 : comment (# gene-model: GENCODE v36) — skip with skip=1
  #   Line 2 : header (gene_id, gene_name, gene_type, unstranded, tpm_unstranded, ...)
  #   Lines 3-6 : STAR summary rows — excluded via SUMMARY_ROWS filter
  #   Lines 7+  : per-gene data
  #
  # Column choices:
  #   unstranded     → raw counts. TCGA used unstranded because samples came
  #                    from many sequencing centers with mixed library prep kits.
  #   tpm_unstranded → TPM. Used by BostonGene TME classifier (run_batch.py).

  cat("  Reading", length(tsv_files), "files...\n")

  counts_list <- vector("list", length(tsv_files))
  tpm_list    <- vector("list", length(tsv_files))

  for (i in seq_along(tsv_files)) {

    tsv_file  <- tsv_files[i]
    file_name <- basename(tsv_file)
    sample_id <- file_to_sample$Sample_ID[file_to_sample$File_Name == file_name]

    if (length(sample_id) == 0 || is.na(sample_id)) {
      cat("  WARNING: No Sample_ID for", file_name, "— skipping\n")
      next
    }

    df <- tryCatch({
      readr::read_tsv(tsv_file, skip = 1, col_names = TRUE, show_col_types = FALSE)
    }, error = function(e) {
      cat("  ERROR reading", file_name, ":", e$message, "\n")
      return(NULL)
    })

    if (is.null(df)) next

    df <- df %>%
      dplyr::filter(!(gene_id %in% SUMMARY_ROWS)) %>%
      # Strip Ensembl version suffix (ENSG00000000003.15 → ENSG00000000003)
      # Version suffix varies across GDC releases — strip for stable IDs
      dplyr::mutate(gene_id = gsub("\\.[0-9]+$", "", gene_id))

    counts_list[[i]] <- df %>%
      dplyr::select(ID = gene_id, !!sample_id := unstranded)

    tpm_list[[i]] <- df %>%
      dplyr::select(ID = gene_id, !!sample_id := tpm_unstranded)

    if (i %% 100 == 0) cat("    Processed", i, "/", length(tsv_files), "\n")
  }

  counts_list <- purrr::compact(counts_list)
  tpm_list    <- purrr::compact(tpm_list)

  if (length(counts_list) == 0) {
    cat("  SKIPPING", proj, "— no valid files after reading\n")
    next
  }

  # ---------------------------------------------------------------------------
  # 2d. Extract gene annotation from first file
  # All TCGA files use GENCODE v36 — gene_id/gene_name mapping is identical
  # across files so reading from the first file is sufficient.
  # ---------------------------------------------------------------------------
  gene_info <- readr::read_tsv(
    tsv_files[1], skip = 1, col_names = TRUE, show_col_types = FALSE
  ) %>%
    dplyr::filter(!(gene_id %in% SUMMARY_ROWS)) %>%
    dplyr::mutate(gene_id = gsub("\\.[0-9]+$", "", gene_id)) %>%
    dplyr::select(ID = gene_id, SYMBOL = gene_name, gene_type)

  # ---------------------------------------------------------------------------
  # 2e. Merge all samples into matrices
  # full_join preserves all genes even if absent in a sample (→ filled with 0)
  # ---------------------------------------------------------------------------
  cat("  Building matrices...\n")

  counts_mat <- purrr::reduce(counts_list, dplyr::full_join, by = "ID")
  tpm_mat    <- purrr::reduce(tpm_list,    dplyr::full_join, by = "ID")

  # Replace NAs from full_join with 0 — gene present in annotation, zero reads
  counts_mat[is.na(counts_mat)] <- 0
  tpm_mat[is.na(tpm_mat)]       <- 0

  # ---------------------------------------------------------------------------
  # 2f. Prepend ID, SYMBOL, gene_type columns
  # ---------------------------------------------------------------------------
  counts_mat <- gene_info %>%
    dplyr::left_join(counts_mat, by = "ID") %>%
    dplyr::select(ID, SYMBOL, gene_type, dplyr::everything())

  tpm_mat <- gene_info %>%
    dplyr::left_join(tpm_mat, by = "ID") %>%
    dplyr::select(ID, SYMBOL, gene_type, dplyr::everything())

  # ---------------------------------------------------------------------------
  # 2g. Remove all-zero genes (not expressed in any sample for this cancer)
  # ---------------------------------------------------------------------------
  nonzero <- rowSums(dplyr::select(counts_mat, -ID, -SYMBOL, -gene_type)) > 0
  counts_mat <- counts_mat[nonzero, ]
  tpm_mat    <- tpm_mat[nonzero, ]

  cat("  Genes retained:", nrow(counts_mat),
      "| Samples:", ncol(counts_mat) - 3, "\n")

  # ---------------------------------------------------------------------------
  # 2h. Save outputs — separate tsv per matrix type + RDS for fast R loading
  # Separate files (not sheets) so filename alone identifies content.
  # RDS for downstream R scripts (Script 03); xlsx for spot-checking.
  # ---------------------------------------------------------------------------
  raw_tsv <- file.path(COUNTS_DIR, paste0(proj, "_Raw_Counts.tsv"))
  tpm_tsv <- file.path(COUNTS_DIR, paste0(proj, "_TPM_Counts.tsv"))
  raw_rds <- file.path(COUNTS_DIR, paste0(proj, "_Raw_Counts.rds"))
  tpm_rds <- file.path(COUNTS_DIR, paste0(proj, "_TPM_Counts.rds"))

  # Save as TSV (tab-delimited)
  readr::write_tsv(counts_mat, raw_tsv)
  readr::write_tsv(tpm_mat,    tpm_tsv)

  # Save as RDS for internal pipeline use
  saveRDS(counts_mat, raw_rds)
  saveRDS(tpm_mat,    tpm_rds)

  cat("  Saved:", basename(raw_tsv), "|", basename(tpm_tsv),
      "|", basename(raw_rds), "|", basename(tpm_rds), "\n")
}

cat("\nDone. Outputs in:", COUNTS_DIR, "\n")
