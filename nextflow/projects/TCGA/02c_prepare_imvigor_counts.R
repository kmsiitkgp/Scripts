# =============================================================================
# 02c_prepare_imvigor_counts.R
#
# Purpose: Reformat IMvigor raw-count xlsx files into the same structure as
#          Script 02a/02b output (ID | SYMBOL | gene_type | sample columns),
#          and derive TPM from raw counts + GTF-based gene lengths, since
#          IMvigor was supplied as raw counts only (no accompanying TPM/FPKM
#          file, unlike TCGA [GDC-computed] or GEO [depositor-supplied]).
#
# TPM note: computed using union-of-exons gene length (see
#           00_create_gene_lengths.R). This is a standard approximation but
#           is NOT methodologically identical to Salmon/GDC effective-length
#           TPM. Expect IMvigor TPM to differ somewhat from TCGA/GEO TPM for
#           the same gene, in a way that is not biological.
#
# Inputs:
#   IMvigor010_Raw_Counts.xlsx / IMvigor210_Raw_Counts.xlsx
#     — columns: gene_id, <sample_1>, <sample_2>, ...  (raw integer counts;
#       "gene_id" column already fixed by hand to contain ENSG IDs, not the
#       gene_name values the original "SYMBOL" header implied)
#   tx2gene_GRCh38.115.csv   — transcript_id, gene_id, gene_name, gene_biotype
#   gene_lengths_GRCh38.115.csv — gene_id, length  (from Script 00)
#
# Outputs (same schema/naming convention as 02a/02b):
#   IMvigor010_Raw_Counts.tsv / .rds
#   IMvigor010_TPM_Counts.tsv / .rds
#   IMvigor210_Raw_Counts.tsv / .rds
#   IMvigor210_TPM_Counts.tsv / .rds
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(openxlsx)
  library(readr)
})

# =============================================================================
# INPUTS / OUTPUTS
# =============================================================================

BASE_DIR    <- "/home/kailasamms/scripts/nextflow/resources/Datasets"

INPUT_DIR   <- file.path(BASE_DIR, "00.Downloads", "IMvigor")
COUNTS_DIR  <- file.path(BASE_DIR, "02.Counts")          # same output dir as 02a/02b

TX2GENE_PATH      <- file.path(BASE_DIR, "tx2gene_GRCh38.115.csv")
GENE_LENGTHS_PATH <- file.path(BASE_DIR, "gene_lengths_GRCh38.115.csv")

dir.create(COUNTS_DIR, recursive = TRUE, showWarnings = FALSE)

cat("INPUT_DIR          :", INPUT_DIR, "\n")
cat("COUNTS_DIR          :", COUNTS_DIR, "\n")
cat("TX2GENE_PATH        :", TX2GENE_PATH, "\n")
cat("GENE_LENGTHS_PATH   :", GENE_LENGTHS_PATH, "\n")

# =============================================================================
# STEP 1: Load reference tables (once, shared across both IMvigor files)
# =============================================================================

cat("\nLoading tx2gene and gene length reference tables...\n")

tx2gene <- readr::read_csv(TX2GENE_PATH, show_col_types = FALSE) %>%
  dplyr::select(gene_id, gene_name, gene_biotype) %>%
  dplyr::distinct()

gene_lengths <- readr::read_csv(GENE_LENGTHS_PATH, show_col_types = FALSE)

cat("  tx2gene: ", nrow(tx2gene), "unique genes\n")
cat("  gene_lengths: ", nrow(gene_lengths), "genes\n")

# =============================================================================
# STEP 2: FIND IMVIGOR RAW COUNT FILES
# =============================================================================

xlsx_files <- list.files(
  path       = INPUT_DIR,
  pattern    = "_Raw_Counts\\.xlsx$",
  full.names = TRUE
)

if (length(xlsx_files) == 0) {
  stop("No *_Raw_Counts.xlsx files found in: ", INPUT_DIR)
}

cat("\nFound", length(xlsx_files), "IMvigor file(s):",
    paste(basename(xlsx_files), collapse = ", "), "\n")

# =============================================================================
# STEP 3: PROCESS EACH FILE
# =============================================================================

for (xlsx_file in sort(xlsx_files)) {

  proj <- sub("_Raw_Counts\\.xlsx$", "", basename(xlsx_file))
  cat("\n", proj, "\n", sep = "")

  # ---------------------------------------------------------------------------
  # 3a. Read raw counts
  # ---------------------------------------------------------------------------
  counts_raw <- openxlsx::read.xlsx(xlsx_file, sheet = 1, colNames = TRUE)

  if (!"gene_id" %in% colnames(counts_raw)) {
    cat("  SKIPPING", proj, "— 'gene_id' column not found\n")
    next
  }

  sample_cols <- setdiff(colnames(counts_raw), "gene_id")
  cat("  Read", nrow(counts_raw), "genes x", length(sample_cols), "samples\n")

  # ---------------------------------------------------------------------------
  # 3b. Join annotation (gene_name, gene_biotype) via tx2gene
  # ---------------------------------------------------------------------------
  counts_annot <- counts_raw %>%
    dplyr::left_join(tx2gene, by = "gene_id") %>%
    dplyr::rename(ID = gene_id, SYMBOL = gene_name, gene_type = gene_biotype)

  n_no_symbol <- sum(is.na(counts_annot$SYMBOL) | counts_annot$SYMBOL == "")
  if (n_no_symbol > 0) {
    cat("  Dropping", n_no_symbol,
        "gene(s) with no resolvable SYMBOL (unannotated / non-coding, no name)\n")
  }
  counts_annot <- counts_annot %>%
    dplyr::filter(!is.na(SYMBOL), SYMBOL != "")

  # ---------------------------------------------------------------------------
  # 3c. Dedup on SYMBOL — highest total raw count wins
  # (consistent with 02a/02b dedup logic)
  # ---------------------------------------------------------------------------
  n_before <- nrow(counts_annot)

  counts_annot <- counts_annot %>%
    dplyr::mutate(
      row_sum = rowSums(dplyr::across(dplyr::all_of(sample_cols)), na.rm = TRUE)
    ) %>%
    dplyr::group_by(SYMBOL) %>%
    dplyr::slice_max(order_by = row_sum, n = 1, with_ties = FALSE) %>%
    dplyr::ungroup()

  n_after <- nrow(counts_annot)
  if (n_after < n_before) {
    cat("  Duplicate SYMBOL collapse:", n_before, "->", n_after,
        "rows (kept highest total-count copy)\n")
  }

  # ---------------------------------------------------------------------------
  # 3d. Remove all-zero genes
  # ---------------------------------------------------------------------------
  nonzero <- counts_annot$row_sum > 0
  n_zero  <- sum(!nonzero)
  if (n_zero > 0) {
    cat("  Removing", n_zero, "all-zero gene(s)\n")
  }
  counts_annot <- counts_annot[nonzero, ]

  counts_final <- counts_annot %>%
    dplyr::select(ID, SYMBOL, gene_type, dplyr::all_of(sample_cols))

  cat("  Final raw count matrix:", nrow(counts_final), "genes x",
      length(sample_cols), "samples\n")

  # ---------------------------------------------------------------------------
  # 3e. Compute TPM using GTF-derived gene lengths
  # ---------------------------------------------------------------------------
  tpm_input <- counts_final %>%
    dplyr::left_join(gene_lengths, by = c("ID" = "gene_id"))

  n_no_length <- sum(is.na(tpm_input$length))
  if (n_no_length > 0) {
    cat("  WARNING:", n_no_length,
        "gene(s) missing GTF length — dropped from both raw and TPM output ",
        "to keep gene sets identical\n")
  }
  tpm_input <- tpm_input %>% dplyr::filter(!is.na(length))

  # rate = count / length (kb), then normalize to sum to 1e6 per sample
  rate_mat <- as.matrix(tpm_input[, sample_cols]) / (tpm_input$length / 1000)
  tpm_mat  <- sweep(rate_mat, 2, colSums(rate_mat, na.rm = TRUE), FUN = "/") * 1e6

  tpm_final <- tpm_input %>%
    dplyr::select(ID, SYMBOL, gene_type)
  tpm_final <- cbind(tpm_final, as.data.frame(tpm_mat))

  # Keep raw counts restricted to the same gene set as TPM (identical IDs
  # across raw/TPM outputs, mirroring 02b's guarantee)
  counts_final <- counts_final %>% dplyr::filter(ID %in% tpm_final$ID)

  cat("  Final TPM matrix:", nrow(tpm_final), "genes x",
      length(sample_cols), "samples\n")

  # ---------------------------------------------------------------------------
  # 3f. Save
  # ---------------------------------------------------------------------------
  raw_tsv <- file.path(COUNTS_DIR, paste0(proj, "_Raw_Counts.tsv"))
  raw_rds <- file.path(COUNTS_DIR, paste0(proj, "_Raw_Counts.rds"))
  tpm_tsv <- file.path(COUNTS_DIR, paste0(proj, "_TPM_Counts.tsv"))
  tpm_rds <- file.path(COUNTS_DIR, paste0(proj, "_TPM_Counts.rds"))

  readr::write_tsv(counts_final, raw_tsv)
  readr::write_tsv(tpm_final,    tpm_tsv)
  saveRDS(counts_final, raw_rds)
  saveRDS(tpm_final,    tpm_rds)

  cat("  Saved:", basename(raw_tsv), "|", basename(tpm_tsv), "|",
      basename(raw_rds), "|", basename(tpm_rds), "\n")
}

cat("\nDone. Outputs in:", COUNTS_DIR, "\n")
