# =============================================================================
# 04_prepare_geo_counts.R
# Reformats GEO count data (Entrez Gene IDs) into the same structure as
# Script 02 output (ID | SYMBOL | gene_type | sample columns), so it can
# feed directly into Script 03 (DESeq2) without any changes.
#
# Run order: 01 → 02 → 03 (TCGA)
#                  04 → 03 (GEO)
#
# INPUTS:
#   1. GEO raw counts file  : rows = Entrez IDs, cols = GSM sample IDs
#   2. GEO TPM file         : same structure as raw counts
#   3. tx2gene CSV          : output of create_tx2gene() — GTF-derived
#                             (transcript_id | gene_id | gene_name | gene_biotype)
#   4. Metadata XLSX        : manually created, must have Sample_ID column
#                             where Sample_ID matches GSM IDs in count files
#
# OUTPUTS (identical structure to Script 02):
#   GSE_Raw_Counts.tsv / .rds
#   GSE_TPM_Counts.tsv / .rds
#
# ENTREZ → SYMBOL MAPPING STRATEGY:
#   Entrez → ENSG         : org.Hs.eg.db (only source for Entrez → ENSG)
#   ENSG  → gene_name     : tx2gene CSV  (GTF-derived, latest Ensembl release)
#   ENSG  → gene_biotype  : tx2gene CSV  (GTF-derived, latest Ensembl release)
#
#   Multi-mapping (many Entrez → same SYMBOL):
#     Dedup is performed on raw counts only, using highest total expression
#     (row sum) per SYMBOL to pick the winning Entrez ID
#     TPM is then filtered to the same winning Entrez IDs — no independent
#     dedup — guaranteeing identical gene sets in both output matrices.
#
#   ENSG present in org.Hs.eg.db but absent from tx2gene:
#     → gene_name = NA → dropped (no reliable symbol from GTF)
#
#   ENSG present in tx2gene but absent from org.Hs.eg.db:
#     → won't appear in entrez map; GEO file won't have these anyway
#     since GEO uses Entrez IDs as row identifiers
# =============================================================================

library(dplyr)
library(readr)
library(tibble)
library(openxlsx)
library(org.Hs.eg.db)
library(AnnotationDbi)

# =============================================================================
# INPUTS / OUTPUTS
# =============================================================================

BASE_DIR     <- "/home/kailasamms/scripts/nextflow/resources/"

# Input directories
DOWNLOADS_DIR <- file.path(BASE_DIR, "Datasets", "00.Downloads", "GEO")  # GEO TSV files
METADATA_DIR  <- file.path(BASE_DIR, "Datasets", "01.Metadata")          # GEO metadata files

# Input files
TX2GENE_CSV   <- file.path(BASE_DIR, "Datasets", "tx2gene_GRCh38.115.csv")  # from create_tx2gene()

# Output directory
COUNTS_DIR    <- file.path(BASE_DIR, "Datasets", "02.Counts")

dir.create(COUNTS_DIR, recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# STEP 1: BUILD ENTREZ → ENSG MAP
# =============================================================================
# Entrez → ENSG from org.Hs.eg.db (only reliable source for this direction)
# ENSG  → gene_name + gene_biotype from tx2gene CSV (GTF-derived, more current)
# =============================================================================

cat("Building Entrez → ENSG map...\n")

# Pull all Entrez → ENSG mappings from org.Hs.eg.db
# One Entrez can map to multiple ENSG — kept here; resolved later per SYMBOL
entrez_to_ensg <- AnnotationDbi::select(
  x       = org.Hs.eg.db,
  keys    = AnnotationDbi::keys(org.Hs.eg.db, keytype = "ENTREZID"),
  columns = c("ENTREZID", "ENSEMBL"),
  keytype = "ENTREZID"
) %>%
  dplyr::rename(entrez_id = ENTREZID, gene_id = ENSEMBL) %>%
  dplyr::filter(!is.na(entrez_id), !is.na(gene_id)) %>%
  dplyr::mutate(entrez_id = as.character(entrez_id))

cat("  Entrez → ENSG pairs from org.Hs.eg.db:", nrow(entrez_to_ensg), "\n")

# Load tx2gene — collapse to gene level (one row per ENSG)
# select(2:4): col2 = gene_id, col3 = gene_name, col4 = gene_biotype
tx2gene <- readr::read_csv(TX2GENE_CSV, show_col_types = FALSE) %>%
  dplyr::select(gene_id = 2, gene_name = 3, gene_biotype = 4) %>%
  dplyr::distinct(gene_id, .keep_all = TRUE) %>%
  dplyr::filter(!is.na(gene_id))

cat("  Unique ENSG in tx2gene:", nrow(tx2gene), "\n")

# Join: Entrez → ENSG (org.Hs.eg.db) + ENSG → gene_name/biotype (tx2gene)
# Left join on gene_id — ENSGs absent from tx2gene get NA gene_name → dropped
entrez_map <- entrez_to_ensg %>%
  dplyr::left_join(tx2gene, by = "gene_id") %>%
  dplyr::filter(!is.na(gene_name))

cat("  Entrez IDs with resolvable ENSG + symbol:", dplyr::n_distinct(entrez_map$entrez_id), "\n")
cat("  ENSG IDs covered:", dplyr::n_distinct(entrez_map$gene_id), "\n")

# =============================================================================
# STEP 2: FIND PER-DATASET GEO FILES
# =============================================================================

cat("Scanning metadata directory...\n")

metadata_files <- list.files(
  path       = METADATA_DIR,
  pattern    = "^GSE[0-9].*_metadata\\.xlsx$",
  full.names = TRUE
)

if (length(metadata_files) == 0) {
  stop("No per-dataset metadata files found in: ", METADATA_DIR)
}

cat("  Found", length(metadata_files), "dataset(s)\n")

# =============================================================================
# STEP 3: PROCESS EACH GEO DATASET
# =============================================================================

for (meta_file in sort(metadata_files)) {

  proj <- sub("_metadata\\.xlsx$", "", basename(meta_file))
  cat("\n", proj, "\n", sep = "")

  # ---------------------------------------------------------------------------
  # 3a. Read metadata
  # ---------------------------------------------------------------------------
  metadata <- openxlsx::read.xlsx(
    xlsxFile = meta_file,
    sheet    = proj,
    colNames = TRUE
  )

  if (!"Sample_ID" %in% colnames(metadata)) {
    cat("  SKIPPING", proj, "— missing Sample_ID\n")
    next
  }

  samples <- metadata %>%
    dplyr::select(Sample_ID) %>%
    dplyr::filter(!is.na(Sample_ID)) %>%
    dplyr::distinct(Sample_ID, .keep_all = TRUE)

  cat("  Samples in metadata:", nrow(samples), "\n")

  # ---------------------------------------------------------------------------
  # 3b. Find TSV files on disk for this dataset
  # ---------------------------------------------------------------------------
  tsv_files <- list.files(
    path       = DOWNLOADS_DIR,
    pattern    = proj,
    full.names = TRUE
  )

  count_file <- grep(pattern = "_Raw_Counts", x = tsv_files, value = TRUE)
  tpm_file   <- grep(pattern = "_TPM_Counts", x = tsv_files, value = TRUE)

  if (length(count_file) == 0 || length(tpm_file) == 0) {
    cat("  SKIPPING", proj, "— raw or TPM file not found\n")
    next
  }

  # ---------------------------------------------------------------------------
  # 3c. Read raw counts and TPM files
  # GEO format: first column = Entrez Gene ID (GeneID), remaining = GSM IDs
  # ---------------------------------------------------------------------------
  cat("  Reading", basename(count_file), "...\n")
  counts_raw <- tryCatch({
    readr::read_tsv(count_file, col_names = TRUE, show_col_types = FALSE) %>%
      dplyr::mutate(GeneID = as.character(GeneID))
  }, error = function(e) {
    cat("  ERROR reading", count_file, ":", e$message, "\n"); return(NULL)
  })

  cat("  Reading", basename(tpm_file), "...\n")
  tpm_raw <- tryCatch({
    readr::read_tsv(tpm_file, col_names = TRUE, show_col_types = FALSE) %>%
      dplyr::mutate(GeneID = as.character(GeneID))
  }, error = function(e) {
    cat("  ERROR reading", tpm_file, ":", e$message, "\n"); return(NULL)
  })

  if (is.null(counts_raw) || is.null(tpm_raw)) next

  # Replace any NAs with 0
  counts_raw[is.na(counts_raw)] <- 0
  tpm_raw[is.na(tpm_raw)]       <- 0

  # ---------------------------------------------------------------------------
  # 3d. Map Entrez → annotation and deduplicate on SYMBOL using raw counts
  #
  # Strategy:
  #   1. Join raw counts to entrez_map → expands to ENSG + SYMBOL + gene_type
  #   2. Compute row_mean across sample columns (raw counts only)
  #   3. group_by(SYMBOL) → slice_max(row_mean) → one row per SYMBOL
  #      ENSG (ID) and gene_type are kept as passengers from the winning row
  #   4. Extract winning entrez_ids from deduped raw counts
  #   5. Filter TPM to those same entrez_ids, then join annotation from
  #      deduped raw counts — no independent dedup on TPM, guarantees
  #      identical gene sets in both output matrices
  # ---------------------------------------------------------------------------

  sample_cols <- setdiff(colnames(counts_raw), "GeneID")

  # Step 1: join raw counts to annotation
  counts_mat <- counts_raw %>%
    dplyr::inner_join(entrez_map, by = c("GeneID" = "entrez_id"))

  cat("  Rows after Entrez → annotation join:", nrow(counts_mat),
      "(from", nrow(counts_raw), "Entrez rows)\n")

  # Step 2: compute row means on sample columns
  counts_mat <- counts_mat %>%
    dplyr::mutate(
      row_sum = rowSums(dplyr::across(dplyr::all_of(sample_cols)), na.rm = TRUE)
    )

  # Step 3: dedup on SYMBOL — highest mean raw count wins
  counts_mat <- counts_mat %>%
    dplyr::group_by(gene_name) %>%
    dplyr::slice_max(order_by = row_sum, n = 1, with_ties = FALSE) %>%
    dplyr::ungroup()

  cat("  Genes retained after SYMBOL dedup:", nrow(counts_mat), "\n")

  # Step 4: extract winning entrez_ids and annotation from deduped raw counts
  keep_entrez <- counts_mat$GeneID

  annotation <- counts_mat %>%
    dplyr::select(GeneID, ID = gene_id, SYMBOL = gene_name, gene_type = gene_biotype)

  # Finalise counts_mat: rename columns, drop working columns, arrange
  counts_mat <- counts_mat %>%
    dplyr::select(
      ID        = gene_id,
      SYMBOL    = gene_name,
      gene_type = gene_biotype,
      dplyr::all_of(sample_cols)
    ) %>%
    dplyr::arrange(ID)

  # Step 5: filter TPM to winning entrez_ids, attach annotation from raw counts
  tpm_mat <- tpm_raw %>%
    dplyr::filter(GeneID %in% keep_entrez) %>%
    dplyr::left_join(annotation, by = "GeneID") %>%
    dplyr::select(ID, SYMBOL, gene_type, dplyr::all_of(sample_cols)) %>%
    dplyr::arrange(ID)

  # Sanity check — gene sets must match
  if (!identical(counts_mat$ID, tpm_mat$ID)) {
    warning("  WARNING: Raw and TPM gene sets differ after mapping — check inputs.")
  }

  # ---------------------------------------------------------------------------
  # 3e. Remove all-zero genes
  # ---------------------------------------------------------------------------
  nonzero    <- rowSums(dplyr::select(counts_mat, dplyr::all_of(sample_cols))) > 0
  counts_mat <- counts_mat[nonzero, ]
  tpm_mat    <- tpm_mat[nonzero, ]

  cat("  Genes retained:", nrow(counts_mat),
      "| Samples:", length(sample_cols), "\n")

  # ---------------------------------------------------------------------------
  # 3f. Save outputs
  # ---------------------------------------------------------------------------
  raw_tsv <- file.path(COUNTS_DIR, paste0(proj, "_Raw_Counts.tsv"))
  tpm_tsv <- file.path(COUNTS_DIR, paste0(proj, "_TPM_Counts.tsv"))
  raw_rds <- file.path(COUNTS_DIR, paste0(proj, "_Raw_Counts.rds"))
  tpm_rds <- file.path(COUNTS_DIR, paste0(proj, "_TPM_Counts.rds"))

  readr::write_tsv(counts_mat, raw_tsv)
  readr::write_tsv(tpm_mat,    tpm_tsv)
  saveRDS(counts_mat,          raw_rds)
  saveRDS(tpm_mat,             tpm_rds)

  cat("  Saved:", basename(raw_tsv), "|", basename(tpm_tsv),
      "|", basename(raw_rds), "|", basename(tpm_rds), "\n")
}

cat("\nDone. Outputs in:", COUNTS_DIR, "\n")
