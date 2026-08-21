# =============================================================================
# 03_create_dds.R
# Builds a DESeq2 DESeqDataSet per cancer and saves normalized count matrices
# for downstream analysis (survival, GSEA, TF activity etc.)
#
# Run order: 01 → 02 → 03
# =============================================================================

library(dplyr)
library(tibble)
library(openxlsx)
library(DESeq2)
library(SummarizedExperiment)

# =============================================================================
# INPUTS / OUTPUTS
# =============================================================================

BASE_DIR       <- "/home/kailasamms/scripts/nextflow/resources/"

# Input directories
COUNTS_DIR     <- file.path(BASE_DIR, "Datasets", "02.Counts")    # output of Script 02
METADATA_DIR   <- file.path(BASE_DIR, "Datasets", "01.Metadata")  # output of Script 01

# Output directory
NORMALIZED_DIR <- file.path(BASE_DIR, "Datasets", "03.Normalized")

dir.create(NORMALIZED_DIR, recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# HELPER: RESTORE SYMBOL AND DEDUPLICATE
# =============================================================================
# DESeq2 uses Ensembl IDs as rownames internally. This function adds SYMBOL
# back to any normalized count matrix after DESeq2 processing.
#
# SYMBOL assignment logic:
#   Has SYMBOL                 → use it (any biotype)
#   No SYMBOL + protein_coding → use Ensembl ID as fallback (keep it)
#   No SYMBOL + other biotype  → drop (unnamed non-coding RNA, pseudogenes etc.)
#
# Duplicate SYMBOL resolution:
#   Multiple Ensembl IDs can map to the same SYMBOL (paralogs, annotation
#   artifacts). Keep the row with highest mean expression across all samples —
#   this is the most likely "real" gene; others are typically pseudogenes or
#   misannotated entries with low signal.

add_symbol_and_dedup <- function(mat, gene_info) {

  row_means <- rowMeans(as.matrix(mat), na.rm = TRUE)

  as.data.frame(mat) %>%
    tibble::rownames_to_column("ID") %>%
    dplyr::mutate(row_mean = row_means) %>%
    dplyr::left_join(gene_info %>% dplyr::select(ID, SYMBOL, gene_type),
                     by = "ID") %>%
    dplyr::mutate(
      SYMBOL = dplyr::case_when(
        !is.na(SYMBOL) & SYMBOL != ""  ~ SYMBOL,
        gene_type == "protein_coding"  ~ ID,
        TRUE                           ~ NA_character_
      )
    ) %>%
    dplyr::filter(!is.na(SYMBOL)) %>%
    dplyr::group_by(SYMBOL) %>%
    dplyr::slice_max(order_by = row_mean, n = 1, with_ties = FALSE) %>%
    dplyr::ungroup() %>%
    dplyr::select(ID, SYMBOL, dplyr::everything(), -gene_type, -row_mean)
}

# =============================================================================
# STEP 1: FIND PER-CANCER COUNT FILES
# =============================================================================

cat("Scanning counts directory...\n")

rds_files <- list.files(
  path       = COUNTS_DIR,
  pattern    = "_Raw_Counts\\.rds$",
  full.names = TRUE
)

if (length(rds_files) == 0) {
  stop("No Raw_Counts.rds files found in: ", COUNTS_DIR)
}

cat("  Found", length(rds_files), "cancer(s)\n")

# =============================================================================
# STEP 2: PROCESS EACH CANCER
# =============================================================================

for (rds_file in sort(rds_files)) {

  proj <- sub("_Raw_Counts\\.rds$", "", basename(rds_file))
  cat("\n", proj, "\n", sep = "")

  # ---------------------------------------------------------------------------
  # Skip if already processed
  # ---------------------------------------------------------------------------
  dds_rds_check <- file.path(NORMALIZED_DIR, paste0(proj, "_dds.rds"))
  if (file.exists(dds_rds_check)) {
    cat("  SKIPPING", proj, "— dds.rds already exists\n")
    next
  }

  # ---------------------------------------------------------------------------
  # 2a. Read raw counts matrix
  # ---------------------------------------------------------------------------
  counts_raw <- readRDS(rds_file)

  # Save gene annotation before dropping these columns for DESeq2
  gene_info <- counts_raw %>%
    dplyr::select(ID, SYMBOL, gene_type) %>%
    dplyr::distinct()

  # DESeq2 needs a pure numeric matrix with Ensembl IDs as rownames.
  # Keep dashes in Sample_IDs — do NOT use make.names() which converts
  # TCGA-3C-AAAU-01A → TCGA.3C.AAAU.01A and breaks joins with metadata.
  counts_mat <- counts_raw %>%
    dplyr::select(-SYMBOL, -gene_type) %>%
    tibble::column_to_rownames("ID") %>%
    as.matrix()

  cat("  Genes:", nrow(counts_mat), "| Samples:", ncol(counts_mat), "\n")

  # ---------------------------------------------------------------------------
  # 2b. Read metadata
  # ---------------------------------------------------------------------------
  meta_file <- file.path(METADATA_DIR, paste0(proj, "_metadata.xlsx"))

  if (!file.exists(meta_file)) {
    cat("  SKIPPING", proj, "— metadata file not found\n")
    next
  }

  metadata <- openxlsx::read.xlsx(
    xlsxFile = meta_file,
    #sheet    = proj,
    colNames = TRUE
  )

  if (!"Sample_ID" %in% colnames(metadata)) {
    cat("  SKIPPING", proj, "— Sample_ID column missing\n")
    next
  }

  metadata <- metadata %>%
    dplyr::filter(!is.na(Sample_ID)) %>%
    as.data.frame()
  rownames(metadata) <- metadata$Sample_ID

  cat("  Metadata samples:", nrow(metadata), "\n")

  # ---------------------------------------------------------------------------
  # 2c. Match samples — keep intersection only
  # ---------------------------------------------------------------------------
  common_samples <- intersect(colnames(counts_mat), rownames(metadata))
  only_counts    <- setdiff(colnames(counts_mat),  rownames(metadata))
  only_meta      <- setdiff(rownames(metadata),     colnames(counts_mat))

  cat("  Common samples:", length(common_samples), "\n")

  if (length(only_counts) > 0)
    cat("  In counts but not metadata:", length(only_counts), "\n")
  if (length(only_meta) > 0)
    cat("  In metadata but not counts:", length(only_meta), "\n")

  if (length(common_samples) == 0) {
    cat("  SKIPPING", proj, "— no common samples\n")
    next
  }

  counts_mat <- counts_mat[, common_samples, drop = FALSE]
  metadata   <- metadata[common_samples, , drop = FALSE]

  # ---------------------------------------------------------------------------
  # 2d. Build DESeqDataSet
  # design = ~1 (intercept only) — we want size-factor normalized expression
  # matrices for ALL samples regardless of Tissue_Type. Tumor vs Normal DEG
  # analysis is done separately in downstream scripts using specific contrasts.
  # ~1 gives unbiased size factor estimation across all sample types.
  # ---------------------------------------------------------------------------
  cat("  Building DESeqDataSet...\n")

  metadata <- metadata %>%
    dplyr::mutate(dplyr::across(where(~ is.character(.) | is.logical(.)),
                                as.factor))

  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = round(counts_mat),  # round() handles any non-integer TPM values
    colData   = metadata,
    design    = ~ 1
  )

  # ---------------------------------------------------------------------------
  # 2e. Pre-filter lowly expressed genes
  # rowSums >= 10 removes absolute zeros and near-zeros to speed up DESeq2.
  # This is intentionally lenient — DESeq2's independent filtering during
  # results() handles biological stringency automatically.
  # ---------------------------------------------------------------------------
  keep <- rowSums(DESeq2::counts(dds)) >= 10
  dds  <- dds[keep, ]
  cat("  Genes after low-count filter:", nrow(dds), "\n")

  # ---------------------------------------------------------------------------
  # 2f. Size factor estimation
  # Default median-of-ratios method fails when every gene has at least one
  # zero-count sample (geometric mean = 0). poscounts uses only positive counts
  # and is robust to sparse data — common in smaller TCGA cohorts.
  # ---------------------------------------------------------------------------
  is_sparse <- all(rowSums(DESeq2::counts(dds) == 0) > 0)

  if (is_sparse) {
    cat("  Sparse data — using poscounts size factor estimation\n")
    dds <- DESeq2::estimateSizeFactors(dds, type = "poscounts")
  }

  # ---------------------------------------------------------------------------
  # 2g. Fit DESeq2 model — select best dispersion fit automatically
  # Parametric fit assumes log-linear dispersion-mean relationship (faster,
  # works well for most cancers). Local (loess) fit is more flexible for
  # heterogeneous datasets. Select by median squared residual — lower = better.
  # ---------------------------------------------------------------------------
  cat("  Fitting parametric dispersion...\n")
  dds_para <- DESeq2::DESeq(
    object                  = dds,
    test                    = "Wald",
    fitType                 = "parametric",
    betaPrior               = FALSE,
    minReplicatesForReplace = 7,
    quiet                   = TRUE
  )

  cat("  Fitting local dispersion...\n")
  dds_local <- DESeq2::DESeq(
    object                  = dds,
    test                    = "Wald",
    fitType                 = "local",
    betaPrior               = FALSE,
    minReplicatesForReplace = 7,
    quiet                   = TRUE
  )

  resid_para  <- median((mcols(dds_para)$dispGeneEst  -
                           mcols(dds_para)$dispFit)^2,  na.rm = TRUE)
  resid_local <- median((mcols(dds_local)$dispGeneEst -
                           mcols(dds_local)$dispFit)^2, na.rm = TRUE)

  if (resid_para <= resid_local) {
    dds      <- dds_para
    fit_type <- "parametric"
  } else {
    dds      <- dds_local
    fit_type <- "local"
  }

  cat("  Dispersion fit:", fit_type,
      "(parametric =", round(resid_para, 6),
      "/ local =", round(resid_local, 6), ")\n")

  # ---------------------------------------------------------------------------
  # 2h. Extract normalized count matrices
  #
  # Norm_Counts     : median-of-ratios normalized, not log transformed.
  #                   Use for: violin/box plots, expression visualization.
  #
  # Log_Norm_Counts : log2(norm + 1).
  #                   Use for: TF activity inference (VIPER, decoupleR).
  #
  # VST_Counts      : variance stabilizing transform (blind = TRUE).
  #                   blind = TRUE is correct here because design = ~1.
  #                   Use for: PCA, clustering, correlation heatmaps.
  # ---------------------------------------------------------------------------
  cat("  Extracting normalized matrices...\n")

  norm_mat    <- DESeq2::counts(dds, normalized = TRUE)
  lognorm_mat <- DESeq2::normTransform(dds) %>% SummarizedExperiment::assay()
  vst_mat     <- DESeq2::vst(dds, blind = TRUE) %>% SummarizedExperiment::assay()

  norm_counts    <- add_symbol_and_dedup(norm_mat,    gene_info)
  lognorm_counts <- add_symbol_and_dedup(lognorm_mat, gene_info)
  vst_counts     <- add_symbol_and_dedup(vst_mat,     gene_info)

  cat("  Genes in output matrices:", nrow(norm_counts), "\n")

  # ---------------------------------------------------------------------------
  # 2i. Save outputs
  # ---------------------------------------------------------------------------
  cat("  Saving outputs...\n")

  dds_rds     <- file.path(NORMALIZED_DIR, paste0(proj, "_dds.rds"))
  norm_tsv    <- file.path(NORMALIZED_DIR, paste0(proj, "_Norm_Counts.tsv"))
  lognorm_tsv <- file.path(NORMALIZED_DIR, paste0(proj, "_Log_Norm_Counts.tsv"))
  vst_tsv     <- file.path(NORMALIZED_DIR, paste0(proj, "_VST_Counts.tsv"))

  saveRDS(dds, dds_rds)

  # Save as TSV (tab-delimited)
  readr::write_tsv(norm_counts,    norm_tsv)
  readr::write_tsv(lognorm_counts, lognorm_tsv)
  readr::write_tsv(vst_counts,     vst_tsv)

  cat("  Saved:", basename(dds_rds), "|",
      basename(norm_tsv), "|",
      basename(lognorm_tsv), "|",
      basename(vst_tsv), "\n")
}

cat("\nDone. Outputs in:", NORMALIZED_DIR, "\n")
