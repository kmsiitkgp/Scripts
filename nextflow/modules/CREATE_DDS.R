#!/usr/bin/env Rscript

# ---- 📦 1. Load Libraries (Always Runs) ----

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tibble)
  library(tximport)
  library(openxlsx)
  library(DESeq2)
})

# ---- 🛠️ 2. Smart Setup (Find & source UTILS.R) ----

get_utils_path <- function() {
  # 1. Windows dev machine
  if (.Platform$OS.type == "windows") {
    return("C:/Users/kailasamms/OneDrive - Cedars-Sinai Health System/Documents/GitHub/Scripts/nextflow/modules/UTILS.R")
  }
  # 2. Interactive Linux / macOS (HPC interactive session)
  if (interactive()) {
    return(file.path(getwd(), "modules", "UTILS.R"))
  }
  # 3. Non-interactive (Nextflow / Rscript)
  initial.options <- commandArgs(trailingOnly = FALSE)
  file_arg        <- grep("--file=", initial.options, value = TRUE)
  if (length(file_arg) == 0) stop("Cannot detect script path for UTILS.R!")
  script_dir      <- dirname(sub("--file=", "", file_arg))
  return(file.path(script_dir, "UTILS.R"))
}

utils_path <- get_utils_path()
if (!file.exists(utils_path)) stop(paste("❌ UTILS.R not found at:", utils_path))
source(utils_path)

# ---- 🧬 3. Function Definition (Always Runs) ----

create_dds <- function(txi, metadata, tx2gene, output_dir, design = 1, expr_mat = NULL) {

  metadata <- load_smart(metadata)
  txi      <- load_smart(txi)
  expr_mat <- load_smart(expr_mat)
  tx2gene  <- load_smart(tx2gene)

  # expr_mat arrives as a dataframe with gene IDs in column 1 and sample counts
  # in the remaining columns. Convert column 1 to rownames so the matrix is
  # purely numeric — required by DESeqDataSetFromMatrix.
  if (!is.null(expr_mat)) {
    rownames(expr_mat) <- make.names(expr_mat[[1]])
    expr_mat           <- expr_mat[, -1, drop = FALSE]
  }

  # ═══════════════════════════════════════════════════════════════════════════
  # SECTION 1: Parse Design Formula
  # ═══════════════════════════════════════════════════════════════════════════
  # Accepts three input forms for maximum flexibility:
  #   Formula object : ~Batch + Group        — passed through unchanged
  #   Character string: "~Batch + Group"     — tilde stripped and rebuilt
  #   Numeric 1 / "1" : intercept-only model — used for QC or SVA runs where
  #                      no biological grouping is needed yet
  #
  # Why not just require a formula object always?
  # Nextflow passes all arguments as strings — the character path handles that
  # case cleanly without requiring the caller to construct an R formula object.

  if (inherits(design, "formula")) {
    design_formula <- design

  } else if ((is.numeric(design) || design == "1") && length(design) == 1) {
    # Intercept-only — valid for exploratory PCA or SVA surrogate variable estimation
    design_formula <- stats::as.formula("~ 1")

  } else if (is.character(design) && length(design) == 1 && nzchar(trimws(design))) {
    # Strip any leading "~" and whitespace before rebuilding — handles both
    # "~Batch+Group" and "Batch+Group" inputs identically
    clean_design   <- sub("^\\s*~\\s*", "", design)
    design_formula <- stats::as.formula(paste0("~", clean_design))

  } else {
    log_error(sample = "", step = "create_dds",
              msg = "`design` must be a formula (e.g. ~Batch+Group), the number 1, or a character string.")
  }

  # ═══════════════════════════════════════════════════════════════════════════
  # SECTION 2: Prepare Metadata
  # ═══════════════════════════════════════════════════════════════════════════

  # ── 2a. Clean and set rownames ───────────────────────────────────────────
  # DESeq2 matches samples between colData and counts by rownames — rownames
  # MUST be set and must match colnames(counts) exactly.
  # Why remove sizeFactor? If a previous dds was saved and metadata re-exported,
  # sizeFactor may be present as a column. Leaving it causes DESeq2 to skip
  # size factor estimation and use the stale values instead.
  # Why as.data.frame()? rownames cannot be set on a tibble — converting first.
  metadata <- metadata %>%
    dplyr::select(-dplyr::any_of("sizeFactor")) %>%
    dplyr::filter(!is.na(Sample_ID)) %>%
    as.data.frame()

  rownames(metadata) <- make.names(metadata$Sample_ID)

  # ── 2b. Ensure Batch column exists ───────────────────────────────────────
  # Batch is expected by downstream QC and batch correction functions. If the
  # metadata doesn't have it (e.g. a single-batch experiment), assign a default
  # of 1 so all samples belong to one batch and downstream code doesn't break.
  if (!"Batch" %in% colnames(metadata)) {
    metadata$Batch <- 1
    log_warn(sample = "", step = "create_dds",
             msg = "No 'Batch' column found. Assigned default value of 1.")
  }

  # ── 2c. Remove samples with NA in design variables ───────────────────────
  # DESeq2 cannot build a model matrix if any sample has NA for a design
  # variable — it will error or silently produce wrong results. Remove these
  # samples early with a clear warning so the user knows what was dropped.
  design_vars <- all.vars(design_formula)
  na_idx      <- apply(X      = metadata[, design_vars, drop = FALSE],
                       MARGIN = 1,
                       FUN    = function(x) any(is.na(x)))

  if (any(na_idx)) {
    na_samples <- rownames(metadata)[na_idx]
    log_warn(sample = "", step = "create_dds",
             msg = glue::glue("Removing {length(na_samples)} sample(s) with NA in design ",
                              "variables: {paste(na_samples, collapse = ', ')}"))
    metadata <- metadata[!na_idx, , drop = FALSE]
  }

  # ── 2d. Coerce columns to factor ─────────────────────────────────────────
  # DESeq2 requires factor columns in colData for categorical design variables.
  # Numeric batch/group IDs (1, 2, 3) would be treated as continuous covariates
  # without this step, producing a wrong model matrix.
  # Why mutate across where(is.character | is.logical)?
  # Logical columns (TRUE/FALSE) and character columns both represent categorical
  # variables and must be factors for DESeq2's model matrix to be correct.
  metadata <- metadata %>%
    dplyr::mutate(Batch     = as.factor(Batch),
                  Sample_ID = as.factor(Sample_ID)) %>%
    dplyr::mutate(dplyr::across(where(~ is.character(.) | is.logical(.)), as.factor))

  # ═══════════════════════════════════════════════════════════════════════════
  # SECTION 3: Load and Validate Count Data
  # ═══════════════════════════════════════════════════════════════════════════
  # Two input sources are supported:
  #   txi     : Salmon quantification imported via tximport — preferred, because
  #             tximport handles transcript-to-gene aggregation and preserves
  #             uncertainty estimates (inferential replicates) that DESeq2 can use
  #   expr_mat: Raw count matrix (e.g. from STAR + featureCounts) — simpler but
  #             loses the transcript-level uncertainty information

  if (!is.null(txi)) {
    log_info(sample = "", step = "create_dds", msg = "Input detected: Salmon (tximport).")
    read_data <- txi$counts

  } else if (!is.null(expr_mat)) {
    log_info(sample = "", step = "create_dds", msg = "Input detected: Raw count matrix.")

    # Sanitise column names — spaces and special characters break DESeq2's
    # colData matching downstream
    colnames(expr_mat)        <- make.names(colnames(expr_mat))
    expr_mat[is.na(expr_mat)] <- 0

    # Coerce to numeric matrix — featureCounts output sometimes has character
    # columns if non-numeric values crept in during file processing
    read_data <- matrix(as.numeric(as.matrix(expr_mat)),
                        nrow     = nrow(expr_mat),
                        ncol     = ncol(expr_mat),
                        dimnames = dimnames(expr_mat))

    if (any(is.na(read_data))) {
      log_error(sample = "", step = "create_dds",
                msg = "`read_data` contains non-numeric values that could not be coerced.")
    }

    # Remove all-zero genes and all-zero samples — DESeq2 cannot fit dispersions
    # for genes with no counts, and zero-count samples indicate failed sequencing
    read_data <- read_data[rowSums(read_data) > 0, , drop = FALSE]
    read_data <- read_data[, colSums(read_data) > 0, drop = FALSE]
    read_data <- read_data[!(is.na(rownames(read_data)) | rownames(read_data) == ""), , drop = FALSE]

  } else {
    log_error(sample = "", step = "create_dds",
              msg = "Neither 'txi' nor 'expr_mat' provided. One is required.")
  }

  # ═══════════════════════════════════════════════════════════════════════════
  # SECTION 4: Align Samples Between Counts and Metadata
  # ═══════════════════════════════════════════════════════════════════════════
  # DESeq2 requires colnames(counts) and rownames(colData) to be identical and
  # in the same order. We subset BOTH to the intersection so extra samples in
  # either don't cause errors — but log any discrepancies clearly.

  log_info(sample = "", step = "create_dds",
           msg = "Aligning samples between counts and metadata...")

  common_samples <- intersect(colnames(read_data), rownames(metadata))

  if (length(common_samples) == 0) {
    log_error(sample = "", step = "create_dds",
              msg = "Zero overlapping samples between counts and metadata. Check Sample_ID formatting.")
  }

  metadata <- metadata[common_samples, , drop = FALSE]

  # For txi, ALL matrices within the list (counts, abundance, length) must be
  # subset identically — subsetting only txi$counts would misalign the others
  if (!is.null(txi)) {
    txi <- lapply(txi, function(x) {
      if (is.matrix(x) || is.data.frame(x)) return(x[, common_samples, drop = FALSE])
      return(x)
    })
  } else {
    read_data <- read_data[, common_samples, drop = FALSE]
  }

  # Verify all design variables are present as columns in the final metadata.
  # Checked here (after subsetting) rather than earlier because subsetting
  # could in theory drop a column — though in practice metadata columns are
  # never dropped by sample subsetting.
  design_vars  <- all.vars(design_formula)
  missing_vars <- setdiff(design_vars, colnames(metadata))

  if (length(missing_vars) > 0) {
    log_error(sample = "", step = "create_dds",
              msg = glue::glue("Design variable(s) missing from metadata: ",
                               "{paste(missing_vars, collapse = ', ')}"))
  }

  # ═══════════════════════════════════════════════════════════════════════════
  # SECTION 5: Build DESeqDataSet
  # ═══════════════════════════════════════════════════════════════════════════
  # DESeqDataSetFromTximport is preferred over DESeqDataSetFromMatrix when txi
  # is available because it internally handles the offset matrix (avgTxLength)
  # that corrects for transcript length bias — important for accurate LFC estimates
  # especially for genes with many/few isoforms.

  if (!is.null(txi)) {
    dds <- DESeq2::DESeqDataSetFromTximport(txi     = txi,
                                            colData = metadata,
                                            design  = design_formula)
  } else {
    # round() is required because featureCounts can produce non-integer counts
    # when scaling is applied — DESeq2 requires integer counts
    dds <- DESeq2::DESeqDataSetFromMatrix(countData = round(read_data),
                                          colData   = metadata,
                                          design    = design_formula)
  }

  # ═══════════════════════════════════════════════════════════════════════════
  # SECTION 6: Model Fitting
  # ═══════════════════════════════════════════════════════════════════════════

  log_info(sample = "", step = "create_dds", msg = "Fitting DESeq2 model...")

  # ── 6a. Sparse data check ────────────────────────────────────────────────
  # DESeq2's default size factor estimation (median-of-ratios) requires at
  # least one non-zero count per gene across all samples to compute the
  # geometric mean. If every gene has at least one zero-count sample (i.e.
  # the data is sparse), this breaks. poscounts uses only positive counts for
  # the geometric mean and is robust to sparse data.
  is_sparse <- FALSE

  if (!is.null(read_data)) {
    is_sparse <- all(rowSums(read_data == 0) > 0)
  } else if (!is.null(txi) && !is.null(txi$counts)) {
    is_sparse <- all(rowSums(txi$counts == 0) > 0)
  }

  if (is_sparse) {
    log_warn(sample = "", step = "create_dds",
             msg = "Sparse count data detected. Using 'poscounts' size factor estimation.")
    dds <- DESeq2::estimateSizeFactors(dds, type = "poscounts")
  } else {
    # Pre-filter lowly expressed genes to reduce noise and improve dispersion
    # estimation convergence. Threshold of 10 total counts is conservative —
    # genes below this are unlikely to yield reliable differential expression calls.
    keep <- rowSums(DESeq2::counts(dds)) >= 10
    dds  <- dds[keep, ]
  }

  # ── 6b. Dispersion model selection ──────────────────────────────────────
  # DESeq2 offers two dispersion fitting strategies:
  #   parametric: assumes a log-linear relationship between dispersion and mean —
  #               works well for well-behaved RNA-seq datasets
  #   local     : fits a local regression (loess) — more flexible, better for
  #               datasets where dispersion doesn't follow a clean log-linear trend
  #               (e.g. very heterogeneous samples, unusual expression distributions)
  #
  # Why fit both and choose?
  # There's no universal rule for which fits better — it depends on the data.
  # We fit both and select by median squared residual (gene-wise dispersion vs
  # fitted dispersion). Lower residual = better fit. This is fully automatic
  # and adds only modest compute time since DESeq is already the slow step.
  dds_para  <- DESeq2::DESeq(object                  = dds,
                              test                    = "Wald",
                              fitType                 = "parametric",
                              betaPrior               = FALSE,
                              minReplicatesForReplace = 7,
                              quiet                   = TRUE)

  dds_local <- DESeq2::DESeq(object                  = dds,
                              test                    = "Wald",
                              fitType                 = "local",
                              betaPrior               = FALSE,
                              minReplicatesForReplace = 7,
                              quiet                   = TRUE)

  residual_para  <- median((mcols(dds_para)$dispGeneEst  - mcols(dds_para)$dispFit)^2,  na.rm = TRUE)
  residual_local <- median((mcols(dds_local)$dispGeneEst - mcols(dds_local)$dispFit)^2, na.rm = TRUE)

  if (residual_para <= residual_local) {
    dds      <- dds_para
    fit_type <- "parametric"
  } else {
    dds      <- dds_local
    fit_type <- "local"
  }

  log_info(sample = "", step = "create_dds",
           msg = glue::glue("Dispersion fit selected: '{fit_type}' ",
                            "(para residual = {round(residual_para, 6)}, ",
                            "local residual = {round(residual_local, 6)})"))

  # ═══════════════════════════════════════════════════════════════════════════
  # SECTION 7: Extract Count Matrices
  # ═══════════════════════════════════════════════════════════════════════════
  # Three count representations are saved, each suited for different uses:
  #
  #   Norm counts     : median-of-ratios normalised, NOT log-transformed
  #                     → use for violin / box plots, expression visualisation
  #   Log norm counts : log2(norm + 1) pseudo-count transform
  #                     → use for TF activity inference (VIPER, decoupleR) which
  #                        expects log-scale expression values; also expression plots
  #   VST blind       : variance-stabilising transform, blind = TRUE (design-unaware)
  #                     → use for QC, PCA, and sample-level correlation heatmaps
  #                        blind = TRUE is correct here because we want QC to be
  #                        unbiased by the design — we're checking data quality,
  #                        not biological signal

  norm_counts_mat <- DESeq2::counts(dds, normalized = TRUE)
  norm_counts     <- process_counts(norm_counts_mat, tx2gene)

  # +1 pseudo-count prevents log2(0) = -Inf for zero-count genes
  log_norm_counts_mat <- DESeq2::normTransform(dds) %>% SummarizedExperiment::assay()
  log_norm_counts     <- process_counts(log_norm_counts_mat, tx2gene)

  vst_blind_mat    <- DESeq2::vst(dds, blind = TRUE) %>% SummarizedExperiment::assay()
  vst_blind_counts <- process_counts(vst_blind_mat, tx2gene)

  # ═══════════════════════════════════════════════════════════════════════════
  # SECTION 8: Save and Return
  # ═══════════════════════════════════════════════════════════════════════════
  # File naming encodes the intended downstream use so outputs are self-documenting
  # when browsing the output directory.

  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

  saveRDS(dds, file = file.path(output_dir, "DESeq2_dds.rds"))

  openxlsx::write.xlsx(x         = list("Norm_Counts"     = norm_counts),
                       file      = file.path(output_dir, "Norm_Counts_ExpressionPlots.xlsx"),
                       overwrite = TRUE)

  openxlsx::write.xlsx(x         = list("Log_Norm_Counts" = log_norm_counts),
                       file      = file.path(output_dir, "Log_Norm_Counts_ExpressionPlots_TFAnalysis.xlsx"),
                       overwrite = TRUE)

  openxlsx::write.xlsx(x         = list("VST_Blind"       = vst_blind_counts),
                       file      = file.path(output_dir, "VST_Blind_Counts_QC_PCA_Correlation.xlsx"),
                       overwrite = TRUE)

  log_info(sample = "", step = "create_dds",
           msg    = glue::glue("dds built with '{fit_type}' dispersion fit. Saved to: {output_dir}"))

  return(invisible(dds))
}

# ---- 🚀 4. Smart Execution (Nextflow Only) ----

if (!interactive()) {
  args <- commandArgs(trailingOnly = TRUE)

  get_arg <- function(idx, default = NULL) {
    if (idx > length(args)) return(default)
    val <- args[idx]
    if (is.na(val) || val == "" || val == "null" || val == "NULL") return(default)
    return(val)
  }

  create_dds(
    txi        = get_arg(1),
    metadata   = get_arg(2),
    tx2gene    = get_arg(3),
    output_dir = get_arg(4, "."),
    design     = get_arg(5, "1"),
    expr_mat   = get_arg(6)
  )
}
