#!/usr/bin/env Rscript

# ---- 📦 1. Load Libraries (Always Runs) ----

suppressPackageStartupMessages({
  library(dplyr)      # Provides the %>% operator
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
    # Assume project root is current working directory
    return(file.path(getwd(), "modules", "UTILS.R"))
  }
  
  # 3. Non-interactive (Nextflow / Rscript)
  initial.options <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("--file=", initial.options, value = TRUE)
  if (length(file_arg) == 0) stop("Cannot detect script path for UTILS.R!")
  script_dir <- dirname(sub("--file=", "", file_arg))
  return(file.path(script_dir, "UTILS.R"))
}

utils_path <- get_utils_path()
if (!file.exists(utils_path)) stop(paste("❌ UTILS.R not found at:", utils_path))
source(utils_path)

# ---- 🧬 3. Function Definition (Always Runs) ----

create_dds <- function(txi, metadata, tx2gene, output_dir, design = 1, expr_mat = NULL) {
  
  # Load the actual data objects using the paths
  metadata <- load_smart(metadata)
  txi      <- load_smart(txi)
  expr_mat <- load_smart(expr_mat)
  tx2gene  <- load_smart(tx2gene) 
  
  # Remove the first column (the IDs) so only numeric counts remain
  if (!is.null(expr_mat)) {
    rownames(expr_mat) <- make.names(expr_mat[[1]])
    expr_mat           <- expr_mat[, -1, drop = FALSE]
  }

  # ---- 🧪 SECTION 1: DESIGN FORMULA ----
  
  # Standardize DESeq2 design formula
  if (inherits(design, "formula")) {
    design_formula <- design
    
  } else if ( (is.numeric(design) || design == "1") && length(design) == 1 ) {
    # Check if it is the number 1 OR the string "1"
    # Intercept-only model (used for QC or SVA)
    design_formula <- stats::as.formula("~ 1")
    
  } else if (is.character(design) && length(design) == 1 && nzchar(trimws(design))) {
    # Clean leading ~ and whitespace (e.g., " ~ Condition " -> "Condition")
    clean_design <- sub("^\\s*~\\s*", "", design)
    design_formula <- stats::as.formula(paste0("~", clean_design))
    
  } else {
    log_error(sample = "", 
              step = "create_dds", 
              msg = "`design` must be a formula (e.g., ~Batch+Group), the number 1, or a character string.")
  }
  
  # ---- 📝 SECTION 2: METADATA ----
  
  # Remove "sizeFactor" column if present and filter out missing Sample_IDs
  metadata <- metadata %>%
    dplyr::select(-any_of("sizeFactor")) %>%
    dplyr::filter(!is.na(Sample_ID)) %>%
    as.data.frame()  # cannot set rownames on tibble. R wont error, just warn.
  
  # Assign valid "Sample_ID" column as row names for DESeq2
  rownames(metadata) <- make.names(metadata$Sample_ID)
  
  # Add default Batch if missing
  if (!"Batch" %in% colnames(metadata)) {
    metadata$Batch <- 1
    log_warn(sample = "", step = "create_dds", msg = "No 'Batch' column. Assigned default '1'.")
  }
  
  # Remove samples with NA in Design Columns
  design_vars <- all.vars(design_formula)
  na_idx <- apply(X = metadata[, design_vars, drop = FALSE], 
                  MARGIN = 1, 
                  FUN = function(x) any(is.na(x)))
  
  if (any(na_idx)) {
    na_samples <- rownames(metadata)[na_idx]
    log_warn(sample = "", 
             step = "create_dds", 
             msg = glue::glue("Removing {length(na_samples)} samples with NA in design: {paste(na_samples, collapse = ', ')}"))
    metadata <- metadata[!na_idx, , drop = FALSE]
  }
  
  # Convert non-numeric columns to factors
  #metadata <- as.data.frame(unclass(metadata), stringsAsFactors = TRUE)
  #metadata[] <- lapply(metadata, as.factor)
  metadata <- metadata %>%
    dplyr::mutate(Batch = as.factor(Batch),  # Force Sample_ID and Batch to factor even if it's numeric (1, 2, 3)
                  Sample_ID = as.factor(Sample_ID)) %>%
    dplyr::mutate(across(where(~ is.character(.) | is.logical(.)), as.factor))
  

  # ---- 🧮 SECTION 3: READ DATA (TXI OR MATRIX) ----
  
  if (!is.null(txi)) {
    # CASE A: Salmon / tximport
    log_info(sample = "", 
             step = "create_dds", 
             msg = "Input detected: Salmon (tximport).")
    
    read_data <- txi$counts
    
  } else if (!is.null(expr_mat)) {
    # CASE B: STAR / FeatureCounts / Matrix
    log_info(sample = "", 
             step = "create_dds", 
             msg = "Input detected: Raw count matrix.")
    
    # Assign valid "Sample_ID" column as column names for DESeq2
    colnames(expr_mat) <- make.names(colnames(expr_mat))
    
    # Replace missing counts with 0
    expr_mat[is.na(expr_mat)] <- 0
    
    # Convert all columns to numeric safely
    read_data <- matrix(as.numeric(as.matrix(expr_mat)),
                        nrow = nrow(expr_mat),
                        ncol = ncol(expr_mat),
                        dimnames = dimnames(expr_mat))
    
    # Detect if any NA still present in read data
    if (any(is.na(read_data))) {
      log_error(sample = "",
                step   = "create_dds",
                msg    = "`read_data` contains non-numeric values that could not be converted.")
    }
    
    # Remove genes with zero counts across all samples
    read_data <- read_data[rowSums(read_data) > 0, , drop = FALSE]
    
    # Remove samples with zero total reads
    read_data <- read_data[, colSums(read_data) > 0, drop = FALSE]
    
    # Remove rows where gene names are NA or empty strings
    invalid_genes <- is.na(rownames(read_data)) | rownames(read_data) == ""
    read_data <- read_data[!invalid_genes, , drop = FALSE]
    
  } else{
    
    log_error(sample = "", 
              step = "create_dds", 
              msg = "Neither 'txi' nor 'expr_mat' provided.")
  }
  
  
  # ---- 🤝 SECTION 4: ALIGNING METADATA, READ DATA & DESIGN ----
  
  log_info(sample = "", 
           step = "create_dds", 
           msg = "Aligning samples between counts and metadata...")
  
  common_samples <- intersect(colnames(read_data), rownames(metadata))
  
  if (length(common_samples) == 0) {
    log_error(sample = "",
              step = "create_dds", 
              msg = "Zero overlapping samples found!")
  }
  
  # Subset Metadata
  metadata <- metadata[common_samples, , drop = FALSE]
  
  # Subset Counts / TXI
  if (!is.null(txi)) {
    txi <- lapply(txi, function(x) {
      if (is.matrix(x) || is.data.frame(x)) return(x[, common_samples, drop = FALSE])
      return(x)
    })
  } else {
    read_data <- read_data[, common_samples, drop = FALSE]
  }
  
  # 🔍 CROSS-CHECK if design variables exist in metadata
  # all.vars() turns "~ Batch + Condition" into c("Batch", "Condition")
  design_vars <- all.vars(design_formula)
  missing_vars <- setdiff(design_vars, colnames(metadata))
  
  if (length(missing_vars) > 0) {
    log_error(sample = "", 
              step = "create_dds", 
              msg = glue::glue("The following variables in your design are MISSING from your metadata: {paste(missing_vars, collapse = ', ')}"))
  }
  
  
  # ---- 🧪 SECTION 5: DDS OBJECT CREATION ----
  
  # Create DESeqDataSet based on input type
  if (!is.null(txi)) {
    dds <- DESeq2::DESeqDataSetFromTximport(txi = txi, 
                                            colData = metadata, 
                                            design = design_formula)
  } else {
    # Rounding protects against non-integer STAR counts if scaled
    dds <- DESeq2::DESeqDataSetFromMatrix(countData = round(read_data), 
                                          colData = metadata, 
                                          design = design_formula)
  }

  # ---- 📉 SECTION 6: MODEL FITTING (DISPERSION & FIT SELECTION) ----

  log_info(sample = "", 
           step = "create_dds", 
           msg = "Fitting DESeq2 model...")
  
  # Auto-detect for 'poscounts' if many zeros exist
  is_scRNA <- FALSE
  
  if (!is.null(read_data)){
    is_scRNA <- all(rowSums(read_data == 0) > 0)
  } else if (!is.null(txi) && !is.null(txi$counts)){
    is_scRNA <- all(rowSums(txi$counts == 0) > 0)
  }
  
  if (is_scRNA) {
    log_warn(sample = "",
             step   = "create_dds",
             msg    = "scRNA-seq–like sparsity detected. Estimating size factors using 'poscounts'.")
    dds <- DESeq2::estimateSizeFactors(dds, type = "poscounts")
  } else {
    # Pre-filter lowly expressed genes to improve sizefactor estimation in next step
    keep <- rowSums(DESeq2::counts(dds)) >= 10
    dds <- dds[keep, ]
  }
  
  # Fit both Parametric and Local to choose the best one
  dds_para <- DESeq2::DESeq(object = dds, 
                            test = "Wald", 
                            fitType = "parametric", 
                            betaPrior = FALSE, 
                            minReplicatesForReplace = 7, 
                            quiet = TRUE)
  
  dds_local <- DESeq2::DESeq(object = dds, 
                             test = "Wald", 
                             fitType = "local", 
                             betaPrior = FALSE, 
                             minReplicatesForReplace = 7, 
                             quiet = TRUE)
  
  # Calculate residuals to find the better fit
  residual_para <- median((mcols(dds_para)$dispGeneEst - mcols(dds_para)$dispFit)^2, na.rm = TRUE)
  residual_local <- median((mcols(dds_local)$dispGeneEst - mcols(dds_local)$dispFit)^2, na.rm = TRUE)
  
  if (residual_para <= residual_local) {
    dds <- dds_para
    fit_type <- "Parametric"
  } else {
    dds <- dds_local
    fit_type <- "Local"
  }
  
  log_info(sample = "", 
           step = "create_dds", 
           msg = glue::glue("DESeq2 model built using '{fit_type}' fit."))
  
  # ---- 📊 COUNT MATRIX EXTRACTION ----
  
  # --- Median-of-Ratios Normalized Counts ---
  # DESeq2 size-factor normalization (median of ratios method).
  # Corrects for library size differences across samples.
  # USE ONLY FOR: violin plots, box plots, etc visualization of raw expression
  norm_counts_mat <- DESeq2::counts(dds, normalized = TRUE)
  norm_counts     <- process_counts(norm_counts_mat, tx2gene)
  
  # --- Log2 Transformed Normalized Counts ---
  # log2(median-of-ratios normalized counts + 1) via DESeq2::normTransform().
  # The +1 pseudo-count prevents log(0) and stabilizes low-count genes.
  # USE ONLY FOR: violin plots, box plots, etc visualization of log expression;
  #               TF activity inference (VIPER, DoRothEA, decoupleR)
  log_norm_counts_mat <- DESeq2::normTransform(dds) %>%
    SummarizedExperiment::assay()
  log_norm_counts     <- process_counts(log_norm_counts_mat, tx2gene)
  
  # --- VST — Blind ---
  # VST computed without awareness of the experimental design (blind = TRUE).
  # Stabilizes variance across the mean-expression range for visualization.
  # Blind = TRUE is appropriate for QC steps where you do not want design
  # information to influence the transformation (e.g., checking for outliers).
  # USE ONLY FOR: QC plots, PCA, sample correlation / heatmaps.
  vst_blind_mat    <- DESeq2::vst(dds, blind = TRUE) %>%
    SummarizedExperiment::assay()
  vst_blind_counts <- process_counts(vst_blind_mat, tx2gene)
  
  # ---- 💾 Save Results ----
  
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  saveRDS(dds, file = file.path(output_dir, "DESeq2_dds.rds"))
  openxlsx::write.xlsx(x         = list("Norm_Counts" = norm_counts),
                       file      = file.path(output_dir, "Norm_Counts_ExpressionPlots.xlsx"),
                       overwrite = TRUE)
  
  openxlsx::write.xlsx(x         = list("Log_Norm_Counts" = log_norm_counts),
                       file      = file.path(output_dir, "Log_Norm_Counts_ExpressionPlots_TFAnalysis.xlsx"),
                       overwrite = TRUE)
  
  openxlsx::write.xlsx(x         = list("VST_Blind" = vst_blind_counts),
                       file      = file.path(output_dir, "VST_Blind_Counts_QC_PCA_Correlation.xlsx"),
                       overwrite = TRUE)
  
  # ---- 🪵 Log Output and Return ----
  
  log_info(sample = "", 
           step   = "create_dds", 
           msg    = "dds created successfully. Saved to: 'DESeq2_dds.rds'")
  
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