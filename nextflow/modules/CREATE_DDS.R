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

# ---- 🛠️ 2. Smart Setup(ONLY runs in Nextflow) ----

if (!interactive()) {
  
  # Source the custom functions from utils.R
  initial.options <- commandArgs(trailingOnly = FALSE)
  script.name <- sub("--file=", "", initial.options[grep("--file=", initial.options)])
  source(file.path(dirname(script.name), "UTILS.R"))
  
  # Capture command line arguments
  args <- commandArgs(trailingOnly = TRUE)
}

# ---- 🧬 3. Function Definition (Always Runs) ----

create_dds <- function(txi_file_path, meta_file_path, output_dir, design = 1, expr_mat_file_path = NULL) {
  
  # Load the actual data objects using the paths
  metadata <- openxlsx::read.xlsx(meta_file_path)
  txi      <- NULL
  expr_mat <- NULL
  
  if (!is.null(txi_file_path)) {
    txi <- readRDS(txi_file_path)
    
    log_info(sample = "", step = "create_dds", msg = "txi loaded.")
  } 
  
  if (!is.null(expr_mat_file_path)) {
    temp_mat <- openxlsx::read.xlsx(expr_mat_file_path)
    
    # Use the first column as rownames
    rownames(temp_mat) <- make.names(temp_mat[[1]])
    
    # Remove the first column (the IDs) so only numeric counts remain
    expr_mat <- temp_mat[, -1, drop = FALSE]
    
    log_info(sample = "", step = "create_dds", msg = "Matrix loaded and Gene IDs assigned to rownames.")
  }
  
  # ---- 🧪 SECTION 1: DESIGN FORMULA ----
  
  # Standardize DESeq2 design formula
  if (inherits(design, "formula")) {
    design_formula <- design
    
  } else if (is.numeric(design) && length(design) == 1 && design == 1) {
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
    dplyr::filter(!is.na(Sample_ID))
  
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
  
  
  # ---- 💾 Save Results ----
  
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  saveRDS(dds, file = file.path(output_dir, "DESeq2_dds.rds"))
  
  # ---- 🪵 Log Output and Return ----
  
  log_info(sample = "", 
           step   = "create_dds", 
           msg    = "dds created successfully. Saved to: 'DESeq2_dds.rds'")
  
  return(invisible(dds))
}

# ---- 🚀 4. Smart Execution (ONLY runs in Nextflow) ----

if (!interactive()) {
  
  get_arg <- function(idx) {
    if (idx > length(args)) return(NULL) # Safety if fewer args provided
    val <- args[idx]
    if (is.na(val) || val == "" || val == "null") return(NULL)
    return(val)
  }
  
  create_dds(
    txi_file_path      = get_arg(1),
    meta_file_path     = get_arg(2),
    output_dir         = get_arg(3),
    design             = get_arg(4), 
    expr_mat_file_path = get_arg(5)
  )
}