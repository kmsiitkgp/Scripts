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

extract_deseq2_results <- function(contrast, dds_rds_path, tx2gene_csv_path, output_dir, batch_vars = NULL, lfc_cutoff = 0, padj_cutoff = 0.1) {
  
  # Load the actual data objects using the paths
  dds  <- readRDS(dds_rds_path)
  
  # For ashr
  set.seed(1234)
  
  # ---- 🧮 Dynamic Contrast Parsing (PARSER METHOD) ----
  
  # Extract design variables (e.g., c("Condition", "Batch"))
  # all.vars is robust; it handles ~Batch+Condition and ~Batch+Condition*Treatment
  design_vars <- all.vars(DESeq2::design(dds))
  
  # Extract the model matrix 
  # This is useful for downstream checking of rank deficiency
  mod_mat <- stats::model.matrix(DESeq2::design(dds), 
                                 data = SummarizedExperiment::colData(dds))
  
  # Create "Groups" column for easy subsetting/contrasts
  if (length(design_vars) > 0) {
    df_groups <- SummarizedExperiment::colData(dds) %>% 
      as.data.frame() %>% 
      tidyr::unite(col = "Groups", dplyr::all_of(design_vars), sep = "_", remove = FALSE)
    
    # Update the dds colData to include this combined Groups column
    SummarizedExperiment::colData(dds)$Groups <- as.factor(df_groups$Groups)
  } else {
    # If design is ~1, Groups is just "All"
    df_groups <- SummarizedExperiment::colData(dds)$Groups <- as.factor("All")
  }
  
  # Define all possible groups that could be compared based on the design
  groups <- unique(df_groups$Groups)
  
  # Map groups to mean model coefficients (Get all possible coefficient vectors)
  group_coef_list <- lapply(groups, function(i) colMeans(as.matrix(mod_mat[df_groups$Groups == i, , drop = FALSE])))
  names(group_coef_list) <- groups
  
  # Recursive function to evaluate the contrast string (e.g., "A-B") as a vector operation
  replace_symbols <- function(node) {
    if (is.symbol(node)) {
      nm <- as.character(node)
      if (nm %in% names(group_coef_list)) {
        return(group_coef_list[[nm]])
      } else {
        # it's an operator like "-" or "+" → return unchanged
        return(node)
      }
    } else if (is.call(node)) {
      return(as.call(lapply(node, replace_symbols)))
    } else {
      return(node)
    }
  }
  
  parsed_expr <- base::parse(text = contrast)[[1]]
  expr_sub <- replace_symbols(parsed_expr)
  contrast_parser <- base::eval(expr_sub)
  
  # ⚠️ Sanity Check: contrast_parser Alignment
  if (length(contrast_parser) != ncol(mod_mat)) {
    log_error(sample = contrast,
              step   = "extract_deseq2_results",
              msg    = glue::glue("Length of contrast_parser ({length(contrast_parser)}) does 
                                   not match the number of columns in the design matrix ({ncol(mod_mat)})."))
  }
  
  if (all(contrast_parser == 0)) {
    log_error(sample = contrast,
              step   = "extract_deseq2_results",
              msg    = "contrast_parser is all zeros — invalid contrast. Check the user-provided contrast string.")
  }
  
  missing_groups <- setdiff(all.vars(parsed_expr), names(group_coef_list))
  if (length(missing_groups) > 0) {
    log_error(sample = contrast,
              step   = "extract_deseq2_results",
              msg    = glue::glue("These groups in contrast do not exist in the design: {paste(missing_groups, collapse=', ')}"))
  }
  
  if (any(is.na(contrast_parser))) {
    log_error(sample = contrast,
              step   = "extract_deseq2_results",
              msg    = "Contrast vector contains NA — possible missing or empty group.")
  }
  
  # ---- 🧹 Zero out user-specified batch covariates ----
  
  if (!is.null(batch_vars)) {
    
    # Find model matrix columns corresponding to nuisance vars
    batch_cols <- grep(paste0("^(", paste(batch_vars, collapse="|"), ")"), 
                          names(contrast_parser), value = TRUE)
    
    if (length(batch_cols) == 0) {
      log_error(sample = contrast,
                step   = "extract_deseq2_results",
                msg    = glue::glue("Nuisance vars {paste(batch_vars, collapse=', ')} 
                                   not found in model matrix columns. 
                                   Are they included in the design formula?"))
    }
    
    # Confounding check — warn if batch is confounded with groups
    for (nv in batch_vars) {
      if (nv %in% colnames(colData(dds))) {
        tab <- table(df_groups$Groups, SummarizedExperiment::colData(dds)[[nv]])
        if (any(rowSums(tab == 0) == ncol(tab) - 1)) {
          log_error(sample = contrast,
                    step   = "extract_deseq2_results",
                    msg    = glue::glue("Nuisance variable '{nv}' appears confounded 
                                       with Groups — results unreliable!"))
        }
      }
    }
    
    log_info(sample = contrast,
             step = "extract_deseq2_results", 
             msg = glue::glue("ℹ️  Zeroing nuisance covariate(s) from contrast vector: '{paste(batch_cols, collapse = ", ")}'"))
    contrast_parser[batch_cols] <- 0
  }
  
  # ---- 🧮 Dynamic Contrast (STANDARD METHOD) ----

  design_vars <- all.vars(design(dds))
  condition_col <- tail(design_vars, 1)
  
  all_contrasts <- DESeq2::resultsNames(dds)
  terms <- all.vars(as.formula(paste0("~", contrast)))
  
  if (length(terms) == 2){
    
    # --- SIMPLE PAIRWISE (e.g., A - B) ---
    contrast_standard <- c(condition_col, terms[1], terms[2])
    
  } else {
    # --- COMPLEX INTERACTION ((A-B) - (C-D)) ---
    contrast_standard <- contrast_parser
  }
  
  # ---- 🧬 Results & LFC Shrinkage ----
  
  res_parser <- DESeq2::results(object               = dds, 
                                contrast             = contrast_parser,
                                lfcThreshold         = lfc_cutoff,
                                altHypothesis        = "greaterAbs",
                                cooksCutoff          = TRUE,
                                independentFiltering = TRUE,
                                alpha                = padj_cutoff,
                                pAdjustMethod        = "BH")
  
  # Use 'ashr' for shrinkage as it is robust for varied effect sizes
  res_parser <- DESeq2::lfcShrink(dds = dds, res = res_parser, type = "ashr", quiet = TRUE)
  #summary(res_parser)
  
  res_standard <- DESeq2::results(object             = dds, 
                                contrast             = contrast_standard,
                                lfcThreshold         = lfc_cutoff,
                                altHypothesis        = "greaterAbs",
                                cooksCutoff          = TRUE,
                                independentFiltering = TRUE,
                                alpha                = padj_cutoff,
                                pAdjustMethod        = "BH")
  
  # Use 'ashr' for shrinkage as it is robust for varied effect sizes
  res_standard <- DESeq2::lfcShrink(dds = dds, res = res_standard, type = "ashr", quiet = TRUE)
  #summary(res_standard)
  
  # ---- 🔍 Detailed Comparison of Results ----
  
  # 1. Convert to dataframes for comparison
  df_p <- as.data.frame(res_parser)
  df_s <- as.data.frame(res_standard)
  
  # 2. Check rownames match
  if (!all(rownames(df_p) == rownames(df_s))) {
    stop("Row names do not match between Parser and Standard results!")
  }
  
  # 3. Check log2FoldChange within tight tolerance
  tolerance <- 1e-8
  lfc_diff <- abs(df_p$log2FoldChange - df_s$log2FoldChange)
  lfc_is_equal <- all(lfc_diff < tolerance, na.rm = TRUE)
  max_lfc_diff <- max(lfc_diff, na.rm = TRUE)
  
  # 4. Check padj within tight tolerance
  padj_diff <- abs(df_p$padj - df_s$padj)
  padj_is_equal <- all(padj_diff < tolerance, na.rm = TRUE)
  max_padj_diff <- max(padj_diff, na.rm = TRUE)
  
  # 5. Compare DEG sets
  deg_p <- rownames(df_p)[which(df_p$padj < padj_cutoff & abs(df_p$log2FoldChange) > lfc_cutoff)]
  deg_s <- rownames(df_s)[which(df_s$padj < padj_cutoff & abs(df_s$log2FoldChange) > lfc_cutoff)]
  
  shared_degs <- intersect(deg_p, deg_s)
  only_parser <- setdiff(deg_p, deg_s)
  only_standard <- setdiff(deg_s, deg_p)
  
  # 6. Print Stringent Comparison Report
  cat("\n--- ⚖️  Parser vs Standard Comparison Report (Stringent) ---\n")
  cat("Numerical LFC Match: ", lfc_is_equal, " (Max Δ = ", formatC(max_lfc_diff, format='e', digits=2), ")\n", sep="")
  cat("Numerical padj Match:", padj_is_equal, " (Max Δ = ", formatC(max_padj_diff, format='e', digits=2), ")\n", sep="")
  cat("\nDEG Set Comparison:\n")
  cat("- Shared DEGs:       ", length(shared_degs), "\n")
  cat("- Only in Parser:    ", length(only_parser), "\n")
  cat("- Only in Standard:  ", length(only_standard), "\n")
  
  if (!lfc_is_equal | !padj_is_equal | length(only_parser) > 0 | length(only_standard) > 0) {
    cat("\n⚠️  Warning: Differences detected! Check p-value or LFC precision.\n")
  } else {
    cat("\n✅ Success: Parser and Standard results are identical (stringent check).\n")
  }
  
  cat("-----------------------------------------------\n")
  
  # For complex interactions, standard gets messy
  # dds_standard <- fit_deseq2_model(expr_mat    = deseq2_data$expr_mat, 
  #                                  txi         = deseq2_data$txi, 
  #                                  metadata    = deseq2_data$metadata,
  #                                  design      = "Cell_Line * Condition * Treatment")
  # 
  # contrast_vec <- setNames(rep(0, length(resultsNames(dds_standard))), 
  #                          resultsNames(dds_standard))
  # 
  # contrast_vec["ConditionWT.Treatment4Gy"] <- -1
  # contrast_vec["Cell_LineARCaPM.ConditionWT.Treatment4Gy"] <- -1
  # 
  # res_standard_ARCaPM <- results(dds_standard, contrast = contrast_vec)
  
  # -------------------TESTING COMPLETE ------------
  
  # ---- 📊 Data Formatting ----

  process_degs <- function(res_obj, tx2gene_path) {

    # Convert to DF and initial cleanup
    df <- as.data.frame(res_obj) %>%
      tibble::rownames_to_column("gene_id")

    # Calculate the minimum non-zero padj safely for this specific dataset
    # We do this outside mutate to avoid the 'padj not found' error
    min_padj <- min(df$padj[df$padj > 0], na.rm = TRUE)

    df %>%
      dplyr::mutate(padj = case_when(is.na(padj) ~ 1,
                                     padj == 0   ~ min_padj,
                                     TRUE        ~ padj)) %>%
      add_annotation(tx2gene_path) %>%
      dplyr::filter(!is.na(SYMBOL)) %>%
      dplyr::group_by(SYMBOL) %>%
      dplyr::slice_min(order_by = padj, n = 1, with_ties = FALSE) %>%
      dplyr::ungroup()
  }

  # Use the function (Much cleaner!)
  DEGs_parser   <- process_degs(res_parser, tx2gene_csv_path)
  DEGs_standard <- process_degs(res_standard, tx2gene_csv_path)

  # ---- 💾 Save Results ----

  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

  # Save both sheets in ONE Excel file for easier comparison
  openxlsx::write.xlsx(x = list("DEGs_Parser" = DEGs_parser, "DEGs_Standard" = DEGs_standard),
                       file = file.path(output_dir, "DEGs_DESeq2.xlsx"),
                       overwrite = TRUE)

  # Keep your RDS files for R-specific backup
  saveRDS(res_parser, file = file.path(output_dir, "res_parser.rds"))
  saveRDS(res_standard, file = file.path(output_dir, "res_standard.rds"))
  
  # ---- 🪵 Log Output and Return ----
  
  log_info(sample = "", 
           step   = "extract_deseq2_results", 
           msg    = "DEG results extracted successfully")
  
  return(invisible(list(df_parser   = DEGs_parser,
                        df_standard = DEGs_standard,
                        res_parser  = res_parser,
                        res_standard = res_standard)))
}

# ---- 🚀 4. Smart Execution (ONLY runs in Nextflow) ----

if (!interactive()) {
  
  # Helper to get args with specific default fallbacks
  get_arg <- function(idx, default) {
    val <- args[idx]
    if (is.na(val) || val == "" || val == "null" || length(args) < idx) {
      return(default)
    }
    return(val)
  }
  
  # Execute with explicit types
  extract_deseq2_results(
    contrast         = get_arg(1, NULL),
    dds_rds_path     = get_arg(2, NULL),
    tx2gene_csv_path = get_arg(3, NULL),
    output_dir       = get_arg(4, "."), 
    # Batch vars should stay NULL if not provided
    batch_vars       = get_arg(5, NULL),
    # Numeric defaults must be coerced
    lfc_cutoff       = as.numeric(get_arg(6, 0)),
    padj_cutoff      = as.numeric(get_arg(7, 0.1))
  )
}