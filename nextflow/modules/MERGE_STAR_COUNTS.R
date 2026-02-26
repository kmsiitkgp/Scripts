#!/usr/bin/env Rscript

# ---- 📦 1. Load Libraries (Always Runs) ----

suppressPackageStartupMessages({
  library(dplyr)      # Provides the %>% operator
  library(readr)
  library(purrr)
  library(openxlsx)
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

merge_counts <- function(counts_dir, output_dir) {
  
  # ---- 🔎 Identify Count Files ----
  
  # Pattern matches STAR (ReadsPerGene.out.tab) files
  count_files <- list.files(path = counts_dir, 
                            pattern = "ReadsPerGene\\.out\\.tab$", 
                            full.names = TRUE)
  
  if (length(count_files) == 0) {
    log_warn(sample = "", step = "merge_counts", msg = "No count files found.")
    return(NULL)
  }
  
  # ---- 🔄 Parse Files and Detect Strandedness ----
  
  # Define metadata rows to exclude (STAR stats)
  special_counters <- c("__no_feature", "__ambiguous", "__too_low_aQual", 
                        "__not_aligned", "__alignment_not_unique", "__assignment_counts",
                        "N_unmapped", "N_multimapping", "N_noFeature", "N_ambiguous")
  
  all_counts <- list()
  
  for (count_file in count_files) {
    
    # Extract Sample ID from filename
    sample_id <- gsub("\\..*$|ReadsPerGene\\.out\\.tab", "", basename(count_file))
    
    # Read a count file and name the columns
    # read.table(file = count_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
    df <- tryCatch({
      readr::read_tsv(file = count_file, 
                      col_names = FALSE, 
                      show_col_types = FALSE) %>%
        dplyr::rename(gene_id    = 1, 
                      unstranded = 2, 
                      forward    = 3, 
                      reverse    = 4)
    }, error = function(e) {
      log_error(sample = sample_id, step = "merge_counts", msg = glue::glue("Error: {e$message}"))
      return(NULL)
    })
    
    if (is.null(df) || ncol(df) < 4) next
    
    # Remove metadata/stats rows
    df <- df %>% dplyr::filter(!(gene_id %in% special_counters))
    
    
    # ---- 🧬 Strandedness Logic ----
    
    fwd_sum <- sum(df$forward, na.rm = TRUE)
    rev_sum <- sum(df$reverse, na.rm = TRUE)
    total   <- fwd_sum + rev_sum
    
    if (total == 0) {
      log_error(sample = sample_id, step = "merge_counts", msg = "CRITICAL: Total counts are zero. Sample excluded.")
      next 
    }
    
    prop_fwd <- fwd_sum / total
    
    # Determine the column to keep
    if (prop_fwd > 0.8) {
      log_info(sample = sample_id, step = "merge_counts", msg = "Detected: Forward Stranded")
      sample_df <- df %>% dplyr::select(ID = gene_id, !!sample_id := forward)
      
    } else if (prop_fwd < 0.2) {
      log_info(sample = sample_id, step = "merge_counts", msg = "Detected: Reverse Stranded")
      sample_df <- df %>% dplyr::select(ID = gene_id, !!sample_id := reverse)
      
    } else {
      # Handle Unstranded (~0.5) and Ambiguous (e.g. 0.7)
      msg <- if(abs(prop_fwd - 0.5) < 0.1) "Detected: Unstranded" else glue::glue("Ambiguous ({round(prop_fwd, 2)}). Using Unstranded.")
      log_info(sample = sample_id, step = "merge_counts", msg = msg)
      sample_df <- df %>% dplyr::select(ID = gene_id, !!sample_id := unstranded)
    }
    
    all_counts[[sample_id]] <- sample_df
  }
  
  if (length(all_counts) == 0) {
    log_error(sample = "", step = "merge_counts", msg = "No valid samples remained after filtering.")
    return(NULL)
  }
  
  # ---- 📊 Build Count Matrix ----
  
  # purrr::reduce handles merging multiple dataframes efficiently
  count_matrix <- all_counts %>% 
    purrr::reduce(dplyr::full_join, by = "ID") 
  
  # Replace any NAs from full_join with 0
  count_matrix[is.na(count_matrix)] <- 0
  
  # Remove genes with 0 counts across all samples
  count_matrix <- count_matrix %>%
    dplyr::filter(rowSums(dplyr::select(., -ID), na.rm = TRUE) > 0)
  
  # Remove samples with 0 counts across all genes
  count_matrix <- count_matrix[, c(TRUE, colSums(count_matrix[,-1]) > 0), drop = FALSE]
  
  # ---- 💾 Save Results ----
  
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  openxlsx::write.xlsx(x = count_matrix, 
                       sheetName = "raw_counts", 
                       file = file.path(output_dir, "STAR_Gene_counts.xlsx"), 
                       overwrite = TRUE)
  
  # ---- 🪵 Log Output and Return ----
  
  log_info(sample = "", 
           step   = "merge_counts", 
           msg    = "Merged counts successfully. Saved to: 'STAR_Gene_counts.xlsx'")
  
  return(invisible(count_matrix))
}

# ---- 🚀 4. Smart Execution (ONLY runs in Nextflow) ----

if (!interactive()) {
  
  get_arg <- function(idx) {
    if (idx > length(args)) return(NULL) # Safety if fewer args provided
    val <- args[idx]
    if (is.na(val) || val == "" || val == "null") return(NULL)
    return(val)
  }
  
  merge_counts(
    counts_dir  = get_arg(1),
    output_dir  = get_arg(2)
  )
}