#!/usr/bin/env Rscript

# ---- 📦 1. Load Libraries (Always Runs) ----

suppressPackageStartupMessages({
  library(dplyr)      # Provides the %>% operator
  library(DESeq2)
  library(SummarizedExperiment)
  library(ggplot2)
  library(tidyr)
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

deep_audit <- function(base_path) {
  
  cat(">>> Starting Deep Audit in:", base_path, "\n")
  
  # --- 1. DYNAMIC SUMMARY DISCOVERY (The Hunter) ---
  found_salmon_summary <- list.files(base_path, pattern = "Salmon_Master_Summary\\.txt$", recursive = TRUE, full.names = TRUE)
  found_star_summary   <- list.files(base_path, pattern = "STAR_Master_Summary\\.txt$", recursive = TRUE, full.names = TRUE)
  
  salmon_file <- if(length(found_salmon_summary) > 0) found_salmon_summary[1] else NULL
  star_file   <- if(length(found_star_summary) > 0)   found_star_summary[1]   else NULL
  
  salmon_stats <- if(!is.null(salmon_file)) {
    cat("   + Found Salmon Summary:", basename(salmon_file), "\n")
    read.delim(salmon_file, check.names=FALSE, stringsAsFactors=FALSE)
  } else NULL
  
  star_stats <- if(!is.null(star_file)) {
    cat("   + Found STAR Summary:", basename(star_file), "\n")
    read.delim(star_file, check.names=FALSE, stringsAsFactors=FALSE)
  } else NULL
  
  # Prepare metadata for matching
  if(!is.null(salmon_stats)) salmon_stats$CleanDir <- basename(dirname(salmon_stats$Directory_ID))
  if(!is.null(star_stats))   star_stats$CleanDir   <- basename(dirname(star_stats$Directory_ID))
  
  potential_species <- c("Xenograft", "Human", "Mouse")
  species_list <- potential_species[dir.exists(file.path(base_path, potential_species))]
  
  all_rows_list <- list()
  
  for(sp in species_list) {
    cat(paste0("\n>>> Processing Folder: ", sp, "\n"))
    sp_dir <- file.path(base_path, sp)
    
    # --- 2. DYNAMIC FILE DISCOVERY (Searching for .rds) ---
    found_dds <- list.files(sp_dir, pattern = "dds\\.rds$", recursive = TRUE, full.names = TRUE)
    found_txi <- list.files(sp_dir, pattern = "txi\\.rds$", recursive = TRUE, full.names = TRUE)
    
    dds_path <- if(length(found_dds) > 0) found_dds[1] else NULL
    txi_path <- if(length(found_txi) > 0) found_txi[1] else NULL
    
    if(is.null(dds_path)) {
      cat("   ! No dds.rds found in", sp, "- skipping.\n")
      next
    }
    
    # Load Primary Objects
    dds <- readRDS(dds_path)
    all_genes <- rownames(dds)
    cnts      <- DESeq2::counts(dds)
    lib_sizes <- colSums(cnts)
    total_c   <- sum(cnts)
    
    # Load TXI Stats
    n_txi_val <- NA; n_zero_txi_val <- NA
    if(!is.null(txi_path)) {
      txi_obj <- readRDS(txi_path)
      counts_mtx <- if(is.list(txi_obj) && "counts" %in% names(txi_obj)) txi_obj$counts else txi_obj
      if(!is.null(counts_mtx)) {
        n_txi_val      <- nrow(counts_mtx)
        n_zero_txi_val <- sum(rowSums(counts_mtx) == 0)
      }
    }
    
    # Species Fractions Logic
    h_idx <- grep("^ENSG", all_genes); m_idx <- grep("^ENSMUS", all_genes)
    h_f <- if(length(h_idx) > 0 && total_c > 0) sum(cnts[h_idx, , drop=FALSE]) / total_c else 0
    m_f <- if(length(m_idx) > 0 && total_c > 0) sum(cnts[m_idx, , drop=FALSE]) / total_c else 0
    if(sp %in% c("Human", "Xenograft") && m_f == 0) h_f <- 1
    if(sp == "Mouse" && h_f == 0) m_f <- 1
    
    s_ord <- identical(colnames(dds), rownames(SummarizedExperiment::colData(dds)))
    
    # --- 3. SAMPLE LOOP & FLAG GENERATION ---
    for(samp in colnames(dds)) {
      s_row  <- if(!is.null(salmon_stats)) salmon_stats[salmon_stats$Sample == samp & salmon_stats$CleanDir == sp, ] else data.frame()
      st_row <- if(!is.null(star_stats))   star_stats[star_stats$Sample == samp & star_stats$CleanDir == sp, ]     else data.frame()
      
      # Data Extraction
      raw_total  <- if(nrow(s_row) > 0) as.numeric(gsub(",", "", as.character(s_row[1, "Total (#)"]))) else NA
      mapped_pct <- if(nrow(st_row) > 0) as.numeric(st_row[1, "Unique (%)"]) else NA
      
      star_assign <- NA
      if(nrow(st_row) > 0) {
        target_col <- grep("Assigned|Feature", colnames(st_row), value = TRUE, ignore.case = TRUE)
        if(length(target_col) > 0) star_assign <- as.numeric(st_row[1, target_col[1]])
      }
      
      # --- IN-LINE FLAG LOGIC ---
      flag_reads  <- if(!is.na(raw_total) && raw_total < 5e6) "[LOW_READS]" else ""
      flag_mapped <- if(!is.na(mapped_pct) && mapped_pct < 50) "[LOW_MAPPING]" else ""
      flag_assign <- if(!is.na(star_assign) && star_assign < 30) "[LOW_ASSIGN]" else ""
      
      all_rows_list[[length(all_rows_list) + 1]] <- data.frame(
        Fastq_Folder                  = sp,
        Sample_ID                     = samp,
        
        # Reads & Flag
        Raw_Total_Reads               = raw_total,
        Flag_Reads                    = flag_reads,
        
        # Mapping & Flag
        STAR_Genome_Mapped_Percent    = mapped_pct,
        Flag_Mapping                  = flag_mapped,
        
        # Assignment & Flag
        STAR_Feature_Assigned_Percent = star_assign,
        Flag_Assignment               = flag_assign,
        
        # Counts & Fractions
        Number_geneID_txi             = n_txi_val,
        Number_zero_count_txi         = n_zero_txi_val,
        Number_geneID_DDS             = length(all_genes),
        Fraction_Human_DDS            = round(h_f, 3),
        Fraction_Mouse_DDS             = round(m_f, 3),
        Lib_Size_Million_DDS          = if(samp %in% names(lib_sizes)) round(lib_sizes[[samp]]/1e6, 2) else NA,
        ColData_Order_OK              = s_ord,
        Salmon_Purity                 = if(sp == "Xenograft") 100 else if(nrow(s_row) > 0) s_row[1, "Purity (%)"] else NA,
        stringsAsFactors = FALSE
      )
    }
  }
  
  final_df <- do.call(rbind, all_rows_list)
  write.csv(final_df, "Sample_Master_Summary.csv", row.names = FALSE)
  cat("\n>>> Audit Complete! Results saved to Deep_Audit_Flagged.csv\n")
  print(final_df)

  # --- Visualization ---
  try({
    plot_df <- summary_stats %>%
      dplyr::filter(Fastq_Folder %in% c("Human", "Mouse")) %>%
      dplyr::select(Sample_ID, Fastq_Folder, Salmon_Purity_Percent) %>%
      tidyr::pivot_wider(names_from = Fastq_Folder, values_from = Salmon_Purity_Percent) %>%
      dplyr::mutate(Unmapped_Other = 100 - (replace_na(Human, 0) + replace_na(Mouse, 0))) %>%
      tidyr::pivot_longer(cols = -Sample_ID, names_to = "Source", values_to = "Percentage")
    
    p <- ggplot(plot_df, aes(x = Sample_ID, y = Percentage, fill = Source)) +
      geom_bar(stat = "identity") +
      scale_fill_manual(values = c("Human"="#2c7bb6", "Mouse"="#d7191c", "Unmapped_Other"="#999999")) +
      theme_minimal() + labs(title="Salmon Composition (Human vs Mouse)", y="Percentage (%)") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    
    ggsave("Salmon_Mapping_Summary.pdf", plot = p, width = 8, height = 5)
  }, silent = TRUE)
  
  return(final_df)
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
  
  deep_audit(
    base_path       = get_arg(1)
  )
}