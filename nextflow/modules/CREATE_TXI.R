#!/usr/bin/env Rscript

# ---- 📦 1. Load Libraries (Always Runs) ----

suppressPackageStartupMessages({
  library(dplyr)      # Provides the %>% operator
  library(readr)
  library(tibble)
  library(tximport)
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

create_txi <- function(species, salmon_dir, output_dir, tx2gene_csv_path = NULL, ensembl_assembly = NULL, ensembl_release = NULL) {
  
  # ---- 🧬 Normalize and Resolve Species ----
  
  # Remove underscores and convert to lowercase (e.g., "Homo_sapiens" -> "homo sapiens")
  clean_species <- gsub(pattern = "_", replacement = " ", x = tolower(species))
  
  # Define the mapping
  if (grepl(pattern = "human|homo|hsap", x = clean_species)) {
    formal_species <- "Homo sapiens"
    
  } else if (grepl(pattern = "mouse|mus|mmus", x = clean_species)) {
    formal_species <- "Mus musculus"
    
  } else if (grepl(pattern = "xeno", x = clean_species)) {
    formal_species <- "Homo sapiens"   # usually we are interested in graft data
    log_warn(sample = species, 
             step = "create_txi", 
             msg = "Xeno detected: Defaulting to Homo sapiens graft.")
    
  } else {
    # Fallback: if it's something else, just use the cleaned string
    formal_species <- clean_species
  }
  
  log_info(sample = species, 
           step = "create_txi", 
           msg = paste("Resolved to:", formal_species))
  
  # ---- 🔍 Fetch or Load Annotation (tx2gene) ----
  
  if (!is.null(tx2gene_csv_path) && file.exists(tx2gene_csv_path)) {
    # CASE 1: User provided a CSV
    log_info(sample = species,
             step = "create_txi",
             msg = glue::glue("✔ FOUND User-provided CSV: {basename(tx2gene_csv_path)}"))
    
    tx2gene <- readr::read_csv(tx2gene_csv_path, show_col_types = FALSE) %>% 
      dplyr::select(1:2)
    
  } else {
    # CASE 2: No CSV, Fetch from AnnotationHub
    log_info(sample = species, 
             step = "create_txi",
             msg = "✖ MISSING CSV. Querying AnnotationHub...")
    
    # Connect to AnnotationHub 
    hub <- AnnotationHub::AnnotationHub()
    
    # Query for required species
    hub_db <- AnnotationHub::query(x = hub, 
                                   pattern = c("EnsDb", formal_species), 
                                   ignore.case = TRUE)
    
    # Convert mcols to dataframe once for easier filtering
    df_hub <- as.data.frame(mcols(hub_db))
    
    # Identify database ID for annotation
    filtered_db <- df_hub
    
    # Filter by Assembly if provided (e.g., "GRCh38")
    if (!is.null(ensembl_assembly)) {
      filtered_db <- filtered_db %>% dplyr::filter(genome == ensembl_assembly)
    }
    
    # Filter by Release if provided (e.g., "113")
    if (!is.null(ensembl_release)) {
      filtered_db <- filtered_db %>% 
        dplyr::filter(grepl(paste("Ensembl", ensembl_release), title))
    }
    
    # Identify Hub ID with Fallback
    if (nrow(filtered_db) > 0) {
      hub_id <- filtered_db %>%
        # Sort to get the most recent record (just in case of duplicates)
        dplyr::arrange(desc(rdatadateadded)) %>%
        head(1) %>%
        rownames()
      
      log_info(sample = species,
               step = "create_txi", 
               msg = glue::glue("Found matching Ensembl Db: {filtered_db[hub_id, 'title']}"))
      
    } else {
      # ABSOLUTE FALLBACK: Get the newest available for this species
      hub_id <- df_hub %>%
        dplyr::arrange(desc(rdatadateadded)) %>%
        head(1) %>%
        rownames()
      
      latest_title <- df_hub[hub_id, "title"]
      log_warn(sample = species, step = "create_txi", 
               msg = glue::glue("Requested version not found. Falling back to latest available: {latest_title}"))
    }
    
    # Download the appropriate Ensembl database
    ensdb <- hub_db[[hub_id]]
    
    # Extract transcript and gene info
    tx2gene <- GenomicFeatures::transcripts(x = ensdb) %>%
      as.data.frame() %>%
      dplyr::select(tx_id, gene_id)
  }
  
  # ---- 📂 Locate Salmon Files ----
  
  quant_files <- list.files(path       = salmon_dir, 
                            pattern    = "quant.sf$", 
                            recursive  = TRUE, 
                            full.names = TRUE)
  
  # Extract folder names as Sample IDs
  # basename(dirname()) is the cleanest way to get the folder name
  names(quant_files) <- make.names(basename(dirname(quant_files)))
  
  # ---- 📥 Import Data with tximport ----
  
  txi <- tximport::tximport(files           = quant_files, 
                            type            = "salmon", 
                            tx2gene         = tx2gene,  
                            ignoreTxVersion = TRUE)
  
  # ---- 💾 Save Results ----
  
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  txi$counts %>%
    as.data.frame() %>%
    tibble::rownames_to_column(var = "GeneID") %>%
    write.table(file.path(output_dir, "Salmon_Gene_Counts.txt"),
                sep = "\t", quote = FALSE, row.names = FALSE)
  
  txi$abundance %>%
    as.data.frame() %>%
    tibble::rownames_to_column(var = "GeneID") %>%
    write.table(file.path(output_dir, "Salmon_TPM_Values.txt"),
                sep = "\t", quote = FALSE, row.names = FALSE)
  
  saveRDS(txi, file = file.path(output_dir, "Salmon_txi.rds"))
  
  # ---- 🪵 Log Output and Return ----
  
  log_info(sample = "", 
           step   = "create_txi", 
           msg    = "txi created successfully. Saved to: 'Salmon_txi.rds'")
  
  return(invisible(txi))
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

  create_txi(
    species          = get_arg(1),
    salmon_dir       = get_arg(2),
    output_dir       = get_arg(3),
    tx2gene_csv_path = get_arg(4), 
    ensembl_assembly = get_arg(5),
    ensembl_release  = get_arg(6)
  )
}