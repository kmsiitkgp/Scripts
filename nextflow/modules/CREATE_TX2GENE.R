#!/usr/bin/env Rscript

# ---- 📦 1. Load Libraries (Always Runs) ----

suppressPackageStartupMessages({
  library(dplyr)      # Provides the %>% operator
  library(readr)
  library(tibble)
  library(AnnotationHub)
  library(rtracklayer)
  library(GenomicFeatures)
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

create_tx2gene <- function(output_dir, gtf_file_path = NULL, species = NULL, ensembl_assembly = NULL, ensembl_release = NULL) {
  
  # ---- MODE A: USER PROVIDED GTF (Highest Priority) ----
  
  if (!is.null(gtf_file_path) && file.exists(gtf_file_path)) {
    log_info(sample = "", step = "create_tx2gene", msg = "MODE: GTF Parsing (Local file)")
    
    # Load GTF
    gtf_data <- rtracklayer::import(gtf_file_path, 
                                    feature.type = "transcript", 
                                    colnames = c("type", "transcript_id", "gene_id", "gene_name", "gene_biotype"))
    
    # Retain only transcripts
    tx <- gtf_data[gtf_data$type == "transcript"]
    
    # Safely extract columns
    get_col <- function(obj, col_name) {
      if (col_name %in% colnames(mcols(obj))) {
        return(mcols(obj)[[col_name]])
      } else {
        return(rep(NA, length(obj))) # Return a vector of NAs the same length as the object
      }
    }
    
    tx2gene <- data.frame(
      transcript_id = get_col(tx, "transcript_id"),
      gene_id       = get_col(tx, "gene_id"),
      gene_name     = get_col(tx, "gene_name"),
      gene_biotype  = get_col(tx, "gene_biotype"),
      stringsAsFactors = FALSE
    )
    
    ### Extract file name
    
    # Get filename: "Mus_musculus.GRCm39.115.gtf"
    base_name <- basename(gtf_file_path)
    
    # Strip extension and species prefix
    # This removes everything before the first dot and the .gtf at the end
    ver_name <- gsub("^[^.]+\\.|\\.gtf$", "", base_name)
    
    csv_name <- paste0("tx2gene_", ver_name, ".csv")
    
  } 
  
  # ---- MODE B: AnnotationHub (Web Query) ----
  
  if (is.null(gtf_file_path) && !is.null(species)) {
    log_info(sample = "", step = "create_tx2gene", msg = "MODE: AnnotationHub (No GTF provided)")
    
    ### Normalize and Resolve Species
    
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
               step = "create_tx2gene", 
               msg = "Xeno detected: Defaulting to Homo sapiens graft.")
      
    } else {
      # Fallback: if it's something else, just use the cleaned string
      formal_species <- clean_species
    }
    
    log_info(sample = species, 
             step = "create_tx2gene", 
             msg = paste("Resolved to:", formal_species))
    
    ### Fetch Annotation
    
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
               step = "create_tx2gene", 
               msg = glue::glue("Found matching Ensembl Db: {filtered_db[hub_id, 'title']}"))
      
    } else {
      # ABSOLUTE FALLBACK: Get the newest available for this species
      hub_id <- df_hub %>%
        dplyr::arrange(desc(rdatadateadded)) %>%
        head(1) %>%
        rownames()
      
      latest_title <- df_hub[hub_id, "title"]
      log_warn(sample = species, step = "create_tx2gene", 
               msg = glue::glue("Requested version not found. Falling back to latest available: {latest_title}"))
    }
    
    # Download the appropriate Ensembl database
    ensdb <- hub_db[[hub_id]]
    
    # Extract transcript and gene info
    tx2gene <- GenomicFeatures::transcripts(ensdb) %>%
      as.data.frame() %>%
      dplyr::select(tx_id, gene_id, gene_name, gene_biotype) %>%
      dplyr::rename(transcript_id = tx_id)
    
    ### Extract file name
    csv_name <- paste0("tx2gene_", species, "_reference.csv")

  } 

  if (!exists("tx2gene")) {
    log_error(sample = "", 
              step = "create_tx2gene", 
              msg = "CRITICAL ERROR: You must provide either a 'gtf_file_path' or a 'species' name.")
  }
  
  # ---- Format tx2gene ----
  
  # Remove exact duplicates
  tx2gene <- unique(tx2gene)
  
  # Retain ONLY rows that have both transcript & gene ID
  tx2gene <- tx2gene[!is.na(tx2gene$transcript_id) & !is.na(tx2gene$gene_id), ]
  
  # If protein-coding gene has gene_id but no gene_name, use gene_id
  # For all other biotypes without gene_name, keep as NA
  tx2gene <- tx2gene %>%
    mutate(gene_name = ifelse((is.na(gene_name) | gene_name == "") & 
                                (gene_biotype == "protein_coding"), gene_id, gene_name))
  
  # ---- 💾 Save Results ----
  
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  csv_file <- file.path(output_dir, csv_name)

  readr::write_csv(x = tx2gene, file = csv_file)

  # ---- 🪵 Log Output and Return ----
  
  log_info(sample = "", 
           step   = "create_tx2gene", 
           msg    = glue::glue("tx2gene created successfully. Saved to: '{csv_file}'"))
  
  return(invisible(tx2gene))
}

# ---- 🚀 4. Smart Execution (ONLY runs in Nextflow) ----

if (!interactive()) {
  
  get_arg <- function(idx) {
    if (idx > length(args)) return(NULL) # Safety if fewer args provided
    val <- args[idx]
    if (is.na(val) || val == "" || val == "null") return(NULL)
    return(val)
  }
  
  create_tx2gene(
    output_dir       = get_arg(1),
    gtf_file_path    = get_arg(2),
    species          = get_arg(3),
    ensembl_assembly = get_arg(4), 
    ensembl_release  = get_arg(5)
  )
}