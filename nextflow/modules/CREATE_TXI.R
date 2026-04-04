#!/usr/bin/env Rscript

# ---- 📦 1. Load Libraries (Always Runs) ----

suppressPackageStartupMessages({
  library(dplyr)
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

create_txi <- function(species, salmon_dir, output_dir, tx2gene_csv_path = NULL,
                       ensembl_assembly = NULL, ensembl_release = NULL) {

  # ═══════════════════════════════════════════════════════════════════════════
  # SECTION 1: Resolve Species Name
  # ═══════════════════════════════════════════════════════════════════════════
  # Normalise the species string to a formal binomial name before any
  # annotation queries. Accepts common names, abbreviations, or underscored
  # strings — e.g. "Homo_sapiens", "human", "HSAP" all resolve to "Homo sapiens".
  # tolower() + gsub("_", " ") is applied first so pattern matching is
  # case-insensitive and underscore-agnostic.

  clean_species <- gsub(pattern = "_", replacement = " ", x = tolower(species))

  if (grepl(pattern = "human|homo|hsap", x = clean_species)) {
    formal_species <- "Homo sapiens"

  } else if (grepl(pattern = "mouse|mus|mmus", x = clean_species)) {
    formal_species <- "Mus musculus"

  } else if (grepl(pattern = "xeno", x = clean_species)) {
    # Xenograft experiments contain both human (graft) and mouse (host) reads.
    # We annotate against human because the graft signal is what we quantify;
    # mouse reads are filtered out at the alignment stage by Salmon.
    formal_species <- "Homo sapiens"
    log_warn(sample = species, step = "create_txi",
             msg = "Xeno detected: defaulting to Homo sapiens graft annotation.")

  } else {
    formal_species <- clean_species    # Unknown species — pass through as-is
  }

  log_info(sample = species, step = "create_txi",
           msg = paste("Species resolved to:", formal_species))

  # ═══════════════════════════════════════════════════════════════════════════
  # SECTION 2: Fetch or Load tx2gene Annotation
  # ═══════════════════════════════════════════════════════════════════════════
  # tximport needs a tx2gene table to aggregate transcript-level quantifications
  # to gene level. Only the first two columns (transcript_id, gene_id) are
  # needed here — extra columns like gene_name and gene_biotype are used by
  # downstream tools (create_dds, process_counts) but not by tximport itself.
  #
  # Two sources, in priority order:
  #   Case 1: User-provided CSV — fastest and guaranteed to match the Salmon index
  #   Case 2: AnnotationHub query — convenient fallback when no CSV is available

  if (!is.null(tx2gene_csv_path) && file.exists(tx2gene_csv_path)) {

    # ── 2a. Load from CSV ─────────────────────────────────────────────────
    log_info(sample = species, step = "create_txi",
             msg = glue::glue("Loading tx2gene from CSV: {basename(tx2gene_csv_path)}"))

    # dplyr::select(1:2) takes the first two columns regardless of their names —
    # robust to minor column naming differences between CSV sources
    tx2gene <- load_smart(tx2gene_csv_path) %>% dplyr::select(1:2)

  } else {

    # ── 2b. Query AnnotationHub ───────────────────────────────────────────
    log_info(sample = species, step = "create_txi",
             msg = "No tx2gene CSV provided. Querying AnnotationHub...")

    hub    <- AnnotationHub::AnnotationHub()
    hub_db <- AnnotationHub::query(x           = hub,
                                   pattern     = c("EnsDb", formal_species),
                                   ignore.case = TRUE)

    df_hub      <- as.data.frame(mcols(hub_db))
    filtered_db <- df_hub

    # Narrow by assembly (e.g. "GRCh38") and/or release (e.g. "113") if provided.
    # Both filters are optional — if neither is provided, all versions for the
    # species are candidates and the most recently added is used.
    if (!is.null(ensembl_assembly)) {
      filtered_db <- filtered_db %>% dplyr::filter(genome == ensembl_assembly)
    }
    if (!is.null(ensembl_release)) {
      filtered_db <- filtered_db %>%
        dplyr::filter(grepl(paste("Ensembl", ensembl_release), title))
    }

    if (nrow(filtered_db) > 0) {
      # Use the most recently added match in case of duplicates for the same version
      hub_id <- filtered_db %>%
        dplyr::arrange(desc(rdatadateadded)) %>%
        dplyr::slice_head(n = 1) %>%
        rownames()

      log_info(sample = species, step = "create_txi",
               msg = glue::glue("Using: {filtered_db[hub_id, 'title']}"))

    } else {
      # Requested assembly/release not found — warn and fall back to latest.
      # The user should verify this version matches the Salmon index used for
      # quantification — a mismatch causes incorrect transcript-to-gene mapping.
      hub_id <- df_hub %>%
        dplyr::arrange(desc(rdatadateadded)) %>%
        dplyr::slice_head(n = 1) %>%
        rownames()

      log_warn(sample = species, step = "create_txi",
               msg = glue::glue("Requested version not found. ",
                                "Falling back to latest: {df_hub[hub_id, 'title']}. ",
                                "Verify this matches your Salmon index!"))
    }

    ensdb <- hub_db[[hub_id]]

    # Extract only tx_id and gene_id — tximport needs nothing else
    tx2gene <- GenomicFeatures::transcripts(x = ensdb) %>%
      as.data.frame() %>%
      dplyr::select(tx_id, gene_id)
  }

  # ═══════════════════════════════════════════════════════════════════════════
  # SECTION 3: Locate Salmon quant.sf Files
  # ═══════════════════════════════════════════════════════════════════════════
  # Salmon writes one quant.sf per sample into its own subdirectory. We
  # recursively search salmon_dir for all quant.sf files and name each by its
  # parent directory, which by convention is the sample ID set during alignment.
  #
  # Why make.names()? Sample IDs with spaces or special characters (e.g. "Sample-1")
  # would cause issues in downstream R operations. make.names() sanitises them
  # to valid R names consistently.

  quant_files <- list.files(path       = salmon_dir,
                            pattern    = "quant.sf$",
                            recursive  = TRUE,
                            full.names = TRUE)

  if (length(quant_files) == 0) {
    log_error(sample = "", step = "create_txi",
              msg = glue::glue("No quant.sf files found under: {salmon_dir}"))
  }

  names(quant_files) <- make.names(basename(dirname(quant_files)))

  log_info(sample = "", step = "create_txi",
           msg = glue::glue("Found {length(quant_files)} quant.sf file(s)."))

  # ═══════════════════════════════════════════════════════════════════════════
  # SECTION 4: Import with tximport
  # ═══════════════════════════════════════════════════════════════════════════
  # tximport aggregates transcript-level quantifications to gene level using
  # the tx2gene mapping. It also computes TPM (abundance) and scales counts to
  # account for transcript length differences across samples.
  #
  # Why ignoreTxVersion = TRUE?
  # Salmon transcript IDs sometimes include version suffixes (e.g. ENST00000456328.2).
  # The tx2gene table from Ensembl/AnnotationHub typically uses unversioned IDs
  # (ENST00000456328). ignoreTxVersion strips the suffix before matching so no
  # transcripts are lost due to version number mismatches.

  txi <- tximport::tximport(files           = quant_files,
                            type            = "salmon",
                            tx2gene         = tx2gene,
                            ignoreTxVersion = TRUE)

  # ═══════════════════════════════════════════════════════════════════════════
  # SECTION 5: Save and Return
  # ═══════════════════════════════════════════════════════════════════════════
  # Three outputs:
  #   Salmon_txi.rds         : full txi object — passed directly to create_dds()
  #   Salmon_Gene_Counts.txt : raw gene-level counts — useful for manual QC
  #   Salmon_TPM_Values.txt  : TPM abundance values — useful for cross-sample
  #                            expression comparisons (not used for DESeq2)

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

  log_info(sample = "", step = "create_txi",
           msg    = glue::glue("txi saved to: {file.path(output_dir, 'Salmon_txi.rds')}"))

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
