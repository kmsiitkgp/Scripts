#!/usr/bin/env Rscript

# ---- 📦 1. Load Libraries (Always Runs) ----

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tibble)
  library(AnnotationHub)
  library(rtracklayer)
  library(GenomicFeatures)
})

# ---- 🛠️ 2. Smart Setup (Find & source UTILS.R) ----

get_modules_dir <- function() {
  # 1. Windows dev machine
  if (.Platform$OS.type == "windows") {
    return("C:/Users/kailasamms/OneDrive - Cedars-Sinai Health System/Documents/GitHub/Scripts/nextflow/modules")
  }
  # 2. Interactive Linux / macOS (HPC interactive session)
  if (interactive()) {
    return(file.path(getwd(), "modules"))
  }
  # 3. Non-interactive (Nextflow / Rscript)
  initial.options <- commandArgs(trailingOnly = FALSE)
  file_arg        <- grep("--file=", initial.options, value = TRUE)
  if (length(file_arg) == 0) stop("Cannot detect modules directory path!")
  modules_dir      <- dirname(sub("--file=", "", file_arg))
  return(modules_dir)
}

find_utils_dir <- function(start_dir) {
  d <- normalizePath(start_dir, mustWork = FALSE)
  repeat {
    if (file.exists(file.path(d, "UTILS.R"))) return(d)
    parent <- dirname(d)
    if (parent == d) stop("❌ UTILS.R not found searching upward from: ", start_dir)
    d <- parent
  }
}

modules_dir <- get_modules_dir()
utils_path  <- file.path(find_utils_dir(modules_dir), "UTILS.R")
source(utils_path)

# ---- 🧬 3. Function Definition (Always Runs) ----

salmon_create_tx2gene <- function(output_dir, ref_gtf_path = NULL, species = NULL,
                           ensembl_assembly = NULL, ensembl_release = NULL) {

  # ═══════════════════════════════════════════════════════════════════════════
  # SECTION 1: Build tx2gene Table  —  Mode A or Mode B
  # ═══════════════════════════════════════════════════════════════════════════
  # Two mutually exclusive input modes:
  #   Mode A (GTF)          : parse a local GTF file — fastest, fully reproducible,
  #                           guaranteed to match the genome index used by Salmon
  #   Mode B (AnnotationHub): query Ensembl EnsDb online — convenient when no
  #                           local GTF is available, but requires internet access
  #
  # Mode A takes priority: if ref_gtf_path is provided AND exists, Mode B is
  # never attempted. This prevents silent fallback to a mismatched annotation.

  # ── 1a. Mode A: Parse local GTF ─────────────────────────────────────────
  # Why filter to feature.type = "transcript" upfront?
  # GTF files contain rows for gene, transcript, exon, CDS, UTR etc. We only
  # need transcript-level rows — filtering at import is faster than loading
  # everything and then subsetting in R.

  if (!is.null(ref_gtf_path) && file.exists(ref_gtf_path)) {
    log_info(sample = "", step = "salmon_create_tx2gene", msg = "MODE: GTF Parsing (local file)")

    gtf_data <- rtracklayer::import(ref_gtf_path,
                                    feature.type = "transcript",
                                    colnames     = c("type", "transcript_id", "gene_id",
                                                     "gene_name", "gene_biotype"))

    tx <- gtf_data[gtf_data$type == "transcript"]

    # Safely extract a metadata column, returning NAs if absent.
    # Why not just mcols(obj)[[col_name]] directly? GTF files from different
    # sources (Ensembl, GENCODE, RefSeq) don't all have the same columns —
    # direct access would throw an error on missing columns. get_col() makes
    # the extraction robust to any GTF format.
    get_col <- function(obj, col_name) {
      if (col_name %in% colnames(mcols(obj))) return(mcols(obj)[[col_name]])
      return(rep(NA, length(obj)))
    }

    tx2gene <- data.frame(
      transcript_id = get_col(tx, "transcript_id"),
      gene_id       = get_col(tx, "gene_id"),
      gene_name     = get_col(tx, "gene_name"),
      gene_biotype  = get_col(tx, "gene_biotype"),
      stringsAsFactors = FALSE
    )

    # Derive output filename from GTF filename so it encodes the genome version.
    # e.g. "Mus_musculus.GRCm39.115.gtf" → "tx2gene_GRCm39.115.csv"
    # The regex strips the leading species prefix (everything up to first dot)
    # and the trailing ".gtf" extension, keeping only the assembly + release.
    base_name <- basename(ref_gtf_path)
    ver_name  <- gsub("^[^.]+\\.|\\.gtf$", "", base_name)
    csv_name  <- paste0("tx2gene_", ver_name, ".csv")
  }

  # ── 1b. Mode B: Query AnnotationHub ─────────────────────────────────────
  # Only runs when no GTF path was provided. Requires internet access.
  # Species string is normalised before querying so common abbreviations
  # ("human", "hsap", "Homo_sapiens") all resolve to the formal binomial name
  # that AnnotationHub expects.

  if (is.null(ref_gtf_path) && !is.null(species)) {
    log_info(sample = "", step = "salmon_create_tx2gene", msg = "MODE: AnnotationHub (no GTF provided)")

    # ── 1b-i. Resolve species name ───────────────────────────────────────
    # tolower() + gsub("_", " ") normalises capitalisation and underscores
    # before grepl pattern matching, so "Homo_sapiens", "HUMAN", and "hsap"
    # all match the "human|homo|hsap" pattern.
    clean_species <- gsub(pattern = "_", replacement = " ", x = tolower(species))

    if (grepl(pattern = "human|homo|hsap", x = clean_species)) {
      formal_species <- "Homo sapiens"

    } else if (grepl(pattern = "mouse|mus|mmus", x = clean_species)) {
      formal_species <- "Mus musculus"

    } else if (grepl(pattern = "xeno", x = clean_species)) {
      # Xenograft experiments contain both human (graft) and mouse (host) reads.
      # We annotate against human because the graft signal is what we want to
      # quantify; mouse reads are typically filtered out at the alignment stage.
      formal_species <- "Homo sapiens"
      log_warn(sample = species, step = "salmon_create_tx2gene",
               msg = "Xeno detected: defaulting to Homo sapiens graft annotation.")

    } else {
      formal_species <- clean_species    # Unknown species — pass through as-is
    }

    log_info(sample = species, step = "salmon_create_tx2gene",
             msg = paste("Resolved to:", formal_species))

    # ── 1b-ii. Query AnnotationHub ───────────────────────────────────────
    # Query returns all EnsDb objects matching the species. We then narrow by
    # assembly and/or release if the user specified them, and fall back to the
    # latest available version if the requested one is not found.
    hub    <- AnnotationHub::AnnotationHub()
    hub_db <- AnnotationHub::query(x           = hub,
                                   pattern     = c("EnsDb", formal_species),
                                   ignore.case = TRUE)

    df_hub      <- as.data.frame(mcols(hub_db))
    filtered_db <- df_hub

    # Narrow by assembly (e.g. "GRCh38") and/or release (e.g. "113") if provided
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

      log_info(sample = species, step = "salmon_create_tx2gene",
               msg = glue::glue("Using: {filtered_db[hub_id, 'title']}"))

    } else {
      # Requested assembly/release not found — warn and fall back to latest.
      # This prevents a hard stop but the user should verify the version matches
      # the genome index used by Salmon upstream.
      hub_id <- df_hub %>%
        dplyr::arrange(desc(rdatadateadded)) %>%
        dplyr::slice_head(n = 1) %>%
        rownames()

      log_warn(sample = species, step = "salmon_create_tx2gene",
               msg = glue::glue("Requested version not found. ",
                                "Falling back to latest: {df_hub[hub_id, 'title']}. ",
                                "Verify this matches your Salmon index!"))
    }

    ensdb <- hub_db[[hub_id]]

    tx2gene <- GenomicFeatures::transcripts(ensdb) %>%
      as.data.frame() %>%
      dplyr::select(tx_id, gene_id, gene_name, gene_biotype) %>%
      dplyr::rename(transcript_id = tx_id)

    csv_name <- paste0("tx2gene_", species, ".csv")
  }

  # Hard stop if neither mode produced a tx2gene table
  if (!exists("tx2gene")) {
    log_error(sample = "", step = "salmon_create_tx2gene",
              msg = "Must provide either 'ref_gtf_path' (local GTF) or 'species' (AnnotationHub).")
  }

  # ═══════════════════════════════════════════════════════════════════════════
  # SECTION 2: Format tx2gene Table
  # ═══════════════════════════════════════════════════════════════════════════

  # ── 2a. Deduplicate ──────────────────────────────────────────────────────
  # GTF files can contain duplicate transcript rows (e.g. when the same
  # transcript appears in multiple biotype annotations). unique() collapses
  # these before any downstream filtering.
  tx2gene <- unique(tx2gene)

  # ── 2b. Remove rows missing transcript_id or gene_id ────────────────────
  # tximport requires both columns to map transcripts to genes — rows with
  # either missing will cause silent mismatches or errors during import.
  tx2gene <- tx2gene[!is.na(tx2gene$transcript_id) & !is.na(tx2gene$gene_id), ]

  # ── 2c. Fill missing gene_name for protein-coding genes ─────────────────
  # Some annotation sources leave gene_name blank for non-canonical transcripts
  # even when the gene has a known symbol. For protein_coding genes only, fall
  # back to gene_id (e.g. "ENSG00000...") so downstream tools always have a
  # non-empty label. Non-coding biotypes are left as NA intentionally — they
  # often don't have meaningful gene symbols and forcing a fallback could
  # mislead interpretation.
  tx2gene <- tx2gene %>%
    dplyr::mutate(gene_name = ifelse(
      (is.na(gene_name) | gene_name == "") & (gene_biotype == "protein_coding"),
      gene_id, gene_name
    ))

  # ═══════════════════════════════════════════════════════════════════════════
  # SECTION 3: Save and Return
  # ═══════════════════════════════════════════════════════════════════════════

  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

  csv_file <- file.path(output_dir, csv_name)
  readr::write_csv(x = tx2gene, file = csv_file)

  log_info(sample = "", step = "salmon_create_tx2gene",
           msg    = glue::glue("tx2gene saved to: '{csv_file}'"))

  return(invisible(tx2gene))
}

# ---- 🚀 4. Smart Execution (Nextflow Only) ----

if (sys.nframe() == 0) {
  args <- commandArgs(trailingOnly = TRUE)

  get_arg <- function(idx, default = NULL) {
    if (idx > length(args)) return(default)
    val <- args[idx]
    if (is.na(val) || val == "" || val == "null" || val == "NULL") return(default)
    return(val)
  }

  salmon_create_tx2gene(
    output_dir       = get_arg(1),
    ref_gtf_path     = get_arg(2),
    species          = get_arg(3),
    ensembl_assembly = get_arg(4),
    ensembl_release  = get_arg(5)
  )
}
