#!/usr/bin/env Rscript

# ---- 📦 1. Load Libraries (Always Runs) ----

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(purrr)
  library(openxlsx)
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

merge_counts <- function(counts_dir, output_dir) {

  # ═══════════════════════════════════════════════════════════════════════════
  # SECTION 1: Identify Count Files
  # ═══════════════════════════════════════════════════════════════════════════
  # STAR writes one ReadsPerGene.out.tab file per sample during the alignment
  # step (--quantMode GeneCounts). Each file has 4 columns:
  #   col 1: gene_id    — Ensembl gene ID
  #   col 2: unstranded — counts ignoring strand
  #   col 3: forward    — counts on the forward strand
  #   col 4: reverse    — counts on the reverse strand
  # We collect all such files across the counts_dir before processing.

  count_files <- list.files(path       = counts_dir,
                            pattern    = "ReadsPerGene\\.out\\.tab$",
                            full.names = TRUE,
                            recursive  = TRUE)

  if (length(count_files) == 0) {
    log_warn(sample = "", step = "merge_counts",
             msg = glue::glue("No ReadsPerGene.out.tab files found under: {counts_dir}"))
    return(NULL)
  }

  log_info(sample = "", step = "merge_counts",
           msg = glue::glue("Found {length(count_files)} count file(s)."))

  # ═══════════════════════════════════════════════════════════════════════════
  # SECTION 2: Parse Files and Detect Strandedness
  # ═══════════════════════════════════════════════════════════════════════════

  # ── 2a. Define summary rows to exclude ──────────────────────────────────
  # STAR prepends 4 summary rows to every ReadsPerGene.out.tab file before the
  # per-gene counts begin. These rows use special gene IDs (N_unmapped etc.)
  # and report alignment statistics, not gene counts. They must be excluded
  # before strandedness detection — including them would inflate forward/reverse
  # totals with non-gene reads and distort the strand proportion calculation.
  # The double-underscore variants (__no_feature etc.) are HTSeq-style labels
  # that appear in some STAR versions or featureCounts outputs.
  special_counters <- c("N_unmapped", "N_multimapping", "N_noFeature", "N_ambiguous",
                        "__no_feature", "__ambiguous", "__too_low_aQual",
                        "__not_aligned", "__alignment_not_unique", "__assignment_counts")

  all_counts <- list()

  for (count_file in count_files) {

    # ── 2b. Extract sample ID from filename ─────────────────────────────────
    # The regex removes the "ReadsPerGene.out.tab" suffix AND any other dot-
    # separated extension that may precede it (e.g. "Sample1.Aligned.out.tab"
    # becomes "Sample1"). The alternation handles both naming conventions.
    sample_id <- gsub("\\..*$|ReadsPerGene\\.out\\.tab", "", basename(count_file))

    # ── 2c. Read the file ────────────────────────────────────────────────────
    # col_names = FALSE because STAR does not write a header row — the first
    # data line is already a summary counter (N_unmapped). We assign names
    # explicitly via rename() to make downstream column references readable.
    df <- tryCatch({
      readr::read_tsv(file           = count_file,
                      col_names      = FALSE,
                      show_col_types = FALSE) %>%
        dplyr::rename(gene_id    = 1,
                      unstranded = 2,
                      forward    = 3,
                      reverse    = 4)
    }, error = function(e) {
      log_error(sample = sample_id, step = "merge_counts",
                msg = glue::glue("Failed to read file: {e$message}"))
      return(NULL)
    })

    if (is.null(df) || ncol(df) < 4) next

    # Remove STAR summary rows — gene-level counts only from here on
    df <- df %>% dplyr::filter(!(gene_id %in% special_counters))

    # ── 2d. Detect strandedness ──────────────────────────────────────────────
    # Strandedness is inferred from the proportion of reads mapping to the
    # forward strand relative to all stranded reads.
    #
    # Why use fwd/(fwd+rev) rather than fwd/(fwd+rev+unstranded)?
    # The unstranded column is NOT a third category of reads — it is the SUM of
    # forward + reverse. Including it in the denominator would double-count all
    # reads and artificially compress the proportion toward 0.5.
    #
    # Threshold rationale:
    #   prop_fwd > 0.8 : ≥80% forward → library is forward-stranded
    #   prop_fwd < 0.2 : ≤20% forward → library is reverse-stranded (≥80% reverse)
    #   0.4–0.6        : ~50/50 split → genuinely unstranded library
    #   0.2–0.4 or 0.6–0.8 : ambiguous — neither clearly stranded nor unstranded;
    #                         falls back to unstranded with a warning so the run
    #                         completes rather than hard-stopping on one sample

    fwd_sum <- sum(df$forward,  na.rm = TRUE)
    rev_sum <- sum(df$reverse,  na.rm = TRUE)
    total   <- fwd_sum + rev_sum

    if (total == 0) {
      log_error(sample = sample_id, step = "merge_counts",
                msg = "Total stranded counts are zero — sample excluded.")
      next
    }

    prop_fwd <- fwd_sum / total

    if (prop_fwd > 0.8) {
      log_info(sample = sample_id, step = "merge_counts",
               msg = glue::glue("Strandedness: Forward ({round(prop_fwd * 100, 1)}% fwd)"))
      sample_df <- df %>% dplyr::select(ID = gene_id, !!sample_id := forward)

    } else if (prop_fwd < 0.2) {
      log_info(sample = sample_id, step = "merge_counts",
               msg = glue::glue("Strandedness: Reverse ({round((1-prop_fwd) * 100, 1)}% rev)"))
      sample_df <- df %>% dplyr::select(ID = gene_id, !!sample_id := reverse)

    } else {
      # abs(prop_fwd - 0.5) < 0.1 catches the genuinely unstranded case
      # (40–60% forward). Values outside this range (0.2–0.4 or 0.6–0.8) are
      # ambiguous — warn the user but proceed with unstranded rather than crash.
      msg <- if (abs(prop_fwd - 0.5) < 0.1) {
        glue::glue("Strandedness: Unstranded ({round(prop_fwd * 100, 1)}% fwd)")
      } else {
        glue::glue("Strandedness: Ambiguous ({round(prop_fwd * 100, 1)}% fwd). ",
                   "Falling back to unstranded — verify library prep protocol.")
      }
      log_info(sample = sample_id, step = "merge_counts", msg = msg)
      sample_df <- df %>% dplyr::select(ID = gene_id, !!sample_id := unstranded)
    }

    all_counts[[sample_id]] <- sample_df
  }

  if (length(all_counts) == 0) {
    log_error(sample = "", step = "merge_counts",
              msg = "No valid samples remained after parsing. Check input files.")
  }

  # ═══════════════════════════════════════════════════════════════════════════
  # SECTION 3: Build Count Matrix
  # ═══════════════════════════════════════════════════════════════════════════

  # ── 3a. Merge per-sample dataframes ─────────────────────────────────────
  # purrr::reduce() applies full_join iteratively across all sample dataframes,
  # joining on the shared "ID" column.
  # Why full_join and not inner_join?
  # inner_join would drop any gene absent from even one sample — in practice
  # this can remove real genes that happen to have zero counts in one sample
  # (e.g. lowly expressed genes). full_join retains all genes and produces NA
  # for absent genes, which we then replace with 0 (correct: gene was present
  # in the annotation but had no detected reads in that sample).
  count_matrix <- purrr::reduce(all_counts, dplyr::full_join, by = "ID")

  # Replace NAs introduced by full_join with 0 — absent gene = zero counts
  count_matrix[is.na(count_matrix)] <- 0

  # ── 3b. Remove all-zero genes and samples ────────────────────────────────
  # All-zero genes: not expressed in any sample — cannot yield DE results and
  # inflate the multiple testing burden. rowSums on all columns except ID.
  count_matrix <- count_matrix %>%
    dplyr::filter(rowSums(dplyr::select(., -ID), na.rm = TRUE) > 0)

  # All-zero samples: indicate a failed sequencing run or sample mislabelling.
  # The c(TRUE, ...) prepends TRUE for the ID column so it is always kept —
  # without this, the logical vector would be one element too short and R would
  # recycle it incorrectly, potentially dropping the ID column.
  count_matrix <- count_matrix[, c(TRUE, colSums(count_matrix[, -1]) > 0), drop = FALSE]

  log_info(sample = "", step = "merge_counts",
           msg = glue::glue("Matrix dimensions: {nrow(count_matrix)} genes × ",
                            "{ncol(count_matrix) - 1} samples."))

  # ═══════════════════════════════════════════════════════════════════════════
  # SECTION 4: Save and Return
  # ═══════════════════════════════════════════════════════════════════════════

  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

  output_file <- file.path(output_dir, "STAR_Gene_counts.xlsx")

  openxlsx::write.xlsx(x         = list("raw_counts" = count_matrix),
                       file      = output_file,
                       overwrite = TRUE)

  log_info(sample = "", step = "merge_counts",
           msg    = glue::glue("Merged count matrix saved to: '{output_file}'"))

  return(invisible(count_matrix))
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

  merge_counts(
    counts_dir = get_arg(1),
    output_dir = get_arg(2)
  )
}
