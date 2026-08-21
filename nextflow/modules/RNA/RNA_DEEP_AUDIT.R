#!/usr/bin/env Rscript

# ---- 📦 1. Load Libraries (Always Runs) ----

suppressPackageStartupMessages({
  library(dplyr)
  library(DESeq2)
  library(SummarizedExperiment)
  library(ggplot2)
  library(tidyr)
  library(glue)       # glue::glue for log messages
  library(grDevices)  # cairo_pdf / dev.off
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

rna_deep_audit <- function(base_path) {

  log_info(sample = "", step = "deep_audit",
           msg = glue::glue("Starting deep audit in: '{base_path}'"))

  # ═══════════════════════════════════════════════════════════════════════════
  # SECTION 1: Load Stats Collector Summaries
  # ═══════════════════════════════════════════════════════════════════════════
  # Both summary files are optional — DEEP_AUDIT can run in de-only mode
  # where no Salmon/STAR outputs exist. Missing files produce NA columns
  # rather than errors so the audit still completes with whatever is available.

  found_salmon_summary <- list.files(base_path,
                                     pattern   = "Salmon_Master_Summary\\.txt$",
                                     recursive = TRUE,
                                     full.names = TRUE)
  found_star_summary   <- list.files(base_path,
                                     pattern   = "STAR_Master_Summary\\.txt$",
                                     recursive = TRUE,
                                     full.names = TRUE)

  salmon_file <- if (length(found_salmon_summary) > 0) found_salmon_summary[1] else NULL
  star_file   <- if (length(found_star_summary)   > 0) found_star_summary[1]   else NULL

  salmon_stats <- if (!is.null(salmon_file)) {
    log_info(sample = "", step = "deep_audit",
             msg = glue::glue("Found Salmon summary: '{basename(salmon_file)}'"))
    read.delim(salmon_file, check.names = FALSE, stringsAsFactors = FALSE)
  } else {
    log_warn(sample = "", step = "deep_audit",
             msg = "No Salmon_Master_Summary.txt found — Salmon columns will be NA.")
    NULL
  }

  star_stats <- if (!is.null(star_file)) {
    log_info(sample = "", step = "deep_audit",
             msg = glue::glue("Found STAR summary: '{basename(star_file)}'"))
    read.delim(star_file, check.names = FALSE, stringsAsFactors = FALSE)
  } else {
    log_warn(sample = "", step = "deep_audit",
             msg = "No STAR_Master_Summary.txt found — STAR columns will be NA.")
    NULL
  }

  # ── CleanDir derivation ───────────────────────────────────────────────────
  # Directory_ID in both summary files encodes the species subfolder as the
  # first path component (e.g. "Human/03.Salmon", "Mouse/gene_counts").
  # dirname() strips the tool subfolder, basename() extracts species name.
  if (!is.null(salmon_stats)) salmon_stats$CleanDir <- basename(dirname(salmon_stats$Directory_ID))
  if (!is.null(star_stats))   star_stats$CleanDir   <- basename(dirname(star_stats$Directory_ID))

  # ═══════════════════════════════════════════════════════════════════════════
  # SECTION 2: Discover Species Folders
  # ═══════════════════════════════════════════════════════════════════════════
  # Priority order: Xenograft first so it appears at the top of the output.
  # Only folders that actually exist in base_path are processed.

  potential_species <- c("Xenograft", "Human", "Mouse")
  species_list      <- potential_species[dir.exists(file.path(base_path, potential_species))]

  if (length(species_list) == 0) {
    log_error(sample = "", step = "deep_audit",
              msg = glue::glue(
                "No species folders (Xenograft / Human / Mouse) found in '{base_path}'. ",
                "Check that base_path points to the correct project directory."
              ))
  }

  log_info(sample = "", step = "deep_audit",
           msg = glue::glue("Species folders found: {paste(species_list, collapse = ', ')}"))

  # ═══════════════════════════════════════════════════════════════════════════
  # SECTION 3: Per-Species Loop — load dds, extract metrics, flag samples
  # ═══════════════════════════════════════════════════════════════════════════

  all_rows_list <- list()

  for (sp in species_list) {

    log_info(sample = "", step = "deep_audit",
             msg = glue::glue("Processing species: '{sp}'"))

    sp_dir <- file.path(base_path, sp)

    # ── 3a. Locate dds and txi RDS files ────────────────────────────────────
    # list.files() is used rather than a fixed path so the audit is robust to
    # minor output directory changes between pipeline versions.
    found_dds <- list.files(sp_dir, pattern = "dds\\.rds$",
                            recursive = TRUE, full.names = TRUE)
    found_txi <- list.files(sp_dir, pattern = "txi\\.rds$",
                            recursive = TRUE, full.names = TRUE)

    dds_path <- if (length(found_dds) > 0) found_dds[1] else NULL
    txi_path <- if (length(found_txi) > 0) found_txi[1] else NULL

    if (is.null(dds_path)) {
      log_warn(sample = "", step = "deep_audit",
               msg = glue::glue("No dds.rds found in '{sp}' — skipping species."))
      next
    }

    # ── 3b. Load DESeqDataSet ────────────────────────────────────────────────
    dds       <- readRDS(dds_path)
    all_genes <- rownames(dds)
    cnts      <- DESeq2::counts(dds)
    lib_sizes <- colSums(cnts)
    total_c   <- sum(cnts)

    # ── 3c. Load tximport object ─────────────────────────────────────────────
    n_txi_val <- NA_integer_; n_zero_txi_val <- NA_integer_
    if (!is.null(txi_path)) {
      txi_obj    <- readRDS(txi_path)
      counts_mtx <- if (is.list(txi_obj) && "counts" %in% names(txi_obj)) txi_obj$counts else txi_obj
      if (!is.null(counts_mtx)) {
        n_txi_val      <- nrow(counts_mtx)
        n_zero_txi_val <- sum(rowSums(counts_mtx) == 0)
      }
    }

    # ── 3d. Species fraction logic ───────────────────────────────────────────
    # Fraction of counts attributable to Human (ENSG*) vs Mouse (ENSMUSG*)
    # gene IDs. For single-species runs one fraction will always be 1.0 — the
    # fallback handles this so fractions are always interpretable.
    h_idx <- grep("^ENSG",   all_genes)
    m_idx <- grep("^ENSMUSG", all_genes)

    h_f <- if (length(h_idx) > 0 && total_c > 0) sum(cnts[h_idx,  , drop = FALSE]) / total_c else 0
    m_f <- if (length(m_idx) > 0 && total_c > 0) sum(cnts[m_idx,  , drop = FALSE]) / total_c else 0

    # Single-species fallback — if no mouse genes found in a Human/Xenograft
    # dds, human fraction is by definition 1.0
    if (sp %in% c("Human", "Xenograft") && m_f == 0) h_f <- 1
    if (sp == "Mouse" && h_f == 0) m_f <- 1

    # ── 3e. ColData order check ──────────────────────────────────────────────
    # DESeq2 requires colnames(dds) and rownames(colData) to match exactly.
    # A mismatch here means sample-level metadata is misaligned with counts —
    # a silent but catastrophic error in differential expression.
    s_ord <- identical(colnames(dds), rownames(SummarizedExperiment::colData(dds)))

    if (!s_ord) {
      log_warn(sample = "", step = "deep_audit",
               msg = glue::glue(
                 "ColData order MISMATCH in '{sp}' dds — sample metadata may be ",
                 "misaligned with count matrix. Verify metadata row order."
               ))
    }

    # ── 3f. Per-sample row construction ─────────────────────────────────────
    for (samp in colnames(dds)) {

      # Match this sample to its row in each stats summary
      s_row  <- if (!is.null(salmon_stats)) {
        salmon_stats[salmon_stats$Sample == samp & salmon_stats$CleanDir == sp, ]
      } else data.frame()

      st_row <- if (!is.null(star_stats)) {
        star_stats[star_stats$Sample == samp & star_stats$CleanDir == sp, ]
      } else data.frame()

      # ── Extract numeric metrics ──────────────────────────────────────────
      raw_total  <- if (nrow(s_row)  > 0) as.numeric(gsub(",", "", as.character(s_row[1,  "Total (#)"])))  else NA_real_
      mapped_pct <- if (nrow(st_row) > 0) as.numeric(st_row[1, "Unique (%)"])                              else NA_real_
      purity     <- if (nrow(s_row)  > 0) as.numeric(s_row[1,  "Purity (%)"])                              else NA_real_

      # STAR feature assignment — column name varies by strand (U/F/R)
      # grep finds whichever "Feature" or "Assigned" column is present first
      star_assign <- NA_real_
      if (nrow(st_row) > 0) {
        target_col <- grep("Feature U|Assigned", colnames(st_row),
                           value = TRUE, ignore.case = TRUE)
        if (length(target_col) > 0) {
          star_assign <- as.numeric(st_row[1, target_col[1]])
        }
      }

      # ── Flag generation ──────────────────────────────────────────────────
      # Thresholds: <5M reads is low depth; <50% mapping suggests QC issue;
      # <30% feature assignment suggests strand mismatch or poor annotation.
      flag_reads  <- if (!is.na(raw_total)   && raw_total   < 5e6) "[LOW_READS]"   else ""
      flag_mapped <- if (!is.na(mapped_pct)  && mapped_pct  < 50)  "[LOW_MAPPING]" else ""
      flag_assign <- if (!is.na(star_assign) && star_assign < 30)  "[LOW_ASSIGN]"  else ""

      all_rows_list[[length(all_rows_list) + 1]] <- data.frame(
        Species                       = sp,
        Sample_ID                     = samp,

        # ── Reads ────────────────────────────────────────────────────────
        Raw_Total_Reads               = raw_total,
        Flag_Reads                    = flag_reads,

        # ── Mapping ──────────────────────────────────────────────────────
        STAR_Unique_Mapped_Pct        = mapped_pct,
        Flag_Mapping                  = flag_mapped,

        # ── Assignment ───────────────────────────────────────────────────
        STAR_Feature_Assigned_Pct     = star_assign,
        Flag_Assignment               = flag_assign,

        # ── Counts ───────────────────────────────────────────────────────
        N_Genes_txi                   = n_txi_val,
        N_Zero_Count_Genes_txi        = n_zero_txi_val,
        N_Genes_dds                   = length(all_genes),
        Lib_Size_Million_dds          = if (samp %in% names(lib_sizes)) round(lib_sizes[[samp]] / 1e6, 2) else NA_real_,

        # ── Species fractions ────────────────────────────────────────────
        Fraction_Human_dds            = round(h_f, 3),
        Fraction_Mouse_dds            = round(m_f, 3),

        # ── QC flags ─────────────────────────────────────────────────────
        ColData_Order_OK              = s_ord,

        # ── Xenograft purity ─────────────────────────────────────────────
        # Purity (%) from salmon_stats_collector: ratio of reads in this
        # species folder to the total reads in the Xenograft (unsorted) run.
        # Meaningful for Human and Mouse folders in Xenograft experiments.
        # NA for standard single-species experiments.
        Salmon_Purity_Pct             = purity,

        stringsAsFactors = FALSE
      )
    }
  }

  # ═══════════════════════════════════════════════════════════════════════════
  # SECTION 4: Assemble and Save Summary Table
  # ═══════════════════════════════════════════════════════════════════════════

  if (length(all_rows_list) == 0) {
    log_error(sample = "", step = "deep_audit",
              msg = paste("No samples processed — no dds.rds found in any species folder.",
                          "Check that DE steps completed successfully before running deep_audit."))
  }

  final_df <- do.call(rbind, all_rows_list)

  write.csv(final_df, "Sample_Master_Summary.csv", row.names = FALSE)

  log_info(sample = "", step = "deep_audit",
           msg = glue::glue("Audit complete. Summary saved to 'Sample_Master_Summary.csv' ",
                            "({nrow(final_df)} samples across {length(species_list)} species)."))

  # ═══════════════════════════════════════════════════════════════════════════
  # SECTION 5: Xenograft Purity Visualization
  # ═══════════════════════════════════════════════════════════════════════════
  # A stacked bar plot showing Human vs Mouse vs unassigned read fractions
  # per sample — only meaningful for Xenograft experiments where both Human
  # and Mouse species folders exist with Salmon_Purity_Pct values.

  is_xenograft <- all(c("Human", "Mouse") %in% unique(final_df$Species))

  if (is_xenograft && any(!is.na(final_df$Salmon_Purity_Pct))) {

    plot_df <- final_df %>%
      dplyr::filter(Species %in% c("Human", "Mouse")) %>%
      dplyr::select(Sample_ID, Species, Salmon_Purity_Pct) %>%
      tidyr::pivot_wider(
        names_from  = Species,
        values_from = Salmon_Purity_Pct
      ) %>%
      dplyr::mutate(
        Human             = tidyr::replace_na(Human, 0),
        Mouse             = tidyr::replace_na(Mouse, 0),
        Unassigned        = pmax(0, 100 - Human - Mouse)
      ) %>%
      tidyr::pivot_longer(
        cols      = c(Human, Mouse, Unassigned),
        names_to  = "Source",
        values_to = "Percentage"
      ) %>%
      dplyr::mutate(Source = factor(Source, levels = c("Human", "Mouse", "Unassigned")))

    p <- ggplot2::ggplot(
      data    = plot_df,
      mapping = ggplot2::aes(x = Sample_ID, y = Percentage, fill = Source)
    ) +
      ggplot2::geom_col(width = 0.75) +
      ggplot2::scale_fill_manual(
        values = c("Human" = "#2c7bb6", "Mouse" = "#d7191c", "Unassigned" = "#999999")
      ) +
      ggplot2::scale_y_continuous(
        limits = c(0, 100),
        breaks = seq(0, 100, by = 25),
        labels = function(x) paste0(x, "%")
      ) +
      ggplot2::coord_cartesian(clip = "off") +
      theme_publication() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
      ggplot2::labs(
        title    = "Xenograft Read Composition per Sample",
        subtitle = "Purity = fraction of reads assigned to each species by Salmon",
        x        = NULL,
        y        = "Percentage of Reads (%)",
        fill     = "Species"
      )

    grDevices::cairo_pdf("Salmon_Mapping_Summary.pdf", width = 10, height = 6)
    print(p)
    grDevices::dev.off()

    log_info(sample = "", step = "deep_audit",
             msg = "Xenograft purity plot saved to 'Salmon_Mapping_Summary.pdf'")

  } else {
    log_info(sample = "", step = "deep_audit",
             msg = "Non-xenograft experiment or no purity data — skipping purity plot.")
  }

  return(invisible(final_df))
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

  rna_deep_audit(
    base_path = get_arg(1)
  )
}
