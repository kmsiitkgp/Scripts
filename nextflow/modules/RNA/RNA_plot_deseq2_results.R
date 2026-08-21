#!/usr/bin/env Rscript

# ---- 📦 1. Load Libraries (Always Runs) ----

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(grid)      # grid.draw for heatmaps
  library(grDevices) # cairo_pdf / dev.off
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

plot_deseq2_results <- function(contrast, res, degs,
                                vst_nonblind    = NULL,
                                metadata        = NULL,
                                col_annotations = NULL,
                                col_cluster_by  = "all",
                                output_dir = NULL) {

  # ═══════════════════════════════════════════════════════════════════════════
  # SECTION 1: Validate Inputs
  # ═══════════════════════════════════════════════════════════════════════════
  # All structural checks fire here before any data is loaded or touched.
  # This follows the same contract as plot_volcano() and plot_heatmap() —
  # fail immediately with a clear message rather than crash deep inside a
  # plotting function with a cryptic error.

  # ── 1a. contrast must be a non-empty string ──────────────────────────────
  if (!is.character(contrast) || length(contrast) != 1 || !nzchar(trimws(contrast))) {
    log_error(sample = "", step = "plot_deseq2_results",
              msg = "'contrast' must be a non-empty string (e.g. 'GroupA - GroupB').")
  }

  # ── 1b. col_annotations must be length >= 1 if provided ──────────────────
  if (!is.null(col_annotations) && length(col_annotations) == 0) {
    log_error(sample = "", step = "plot_deseq2_results",
              msg = "'col_annotations' must be NULL or a non-empty character vector.")
  }

  # ── 1c. vst_nonblind and metadata must both be provided or both NULL ──────
  # vst_nonblind is only used for the heatmap, which also requires metadata
  # to derive sample columns. Providing one without the other is always a
  # mistake — fail early rather than silently skip the heatmap.
  if (xor(is.null(vst_nonblind), is.null(metadata))) {
    log_error(sample = "", step = "plot_deseq2_results",
              msg = paste("'vst_nonblind' and 'metadata' must both be provided or both NULL.",
                          "The heatmap requires both to subset samples for this contrast."))
  }

  # ═══════════════════════════════════════════════════════════════════════════
  # SECTION 2: Load Data
  # ═══════════════════════════════════════════════════════════════════════════

  res      <- load_smart(input_path = res)
  degs_df  <- load_smart(input_path = degs)
  expr_mat <- load_smart(input_path = vst_nonblind)
  metadata <- load_smart(input_path = metadata)

  # ═══════════════════════════════════════════════════════════════════════════
  # SECTION 3: Prep expr_mat
  # ═══════════════════════════════════════════════════════════════════════════
  # expr_mat comes in as a data.frame with SYMBOL + gene_id + gene_biotype cols.
  # We need it subsetted to only the samples relevant to this contrast.
  # expr_mat is only used for the DEG heatmap — skipped if not provided.
  # Note: column_to_rownames() is not called here — plot_heatmap() now accepts
  # a data.frame directly with row_id = "SYMBOL".

  if (!is.null(expr_mat) && !is.null(metadata)) {

    # Parse contrast string ("GroupA - GroupB") into group names
    # all.vars() handles any valid R expression, not just " - " format
    target_groups <- all.vars(base::parse(text = contrast))

    # Find which metadata column contains these group labels
    # (could be "Group", "Condition", "Treatment" etc — user-agnostic)
    group_col <- names(metadata)[
      sapply(metadata, function(x) any(target_groups %in% x))
    ][1]

    if (is.na(group_col)) {
      log_error(sample = "", step = "plot_deseq2_results",
                msg = "Could not find contrast groups in any metadata column.")
    }

    relevant_samples <- metadata %>%
      dplyr::filter(.data[[group_col]] %in% target_groups) %>%
      dplyr::pull(Sample_ID)

    # Subset to contrast-relevant samples — keep SYMBOL + sample columns only
    expr_mat <- expr_mat %>%
      dplyr::select(dplyr::any_of(c("SYMBOL", relevant_samples)))

    metadata <- metadata %>% dplyr::filter(Sample_ID %in% relevant_samples)
  }

  # ═══════════════════════════════════════════════════════════════════════════
  # SECTION 4: Set Up Output Path
  # ═══════════════════════════════════════════════════════════════════════════

  # Set global ggplot theme — plot_ma() and plot_volcano() apply
  # theme_publication() internally. theme_set() here ensures any additional
  # ggplot calls in this function inherit the same theme automatically.
  ggplot2::theme_set(theme_publication())

  # Plots are saved in a contrast-specific subfolder — consistent with
  # plot_pathways() and plot_tfs() which use the same structure.

  safe_contrast <- gsub("[^[:alnum:]_-]", "_", contrast)
  output_dir    <- file.path(output_dir %||% ".", safe_contrast)
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

  # ═══════════════════════════════════════════════════════════════════════════
  # SECTION 5: MA Plot
  # ═══════════════════════════════════════════════════════════════════════════
  # res is a DESeqResults object from DESeq2::results(dds, contrast = ...).
  # plot_ma() accepts DESeqResults or data.frame — dds is not a valid input.

  plot_ma(res        = res,
          output_dir = output_dir)

  # ═══════════════════════════════════════════════════════════════════════════
  # SECTION 6: Volcano Plot
  # ═══════════════════════════════════════════════════════════════════════════

  plot_volcano(df         = degs_df,
               x_col      = "log2FoldChange",
               y_col      = "padj",
               label_col  = "SYMBOL",
               output_dir = output_dir)

  # ═══════════════════════════════════════════════════════════════════════════
  # SECTION 7: DEG Heatmap
  # ═══════════════════════════════════════════════════════════════════════════
  # Skipped entirely if vst_nonblind or metadata were not provided.
  # sig_genes is intersected with expr_mat$SYMBOL rather than rownames() since
  # expr_mat is now a data.frame with a SYMBOL column, not a matrix.

  if (!is.null(expr_mat) && !is.null(metadata)) {

    sig_genes <- degs_df %>%
      dplyr::filter(padj <= 0.05) %>%
      dplyr::pull(SYMBOL) %>%
      base::intersect(expr_mat$SYMBOL)

    if (length(sig_genes) >= 2) {

      # Build metadata_col for annotation sidebar — subset to annotation columns
      # + Sample_ID only. NULL when no col_annotations requested.
      metadata_col <- NULL
      if (!is.null(col_annotations)) {
        metadata_col <- metadata %>%
          dplyr::select(dplyr::any_of(col_annotations), Sample_ID)
      }

      # Why save and restore par()? plot_heatmap() sets par() for font family
      # but does not restore it afterward. par(old_par) prevents the font
      # change from leaking into any plots drawn after this heatmap.
      old_par <- par(family = "Helvetica")
      ph <- plot_heatmap(
        df               = expr_mat %>% dplyr::filter(SYMBOL %in% sig_genes),
        row_id           = "SYMBOL",
        col_id           = "Sample_ID",
        data_type        = "counts_log",
        label_genes      = NULL,
        metadata_col     = metadata_col,
        metadata_row     = NULL,
        col_annotations  = col_annotations,
        row_annotations  = NULL,
        col_gap_by       = NULL,
        row_gap_by       = NULL,
        col_cluster_by   = col_cluster_by,
        row_cluster_by   = "all",
        plot_title       = contrast,
        heatmap_palette  = "rdbu",
        border_color     = NA,
        show_expr_legend = TRUE
      )
      par(old_par)

      if (!is.null(ph)) {
        grDevices::cairo_pdf(filename = file.path(output_dir, "Heatmap.pdf"),
                             width    = 10,
                             height   = 11.5,
                             onefile  = TRUE)
        grid::grid.draw(x = ph$ph$gtable)
        grDevices::dev.off()
      }

    } else {
      log_warn(sample = "", step = "plot_deseq2_results",
               msg = glue::glue(
                 "Fewer than 2 significant genes (padj <= 0.05) found in expr_mat ",
                 "for contrast '{contrast}'. Skipping heatmap."
               ))
    }
  }

  # ═══════════════════════════════════════════════════════════════════════════
  # SECTION 8: Log and Return
  # ═══════════════════════════════════════════════════════════════════════════

  log_info(sample = "", step = "plot_deseq2_results",
           msg    = glue::glue("DESeq2 plots saved to: '{output_dir}'"))

  return(invisible(NULL))
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

  # Split comma-separated col_annotations string into a character vector
  raw_col_ann  <- get_arg(6)
  col_ann_list <- if (!is.null(raw_col_ann)) trimws(strsplit(raw_col_ann, ",")[[1]]) else NULL

  plot_deseq2_results(
    contrast        = get_arg(1),
    res             = get_arg(2),
    degs            = get_arg(3),
    vst_nonblind    = get_arg(4),
    metadata        = get_arg(5),
    col_annotations = col_ann_list,
    col_cluster_by  = get_arg(7, default = "all"),
    output_dir      = get_arg(8, default = ".")
  )
}
