#!/usr/bin/env Rscript

# ---- 📦 1. Load Libraries (Always Runs) ----

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(grid)      # grid.draw for heatmaps
  library(grDevices) # cairo_pdf / dev.off
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

plot_deseq2_qc <- function(dds        = NULL,
                           vst_blind  = NULL,
                           metadata   = NULL,
                           output_dir = NULL) {

  # ═══════════════════════════════════════════════════════════════════════════
  # SECTION 1: Validate Inputs
  # ═══════════════════════════════════════════════════════════════════════════
  # All structural checks fire here before any data is loaded or touched.
  # This follows the same contract as plot_volcano() and plot_heatmap() —
  # fail immediately with a clear message rather than crash deep inside a
  # plotting function with a cryptic error.

  # ── 1a. vst_blind and metadata must both be provided or both NULL ──────
  # vst_blind is only used for the pca, which also requires metadata
  # to derive sample groups. Providing one without the other is always a
  # mistake — fail early rather than silently skip the pca.
  if (xor(is.null(vst_blind), is.null(metadata))) {
    log_error(sample = "", step = "plot_deseq2_results",
              msg = paste("'vst_blind' and 'metadata' must both be provided or both NULL.",
                          "PCA requires both to group samples."))
  }

  # ═══════════════════════════════════════════════════════════════════════════
  # SECTION 2: Load Data
  # ═══════════════════════════════════════════════════════════════════════════

  dds      <- load_smart(input_path = dds)
  expr_mat <- load_smart(input_path = vst_blind)
  metadata <- load_smart(input_path = metadata)

  # ═══════════════════════════════════════════════════════════════════════════
  # SECTION 4: Set Up Output Path
  # ═══════════════════════════════════════════════════════════════════════════
  
  # Set global ggplot theme — plot_ma() and plot_volcano() apply
  # theme_publication() internally. theme_set() here ensures any additional
  # ggplot calls in this function inherit the same theme automatically.
  ggplot2::theme_set(theme_publication())

  output_path   <- file.path(output_dir)

  if (!dir.exists(output_path)) dir.create(output_path, recursive = TRUE)

  # ═══════════════════════════════════════════════════════════════════════════
  # SECTION 5: Dispersion Plot
  # ═══════════════════════════════════════════════════════════════════════════
  # res is a DESeqResults object from DESeq2::results(dds, contrast = ...).
  # plot_ma() accepts DESeqResults or data.frame — dds is not a valid input.
  plot_dispersion(dds        = dds,
                  output_dir = output_path)

  # ═══════════════════════════════════════════════════════════════════════════
  # SECTION 6: PCA Plot
  # ═══════════════════════════════════════════════════════════════════════════

  plot_pca(df         = expr_mat,
           metadata   = metadata,
           output_dir = output_path)

  # ═══════════════════════════════════════════════════════════════════════════
  # SECTION 7: Log and Return
  # ═══════════════════════════════════════════════════════════════════════════

  log_info(sample = "", step = "plot_deseq2_qc",
           msg    = glue::glue("DESeq2 QC plots saved to: '{output_path}'"))

  return(invisible(NULL))
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

  plot_deseq2_qc(
    dds             = get_arg(1),
    vst_blind       = get_arg(2),
    metadata        = get_arg(3),
    output_dir      = get_arg(4, ".")
  )
}
