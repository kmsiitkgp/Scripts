# ---- 📦 1. Load Libraries (Always Runs) ----

suppressPackageStartupMessages({
  library(dplyr)
  library(glue)
  library(ggplot2)
  library(Seurat)
  library(SeuratObject)
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
  modules_dir     <- dirname(sub("--file=", "", file_arg))
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

integrate <- function(seurat_rds, output_dir = NULL) {

  set.seed(1234)

  # ═══════════════════════════════════════════════════════════════════════════
  # SECTION 1: Load Seurat Object
  # ═══════════════════════════════════════════════════════════════════════════

  sct_seurat <- load_smart(input_path = seurat_rds)

  log_info(sample = "", step = "INTEGRATE",
           msg = glue("Loaded: {ncol(sct_seurat)} cells across ",
                      "{length(unique(sct_seurat@meta.data$Sample))} samples."))

  # ═══════════════════════════════════════════════════════════════════════════
  # SECTION 2: Validate Prerequisites
  # ═══════════════════════════════════════════════════════════════════════════
  # sct_pca must exist — Harmony operates on PCA embeddings, not raw counts.
  # SCT assay must exist — integration always runs on SCT, never raw RNA.
  # RNA assay must exist — required by FindAllMarkers in downstream steps.

  if (!"sct_pca" %in% names(sct_seurat@reductions)) {
    log_error(sample = "", step = "INTEGRATE",
              msg = "Reduction 'sct_pca' missing. Run SCT_TRANSFORM first.")
  }

  if (!"SCT" %in% Seurat::Assays(sct_seurat)) {
    log_error(sample = "", step = "INTEGRATE",
              msg = "SCT assay missing. Run SCT_TRANSFORM first.")
  }

  if (!"RNA" %in% Seurat::Assays(sct_seurat)) {
    log_error(sample = "", step = "INTEGRATE",
              msg = "RNA assay missing — required for FindAllMarkers downstream.")
  }

  # ═══════════════════════════════════════════════════════════════════════════
  # SECTION 3: Set Parameters
  # ═══════════════════════════════════════════════════════════════════════════
  #
  # INTEGRATION METHOD — Harmony only:
  #   Harmony corrects batch effects at the PCA embedding level — operates on
  #   sct_pca (n_cells × n_dims), which is tiny regardless of cell count.
  #   This makes it memory-safe at 500k+ cells where anchor-based methods fail.
  #
  #   RPCA was dropped: anchor-based integration builds pairwise cross-sample
  #   matrices that cause R long vector crashes (> 2^31 elements) at this scale,
  #   even with sketch-based cell count reduction. Reference-based RPCA was
  #   explored but layer ordering after split() is not guaranteed in Seurat v5,
  #   making safe index mapping unreliable.
  #
  #   CCA was excluded: known to over-integrate heterogeneous datasets,
  #   collapsing real biological variation between conditions.
  #
  #   JointPCA was excluded: less battle-tested than Harmony in the literature.
  #
  #   scVI was excluded: requires Python/conda and GPU for meaningful speed
  #   advantage. Marginal accuracy gain does not justify pipeline complexity.
  #
  # n_dims = 30:
  #   Standard for datasets of this scale. Captures major axes of variation
  #   without including noise-dominated higher PCs.

  n_dims   <- 30
  n_cells  <- ncol(sct_seurat)
  max_dims <- min(n_dims, ncol(Seurat::Embeddings(sct_seurat, reduction = "sct_pca")))

  log_info(sample = "", step = "INTEGRATE",
           msg = glue("Total cells: {n_cells}. Using {max_dims} PCA dims for Harmony."))

  # ═══════════════════════════════════════════════════════════════════════════
  # SECTION 4: Harmony Integration
  # ═══════════════════════════════════════════════════════════════════════════
  # group.by.vars = "Sample" targets sample-of-origin technical variation only.
  # Biological variables (Condition, Age, etc.) are intentionally not corrected.
  # IntegrateLayers is the Seurat v5 interface — handles SCT correctly without
  # requiring manual layer splitting.
  #
  # We extract only the integ_harmony reduction from the result and attach it
  # to our object — this avoids carrying duplicate large objects in memory.

  integrated_seurat <- sct_seurat
  Seurat::DefaultAssay(integrated_seurat) <- "SCT"

  # 50 * 1024^3 bytes ≈ 53.7GB — large enough for typical merged Seurat objects 
  # at this stage without disabling the safety check entirely.
  options(future.globals.maxSize = 50 * 1024^3)

  log_info(sample = "", step = "INTEGRATE",
           msg = glue("Running Harmony on {n_cells} cells → 'integ_harmony'."))

  harmony_result <- Seurat::IntegrateLayers(
    object               = integrated_seurat,
    method               = "HarmonyIntegration",
    assay                = "SCT",
    normalization.method = "SCT",
    orig.reduction       = "sct_pca",
    new.reduction        = "integ_harmony",
    group.by.vars        = "Sample",
    dims                 = 1:max_dims,
    verbose              = FALSE
  )

  integrated_seurat[["integ_harmony"]] <- harmony_result[["integ_harmony"]]
  rm(harmony_result); gc()

  log_info(sample = "", step = "INTEGRATE",
           msg = "Harmony complete → 'integ_harmony' added.")

  # ═══════════════════════════════════════════════════════════════════════════
  # SECTION 5: Post-Integration UMAP Plots
  # ═══════════════════════════════════════════════════════════════════════════
  # We run a quick UMAP here on the Harmony embedding purely for QC visualization
  # — to verify that Harmony mixed samples correctly without collapsing
  # biological variation. This UMAP is NOT used downstream; CLUSTER_AND_FINDMARKERS
  # runs its own UMAP after clustering on the clean dataset.
  #
  # Plots saved:
  #   UMAP_Sample  — cells colored by sample, split by sample
  #                  Good integration: samples overlap. Bad: samples separate.
  #   UMAP_Phase   — cells colored by cell cycle phase
  #                  Verifies cell cycle regression worked in SCT_TRANSFORM.
  #   UMAP_QC      — nUMIs, nGenes, MitoRatio, cell cycle scores
  #                  Verifies no technical gradients remain after integration.

  if (!is.null(output_dir)) {

    log_info(sample = "", step = "INTEGRATE",
             msg = "Running QC UMAP for integration validation.")

    integrated_seurat <- Seurat::RunUMAP(
      object         = integrated_seurat,
      reduction      = "integ_harmony",
      dims           = 1:max_dims,
      n.neighbors    = 30L,
      reduction.name = "umap_harmony_qc",
      return.model   = FALSE,
      verbose        = FALSE
    )

    plot_seurat(integrated_seurat,
                features     = "Sample",
                feature_type = "metadata",
                reduction    = "umap_harmony_qc",
                filename     = "UMAP_Integration_Sample",
                split_col    = "Sample",
                output_dir   = output_dir)

    plot_seurat(integrated_seurat,
                features     = "Phase",
                feature_type = "metadata",
                reduction    = "umap_harmony_qc",
                filename     = "UMAP_Integration_Phase",
                split_col    = "Phase",
                output_dir   = output_dir)

    plot_seurat(integrated_seurat,
                features     = c("nUMIs", "nGenes", "MitoRatio",
                                 "S.Score", "G2M.Score", "CC.Score"),
                feature_type = "metadata",
                reduction    = "umap_harmony_qc",
                filename     = "UMAP_Integration_QC_Metrics",
                output_dir   = output_dir)

    log_info(sample = "", step = "INTEGRATE",
             msg = "Integration QC plots saved.")
  }

  # ═══════════════════════════════════════════════════════════════════════════
  # SECTION 6: Store Parameters in @misc
  # ═══════════════════════════════════════════════════════════════════════════
  # Downstream modules (CLUSTER_AND_FINDMARKERS, CALC_MODULE_SCORES,
  # ANNOTATE_CLUSTERS) read integration_params from @misc to access the
  # harmony reduction name and dims without hardcoding them.

  integrated_seurat@misc$integration_params <- list(
    method               = "Harmony",
    n_dims               = max_dims,
    normalization_method = "SCT",
    group_by_vars        = "Sample",
    harmony_reduction    = "integ_harmony"
  )

  # ═══════════════════════════════════════════════════════════════════════════
  # SECTION 7: Save and Return
  # ═══════════════════════════════════════════════════════════════════════════

  if (!is.null(output_dir)) {

    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

    rds_out <- file.path(output_dir, "integrated_seurat.rds")
    saveRDS(object = integrated_seurat, file = rds_out)

    log_info(sample = "", step = "INTEGRATE",
             msg = glue("Saved: '{rds_out}' ({ncol(integrated_seurat)} cells)."))

  } else {
    log_info(sample = "", step = "INTEGRATE",
             msg = "No output_dir provided — skipping save.")
  }

  log_info(sample = "", step = "INTEGRATE",
           msg = "Integration completed successfully.")

  return(invisible(integrated_seurat))
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

  integrate(
    seurat_rds = get_arg(1),
    output_dir = get_arg(2, default = ".")
  )
}
