# ---- 📦 1. Load Libraries (Always Runs) ----

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tibble)
  library(readxl)
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

run_sctransform <- function(seurat_rds, cellcycle_marker_xlsx, output_dir = NULL) {

  # Set seed for reproducible stochastic processes (PCA, UMAP)
  set.seed(1234)

  # ═══════════════════════════════════════════════════════════════════════════
  # SECTION 1: Load Seurat Object
  # ═══════════════════════════════════════════════════════════════════════════

  filtered_seurat <- load_smart(input_path = seurat_rds)

  # Automatically detects "RNA" or "Spatial"
  assay <- Seurat::DefaultAssay(filtered_seurat)
  is_spatial <- length(filtered_seurat@images) > 0

  log_info(sample = "", step = "run_sctransform",
           msg = glue::glue("Loaded: {ncol(filtered_seurat)} cells across ",
                            "{length(unique(filtered_seurat@meta.data$Sample))} samples."))

  # ═══════════════════════════════════════════════════════════════════════════
  # SECTION 2: Load Cell Cycle Genes
  # ═══════════════════════════════════════════════════════════════════════════
  # Pull both Mouse_Gene and Human_Gene columns for both phases.
  # Species filtering happens automatically via intersect with present features —
  # human genes won't exist in a mouse object and vice versa.

  cc_genes <- tryCatch({
    readxl::read_excel(cellcycle_marker_xlsx)
  }, error = function(e) {
    log_error(sample = "", step = "run_sctransform",
              msg = glue::glue("Failed to load cell cycle genes: '{cellcycle_marker_xlsx}' — {e$message}"))
  })

  # Validate required columns exist
  required_cols <- c("Phase", "Mouse_Gene", "Human_Gene")
  missing_cols  <- setdiff(required_cols, colnames(cc_genes))
  if (length(missing_cols) > 0) {
    log_error(sample = "", step = "run_sctransform",
              msg = glue::glue("Missing columns in cellcycle_marker_xlsx: ",
                               "{paste(missing_cols, collapse = ', ')}"))
  }

  # Pull S and G2M genes from both species columns — deduplicate ignoring case
  s_genes <- cc_genes %>%
    dplyr::filter(Phase == "S") %>%
    dplyr::select(Mouse_Gene, Human_Gene) %>%
    unlist() %>%
    na.omit() %>%
    unique()

  g2m_genes <- cc_genes %>%
    dplyr::filter(Phase == "G2/M") %>%
    dplyr::select(Mouse_Gene, Human_Gene) %>%
    unlist() %>%
    na.omit() %>%
    unique()

  # ═══════════════════════════════════════════════════════════════════════════
  # SECTION 3: Intersect Cell Cycle Genes with Present Features
  # ═══════════════════════════════════════════════════════════════════════════
  # SeuratObject::Features() returns union of all feature names at the assay
  # level. Intersect removes species-mismatched genes and any genes absent
  # from this dataset — prevents CellCycleScoring errors on missing features.

  present_features <- SeuratObject::Features(x = filtered_seurat,
                                             layer = "counts",
                                             simplify = TRUE)

  s_genes_all   <- s_genes    # retain full list for setdiff on variable features later
  g2m_genes_all <- g2m_genes

  s_genes   <- intersect(s_genes,   present_features)
  g2m_genes <- intersect(g2m_genes, present_features)

  if (length(s_genes) == 0) {
    log_error(sample = "", step = "run_sctransform",
              msg = glue::glue("Zero S-phase genes found in '{assay}' assay. ",
                               "Cannot perform cell cycle scoring."))
  }

  if (length(g2m_genes) == 0) {
    log_error(sample = "", step = "run_sctransform",
              msg = glue::glue("Zero G2M-phase genes found in '{assay}' assay. ",
                               "Cannot perform cell cycle scoring."))
  }

  log_info(sample = "", step = "run_sctransform",
           msg = glue::glue("{length(s_genes)} S-phase genes and ",
                            "{length(g2m_genes)} G2M-phase genes found and will be used for scoring."))

  # ═══════════════════════════════════════════════════════════════════════════
  # SECTION 4: Normalize Data (required for CellCycleScoring)
  # ═══════════════════════════════════════════════════════════════════════════
  # LogNormalize required before CellCycleScoring — scores are computed on
  # log-normalized data, not raw counts. SCTransform will renormalize
  # separately in Section 7.

  filtered_seurat <- Seurat::NormalizeData(object             = filtered_seurat,
                                           assay              = assay,
                                           normalization.method = "LogNormalize",
                                           scale.factor       = 10000,
                                           margin             = 1,
                                           verbose            = FALSE)

  # ═══════════════════════════════════════════════════════════════════════════
  # SECTION 5: Join Layers (required for CellCycleScoring)
  # ═══════════════════════════════════════════════════════════════════════════
  # CellCycleScoring requires a single data layer. Per-sample layers from
  # the merge step must be joined before scoring.

  filtered_seurat@assays[[assay]] <- SeuratObject::JoinLayers(
    object = filtered_seurat@assays[[assay]]
  )

  # ═══════════════════════════════════════════════════════════════════════════
  # SECTION 6: Score Cell Cycle Phases
  # ═══════════════════════════════════════════════════════════════════════════
  # Adds S.Score, G2M.Score, and Phase columns to @meta.data.
  # Scores are retained for annotation and visualization — NOT regressed out.
  # Cell cycle variance is handled downstream by removing cc genes from
  # variable features (Section 9), which is more transparent than regression.

  filtered_seurat <- Seurat::CellCycleScoring(object       = filtered_seurat,
                                              s.features   = s_genes,
                                              g2m.features = g2m_genes,
                                              ctrl         = NULL)

  # CC.Score = difference between G2M and S scores.
  # Positive values → G2M-dominant; negative → S-dominant; near zero → G1/quiescent.
  # Stored for annotation — not used for regression.
  filtered_seurat$CC.Score <- filtered_seurat$G2M.Score - filtered_seurat$S.Score

  log_info(sample = "", step = "run_sctransform",
           msg = "Cell cycle scoring complete.")

  # ═══════════════════════════════════════════════════════════════════════════
  # SECTION 7: Remove Small Samples Before SCTransform
  # ═══════════════════════════════════════════════════════════════════════════
  # SCTransform models each sample separately — too few cells produces
  # unreliable variance estimates. Minimum 50 cells matches the threshold
  # established in filter_singlets() for doublet detection reliability.
  # Removed samples still appear in QC plots from MERGE_FILTERED.

  cell_counts   <- table(filtered_seurat@meta.data$Sample)
  small_samples <- names(cell_counts[cell_counts < 50])

  if (length(small_samples) > 0) {
    log_warn(sample = "", step = "run_sctransform",
             msg = glue::glue("Removing {length(small_samples)} sample(s) with < 50 cells: ",
                              "{paste(small_samples, collapse = ', ')}"))
    
    # Capture barcodes about to be dropped as a named vector
    dropped_barcodes <- colnames(filtered_seurat)[filtered_seurat@meta.data$Sample %in% small_samples]
    new_excluded <- stats::setNames(rep("SCT_SampleDropped", length(dropped_barcodes)), dropped_barcodes)

    filtered_seurat <- subset(filtered_seurat,
                              subset = Sample %in% names(cell_counts[cell_counts >= 50]))
  } else {
    new_excluded <- NULL
  }
  
  log_info(sample = "", step = "run_sctransform",
           msg = glue::glue("{ncol(filtered_seurat)} cells remaining across ",
                            "{length(unique(filtered_seurat@meta.data$Sample))} samples."))

  # ═══════════════════════════════════════════════════════════════════════════
  # SECTION 8: Split Layers by Sample (batch-aware SCTransform)
  # ═══════════════════════════════════════════════════════════════════════════
  # SCTransform must model each batch (sample) separately to avoid confounding
  # biological variance with technical batch effects. All cells within the
  # same sample are analyzed together.

  if (!is_spatial) {
    filtered_seurat@assays[[assay]] <- base::split(
      x = filtered_seurat@assays[[assay]],
      f = filtered_seurat@meta.data[["Sample"]]
    )
  }

  # ═══════════════════════════════════════════════════════════════════════════
  # SECTION 9: Run SCTransform
  # ═══════════════════════════════════════════════════════════════════════════
  # Regress MitoRatio only — technical covariate for damaged/dying cells.
  # Cell cycle is handled via setdiff on variable features (Section 9).
  # SCTransform internally models sequencing depth via NB regression so
  # nUMIs and nGenes do not need to be regressed separately.
  # https://github.com/satijalab/seurat/issues/7342

  # Dynamically determine regression targets
  regress_targets <- NULL
  if ("MitoRatio" %in% colnames(filtered_seurat@meta.data)) {
    regress_targets <- "MitoRatio"
  }

  # future.globals.maxSize raised to handle large objects passed to parallel
  # workers — default 500MB is too small for datasets with 500K+ cells.
  # 500 * 1024^3 bytes ≈ 537GB — large enough for typical merged Seurat objects 
  # at this stage without disabling the safety check entirely.
  options(future.globals.maxSize = 500 * 1024^3)

  log_info(sample = "", step = "run_sctransform",
           msg = glue::glue("Running SCTransform on assay: '{assay}'."))
  sct_seurat <- Seurat::SCTransform(object    = filtered_seurat,
                        assay                 = assay,
                        new.assay.name        = "SCT",
                        do.correct.umi        = TRUE,
                        ncells                = 5000,
                        variable.features.n   = 3000,
                        vars.to.regress       = regress_targets,
                        do.scale              = FALSE,
                        do.center             = TRUE,
                        vst.flavor            = "v2",
                        return.only.var.genes = TRUE,
                        verbose               = TRUE)

  log_info(sample = "", step = "run_sctransform",
           msg = glue::glue("SCTransform complete. ",
                            "{length(Seurat::VariableFeatures(sct_seurat, assay = 'SCT'))} ",
                            "variable features selected."))

  # ═══════════════════════════════════════════════════════════════════════════
  # SECTION 10: Rejoin RNA Layers
  # ═══════════════════════════════════════════════════════════════════════════
  # RNA layers were split in Section 7 for SCTransform only. Rejoin immediately
  # after — no downstream step requires RNA to remain split.
  if (!is_spatial) {
    sct_seurat@assays[[assay]] <- SeuratObject::JoinLayers(
      object = sct_seurat@assays[[assay]]
    )
  }

  # ═══════════════════════════════════════════════════════════════════════════
  # SECTION 11: Prepare SCT Assay for Differential Expression
  # ═══════════════════════════════════════════════════════════════════════════
  # PrepSCTFindMarkers corrects the SCT counts across samples to account for
  # differences in sequencing depth — required ONLY if FindAllMarkers/FindMarkers
  # use SCT assay.

  # sct_seurat <- quiet_msg(
    # Seurat::PrepSCTFindMarkers(object  = sct_seurat,
                               # assay   = "SCT",
                               # verbose = FALSE)
  # )

  # ═══════════════════════════════════════════════════════════════════════════
  # SECTION 12: Filter Variable Features
  # ═══════════════════════════════════════════════════════════════════════════
  # Remove genes that should not drive PCA or clustering:
  #   Ribosomal   : ^Rps / ^Rpl — technical noise, already excluded from regression
  #   RIKEN       : ...Rik — mouse uncharacterized transcripts, no biological label
  #   Mitochondrial: ^Mt- — already regressed via MitoRatio, individual genes redundant
  #   Predicted   : ^Gm[0-9]+ — mouse predicted genes, uncharacterized
  #   Cell cycle  : s_genes + g2m_genes — removed via setdiff to prevent phase-driven
  #                 PCA axes while preserving scores in metadata for annotation

  var_f <- Seurat::VariableFeatures(sct_seurat, assay = "SCT")

  # Remove uninformative gene classes via regex
  var_f <- var_f[!grepl(
    pattern = "^[Rr][Pp][SsLl]|[Rr][Ii][Kk]$|^[Mm][Tt]-|^G[Mm][0-9.]+$",
    x       = var_f
  )]

  # Remove cell cycle genes — use full species-combined list before intersect
  # so that any cc gene present in the object is excluded regardless of which
  # species column it came from
  var_f <- setdiff(var_f, c(s_genes_all, g2m_genes_all))

  # Sort Variable Features Deterministically
  var_f <- sort(var_f)
  Seurat::VariableFeatures(sct_seurat, assay = "SCT") <- var_f

  log_info(sample = "", step = "run_sctransform",
           msg = glue::glue("Final variable features after filtering: {length(var_f)}."))

  # ═══════════════════════════════════════════════════════════════════════════
  # SECTION 13: Run PCA on SCT Assay
  # ═══════════════════════════════════════════════════════════════════════════
  # PCA on SCT residuals — used as input to IntegrateLayers in INTEGRATE_SCT.
  # reduction.key set explicitly for consistent naming across integration methods.

  sct_seurat <- Seurat::RunPCA(object         = sct_seurat,
                               assay          = "SCT",
                               features       = Seurat::VariableFeatures(sct_seurat, assay = "SCT"),
                               reduction.name = "sct_pca",
                               reduction.key  = "sctpca_",
                               verbose        = FALSE)

  log_info(sample = "", step = "run_sctransform",
           msg = "PCA on SCT assay complete.")

  # ═══════════════════════════════════════════════════════════════════════════
  # SECTION 14: Store SCT Parameters in @misc
  # ═══════════════════════════════════════════════════════════════════════════
  # Stored at object level — downstream steps read from here rather than
  # requiring these values to be re-passed as CLI arguments.

  sct_seurat@misc$sct_params <- list(
    assay               = assay,
    vars_to_regress     = "MitoRatio",
    vst_flavor          = "v2",
    n_variable_features = 3000,
    n_cells             = 5000,
    cc_genes_removed    = TRUE    # documents that setdiff was applied in Section 11
  )
  
  # Combine inherited ledger vector + this stage's new drops vector cleanly
  prior_excluded <- filtered_seurat@misc$excluded_barcodes # could be NULL, a vector, etc.
  sct_seurat@misc$excluded_barcodes <- c(prior_excluded, new_excluded)

  # ═══════════════════════════════════════════════════════════════════════════
  # SECTION 15: Save and Return
  # ═══════════════════════════════════════════════════════════════════════════

  if (!is.null(output_dir)) {

    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

    rds_out <- file.path(output_dir, "sct_seurat.rds")
    saveRDS(object = sct_seurat, file = rds_out)

    log_info(sample = "", step = "run_sctransform",
             msg = glue::glue("SCT Seurat object saved to: '{rds_out}'"))

  } else {
    log_info(sample = "", step = "run_sctransform",
             msg = "No output_dir provided — skipping save.")
  }

  log_info(sample = "", step = "run_sctransform",
           msg = "SCTransform completed successfully.")

  return(invisible(sct_seurat))
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

  run_sctransform(
    seurat_rds            = get_arg(1),
    cellcycle_marker_xlsx = get_arg(2),
    output_dir            = get_arg(3, default = ".")
  )
}
