# ---- 📦 1. Load Libraries (Always Runs) ----

suppressPackageStartupMessages({
  library(dplyr)
  library(glue)
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

subset_cells <- function(seurat_rds,
                         filter_col  = "Lineage",
                         filter_values,
                         label       = NULL,
                         output_dir  = NULL) {

  # ═══════════════════════════════════════════════════════════════════════════
  # SECTION 1: Load Seurat Object
  # ═══════════════════════════════════════════════════════════════════════════
  # filter_col drives which metadata column this call filters on — "Lineage" for
  # the Whole -> {label}_Initial subset (this pipeline's original/default use),
  # "CellType" for a {label}_Initial -> {label}_Final subset (second-pass
  # mixed-population removal: filter_values = "Retain" only). Any label's own
  # annotated object works as the source as long as filter_col is populated on it.

  seurat_obj <- load_smart(input_path = seurat_rds)

  log_info(sample = "", step = "subset_cells",
           msg = glue("Loaded: {ncol(seurat_obj)} cells, {nrow(seurat_obj)} features. ",
                      "Target label: '{label %||% 'unspecified'}'. Filter column: '{filter_col}'."))

  if (!filter_col %in% colnames(seurat_obj@meta.data)) {
    log_error(sample = "", step = "subset_cells",
              msg = glue("No {filter_col} column in metadata. Run annotation (auto or manual) first."))
  }

  # ═══════════════════════════════════════════════════════════════════════════
  # SECTION 2: Strip Non-RNA Assays and Reductions BEFORE Subsetting
  # ═══════════════════════════════════════════════════════════════════════════
  # MUST happen before base::subset() below, not after — subsetting cells
  # while multiple assays are present (especially SCT, which stores per-
  # object model parameters that do not reliably survive a cell-subset
  # operation) risks an error or a corrupted SCT model on the result.
  # Setting RNA as the only assay first makes subset() operate on plain
  # counts only — safe regardless of what the source object's SCT/other
  # assay state looked like.
  #
  # Reductions (pca, sct_pca, integ_harmony, umap_harmony, umap_harmony_qc,
  # ...) and neighbor graphs are dropped here too, for the same underlying
  # reason as the assay strip: they were fit using variance dominated by
  # between-lineage differences (PCA on "all animals" doesn't resolve cat
  # breeds) and are dead weight/stale baggage once this branch gets its own
  # fresh SCT_TRANSFORM_SUB -> INTEGRATE_SUB -> CLUSTER_AND_FINDMARKERS_SUB
  # pass — stripping them now avoids anyone downstream accidentally reading
  # or plotting a stale inherited embedding by mistake.

  SeuratObject::DefaultAssay(seurat_obj) <- "RNA"
  for (a in SeuratObject::Assays(seurat_obj)) {
    if (a != "RNA") seurat_obj[[a]] <- NULL
  }

  for (r in SeuratObject::Reductions(seurat_obj)) {
    seurat_obj[[r]] <- NULL
  }

  seurat_obj@graphs <- list()

  log_info(sample = "", step = "subset_cells",
           msg = "Stripped non-RNA assays, all reductions, all graphs before subsetting.")

  # ═══════════════════════════════════════════════════════════════════════════
  # SECTION 3: Validate and Apply Filter (on filter_col)
  # ═══════════════════════════════════════════════════════════════════════════
  # Two distinct failure modes, handled differently:
  #   - ALL requested values missing -> nothing to build -> hard error.
  #   - SOME requested values missing -> warn and proceed with whichever
  #     requested value(s) are actually present (e.g. a merged filter like
  #     c("T Cell","B Cell","Myeloid") where "Myeloid" genuinely isn't present
  #     in this tissue is a legitimate partial run, not a typo to block on).
  # matched_types (the intersection) is used for the actual subset — not the
  # original filter_values — so the code says what it means rather than
  # relying on %in%'s silent-drop behaviour to do the right thing implicitly.

  available_types <- sort(unique(seurat_obj[[filter_col]][, 1]))
  matched_types    <- intersect(filter_values, available_types)
  unknown_types    <- setdiff(filter_values, available_types)

  if (length(matched_types) == 0) {
    log_warn(sample = "", step = "subset_cells",
             msg = glue("SKIP: none of the requested filter_values ",
                       "({paste(filter_values, collapse=', ')}) are present in {filter_col} column ",
                       "for label '{label %||% 'unspecified'}'. Available: ",
                       "{paste(available_types, collapse=', ')}. No output written for this label ",
                       "(this is expected when a value is genuinely absent from this dataset)."))
    return(invisible(NULL))
  }

  if (length(unknown_types) > 0) {
    log_info(sample = "", step = "subset_cells",
             msg = glue("WARNING: {length(unknown_types)} requested {filter_col} value(s) not found ",
                       "and will be skipped: {paste(unknown_types, collapse=', ')}. ",
                       "Proceeding with: {paste(matched_types, collapse=', ')}."))
  }

  cells_keep     <- colnames(seurat_obj)[seurat_obj@meta.data[[filter_col]] %in% matched_types]
  
  # Capture barcodes about to be dropped
  dropped_barcodes <- setdiff(colnames(seurat_obj), cells_keep)
  
  new_excluded <- NULL
  if (length(dropped_barcodes) > 0) {
    
    # Context-aware reason assignment
    reasons <- if (filter_col == "Lineage") {
      # First Pass (Whole -> Lineage Initial): auto-junk/sparse cells are already
      # gone by this point (removed inside cluster_and_findmarkers() upstream) --
      # anything dropped here was a real, clustered cell whose Lineage call just
      # isn't one of the requested lineages for this branch.
      rep("Lineage_Excluded", length(dropped_barcodes))
      
    } else if (filter_col == "CellType") {
      # Second Pass ({label}_Initial -> {label}_Final or similar):
      # Cells dropped here during CellType filtering are contaminants/manual excludes
      raw_vals <- seurat_obj@meta.data[dropped_barcodes, filter_col]
      ifelse(raw_vals == "Exclude", "Contaminant", paste0("Dropped_", raw_vals))
      
    } else {
      # Generic fallback for any other custom filter columns
      rep(paste0("Dropped_", filter_col), length(dropped_barcodes))
    }
    
    new_excluded <- stats::setNames(reasons, dropped_barcodes)
  }

  # Carry forward any inherited ledger vector, then combine with this stage's drops
  prior_excluded    <- seurat_obj@misc$excluded_barcodes   # NULL on first run, or a named vector
  combined_excluded <- c(prior_excluded, new_excluded)

  seurat_subset  <- subset(seurat_obj, cells = cells_keep)

  log_info(sample = "", step = "subset_cells",
           msg = glue("Subset to {filter_col} in [{paste(matched_types, collapse=', ')}]: ",
                      "{ncol(seurat_subset)} / {ncol(seurat_obj)} cells retained ",
                      "({round(ncol(seurat_subset)/ncol(seurat_obj)*100, 1)}%)."))

  if (ncol(seurat_subset) == 0) {
    log_error(sample = "", step = "subset_cells",
              msg = "Subset contains 0 cells despite matched_types being non-empty — unexpected, check Lineage values above.")
  }

  rm(seurat_obj); gc()

  # ═══════════════════════════════════════════════════════════════════════════
  # SECTION 4: Strip Remaining Whole-Dataset-Derived Metadata/Provenance
  # ═══════════════════════════════════════════════════════════════════════════
  # Assays/reductions/graphs are already gone (Section 2, before subsetting).
  # What's left: metadata COLUMNS (Lineage, cluster assignments, module
  # scores, etc.) and @misc provenance from the Whole pipeline run this
  # object came from. Every label ("Whole", "Epithelial", ...) must enter
  # SCT_TRANSFORM with an identical contract: raw RNA counts + basic QC
  # metadata only, nothing computed.
  #
  # Column lists are read from @misc *_params BEFORE wiping @misc — reusing
  # the exact provenance each upstream step already recorded, rather than
  # guessing at column names via regex.

  cols_to_drop <- unique(c(
    seurat_subset@misc$clustering_params$cluster_cols,
    seurat_subset@misc$scoring_params$seurat_module_cols,
    seurat_subset@misc$scoring_params$ucell_module_cols,
    seurat_subset@misc$annotation_params$consensus_cols,
    grep("_MeanProb$|_Entropy$", colnames(seurat_subset@meta.data), value = TRUE),
    filter_col, "Stability_Score"
  ))
  cols_to_drop <- intersect(cols_to_drop, colnames(seurat_subset@meta.data))

  if (length(cols_to_drop) > 0) {
    seurat_subset@meta.data[, cols_to_drop] <- NULL
  }

  # @misc — wipe all upstream provenance, replace with subset provenance only. Dropping
  # filter_col here (whichever column this call filtered on — Lineage or CellType) means
  # the NEXT annotation pass on this object starts from a clean slate rather than
  # inheriting a stale Retain/Exclude or Lineage call from the object it was built from.
  seurat_subset@misc <- list(
    subset_params = list(
      source_rds        = seurat_rds,
      label             = label,
      filter_col        = filter_col,        # which column was filtered on
      filter_values     = filter_values,     # as requested
      matched_types     = matched_types,     # as actually applied
      unknown_types     = unknown_types,     # requested but not found, if any
      n_cells_subset    = ncol(seurat_subset)
    ),
    excluded_barcodes = combined_excluded
  )

  log_info(sample = "", step = "subset_cells",
           msg = glue("Stripped {length(cols_to_drop)} metadata column(s) and all @misc ",
                      "provenance except subset_params. Object now matches the same ",
                      "entry-point contract as MERGE_AND_PLOT_QC's filtered_seurat.rds."))

  # ═══════════════════════════════════════════════════════════════════════════
  # SECTION 5: Save
  # ═══════════════════════════════════════════════════════════════════════════
  # Filename is always "filtered_seurat.rds" regardless of label — identical
  # entry-point name across every label's folder (04.Seurat/{label}/). Which
  # folder it lands in is decided by the .nf process's publishDir, not here.

  if (!is.null(output_dir)) {
    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

    rds_out <- file.path(output_dir, "filtered_seurat.rds")
    saveRDS(object = seurat_subset, file = rds_out)

    log_info(sample = "", step = "subset_cells",
             msg = glue("Saved: '{rds_out}' ({ncol(seurat_subset)} cells)."))
  } else {
    log_info(sample = "", step = "subset_cells", msg = "No output_dir provided — skipping save.")
  }

  return(invisible(seurat_subset))
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

  # filter_values passed as comma-separated string, e.g. "Epithelial" or
  # "T Cell,B Cell,Myeloid" for a multi-type merged subset (filter_col="Lineage"),
  # or "Retain" for a second-pass mixed-population cleanup (filter_col="CellType")
  filter_values_raw <- get_arg(3)
  if (is.null(filter_values_raw)) {
    stop("❌ ERROR: filter_values (arg 3) is required, e.g. \"Epithelial\" or \"Retain\"")
  }
  filter_values_list <- trimws(strsplit(filter_values_raw, ",")[[1]])

  subset_cells(
    seurat_rds     = get_arg(1),
    filter_col     = get_arg(2, default = "Lineage"),
    filter_values  = filter_values_list,
    label          = get_arg(4),
    output_dir     = get_arg(5, default = ".")
  )
}
