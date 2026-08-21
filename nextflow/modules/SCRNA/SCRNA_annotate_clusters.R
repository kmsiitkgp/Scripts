# ---- 📦 1. Load Libraries (Always Runs) ----

suppressPackageStartupMessages({
  library(dplyr)
  library(glue)
  library(openxlsx)
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

# =============================================================================
# FUNCTION: annotate_clusters
# =============================================================================
# Manual cluster annotation -- the only annotation path (no automatic path
# exists any more). Which metadata column gets written is driven by `label`:
#   - label == "Whole"   -> writes/validates/logs against "Lineage"
#   - any other label    -> writes/validates/logs against "CellType"
# This is a structural distinction, not cosmetic: for "Whole", the object
# ends up with a Lineage column (and NO CellType column at all -- see
# Section 3); for a subset branch, it ends up with a CellType column, exactly
# as before. SCRNA_ADD_CELLTYPE reads $Lineage directly off the Whole object
# and $CellType off each subset object -- no renaming happens there any more.
# Candidate-marker-for-new-cell-types generation has moved to
# SCRNA_ADD_CELLTYPE (runs once, on the final assembled object) and is no
# longer part of this script.
# =============================================================================

annotate_clusters <- function(seurat_rds,
                              label,
                              annotation_xlsx,
                              cluster_col,
                              exclude_label = "Exclude",
                              output_dir    = NULL) {
  
  # target_col: the metadata column this run annotates. "Whole" gets the
  # coarse Lineage call; every other label gets the fine-grained CellType
  # call. Used throughout for validation, logging, plotting, and the final
  # metadata column written.
  target_col <- if (identical(label, "Whole")) "Lineage" else "CellType"
  
  # ═══════════════════════════════════════════════════════════════════════════
  # SECTION 1: Load Seurat Object
  # ═══════════════════════════════════════════════════════════════════════════
  # Expects clustered_seurat.rds produced by CLUSTER_AND_FINDMARKERS — junk/
  # sparse cells already physically removed upstream, so every barcode here
  # needs a template label.
  
  seurat_obj <- load_smart(input_path = seurat_rds)
  
  log_info(sample = "", step = "annotate_clusters",
           msg = glue("Loaded: {ncol(seurat_obj)} cells, ",
                      "{nrow(seurat_obj)} features."))
  
  if (!cluster_col %in% colnames(seurat_obj@meta.data)) {
    log_error(sample = "", step = "annotate_clusters",
              msg = glue("cluster_col '{cluster_col}' not found in metadata. ",
                         "Available: {paste(grep('^harmony_res', colnames(seurat_obj@meta.data), value = TRUE), collapse=', ')}"))
  }
  
  # ═══════════════════════════════════════════════════════════════════════════
  # SECTION 2: Load and Validate Filled-In Annotation File
  # ═══════════════════════════════════════════════════════════════════════════
  # Reads the "Annotation" sheet of the template produced by
  # generate_annotation_template() (in CLUSTER_AND_FINDMARKERS), after the
  # user has filled in the Lineage (Whole) or CellType (subset) column by
  # hand.
  #
  # Two-directional validation:
  #   1. Every cluster PRESENT IN THE OBJECT must have a non-blank label
  #      in the file -> hard error listing exactly which clusters are missing.
  #      A silently-NA label downstream is a much worse failure mode than
  #      stopping here.
  #   2. Any cluster_id in the file but NOT in the object (stray leftover
  #      row, e.g. copy-paste error) -> warning only, ignored.
  
  if (!file.exists(annotation_xlsx)) {
    log_error(sample = "", step = "annotate_clusters",
              msg = glue("annotation_xlsx not found: '{annotation_xlsx}'"))
  }
  
  anno_df <- openxlsx::read.xlsx(annotation_xlsx, sheet = "Annotation")
  
  # Check if optional SubType column is present
  has_subtype <- "SubType" %in% colnames(anno_df)
  
  required_cols <- c("cluster_id", target_col)
  missing_cols  <- setdiff(required_cols, colnames(anno_df))
  if (length(missing_cols) > 0) {
    log_error(sample = "", step = "annotate_clusters",
              msg = glue("annotation_xlsx missing required column(s): ",
                         "{paste(missing_cols, collapse=', ')}"))
  }
  
  anno_df$cluster_id           <- as.character(anno_df$cluster_id)
  anno_df[[target_col]]        <- trimws(as.character(anno_df[[target_col]]))
  if (has_subtype) {
    anno_df$SubType <- trimws(as.character(anno_df$SubType))
    # Replace empty strings or NA with NA or a placeholder if needed
    anno_df$SubType[anno_df$SubType == ""] <- NA
  }
  
  object_clusters <- unique(as.character(seurat_obj@meta.data[[cluster_col]]))
  file_clusters    <- unique(anno_df$cluster_id)
  
  # Direction 1: object clusters missing a label, or blank/NA label
  labeled_lookup   <- stats::setNames(anno_df[[target_col]], anno_df$cluster_id)
  # Optional SubType lookup setup
  if (has_subtype) {
    subtype_lookup <- stats::setNames(anno_df$SubType, anno_df$cluster_id)
  }
  missing_labels   <- object_clusters[
    is.na(labeled_lookup[object_clusters]) | labeled_lookup[object_clusters] == ""
  ]
  if (length(missing_labels) > 0) {
    log_error(sample = "", step = "annotate_clusters",
              msg = glue("{length(missing_labels)} cluster(s) in '{cluster_col}' have no ",
                         "{target_col} in annotation_xlsx: {paste(sort(missing_labels), collapse=', ')}. ",
                         "Fill in every row of the Annotation sheet before rerunning."))
  }
  
  # Direction 2: stray rows in file not present in the object
  stray_rows <- setdiff(file_clusters, object_clusters)
  if (length(stray_rows) > 0) {
    log_info(sample = "", step = "annotate_clusters",
             msg = glue("WARNING: {length(stray_rows)} cluster_id(s) in annotation_xlsx ",
                        "not present in '{cluster_col}' — ignored: {paste(stray_rows, collapse=', ')}"))
  }
  
  log_info(sample = "", step = "annotate_clusters",
           msg = glue("Validated: all {length(object_clusters)} clusters in '{cluster_col}' ",
                      "have a {target_col} label."))
  
  # ═══════════════════════════════════════════════════════════════════════════
  # SECTION 3: Map Label onto Cells
  # ═══════════════════════════════════════════════════════════════════════════
  # Junk/sparse cells are physically removed upstream by cluster_and_findmarkers()
  # before this object is ever saved, so there's no pre-filled auto-Exclude to
  # preserve here -- every barcode present just takes its cluster's template label.
  # FINAL label goes to target_col ("Lineage" for Whole, "CellType" for subsets).
  
  # No pre-flagged auto-Exclude barcodes to preserve anymore -- junk/sparse cells
  # are physically removed by cluster_and_findmarkers() before this object is ever
  # saved, so every barcode present here just takes its cluster's template label.
  cell_clusters   <- as.character(seurat_obj@meta.data[[cluster_col]])
  final_labels    <- unname(labeled_lookup[cell_clusters])
  
  seurat_obj@meta.data[[target_col]] <- final_labels
  
  # Map SubType if present in the spreadsheet
  if (has_subtype) {
    final_subtypes <- unname(subtype_lookup[cell_clusters])
    seurat_obj@meta.data$SubType <- final_subtypes
    
    log_info(sample = "", step = "annotate_clusters",
             msg = glue("{ncol(seurat_obj)} barcode(s) assigned SubType from template."))
  }
  
  # Full confidence by construction — no consensus vote to derive a graded
  # score from, unlike an automatic path. This applies equally to
  # template-assigned cells AND pre-flagged auto-Exclude cells: the auto-flag
  # is barcode-level, cross-resolution-consistent evidence (see
  # identify_junk_and_sparse's intersect logic) — at least as reliable as a
  # single cluster-wide manual label, not less.
  seurat_obj$Stability_Score <- 1
  
  # ═══════════════════════════════════════════════════════════════════════════
  # SECTION 4: Handle Excluded Clusters
  # ═══════════════════════════════════════════════════════════════════════════
  # Clusters labeled with exclude_label (default "Exclude") are junk/mixed
  # populations the user identified by eye that automatic junk-detection
  # missed — on top of any individual barcodes already pre-flagged Exclude
  # by CLUSTER_AND_FINDMARKERS's automated detection (Section 3), which may
  # sit inside a cluster the template otherwise labels as a real cell type.
  
  is_excluded             <- seurat_obj@meta.data[[target_col]] == exclude_label
  n_excluded_cells        <- sum(is_excluded)
  # Clusters the TEMPLATE itself labeled Exclude (a coarser signal than
  # is_excluded, which also includes individual auto-flagged barcodes sitting
  # in clusters the template labeled as a real cell type).
  template_excluded_clusters <- sort(unique(names(labeled_lookup)[labeled_lookup == exclude_label]))
  
  log_info(sample = "", step = "annotate_clusters",
           msg = glue("Exclude total: {n_excluded_cells} cells ({round(mean(is_excluded)*100, 1)}%) ",
                      "across {length(template_excluded_clusters)} template-excluded cluster(s) — ",
                      "all stay in this object; physical removal happens downstream in SCRNA_ADD_CELLTYPE."))
  
  # ═══════════════════════════════════════════════════════════════════════════
  # SECTION 5: Store Annotation Parameters in @misc
  # ═══════════════════════════════════════════════════════════════════════════
  # celltype_levels excludes exclude_label — downstream consumers should
  # never see "Exclude" treated as a real cell type.
  
  used_types <- sort(unique(anno_df[[target_col]][anno_df[[target_col]] != exclude_label & anno_df[[target_col]] != ""]))
  
  seurat_obj@misc$annotation_params <- list(
    method            = "manual",
    label_col         = target_col,
    cluster_col       = cluster_col,
    annotation_xlsx   = annotation_xlsx,
    exclude_label     = exclude_label,
    excluded_clusters = template_excluded_clusters,
    celltype_levels   = used_types
  )
  
  # ═══════════════════════════════════════════════════════════════════════════
  # SECTION 6: UMAP Plot
  # ═══════════════════════════════════════════════════════════════════════════
  # All cells, including Exclude-labeled ones — lets the user visually
  # confirm the excluded cells are where they expect before they're
  # eventually removed downstream in SCRNA_ADD_CELLTYPE.
  
  if (!is.null(output_dir)) {
    plot_seurat(seurat_obj,
                features     = target_col,
                feature_type = "metadata",
                reduction    = "umap_harmony",
                filename     = glue("UMAP_Annotation_{target_col}"),
                output_dir   = output_dir)

    if ("SubType" %in% colnames(seurat_obj@meta.data)) {
      plot_seurat(seurat_obj,
                  features     = "SubType",
                  feature_type = "metadata",
                  reduction    = "umap_harmony",
                  filename     = "UMAP_Annotation_SubType",
                  output_dir   = output_dir)
    }

    log_info(sample = "", step = "annotate_clusters",
             msg = "Annotation UMAP plot saved.")
  }
  
  # ═══════════════════════════════════════════════════════════════════════════
  # SECTION 7: Save Seurat Object
  # ═══════════════════════════════════════════════════════════════════════════
  # ONE object saved — all cells, including every Exclude-labeled one
  # (manual template calls and auto-flagged barcodes alike). No clean
  # variant generated here: CellType isn't truly final until SCRNA_ADD_CELLTYPE
  # has assembled Whole + every subset branch's own CellType calls (a
  # subset's fine-grained CellType, including its own Exclude calls,
  # overrides Whole's coarse Lineage for its barcodes) — filtering out
  # Exclude before that point would be premature. Clean is computed exactly
  # once, there, after CellType is genuinely final.
  
  if (!is.null(output_dir)) {
    
    rds_out <- file.path(output_dir, "annotated_seurat.rds")
    saveRDS(object = seurat_obj, file = rds_out)
    log_info(sample = "", step = "annotate_clusters",
             msg = glue("Annotated object saved: '{rds_out}' ",
                        "({ncol(seurat_obj)} cells, includes Exclude-labeled cells)."))
    
    # -------------------------------------------------------------------------
    # Lightweight barcode/metadata CSV — written on EVERY call (Whole,
    # {label}_Initial, {label}_Final alike), not gated on label, so this stays
    # one simple unconditional block rather than special-cased per branch.
    # Only the "_Final" branches' CSVs are actually consumed downstream (by
    # SCRNA_PREP_VELOCITY_INPUT, via lineage_ch) — Whole/_Initial CSVs just sit
    # unused in their own folders, same pattern as *.pdf being emitted
    # unconditionally regardless of whether every label's plots get looked at.
    #
    # target_col is written generically as "CellType" in the CSV output
    # regardless of whether the source column was actually "Lineage" (Whole)
    # or "CellType" (every other label) — downstream Python code expects one
    # consistent column name, and doesn't care which pass produced it.
    
    #Filter out excluded cells for the exported CSV if desired
    barcodes_to_keep <- seurat_obj@meta.data[[target_col]] != exclude_label
    n_total_csv      <- ncol(seurat_obj)
    n_kept_csv       <- sum(barcodes_to_keep)
    n_dropped_csv    <- n_total_csv - n_kept_csv
    
    csv_out <- file.path(output_dir, "barcodes.csv")
    barcodes_df <- data.frame(
      barcode  = colnames(seurat_obj)[barcodes_to_keep],
      CellType = seurat_obj@meta.data[[target_col]][barcodes_to_keep],
      cluster  = seurat_obj@meta.data[[cluster_col]][barcodes_to_keep],
      row.names = NULL
    )
    if ("umap_harmony" %in% names(seurat_obj@reductions)) {
      umap_coords <- Embeddings(seurat_obj, "umap_harmony")[barcodes_to_keep, , drop = FALSE]
      barcodes_df$UMAP_1 <- umap_coords[, 1]
      barcodes_df$UMAP_2 <- umap_coords[, 2]
    } else {
      log_info(sample = "", step = "annotate_clusters",
               msg = "No 'umap_harmony' reduction found — barcodes.csv written without UMAP columns.")
    }
    write.csv(barcodes_df, csv_out, row.names = FALSE)
    log_info(sample = "", step = "annotate_clusters",
             msg = glue("Metadata CSV saved: '{csv_out}' ({n_kept_csv} barcodes retained, ",
                        "{n_dropped_csv} '{exclude_label}' cells excluded from CSV)."))
    
  } else {
    log_info(sample = "", step = "annotate_clusters",
             msg = "No output_dir provided — skipping save.")
  }
  
  # ═══════════════════════════════════════════════════════════════════════════
  # SECTION 8: Log Summary and Return
  # ═══════════════════════════════════════════════════════════════════════════
  
  label_table <- sort(table(seurat_obj@meta.data[[target_col]]), decreasing = TRUE)
  log_info(sample = "", step = "annotate_clusters", msg = glue("{target_col} distribution:"))
  for (ct in names(label_table)) {
    pct <- round(label_table[ct] / ncol(seurat_obj) * 100, 1)
    log_info(sample = "", step = "annotate_clusters",
             msg = glue("  {ct}: {label_table[ct]} cells ({pct}%)"))
  }
  
  log_info(sample = "", step = "annotate_clusters",
           msg = glue("Manual annotation completed. New metadata: {target_col}, Stability_Score."))
  
  return(invisible(seurat_obj))
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
  
  annotate_clusters(
    seurat_rds      = get_arg(1),
    label           = get_arg(2),
    annotation_xlsx = get_arg(3),
    cluster_col     = get_arg(4),
    exclude_label   = get_arg(5, default = "Exclude"),
    output_dir      = get_arg(6, default = ".")
  )
}
